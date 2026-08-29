import ROOT
import subprocess
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as mh
from datacard.ggHdatacardmaker import main as make_datacard
from ggHparameters import signal_path, bkg_path, lumi, order_fit


era_years = {
    "Run2": ["2017", "2018"],
    "Run3": ["2022", "2023", "2024"],
}
era_years["Run2Run3"] = era_years["Run2"] + era_years["Run3"]
sig_eras = {
    "2022": ["2022preEE", "2022postEE"],
    "2023": ["2023preBPix", "2023postBPix"],
    "2024": ["2024"],
}
era_aliases = {"run2": "Run2", "run3": "Run3",
                "run2run3": "Run2Run3", "run3run2": "Run2Run3",
                "run2+run3": "Run2Run3"}
data_years = {"2017", "2018", "2022", "2023", "2024"}


def resolve_era(token):
    key = era_aliases.get(token.strip().lower())
    if key:
        return key, list(era_years[key])
    if token in data_years:
        return token, [token]
    raise ValueError(f"unsupported year or era '{token}'")


def era_lumi(years):
    return sum(lumi[y] for y in years)


quantile_keys = {0.025: "exp-2", 0.16: "exp-1",0.5: "exp0", 0.84: "exp+1", 0.975: "exp+2"}


def get_limits(root_file):
    f = ROOT.TFile.Open(root_file)
    tree = f.Get("limit") if f else None
    out = {}
    if not tree:
        print(f"[get_limits] WARNING: no 'limit' tree in {root_file}; skipping point.")
        if f:
            f.Close()
        return out
    for entry in tree:
        for q, key in quantile_keys.items():
            if abs(entry.quantileExpected - q) < 1e-3:
                out[key] = entry.limit
    f.Close()
    return out


def save_results(results, filename):
    with open(filename, "w") as fout:
        json.dump(results, fout, indent=2)
    print(f"[save_results] wrote {filename}")


def load_results(filename):
    with open(filename) as fin:
        return json.load(fin)


def combine_cards(labelled_cards, out_txt):
    args = [f"{label}={card}" for label, card in labelled_cards]
    with open(out_txt, "w") as f:
        subprocess.run(["combineCards.py", *args], stdout=f, check=True)
    return out_txt


def scan_mass_lifetime(masses, lifetimes, era, years, bins, finalstate="4g", physics="ggH", order_fit=order_fit, results_json=None):
    ROOT.gROOT.SetBatch(True)
    if results_json is None:
        results_json = f"limits_UL_vs_mass_{era}.json"
    results = {era: {}}
    for mass in masses:
        for ctau in lifetimes:
            labelled_cards = []
            for year in years:
                signal_years = sig_eras.get(year, [year])
                signals = [signal_path(mass, ctau, signal_year) for signal_year in signal_years]
                bkg = bkg_path(year)
                make_datacard(paths=signals+[bkg], isMC=[1]*len(signals)+[0], trees=["ggH4g"]*(len(signals)+1), var=f"best_4g_corr_mass_m{mass}", period=year, bins=bins, lifetime=ctau, mass=mass, finalstate=finalstate, physics=physics, order_fit=order_fit, signal_lumis=[lumi[signal_year] for signal_year in signal_years])
                card = f"datacard_{physics}_{finalstate}_m{mass}_ct{ctau}_{year}.txt"
                labelled_cards.append((year, card))

            combined_txt = f"datacard_{physics}_{finalstate}_m{mass}_ct{ctau}_{era}.txt"
            combined_root = f"datacard_{physics}_{finalstate}_m{mass}_ct{ctau}_{era}.root"
            combine_cards(labelled_cards, combined_txt)
            subprocess.run(["text2workspace.py", combined_txt, "-o", combined_root], check=True)
            subprocess.run(["combine", "-M", "AsymptoticLimits", combined_root, "-m", "125", "--run", "blind"], check=True)

            result_root = f"higgsCombineTest_m{mass}_ct{ctau}_{era}.AsymptoticLimits.mH125.root"
            subprocess.run(["mv", "higgsCombineTest.AsymptoticLimits.mH125.root", result_root], check=True)

            lims = get_limits(result_root)
            if lims:
                results[era].setdefault(str(ctau), {})[str(mass)] = lims
            print(f"m={mass} ct={ctau} era={era}: expected UL on r = {lims.get('exp0')}")

    save_results(results, results_json)
    return results


def panel_arrays(panel):
    masses = sorted(int(m) for m in panel)
    needed = ["exp-2", "exp-1", "exp0", "exp+1", "exp+2"]
    masses = [m for m in masses if all(k in panel[str(m)] for k in needed)]
    arrays = {k: [panel[str(m)][k] for m in masses] for k in needed}
    return masses, arrays


def plot_UL_vs_mass(results, era, ctaus, total_lumi=None, br_scale=1.0, ytitle=r"95% CL Upper Limit on $\mathcal{B}(H\rightarrow\Phi\Phi)$", extra_labels=(r"$\mathcal{B}(\Phi\rightarrow\gamma\gamma)=1$",), yrange=None, outname=None, center_of_mass=13):
    ROOT.gROOT.SetBatch(True)
    mh.style.use("CMS")

    if total_lumi is None:
        total_lumi = lumi.get(era)
    year_res = results[era] if era in results else results
    ctaus = [ct for ct in ctaus if str(ct) in year_res]
    n = len(ctaus)
    if n == 0:
        print("[plot_UL_vs_mass] no lifetimes with data to plot.")
        return

    panels = []
    for ct in ctaus:
        masses, arrays = panel_arrays(year_res[str(ct)])
        arrays = {k: [v * br_scale for v in vals] for k, vals in arrays.items()}
        panels.append((ct, masses, arrays))

    if not any(masses for _, masses, _ in panels):
        print("[plot_UL_vs_mass] no complete points to plot.")
        return

    band2sigma_color = "#85d2fb"
    band1sigma_color = "#ffde9c"
    fig, axes = plt.subplots(1, n, sharey=True, squeeze=False, figsize=(7.5 * n, 20))
    axes = axes[0]
    plt.subplots_adjust(wspace=0)

    for i, (ct, masses, arrays) in enumerate(panels):
        ax = axes[i]
        if masses:
            ax.fill_between(masses, arrays["exp-2"], arrays["exp+2"], color=band2sigma_color, label=r"$\pm 2\sigma$")
            ax.fill_between(masses, arrays["exp-1"], arrays["exp+1"], color=band1sigma_color, label=r"$\pm 1\sigma$")
            ax.plot(masses, arrays["exp0"], linestyle="--", color="black", linewidth=4, label="Expected")

        ax.text(0.5, 0.02, rf"$c\tau = {ct}\ \mathrm{{mm}}$", transform=ax.transAxes, fontsize=40, ha="center", va="center", bbox=dict(facecolor="white", edgecolor="none", alpha=0.5))
        ax.tick_params(axis="both", which="major", labelsize=40)
        ax.margins(x=0)
        ax.minorticks_on()

    axes[0].legend(loc="upper right", fontsize=40)
    axes[-1].set_xlabel(r"$m_\Phi$ (GeV)", fontsize=40)
    axes[0].set_ylabel(ytitle, fontsize=40)
    axes[0].set_yscale("log")
    if yrange is not None:
        axes[0].set_ylim(*yrange)

    mh.cms.label("Work in Progress", data=True, rlabel="", ax=axes[0], loc=0, fontsize=60)
    mh.cms.label(None, exp="", data=True, llabel="", ax=axes[-1], loc=0, lumi=f"{total_lumi / 1000:.2f}", com=center_of_mass, fontsize=40)
    for j, label in enumerate(extra_labels):
        axes[0].text(0.5, 0.70 - 0.06 * j, label, transform=axes[0].transAxes, fontsize=40, ha="center", va="center", bbox=dict(facecolor="white", edgecolor="none", alpha=0.5))

    if outname is None:
        outname = f"UL_vs_mass_{era}"
    fig.savefig(f"{outname}.png", dpi=400, bbox_inches="tight")
    fig.savefig(f"{outname}.pdf", bbox_inches="tight")
    return fig
