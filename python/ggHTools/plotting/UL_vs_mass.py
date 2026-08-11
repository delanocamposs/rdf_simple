import ROOT
import subprocess
import json
import os
import sys
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as mh
from datacard.ggHdatacardmaker import main as make_datacard
from ggHparameters import signal_path, bkg_path, lumi


# Run3 only has per-subera signal MC (no aggregated 2022/2023 signal files), so
# the era->years map uses the suberas that actually have datacards.
ERA_YEARS = {
    "Run2": ["2017", "2018"],
    "Run3": ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"],
}
ERA_YEARS["Run2Run3"] = ERA_YEARS["Run2"] + ERA_YEARS["Run3"]
_ERA_ALIASES = {"run2": "Run2", "run3": "Run3",
                "run2run3": "Run2Run3", "run3run2": "Run2Run3",
                "run2+run3": "Run2Run3"}


def resolve_era(token):
    """Map a CLI token to (era_label, [years]).

    Accepts a single year tag (e.g. "2018", "2022preEE") or an era keyword
    (Run2, Run3, Run2Run3 / Run2+Run3, case-insensitive).
    """
    key = _ERA_ALIASES.get(token.strip().lower())
    if key:
        return key, list(ERA_YEARS[key])
    return token, [token]


def era_lumi(years):
    return sum(lumi[y] for y in years)


# AsymptoticLimits stores one tree entry per quantile. Map the quantileExpected
# value combine writes to the keys we use everywhere downstream.
_QUANTILE_KEYS = {0.025: "exp-2", 0.16: "exp-1",
                  0.5: "exp0", 0.84: "exp+1", 0.975: "exp+2"}


def get_limits(root_file):
    """Return all 95% CL upper limits on r from an AsymptoticLimits root file.

    Returns a dict like {"exp-2": .., "exp-1": .., "exp0": ..,
    "exp+1": .., "exp+2": ..}. Missing quantiles are simply absent.
    """
    f = ROOT.TFile.Open(root_file)
    tree = f.Get("limit") if f else None
    out = {}
    if not tree:
        print(f"[get_limits] WARNING: no 'limit' tree in {root_file}; skipping point.")
        if f:
            f.Close()
        return out
    for entry in tree:
        for q, key in _QUANTILE_KEYS.items():
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


def _combine_cards(labelled_cards, out_txt):
    """combineCards.py wrapper. labelled_cards is a list of (label, txt) pairs;
    explicit labels keep channel/nuisance names unique across years."""
    args = [f"{label}={card}" for label, card in labelled_cards]
    with open(out_txt, "w") as f:
        subprocess.run(["combineCards.py", *args], stdout=f, check=True)
    return out_txt


def scan_mass_lifetime(masses, lifetimes, era, years, categories, bins, finalstate="4g",
                       physics="ggH", order_fit=4, results_json=None,
                       photon_id="custom"):
    """Build per-(year, category) datacards, statistically combine them across the
    whole era, and run AsymptoticLimits over the mass/lifetime grid.

    For a single-year era this combines only the categories; for Run2/Run3/Run2Run3
    it combines every (year, category) channel into one card per (mass, ctau).

    Returns (and saves to JSON) a nested dict:
        results[era][str(ctau)][str(mass)] = {"exp0":.., "exp-1":.., ...}   # raw r
    Raw r limits are stored; BR scaling is applied at plot time via br_scale.
    """
    ROOT.gROOT.SetBatch(True)
    if results_json is None:
        results_json = f"limits_UL_vs_mass_{era}.json"
    cats = list(categories)
    if "none" in cats and len(cats) != 1:
        raise ValueError("the inclusive 'none' category cannot be combined with exclusive categories")
    results = {era: {}}
    for mass in masses:
        for ctau in lifetimes:
            labelled_cards = []
            for year in years:
                sig = signal_path(mass, ctau, year)
                bkg = bkg_path(year)
                for cat in cats:
                    make_datacard(paths=[sig, bkg], isMC=[1, 0], trees=["ggH4g", "ggH4g"],
                                  var=f"best_4g_corr_mass_m{mass}", categories=[cat], period=year,
                                  bins=bins, lifetime=ctau, mass=mass, finalstate=finalstate,
                                  physics=physics, order_fit=order_fit, photon_id=photon_id)
                    card = f"datacard_{physics}_{finalstate}_m{mass}_ct{ctau}_{cat}_{year}.txt"
                    labelled_cards.append((f"{cat}_{year}", card))

            combined_txt = f"datacard_{physics}_{finalstate}_m{mass}_ct{ctau}_{era}.txt"
            combined_root = f"datacard_{physics}_{finalstate}_m{mass}_ct{ctau}_{era}.root"
            _combine_cards(labelled_cards, combined_txt)
            subprocess.run(["text2workspace.py", combined_txt, "-o", combined_root], check=True)
            subprocess.run(["combine", "-M", "AsymptoticLimits", combined_root, "-m", "125", "--run", "expected"], check=True)

            result_root = f"higgsCombineTest_m{mass}_ct{ctau}_{era}.AsymptoticLimits.mH125.root"
            subprocess.run(["mv", "higgsCombineTest.AsymptoticLimits.mH125.root", result_root], check=True)

            lims = get_limits(result_root)
            if lims:
                results[era].setdefault(str(ctau), {})[str(mass)] = lims
            print(f"m={mass} ct={ctau} era={era}: expected UL on r = {lims.get('exp0')}")

    save_results(results, results_json)
    return results


def _panel_arrays(panel):
    masses = sorted(int(m) for m in panel)
    needed = ["exp-2", "exp-1", "exp0", "exp+1", "exp+2"]
    masses = [m for m in masses if all(k in panel[str(m)] for k in needed)]
    arrays = {k: [panel[str(m)][k] for m in masses] for k in needed}
    return masses, arrays


def plot_UL_vs_mass(results, era, ctaus, total_lumi=None, br_scale=1.0,
                    ytitle=r"95% upper limit on $\mathcal{B}(H\rightarrow\Phi\Phi)$",
                    extra_labels=(r"$\mathcal{B}(\Phi\rightarrow\gamma\gamma)=1$",),
                    yrange=None, outname=None):
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
        masses, arrays = _panel_arrays(year_res[str(ct)])
        arrays = {k: [v * br_scale for v in vals] for k, vals in arrays.items()}
        panels.append((ct, masses, arrays))

    if not any(masses for _, masses, _ in panels):
        print("[plot_UL_vs_mass] no complete points to plot.")
        return

    band2sigma_color = "#85d2fb"
    band1sigma_color = "#ffde9c"
    fig, axes = plt.subplots(1, n, sharey=True, squeeze=False,
                             figsize=(7.5 * n, 20))
    axes = axes[0]
    plt.subplots_adjust(wspace=0)

    for i, (ct, masses, arrays) in enumerate(panels):
        ax = axes[i]
        if masses:
            ax.fill_between(masses, arrays["exp-2"], arrays["exp+2"],
                            color=band2sigma_color, label=r"$\pm 2\sigma$")
            ax.fill_between(masses, arrays["exp-1"], arrays["exp+1"],
                            color=band1sigma_color, label=r"$\pm 1\sigma$")
            ax.plot(masses, arrays["exp0"], linestyle="--", color="black",
                    linewidth=4, label="Expected")

        ax.text(0.5, 0.02, rf"$c\tau = {ct}\ \mathrm{{mm}}$",
                transform=ax.transAxes, fontsize=40, ha="center", va="center",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.5))
        ax.tick_params(axis="both", which="major", labelsize=40)
        ax.margins(x=0)
        ax.minorticks_on()

    axes[0].legend(loc="upper right", fontsize=40)
    axes[-1].set_xlabel(r"$m_\Phi$ (GeV)", fontsize=40)
    axes[0].set_ylabel(ytitle, fontsize=40)
    axes[0].set_yscale("log")
    if yrange is not None:
        axes[0].set_ylim(*yrange)

    mh.cms.label("Preliminary", data=True, rlabel="", ax=axes[0], loc=0,
                 fontsize=60)
    mh.cms.label(None, exp="", data=True, llabel="", ax=axes[-1], loc=0,
                 lumi=f"{total_lumi / 1000:.2f}", com=13, fontsize=40)
    for j, label in enumerate(extra_labels):
        axes[0].text(0.5, 0.70 - 0.06 * j, label,
                     transform=axes[0].transAxes, fontsize=40,
                     ha="center", va="center",
                     bbox=dict(facecolor="white", edgecolor="none", alpha=0.5))

    if outname is None:
        outname = f"UL_vs_mass_{era}"
    fig.savefig(f"{outname}.png", dpi=400, bbox_inches="tight")
    fig.savefig(f"{outname}.pdf", bbox_inches="tight")
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "UL vs mass limits. Scans all mass/ctau/category points and builds the "
        "combined-dataset limit for a single year or a whole era.")
    parser.add_argument("-y", "--year", type=str,
                        help="single year (e.g. 2018, 2022preEE) or era: Run2, Run3, Run2Run3")
    parser.add_argument("-process_run2", "--process_run2", dest="process_run2", type=int,
                        help="1=run full Run2 (2017+2018) combined")
    parser.add_argument("-process_run3", "--process_run3", dest="process_run3", type=int,
                        help="1=run full Run3 (2022/2023/2024 suberas) combined")
    parser.add_argument("-process_run2run3", "--process_run2run3", dest="process_run2run3", type=int,
                        help="1=run combined Run2+Run3")
    parser.add_argument("--rescan", action="store_true",
                        help="force re-running combine even if a cached json exists")
    parser.add_argument("--inclusive", action="store_true",
                        help="use one inclusive 'none' category instead of prompt/asym/displaced")
    parser.add_argument("--photon-id", choices=["custom", "LooseEGM", "MediumEGM", "TightEGM"],
                        default="custom", help="photon ID applied to signal and sideband data")
    args = parser.parse_args()

    if args.process_run2:
        era_token = "Run2"
    elif args.process_run3:
        era_token = "Run3"
    elif args.process_run2run3:
        era_token = "Run2Run3"
    elif args.year:
        era_token = args.year
    else:
        print("\033[91mERROR: specify -y/--year <year|Run2|Run3|Run2Run3>, or "
              "-process_run2/-process_run3/-process_run2run3 1.\033[0m")
        sys.exit(1)

    era, years = resolve_era(era_token)

    masses = [15, 20, 30, 40, 50, 55]
    lifetimes = [0, 10, 20, 50, 100, 1000]
    categories = ["none"] if args.inclusive else ["prompt", "asym", "displaced"]
    scan_label = f"{era}_inclusive" if args.inclusive else era
    if args.photon_id != "custom":
        scan_label = f"{scan_label}_{args.photon_id}"
    results_json = f"limits_UL_vs_mass_{scan_label}.json"

    # Re-use existing limits if present; otherwise run the (slow) combine scan.
    if os.path.exists(results_json) and not args.rescan:
        print(f"loading cached limits from {results_json} (pass --rescan to re-run)")
        results = load_results(results_json)
    else:
        results = scan_mass_lifetime(masses, lifetimes, scan_label, years, categories,
                                     bins=[30, 110, 140], results_json=results_json,
                                     photon_id=args.photon_id)

    # br_scale: r->BR conversion applied to every y value (= reference BR=1e-4).
    id_label = "Custom photon ID" if args.photon_id == "custom" else args.photon_id.replace("EGM", " EGM photon ID")
    plot_UL_vs_mass(results, scan_label, lifetimes, total_lumi=era_lumi(years),
                    br_scale=1e-4,
                    extra_labels=(r"$\mathcal{B}(\Phi\rightarrow\gamma\gamma)=1$", id_label))
