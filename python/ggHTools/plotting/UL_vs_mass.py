import ROOT
import subprocess
import json
import os
import sys
import argparse
from datacard.ggHdatacardmaker import main as make_datacard
from ggHparameters import signal_path, bkg_path, lumi
from plotting.style import tdrstyle


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
                       physics="ggH", order_fit=4, results_json=None):
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
    cats = [c for c in categories if c != "none"]
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
                                  physics=physics, order_fit=order_fit)
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


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _band_graph(masses, lo, hi, color):
    """Closed-polygon TGraph for a Brazil band between lo[] and hi[] vs masses[]."""
    n = len(masses)
    g = ROOT.TGraph(2 * n)
    for i in range(n):
        g.SetPoint(i, masses[i], hi[i])
    for i in range(n):
        g.SetPoint(n + i, masses[n - 1 - i], lo[n - 1 - i])
    g.SetFillColor(color)
    g.SetLineColor(color)
    return g


def _line_graph(masses, vals, style, width=2):
    g = ROOT.TGraph(len(masses))
    for i, (x, y) in enumerate(zip(masses, vals)):
        g.SetPoint(i, x, y)
    g.SetLineColor(ROOT.kBlack)
    g.SetLineStyle(style)
    g.SetLineWidth(width)
    g.SetMarkerSize(0)
    return g


def _panel_arrays(panel):
    masses = sorted(int(m) for m in panel)
    needed = ["exp-2", "exp-1", "exp0", "exp+1", "exp+2"]
    masses = [m for m in masses if all(k in panel[str(m)] for k in needed)]
    arrays = {k: [panel[str(m)][k] for m in masses] for k in needed}
    return masses, arrays


def plot_UL_vs_mass(results, era, ctaus, total_lumi=None, br_scale=1.0,
                    ytitle="95% upper limit on BR(H#rightarrow#Phi#Phi)",
                    extra_labels=("BR(#Phi#rightarrow#gamma#gamma) = 1",),
                    yrange=None, outname=None):
    ROOT.gROOT.SetBatch(True)
    tdrstyle.setTDRStyle()
    ROOT.gStyle.SetOptStat(0)
    # lighter, less distracting grid lines
    ROOT.gStyle.SetGridColor(ROOT.kGray)
    ROOT.gStyle.SetGridStyle(3)
    ROOT.gStyle.SetGridWidth(1)

    if total_lumi is None:
        total_lumi = lumi.get(era)
    year_res = results[era] if era in results else results
    ctaus = [ct for ct in ctaus if str(ct) in year_res]
    n = len(ctaus)
    if n == 0:
        print("[plot_UL_vs_mass] no lifetimes with data to plot.")
        return

    # gather scaled arrays per panel + global ranges
    panels, all_y, all_m = [], [], []
    for ct in ctaus:
        masses, arrays = _panel_arrays(year_res[str(ct)])
        arrays = {k: [v * br_scale for v in vals] for k, vals in arrays.items()}
        panels.append((ct, masses, arrays))
        for vals in arrays.values():
            all_y += vals
        all_m += masses

    if not all_y:
        print("[plot_UL_vs_mass] no complete points to plot.")
        return

    ymin,ymax =(min(all_y) * 0.5, max(all_y) * 2.0) if yrange is None else yrange
    xlo,xhi= min(all_m), max(all_m)
    lm,rm,tm,bm = 0.23, 0.035, 0.07, 0.14
    denom = (1.0 / (1 - lm)) + max(n - 2, 0) + (1.0 / (1 - rm) if n > 1 else 0)
    if n == 1:
        denom = 1.0 / (1 - lm - rm)
    fw = 1.0 / denom
    widths = []
    for i in range(n):
        if n == 1:
            widths.append(fw / (1 - lm - rm))
        elif i == 0:
            widths.append(fw / (1 - lm))
        elif i == n - 1:
            widths.append(fw / (1 - rm))
        else:
            widths.append(fw)
    edges = [0.0]
    for w in widths:
        edges.append(edges[-1] + w)
    edges = [e / edges[-1] for e in edges]

    canv = ROOT.TCanvas("UL_vs_mass", "UL_vs_mass", 260 * max(n, 4), 600)
    keep = []
    green, yellow = ROOT.kGreen + 1, ROOT.kOrange

    for i, (ct, masses, arrays) in enumerate(panels):
        canv.cd()
        pad = ROOT.TPad(f"pad{i}", "", edges[i], 0.0, edges[i + 1], 1.0)
        pad.SetLeftMargin(lm if (i == 0 or n == 1) else 0.0)
        pad.SetRightMargin(rm if (i == n - 1 or n == 1) else 0.0)
        pad.SetTopMargin(tm)
        pad.SetBottomMargin(bm)
        pad.SetLogy()
        pad.SetTicks(1, 1)
        pad.SetGridx()
        pad.SetGridy()
        pad.Draw()
        keep.append(pad)
        pad.cd()

        frame = pad.DrawFrame(xlo, ymin, xhi, ymax)
        frame.GetXaxis().SetNdivisions(505)
        frame.GetXaxis().SetTitleSize(0.06)
        frame.GetXaxis().SetLabelSize(0.055)
        frame.GetXaxis().SetLabelOffset(-0.005)
        frame.GetYaxis().SetTitle(ytitle if i == 0 else "")
        frame.GetYaxis().SetTitleSize(0.062)
        frame.GetYaxis().SetTitleOffset(1.8)
        if i != 0 and n > 1:
            frame.GetYaxis().SetLabelSize(0.0)
        else:
            frame.GetYaxis().SetLabelSize(0.045)
        keep.append(frame)

        if masses:
            g2 = _band_graph(masses, arrays["exp-2"], arrays["exp+2"], yellow)
            g1 = _band_graph(masses, arrays["exp-1"], arrays["exp+1"], green)
            g2.Draw("F same")
            g1.Draw("F same")
            gexp = _line_graph(masses, arrays["exp0"], 2)  # dashed expected
            gexp.Draw("L same")
            keep += [g2, g1, gexp]

        pad.RedrawAxis()
        pad.RedrawAxis("g")  # grid lines on top of the bands

        cx = pad.GetLeftMargin() + 0.5 * (1 - pad.GetLeftMargin() - pad.GetRightMargin())
        tl = ROOT.TLatex()
        tl.SetNDC()
        tl.SetTextFont(42)
        tl.SetTextAlign(22)
        tl.SetTextSize(0.07)
        tl.DrawLatex(cx, bm + 0.04, f"c#tau = {ct} mm")
        keep.append(tl)

        # legend + extra labels only in the first panel
        if i == 0 and masses:
            leg = ROOT.TLegend(lm + 0.04, 0.65, lm + 1.1, 0.90)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            leg.SetTextSize(0.06)
            leg.AddEntry(g1, "#pm1#sigma", "f")
            leg.AddEntry(g2, "#pm2#sigma", "f")
            leg.AddEntry(gexp, "Expected", "l")
            leg.Draw()
            keep.append(leg)

            lab = ROOT.TLatex()
            lab.SetNDC()
            lab.SetTextFont(42)
            lab.SetTextSize(0.05)
            for j, txt in enumerate(extra_labels):
                lab.DrawLatex(lm + 0.04, 0.60 - 0.06 * j, txt)
            keep.append(lab)

    # header drawn on the canvas, spanning full width
    canv.cd()
    head = ROOT.TLatex()
    head.SetNDC()
    head.SetTextAlign(13)
    cms_size = 0.05
    head.SetTextFont(61)
    head.SetTextSize(cms_size)
    head.DrawLatex(0.04, 0.985, "CMS")
    # place "Preliminary" right after "CMS": its NDC width scales with the (wide)
    # canvas aspect ratio, so a fixed NDC offset would leave a large pixel gap.
    cms_w = cms_size * (canv.GetWh() / canv.GetWw()) * 1.8
    head.SetTextFont(52)
    head.SetTextSize(0.038)
    head.DrawLatex(0.04 + cms_w + 0.006, 0.983, "Preliminary")
    head.SetTextFont(42)
    head.SetTextAlign(33)
    head.SetTextSize(0.038)
    head.DrawLatex(1 - rm, 0.985, f"{total_lumi / 1000:.2f} fb^{{-1}} (13 TeV)")
    # global x-axis title near the bottom-right
    head.SetTextAlign(31)
    head.SetTextSize(0.045)
    head.DrawLatex(1 - rm, 0.035, "m_{#Phi} (GeV)")
    keep.append(head)

    canv.Update()
    if outname is None:
        outname = f"UL_vs_mass_{era}"
    canv.SaveAs(f"{outname}.png")
    canv.SaveAs(f"{outname}.pdf")
    return canv


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
    categories = ["prompt", "asym", "displaced"]
    results_json = f"limits_UL_vs_mass_{era}.json"

    # Re-use existing limits if present; otherwise run the (slow) combine scan.
    if os.path.exists(results_json) and not args.rescan:
        print(f"loading cached limits from {results_json} (pass --rescan to re-run)")
        results = load_results(results_json)
    else:
        results = scan_mass_lifetime(masses, lifetimes, era, years, categories,
                                     bins=[30, 110, 140], results_json=results_json)

    # br_scale: r->BR conversion applied to every y value (= reference BR=1e-4).
    plot_UL_vs_mass(results, era, lifetimes, total_lumi=era_lumi(years), br_scale=1e-4)
