import ROOT
import os
import argparse
from array import array
from plotting.style.ggHcmsstyle import CMSstyle
from ggHparameters import lumi, signal_path, bins
import ggHcuts as cuts
ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT()

all_masses = [15, 20, 30, 40, 50, 55]
all_ctaus = [0, 10, 20, 50, 100, 1000]

ctau_colors = {0: ROOT.kOrange,
               10: ROOT.kRed,
               20: ROOT.kViolet,
               50: ROOT.kAzure,
               100: ROOT.kBlack,
               1000: ROOT.kGreen}


def signal_hist(mass, ctau, year, category, hist_bins):
    path = signal_path(mass, ctau, year)
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    sumw = 0.0
    with ROOT.TFile.Open(path) as f:
        runs_tree = f.Get("Runs")
        for entry in runs_tree:
            sumw += entry.genEventSumw
    if sumw == 0:
        raise ValueError("sum of weights is 0")

    cats = cuts.categories(mass)
    if category not in cats:
        raise ValueError(f"category must be one of {list(cats)}, got '{category}'")

    var = f"best_4g_corr_mass_m{mass}"
    df = ROOT.RDataFrame("ggH4g", path)
    df = df.Define("event_weight", cuts.mc_weight(sumw))
    df = df.Filter(cuts.combine(cuts.trigger(), cuts.dxy_valid(mass),
                                cuts.preselection(mass), cuts.full_id(mass),
                                cuts.pileup()))
    df = df.Filter(cats[category])

    h = df.Histo1D((f"sig_ct{ctau}", f"sig_ct{ctau};m_{{4#gamma}}^{{corr}} [GeV];Events",
                    hist_bins[0], hist_bins[1], hist_bins[2]), var, "event_weight")
    raw = df.Count().GetValue()
    out = h.GetValue().Clone(f"m4g_ct{ctau}")
    out.SetDirectory(0)
    out.Scale(lumi[year])
    return out, raw


def run(mass, year, category, hist_bins, normalize, ctaus):
    hists = {}
    for ct in ctaus:
        try:
            h, raw = signal_hist(mass, ct, year, category, hist_bins)
        except Exception as e:
            print(f"  SKIP ct{ct}: {e}")
            continue
        if h.Integral() <= 0:
            print(f"  SKIP ct{ct}: empty in range")
            continue
        hists[ct] = (h, raw)

    if not hists:
        raise RuntimeError(f"no samples found for m{mass} {year} {category}")

    print(f"m{mass} {year} category={category}  "
          f"range [{hist_bins[1]}, {hist_bins[2]}] GeV")
    probs = array('d', [0.16, 0.50, 0.84])
    for ct, (h, raw) in hists.items():
        quants = array('d', [0.0, 0.0, 0.0])
        h.GetQuantiles(3, quants, probs)
        print(f"  ct={ct:>5} mm  N_raw={raw:>7}  N_wgt={h.Integral():8.3f}  "
              f"mean={h.GetMean():7.3f}  RMS={h.GetRMS():6.3f}  "
              f"median={quants[1]:7.3f}  q16={quants[0]:7.3f}  q84={quants[2]:7.3f}")

    for ct, (h, _) in hists.items():
        if normalize:
            h.Scale(1.0/h.Integral())
        h.SetLineColor(ctau_colors.get(ct, ROOT.kGray+2))
        h.SetLineWidth(3)
        h.SetFillStyle(0)
        h.SetMarkerSize(0)
        h.SetStats(0)
        for axis in (h.GetXaxis(), h.GetYaxis()):
            axis.SetTitleSize(0.05)
            axis.SetLabelSize(0.04)
        h.GetXaxis().SetTitleOffset(1.0)
        h.GetYaxis().SetTitleOffset(1.25)

    c = ROOT.TCanvas("", "", 800, 800)
    right = 0.08
    left = 0.14
    up = 0.08
    down = 0
    l = ROOT.TLegend(0.56+right-left, 0.5+up-down, 0.86+right-left, 0.8+up-down)

    per_bin = round((hist_bins[2]-hist_bins[1])/hist_bins[0], 2)
    ytitle = f"a.u. / {per_bin} GeV" if normalize else f"Events / {per_bin} GeV"
    ymax = 1.6*max(h.GetMaximum() for h, _ in hists.values())

    first = True
    for ct, (h, _) in hists.items():
        h.GetYaxis().SetTitle(ytitle)
        h.SetMaximum(ymax)
        h.SetMinimum(0)
        h.Draw("HIST" if first else "HIST,SAME")
        first = False

    ROOT.gPad.RedrawAxis()
    c.Update()
    l.SetTextSize(0.035)
    CMSstyle(c, l, year, lumi[year], [f"H #rightarrow #phi#phi #rightarrow 4#gamma",
                                      f"m_{{#phi}} = {mass} #font[42]{{GeV}}",
                                      f"category: {category}",
                                      "area normalized" if normalize else "L-scaled yields"])
    for ct, (h, _) in hists.items():
        l.AddEntry(h, f"c#tau = {ct} mm", "l")
    l.Draw("SAME")

    c.Update()
    tag = "shape" if normalize else "yield"
    c.SaveAs(f"m4g_allctau_{tag}_{year}_m{mass}_{category}_corr.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Corrected 4#gamma mass overlaid for all lifetimes")
    parser.add_argument("-m", "--mass", type=str, help="mass of sample")
    parser.add_argument("-y", "--year", type=str, default="2018", help="year of MC")
    parser.add_argument("-c", "--category", type=str, default="none",
                        help="displaced | asym | prompt | none (none = no lxy split)")
    parser.add_argument("-b", "--bins", type=float, nargs=3, default=None,
                        metavar=("NBINS", "LOW", "HIGH"),
                        help=f"histogram binning (default {bins})")
    parser.add_argument("-ct", "--ctaus", type=int, nargs="+", default=all_ctaus,
                        help="lifetimes to overlay")
    parser.add_argument("--absolute", action="store_true",
                        help="draw lumi-scaled yields instead of area-normalized shapes")
    parser.add_argument("-a", "--all", action="store_true",
                        help="loop over all mass points for the given year")
    args = parser.parse_args()

    hist_bins = bins if args.bins is None else [int(args.bins[0]), args.bins[1], args.bins[2]]

    if args.all:
        for m in all_masses:
            print(f"\n=== m{m} {args.year} {args.category} ===")
            try:
                run(str(m), args.year, args.category, hist_bins,
                    not args.absolute, args.ctaus)
            except Exception as e:
                print(f"  SKIP m{m}: {e}")
    else:
        if not args.mass:
            parser.error("-m/--mass is required unless --all is given")
        run(args.mass, args.year, args.category, hist_bins,
            not args.absolute, args.ctaus)
