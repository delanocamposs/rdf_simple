import ROOT 
import os 
import sys 
import argparse, subprocess
from plotting.style.ggHcmsstyle import CMSstyle
from ggHparameters import lumi, signal_window, signal_path, bkg_path
import ggHcuts as cuts
ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT()


def run(mass, ctau, year):
    c=ROOT.TCanvas("", "", 800, 800)
    right=0.08
    left=0.14
    up=0.08
    down=0
    l=ROOT.TLegend(0.5+right-left,0.67+up-down, 0.95+right-left, 0.8+up-down)
    bins=[55, 60, 180]
    var=f"best_4g_corr_mass_m{mass}"
    
    sig=signal_path(mass, ctau, year)
    bkg=bkg_path(year)
    
    sig_open=ROOT.TFile.Open(sig)
    sumw=0.0
    signal_df=ROOT.RDataFrame("ggH4g", sig)
    with sig_open as f:
        runs_tree = f.Get("Runs")
        for entry in runs_tree:
            sumw += entry.genEventSumw
    weight_formula = cuts.mc_weight(sumw)
    signal_df=signal_df.Define("event_weight", weight_formula)
    
    data_ID_df=ROOT.RDataFrame("ggH4g", bkg)
    data_pre_df=ROOT.RDataFrame("ggH4g", bkg)
    
    signal_df=signal_df.Filter(cuts.combine(cuts.trigger(), cuts.dxy_valid(mass), cuts.full_id(mass), cuts.pileup()))
    data_ID_df=data_ID_df.Filter(cuts.combine(cuts.preselection(mass), cuts.trigger(), cuts.dxy_valid(mass), cuts.full_id(mass), cuts.blind(mass)))
    data_pre_df=data_pre_df.Filter(cuts.combine(cuts.trigger(), cuts.dxy_valid(mass), cuts.preselection(mass), cuts.blind(mass)))
    
    signal_histo=signal_df.Histo1D(("hist1_1", f"hist1_1;{var};Events", bins[0], bins[1], bins[2]), f"best_4g_corr_mass_m{mass}", "event_weight")
    signal_histo.Scale(lumi[year])
    
    data_ID_histo=data_ID_df.Histo1D(("hist2", f"hist2;4#gamma mass;Events", bins[0], bins[1], bins[2]), f"{var}")
    data_pre_histo=data_pre_df.Histo1D(("hist3", f"hist3;4#gamma mass;Events", bins[0], bins[1], bins[2]), f"{var}")
    lower_sb=f"{var}>={bins[1]} && {var}<{signal_window[0]}"
    N_data_ID=data_ID_df.Filter(lower_sb).Count().GetValue()
    N_data_pre=data_pre_df.Filter(lower_sb).Count().GetValue()
    N_sig=signal_histo.Integral()
    print(f"lower sideband [{bins[1]}, {signal_window[0]}) GeV")
    print("N_ID", N_data_ID)
    print("N_pre", N_data_pre)
    print("N_sig", N_sig)
    
    data_ID_color=ROOT.kOrange+8
    data_pre_color=ROOT.kPink+1
    signal_color=ROOT.kAzure-4
    data_ID_histo.SetFillStyle(1001)
    data_ID_histo.SetLineColor(ROOT.kBlack)
    data_ID_histo.SetLineStyle(1)
    data_ID_histo.SetFillColor(data_ID_color)
    data_pre_histo.SetFillStyle(1001)
    data_pre_histo.SetLineColor(ROOT.kBlack)
    data_pre_histo.SetLineStyle(1)
    data_pre_histo.SetFillColor(data_pre_color)
    signal_histo.SetFillStyle(1001)
    signal_histo.SetFillColor(signal_color)
    signal_histo.SetLineColor(ROOT.kBlack)
    signal_histo.SetLineStyle(1)
    per_bin=round((bins[2]-bins[1])/bins[0], 2)
    
    data_ID_clone=data_ID_histo.Clone("data_ID_clone")
    data_pre_clone=data_pre_histo.Clone("data_pre_clone")
    signal_histo_clone=signal_histo.Clone("signal_histo_clone")
    signal_histo_clone.SetFillColor(signal_color)
    data_ID_clone.SetFillColor(data_ID_color)
    data_pre_clone.SetFillColor(data_pre_color)
    data_ID_clone.SetFillStyle(1001)
    data_pre_clone.SetFillStyle(1001)
    signal_histo_clone.SetFillStyle(1001)
    data_ID_clone.SetLineColor(ROOT.kBlack)
    data_pre_clone.SetLineColor(ROOT.kBlack)
    signal_histo_clone.SetLineColor(ROOT.kBlack)
    
    data_ID_histo.SetLineColor(ROOT.kBlack)
    data_pre_histo.SetLineColor(ROOT.kBlack)
    signal_histo.SetLineColor(ROOT.kBlack)
    data_ID_histo.SetLineWidth(2)
    data_pre_histo.SetLineWidth(2)
    signal_histo.SetLineWidth(2)
    ymax=100*max(data_ID_histo.GetMaximum(), data_pre_histo.GetMaximum(), signal_histo.GetMaximum())
    data_pre_histo.SetMaximum(ymax)
    data_pre_histo.SetMinimum(0.01)
    
    c.SetLogy()
    data_pre_histo.Draw("HIST")
    data_ID_histo.Draw("HIST,SAME")
    signal_histo.Draw("HIST,SAME")
    
    ROOT.gPad.RedrawAxis()  
    ROOT.gPad.RedrawAxis("g") 
    
    c.Update()
    l.SetTextSize(0.05)
    CMSstyle(c, l,year, lumi[year],[f"H #rightarrow #phi#phi #rightarrow 4#gamma", f"c#tau = {ctau} mm", f"m_{{#phi}} = {mass} #font[42]{{GeV}}"])
    l.AddEntry(signal_histo_clone, "Signal MC", "f")
    l.AddEntry(data_ID_clone, "ID Data (blinded)", "f")
    l.AddEntry(data_pre_clone, "pre Data (blinded)", "f")
    l.Draw("SAME")
    
    c.Update()
    c.SaveAs(f"ID_pre_sideband_{year}_ct{ctau}_m{mass}.png")
    SF=N_data_ID/N_data_pre
    print(f"SF={SF}") 


if __name__=="__main__":
    parser = argparse.ArgumentParser("Running sideband plot")
    parser.add_argument("-m","--mass", type=str, help="mass of sample")
    parser.add_argument("-ct","--ctau", type=str, help="lifetime of sample")
    parser.add_argument("-y","--year", type=str, help="year of MC and data")
    args = parser.parse_args()
    mass=args.mass
    lifetime=args.ctau
    year=args.year
    run(mass, lifetime, year)
    
