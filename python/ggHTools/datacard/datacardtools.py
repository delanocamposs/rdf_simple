import ROOT 
import uproot
import awkward as ak
from ROOT import RooFit, RooArgList
import sys
import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.special import comb
from scipy.integrate import simps
from scipy.stats import beta
from datacard.ggHfitter import fitBKG
from ggHparameters import fit_bins
import ggHcuts as cuts
ROOT.gROOT.SetBatch(False)



def sig_bkg_histos(files, isMC, trees, mass, lifetime, selections, var, output_names, bins, year, histo_names=[], bkg_weight=True, lumi_scaling=1):
    '''
    files = path to the root data/MC files to construct the histos from
    isMC = is the file MC or not (matters because it gets weighted and scaled correctly if it is)
    trees = the trees inside each file that contain the event data 
    selections = cuts to be applied to the files provided 
    var = the independetn variable in the histograms 
    output_names = the names of the new files which hold the new histos 
    histo_names = new names of the new histos made inside the new root output file 
    bkg_weight = the weight we scale the background histogram by (derived SF from looking at blinded data compared to background)
    bins = the binning of the new histograms

    EXAMPLE:

    sig_bkg_histos(["path/to/MC", "path/to/bkg"], [1, 0], ["ggH4g", "ggH4g"], ["best_4g_phi1_dxy_m{mass}>50", "best_4g_phi2_dxy_m{mass}<50"], "best_4g_cor_mass_m{mass}", ["output1.root", "output2.root"], [["hist1", "hist2"], ["hist3", "hist4"]])
    '''

    outputs=[]
    sumw_dict={}
    if len(output_names) != len(selections):
        print(f"ERROR: the number of selections must match the number of desired saved output file names. one saved file for each selection")
        return
    selection_num = len(selections)
    file_num = len(files)
    if len(histo_names) == 0:
        histo_names = [[] for _ in range(selection_num)]
        for i in range(selection_num):
            for j in range(file_num):
                histo_names[i].append(f"hist_{i}_{j}")
    histo_obj = [[] for _ in range(selection_num)]

    for i in range(len(files)):
        filepath_i = files[i]
        if isMC[i]:
            sumw=0.0
            with ROOT.TFile.Open(filepath_i) as f:
                runs_tree = f.Get("Runs")
                if runs_tree:
                    for entry in runs_tree:
                        sumw += entry.genEventSumw
            if sumw == 0:
                print("sum of weights is 0")
                sumw = 1.0
            sumw_dict[filepath_i]=sumw
    for j in range(len(selections)):
        selection_j = selections[j]
        output_file_j = ROOT.TFile(f"{output_names[j]}", "RECREATE")
        outputs.append(output_file_j)
        for k in range(len(files)):
            rdf_j_k=ROOT.RDataFrame(trees[k], files[k])
            rdf_j_k = rdf_j_k.Filter(cuts.combine(cuts.trigger(), cuts.dxy_valid(mass))).Filter(selection_j)
            if isMC[k]:
                weight_formula_k = cuts.mc_weight(sumw_dict[files[k]])
                rdf_j_k = rdf_j_k.Filter(cuts.combine(cuts.preselection(mass), cuts.full_id(mass), cuts.pileup())).Define("event_weight", weight_formula_k)
                hist_j_k = rdf_j_k.Histo1D((f"{histo_names[j][k]}", f"{j}_{k};{var};Events", bins[0], bins[1], bins[2]), f"{var}", "event_weight")
                hist_j_k.Scale(lumi_scaling)
                hist_j_k.Write()
                histo_obj[j].append(hist_j_k)
            else:
                rdf_j_k = rdf_j_k.Filter(cuts.combine(cuts.preselection(mass), cuts.full_id(mass), cuts.sidebands(mass)))
                hist_j_k = rdf_j_k.Histo1D((f"{histo_names[j][k]}", f"{j}_{k};{var};Events", fit_bins[0], fit_bins[1], fit_bins[2]), f"{var}")
                hist_j_k.Scale(lumi_scaling)
                hist_j_k.Write()
                histo_obj[j].append(hist_j_k)
        output_file_j.Close()
    return outputs, output_names, histo_names, histo_obj

def extract_JSON(root_filename, workspace_name, json_filename):
    '''
    creates a JSON file in the format the datcard needs from a root file, the name of the 
    RooWorkspace in the file.
    '''
    f = ROOT.TFile(root_filename)
    ws = f.Get(workspace_name)
    if not ws:
        print(f"Error: Workspace '{workspace_name}' not found in {root_filename}")
        return
    params = {}
    all_vars = ws.allVars()
    it = all_vars.createIterator()
    var = it.Next()
    while var:
        if var.InheritsFrom("RooRealVar"):
            params[var.GetName()] = {"value": var.getVal(), "error": var.getError()}
        var = it.Next()
    with open(json_filename, "w") as json_file:
        json.dump(params, json_file, indent=4)
    f.Close()


def generate_data_hist(file, bins_num, norm, output_name):
    f=ROOT.TFile.Open(file)
    w=f.Get("w")
    pdf =w.pdf("model")
    x=w.var("mass")
    h_pdf =pdf.createHistogram("h_pdf", x, ROOT.RooFit.Binning(bins_num))
    h_pdf1 =pdf.createHistogram("h_pdf1", x, ROOT.RooFit.Binning(bins_num))
    h_pdf.Scale(norm)
    h_pdf1.Scale(norm)
    total_expected = float(h_pdf.Integral())
    if total_expected < 1:
        print("[generate_data_hist] WARNING: expected total events is < 1. A fully empty toy histogram is likely and can be statistically consistent.")
    output_file = ROOT.TFile(output_name, "RECREATE")
    h_pdf1.Write()
    for i in range(1, h_pdf.GetNbinsX()+1):
        density=h_pdf.GetBinContent(i)
        bw=h_pdf.GetXaxis().GetBinWidth(i)
        h_pdf.SetBinContent(i, np.random.poisson(density))
    h_pdf.Write()
    output_file.Close()
    return output_name


def clopper_pearson(X, n, alpha=0.05):
    "by default i am doing 95% CL, change alpha to change this: CL=1-alpha"
    lower=beta.ppf(alpha/2, X,n-X+1)
    upper=beta.ppf(1-(alpha/2),X+1,n-X)
    return lower,upper
    

