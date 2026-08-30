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



def sig_bkg_histos(files, isMC, trees, mass, var, output_name, bins, photon_id,histo_names=None, lumi_scaling=1, file_scalings=None):
    if not (len(files) == len(isMC) == len(trees)):
        raise ValueError("files, isMC, and trees must have the same length")
    if file_scalings is None:
        file_scalings = [lumi_scaling] * len(files)
    if len(file_scalings) != len(files):
        raise ValueError("file_scalings must have one value per input file")

    sumw_dict={}
    file_num = len(files)
    if histo_names is None:
        histo_names = [f"hist_{i}" for i in range(file_num)]
    if len(histo_names) != file_num:
        raise ValueError("histo_names must have one value per input file")
    histo_obj = []

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
    output_file = ROOT.TFile(output_name, "RECREATE")
    for i in range(file_num):
        dataframe = ROOT.RDataFrame(trees[i], files[i])

        if isMC[i]:
            weight_formula = cuts.mc_weight(sumw_dict[files[i]])
            dataframe = dataframe.Filter(cuts.signal_selection(mass, photon_id))
            dataframe=dataframe.Define("event_weight", weight_formula)
            histogram = dataframe.Histo1D((histo_names[i], f"{i};{var};Events", bins[0], bins[1], bins[2]),var,"event_weight")
        else:
            dataframe = dataframe.Filter(cuts.background_selection(mass, photon_id))
            histogram = dataframe.Histo1D((histo_names[i], f"{i};{var};Events", fit_bins[0], fit_bins[1], fit_bins[2]),var)

        histogram.Scale(file_scalings[i])
        histogram.Write()
        histo_obj.append(histogram)
    output_file.Close()
    return output_name, histo_names, histo_obj

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
    
