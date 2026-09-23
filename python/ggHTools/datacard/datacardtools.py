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

#### function that just builds signal and backgroudn histgorams with all cuts applied and returns the object

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

def preselected_counts(path, tree, mass, photon_id):

    dataframe=ROOT.RDataFrame(tree, path)
    base=dataframe.Filter(cuts.preselected(mass)).Filter(cuts.fails_photon_id(mass, photon_id))
    sideband=base.Filter(cuts.sidebands(mass)).Count()
    signal_region=base.Filter(cuts.signal_region(mass)).Count()
    return int(sideband.GetValue()), int(signal_region.GetValue())


def extract_JSON(root_filename, workspace_name, json_filename):

### makes a json file by extracting initial parameters from rooworkspace

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


def expected_hist(file, name, bins_num, norm, pdf_name="model", var_name="mass"):

    f=ROOT.TFile.Open(file)
    w=f.Get("w")
    if not w:
        raise ValueError("no workspace 'w' in {}".format(file))
    pdf=w.pdf(pdf_name)
    x=w.var(var_name)
    if not pdf or not x:
        raise ValueError("workspace in {} is missing pdf '{}' or variable '{}'".format(file, pdf_name, var_name))
    histogram=pdf.createHistogram(name, x, ROOT.RooFit.Binning(bins_num))
    histogram.SetDirectory(0)
    integral=float(histogram.Integral())
    if integral<=0.0:
        raise ValueError("pdf '{}' in {} integrates to {} over the histogram range".format(pdf_name, file, integral))
    histogram.Scale(norm/integral)
    f.Close()
    return histogram


def generate_data_hist(file, bins_num, norm, output_name, signal_file=None, signal_norm=0.0, seed=None):

    expected=expected_hist(file, "h_pdf1", bins_num, norm)
    if signal_file is not None and signal_norm>0.0:
        signal_expected=expected_hist(signal_file, "h_sig", bins_num, signal_norm)
        expected.Add(signal_expected)
    toy=expected.Clone("h_pdf")
    toy.SetDirectory(0)
    total_expected=float(expected.Integral())
    if total_expected<1:
        print("[generate_data_hist] WARNING: expected total events is < 1. A fully empty toy histogram is likely and can be statistically consistent.")
    rng=np.random.default_rng(seed)
    for i in range(1, toy.GetNbinsX()+1):
        toy.SetBinContent(i, float(rng.poisson(expected.GetBinContent(i))))
        toy.SetBinError(i, 0.0)
    output_file=ROOT.TFile(output_name, "RECREATE")
    output_file.cd()
    expected.Write("h_pdf1__mass")
    toy.Write("h_pdf__mass")
    output_file.Close()
    return output_name, total_expected, float(toy.Integral())


def fitted_poi(fitdiagnostics_file, poi="r", fit="fit_s"):

    f=ROOT.TFile.Open(fitdiagnostics_file)
    if not f or f.IsZombie():
        return None
    result=f.Get(fit)
    if not result:
        f.Close()
        return None
    parameter=result.floatParsFinal().find(poi)
    if not parameter:
        f.Close()
        return None
    values=(parameter.getVal(), parameter.getErrorLo(), parameter.getErrorHi(), result.status(), result.covQual())
    f.Close()
    return values


def clopper_pearson(X, n, alpha=0.05):

# clopper pearson interval, default at 95% (CL=1-alpha)

    lower=beta.ppf(alpha/2, X,n-X+1)
    upper=beta.ppf(1-(alpha/2),X+1,n-X)
    return lower,upper
    
