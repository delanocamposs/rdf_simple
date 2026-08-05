from datacard import datacardtools
from datacard import ggHfitter
from datacard.ggHdatacardworkspace import DatacardWorkspace
from ggHparameters import order_fit, order_gen, lxy1, lxy2, smear_resolution, lumi, xsec_unc, pdf_alphas_unc, lumi_unc
import json 
import os 
import ROOT
import argparse, subprocess
import numpy as np
import shutil
import glob

ROOT.gROOT.SetBatch(False)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
ROOT.gErrorIgnoreLevel = ROOT.kError

def cleanup(year, finalstate, physics, cat, mass, lifetime):
    subprocess.run(["rm", f"cache.root"])
    subprocess.run(["mkdir", f"m{mass}_ct{lifetime}_{cat}_{year}_{finalstate}_{physics}"])
    subprocess.run(["mv", f"fit_bkg_m{mass}_ct{lifetime}_{cat}_{year}_fit.root", f"m{mass}_ct{lifetime}_{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"fit_bkg_m{mass}_ct{lifetime}_{cat}_{year}_gen.root", f"m{mass}_ct{lifetime}_{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"fit_sig_m{mass}_ct{lifetime}_{cat}_{year}.root", f"m{mass}_ct{lifetime}_{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"bkg_parameters_m{mass}_ct{lifetime}_{cat}_{year}_fit.json", f"m{mass}_ct{lifetime}_{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"bkg_parameters_m{mass}_ct{lifetime}_{cat}_{year}_gen.json", f"m{mass}_ct{lifetime}_{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"sig_parameters_m{mass}_ct{lifetime}_{cat}_{year}.json", f"m{mass}_ct{lifetime}_{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"data_obs_m{mass}_ct{lifetime}_{cat}_{year}.root", f"m{mass}_ct{lifetime}_{cat}_{year}_{finalstate}_{physics}/"])
#    subprocess.run(["mv", f"datacardInputs_{physics}_{finalstate}_m{mass}_ct{lifetime}_{cat}_{year}.root", f"{cat}_{year}_{finalstate}_{physics}/"])
    subprocess.run(["mv", f"rate_histos_m{mass}_ct{lifetime}_{cat}_{year}.root", f"m{mass}_ct{lifetime}_{cat}_{year}_{finalstate}_{physics}/"])
#    subprocess.run(["mv", f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{cat}_{year}.txt", f"{cat}_{year}_{finalstate}_{physics}/"])
    
def main(paths, isMC, trees, var, categories, period, bins, lifetime, mass, finalstate="4g", physics="ggH", bkg_weight=True,order_fit=order_fit, order_gen=order_gen,lxy1=lxy1, lxy2=lxy2, lumi_scaling=1):
    ROOT.gROOT.SetBatch(True)
    year=period

    teal = "\033[38;5;44m"
    reset = "\033[0m"
    print(f"{teal}processing datacard for ct={lifetime} mm, mass={mass} GeV, year={year} in categories: {categories}{reset}")

    cat_dict = {"displaced" : {"cut" : f"(best_4g_phi1_dxy_m{mass}>{lxy1})&&(best_4g_phi2_dxy_m{mass}>{lxy2})", "file" : ""},
                    "asym" : {"cut" : f"(best_4g_phi1_dxy_m{mass}>{lxy1})&&(best_4g_phi2_dxy_m{mass}<{lxy2})||(best_4g_phi1_dxy_m{mass}<{lxy1})&&(best_4g_phi2_dxy_m{mass}>{lxy2})", "file" : ""},
                    "prompt" : {"cut" : f"(best_4g_phi1_dxy_m{mass}<{lxy1})&&(best_4g_phi2_dxy_m{mass}<{lxy2})", "file" : ""},
                    "none" : {"cut" : f"(best_4g_phi1_dxy_m{mass}<{lxy1})&&(best_4g_phi2_dxy_m{mass}<{lxy2})||(best_4g_phi1_dxy_m{mass}<{lxy1})&&(best_4g_phi2_dxy_m{mass}>{lxy2})||(best_4g_phi1_dxy_m{mass}>{lxy1})&&(best_4g_phi2_dxy_m{mass}<{lxy2})||(best_4g_phi1_dxy_m{mass}>{lxy1})&&(best_4g_phi2_dxy_m{mass}>{lxy2})", "file" : ""}}


    output_names = ["rate_histos_m{}_ct{}_{}_{}.root".format(mass, lifetime, cat, year) for cat in categories]
    histo_names = []

    for i in range(len(categories)):
        pair = ["hist{}".format(2*i+1), "hist{}".format(2*i+2)]
        histo_names.append(pair)

    xsec_quad_up = np.sqrt((xsec_unc["ggH"][0])**2+(xsec_unc["VBF"][0])**2)
    xsec_quad_down = np.sqrt((xsec_unc["ggH"][1])**2+(xsec_unc["VBF"][1])**2)
    PDF_alphas_unc = np.sqrt((pdf_alphas_unc["ggH"][0])**2 + (pdf_alphas_unc["VBF"][0])**2)
        
    signal_samples={}
    for j in isMC:
        if isMC[j]==1:
            signal_samples[paths[j]]=trees[j]

    selections = []
    for cat in cat_dict:    
        if cat in categories:
            selections.append(cat_dict[cat]["cut"])

    selections.reverse()
    #print(selections)
    
    th1d_files, th1d_filenames, th1d_histos, th1d_histo_obj = datacardtools.sig_bkg_histos(paths, isMC, trees,mass,lifetime, selections, var, output_names, bins, year, histo_names, bkg_weight, lumi_scaling=lumi_scaling)

    i=0
    for cat in categories:
        N_sb=th1d_histo_obj[i][1].Integral()
        dcm_cat_year = DatacardWorkspace(finalstate, cat, period, lifetime, mass, lumi[year], physics)

        #we need to fit bkg twice because one fit is used to generate the data and one fit is used for combine fit
        #two jsons with parameters: fit and gen. gen=used to generate data. fit=used in combine.
        ratio = ggHfitter.fitBKG(f"{th1d_filenames[i]}", f"{th1d_histos[i][1]}", f"fit_bkg_m{mass}_ct{lifetime}_{cat}_{year}_fit.root", order=order_fit)
        datacardtools.extract_JSON(f"fit_bkg_m{mass}_ct{lifetime}_{cat}_{year}_fit.root", "w", f"bkg_parameters_m{mass}_ct{lifetime}_{cat}_{year}_fit.json")

        ggHfitter.fitBKG(f"{th1d_filenames[i]}", f"{th1d_histos[i][1]}", f"fit_bkg_m{mass}_ct{lifetime}_{cat}_{year}_gen.root", order=order_gen)
        datacardtools.extract_JSON(f"fit_bkg_m{mass}_ct{lifetime}_{cat}_{year}_gen.root", "w", f"bkg_parameters_m{mass}_ct{lifetime}_{cat}_{year}_gen.json")

        bkg_rate=N_sb*ratio

        ggHfitter.fitSIG(f"{th1d_filenames[i]}", f"{th1d_histos[i][0]}", f"fit_sig_m{mass}_ct{lifetime}_{cat}_{year}.root")
        datacardtools.extract_JSON(f"fit_sig_m{mass}_ct{lifetime}_{cat}_{year}.root", "w", f"sig_parameters_m{mass}_ct{lifetime}_{cat}_{year}.json")

        dcm_cat_year.addDCB("signal", "mass", f"sig_parameters_m{mass}_ct{lifetime}_{cat}_{year}.json", resolution={f"nuisance_smear_m{mass}_ct{lifetime}_{cat}_{year}":str(smear_resolution)})
        dcm_cat_year.addBernstein("background", "mass", f"bkg_parameters_m{mass}_ct{lifetime}_{cat}_{year}_fit.json")

        dcm_cat_year.addSystematic(name=f"nuisance_smear_m{mass}_ct{lifetime}_{cat}_{year}", kind = "param", values=[0.0, 1.0])
        dcm_cat_year.addSystematic(name=f"xsec_unc_m{mass}_ct{lifetime}_{cat}_{year}", kind = "lnN", values={"signal":f"{1+xsec_quad_up}/{1-xsec_quad_down}"})
        dcm_cat_year.addSystematic(name=f"lumi_unc_m{mass}_ct{lifetime}_{cat}_{year}", kind = "lnN", values={"signal":"{}".format(1+lumi_unc[year])})
        dcm_cat_year.addSystematic(name=f"PDF_alphas_unc_m{mass}_ct{lifetime}_{cat}_{year}", kind = "lnN", values={"signal":"{}".format(1+PDF_alphas_unc)})
        dcm_cat_year.addSystematic(name=f"bkg_rate_m{mass}_ct{lifetime}_{cat}_{year}", kind = "rateParam", values=[dcm_cat_year.tag, "background", "1", "[0,10]"])

        dcm_cat_year.addFixedYield(name="background", ID=1, value=bkg_rate)
        dcm_cat_year.addFixedYieldFromFile(name="signal", ID=0, filename=th1d_filenames[i], histoName=th1d_histos[i][0], lumi=True)

        workspace_file_cat_year = "datacardInputs_"+dcm_cat_year.tag+".root"

        data_cat_year_name = datacardtools.generate_data_hist(f"fit_bkg_m{mass}_ct{lifetime}_{cat}_{year}_gen.root", bins_num=bins[0], norm=bkg_rate, output_name=f"data_obs_m{mass}_ct{lifetime}_{cat}_{year}.root")

        dcm_cat_year.importBinnedData(data_cat_year_name, "h_pdf__mass", ["mass"])
            
        dcm_cat_year.makeCard()

        cleanup(year, finalstate, physics, cat, mass, lifetime)

        i+=1

