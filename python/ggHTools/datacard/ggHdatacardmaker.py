from datacard import datacardtools
from datacard import ggHfitter
from datacard.ggHdatacardworkspace import DatacardWorkspace
from ggHparameters import order_fit, order_gen, smear_resolution, lumi, xsec_unc, pdf_alphas_unc, lumi_unc
import ROOT
import subprocess
import numpy as np

ROOT.gROOT.SetBatch(False)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
ROOT.gErrorIgnoreLevel = ROOT.kError

# helper function to get the correct photon id by year (loose=run2, medium=run3)
def recommended_photon_id(period):
    period = str(period)
    if period in ["2016", "2017", "2018", "Run2"]:
        return "LooseEGM"
    if period in ["2022", "2023", "2024", "2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "Run3"]:
        return "MediumEGM"
    raise ValueError(f"no recommended photon ID is configured for period '{period}'")

def cleanup(year, finalstate, physics, mass, lifetime):
    output_dir = f"m{mass}_ct{lifetime}_{year}_{finalstate}_{physics}"
    subprocess.run(["rm", f"cache.root"])
    subprocess.run(["mkdir", output_dir])
    subprocess.run(["mv", f"fit_bkg_m{mass}_ct{lifetime}_{year}_fit.root", f"{output_dir}/"])
    subprocess.run(["mv", f"fit_bkg_m{mass}_ct{lifetime}_{year}_gen.root", f"{output_dir}/"])
    subprocess.run(["mv", f"fit_sig_m{mass}_ct{lifetime}_{year}.root", f"{output_dir}/"])
    subprocess.run(["mv", f"bkg_parameters_m{mass}_ct{lifetime}_{year}_fit.json", f"{output_dir}/"])
    subprocess.run(["mv", f"bkg_parameters_m{mass}_ct{lifetime}_{year}_gen.json", f"{output_dir}/"])
    subprocess.run(["mv", f"sig_parameters_m{mass}_ct{lifetime}_{year}.json", f"{output_dir}/"])
    subprocess.run(["mv", f"data_obs_m{mass}_ct{lifetime}_{year}.root", f"{output_dir}/"])
    subprocess.run(["mv", f"rate_histos_m{mass}_ct{lifetime}_{year}.root", f"{output_dir}/"])
    
def main(paths, isMC, trees, var, period, bins, lifetime, mass,finalstate="4g", physics="ggH", order_fit=order_fit,order_gen=order_gen, lumi_scaling=1, signal_lumis=None):
    ROOT.gROOT.SetBatch(True)
    year=period
    photon_id = recommended_photon_id(period)

    teal = "\033[38;5;44m"
    reset = "\033[0m"
    print(f"{teal}processing datacard for ct={lifetime} mm, mass={mass} GeV, year={year}, photon ID={photon_id}{reset}")

    output_name = f"rate_histos_m{mass}_ct{lifetime}_{year}.root"
    histo_names = [f"hist_{i}" for i in range(len(paths))]

    xsec_quad_up = np.sqrt((xsec_unc["ggH"][0])**2+(xsec_unc["VBF"][0])**2)
    xsec_quad_down = np.sqrt((xsec_unc["ggH"][1])**2+(xsec_unc["VBF"][1])**2)
    PDF_alphas_unc = np.sqrt((pdf_alphas_unc["ggH"][0])**2 + (pdf_alphas_unc["VBF"][0])**2)
        
    signal_indices = [i for i, flag in enumerate(isMC) if flag]
    background_indices = [i for i, flag in enumerate(isMC) if not flag]
    if not signal_indices:
        raise ValueError("makes no sense. >=1 signal input is required")
    if len(background_indices) != 1:
        raise ValueError("onl 1 data/background input is required")
    background_index = background_indices[0]
    if signal_lumis is not None and len(signal_lumis) != len(signal_indices):
        raise ValueError("signal_lumis must  only have one value per signal input")
    file_scalings = None
    if signal_lumis is not None:
        file_scalings = []
        signal_number = 0
        for flag in isMC:
            if flag:
                file_scalings.append(signal_lumis[signal_number])
                signal_number += 1
            else:
                file_scalings.append(1.0)

#generates the signal and background histograms
    th1d_filename, th1d_histos, th1d_histo_obj = datacardtools.sig_bkg_histos(paths, isMC, trees, mass, var, output_name, bins, photon_id, histo_names=histo_names, lumi_scaling=lumi_scaling, file_scalings=file_scalings)

    N_sb_raw=float(th1d_histo_obj[background_index].Integral())
    N_sb=int(round(N_sb_raw))
    if not np.isclose(N_sb_raw, N_sb, rtol=0.0, atol=1e-9):
        raise ValueError("background sideband yield must be an unweighted integer count, found {}".format(N_sb_raw))
    dcm_year = DatacardWorkspace(finalstate, period, lifetime, mass, lumi[year], physics)

    signal_hist_name = "signal_combined"
    signal_hist = th1d_histo_obj[signal_indices[0]].GetValue().Clone(signal_hist_name)
    signal_hist.SetDirectory(0)
    for signal_index in signal_indices[1:]:
        signal_hist.Add(th1d_histo_obj[signal_index].GetValue())
    rate_file = ROOT.TFile(th1d_filename, "UPDATE")
    rate_file.cd()
    signal_hist.Write(signal_hist_name, ROOT.TObject.kOverwrite)
    rate_file.Close()

# i have to fit bkg twice. 1: to generate fake data 2: actual combine bkg fit
    ratio = ggHfitter.fitBKG(th1d_filename, th1d_histos[background_index], f"fit_bkg_m{mass}_ct{lifetime}_{year}_fit.root", order=order_fit)
    if not np.isfinite(ratio) or ratio <= 0.0:
        raise ValueError("fitted SR/sideband extrapolation ratio must be finite and positive, found {}".format(ratio))
    datacardtools.extract_JSON(f"fit_bkg_m{mass}_ct{lifetime}_{year}_fit.root", "w", f"bkg_parameters_m{mass}_ct{lifetime}_{year}_fit.json")
    ggHfitter.fitBKG(th1d_filename, th1d_histos[background_index], f"fit_bkg_m{mass}_ct{lifetime}_{year}_gen.root", order=order_gen)
    datacardtools.extract_JSON(f"fit_bkg_m{mass}_ct{lifetime}_{year}_gen.root", "w", f"bkg_parameters_m{mass}_ct{lifetime}_{year}_gen.json")

    bkg_rate=N_sb*ratio

    ggHfitter.fitSIG(th1d_filename, signal_hist_name, f"fit_sig_m{mass}_ct{lifetime}_{year}.root")
    datacardtools.extract_JSON(f"fit_sig_m{mass}_ct{lifetime}_{year}.root", "w", f"sig_parameters_m{mass}_ct{lifetime}_{year}.json")

#import the DCB and Bernstein polynomial into the workspace
    dcm_year.addDCB("signal", "mass", f"sig_parameters_m{mass}_ct{lifetime}_{year}.json", resolution={f"nuisance_smear_m{mass}_ct{lifetime}_{year}":str(smear_resolution)})
    dcm_year.addBernstein("background", "mass", f"bkg_parameters_m{mass}_ct{lifetime}_{year}_fit.json")

#add systematics
    dcm_year.addSystematic(name=f"nuisance_smear_m{mass}_ct{lifetime}_{year}", kind = "param", values=[0.0, 1.0])
    dcm_year.addSystematic(name=f"xsec_unc_m{mass}_ct{lifetime}_{year}", kind = "lnN", values={"signal":f"{1-xsec_quad_down}/{1+xsec_quad_up}"})
    dcm_year.addSystematic(name=f"lumi_unc_m{mass}_ct{lifetime}_{year}", kind = "lnN", values={"signal":"{}".format(1+lumi_unc[year])})
    dcm_year.addSystematic(name=f"PDF_alphas_unc_m{mass}_ct{lifetime}_{year}", kind = "lnN", values={"signal":"{}".format(1+PDF_alphas_unc)})
    #dcm_year.addSystematic(name=f"bkg_rate_m{mass}_ct{lifetime}_{year}", kind = "rateParam", values=[dcm_year.tag, "background", "1", "[0,10]"])
    dcm_year.addGmN(name=f"bkg_sideband_stat_m{mass}_ct{lifetime}_{year}", count=N_sb, values={"background":ratio})

#yields are added by integrating the histograms for sig and bkg
    dcm_year.addFixedYield(name="background", ID=1, value=bkg_rate)
    dcm_year.addFixedYieldFromFile(name="signal", ID=0, filename=th1d_filename, histoName=signal_hist_name, lumi=(signal_lumis is None))

#this part is just for creating fake data inside the SR (see datacardtools.py for it_
    data_year_name = datacardtools.generate_data_hist(f"fit_bkg_m{mass}_ct{lifetime}_{year}_gen.root", bins_num=bins[0], norm=bkg_rate, output_name=f"data_obs_m{mass}_ct{lifetime}_{year}.root")

# add in the fake data to the card, make the card and clean up all the file junk
    dcm_year.importBinnedData(data_year_name, "h_pdf__mass", ["mass"])
    dcm_year.makeCard()
    cleanup(year, finalstate, physics, mass, lifetime)
