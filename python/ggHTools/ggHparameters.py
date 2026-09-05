order_fit = 4
order_gen = 4

#obviously this will need to change if someone else is runing ggH
eos_base = "/eos/uscms/store/user/dacampos/analysis"

signal_masses = [15, 20, 30, 40, 50, 55]
signal_ctaus = [0, 10, 20, 50, 100, 1000]
run2_signal_years = ["2017", "2018"]
run3_signal_years=["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]
run2_data_years=["2017","2018"]
run3_data_years=["2022","2023","2024"]
combined_data_years=("Run2","Run3")
year_config={
    "2017":{"data_years":["2017"],"signal_years":["2017"],"period":4},
    "2018":{"data_years":["2018"],"signal_years":["2018"],"period":4},
    "Run2":{"data_years":run2_data_years,"signal_years":run2_signal_years,"period":4},
    "2022":{"data_years":["2022"],"signal_years":["2022preEE","2022postEE"],"period":5},
    "2023":{"data_years":["2023"],"signal_years":["2023preBPix","2023postBPix"],"period":5},
    "2024":{"data_years":["2024"],"signal_years":["2024"],"period":5},
    "Run3":{"data_years":run3_data_years,"signal_years":run3_signal_years,"period":5},
}
# helper function to get the correct photon id by year (loose=run2, medium=run3)
def recommended_photon_id(period):
    period = str(period)
    if period in ["2016", "2017", "2018", "Run2"]:
        return "LooseEGM"
    if period in ["2022", "2023", "2024", "2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "Run3"]:
        return "MediumEGM"
    raise ValueError(f"no recommended photon ID is configured for period '{period}'")

delta_r_cut = 0.3

def signal_path(mass, ctau, year):
    return f"{eos_base}/signal_newTrig/ggH4g_M{mass}_ctau{ctau}_{year}_0_ggH4g_M{mass}_ctau{ctau}_{year}_ggH4g.root"

def bkg_path(year):
    directory=year if year in combined_data_years else f"EGamma_{year}"
    return f"{eos_base}/data_newTrig/{directory}/EGamma_{year}_ggH4g_all.root"

signal_window = (110, 140)
fit_window=(100,220)

lower_sb=(fit_window[0],signal_window[0])
upper_sb=(signal_window[1],fit_window[1])

analysis_tree="ggH4g"


n_bins = 30
bins = [n_bins, signal_window[0], signal_window[1]]

bin_width = (signal_window[1] - signal_window[0]) / n_bins
n_fit_bins = int(round((fit_window[1] - fit_window[0]) / bin_width))
fit_bins = [n_fit_bins, fit_window[0], fit_window[1]]

summary_bin_width = 5.0

lxy1=50
lxy2 = 50

dxy_min=-20 #cm
dxy_max= 110 #cm

signal_xsec = 52.143
BR = 1e-4

smear_resolution = 0.264

##luminosities in pb^-1 from: https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2
##all values computed with brilcalc using goldenJSON for each year
lumi = {"2017": 27461.448215758, # we only use eraas D-F, no longer B,C. not the full 41.48 fb^-1 for 2017
        "2018": 59557.110211607,
        "Run2": 27461.448215758 + 59557.110211607,
        "2022preEE": 8078.592654770,
        "2022postEE": 26671.382348627,
        "2022": 8078.592654770 + 26671.382348627,
        "2023preBPix": 18579.703603083,
        "2023postBPix": 9675.416900571,
        "2023": 18579.703603083 + 9675.416900571,
        "2024": 109143.514819334,
        "Run3": (8078.592654770 + 26671.382348627) + (18579.703603083 + 9675.416900571) + 109143.514819334}

#######IMPORTANT#######
#since the cross section we use to scale the MC is both ggH+VBF combined, these xsec and lumi unc are derived from adding the unc from
#each process in quadrature

#obtained xsec and PDF from:https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV#gluon_gluon_Fusion_Process
xsec_unc = {"ggH": [0.046, -0.067],
            "VBF": [0.004, -0.003]}

pdf_alphas_unc = {"ggH": [0.032, 0.032],
                  "VBF": [0.021, 0.021]}

#https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun3
#https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2
lumi_unc = {"2017": 0.008072175,
            "2018": 0.008409014,
            "Run2": 0.011656,
            "2022preEE": 0.0136,  #assuming suberas use the full-year value
            "2022postEE": 0.0136,
            "2022": 0.0136,
            "2023preBPix": 0.0115,  #assuming suberas use the full-year value
            "2023postBPix": 0.0115,
            "2023": 0.0115,
            "2024": 0.0161,
            "Run3": 0.011947}

dcb_mean = (125, 120, 130)
dcb_sigma = (2, 0.1, 20)
dcb_alpha1 = (2, 1, 20)
dcb_n1 = (2, 1, 50)
dcb_alpha2 = (2, 1, 20)
dcb_n2 = (2, 1, 50)

bernstein_coeff = (0, 10000)
bernstein_coeff_card = (0, 10000)
