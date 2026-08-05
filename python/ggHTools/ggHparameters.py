order_fit = 4
order_gen = 4

eos_base = "/eos/uscms/store/user/dacampos/analysis"

def signal_path(mass, ctau, year):
    return f"{eos_base}/signal/ggH4g_M{mass}_ctau{ctau}_{year}_0_ggH4g_M{mass}_ctau{ctau}_{year}_ggH4g.root"

def bkg_path(year):
    return f"{eos_base}/data/EGamma_{year}_updated/EGamma_{year}_all_ggH4g.root"

signal_window = (110, 140)

fit_window = (70, 180)
lower_sb = (70, 110)
upper_sb = (140, 180)

n_bins = 30
bins = [n_bins, signal_window[0], signal_window[1]]

bin_width = (signal_window[1] - signal_window[0]) / n_bins
n_fit_bins = int(round((fit_window[1] - fit_window[0]) / bin_width))
fit_bins = [n_fit_bins, fit_window[0], fit_window[1]]

lxy1 = 50
lxy2 = 50

dxy_min = -20

signal_xsec = 52.143
BR = 1e-4

smear_resolution = 0.264

##luminosities in pb^-1 from: https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2
lumi = {"2017": 41480,
        "2018": 59830,
        "Run2": 41480 + 59830,
        "2022preEE": 7990,
        "2022postEE": 26680,
        "2022": 7990 + 26680,
        "2023preBPix": 17960,
        "2023postBPix": 9680,
        "2023": 17960 + 9680,
        #2024 = 110.11 +/- 1.77 fb^-1: https://indico.cern.ch/event/1617597/contributions/6820152/attachments/3184909/5675974/Run%202%20and%20Run%203%20combination.pdf
        "2024": 110110,
        "Run3": (7990 + 26680) + (17960 + 9680) + 110110}

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
