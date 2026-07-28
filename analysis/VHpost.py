import os
import ROOT
from pathlib import Path
import array
import pandas as pd
import subprocess
from scipy.stats import beta
import matplotlib
matplotlib.use('TkAgg')
from common.plotter import *
from common.datacard_maker import cnc_datacard_maker
import numpy as np
import pickle        
import matplotlib.colors as mpcolors 

ROOT.gInterpreter.Declare('#include "common/chelpers.h"')
ROOT.gInterpreter.Declare('#include "common/signalEfficiency.h"')
ROOT.gInterpreter.Declare('#include "common/scaleFactors.h"')

ROOT.gInterpreter.Declare('#include "common/vhFakeRates.h"')        
ROOT.gInterpreter.Declare('#include "analysis/ddp_vertex.h"')        
ROOT.ROOT.EnableImplicitMT()

analysis_status='Supplementary'

masses = [15,20,30,40,50,55]
lifetimes = [0,10,20,50,100,1000]
signal_colors={0:'#d31a17',
               10:'#005a9a',
               20:'chocolate',
               50:'darkorange',
               100:'darkgreen',
               1000:'#282d65'}
signal_colors_m={15:'#d31a17',
               20:'chocolate',
               30:'#005a9a',                 
               40:'darkorange',
               50:'#282d65',
               55:'darkgreen'}

analyses  = ['wmn2g','wen2g','zmm2g','zee2g']
signals=['ZH','ggZH','WH','ttH']

lumifb = {'2018': "60",
          '2017': '41',
          '2016': '36',
          'Run2': '138'}

center_of_mass = {'2018': 13,
                  '2017': 13,
                  '2016': 13,
                  'Run2': 13}

#Colors:
cmscolors = ["#218fcf", "#ffde9c",  "#98c565","#e5e9ee", 
             "#e62d1a", "#282d65", "#005a9a", "#85d2fb", "#3cbee9"]
mh.style.use(["CMS", {"axes.prop_cycle": cycler("color", cmscolors)}])


#FONTS
#plt.rcParams['text.usetex'] = True
# Optional: Ensure the rest of the text matches your math style
#plt.rcParams['font.family'] = 'serif'


# Standard Analysis Cuts
cuts = {}
# Tight + veto lepton cuts
cuts['W'] = {'MU': "(Muon_isTrigger[W_l1_idx] && Muon_nloose==1 && Electron_nloose==0)",
             'ELE': "(Electron_isTrigger[W_l1_idx] && Muon_nloose==0 && Electron_nloose==1)"}
            
cuts['Z'] = {'MU': "((Muon_isTrigger[Z_idx[0]] || Muon_isTrigger[Z_idx[1]]) && Muon_nloose==2 && Electron_nloose==0)",
             'ELE': "((Electron_isTrigger[Z_idx[0]] || Electron_isTrigger[Z_idx[1]]) && Electron_nloose==2 && Muon_nloose==0)"}
# Photon pt and kinematic cuts
cuts['pt'] = {}
cuts['photons'] = {}
for m in masses:
    cuts['pt'][m] = "(Photon_pt[best_2g_idx1_m{m}] > 35 && Photon_pt[best_2g_idx2_m{m}] > 25)".format(m=m)
    cuts['photons'][m] = "(best_2g_pt_m{m} > 20 && best_2g_raw_mass_m{m} > 4)".format(m=m)

# w->e+nu misID cuts
cuts['misID'] = {}
for v in ['W', 'Z']:
    cuts['misID'][v] = {}
    for l in ['ELE', 'MU']:
        cuts['misID'][v][l] = {}
        for m in masses:
            if (v == 'W' and l == 'ELE'):
                cuts['misID'][v][l][m] = "(abs(best_2g_misID1_m{m}-91) > 10 && abs(best_2g_misID2_m{m}-91) > 10 && abs(best_2g_misID3_m{m} - 91) > 15)".format(m=m)
            else:
                cuts['misID'][v][l][m] = "(1)"
# Z->ll FSR cuts
#cuts['fsr'] = {}
#for v in ['W', 'Z']:
#    cuts['fsr'][v] = {}
#    for m in masses:
#        cuts['fsr'][v][m] = "(1)"
#        if (v=='W'):
#            cuts['fsr'][v][m] = "(1)"
#        else:
#            cuts['fsr'][v][m] = "(Photon_isFSR[best_2g_idx1_m{m}]==0 && Photon_isFSR[best_2g_idx2_m{m}]==0)".format(m=m)

            
#redefine the cuts here based on the analysis for easier parsing
for analysis in ['wmn2g','wen2g','zmm2g','zee2g']:
    cuts[analysis]={}
    if analysis=='wmn2g':
        v='W'
        l='MU'
    elif analysis=='zmm2g':
        v='Z'
        l='MU'
    elif analysis=='wen2g':
        v='W'
        l='ELE'
    elif analysis=='zee2g':
        v='Z'
        l='ELE'
    for m in masses:
        cuts[analysis][m] = {} 
        cuts[analysis][m]['presr']='&&'.join([cuts[v][l],
                                           cuts['pt'][m],
                                           cuts['photons'][m],
                                           cuts['misID'][v][l][m],
                                           f"((Photon_passCutBasedID[best_2g_idx1_m{m}]+Photon_passCutBasedID[best_2g_idx2_m{m}])==2)"])
        cuts[analysis][m]['precr']='&&'.join([cuts[v][l],
                                              cuts['pt'][m],
                                              cuts['photons'][m],
                                              cuts['misID'][v][l][m],
                                              f"(Photon_passCutBasedID[best_2g_idx1_m{m}]>0)"])
        cuts[analysis][m]['noid']='&&'.join([cuts[v][l],
                                              cuts['pt'][m],
                                              cuts['photons'][m],
                                             cuts['misID'][v][l][m]])
        
        cuts[analysis][m]['precr_abcd']='&&'.join([cuts[v][l],
                                              cuts['pt'][m],
                                              cuts['photons'][m],
                                              cuts['misID'][v][l][m],
                                              f"((Photon_passCutBasedID[best_2g_idx1_m{m}]+Photon_passCutBasedID[best_2g_idx2_m{m}])==1)"])


        cuts[analysis][m]['sr']='&&'.join([cuts[analysis][m]['presr'],f'(best_2g_dxy_m{m}>-10)'])
        cuts[analysis][m]['cr']='&&'.join([cuts[analysis][m]['precr'],f'(best_2g_dxy_m{m}>-10)'])
        cuts[analysis][m]['cr_abcd']='&&'.join([cuts[analysis][m]['precr_abcd'],f'(best_2g_dxy_m{m}>-10)'])
        cuts[analysis][m]['ssb']='&&'.join([cuts[analysis][m]['presr'],f'(best_2g_raw_mass_m{m}>62.5)'])
        cuts[analysis][m]['csb_abcd']='&&'.join([cuts[analysis][m]['precr_abcd'],f'(best_2g_raw_mass_m{m}>62.5)'])                  
# Cuts for signal efficiencies:
effCuts = {}
effCuts['HLT'] = {'ELE': {'2018': "(HLT_Ele32_WPTight_Gsf)",
                          '2017': "(HLT_Ele32_WPTight_Gsf_L1DoubleEG)",
                          '2016': "(HLT_Ele27_WPTight_Gsf)"},
                  'MU': {'2018': "(HLT_IsoMu24)",
                         '2017': "(HLT_IsoMu27)",
                         '2016': "(HLT_IsoMu24)"}
              }

# Cuts for W->e+nu 2018 HEM
cutsHEM = "(Electron_eta[W_l1_idx]>-1.3 || Electron_eta[W_l1_idx]<-3.0) && (Electron_phi[W_l1_idx]>-0.85 || Electron_phi[W_l1_idx]<-1.57)"








binning={'wen2g': {7: [((4.0, 6.0), [-10.0, 40.0, 70.0, 110.0]), ((6.0, 8.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((8.0, 10.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((10.0, 12.0), [-10.0, 20.0, 110.0])], 15: [((4.0, 8.0), [-10.0, 40.0, 70.0, 110.0]), ((8.0, 12.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((12.0, 16.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((16.0, 20.0), [-10.0, 20.0, 110.0])], 20: [((4.0, 9.25), [-10.0, 40.0, 70.0, 110.0]), ((9.25, 14.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((14.5, 19.75), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((19.75, 25.0), [-10.0, 20.0, 110.0])], 30: [((4.0, 11.75), [-10.0, 40.0, 70.0, 110.0]), ((11.75, 19.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((19.5, 27.25), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((27.25, 35.0), [-10.0, 20.0, 110.0])], 40: [((4.0, 14.25), [-10.0, 40.0, 110.0]), ((14.25, 24.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((24.5, 34.75), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((34.75, 45.0), [-10.0, 20.0, 110.0])], 50: [((4.0, 16.75), [-10.0, 40.0, 70.0, 110.0]), ((16.75, 29.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((29.5, 42.25), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((42.25, 55.0), [-10.0, 20.0, 110.0])], 55: [((4.0, 18.0), [-10.0, 40.0, 70.0, 110.0]), ((18.0, 32.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((32.0, 46.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((46.0, 60.0), [-10.0, 20.0, 110.0])]}, 'wmn2g': {7: [((4.0, 6.0), [-10.0, 40.0, 70.0, 110.0]), ((6.0, 8.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((8.0, 10.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((10.0, 12.0), [-10.0, 20.0, 110.0])], 15: [((4.0, 8.0), [-10.0, 40.0, 70.0, 110.0]), ((8.0, 12.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((12.0, 16.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((16.0, 20.0), [-10.0, 20.0, 110.0])], 20: [((4.0, 9.25), [-10.0, 40.0, 70.0, 110.0]), ((9.25, 14.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((14.5, 19.75), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((19.75, 25.0), [-10.0, 20.0, 110.0])], 30: [((4.0, 11.75), [-10.0, 15.0, 70.0, 110.0]), ((11.75, 19.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((19.5, 27.25), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((27.25, 35.0), [-10.0, 20.0, 110.0])], 40: [((4.0, 14.25), [-10.0, 40.0, 70.0, 110.0]), ((14.25, 24.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((24.5, 34.75), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((34.75, 45.0), [-10.0, 20.0, 110.0])], 50: [((4.0, 16.75), [-10.0, 40.0, 70.0, 110.0]), ((16.75, 29.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((29.5, 42.25), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((42.25, 55.0), [-10.0, 20.0, 110.0])], 55: [((4.0, 18.0), [-10.0, 40.0, 70.0, 110.0]), ((18.0, 32.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((32.0, 46.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((46.0, 60.0), [-10.0, 20.0, 110.0])]}, 'zee2g': {7: [((4.0, 8.0), [-10.0, 26.0, 62.0, 110.0]), ((8.0, 12.0), [-10.0, 110.0])], 15: [((4.0, 12.0), [-10.0, 26.0, 62.0, 110.0]), ((12.0, 20.0), [-10.0, 110.0])], 20: [((4.0, 14.5), [-10.0, 110.0]), ((14.5, 25.0), [-10.0, 18.0, 45.0, 72.0, 110.0])], 30: [((4.0, 11.75), [-10.0, 40.0, 70.0, 110.0]), ((11.75, 19.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((19.5, 27.25), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((27.25, 35.0), [-10.0, 20.0, 110.0])], 40: [((4.0, 14.25), [-10.0, 40.0, 110.0]), ((14.25, 24.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((24.5, 34.75), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((34.75, 45.0), [-10.0, 10.0, 30.0, 110.0])], 50: [((4.0, 16.75), [-10.0, 40.0, 110.0]), ((16.75, 29.5), [-10.0, 18.0, 72.0, 110.0]), ((29.5, 42.25), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((42.25, 55.0), [-10.0, 20.0, 110.0])], 55: [((4.0, 18.0), [-10.0, 40.0, 110.0]), ((18.0, 32.0), [-10.0, 18.0, 72.0, 110.0]), ((32.0, 46.0), [-10.0, 18.0, 45.0, 110.0]), ((46.0, 60.0), [-10.0, 20.0, 110.0])]}, 'zmm2g': {7: [((4.0, 6.0), [-10.0, 40.0, 70.0, 110.0]), ((6.0, 8.0), [-10.0, 45.0, 72.0, 110.0]), ((8.0, 10.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((10.0, 12.0), [-10.0, 20.0, 110.0])], 15: [((4.0, 8.0), [-10.0, 40.0, 70.0, 110.0]), ((8.0, 12.0), [-10.0, 45.0, 72.0, 110.0]), ((12.0, 16.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((16.0, 20.0), [-10.0, 20.0, 110.0])], 20: [((4.0, 9.25), [-10.0, 40.0, 70.0, 110.0]), ((9.25, 14.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((14.5, 19.75), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((19.75, 25.0), [-10.0, 20.0, 110.0])], 30: [((4.0, 11.75), [-10.0, 40.0, 110.0]), ((11.75, 19.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((19.5, 27.25), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((27.25, 35.0), [-10.0, 20.0, 110.0])], 40: [((4.0, 14.25), [-10.0, 40.0, 70.0, 110.0]), ((14.25, 24.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((24.5, 34.75), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((34.75, 45.0), [-10.0, 10.0, 30.0, 110.0])], 50: [((4.0, 16.75), [-10.0, 40.0, 70.0, 110.0]), ((16.75, 29.5), [-10.0, 18.0, 72.0, 110.0]), ((29.5, 42.25), [-10.0, 45.0, 72.0, 110.0]), ((42.25, 55.0), [-10.0, 20.0, 110.0])], 55: [((4.0, 18.0), [-10.0, 40.0, 70.0, 110.0]), ((18.0, 32.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((32.0, 46.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((46.0, 60.0), [-10.0, 20.0, 110.0])]}}


binning1d={'wen2g':np.array([-10.,10.,30.,50.,75.,110.]),'wmn2g':np.array([-10.,10.,30.,50.,75.,110.]),'zmm2g':np.array([-10.,30.,75.,110.]),'zee2g':np.array([-10.,30.,75.,110.])}
binning1d_sb={'wen2g':np.array([-160.,-130.,-110.,-90.,-70.,-50.,-30.,-10.]),'wmn2g':np.array([-160.,-130.,-110.,-90.,-70.,-50.,-30.,-10.]),'zmm2g':np.array([-160.,-110.,-60.,-10.]),'zee2g':np.array([-160.,-110.,-60.,-10.])}
binning_mass={'wen2g':np.array([4.,12.,20.,28.,36.,44.,52.,60.,68.]),'wmn2g':np.array([4.,12.,20.,28.,36.,44.,52.,60.,68.]),'zmm2g':np.array([4.,20.,44.,68.]),'zee2g':np.array([4.,20.,44.,68.])}


# Scale factor stuff 
leptonSF = {}
leptonSF['wmn2g'] = "*".join(['Muon_recoSF_val[W_l1_idx]', 'Muon_idSF_val[W_l1_idx]', 'Muon_isoSF_val[W_l1_idx]', 'Muon_trigSF_val[W_l1_idx]'])
leptonSF['wen2g'] = "*".join(['Electron_recoSF_val[W_l1_idx]', 'Electron_idSF_val[W_l1_idx]', 'Electron_trigSF_val[W_l1_idx]'])
leptonSF['zee2g'] = "*".join(['Electron_recoSF_val[Z_idx[0]]', 'Electron_recoSF_val[Z_idx[1]]', 'Electron_idSF_val[Z_idx[0]]', 'Electron_idSF_val[Z_idx[1]]', 'Electron_trigSF_val[Z_idx[0]]'])
leptonSF['zmm2g'] = "*".join(['Muon_recoSF_val[Z_idx[0]]', 'Muon_recoSF_val[Z_idx[1]]', 'Muon_idSF_val[Z_idx[0]]', 'Muon_idSF_val[Z_idx[1]]', 'Muon_isoSF_val[Z_idx[0]]', 'Muon_isoSF_val[Z_idx[1]]', 'Muon_trigSF_val[Z_idx[0]]'])

    
photonSF={}
for m in masses:
    photonSF[m] = "*".join(['Photon_idSF_val[best_2g_idx1_m{}]'.format(m), 'Photon_idSF_val[best_2g_idx2_m{}]'.format(m), 'Photon_pixSF_val[best_2g_idx1_m{}]'.format(m), 'Photon_pixSF_val[best_2g_idx2_m{}]'.format(m)])




def convertOldCustomBins(customBins):
    from ctypes import c_int, byref
    d={}
    for v in ['W','Z']:
        for l in ['ELE','MU']:
            if v =='W':
                if l=='ELE':
                    ana='wen2g'
                else:    
                    ana='wmn2g'
            else: 
                if l=='ELE':
                    ana='zee2g'
                else:    
                    ana='zmm2g'
            d[ana]={}        
            for mass in [7,15,20,30,40,50,55]:
                #create a fake histogram
                if v == 'Z' and mass == 7 and l == 'ELE':
                    binsM = 2
                elif v == 'Z' and mass == 15 and l == 'ELE':
                    binsM = 2
                elif v == 'Z' and mass ==20 and l == 'ELE':
                    binsM = 2
                else:
                    binsM =4
                models=[]    
                h=ROOT.TH2D("test","test",binsM,4,mass+5,110,-10,100)
                hU = ROOT.TH1D("test2","test2",binsM*110,0,binsM*110)
                #loop on the mass bins
                for i in range(1,h.GetNbinsX()+1):
                    xedges=(h.GetXaxis().GetBinLowEdge(i),h.GetXaxis().GetBinUpEdge(i))
                    #now we need to identify the y bins associated with this xbin
                    yedges=[]
                    for b in customBins[v][l][mass][:-1]:
                        binglobal = hU.GetXaxis().FindBin(b)
                        bx = int((binglobal-1)/110)+1
                        if bx!=i:
                            continue
                        by = int((binglobal-1)%110)+1
                        yedges.append(h.GetYaxis().GetBinLowEdge(by))
                    yedges.append(110.)
                    models.append((xedges,yedges))
                d[ana][mass]=models
    print(d)            
                        
                    
                               


# Histogram methods for data cards
def unfoldTH2(hist):
    binsX = hist.GetNbinsX()
    binsY = hist.GetNbinsY()
    hOut = ROOT.TH1D("hOut", "", binsX*binsY, 0, binsX*binsY)
    for x in range(binsX):
        for y in range(binsY):
            hOut.SetBinContent(1+y+x*binsY, hist.GetBinContent(x+1, y+1))
            hOut.SetBinError(1+y+x*binsY, hist.GetBinError(x+1, y+1))
    return hOut


def isValidFile(fname, tree):
    f = ROOT.TFile.Open(fname)  # ✅ works for root:// and local paths
    if not f or f.IsZombie():
        return False
    try:
        t = f.Get(tree)
        if not t:
            f.Close()
            return False
        t.GetEntry(0)  # test if it’s readable
        _ = t.Photon_pt  # test access
    except Exception:
        f.Close()
        return False
    f.Close()
    return True







def getFiles(query,sampleDir,sampleType,era,prod):
    if '/store' in sampleDir:
        ser = pd.Series(subprocess.check_output(['xrdfs', 'root://cmseos.fnal.gov', 'ls', f"{sampleDir}/{sampleType}{era}_{prod}/"], text=True).split("\n"))
        return(list('root://cmseos.fnal.gov/' + ser[ser.str.contains(query)]))
    else:
        # Your original variables setup
        target_dir = f"{sampleDir}/{sampleType}{era}_{prod}/"
        match_str = f"{query}"

        files = []

        # If the directory doesn't exist, avoid scanning entirely
        if os.path.exists(target_dir):
            # 'with os.scandir' forces the OS to immediately release 
            # the directory descriptor the moment this block finishes!
            with os.scandir(target_dir) as entries:
                for entry in entries:
                    # Check if it's a file, ends with .root, and contains your query
                    if entry.is_file() and entry.name.endswith(".root") and match_str in entry.name:
                        files.append(entry.path) # entry.path is already a plain string!
                        
        return files


def getPlotter(sample,sampleDir,sampleType,eras,prod,analysis):
    plotters=[]
    for era in eras:
        #if we are doing data overrtide the search
        files=[]
        if sampleType=='DATA':
            if analysis in ['wen2g','zee2g','wegamma']:
                if era=='2018':
                    files = getFiles('EGamma',sampleDir,sampleType,era,prod)
                else:
                    files = getFiles('SingleElectron',sampleDir,sampleType,era,prod)
            else:
                files = getFiles('SingleMuon',sampleDir,sampleType,era,prod)

                
        else: #MC means query
            files = getFiles(sample,sampleDir,sampleType,era,prod)
        for f in files:
            plotters.append(rdf_plotter(f,isMC=(sampleType=='MC'), tree=analysis, report = "Report_" + analysis))
            #scale with the luminosity
            if sampleType=='MC':
                plotters[-1].addCorrectionFactor(lumifb[era], "flat")
                plotters[-1].addCorrectionFactor('1000', "flat") #to conevrt to pb-1
            
            #Deal with the HEM cuts
            if era == '2018' and analysis == 'wen2g' and sampleType=='DATA':
                plotters[-1].defaultCuts = "((run>=319077&&(Electron_eta[W_l1_idx]>-1.3||Electron_eta[W_l1_idx]<-3.0)&&(Electron_phi[W_l1_idx]>-0.87||Electron_phi[W_l1_idx]<-1.57))||(run<319077))"
            elif era=='2018' and  analysis == 'wen2g' and sampleType=='MC':   
                plotters[-1].addCorrectionFactor(str(21080.0/59830), "flat")
                plotters.append(rdf_plotter(f, True, tree = analysis, defaultCuts = cutsHEM))
                plotters[-1].addCorrectionFactor(lumifb[era], "flat")
                plotters[-1].addCorrectionFactor('1000', "flat") #to conevrt to pb-1                
                plotters[-1].addCorrectionFactor(str(38750./59830), "flat")
    p = merged_plotter(plotters)

    return p


def calculate_fake_rate(sampleDir,prod,eras=['2016','2017','2018'],ana='wmugamma',ptbins=[25.,30.,35.,40.,50.,60.,80.,150.],etabins=[-2.5,-2.0,-1.57,-1.44,-0.8,0.8,1.44,1.57,2.0,2.5],arrayName="fake_rate_VH",doMCClosure=False,outdir='.',file_extension='png'):
    def clopper_pearson(k, n, confidence=0.68):
        """
        Compute the Clopper-Pearson interval for a binomial proportion.
        
        k: number of successes (array-like)
        n: number of trials (array-like)
        confidence: confidence level (0.95 for 95%)
        """
        alpha = 1 - confidence
    
        # Lower bound: alpha/2 quantile of Beta(k, n-k+1)
        # We use np.maximum to handle the edge case where k=0
        lower = beta.ppf(alpha / 2, k, n - k + 1)
        lower = np.nan_to_num(lower, nan=0.0)
        
        # Upper bound: 1-alpha/2 quantile of Beta(k+1, n-k)
        # We use np.nan_to_num to handle the edge case where k=n
        upper = beta.ppf(1 - alpha / 2, k + 1, n - k)
        upper = np.nan_to_num(upper, nan=1.0)
    
        return lower, upper

    if doMCClosure:
        wjets=getPlotter('WJetsToLNu_HT',sampleDir,'MC',eras,prod,ana)
        vjets=getPlotter('DYJetsToLL_M50_LO',sampleDir,'MC',eras,prod,ana)
        tt=getPlotter('TTJets',sampleDir,'MC',eras,prod,ana)
        fr = merged_plotter([wjets,vjets,tt])
    else:    
        fr = getPlotter('nothing',sampleDir,'DATA',eras,prod,ana)

    if ana=='wmugamma':
        lepton='MU'
        cutsMisID='(Photon_pt>0)'#dummy
    else:    
        lepton='ELE'
        fr.define("Photon_EGMass","invMassEG(Electron_pt[W_l1_idx],Electron_eta[W_l1_idx],Electron_phi[W_l1_idx],Electron_mass[W_l1_idx],Photon_pt, Photon_eta,Photon_phi,Photon_mass)")
        cutsMisID='(abs(Photon_EGMass-91)>10.)'

    cuts_denom = "&&".join([cuts['W'][lepton],"(Sum(Photon_preselection)==1)",cutsMisID,'(Photon_preselection==1)'])
    cuts_num = '&&'.join([cuts['W'][lepton],"(Sum(Photon_preselection)==1)",cutsMisID,'((Photon_preselection*Photon_cutBased)>0)'])

    xedges,yedges,denominator,w2_denom = fr.array2d('Photon_pt','Photon_eta',cuts_denom,('denom','denom',len(ptbins)-1,np.array(ptbins),len(etabins)-1,np.array(etabins)))
    xedges,yedges,numerator,w2_num = fr.array2d('Photon_pt','Photon_eta',cuts_num,('num','num',len(ptbins)-1,np.array(ptbins),len(etabins)-1,np.array(etabins)))

    #plot the numerator and denominator separately
    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(xedges, yedges, denominator, cmap='plasma', edgecolors='white', linewidth=0.5)
    plt.colorbar(mesh, label='fake rate denominator')
    plt.savefig(f'{outdir}/{arrayName}_{lepton}_denominator.{file_extension}', dpi=400, bbox_inches='tight')
    plt.close()
    
    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(xedges, yedges, denominator, cmap='plasma', edgecolors='white', linewidth=0.5)
    plt.colorbar(mesh, label='fake rate numerator')
    plt.savefig(f'{outdir}/{arrayName}_{lepton}_numerator.{file_extension}', dpi=400, bbox_inches='tight')
    plt.close()
  
    
    result = np.zeros_like(numerator)
    fake_rate = np.divide(numerator,denominator,out=result,where=(denominator !=0))
        
    fake_rate_down,fake_rate_up = clopper_pearson(numerator,denominator)
    #make a fake rate plot               
    fig, ax = plt.subplots(figsize=(20,10))
    mesh = ax.pcolormesh(xedges, yedges, fake_rate, cmap='rainbow', edgecolors='white', linewidth=0.0)
    cbar=plt.colorbar(mesh)
    if ana=='wmugamma':
        cbar.set_label(r'${\mathcal{f}}_\gamma (W\rightarrow \mu\nu)$',size=20)
    else:
        cbar.set_label(r'${\mathcal{f}}_\gamma (W\rightarrow e\nu)$',size=20)
        
    # 4. Add labels and uncertainties (The "ROOT TEXT" part)
    # Calculate centers for text placement
    x_centers = (xedges[:-1] + xedges[1:]) / 2
    y_centers = (yedges[:-1] + yedges[1:]) / 2

    #put the values on the plot but also make a string
    st=f'std::vector<std::vector<std::vector<float> > > {arrayName}_{lepton}_vals = {{'
    for i, x_val in enumerate(x_centers):
        st=st+'{'
        for j, y_val in enumerate(y_centers):
            st=st+'{'+','.join([str(fake_rate[j, i]),str(fake_rate_up[j, i]),str(fake_rate_down[j, i])])+'}'
            if j!=(len(y_centers)-1):
                st=st+','
            val = fake_rate[j, i] # Note: pcolormesh uses (row, col) which is (y, x)
            errUp = fake_rate_up[j, i]-val
            errDwn = val-fake_rate_down[j, i]                    
            ax.set_xlabel(r"$\gamma p_{T}$ (GeV)",fontsize=20)
            ax.set_ylabel(r"$\gamma \eta$",fontsize=20)
            ax.tick_params(axis='both', which='major', labelsize=20)
            
        st=st+'}'
        if i!=(len(x_centers)-1):
            st=st+','
    st=st+'};'
    if doMCClosure:
        mh.cms.label(analysis_status,data=False, ax=ax, loc=0)
    else:
        mh.cms.label(analysis_status, data=True, lumi=(lumifb[eras[0]] if len(eras)==1 else lumifb['Run2']), ax=ax, loc=0)
        #print it in C++ format
        #note that we remove the last edge on how the code is defined to work
    st=st+'\n'+f'std::vector<float> {arrayName}_{lepton}_xbins = {{'+','.join([str(x) for x in xedges[:-1]])+'};\n'+f'std::vector<float> {arrayName}_{lepton}_ybins = {{'+','.join([str(y) for y in yedges[:-1]])+'};\n'
    plt.savefig(f'{outdir}/{arrayName}_{lepton}.{file_extension}', dpi=400, bbox_inches='tight')
    with open(f'{outdir}/{arrayName}_{lepton}.pickle', "wb") as file:
        pickle.dump(fig, file)
    
    return st


    



def getSignalPlotter(sampleDir,prod,eras,analysis,mass,lifetime,signals=['ZH','ggZH','WH','ttH'],brgg=0.5):
    xsecs = {'ggZ': "0.1057", #units in pb
             'Z': "0.8696",
             'Wplus': "0.84",
             'Wminus': "0.5328",
             'tt':"0.5071"} 

    #source https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV
    #asymetric
    xsecUnc = {'WH'  : [0.005, -0.007],
               'ZH'  : [0.038, -0.031],
               'ggZH': [0.251, -0.189],
               'ttH' : [0.058, -0.092]}

    #symmetric
    pdfUnc = {'WH'  : 1.019,
              'ZH'  : 1.016,
              'ggZH': 1.024,
              'ttH' : 1.036}


    #BRs from PDG
    BRs = {"Z": "(.03363+.03366+.033696)",  #e⁺e⁻, μ⁺μ⁻, τ⁺τ⁻
           "W": "(.1063+.1071+.1138)",      #μν, eν, τν
           "ttSemiLeptonic": "2*(1 - (.1110 + .1140 + .107))*(.1110 + .1140 + .107)" #tt̄->W⁺W⁻bb̄(100%) -> W->(hadronic)W->(leptonic)bb̄->(jets) * 2
           }

    V=[]
    decays=[]
    for sig in signals:
        if sig == 'ttH' or sig == 'ttH2G2Q' or sig == 'ttH4G':
            V.extend(['tt'])
        elif sig == 'WH' or sig == 'WH2G2Q' or sig == 'WH4G':
            V.extend(['Wplus','Wminus'])
        elif sig == 'ZH' or sig == 'ZH2G2Q' or sig == 'ZH4G':
            V.extend(['Z'])
        elif sig == 'ggZH' or sig == 'ggZH2G2Q' or sig == 'ggZH4G':
            V.extend(['ggZ'])
        if '2G2Q' in sig:
            decays.extend(['2G2Q'])
        elif '4G' in sig:
            decays.extend(['4G'])
        else :
            decays.extend(['4G','2G2Q'])
            
    decays=list(set(decays))        
    plotters=[]
    for era in eras:
        for sig in V:
            for br in decays:
                fs=getFiles(f"{sig}H{br}_M{mass}_ctau{lifetime}_{era}",sampleDir,"MC",era,prod)
                if len(fs)==0:                    
                    print(f"WARNING! NO FILE FOUND matching pattern: {sig}H{br}_M{mass}_ctau{lifetime}_{era}")
                for f in fs:
                    plotters.append(rdf_plotter(f, tree=analysis,isMC=True, report = "Report_" + analysis))
                    plotters[-1].addCorrectionFactor(lumifb[era], "flat")                
                    plotters[-1].addCorrectionFactor('1000', "flat") #to conevrt to pb-1
#                    print(f"{sig}H{br}_M{mass}_ctau{lifetime}_{era}")
                    plotters[-1].define('nPhiGamma', "nSpecificGenParticles(GenPart_pdgId,GenPart_genPartIdxMother,GenPart_status,22,9000006,1)") #To distinguis events based on the branching ratio
                    if br=='2G2Q':
                        plotters[-1].addCorrectionFactor('nPhiGamma==2', "flat")
                        w=brgg*(1-brgg)/(0.5*0.5)
                        weight = f"({w})"
                    else:
                        plotters[-1].addCorrectionFactor('nPhiGamma==4', "flat")
                        w=brgg*brgg
                        weight = f"({w})"
                        

                    weight +="*"+xsecs[sig]
                    if sig in ['Z','ggZ']:
                        weight+="*"+BRs['Z']
                    elif sig in ['Wplus','Wminus']:
                        weight+="*"+BRs['W']
                    elif sig in ['tt']:
                        weight+="*"+BRs['ttSemiLeptonic']
                    plotters[-1].addCorrectionFactor(weight, "flat")
                    
                    #Deal with the HEM cuts
                    if era=='2018' and  analysis == 'wen2g':   
                        plotters[-1].addCorrectionFactor(str(21080.0/59830), "flat")
                        plotters.append(rdf_plotter(f, isMC=True, tree = analysis, defaultCuts = cutsHEM))
                        plotters[-1].addCorrectionFactor(lumifb[era], "flat")
                        plotters[-1].addCorrectionFactor('1000', "flat") #to conevrt to pb-1                             
                        plotters[-1].addCorrectionFactor(str(38750./59830), "flat")
                        plotters[-1].addCorrectionFactor(weight, "flat")
                        plotters[-1].define('nPhiGamma', "nSpecificGenParticles(GenPart_pdgId,GenPart_genPartIdxMother,GenPart_status,22,9000006,1)") #To distinguis events based on the branching ratio
                        if br=='2G2Q':
                            plotters[-1].addCorrectionFactor('nPhiGamma==2', "flat")                
                        else:
                            plotters[-1].addCorrectionFactor('nPhiGamma==4', "flat")                                
                    
    p=None
    if len(plotters)>0:                
        p = merged_plotter(plotters)
    return p

def getAnalysis(sampleDir,prod,ana,era='Run2',masses=masses,lifetimes=lifetimes,signals=['ZH','ggZH','WH','ttH'],brphiphi=0.01,brgg=0.5,background_method="fakerate"):
    analysis={}
    if ana in ['wmn2g','zmm2g']:
        lepton='MU'
    else:    
        lepton='ELE'
    
    if era=='Run2':
        eras=['2016','2017','2018']
    else:
        eras=[era]
    analysis['data']=getPlotter('nothing',sampleDir,'DATA',eras,prod,ana)        
    #now create a background plotter with the sideband method we are talking about
    analysis['bkg']={}
    for m in masses:
        if background_method=="abcd":
            analysis['bkg'][m]=abcd_plotter(cuts[ana][m]['sr'],
                                                  cuts[ana][m]['cr_abcd'],
                                                  cuts[ana][m]['ssb'],
                                                  cuts[ana][m]['csb_abcd'],analysis['data'].plotters)
        elif background_method=="fakerate":
            if era!='Run2':
                st = f'fake_rate(Photon_pt[best_2g_idx1_m{m}],Photon_eta[best_2g_idx1_m{m}],Photon_pt[best_2g_idx2_m{m}],Photon_eta[best_2g_idx2_m{m}],(Photon_cutBased[best_2g_idx1_m{m}]>0),(Photon_cutBased[best_2g_idx2_m{m}]>0),fake_rate_{era}_{lepton}_vals,fake_rate_{era}_{lepton}_xbins,fake_rate_{era}_{lepton}_ybins)'
            else:#assume run2
                st = f'fake_rate(Photon_pt[best_2g_idx1_m{m}],Photon_eta[best_2g_idx1_m{m}],Photon_pt[best_2g_idx2_m{m}],Photon_eta[best_2g_idx2_m{m}],(Photon_cutBased[best_2g_idx1_m{m}]>0),(Photon_cutBased[best_2g_idx2_m{m}]>0),fake_rate_{lepton}_vals,fake_rate_{lepton}_xbins,fake_rate_{lepton}_ybins)'
                
            analysis['bkg'][m]=fakerate_plotter(cuts[ana][m]['presr'],
                                                  cuts[ana][m]['precr'],
                                                  analysis['data'].plotters,
                                                  'fakeRate',
                                                  st
                                                  )                                                  
            #define systematics
            analysis['bkg'][m].define("fakeRate_val","fakeRate[0]")           
            analysis['bkg'][m].define("fakeRate_up","fakeRate[1]")
            analysis['bkg'][m].define("fakeRate_down","fakeRate[2]")
        if background_method!="":
    
            analysis['bkg'][m].setFillProperties(1001, ROOT.kAzure+5)
            analysis['bkg'][m].setLineProperties(1, ROOT.kAzure+5, 3)        


    #For these MC sampleswe do not apply photon scale factors since they are not used for data MC/comparison
    #We do it for the signals though
    analysis['wjets']=getPlotter('WJetsToLNu_HT',sampleDir,'MC',eras,prod,ana)
    analysis['wjets'].addCorrectionFactor(leptonSF[ana],'flat')
    
    analysis['zjets']=getPlotter('DYJetsToLL_M50_LO',sampleDir,'MC',eras,prod,ana)
    analysis['zjets'].addCorrectionFactor(leptonSF[ana],'flat')
    if ana in ['wmn2g','wenu2g']:
        analysis['tt']=getPlotter('TTJets',sampleDir,'MC',eras,prod,ana)
        analysis['tt'].addCorrectionFactor(leptonSF[ana],'flat')
    else:
        analysis['tt']=getPlotter('TTJets_DiLept',sampleDir,'MC',eras,prod,ana)
        analysis['tt'].addCorrectionFactor(leptonSF[ana],'flat')
        
    analysis['tt'].setFillProperties(1001, ROOT.kAzure-2)
    analysis['tt'].setLineProperties(1, ROOT.kAzure-2, 3)

    analysis['wg']=getPlotter('WGToLNuG',sampleDir,'MC',eras,prod,ana)
    analysis['wg'].addCorrectionFactor(leptonSF[ana],'flat')
    
    analysis['wgg']=getPlotter('WGG',sampleDir,'MC',eras,prod,ana)
    analysis['wgg'].addCorrectionFactor(leptonSF[ana],'flat')

    
    analysis['signal']={}
    for m in masses:
        analysis['signal'][m]={}
        i=0
        for ct in lifetimes:
            analysis['signal'][m][ct]={}
            p=getSignalPlotter(sampleDir,prod,eras,ana,m,ct,signals,brgg)
            if p!=None:
                analysis['signal'][m][ct]['sum']=p
                analysis['signal'][m][ct]['sum'].addCorrectionFactor(leptonSF[ana],'flat')
                analysis['signal'][m][ct]['sum'].addCorrectionFactor(photonSF[m],'flat')            
                analysis['signal'][m][ct]['sum'].addCorrectionFactor(str(brphiphi),'flat')
            for signal in signals:
                p=getSignalPlotter(sampleDir,prod,eras,ana,m,ct,[signal],brgg)
                if p!=None:
                    analysis['signal'][m][ct][signal]=p
                    analysis['signal'][m][ct][signal].addCorrectionFactor(leptonSF[ana],'flat')
                    analysis['signal'][m][ct][signal].addCorrectionFactor(photonSF[m],'flat')            
                    analysis['signal'][m][ct][signal].addCorrectionFactor(str(brphiphi),'flat')
                
    return analysis
            




def getLimitData(outputDir='VHresults',prefix='asymptotics',lifetimes=lifetimes,brgg=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]):
        target_dir = f"{outputDir}/limits/"
        limits={'br':np.array([]),
                'm':np.array([]),
                'ctau':np.array([]),
                'limit':np.array([]),
                'quantile':np.array([])}
        
        for br in brgg:
            for ctau in lifetimes:
                rdf = ROOT.RDataFrame('limit',f"{target_dir}/{prefix}_ctau{ctau}_bgg{br}.root")
                npy=rdf.AsNumpy(columns=['mh','limit','quantileExpected'])
                npy['CL']=np.array([],dtype=int)
                npy['br']=np.full(len(npy['mh']),float(br))
                for p in npy['quantileExpected']:
                    if p<0:
                        npy['CL']=np.append(npy['CL'],-1)
                    elif p>0.024 and p<0.026:
                        npy['CL']=np.append(npy['CL'],25)
                    elif p>0.15 and p<0.17:
                       npy['CL']=np.append(npy['CL'],160)
 
                    elif p>0.49 and p<0.51:
                       npy['CL']=np.append(npy['CL'],500)
 
                    elif p>0.83 and p<0.85:
                        npy['CL']=np.append(npy['CL'],840)
 
                    elif p>0.974 and p<0.976:
                       npy['CL']=np.append(npy['CL'],975)
  
                        
                npy['ctau']=np.full(len(npy['mh']),float(ctau))
                limits['br']=np.concatenate([limits['br'],npy['br']])
                limits['m']=np.concatenate([limits['m'],npy['mh']])
                limits['ctau']=np.concatenate([limits['ctau'],npy['ctau']])
                limits['limit']=np.concatenate([limits['limit'],npy['limit']])
                limits['quantile']=np.concatenate([limits['quantile'],npy['CL']])
        df = pd.DataFrame(limits).sort_values(by=['m','ctau','br','quantile'])

        return df
                
                                  

        




def runAction(sampleDir,prod,action='fakerate_closure',masses=masses,outputDir='VHresults',era='Run2',analyses=analyses,signals=['ZH','ggZH','WH','ttH'],lifetimes=lifetimes,brphiphi=0.01,blinded=False,file_extension='png'):



    if era=='Run2':
        eras=['2016','2017','2018']
    else:
        eras=[era]


    #ACTION: Kinematic Fit plots
    if action=="kinfit_plots":
        for m in masses:
            ana='wmn2g'
            analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,signals=signals,lifetimes=lifetimes)
            stack=mplhep_plotter(com=center_of_mass[era],data=False,lumi=None,label=analysis_status)
            stack.stack=False
            signal_plotters=[]
            for ctau in lifetimes:
                analysis['signal'][m][ctau]['sum'].define("deltaLXY",f"best_2g_dxy_m{m}-genLxy(GenPart_vx[GenPart_isSignal], GenPart_vy[GenPart_isSignal])")
                signal_plotters.append(analysis['signal'][m][ctau]['sum'])
                stack.add_plotter(analysis['signal'][m][ctau]['sum'],label=r'$c\tau=$'+f"{ctau} mm",typeP='signal',error_mode='w2',color=signal_colors[ctau])                

            #make the resolution plot vs gen Lxy
#            merged=merged_plotter(signal_plotters)
#            xedges,yedges,delta2d,w2_2d = merged.array2d('genLXY','deltaLXY',cuts[ana][m]['sr'],('delta','delta',22,0.,110.,24,-8.,8.))
#            x_centers = (xedges[:-1] + xedges[1:]) / 2
#            y_centers = (yedges[:-1] + yedges[1:]) / 2
#            y_centers_sq = y_centers ** 2
#            numerator = np.sum(delta2d.T* y_centers_sq, axis=1)
            # Denominator: Total counts in each x bin#
#            denominator = np.sum(delta2d.T, axis=1)

            # Suppress division-by-zero warnings for empty x-bins
#            with np.errstate(divide='ignore', invalid='ignore'):
#                mean_y_sq = numerator / denominator
#                y_rms = np.sqrt(mean_y_sq)
#            fig, ax = plt.subplots()
#            mesh = ax.pcolormesh(xedges, yedges, delta2d, cmap='PuBu',norm=colors.LogNorm(),shading='flat')
#            plt.plot(x_centers,y_rms,color='k')
#            plt.savefig(f'{outputDir}/kinfit_delta2d.{file_extension}', dpi=400, bbox_inches='tight')
#            plt.close()
            
                
            dlxy,ax=stack.hist1d("deltaLXY",cuts[ana][m]['sr'],model=('a','a',30,-15,15),alpha=1,xlabel=r"$\Delta L_{xy}$ ",xunits="cm",show=False,legend_loc='upper left',signalstep=True,top_space=0.3)

            ax.text(0.05, 0.65, rf"$m_\Phi={m}\ \mathrm{{GeV}}$" , transform=ax.transAxes,
                     fontsize=30, ha='left', va='center', 
                     bbox=dict(facecolor='white',edgecolor='none', alpha=1))

            plt.savefig(f'{outputDir}/kinfit_deltaLXY_m{m}.{file_extension}', dpi=400, bbox_inches='tight')
            plt.close()

            with open(f'{outputDir}/kinfit_deltaLXY_m{m}.pickle', "wb") as file:
                pickle.dump(dlxy, file)
            mgg,ax=stack.hist1d(f"best_2g_raw_mass_m{m}",cuts[ana][m]['sr'],model=('a','a',100,8,m+2),alpha=1,xlabel=r"$m_{\gamma\gamma}$",xunits="GeV",show=False,legend_loc='upper left',signalstep=True,top_space=0.3)

            ax.text(0.05, 0.65, rf"$m_\Phi={m}\ \mathrm{{GeV}}$" , transform=ax.transAxes,
                    fontsize=30, ha='left', va='center', 
                    bbox=dict(facecolor='white',edgecolor='none', alpha=1))

            with open(f'{outputDir}/kinfit_mass_m{m}.pickle', "wb") as file:
                pickle.dump(mgg, file)            
            plt.savefig(f'{outputDir}/kinfit_mass_m{m}.{file_extension}', dpi=400, bbox_inches='tight')
            dxy,ax=stack.hist1d(f"best_2g_dxy_m{m}",cuts[ana][m]['sr'],model=('a','a',90,-10,80),alpha=1,xlabel=r"$d_{xy}$",xunits="cm",show=False,signalstep=True,top_space=0.3)
            ax.text(0.75, 0.65, rf"$m_\Phi={m}\ \mathrm{{GeV}}$" , transform=ax.transAxes,
                    fontsize=30, ha='left', va='center', 
                    bbox=dict(facecolor='white',edgecolor='none', alpha=1))
            
            plt.savefig(f'{outputDir}/kinfit_dxy{m}.{file_extension}', dpi=400, bbox_inches='tight')
            plt.close()
            
            with open(f'{outputDir}/kinfit_dxy{m}.pickle', "wb") as file:
                pickle.dump(dxy, file)            
        for ctau in lifetimes:
            ana='wmn2g'
            analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,signals=signals,lifetimes=lifetimes)
            stack=mplhep_plotter(com=center_of_mass[era],data=False,lumi=None,label=analysis_status)
            stack.stack=False
            
            for m in masses:
                analysis['signal'][m][ctau]['sum'].define("deltaLXY",f"best_2g_dxy_m{m}-genLxy(GenPart_vx[GenPart_isSignal], GenPart_vy[GenPart_isSignal])")
                stack.add_plotter(analysis['signal'][m][ctau]['sum'],label=rf"$m_{{\Phi}}= {m}\ \mathrm{{GeV}}$",typeP='signal',error_mode='w2',color=signal_colors_m[m])                

            dlxy,ax=stack.hist1d("deltaLXY",cuts[ana][m]['sr'],model=('a','a',30,-15,15),alpha=1,xlabel=r"$\Delta L_{xy}$ ",xunits="cm",show=False,legend_loc='upper left',signalstep=True,top_space=0.3)

            ax.text(0.05, 0.65, rf"$c\tau={ctau}\ \mathrm{{mm}}$" , transform=ax.transAxes,
                     fontsize=30, ha='left', va='center', 
                     bbox=dict(facecolor='white',edgecolor='none', alpha=1))

            plt.savefig(f'{outputDir}/kinfit_deltaLXY_ct{ctau}.{file_extension}', dpi=400, bbox_inches='tight')
            plt.close()

            with open(f'{outputDir}/kinfit_deltaLXY_ct{ctau}.pickle', "wb") as file:
                pickle.dump(dlxy, file)
            mgg,ax=stack.hist1d(f"best_2g_raw_mass_m{m}",cuts[ana][m]['sr'],model=('a','a',100,8,m+2),alpha=1,xlabel=r"$m_{\gamma\gamma}$",xunits="GeV",show=False,legend_loc='upper left',signalstep=True,top_space=0.3)

            ax.text(0.05, 0.65, rf"$c\tau={ctau}\ \mathrm{{mm}}$" , transform=ax.transAxes,
                     fontsize=30, ha='left', va='center', 
                     bbox=dict(facecolor='white',edgecolor='none', alpha=1))

            with open(f'{outputDir}/kinfit_mass_ct{ctau}.pickle', "wb") as file:
                pickle.dump(mgg, file)            
            plt.savefig(f'{outputDir}/kinfit_mass_ct{ctau}.{file_extension}', dpi=400, bbox_inches='tight')
            dxy,ax=stack.hist1d(f"best_2g_dxy_m{m}",cuts[ana][m]['sr'],model=('a','a',90,-10,80),alpha=1,xlabel=r"$d_{xy}$",xunits="cm",show=False,signalstep=True,top_space=0.3)
            ax.text(0.75, 0.65, rf"$c\tau={ctau}\ \mathrm{{mm}}$" , transform=ax.transAxes,
                    fontsize=30, ha='left', va='center', 
                    bbox=dict(facecolor='white',edgecolor='none', alpha=1))
            
            plt.savefig(f'{outputDir}/kinfit_dxy_ct{ctau}.{file_extension}', dpi=400, bbox_inches='tight')
            plt.close()
            
            with open(f'{outputDir}/kinfit_dxy_ct{ctau}.pickle', "wb") as file:
                pickle.dump(dxy, file)            

    #ID vs LXY            
    if action=="id_vs_lxy":
        print("To run this you need to merge all W signals together to a file and pint the code to it.")
        rdf=ROOT.RDataFrame('wmn2g','/tank/ddp/DDP/allW_signals.root')
        rdf=rdf.Define("genLXY","genLxy(GenPart_vx[GenPart_isSignal], GenPart_vy[GenPart_isSignal])")
        rdf=rdf.Filter(cuts['wmn2g'][20]['precr'])
        denom=rdf.Histo1D(("denom","denom",20,0.,100.),"genLXY")
        rdf2=rdf.Filter(cuts['wmn2g'][20]['presr'])        
        num=rdf2.Histo1D(("denom","denom",20,0.,100.),"genLXY")
        g=ROOT.TGraphAsymmErrors()
        g.Divide(num.GetValue(),denom.GetValue())
        f=ROOT.TFile("effvslxy.root","RECREATE")
        f.cd()
        g.Write("eff")
        denom.Write("denom")
        num.Write("num")
        f.Close()
        


                
    #ACTION: Electron mis-id  plots
    if action=="electron_misID_plots":
        for m in masses:
            ana='wen2g'
            analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,signals=['WH'],lifetimes=lifetimes)
            stack=mplhep_plotter(com=center_of_mass[era],data=False,lumi=None)
            stack.stack=False

            for ctau in lifetimes:
                stack.add_plotter(analysis['signal'][m][ctau]['sum'],label=r'$c\tau=$'+f"{ctau} mm",typeP='signal',error_mode='w2',color=signal_colors[ctau])                
            stack.add_plotter(analysis['zjets'],label='Z+jets',typeP='background',error_mode='w2')
            stack.add_plotter(analysis['wjets'],label='W+jets',typeP='background',error_mode='w2')
            
            stack.add_plotter(analysis['tt'],label=r'$t\bar{t}$+jets',typeP='background',error_mode='w2')
            cutsSpecial='&&'.join([cuts['W']['ELE'],
                            cuts['pt'][m],
                            cuts['photons'][m],
#                            cuts['misID']['W']['ELE'][m],
                            f"((Photon_passCutBasedID[best_2g_idx1_m{m}]+Photon_passCutBasedID[best_2g_idx2_m{m}])==2)"])
            
            fig,ax=stack.hist1d(f"best_2g_misID1_m{m}",cutsSpecial,model=('a','a',50,0,200),alpha=0.25,xlabel=r"$m_{\text{tag1}}$ ",xunits="GeV",show=False,legend_loc='upper left')
            ax.text(0.05, 0.55, rf"$m_\Phi={m}\ \mathrm{{GeV}}$" , transform=ax.transAxes,
                    fontsize=20, ha='left', va='center', 
                    bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
            
            plt.savefig(f'{outputDir}/electron_misID_tag1_m{m}.{file_extension}', dpi=400, bbox_inches='tight')
            with open(f'{outputDir}/electron_misID_tag1_m{m}.pickle', "wb") as file:
                pickle.dump(fig, file)            

            
            fig,ax=stack.hist1d(f"best_2g_misID2_m{m}",cutsSpecial,model=('a','a',50,0,200),alpha=0.25,xlabel=r"$m_{\text{tag2}}$ ",xunits="GeV",show=False,legend_loc='upper left')
            ax.text(0.05, 0.55, rf"$m_\Phi={m}\ \mathrm{{GeV}}$" , transform=ax.transAxes,
                    fontsize=20, ha='left', va='center', 
                    bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
            
            plt.savefig(f'{outputDir}/electron_misID_tag2_m{m}.{file_extension}', dpi=400, bbox_inches='tight')
            with open(f'{outputDir}/electron_misID_tag2_m{m}.pickle', "wb") as file:
                pickle.dump(fig, file)            

            fig,ax=stack.hist1d(f"best_2g_misID3_m{m}",cutsSpecial,model=('a','a',50,0,200),alpha=0.25,xlabel=r"$m_{\text{tag3}}$ ",xunits="GeV",show=False,legend_loc='upper left')
            ax.text(0.05, 0.55, rf"$m_\Phi={m}\ \mathrm{{GeV}}$" , transform=ax.transAxes,
                    fontsize=20, ha='left', va='center', 
                    bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
            
            plt.savefig(f'{outputDir}/electron_misID_tag3_m{m}.{file_extension}', dpi=400, bbox_inches='tight')               
            with open(f'{outputDir}/electron_misID_tag3_m{m}.pickle', "wb") as file:
                pickle.dump(fig, file)            

        
    #ACTION: Low pt Photon Background
    if action=="lowpt_photon_background":
        for m in masses:
            ana='wmn2g'
            analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,signals=['WH'],lifetimes=lifetimes)
            stack=mplhep_plotter(com=center_of_mass[era],data=False,lumi=None)
            stack.stack=False

            for ctau in lifetimes:
                stack.add_plotter(analysis['signal'][m][ctau]['sum'],label=r'$c\tau=$'+f"{ctau} mm",typeP='signal',error_mode='w2',color=signal_colors[ctau])                
            stack.add_plotter(analysis['wjets'],label='W+jets',typeP='background',error_mode='w2')
            stack.add_plotter(analysis['zjets'],label='Z+jets',typeP='background',error_mode='w2')
            stack.add_plotter(analysis['tt'],label=r'$t\bar{t}$+jets',typeP='background',error_mode='w2')
            cutsSpecial='&&'.join([cuts['W']['MU'],
#                            cuts['pt'][m],
                            cuts['photons'][m],
                            cuts['misID']['W']['MU'][m],
                            f"((Photon_passCutBasedID[best_2g_idx1_m{m}]+Photon_passCutBasedID[best_2g_idx2_m{m}])==2)"])
            stack.define("pt1",f"Photon_pt[best_2g_idx1_m{m}]")
            stack.define("pt2",f"Photon_pt[best_2g_idx2_m{m}]")
            
            fig,ax=stack.hist1d("pt1",cutsSpecial,model=('a','a',20,20,100),alpha=0.5,xlabel=r"$pt_{\gamma_1}$ ",xunits="GeV",show=False,legend_loc='upper right',signalstep=True)
            ax.text(0.75, 0.55, rf"$m_\Phi={m}\ \mathrm{{GeV}}$" , transform=ax.transAxes,
                    fontsize=20, ha='left', va='center', 
                    bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))

            plt.savefig(f'{outputDir}/lowpt_photon_background1_m{m}.{file_extension}', dpi=400, bbox_inches='tight')
            with open(f'{outputDir}/lowpt_photon_background1_m{m}.pickle', "wb") as file:
                pickle.dump(fig, file)            
            
            fig,ax=stack.hist1d("pt2",cutsSpecial,model=('a','a',20,20,100),alpha=0.5,xlabel=r"$pt_{\gamma_2}$ ",xunits="GeV",show=False,legend_loc='upper right',signalstep=True)
            ax.text(0.75, 0.55, rf"$m_\Phi={m}\ \mathrm{{GeV}}$" , transform=ax.transAxes,
                    fontsize=20, ha='left', va='center', 
                    bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
            
            plt.savefig(f'{outputDir}/lowpt_photon_background2_m{m}.{file_extension}', dpi=400, bbox_inches='tight')
            with open(f'{outputDir}/lowpt_photon_background2_m{m}.pickle', "wb") as file:
                pickle.dump(fig, file)            

        

        
    #ACTION: Calculate fake rates
    elif action=="fakerate_calc":
        print("Running Calculation of Fake rates")
        with open('common/vhFakeRates.h', "w") as file:
            file.write("#ifndef FAKERATES\n")
            file.write("#define FAKERATES\n")
            for ana in ['wmugamma','wegamma']:
                for e in eras:
                    fr=calculate_fake_rate(sampleDir,prod,[e],ana=ana,arrayName=f"fake_rate_{e}",outdir=outputDir,file_extension=file_extension)
                    file.write(fr)
                #write the average all run fake rate for studies
                fr=calculate_fake_rate(sampleDir,prod,eras,ana=ana,arrayName="fake_rate",outdir=outputDir,file_extension=file_extension)
                file.write(fr)
                #write the average all run MC fake rate for studies
                fr=calculate_fake_rate(sampleDir,prod,eras,ana=ana,arrayName="fake_rate_MC",outdir=outputDir,doMCClosure=True,file_extension=file_extension)
                file.write(fr)
            file.write("#endif\n")
    


    #ACTION: fake rate MC Closure
    elif action=="fakerate_closure":
        print("Running MC Closure of Fake rates")
        myTexts = {
            'wmn2g':r'$W\rightarrow \mu \nu$',
            'wen2g':r'$W\rightarrow e \nu$',
            'zmm2g':r'$Z\rightarrow \mu \mu$',
            'zee2g':r'$Z\rightarrow e e$',

            }
        mh.style.use('CMS')

        
        for m in masses:    
            #create the common axes            
            fig,ax = plt.subplots(nrows=2,ncols=len(analyses),sharex="col",sharey="row",figsize=(20,12),gridspec_kw={'height_ratios':[5,1],
                                                                                                                    'width_ratios':list(np.full(len(analyses),1))})
            plt.subplots_adjust(wspace=0,hspace=0)        
        
            for i,ana in enumerate(analyses):

                if ana in ['wmn2g','zmm2g']:
                    lepton='MU'
                else:
                    lepton='ELE'
                #create a plotter that has all MC as data
                analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,signals=[],lifetimes=[])
                plotters=[analysis['wjets'],analysis['zjets'],analysis['tt']]
                #create a fake rate MC plotter
                print(f"Running {ana} m={m} GeV") 
                bkg=fakerate_plotter(cuts[ana][m]['sr'],
                                     cuts[ana][m]['cr'],
                                     plotters,
                                     'fakeRate_MC',
                                     f'fake_rate(Photon_pt[best_2g_idx1_m{m}],Photon_eta[best_2g_idx1_m{m}],Photon_pt[best_2g_idx2_m{m}],Photon_eta[best_2g_idx2_m{m}],(Photon_cutBased[best_2g_idx1_m{m}]>0),(Photon_cutBased[best_2g_idx2_m{m}]>0),fake_rate_MC_{lepton}_vals,fake_rate_MC_{lepton}_xbins,fake_rate_MC_{lepton}_ybins)')                                                  
                #make a stack plotter and plot the stack
                stack=mplhep_plotter(com=center_of_mass[era],data=False,lumi=None)
                stack.add_plotter(bkg,label='Estimated by method',typeP='data',error_mode='poisson_bootstrap')               
                stack.add_plotter(analysis['wjets'],label='W+jets',typeP='background',error_mode='w2')
                stack.add_plotter(analysis['zjets'],label='DY+jets',typeP='background',error_mode='w2')
                stack.add_plotter(analysis['tt'],label=r'$t\bar{t}$+jets',typeP='background',error_mode='w2')
                stack.hist1d(f"best_2g_dxy_m{m}",cuts[ana][m]['sr'],('a','a',len(binning1d[ana])-1,binning1d[ana]),alpha=1.0,xlabel=r"$L_{xy}$",xunits="cm",legend_loc='upper right',show=False,ax=ax[0,i],dndx=False)
                stack.pull1d(f"best_2g_dxy_m{m}",cuts[ana][m]['sr'],('a','a',len(binning1d[ana])-1,binning1d[ana]),alpha=1.0,xlabel=r"$L_{xy}$",xunits="cm",legend_loc='upper right',show=False,ax=ax[1,i])
                
                analysis=None
                ax[0,i].text(0.5, 0.95, myTexts[ana] , transform=ax[0,i].transAxes,
                             fontsize=35, ha='center', va='center', 
                             bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
                ax[0,i].margins(x=0,y=0)
                ax[1,i].margins(x=0,y=0)
                
                if i==(len(analyses)-1):
                    ax[0,-1].legend(loc='upper left', bbox_to_anchor=(0.1,0.9),fontsize=35)
                    
                ax[0,i].tick_params(axis='both', which='major', labelsize=27.5)
                ax[1,i].tick_params(axis='both', which='major', labelsize=27.5)
                    #then stack backgrounds and then draw band
            lo, hi = ax[0,0].get_ylim()
            ax[0,0].set_ylim(lo, hi + (hi - lo) * 1.5)
            ax[1,0].set_ylim(-4,4)
            ax[1,0].set_yticks([-2,0,2])
            
            ax[0,0].text(0.45, 0.8, rf"$m_\Phi={m}\ \mathrm{{GeV}}$" , transform=ax[0,0].transAxes,
                          fontsize=35, ha='left', va='center', 
                          bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))

                    #Labels
            mh.cms.label(analysis_status,data=False,rlabel="", ax=ax[0,0], loc=0,fontsize=45)
            mh.cms.label(None,exp='',data=False,llabel="", ax=ax[0,-1], loc=0,com=center_of_mass[era],fontsize=35)
            ax[1,-1].set_xlabel(r"$L_{xy}$ (cm)",fontsize=35)
            ax[0,0].set_ylabel("Events/bin",fontsize=35)
            ax[1,0].set_ylabel("Pull",fontsize=35)
            fig.align_ylabels([ax[0,0], ax[1,0]])
            plt.savefig(f'{outputDir}/fakerate_closure_{m}_{era}.{file_extension}', dpi=400, bbox_inches='tight')
            with open(f'{outputDir}/fakerate_closure_{m}_{era}.pickle', "wb") as file:
                pickle.dump(fig, file)




        
#        print("Running MC Closure of Fake rates")
#            
#        for ana in analyses:
#            if ana in ['wmn2g','zmm2g']:
#                lepton='MU'
#            else:
#                lepton='ELE'
            #create a plotter that has all MC as data
#            analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,signals=[],lifetimes=[])
#            plotters=[analysis['wjets'],analysis['zjets'],analysis['tt']]
#            #create a fake rate MC plotter
#            for m in masses:
#                print(f"Running {ana} m={m} GeV") 
#                bkg=fakerate_plotter(cuts[ana][m]['sr'],
#                                            cuts[ana][m]['cr'],
#                                            plotters,
#                                            'fakeRate_MC',
#                                            f'fake_rate(Photon_pt[best_2g_idx1_m{m}],Photon_eta[best_2g_idx1_m{m}],Photon_pt[best_2g_idx2_m{m}],Photon_eta[best_2g_idx2_m{m}],(Photon_cutBased[best_2g_idx1_m{m}]>0)#,(Photon_cutBased[best_2g_idx2_m{m}]>0),fake_rate_MC_{lepton}_vals,fake_rate_MC_{lepton}_xbins,fake_rate_MC_{lepton}_ybins)')                                                  
                #make a stack plotter and plot the stack
#                stack=mplhep_plotter(com=center_of_mass[era],data=False,lumi=None)
#                stack.add_plotter(bkg,label='Fake rate estimate',typeP='data',error_mode='poisson_bootstrap')               
#                stack.add_plotter(analysis['wjets'],label='W+jets',typeP='background',error_mode='w2')
#                stack.add_plotter(analysis['zjets'],label='DY+jets',typeP='background',error_mode='w2')
#                stack.add_plotter(analysis['tt'],label=r'$t\bar{t}$+jets',typeP='background',error_mode='w2')
                #draw a plot
#                stack.unrolledCustom(f"best_2g_raw_mass_m{m}",f"best_2g_dxy_m{m}",cuts[ana][m]['sr'],binning[ana][m],alpha=1.0,xlabel=r"$d_{xy}$",xunits="cm",show=False,ylabel=r'$m_{\gamma\gamma}$',yunits='GeV',textx=0.7)
#                plt.savefig(f'{outputDir}/fakerate_closure_{ana}_{m}.{file_extension}', dpi=400, bbox_inches='tight')
#                stack=None
#                fr_plotter=None
#            analysis=None
#            plotters=None
                


    #ACTION: ABCD MC Closure               
    elif action=="abcd_closure":
        print("Running MC Closure of ABCD Method")
        
        for ana in analyses:
            #create a plotter that has all MC as data
            analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,signals=[],lifetimes=[])
            if ana in ['wmn2g','wen2g']:
                plotters=[analysis['wjets'],analysis['zjets'],analysis['tt']]
            else:
                plotters=[analysis['zjets'],analysis['tt']]                
            for m in masses:
                print(f"Running {ana} m={m} GeV")
                #create a ABCD MC plotter
                bkg=abcd_plotter(cuts[ana][m]['sr'],
                                 cuts[ana][m]['cr_abcd'],
                                 cuts[ana][m]['ssb'],
                                 cuts[ana][m]['csb_abcd'],plotters)
                #make a stack plotter and plot the stack
                stack=mplhep_plotter(com=center_of_mass[era],data=False,lumi=None)
                stack.add_plotter(bkg,label='ABCD estimate',typeP='data',error_mode='poisson_bootstrap')              #                stack.add_plotter(analysis['data'],label='Data',typeP='data',error_mode='poisson',color='red')               
                stack.add_plotter(analysis['wjets'],label='W+jets',typeP='background',error_mode='w2')
                stack.add_plotter(analysis['zjets'],label='DY+jets',typeP='background',error_mode='w2')
                stack.add_plotter(analysis['tt'],label=r'$t\bar{t}$+jets',typeP='background',error_mode='w2')
                #draw a plot
                stack.unrolledCustom(f"best_2g_raw_mass_m{m}",f"best_2g_dxy_m{m}",cuts[ana][m]['sr'],binning[ana][m],alpha=1.0,xlabel=r"$d_{xy}$",xunits="cm",show=False)
                plt.savefig(f'{outputDir}/abcd_closure_{ana}_{m}.{file_extension}', dpi=400, bbox_inches='tight')
                stack=None
                fr_plotter=None
            analysis=None
            plotters=None
                
    #ACTION: data_vs_background
    elif action=="data_vs_background":
        print("Make data vs background plots in control regions")
        for ana in analyses:
            analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,brphiphi=brphiphi,signals=signals,lifetimes=[100])
            for m in masses:
                stack=mplhep_plotter(label=analysis_status,data=True,lumi=lumifb[era],com=center_of_mass[era])
                stack.add_plotter(analysis['bkg'][m],label="Background",typeP='background',error_mode='w2')
                stack.add_plotter(analysis['data'],label="Data",typeP='data',error_mode='poisson')


                if ana in ['wmn2g','wen2g']:
                    stack.hist1d('W_mt',cuts[ana][m]['presr']+f"*(best_2g_raw_mass_m{m}>65)",("a","a",20,0,150),xlabel=r"$M_T$",xunits="GeV",show=False)
                    plt.savefig(f'{outputDir}/data_vs_bkg_MT_{ana}_{m}_{era}.{file_extension}', dpi=400, bbox_inches='tight')
#                    stack.hist1d(f'best_2g_dxy_m{m}',cuts[ana][m]['presr'],("a","a",25,-160,-10),xlabel=r"$L_{xy}$",xunits="cm",show=False)
#                    plt.savefig(f'{outputDir}/data_vs_bkg_lxy_{ana}_{m}_{era}.{file_extension}', dpi=400, bbox_inches='tight')
                    
                else:    
                    stack.hist1d('Z_mass',cuts[ana][20]['presr']+f"*(best_2g_raw_mass_m{m}>65)",("a","a",20,0,120),xlabel=r"$m_{\ell\ell}$",xunits="GeV",show=False)
                    plt.savefig(f'{outputDir}/data_vs_bkg_MLL_{ana}_{m}_{era}.{file_extension}', dpi=400, bbox_inches='tight')
#                    stack.hist1d(f'best_2g_dxy_m{m}',cuts[ana][m]['presr'],("a","a",5,-160,-10),xlabel=r"$L_{xy}$",xunits="cm",show=False)
#                    plt.savefig(f'{outputDir}/data_vs_bkg_lxy_{ana}_{m}_{era}.{file_extension}', dpi=400, bbox_inches='tight')
            
                
                print(f"Running {ana} m={m} GeV")
                #draw a plot


    #ACTION: signal_contamination
    elif action=="signal_contamination":
        mySignals={
            'wmn2g':['WH','ttH'],
            'wen2g':['WH','ttH'],
            'zmm2g':['ZH','ggZH'],
            'zee2g':['ZH','ggZH']}
        myTexts = {
            'wmn2g':r'$W\rightarrow \mu \nu$',
            'wen2g':r'$W\rightarrow e \nu$',
            'zmm2g':r'$Z\rightarrow \mu \mu$',
            'zee2g':r'$Z\rightarrow e e$',

            }
        mh.style.use('CMS')

        
        for m in masses:    
            #create the common axes            
            fig,ax = plt.subplots(nrows=1,ncols=len(analyses),sharex="col",sharey="row",figsize=(25,10))
            plt.subplots_adjust(wspace=0)        
        
            for i,ana in enumerate(analyses):
                #create a plotter that has all MC as data
                analysis=getAnalysis(sampleDir,prod,ana,background_method='',era=era,brphiphi=brphiphi,signals=mySignals[ana],lifetimes=lifetimes,brgg=0.5)
                print(f"Running {i} {ana} m={m} GeV")
                stack=mplhep_plotter(label=analysis_status,data=True,lumi=lumifb[era],com=center_of_mass[era],capsize=0)
                stack.add_plotter(analysis['data'],label="Data",typeP='data',error_mode='poisson')               
                for ctau in lifetimes:
                    stack.add_plotter(analysis['signal'][m][ctau]['sum'],label=r" $H \rightarrow\Phi\Phi,  c\tau =$ "+f"{ctau} mm",typeP='signal',error_mode='w2',color=signal_colors[ctau])
                #draw a plot
                stack.hist1d(f"best_2g_dxy_m{m}",cuts[ana][m]['cr'],('a','a',len(binning1d[ana])-1,binning1d[ana]),alpha=1.0,xlabel=r"$L_{xy}$",xunits="cm",legend_loc='upper right',show=False,ax=ax[i],dndx=False)
                
                analysis=None
                ax[i].text(0.5, 0.95, myTexts[ana] , transform=ax[i].transAxes,
                           fontsize=20, ha='center', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
                ax[i].margins(x=0,y=0)
                if i==(len(analyses)-1):
                    ax[i].legend(loc='upper left', bbox_to_anchor=(0.01,0.9))
                ax[i].tick_params(axis='both', which='major', labelsize=18)
            #then stack backgrounds and then draw band
            lo, hi = ax[0].get_ylim()
            ax[0].set_ylim(lo, hi + (hi - lo) * 0.25)

            ax[-1].text(0.1, 0.55, rf"$m_\Phi={m}\ \mathrm{{GeV}}$" , transform=ax[-1].transAxes,
                           fontsize=20, ha='left', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))

            ax[-1].text(0.1, 0.5, r"$\mathcal{B} (H \rightarrow \Phi\Phi)=0.01$" , transform=ax[-1].transAxes,
                           fontsize=20, ha='left', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
            ax[-1].text(0.1, 0.45, r"$\mathcal{B}(\Phi \rightarrow \gamma\gamma)=0.5$" , transform=ax[-1].transAxes,
                           fontsize=20, ha='left', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
            ax[-1].text(0.1, 0.4, "Bkg. extrapolation region" , transform=ax[-1].transAxes,
                           fontsize=20, ha='left', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
            
            #Labels
            mh.cms.label(analysis_status,data=True,rlabel="", ax=ax[0], loc=0)
            mh.cms.label(None,exp='',data=True,llabel="", ax=ax[-1], loc=0,lumi=lumifb[era],com=center_of_mass[era])
            ax[-1].set_xlabel(r"$L_{xy}$ (cm)",fontsize=20)
            ax[0].set_ylabel("Events",fontsize=20)

            plt.savefig(f'{outputDir}/signalContam_1D_{m}_{era}.{file_extension}', dpi=400, bbox_inches='tight')
            with open(f'{outputDir}/signalContam_1D_{m}_{era}.pickle', "wb") as file:
                pickle.dump(fig, file)




    #ACTION: step by step              
    elif action=="step_by_step":
        mySignals={
            'wmn2g':['WH','ttH'],
            'wen2g':['WH','ttH'],
            'zmm2g':['ZH','ggZH'],
            'zee2g':['ZH','ggZH']}
        
        for ana in analyses:
            if ana=='wmn2g':
                v='W'
                l='MU'
            elif ana=='wen2g':
                v='W'
                l='ELE'
            elif ana=='zmm2g':
                v='Z'
                l='MU'
            elif ana=='zee2g':
                v='Z'
                l='ELE'
            cutDescriptors=['Preselection & HLT','W/Z reco',r'$\gamma \ p_{T}$ cuts',r'$p_{T}>20$ GeV, m>4 GeV',r'$e\gamma$ mis-ID','final photon ID',r'$L_{xy}>-10$ cm',]
            colors=['dimgrey','lightcoral','chocolate','yellowgreen','turquoise','deepskyblue','darkviolet','magenta']
            analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,brphiphi=brphiphi,signals=mySignals[ana],lifetimes=lifetimes)
            mh.style.use('CMS')            
            fig,ax = plt.subplots(1,len(lifetimes),sharey=True,figsize=(10*len(lifetimes),10))
            plt.subplots_adjust(wspace=0)        
            for N,ctau in enumerate(lifetimes):
                for i,cutDesc in enumerate(cutDescriptors):
                    efficiency=[]
                    for m in masses:
                        myCuts=['1',cuts[v][l],cuts['pt'][m],cuts['photons'][m],cuts['misID'][v][l][m],f"((Photon_passCutBasedID[best_2g_idx1_m{m}]+Photon_passCutBasedID[best_2g_idx2_m{m}])==2)",f'(best_2g_dxy_m{m}>-10)']
                        analysis['signal'][m][ctau]['sum'].define("dummyDCVar",'1.0')
                        edges,data,w2=analysis['signal'][m][ctau]['sum'].array1d('dummyDCVar',"*".join(myCuts[0:i+1]),('a','a',2,0,2),error_mode='w2')
                        rate=float(np.sum(data))
                        error=np.sqrt(w2)
                        lower=np.sum(error[0,:])
                        upper=np.sum(error[1,:])
                        err=0.5*(lower+upper)
                        efficiency.append([rate,err])
                    ax[N].errorbar(masses,[s[0] for s in efficiency],[s[1] for s in efficiency],label=cutDesc,fmt='-o',color=colors[i])
                    ax[N].text(0.1, 0.95, r'$c\tau=$'+f'{ctau} mm' , transform=ax[N].transAxes,
                       fontsize=20, ha='center', va='center', 
                       bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))

            mh.cms.label(rlabel="", ax=ax[0], loc=0)
            mh.cms.label(None,exp='',llabel="", ax=ax[-1], loc=0,lumi=None,com=None,rlabel='13 TeV')
            ax[-1].set_xlabel(r"$m_{\Phi}$ (GeV)")
            ax[0].set_ylabel(r"Events for $\mathcal{BR}(H \rightarrow \Phi\Phi)$=0.01")
            ax[-1].legend(loc='upper right')
            limits=ax[0].get_ylim()
            plt.savefig(f'{outputDir}/step_by_step_{ana}_{era}_signal.{file_extension}', dpi=400, bbox_inches='tight')
                
            #Now background
            samp=['W+jets','DY+jets','tt+Jets']
            samps = ['wjets','zjets','tt']
            if ana in ['zee2g','zmm2g']:
                samp=['DY+jets','tt+Jets']
                samps=['zjets','tt']
            fig,ax = plt.subplots(1,len(samp),sharey=True,figsize=(10*len(samp),10))
            plt.subplots_adjust(wspace=0)

            for N,bkg in enumerate(samps):
                for i,cutDesc in enumerate(cutDescriptors):
                    efficiency=[]
                    analysis[bkg].define("dummyDCVar",'1.0')                    
                    for m in masses:
                        myCuts=['1',cuts[v][l],cuts['pt'][m],cuts['photons'][m],cuts['misID'][v][l][m],f"((Photon_passCutBasedID[best_2g_idx1_m{m}]+Photon_passCutBasedID[best_2g_idx2_m{m}])==2)",f'(best_2g_dxy_m{m}>-10)']
                        edges,data,w2=analysis[bkg].array1d('dummyDCVar',"*".join(myCuts[0:i+1]),('a','a',2,0,2),error_mode='w2')
                        rate=float(np.sum(data))
                        error=np.sqrt(w2)
                        lower=np.sum(error[0,:])
                        upper=np.sum(error[1,:])
                        err=0.5*(lower+upper)
                        efficiency.append([rate,err])
                    ax[N].errorbar(masses,[s[0] for s in efficiency],[s[1] for s in efficiency],label=cutDesc,fmt='-o',color=colors[i])
                    ax[N].text(0.1, 0.95, samp[N] , transform=ax[N].transAxes,
                       fontsize=20, ha='center', va='center', 
                       bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))

            mh.cms.label(rlabel="", ax=ax[0], loc=0)
            mh.cms.label(None,exp='',llabel="", ax=ax[-1], loc=0,lumi=None,com=None,rlabel='13 TeV')
            ax[-1].set_xlabel(r"$m_{\Phi}$ hypothesis (GeV)")
            ax[0].set_ylabel(r"Events")
            ax[-1].legend(loc='upper right')
#            ax[0].legend(loc='upper left')
#            limits=ax[0].get_ylim()
            ax[0].set_ylim(1e-1,1e+9)
            if ana in ['zee2g','zmm2g']:
                ax[0].set_ylim(1e-3,1e+7)
            
            ax[0].set_yscale('log')
            plt.savefig(f'{outputDir}/step_by_step_{ana}_{era}_bkg.{file_extension}', dpi=400, bbox_inches='tight')


                    
                    
                
        


            
    #ACTION: Final Plots              
    elif action=="final_plots":
        print("Make final Plots")
        mySignals={
            'wmn2g':['WH','ttH'],
            'wen2g':['WH','ttH'],
            'zmm2g':['ZH','ggZH'],
            'zee2g':['ZH','ggZH']}
            
        
        for ana in analyses:
            #create a plotter that has all MC as data
            analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,brphiphi=brphiphi,signals=mySignals[ana],lifetimes=lifetimes)
            for m in masses:
                print(f"Running {ana} m={m} GeV")
                stack=mplhep_plotter(label=analysis_status,data=True,lumi=lumifb[era],com=center_of_mass[era])
                stack.add_plotter(analysis['bkg'][m],label='Background',typeP='background',error_mode='poisson_bootstrap')
                for ctau in lifetimes:
                    stack.add_plotter(analysis['signal'][m][ctau]['sum'],label=r'$m_{\phi}$='+f"{m} GeV,"+r" $c\tau =$ "+f"{ctau} mm",typeP='signal',error_mode='w2',color=signal_colors[ctau])
                if blinded==False:
                    stack.add_plotter(analysis['data'],label="Data",typeP='data',error_mode='poisson')               
                #draw a plot
                stack.unrolledCustom(f"best_2g_raw_mass_m{m}",f"best_2g_dxy_m{m}",cuts[ana][m]['sr'],binning[ana][m],alpha=1.0,xlabel=r"$d_{xy}$",xunits="cm",ylabel=r'$m_{\gamma\gamma}$',yunits='GeV',textx=0.5,texty=0.95,show=False,legend_ax=0,legend_loc='upper left',legend_bbox=(0.01,0.9))
                if blinded:
                    plt.savefig(f'{outputDir}/blinded_prefit_{ana}_{m}_{era}.{file_extension}', dpi=400, bbox_inches='tight')
                else:
                    plt.savefig(f'{outputDir}/prefit_{ana}_{m}_{era}.{file_extension}', dpi=400, bbox_inches='tight')
                stack=None
            analysis=None

    elif action=="background_plots_negativelxy":
        print("Make 1D final Plots")
        myTexts = {
            'wmn2g':r'$W\rightarrow \mu \nu$',
            'wen2g':r'$W\rightarrow e \nu$',
            'zmm2g':r'$Z\rightarrow \mu \mu$',
            'zee2g':r'$Z\rightarrow e e$',

            }
        for m in masses:

            
            fig,ax = plt.subplots(nrows=2,ncols=len(analyses),sharex="col",sharey="row",figsize=(25,12),gridspec_kw={'height_ratios':[5,1],
                                                                                                                    'width_ratios':list(np.full(len(analyses),1))})
            plt.subplots_adjust(wspace=0,hspace=0)        
        
            for i,ana in enumerate(analyses):

                #create a plotter that has all MC as data
                analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,brphiphi=brphiphi,signals=[],lifetimes=[],brgg=0.5)
                print(f"Running {i} {ana} m={m} GeV")
                stack=mplhep_plotter(label=analysis_status,data=True,lumi=lumifb[era],com=center_of_mass[era],capsize=0)
                stack.add_plotter(analysis['bkg'][m],label='Background',typeP='background',error_mode='poisson_bootstrap')
                stack.add_plotter(analysis['data'],label="Data",typeP='data',error_mode='poisson')               
                #draw a plot
                stack.hist1d(f"best_2g_dxy_m{m}",cuts[ana][m]['presr']+f"*(best_2g_raw_mass_m{m}>65)",('a','a',len(binning1d_sb[ana])-1,binning1d_sb[ana]),alpha=1.0,xlabel=r"$L_{xy}$",xunits="cm",legend_loc='upper right',show=False,ax=ax[0,i],dndx=False)
                stack.pull1d(f"best_2g_dxy_m{m}",cuts[ana][m]['presr']+f"*(best_2g_raw_mass_m{m}>65)",('a','a',len(binning1d_sb[ana])-1,binning1d_sb[ana]),alpha=1.0,xlabel=r"$L_{xy}$",xunits="cm",legend_loc='upper right',show=False,ax=ax[1,i])
                
                analysis=None
                ax[0,i].text(0.5, 0.95, myTexts[ana] , transform=ax[0,i].transAxes,
                           fontsize=20, ha='center', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
                ax[0,i].margins(x=0,y=0)
                ax[1,i].margins(x=0,y=0)
                
                if i==(len(analyses)-1):
                    ax[0,i].legend(loc='upper left', bbox_to_anchor=(0.01,0.9))
                ax[0,i].tick_params(axis='both', which='major', labelsize=18)
                ax[1,i].tick_params(axis='both', which='major', labelsize=18)
            #then stack backgrounds and then draw band
            lo, hi = ax[0,0].get_ylim()
            ax[0,0].set_ylim(lo, hi + (hi - lo) * 0.25)
            ax[1,0].set_ylim(-3,3)



            ax[0,-1].text(0.1, 0.55, rf"$m_\Phi={m}\ \mathrm{{GeV}}$" , transform=ax[0,-1].transAxes,
                           fontsize=20, ha='left', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))

            #Labels
            mh.cms.label(analysis_status,data=True,rlabel="", ax=ax[0,0], loc=0)
            mh.cms.label(None,exp='',data=True,llabel="", ax=ax[0,-1], loc=0,lumi=lumifb[era],com=center_of_mass[era])
            ax[1,-1].set_xlabel(r"$L_{xy}$ (cm)",fontsize=20)
            ax[0,0].set_ylabel("< Events / cm >",fontsize=20)
            ax[1,0].set_ylabel(r"$(\mathrm{Data}-\mathrm{Pred.})/\sigma$",fontsize=20)
            fig.align_ylabels([ax[0,0], ax[1,0]])
            plt.savefig(f'{outputDir}/sideband_lxy_1D_{m}_{era}.{file_extension}', dpi=400, bbox_inches='tight')
            with open(f'{outputDir}/sideband_lxy_1D_{m}_{era}.pickle', "wb") as file:
                pickle.dump(fig, file)
            


    elif action=="background_plots_negativelxy_x2":
        print("Make 1D final Plots")
        myTexts = {
            'wmn2g':r'$W\rightarrow \mu \nu$',
            'wen2g':r'$W\rightarrow e \nu$',
            'zmm2g':r'$Z\rightarrow \mu \mu$',
            'zee2g':r'$Z\rightarrow e e$',

            }
        for m in masses:
            
            fig,ax = plt.subplots(nrows=2,ncols=len(analyses),sharex="col",sharey="row",figsize=(20,12),gridspec_kw={'height_ratios':[5,1],
                                                                                                                    'width_ratios':list(np.full(len(analyses),1))})
            plt.subplots_adjust(wspace=0,hspace=0)        
        
            for i,ana in enumerate(analyses):

                #create a plotter that has all MC as data
                analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,brphiphi=0.0,signals=[],lifetimes=[],brgg=0.0)
                print(f"Running {i} {ana} m={m} GeV")
                stack=mplhep_plotter(label=analysis_status,data=True,lumi=lumifb[era],com=center_of_mass[era],capsize=0)
                stack.add_plotter(analysis['bkg'][m],label='Background',typeP='background',error_mode='poisson_bootstrap')
                stack.add_plotter(analysis['data'],label="Data",typeP='data',error_mode='poisson')               
                
                #draw a plot
                stack.hist1d(f"best_2g_dxy_m{m}",cuts[ana][m]['presr']+f"*(best_2g_raw_mass_m{m}>62.5)",('a','a',len(binning1d_sb[ana])-1,binning1d_sb[ana]),alpha=1.0,xlabel=r"$L_{xy}$",xunits="cm",legend_loc='upper right',show=False,ax=ax[0,i],dndx=False)
                stack.pull1d(f"best_2g_dxy_m{m}",cuts[ana][m]['presr']+f"*(best_2g_raw_mass_m{m}>62.5)",('a','a',len(binning1d_sb[ana])-1,binning1d_sb[ana]),alpha=1.0,xlabel=r"$L_{xy}$",xunits="cm",legend_loc='upper right',show=False,ax=ax[1,i])


                analysis=None
                ax[0,i].text(0.5, 0.95, myTexts[ana] , transform=ax[0,i].transAxes,
                           fontsize=35, ha='center', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
                ax[0,i].margins(x=0,y=0)
                ax[1,i].margins(x=0,y=0)
                
                if i==(len(analyses)-1):
                    handles, labels = ax[0,-1].get_legend_handles_labels()
                    hl_dict = dict(zip(labels, handles))
                    ordered_handles=[hl_dict['Data']]
                    ordered_labels=['Data']
                    for l,h in hl_dict.items():
                        if l!='Data':
                            ordered_labels.append(l)
                            ordered_handles.append(h)
                    ax[0,-1].legend(ordered_handles,ordered_labels,loc='upper left', bbox_to_anchor=(0.1,0.9),fontsize=35)
                ax[0,i].tick_params(axis='both', which='major', labelsize=27.5)
                ax[1,i].tick_params(axis='both', which='major', labelsize=27.5)
            #then stack backgrounds and then draw band
            lo, hi = ax[0,0].get_ylim()
            ax[0,0].set_ylim(lo, hi + (hi - lo) * (1.0 if 'wmn2g' in analyses else 2.0))
            ax[1,0].set_ylim(-4,4)
            ax[1,0].set_yticks([-2,0,2])



            ax[0,0].text(0.45, 0.8, rf"$m_\Phi={m}\ \mathrm{{GeV}}$" , transform=ax[0,0].transAxes,
                           fontsize=35, ha='left', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))

            #Labels
            mh.cms.label(analysis_status,data=True,rlabel="", ax=ax[0,0], loc=0,fontsize=45)
            mh.cms.label(None,exp='',data=True,llabel="", ax=ax[0,-1], loc=0,lumi=lumifb[era],com=center_of_mass[era],fontsize=35)
            ax[1,-1].set_xlabel(r"$L_{xy}$ (cm)",fontsize=35)
            ax[0,0].set_ylabel("Events / bin ",fontsize=35)
#            ax[1,0].set_ylabel(r"$(\mathrm{Data}-\mathrm{Pred.})/\sigma$",fontsize=24)
            ax[1,0].set_ylabel("Pull",fontsize=35,loc='center')
            fig.align_ylabels([ax[0,0], ax[1,0]])
            anastr = '_'.join(analyses)
            plt.savefig(f'{outputDir}/sideband_lxy_1D_{m}_{era}_{anastr}.{file_extension}', dpi=400, bbox_inches='tight')
            with open(f'{outputDir}/sideband_lxy_1D_{m}_{era}_{anastr}.pickle', "wb") as file:
                pickle.dump(fig, file)

                
                
        #ACTION: Final Plots 1D             
    elif action=="final_mgg_plot":
        brgg=0.5
        print("Make 1D final Plots")
        mySignals={
            'wmn2g':['WH','ttH'],
            'wen2g':['WH','ttH'],
            'zmm2g':['ZH','ggZH'],
            'zee2g':['ZH','ggZH']}
        myTexts = {
            'wmn2g':r'$W\rightarrow \mu \nu$',
            'wen2g':r'$W\rightarrow e \nu$',
            'zmm2g':r'$Z\rightarrow \mu \mu$',
            'zee2g':r'$Z\rightarrow e e$',

            }
        mh.style.use('CMS')       
        for m in masses:    
            fig,ax = plt.subplots(nrows=2,ncols=len(analyses),sharex="col",sharey="row",figsize=(20,12),gridspec_kw={'height_ratios':[5,1],
                                                                                                                    'width_ratios':list(np.full(len(analyses),1))})
            plt.subplots_adjust(wspace=0,hspace=0)        
        
            for i,ana in enumerate(analyses):

                #create a plotter that has all MC as data
                analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,brphiphi=brphiphi,signals=mySignals[ana],lifetimes=lifetimes,brgg=brgg)
                print(f"Running {i} {ana} m={m} GeV")
                stack=mplhep_plotter(label=analysis_status,data=True,lumi=lumifb[era],com=center_of_mass[era],capsize=0)
                stack.add_plotter(analysis['bkg'][m],label='Background',typeP='background',error_mode='poisson_bootstrap')
                for ctau in lifetimes:
                    stack.add_plotter(analysis['signal'][m][ctau]['sum'],label=r" $H \rightarrow\Phi\Phi,  c\tau =$ "+f"{ctau} mm",typeP='signal',error_mode='w2',color=signal_colors[ctau])
                if blinded==False:
                    stack.add_plotter(analysis['data'],label="Data",typeP='data',error_mode='poisson')               
                #draw a plot
                stack.hist1d(f"best_2g_raw_mass_m{m}",cuts[ana][m]['presr'],('a','a',len(binning_mass[ana])-1,binning_mass[ana]),alpha=1.0,xlabel=r"$m_{\gamma \gamma}$",xunits="GeV",legend_loc='upper right',show=False,ax=ax[0,i],dndx=True)
                stack.pull1d(f"best_2g_raw_mass_m{m}",cuts[ana][m]['presr'],('a','a',len(binning_mass[ana])-1,binning_mass[ana]),alpha=1.0,xlabel=r"$m_{\gamma \gamma}$",xunits="GeV",legend_loc='upper right',show=False,ax=ax[1,i])
                analysis=None
                ax[0,i].text(0.5, 0.95, myTexts[ana] , transform=ax[0,i].transAxes,
                           fontsize=35, ha='center', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
                ax[0,i].margins(x=0,y=0)
                ax[1,i].margins(x=0,y=0)
                
                if i==(len(analyses)-1):
                    handles, labels = ax[0,-1].get_legend_handles_labels()
                    hl_dict = dict(zip(labels, handles))
                    ordered_handles=[hl_dict['Data']]
                    ordered_labels=['Data']
                    for l,h in hl_dict.items():
                        if l!='Data':
                            ordered_labels.append(l)
                            ordered_handles.append(h)
                    ax[0,-1].legend(ordered_handles,ordered_labels,loc='upper left', bbox_to_anchor=(0.1,0.9),fontsize=35)
                ax[0,i].tick_params(axis='both', which='major', labelsize=27.5)
                ax[1,i].tick_params(axis='both', which='major', labelsize=27.5)
            #then stack backgrounds and then draw band
            lo, hi = ax[0,0].get_ylim()
            ax[0,0].set_ylim(lo, hi + (hi - lo) * (1.0 if 'wmn2g' in analyses else 2.0))
            ax[1,0].set_ylim(-4,4)
            ax[1,0].set_yticks([-2,0,2])



            ax[0,0].text(0.45, 0.8, rf"$m_\Phi={m}\ \mathrm{{GeV}}$" , transform=ax[0,0].transAxes,
                           fontsize=35, ha='left', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))

            ax[0,0].text(0.45, 0.7, r"$\mathcal{B} (H \rightarrow \Phi\Phi)=0.01$" , transform=ax[0,0].transAxes,
                           fontsize=35, ha='left', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
            ax[0,0].text(0.45, 0.6, r"$\mathcal{B}(\Phi \rightarrow \gamma\gamma)=0.5$" , transform=ax[0,0].transAxes,
                           fontsize=35, ha='left', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
            
            #Labels
            mh.cms.label(analysis_status,data=True,rlabel="", ax=ax[0,0], loc=0,fontsize=45)
            mh.cms.label(None,exp='',data=True,llabel="", ax=ax[0,-1], loc=0,lumi=lumifb[era],com=center_of_mass[era],fontsize=35)
            ax[1,-1].set_xlabel(r"$m_{\gamma\gamma}$ (GeV)",fontsize=35)
            ax[0,0].set_ylabel("Events / bin ",fontsize=35)
#            ax[1,0].set_ylabel(r"$(\mathrm{Data}-\mathrm{Pred.})/\sigma$",fontsize=24)
            ax[1,0].set_ylabel("Pull",fontsize=35,loc='center')
            fig.align_ylabels([ax[0,0], ax[1,0]])
            anastr = '_'.join(analyses)
            plt.savefig(f'{outputDir}/prefit_mgg_{m}_{era}_{anastr}.{file_extension}', dpi=400, bbox_inches='tight')
            with open(f'{outputDir}/prefit_mgg_{m}_{era}_{anastr}.pickle', "wb") as file:
                pickle.dump(fig, file)




        #ACTION: Final Plots 1D             
    elif action=="final_plots_1d":
        brgg=0.5
        print("Make 1D final Plots")
        mySignals={
            'wmn2g':['WH','ttH'],
            'wen2g':['WH','ttH'],
            'zmm2g':['ZH','ggZH'],
            'zee2g':['ZH','ggZH']}
        myTexts = {
            'wmn2g':r'$W\rightarrow \mu \nu$',
            'wen2g':r'$W\rightarrow e \nu$',
            'zmm2g':r'$Z\rightarrow \mu \mu$',
            'zee2g':r'$Z\rightarrow e e$',

            }
        mh.style.use('CMS')

        
        
        
        for m in masses:    
            #create the common axes
            
            fig,ax = plt.subplots(nrows=2,ncols=len(analyses),sharex="col",sharey="row",figsize=(25,12),gridspec_kw={'height_ratios':[5,1],
                                                                                                                    'width_ratios':list(np.full(len(analyses),1))})
            plt.subplots_adjust(wspace=0,hspace=0)        
        
            for i,ana in enumerate(analyses):

                #create a plotter that has all MC as data
                analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,brphiphi=brphiphi,signals=mySignals[ana],lifetimes=lifetimes,brgg=brgg)
                print(f"Running {i} {ana} m={m} GeV")
                stack=mplhep_plotter(label=analysis_status,data=True,lumi=lumifb[era],com=center_of_mass[era],capsize=0)
                stack.add_plotter(analysis['bkg'][m],label='Background',typeP='background',error_mode='poisson_bootstrap')
                for ctau in lifetimes:
                    stack.add_plotter(analysis['signal'][m][ctau]['sum'],label=r" $H \rightarrow\Phi\Phi,  c\tau =$ "+f"{ctau} mm",typeP='signal',error_mode='w2',color=signal_colors[ctau])
                if blinded==False:
                    stack.add_plotter(analysis['data'],label="Data",typeP='data',error_mode='poisson')               
                #draw a plot
                stack.hist1d(f"best_2g_dxy_m{m}",cuts[ana][m]['sr'],('a','a',len(binning1d[ana])-1,binning1d[ana]),alpha=1.0,xlabel=r"$L_{xy}$",xunits="cm",legend_loc='upper right',show=False,ax=ax[0,i],dndx=True)
                stack.pull1d(f"best_2g_dxy_m{m}",cuts[ana][m]['sr'],('a','a',len(binning1d[ana])-1,binning1d[ana]),alpha=1.0,xlabel=r"$L_{xy}$",xunits="cm",legend_loc='upper right',show=False,ax=ax[1,i])
                
                analysis=None
                ax[0,i].text(0.5, 0.95, myTexts[ana] , transform=ax[0,i].transAxes,
                           fontsize=20, ha='center', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
                ax[0,i].margins(x=0,y=0)
                ax[1,i].margins(x=0,y=0)
                
                if i==(len(analyses)-1):
                    ax[0,i].legend(loc='upper left', bbox_to_anchor=(0.01,0.9))
                ax[0,i].tick_params(axis='both', which='major', labelsize=18)
                ax[1,i].tick_params(axis='both', which='major', labelsize=18)
            #then stack backgrounds and then draw band
            lo, hi = ax[0,0].get_ylim()
            ax[0,0].set_ylim(lo, hi + (hi - lo) * 0.25)
            ax[1,0].set_ylim(-3,3)



            ax[0,-1].text(0.1, 0.55, rf"$m_\Phi={m}\ \mathrm{{GeV}}$" , transform=ax[0,-1].transAxes,
                           fontsize=20, ha='left', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))

            ax[0,-1].text(0.1, 0.5, r"$\mathcal{B} (H \rightarrow \Phi\Phi)=0.01$" , transform=ax[0,-1].transAxes,
                           fontsize=20, ha='left', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
            ax[0,-1].text(0.1, 0.45, r"$\mathcal{B}(\Phi \rightarrow \gamma\gamma)=0.5$" , transform=ax[0,-1].transAxes,
                           fontsize=20, ha='left', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
            
            #Labels
            mh.cms.label(analysis_status,data=True,rlabel="", ax=ax[0,0], loc=0)
            mh.cms.label(None,exp='',data=True,llabel="", ax=ax[0,-1], loc=0,lumi=lumifb[era],com=center_of_mass[era])
            ax[1,-1].set_xlabel(r"$L_{xy}$ (cm)",fontsize=20)
            ax[0,0].set_ylabel("< Events / cm >",fontsize=20)
            ax[1,0].set_ylabel(r"$(\mathrm{Data}-\mathrm{Pred.})/\sigma$",fontsize=20)
            fig.align_ylabels([ax[0,0], ax[1,0]])
            if blinded:
                plt.savefig(f'{outputDir}/blinded_prefit_1D_{m}_{era}.{file_extension}', dpi=400, bbox_inches='tight')
            else:
                plt.savefig(f'{outputDir}/prefit_1D_{m}_{era}.{file_extension}', dpi=400, bbox_inches='tight')
            with open(f'{outputDir}/prefit_1D_{m}_{era}.pickle', "wb") as file:
                pickle.dump(fig, file)

    elif action=="final_plots_1d_x2":
        brgg=0.5
        print("Make 1D final Plots")
        mySignals={
            'wmn2g':['WH','ttH'],
            'wen2g':['WH','ttH'],
            'zmm2g':['ZH','ggZH'],
            'zee2g':['ZH','ggZH']}
        myTexts = {
            'wmn2g':r'$W\rightarrow \mu \nu$',
            'wen2g':r'$W\rightarrow e \nu$',
            'zmm2g':r'$Z\rightarrow \mu \mu$',
            'zee2g':r'$Z\rightarrow e e$',

            }
        mh.style.use('CMS')

        
        
        
        for m in masses:    
            #create the common axes
            
            fig,ax = plt.subplots(nrows=2,ncols=len(analyses),sharex="col",sharey="row",figsize=(20,12),gridspec_kw={'height_ratios':[5,1],
                                                                                                                    'width_ratios':list(np.full(len(analyses),1))})
            plt.subplots_adjust(wspace=0,hspace=0)        
        
            for i,ana in enumerate(analyses):

                #create a plotter that has all MC as data
                analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,brphiphi=brphiphi,signals=mySignals[ana],lifetimes=lifetimes,brgg=brgg)
                print(f"Running {i} {ana} m={m} GeV")
                stack=mplhep_plotter(label=analysis_status,data=True,lumi=lumifb[era],com=center_of_mass[era],capsize=0)
                stack.add_plotter(analysis['bkg'][m],label='Background',typeP='background',error_mode='poisson_bootstrap')
                for ctau in lifetimes:
                    stack.add_plotter(analysis['signal'][m][ctau]['sum'],label=r" $H \rightarrow\Phi\Phi,  c\tau =$ "+f"{ctau} mm",typeP='signal',error_mode='w2',color=signal_colors[ctau])
                if blinded==False:
                    stack.add_plotter(analysis['data'],label="Data",typeP='data',error_mode='poisson')               
                #draw a plot
                stack.hist1d(f"best_2g_dxy_m{m}",cuts[ana][m]['sr'],('a','a',len(binning1d[ana])-1,binning1d[ana]),alpha=1.0,xlabel=r"$L_{xy}$",xunits="cm",legend_loc='upper right',show=False,ax=ax[0,i],dndx=False)
                stack.pull1d(f"best_2g_dxy_m{m}",cuts[ana][m]['sr'],('a','a',len(binning1d[ana])-1,binning1d[ana]),alpha=1.0,xlabel=r"$L_{xy}$",xunits="cm",legend_loc='upper right',show=False,ax=ax[1,i])
                
                analysis=None
                ax[0,i].text(0.5, 0.95, myTexts[ana] , transform=ax[0,i].transAxes,
                           fontsize=35, ha='center', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
                ax[0,i].margins(x=0,y=0)
                ax[1,i].margins(x=0,y=0)
                
                if i==(len(analyses)-1):
                    handles, labels = ax[0,-1].get_legend_handles_labels()
                    hl_dict = dict(zip(labels, handles))
                    ordered_handles=[hl_dict['Data']]
                    ordered_labels=['Data']
                    for l,h in hl_dict.items():
                        if l!='Data':
                            ordered_labels.append(l)
                            ordered_handles.append(h)
                    ax[0,-1].legend(ordered_handles,ordered_labels,loc='upper left', bbox_to_anchor=(0.1,0.9),fontsize=35)
                ax[0,i].tick_params(axis='both', which='major', labelsize=27.5)
                ax[1,i].tick_params(axis='both', which='major', labelsize=27.5)
            #then stack backgrounds and then draw band
            lo, hi = ax[0,0].get_ylim()
            ax[0,0].set_ylim(lo, hi + (hi - lo) * (1.0 if 'wmn2g' in analyses else 2.0))
            ax[1,0].set_ylim(-4,4)
            ax[1,0].set_yticks([-2,0,2])



            ax[0,0].text(0.45, 0.8, rf"$m_\Phi={m}\ \mathrm{{GeV}}$" , transform=ax[0,0].transAxes,
                           fontsize=35, ha='left', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))

            ax[0,0].text(0.45, 0.7, r"$\mathcal{B} (H \rightarrow \Phi\Phi)=0.01$" , transform=ax[0,0].transAxes,
                           fontsize=35, ha='left', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
            ax[0,0].text(0.45, 0.6, r"$\mathcal{B}(\Phi \rightarrow \gamma\gamma)=0.5$" , transform=ax[0,0].transAxes,
                           fontsize=35, ha='left', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
            
            #Labels
            mh.cms.label(analysis_status,data=True,rlabel="", ax=ax[0,0], loc=0,fontsize=45)
            mh.cms.label(None,exp='',data=True,llabel="", ax=ax[0,-1], loc=0,lumi=lumifb[era],com=center_of_mass[era],fontsize=35)
            ax[1,-1].set_xlabel(r"$L_{xy}$ (cm)",fontsize=35)
            ax[0,0].set_ylabel("Events / bin ",fontsize=35)
#            ax[1,0].set_ylabel(r"$(\mathrm{Data}-\mathrm{Pred.})/\sigma$",fontsize=24)
            ax[1,0].set_ylabel("Pull",fontsize=35,loc='center')
            fig.align_ylabels([ax[0,0], ax[1,0]])
            anastr = '_'.join(analyses)
            if blinded:
                plt.savefig(f'{outputDir}/blinded_prefit_1D_{m}_{era}_{anastr}.{file_extension}', dpi=400, bbox_inches='tight')
            else:
                plt.savefig(f'{outputDir}/prefit_1D_{m}_{era}_{anastr}.{file_extension}', dpi=400, bbox_inches='tight')
            with open(f'{outputDir}/prefit_1D_{m}_{era}_{anastr}.pickle', "wb") as file:
                pickle.dump(fig, file)
            
    elif action=="check_files":
        mySignals={
            'wmn2g':['WH','ttH'],
            'wen2g':['WH','ttH'],
            'zmm2g':['ZH','ggZH'],
            'zee2g':['ZH','ggZH']}       
        for i,ana in enumerate(analyses):
            analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era='Run2',brphiphi=brphiphi,signals=mySignals[ana],lifetimes=lifetimes)


    #ACTION: Make datacards 1D            
    elif action=="make_datacards_1d":
        import copy
        print("Make Datacards 1D")
        lumiUnc = {'2018': 1.025,
                   '2017': 1.023,
                   '2016': 1.012}
        xsecUnc = {'WH'  : [1+0.005, 1-0.007],
                   'ZH'  : [1+0.038, 1-0.031],
                   'ggZH': [1+0.251, 1-0.189],
                   'ttH' : [1+0.058, 1-0.092]
                   }

        #symmetric
        pdfUnc = {'WH'  : 1.019,
                  'ZH'  : 1.016,
                  'ggZH': 1.024,
                  'ttH' : 1.036,
                  }
        mySignals={
            'wmn2g':['WH2G2Q','ttH2G2Q','WH4G','ttH4G'],
            'wen2g':['WH2G2Q','ttH2G2Q','WH4G','ttH4G'],
            'zmm2g':['ZH2G2Q','ggZH2G2Q','ZH4G','ggZH4G'],
            'zee2g':['ZH2G2Q','ggZH2G2Q','ZH4G','ggZH4G']}
            
        for ana in analyses:
            for e in eras:
                for m in masses:
                    analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=e,brphiphi=brphiphi,masses=[m],signals=mySignals[ana],lifetimes=lifetimes)                        
                    for ibin,binval in enumerate(binning1d[ana][:-1]):
                        dxy_min=binning1d[ana][ibin]
                        dxy_max=binning1d[ana][ibin+1]
                        print(f"Bin  {dxy_min}<= dxy <={dxy_max}")                                
                        cutstring = cuts[ana][m]['sr']+f"*(best_2g_dxy_m{m}>={dxy_min}&&best_2g_dxy_m{m}<{dxy_max})"
                        print(f"Making Datacards for {ana} in era {e} for m={m} GeV")
                        dcmBase = cnc_datacard_maker(outDir=outputDir,binname="abstract",cuts=cutstring)
                        dcmBase.add('data','data',analysis['data'],{})
                        dcmBase.add('bkg','background',analysis['bkg'][m],uncertainties={
                            f"CMS_VHDDP_{ana}_{e}_bin_{ibin}_CR_stats":{'type':'statAsym'},
                            f"CMS_VHDDP_fakeRateUnc_{e}":{'type':'weightAsymm','weightUp':'fakeRate_up','weightDown':'fakeRate_down','weightOrig':'fakeRate_val'},
                            f"CMS_VHDDP_{ana}_{e}_bin_{ibin}_bkg_lowstat":{'type':'zeroRate','value':0.18}

                        },error_mode='poisson_bootstrap')
                    
                        for ctau in lifetimes:
                            dcm=copy.deepcopy(dcmBase)
                            dcm.binname=f"{ana}_m{m}_ctau{ctau}_era{e}_bin{ibin}"
                            print(f"ctau={ctau} mm")
                            for signal in mySignals[ana]:
                                if not (signal in analysis['signal'][m][ctau].keys()):
                                    print(f"Signal not found {signal} skipping")
                                    continue
                                # add a parameter for branching fraction
                                if '2G2Q' in signal:
                                    dcm.rateParameters.append({'name':f"CMS_VHDDP_{ana}_{e}_bin_{ibin}_{signal}_scale",'bin':f"{dcm.binname}",'signal':f'{signal}', 'formula':'(@0*(1.0-@0)/(0.5*0.5)) brgg'})
                                    theory=signal.split('2G2Q')[0]    
                                    
                                elif '4G' in signal:
                                    dcm.rateParameters.append({'name':f"CMS_VHDDP_{ana}_{e}_bin_{ibin}_{signal}_scale",'bin':f"{dcm.binname}",'signal':f'{signal}', 'formula':'(@0*@0)/(0.5*0.5) brgg'})

                                    theory=signal.split('4G')[0]    
                                    
                                signalUncertainties={
                                    f'CMS_lumi_{e}':{'type':'adhoc','kind':'lnN','value':lumiUnc[e]},                                        
                                    f'CMS_{theory}_xsec':{'type':'adhoc','kind':'lnN','value':f"{xsecUnc[theory][1]}/{xsecUnc[theory][0]}"},
                                    f"CMS_VHDDP_{signal}_{ana}_{e}_bin_{ibin}_stats":{'type':'statSym'},                                    
                                    'CMS_pdf':{'type':'adhoc','kind':'lnN','value':pdfUnc[theory]}
                                }
                                #add lepton ID SF, some string manipulations in place
                                leptonSFs = leptonSF[ana]
                                individuals = list(set([x.split('SF_val')[0] for x in leptonSFs.split("*")]))
                                for unc in individuals:
                                    l = unc.replace('Muon','mu').replace('Electron','ele')
                                    signalUncertainties[f'CMS_{l}_{e}']={'type':'weightAsymm','weightUp':leptonSFs.replace(unc+'SF_val',unc+'SF_up'),'weightDown':leptonSFs.replace(unc+'SF_val',unc+'SF_down'),'weightOrig':leptonSFs}

                                #Same for photons
                                photonSFs = photonSF[m]
                                individuals = list(set([x.split('SF_val')[0] for x in photonSFs.split("*")]))

                                for unc in individuals:
                                    l = unc.replace('Photon','gamma')
                                    signalUncertainties[f'CMS_{l}_{e}']={'type':'weightAsymm','weightUp':photonSFs.replace(unc+'SF_val',unc+'SF_up'),'weightDown':photonSFs.replace(unc+'SF_val',unc+'SF_down'),'weightOrig':photonSFs}

                                #now evaluate the systematics because of photon scale and resolution
                                #first define the new photon momenta
                                analysis['signal'][m][ctau][signal].define("Photon_ptSmearUp", "ptSmearUp(Photon_pt, Photon_eta, Photon_phi, Photon_dEsigmaUp)")
                                analysis['signal'][m][ctau][signal].define("Photon_ptSmearDown", "ptSmearDown(Photon_pt, Photon_eta, Photon_phi, Photon_dEsigmaDown)")
                                analysis['signal'][m][ctau][signal].define("Photon_ptScaleUp", "Photon_pt*photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaledown_{era}_val, PHO_scaledown_{era}_bins, 1)".format(era = "2016preVFP" if e=="2016" else e))
                                analysis['signal'][m][ctau][signal].define("Photon_ptScaleDown", "Photon_pt*photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaleup_{era}_val, PHO_scaledown_{era}_bins, 1)".format(era = "2016preVFP" if e=="2016" else e))
                                for unc in ['Scale','Smear']:
                                    for direction in ['Up','Down']:
                                        #calculate the mass and the dxy with new pt
                                        analysis['signal'][m][ctau][signal].define(f"best_2g_{unc}{direction}_m{m}_info", f"kinfit_systematics(Photon_pt{unc}{direction}, Photon_eta, Photon_phi, Photon_isScEtaEB, Photon_isScEtaEE, best_2g_idx1_m{m}, best_2g_idx2_m{m}, {m})")
                                        analysis['signal'][m][ctau][signal].define(f"best_2g_dxy_{unc}{direction}_m{m}", f"best_2g_{unc}{direction}_m{m}_info[0]")
                                        analysis['signal'][m][ctau][signal].define(f"best_2g_raw_mass_{unc}{direction}_m{m}", f"best_2g_{unc}{direction}_m{m}_info[1]")
                                    systName=unc.replace('Scale','scale').replace('Smear','res')     
                                    signalUncertainties[f'CMS_gamma_{systName}_{e}']={'type':'replication','originals':['Photon_pt',f'best_2g_raw_mass_m{m}',f'best_2g_dxy_m{m}'],
                                                                                      'replacementsUp':[f"Photon_pt{unc}Up",f"best_2g_raw_mass_{unc}Up_m{m}",f"best_2g_dxy_{unc}Up_m{m}"],
                                                                                      'replacementsDown':[f"Photon_pt{unc}Down",f"best_2g_raw_mass_{unc}Down_m{m}",f"best_2g_dxy_{unc}Down_m{m}"]}
                                dcm.add(signal,'signal',analysis['signal'][m][ctau][signal],uncertainties=signalUncertainties)
                                    
                                    
                            #write only reasonable cards
                            write_card_signal = sum([(t=='signal') for s,t in  dcm.types.items()])
                            write_card_bkg = sum([(t=='background') for s,t in  dcm.types.items()])
                            dcm.parameters.append(f"brgg extArg 0.5 [0,1]")                                
                            if write_card_signal>0:
                                dcm.write()
                            del dcm
                            
                    analysis=None

        
    #ACTION: Make datacards             
    elif action=="make_datacards":
        print("Make Datacards")
        lumiUnc = {'2018': 1.025,
                   '2017': 1.023,
                   '2016': 1.012}
        xsecUnc = {'WH'  : [1+0.005, 1-0.007],
                   'ZH'  : [1+0.038, 1-0.031],
                   'ggZH': [1+0.251, 1-0.189],
                   'ttH' : [1+0.058, 1-0.092]}

        #symmetric
        pdfUnc = {'WH'  : 1.019,
                  'ZH'  : 1.016,
                  'ggZH': 1.024,
                  'ttH' : 1.036}
        mySignals={
            'wmn2g':['WH','ttH'],
            'wen2g':['WH','ttH'],
            'zmm2g':['ZH','ggZH'],
            'zee2g':['ZH','ggZH']}
            
        for ana in analyses:
            for e in eras:
                for m in masses:                    
                    for ctau in lifetimes:
                        analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=e,brphiphi=brphiphi,masses=[m],signals=signals,lifetimes=[ctau])                        
                        print(f"Making Datacards for {ana} in era {e} for m={m} GeV and ctau={ctau} mm")
                        for ibinx,bin_setup in enumerate(binning[ana][m]):
                            mass_min = bin_setup[0][0]
                            mass_max=  bin_setup[0][1]
                            for ibiny,biny in enumerate(bin_setup[1][:-1]):
                                dxy_min=bin_setup[1][ibiny]
                                dxy_max=bin_setup[1][ibiny+1]
                                print(f"Bin {mass_min}<=m<{mass_max} {dxy_min}<= dxy <={dxy_max}")
                                
                                cutstring = cuts[ana][m]['sr']+f"*(best_2g_raw_mass_m{m}>={mass_min}&&best_2g_raw_mass_m{m}<{mass_max})*(best_2g_dxy_m{m}>={dxy_min}&&best_2g_dxy_m{m}<{dxy_max})"

                                dcm = cnc_datacard_maker(outDir=outputDir,binname=f"{ana}_m{m}_ctau{ctau}_era{e}_binm_{ibinx}_bindxy_{ibiny}",cuts=cutstring)
                                dcm.add('data','data',analysis['data'],{})

                                for signal in mySignals[ana]:
                                    if not (signal in analysis['signal'][m][ctau].keys()):
                                        print(f"Signal not found {signal} skipping")
                                        continue
                                    signalUncertainties={
                                        f'CMS_lumi_{e}':{'type':'adhoc','kind':'lnN','value':lumiUnc[e]},                                        
                                        f'CMS_{signal}_xsec':{'type':'adhoc','kind':'lnN','value':f"{xsecUnc[signal][1]}/{xsecUnc[signal][0]}"},
                                        'CMS_pdf':{'type':'adhoc','kind':'lnN','value':pdfUnc[signal]}
                                    }
                                    #add lepton ID SF, some string manipulations in place
                                    leptonSFs = leptonSF[ana]
                                    individuals = list(set([x.split('SF_val')[0] for x in leptonSFs.split("*")]))
                                    for unc in individuals:
                                        l = unc.replace('Muon','mu').replace('Electron','ele')
                                        signalUncertainties[f'CMS_{l}_{e}']={'type':'weightAsymm','weightUp':leptonSFs.replace(unc+'SF_val',unc+'SF_up'),'weightDown':leptonSFs.replace(unc+'SF_val',unc+'SF_down'),'weightOrig':leptonSFs}

                                    #Same for photons
                                    photonSFs = photonSF[m]
                                    individuals = list(set([x.split('SF_val')[0] for x in photonSFs.split("*")]))

                                    for unc in individuals:
                                        l = unc.replace('Photon','gamma')
                                        signalUncertainties[f'CMS_{l}_{e}']={'type':'weightAsymm','weightUp':photonSFs.replace(unc+'SF_val',unc+'SF_up'),'weightDown':photonSFs.replace(unc+'SF_val',unc+'SF_down'),'weightOrig':photonSFs}

                                    #now evaluate the systematics because of photon scale and resolution
                                    #first define the new photon momenta
                                    analysis['signal'][m][ctau][signal].define("Photon_ptSmearUp", "ptSmearUp(Photon_pt, Photon_eta, Photon_phi, Photon_dEsigmaUp)")
                                    analysis['signal'][m][ctau][signal].define("Photon_ptSmearDown", "ptSmearDown(Photon_pt, Photon_eta, Photon_phi, Photon_dEsigmaDown)")
                                    analysis['signal'][m][ctau][signal].define("Photon_ptScaleUp", "Photon_pt*photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaledown_{era}_val, PHO_scaledown_{era}_bins, 1)".format(era = "2016preVFP" if e=="2016" else e))
                                    analysis['signal'][m][ctau][signal].define("Photon_ptScaleDown", "Photon_pt*photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaleup_{era}_val, PHO_scaledown_{era}_bins, 1)".format(era = "2016preVFP" if e=="2016" else e))
                                    for unc in ['Scale','Smear']:
                                        for direction in ['Up','Down']:
                                            #calculate the mass and the dxy with new pt
                                            analysis['signal'][m][ctau][signal].define(f"best_2g_{unc}{direction}_m{m}_info", f"kinfit_systematics(Photon_pt{unc}{direction}, Photon_eta, Photon_phi, Photon_isScEtaEB, Photon_isScEtaEE, best_2g_idx1_m{m}, best_2g_idx2_m{m}, {m})")
                                            analysis['signal'][m][ctau][signal].define(f"best_2g_dxy_{unc}{direction}_m{m}", f"best_2g_{unc}{direction}_m{m}_info[0]")
                                            analysis['signal'][m][ctau][signal].define(f"best_2g_raw_mass_{unc}{direction}_m{m}", f"best_2g_{unc}{direction}_m{m}_info[1]")
                                        systName=unc.replace('Scale','scale').replace('Smear','res')     
                                        signalUncertainties[f'CMS_gamma_{systName}_{e}']={'type':'replication','originals':['Photon_pt',f'best_2g_raw_mass_m{m}',f'best_2g_dxy_m{m}'],
                                                                                      'replacementsUp':[f"Photon_pt{unc}Up",f"best_2g_raw_mass_{unc}Up_m{m}",f"best_2g_dxy_{unc}Up_m{m}"],
                                                                                      'replacementsDown':[f"Photon_pt{unc}Down",f"best_2g_raw_mass_{unc}Down_m{m}",f"best_2g_dxy_{unc}Down_m{m}"]}
                                    dcm.add(signal,'signal',analysis['signal'][m][ctau][signal],uncertainties=signalUncertainties)
                                    
                                    
                                dcm.add('bkg','background',analysis['bkg'][m],uncertainties={
                                    f"CMS_DDP_{ana}_{e}_binm_{ibinx}_bindxy_{ibiny}_CR_stats":{'type':'statAsym'},
                                    f"CMS_DDP_fakeRateUnc_{e}":{'type':'weightAsymm','weightUp':'fakeRate_up','weightDown':'fakeRate_down','weightOrig':'fakeRate_val'},
                                    f"CMS_DDP_{ana}_{e}_binm_{ibinx}_bindxy_{ibiny}_bkg_lowstat":{'type':'zeroRate','value':0.18}
                                    
                                },error_mode='poisson_bootstrap')
                                #write only reasonable cards
                                write_card_signal = sum([(t=='signal') for s,t in  dcm.types.items()])
                                write_card_bkg = sum([(t=='background') for s,t in  dcm.types.items()])
                                
                                if write_card_signal>0:
                                    dcm.write()
                                    
                        analysis=None


    elif action=="limit_plots":
        band2sigma_color='#85d2fb'  
        band1sigma_color='#ffde9c'
        brgg=[0.01,0.02,0.05,0.1,0.5,1]
        combined=getLimitData(outputDir='VHresults',prefix='asymptotics',lifetimes=lifetimes,brgg=brgg)
        for br in brgg:
            print(f"Unrolled Brazilian flag as a function of m and lifetime for BR={br}")
        
            plotter=limit_plotter(label=analysis_status,lumi=lumifb['Run2'],data=True,com=13,text=None,scale=brphiphi)
            fig,ax = plt.subplots(1,len(lifetimes),sharey=True,figsize=(7.5*len(lifetimes),20))
            plt.subplots_adjust(wspace=0)
        
            for i,ctau in enumerate(lifetimes):
                xdim = int(i%3)
                ydim = int(i/3)
            
                plotter.brazilian_flag(combined[((combined['ctau']==ctau)& (combined['br']==br))],'m',band2sigma_color=band2sigma_color,band1sigma_color=band1sigma_color,ax=ax[i],show=False,quiet=True)

                #            plotter.expected_vs_observed(combined[((combined['ctau']==ctau)& (combined['br']==1.0))],'m',col='grey',ax=ax[i],show=False,quiet=True)
                #            plotter.expected_vs_observed(combined[((combined['ctau']==ctau)& (combined['br']==0.1))],'m',col='grey',ax=ax[i],show=False,quiet=True)
            
                ax[i].text(0.5, 0.02,rf"$c\tau = {ctau}\ \mathrm{{mm}}$", transform=ax[i].transAxes,
                           fontsize=40, ha='center', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
                ax[i].tick_params(axis='both', which='major', labelsize=40)
            ax[0].legend(loc='upper right',fontsize=40)
            ax[-1].set_xlabel(r"$m_\Phi$ (GeV)",fontsize=40)
            ax[0].set_ylabel(r'95% upper limit on $\mathcal{B}(H\rightarrow\Phi\Phi)$',fontsize=40)
            mh.cms.label(analysis_status,data=True,rlabel="", ax=ax[0], loc=0,fontsize=60)
            mh.cms.label(None,exp='',data=True,llabel="", ax=ax[-1], loc=0,lumi=lumifb['Run2'],com=13,fontsize=40)
            ax[0].set_yscale('log')
#            ax[0].set_ylim(0.0005,10)
            ax[0].text(0.5,0.7,rf'$\mathcal{{B}}(\Phi\rightarrow\gamma\gamma)={br}$',transform=ax[0].transAxes,
                   fontsize=40, ha='center', va='center', 
                   bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))

            brstr=str(br).replace('.','_')
            
            fig.savefig(f'{outputDir}/combined_limits_unrolled_m_ctau_br{brstr}.{file_extension}', dpi=400, bbox_inches='tight')
            with open(f'{outputDir}/combined_limits_unrolled_m_ctau_br{brstr}.pickle', "wb") as file:
                pickle.dump(fig, file)
            plt.close()
        for ctau in lifetimes:
            print(f"Unrolled Brazilian flag as a function of BR->gg and mtau for ctau={ctau} ")        
            brToPlot=[0.01,0.02,0.05,0.1,0.5,1.0]
            plotter=limit_plotter(label=analysis_status,lumi=lumifb['Run2'],data=True,com=13,text=None,scale=brphiphi)            
            fig,ax = plt.subplots(1,len(brToPlot),sharey=True,figsize=(7.5*len(brToPlot),20))
            plt.subplots_adjust(wspace=0)
        
            for i,br in enumerate(brToPlot):
                xdim = int(i%3)
                ydim = int(i/3)
            
                plotter.brazilian_flag(combined[((combined['br']==br)& (combined['ctau']==ctau))],'m',band2sigma_color=band2sigma_color,band1sigma_color=band1sigma_color,ax=ax[i],show=False,quiet=True)
                ax[i].text(0.5, 0.02,rf"$\mathcal{{B}}(\Phi\rightarrow\gamma\gamma)={br}$", transform=ax[i].transAxes,
                       fontsize=40, ha='center', va='center', 
                       bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
                ax[i].tick_params(axis='both', which='major', labelsize=40)
            ax[-1].legend(loc='upper right',fontsize=40)
            ax[-1].set_xlabel(r"$m_\Phi$ (GeV)",fontsize=40)
            ax[0].set_ylabel(r'95% upper limit on $\mathcal{B}(H\rightarrow\Phi\Phi)$',fontsize=40)
            mh.cms.label(analysis_status,data=True,rlabel="", ax=ax[0], loc=0,fontsize=60)
            mh.cms.label(None,exp='',data=True,llabel="", ax=ax[-1], loc=0,lumi=lumifb['Run2'],com=13,fontsize=40)
            ax[0].set_yscale('log')
            ax[0].set_ylim(0.0005,10)
            
            ax[-1].text(0.5,0.7,rf'$c\tau = {ctau}\ \mathrm{{mm}}$',transform=ax[-1].transAxes,
                        fontsize=40, ha='center', va='center', 
                        bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
   
            fig.savefig(f'{outputDir}/combined_limits_unrolled_m_br_ctau{ctau}.{file_extension}', dpi=400, bbox_inches='tight')
            with open(f'{outputDir}/combined_limits_unrolled_m_br_ctau{ctau}.pickle', "wb") as file:
                pickle.dump(fig, file)
            plt.close()

        #expected vs observed for each analysis         
        brgg=[0.01,0.02,0.05,0.1,0.5,1]       
        combined=getLimitData(outputDir='VHresults',prefix='asymptotics',lifetimes=lifetimes,brgg=brgg)       
        wmn=getLimitData(outputDir='VHresults',prefix='asymptotics_wmn2g',lifetimes=lifetimes,brgg=brgg)
        wen=getLimitData(outputDir='VHresults',prefix='asymptotics_wen2g',lifetimes=lifetimes,brgg=brgg)
        zmm=getLimitData(outputDir='VHresults',prefix='asymptotics_zmm2g',lifetimes=lifetimes,brgg=brgg)
        zee=getLimitData(outputDir='VHresults',prefix='asymptotics_zee2g',lifetimes=lifetimes,brgg=brgg)
        for br in brgg:
            print("Unrolled overlayed expected vs observed")
        
            plotter=limit_plotter(label="Preliminary",lumi=lumifb['Run2'],data=True,com=13,text=None,scale=brphiphi)
            fig,ax = plt.subplots(1,len(lifetimes),sharey=True,figsize=(7.5*len(lifetimes),20))
            plt.subplots_adjust(wspace=0)
        
            for i,ctau in enumerate(lifetimes):
                xdim = int(i%3)
                ydim = int(i/3)
            
                #plotter.brazilian_flag(combined[((combined['ctau']==ctau)& (combined['br']==br))],'m',band2sigma_color=band2sigma_color,band1sigma_color=band1sigma_color,ax=ax[i],show=False,quiet=True)
                plotter.expected_vs_observed(wmn[((wmn['ctau']==ctau)& (wmn['br']==br))],'m',col='#e62d1a',ax=ax[i],show=False,quiet=True,label_obs = r"$W\rightarrow \mu\nu$")
                plotter.expected_vs_observed(wen[((wen['ctau']==ctau)& (wen['br']==br))],'m',col='#282d65',ax=ax[i],show=False,quiet=True,label_obs = r"$W\rightarrow e\nu$")
                plotter.expected_vs_observed(zmm[((zmm['ctau']==ctau)& (zmm['br']==br))],'m',col='#98c565',ax=ax[i],show=False,quiet=True,label_obs = r"$Z\rightarrow \mu\mu$")
                plotter.expected_vs_observed(zee[((zee['ctau']==ctau)& (zee['br']==br))],'m',col='#636470',ax=ax[i],show=False,quiet=True,label_obs = r"$Z\rightarrow ee$")
                plotter.expected_vs_observed(combined[((combined['ctau']==ctau)& (combined['br']==br))],'m',col='black',ax=ax[i],show=False,quiet=True,label_obs = "Combined",label_exp = "Expected")

            
                ax[i].text(0.5, 0.02,rf"$c\tau = {ctau}\ \mathrm{{mm}}$", transform=ax[i].transAxes,
                           fontsize=40, ha='center', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
                ax[i].tick_params(axis='both', which='major', labelsize=40)
            ax[0].legend(loc='upper right',fontsize=40)
            ax[-1].set_xlabel(r"$m_\Phi$ (GeV)",fontsize=40)
            ax[0].set_ylabel(r'95% upper limit on $\mathcal{B}(H\rightarrow\Phi\Phi)$',fontsize=40)
            mh.cms.label(analysis_status,data=True,rlabel="", ax=ax[0], loc=0,fontsize=60)
            mh.cms.label(None,exp='',data=True,llabel="", ax=ax[-1], loc=0,lumi=lumifb['Run2'],com=13,fontsize=40)
            ax[0].set_yscale('log')
            ax[0].set_ylim(0.0005,10)
            ax[0].text(0.5,0.55,rf'$\mathcal{{B}}(\Phi\rightarrow\gamma\gamma)={br}$',transform=ax[0].transAxes,
                   fontsize=40, ha='center', va='center', 
                   bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))

    
            fig.savefig(f'{outputDir}/obs_vs_exp_limits_unrolled_m_ctau_br{br}.{file_extension}', dpi=400, bbox_inches='tight')
            with open(f'{outputDir}/obs_vs_exp_limits_unrolled_m_ctau_br{br}.pickle', "wb") as file:
                pickle.dump(fig, file)
            plt.close()
                        


        #expected vs observed for each year         
        brgg=[0.01,0.02,0.05,0.1,0.5,1]
        combined=getLimitData(outputDir='VHresults',prefix='asymptotics',lifetimes=lifetimes,brgg=brgg)       
        era2016=getLimitData(outputDir='VHresults',prefix='asymptotics_2016',lifetimes=lifetimes,brgg=brgg)
        era2017=getLimitData(outputDir='VHresults',prefix='asymptotics_2017',lifetimes=lifetimes,brgg=brgg)
        era2018=getLimitData(outputDir='VHresults',prefix='asymptotics_2018',lifetimes=lifetimes,brgg=brgg)
        for br in brgg:
            print("Unrolled overlayed expected vs observed")
        
            plotter=limit_plotter(label="Preliminary",lumi=lumifb['Run2'],data=True,com=13,text=None,scale=brphiphi)
            fig,ax = plt.subplots(1,len(lifetimes),sharey=True,figsize=(7.5*len(lifetimes),20))
            plt.subplots_adjust(wspace=0)
        
            for i,ctau in enumerate(lifetimes):
                xdim = int(i%3)
                ydim = int(i/3)
            

                plotter.expected_vs_observed(era2016[((era2016['ctau']==ctau)& (era2016['br']==br))],'m',col='#e62d1a',ax=ax[i],show=False,quiet=True,label_obs = "2016")
                plotter.expected_vs_observed(era2017[((era2017['ctau']==ctau)& (era2017['br']==br))],'m',col='#282d65',ax=ax[i],show=False,quiet=True,label_obs = "2017")
                plotter.expected_vs_observed(era2018[((era2018['ctau']==ctau)& (era2018['br']==br))],'m',col='#98c565',ax=ax[i],show=False,quiet=True,label_obs = "2018")
                plotter.expected_vs_observed(combined[((combined['ctau']==ctau)& (combined['br']==br))],'m',col='black',ax=ax[i],show=False,quiet=True,label_obs = "Combined",label_exp = "Expected")

            
                ax[i].text(0.5, 0.02,rf"$c\tau = {ctau}\ \mathrm{{mm}}$", transform=ax[i].transAxes,
                           fontsize=40, ha='center', va='center', 
                           bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))
                ax[i].tick_params(axis='both', which='major', labelsize=40)
            ax[0].legend(loc='upper right',fontsize=40)
            ax[-1].set_xlabel(r"$m_\Phi$ (GeV)",fontsize=40)
            ax[0].set_ylabel(r'95% upper limit on $\mathcal{B}(H\rightarrow\Phi\Phi)$',fontsize=40)
            mh.cms.label(analysis_status,data=True,rlabel="", ax=ax[0], loc=0,fontsize=60)
            mh.cms.label(None,exp='',data=True,llabel="", ax=ax[-1], loc=0,lumi=lumifb['Run2'],com=13,fontsize=40)
            ax[0].set_yscale('log')
            ax[0].set_ylim(0.0005,10)
            
            ax[0].text(0.5,0.55,rf'$\mathcal{{B}}(\Phi\rightarrow\gamma\gamma)={br}$',transform=ax[0].transAxes,
                   fontsize=40, ha='center', va='center', 
                   bbox=dict(facecolor='white',edgecolor='none', alpha=0.5))

    
            fig.savefig(f'{outputDir}/obs_vs_exp_year_limits_unrolled_m_ctau_br{br}.{file_extension}', dpi=400, bbox_inches='tight')
            with open(f'{outputDir}/obs_vs_exp_year_limits_unrolled_m_ctau_br{br}.pickle', "wb") as file:
                pickle.dump(fig, file)
            plt.close()
            
