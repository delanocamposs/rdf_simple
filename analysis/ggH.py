import ROOT
import random
import numpy as np
ROOT.gInterpreter.Declare('#include "analysis/ddp_vertex.h"')
ROOT.gInterpreter.Declare('#include "common/scaleFactors.h"')

opts = ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"
opts.fOverwriteIfExists = True

from common.pyhelpers import load_meta_data


#cols = "best_3g.*|best_4g.*|sample_.*|^Photon_.*|^Muon_.*|^Z_.*|Weight.*|^Gen.*|^weight.*|^TrigObj_.*|^event.*|^Electron_.*|^Pileup_.*|^run.*"

iso = {'2016': 'Photon_pfRelIso03_all',
       '2017': 'Photon_pfRelIso03_all',
       '2018': 'Photon_pfRelIso03_all',
       '2022': 'Photon_pfRelIso03_all_quadratic',
       '2023': 'Photon_pfRelIso03_all_quadratic',
       '2024': 'Photon_pfRelIso03_all_quadratic'}
iso_scouting={
        "2022":"",
        "2023":"",
        "2024":"ScoutingPhoton_ecalIso"}

def get_ID_val(var, ID_type):
    cuts_ID={
            "EGM":{"endcap":{"iso":[0.1], 
                             "hoe":[0.0590], 
                             "sieie":[0.0272]}, 
                   "barrel":{"iso":[0.1], 
                             "hoe":[0.04596], 
                             "sieie":[0.0106]}},
            "custom":{"endcap":{"iso":[0.3], 
                                "hoe":[0.2], 
                                "sieie":[0.045]}, 
                      "barrel":{"iso":[0.25], 
                                "hoe":[0.3], 
                                "sieie":[0.035]}}
            }

    return cuts_ID[ID_type]['barrel'][var][0], cuts_ID[ID_type]['endcap'][var][0]
    
# Common Object ID:
def muonAna(dataframe):

    # Common Muon ID definitions (No isolation)
    muons = dataframe.Define("loose_muon", "Muon_looseId==1&&abs(Muon_eta)<2.4&&abs(Muon_dxy)<0.2&&abs(Muon_dz)<0.5&&Muon_pt>10&&Muon_pfIsoId>1")
    muons = muons.Define("tight_muon", "loose_muon&&Muon_tightId&&Muon_pfIsoId>3")
    muons = muons.Define("veto_muon", "Muon_pt>5&&abs(Muon_eta)<2.4&&abs(Muon_dxy)<0.2&&abs(Muon_dz)<0.5&&(loose_muon==0)&&(tight_muon==0)")
    muons = muons.Define("Muon_nloose", "Sum(loose_muon)")
    muons = muons.Define("Muon_ntight", "Sum(tight_muon)")
    muons = muons.Define("Muon_nveto", "Sum(veto_muon)")
    return muons

def electronAna(dataframe):

    # Common Electron ID definitions
    electrons = dataframe.Define("loose_electron", "Electron_pt>15&&abs(Electron_eta)<2.5&&(abs(Electron_eta)>1.57||abs(Electron_eta)<1.44)&&abs(Electron_dxy)<0.2&&abs(Electron_dz)<0.2&&Electron_lostHits<2&&Electron_convVeto&&Electron_cutBased>0")
    electrons = electrons.Define("tight_electron", "loose_electron&&Electron_cutBased>3")
    electrons = electrons.Define("veto_electron", "Electron_pt>55&&abs(Electron_eta)<2.5&&(abs(Electron_eta)>1.57||abs(Electron_eta)<1.44)&&abs(Electron_dxy)<0.2&&abs(Electron_dz)<0.2&&Electron_lostHits<2&&Electron_convVeto&&(tight_electron==0)&&(loose_electron==0)")
    electrons = electrons.Define("Electron_nloose", "Sum(loose_electron)")
    electrons = electrons.Define("Electron_ntight", "Sum(tight_electron)")
    electrons = electrons.Define("Electron_nveto", "Sum(veto_electron)")
    return electrons

def photonAnaScouting(dataframe):
    # Overlap with loose leptons
    #photons = dataframe.Define("Photon_muOverlap", "overlapClean(Photon_phi, Photon_eta, Muon_phi[loose_muon], Muon_eta[loose_muon])")
    #photons = photons.Define("Photon_eleOverlap", "overlapClean(Photon_phi, Photon_eta, Electron_phi[loose_electron], Electron_eta[loose_electron])")
    #photons = photons.Define("Photon_overlap", "Photon_muOverlap||Photon_eleOverlap")
    photons=dataframe.Define("Photon_isScEtaEE", "(abs(ScoutingPhoton_eta)>1.3)&&(abs(ScoutingPhoton_eta)<2.5)")
    photons=photons.Define("Photon_isScEtaEB", "abs(ScoutingPhoton_eta)<1.29")
    photons = photons.Define("Photon_preselection", "ScoutingPhoton_pt>20&&abs(ScoutingPhoton_eta)<2.5&&(abs(ScoutingPhoton_eta)>1.57||abs(ScoutingPhoton_eta)<1.44)&&(Photon_isScEtaEE||Photon_isScEtaEB)")
    #photons = photons.Define("Photon_rho", "fixedGridRhoFastjetAll")
    #photons=photons.Define("Photon_PassPhIso" , "passPhIso(Photon_vidNestedWPBitmap)")
    sieie1, sieie2 = get_ID_val("sieie", "custom")
    hoe1, hoe2 = get_ID_val("hoe", "custom")
    photons = photons.Define("Photon_IdNoIso",f"((Photon_isScEtaEB&&ScoutingPhoton_hOverE<{hoe1}&&ScoutingPhoton_sigmaIetaIeta<{sieie1})||(Photon_isScEtaEE&&ScoutingPhoton_hOverE<{hoe2}&&ScoutingPhoton_sigmaIetaIeta<{sieie2}))")
    return photons    
    
def photonAna(dataframe):
    # Overlap with loose leptons
    #photons = dataframe.Define("Photon_muOverlap", "overlapClean(Photon_phi, Photon_eta, Muon_phi[loose_muon], Muon_eta[loose_muon])")
    #photons = photons.Define("Photon_eleOverlap", "overlapClean(Photon_phi, Photon_eta, Electron_phi[loose_electron], Electron_eta[loose_electron])")
    #photons = photons.Define("Photon_overlap", "Photon_muOverlap||Photon_eleOverlap")

    photons = dataframe.Define("Photon_preselection", "Photon_pt>20&&!Photon_pixelSeed&&abs(Photon_eta)<2.5&&(abs(Photon_eta)>1.57||abs(Photon_eta)<1.44)&&(Photon_isScEtaEE||Photon_isScEtaEB)")
    #photons = photons.Define("Photon_rho", "fixedGridRhoFastjetAll")
    photons=photons.Define("Photon_PassPhIso" , "passPhIso(Photon_vidNestedWPBitmap)")
    sieie1, sieie2 = get_ID_val("sieie", "custom")
    hoe1, hoe2 = get_ID_val("hoe", "custom")
    photons = photons.Define("Photon_IdNoIso",f"((Photon_isScEtaEB&&Photon_hoe<{hoe1}&&Photon_sieie<{sieie1})||(Photon_isScEtaEE&&Photon_hoe<{hoe2}&&Photon_sieie<{sieie2}))")
    return photons    

def save_report(df, report_name, sample, opts, actions):
        report = ROOT.RDataFrame(1)  # Create a dummy dataframe with one entry
        r = df.Report()
        for cut in r:
            # Define new columns for total and passing entries for each cut
            report = report.Define(f"report_{cut.GetName()}_all", f"{cut.GetAll()}")
            report = report.Define(f"report_{cut.GetName()}_pass", f"{cut.GetPass()}")
        # Append the snapshot action to the actions list
        actions.append(report.Snapshot(report_name, f"{sample}.root", "", opts))

def ggH(data,phi_mass,sample):
    def emulate_scouting_trigger(df):
        ROOT.gInterpreter.Declare(r'''
        #include "ROOT/RVec.hxx"
        #include <cmath>
        using namespace ROOT;
        using namespace ROOT::VecOps;
        
        bool PassL1DoubleEG(const RVec<float>& pt,
                            const RVec<float>& eta,
                            const RVec<float>& phi,
                            float thrLeading, float thrSubleading) {
            const size_t n = pt.size();
            if (n < 2) return false;
            if (eta.size() != n || phi.size() != n) return false;
            for (size_t i = 0; i < n; ++i) {
                if (std::abs(eta[i]) > 1.2f) continue;
                for (size_t j = i + 1; j < n; ++j) {
                    if (std::abs(eta[j]) > 1.2f) continue;
                    const bool passEt =
                        (pt[i] >= thrLeading && pt[j] >= thrSubleading) ||
                        (pt[j] >= thrLeading && pt[i] >= thrSubleading);
                    if (!passEt) continue;
                    const float dr = DeltaR(eta[i], eta[j], phi[i], phi[j]);
                    if (dr < 0.6f) return true;
                }
            }
            return false;
        }
        ''')
        df=df.Define("Pass_L1_DoubleEG15_11", "PassL1DoubleEG(Photon_pt, Photon_eta, Photon_phi, 15.0f, 11.0f)")
        df=df.Define("Pass_L1_DoubleEG16_11", "PassL1DoubleEG(Photon_pt, Photon_eta, Photon_phi, 16.0f, 11.0f)")
        df=df.Define("Pass_L1_DoubleEG17_11", "PassL1DoubleEG(Photon_pt, Photon_eta, Photon_phi, 17.0f, 11.0f)")
        df=df.Define("Pass_L1_DoubleEG_OR",   "Pass_L1_DoubleEG15_11 || Pass_L1_DoubleEG16_11 || Pass_L1_DoubleEG17_11")
        return df


    cols = "best_.*|sample_.*|^Photon_.*|^Electron_.*|Weight.*|^Gen.*|^weight.*|^TrigObj_.*|^event.*|^Pileup_.*|^run.*|gen.*|.*LHE.*|^PV.*|luminosity|Block|genWeight|HLT_passed|sorted_photon_pt|Pass_L1_DoubleEG15_11|Pass_L1_DoubleEG16_11|Pass_L1_DoubleEG17_11|Pass_L1_DoubleEG_OR"
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    ggH=dataframe["Events"]


    if data["isMC"]:
        ggH = ggH.Define("Pileup_weight", "getPUweight(Pileup_nPU, puWeight_{}, sample_isMC)".format(data["era"]))

    

    #Filter out muons above 10Gev and electrons above 15GeV
    #ggH=electronAna(ggH)
    #ggH=muonAna(ggH)
    
    #ggH=ggH.Filter("Sum(loose_muon==1)==0",'muon_veto')
    #ggH=ggH.Filter("Sum(loose_electron==1)==0",'electron_veto')

    #Add photon preselection and common photon ID definitions 
    scouting=0
    if dataframe['isScouting']==1:
        scouting=1
        ggH=photonAnaScouting(ggH)
        ggH4g=ggH.Filter('nScoutingPhoton>3','at_least_4_photons')
    else:
        ggH=photonAna(ggH)
        ggH4g=ggH.Filter('nPhoton>3','at_least_4_photons')
    
    ggH4g=emulate_scouting_trigger(ggH4g)
    ggH4g=ggH4g.Filter('Sum(Photon_preselection==1)>3','at_least_3_preselected_photons')

    def four_gamma(df, mass, scouting):
        if scouting==0:
            df=df.Define('raw_best_4g_m{}'.format(mass),"best_4gamma(Photon_pt,Photon_eta,Photon_phi,Photon_isScEtaEB,Photon_isScEtaEE,Photon_preselection,Photon_IdNoIso,Photon_corrIso_m{},{})".format(mass,float(mass)))
        else:
            df=df.Define('raw_best_4g_m{}'.format(mass),"best_4gamma(ScoutingPhoton_pt,ScoutingPhoton_eta,ScoutingPhoton_phi,Photon_isScEtaEB,Photon_isScEtaEE,Photon_preselection,Photon_IdNoIso,Photon_corrIso_m{},{})".format(mass,float(mass)))
        return df

    def isolation_vars(df, mass, scouting):
        if scouting==0:
            #pass isolation criteria using bitmap by EGM here: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
            df=df.Define("Photon_passBitMap_loose_iso_gamma1_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx1_m{m}]==1".format(m=mass))
            df=df.Define("Photon_passBitMap_loose_iso_gamma2_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx2_m{m}]==1".format(m=mass))
            df=df.Define("Photon_passBitMap_loose_iso_gamma3_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx3_m{m}]==1".format(m=mass))
            df=df.Define("Photon_passBitMap_loose_iso_gamma4_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx4_m{m}]==1".format(m=mass))
            df=df.Define("Photon_passBitMap_medium_iso_gamma1_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx1_m{m}]==2".format(m=mass))
            df=df.Define("Photon_passBitMap_medium_iso_gamma2_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx2_m{m}]==2".format(m=mass))
            df=df.Define("Photon_passBitMap_medium_iso_gamma3_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx3_m{m}]==2".format(m=mass))
            df=df.Define("Photon_passBitMap_medium_iso_gamma4_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx4_m{m}]==2".format(m=mass))
            df=df.Define("Photon_passBitMap_tight_iso_gamma1_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx1_m{m}]==3".format(m=mass))
            df=df.Define("Photon_passBitMap_tight_iso_gamma2_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx2_m{m}]==3".format(m=mass))
            df=df.Define("Photon_passBitMap_tight_iso_gamma3_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx3_m{m}]==3".format(m=mass))
            df=df.Define("Photon_passBitMap_tight_iso_gamma4_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx4_m{m}]==3".format(m=mass))
 
            #EVENT-LEVEL BOOLEAN IF BEST 4 PHOTONS PASS LOOSE EGM ISOLATION. THIS IS THE RECOMMENDED METHOD TO ISOLATION ID.
            df=df.Define("best_4g_passBitMap_loose_iso_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx1_m{m}]==1 && Photon_PassPhIso[best_4g_idx2_m{m}]==1 && Photon_PassPhIso[best_4g_idx3_m{m}]==1 && Photon_PassPhIso[best_4g_idx4_m{m}]==1".format(m=mass))
            df=df.Define("best_4g_passBitMap_medium_iso_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx2_m{m}]==2 && Photon_PassPhIso[best_4g_idx2_m{m}]==2 && Photon_PassPhIso[best_4g_idx3_m{m}]==2 && Photon_PassPhIso[best_4g_idx4_m{m}]==2".format(m=mass))
            df=df.Define("best_4g_passBitMap_tight_iso_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx3_m{m}]==3 && Photon_PassPhIso[best_4g_idx2_m{m}]==3 && Photon_PassPhIso[best_4g_idx3_m{m}]==3 && Photon_PassPhIso[best_4g_idx4_m{m}]==3".format(m=mass))
            return df
        else:
            return df

        #the relative isolation values for each of the 4 best photons. i
        # ANYTHING WITHIN THE EQUAL SIGNS NOT RECOMMENDED TO USE BY EGM - ONLY FOR QUICK STUDIES/TESTS.
#==========================================================================================================================
#==========================================================================================================================
        #df=df.Define("Photon_RelIso_gamma1_m{m}".format(m=mass), "Photon_pfRelIso03_all_quadratic[best_4g_idx1_m{m}]".format(m=mass))
        #df=df.Define("Photon_RelIso_gamma2_m{m}".format(m=mass), "Photon_pfRelIso03_all_quadratic[best_4g_idx2_m{m}]".format(m=mass))
        #df=df.Define("Photon_RelIso_gamma3_m{m}".format(m=mass), "Photon_pfRelIso03_all_quadratic[best_4g_idx3_m{m}]".format(m=mass))
        #df=df.Define("Photon_RelIso_gamma4_m{m}".format(m=mass), "Photon_pfRelIso03_all_quadratic[best_4g_idx4_m{m}]".format(m=mass))

        #sanity checks to be sure the boolean logic on the isolation with the photons is correct. checking the best 4 photons individually. 
        #i check the n+p iso value, the n+p cut value, if it passes for all 4 then repeat for the charged. crucial to make sure it works.
        #df=df.Define("Photon_valIso_neutral_plus_photon_gamma1_m{m}".format(m=mass), "Photon_pfRelIso03_all_quadratic[best_4g_idx1_m{m}]-Photon_pfRelIso03_chg[best_4g_idx1_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutEGMIso_neutral_plus_photon_gamma1_m{m}".format(m=mass), "Photon_looseEGMIso_neutral_plus_photon_cut[best_4g_idx1_m{m}]".format(m=mass))
        #df=df.Define("Photon_passEGMIso_neutral_plus_photon_gamma1_m{m}".format(m=mass), "Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx1_m{m}]".format(m=mass))
        #df=df.Define("Photon_valIso_chg_gamma1_m{m}".format(m=mass), "Photon_pfRelIso03_chg[best_4g_idx1_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutEGMIso_chg_gamma1_m{m}".format(m=mass), "Photon_looseEGMIso_chg_cut[best_4g_idx1_m{m}]".format(m=mass))
        #df=df.Define("Photon_passEGMIso_chg_gamma1_m{m}".format(m=mass), "Photon_passlooseEGMIso_charged[best_4g_idx1_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutcustomIso_neutral_plus_photon_gamma1_m{m}".format(m=mass), "Photon_loosecustomIso_neutral_plus_photon_cut[best_4g_idx1_m{m}]".format(m=mass))
        #df=df.Define("Photon_passcustomIso_neutral_plus_photon_gamma1_m{m}".format(m=mass), "Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx1_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutcustomIso_chg_gamma1_m{m}".format(m=mass), "Photon_loosecustomIso_chg_cut[best_4g_idx1_m{m}]".format(m=mass))
        #df=df.Define("Photon_passcustomIso_chg_gamma1_m{m}".format(m=mass), "Photon_passloosecustomIso_charged[best_4g_idx1_m{m}]".format(m=mass))

        #df=df.Define("Photon_valIso_neutral_plus_photon_gamma2_m{m}".format(m=mass), "Photon_pfRelIso03_all_quadratic[best_4g_idx2_m{m}]-Photon_pfRelIso03_chg[best_4g_idx2_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutEGMIso_neutral_plus_photon_gamma2_m{m}".format(m=mass), "Photon_looseEGMIso_neutral_plus_photon_cut[best_4g_idx2_m{m}]".format(m=mass))
        #df=df.Define("Photon_passEGMIso_neutral_plus_photon_gamma2_m{m}".format(m=mass), "Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx2_m{m}]".format(m=mass))
        #df=df.Define("Photon_valIso_chg_gamma2_m{m}".format(m=mass), "Photon_pfRelIso03_chg[best_4g_idx2_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutEGMIso_chg_gamma2_m{m}".format(m=mass), "Photon_looseEGMIso_chg_cut[best_4g_idx2_m{m}]".format(m=mass))
        #df=df.Define("Photon_passEGMIso_chg_gamma2_m{m}".format(m=mass), "Photon_passlooseEGMIso_charged[best_4g_idx2_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutcustomIso_neutral_plus_photon_gamma2_m{m}".format(m=mass), "Photon_loosecustomIso_neutral_plus_photon_cut[best_4g_idx2_m{m}]".format(m=mass))
        #df=df.Define("Photon_passcustomIso_neutral_plus_photon_gamma2_m{m}".format(m=mass), "Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx2_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutcustomIso_chg_gamma2_m{m}".format(m=mass), "Photon_loosecustomIso_chg_cut[best_4g_idx2_m{m}]".format(m=mass))
        #df=df.Define("Photon_passcustomIso_chg_gamma2_m{m}".format(m=mass), "Photon_passloosecustomIso_charged[best_4g_idx2_m{m}]".format(m=mass))

        #df=df.Define("Photon_valIso_neutral_plus_photon_gamma3_m{m}".format(m=mass), "Photon_pfRelIso03_all_quadratic[best_4g_idx3_m{m}]-Photon_pfRelIso03_chg[best_4g_idx3_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutEGMIso_neutral_plus_photon_gamma3_m{m}".format(m=mass), "Photon_looseEGMIso_neutral_plus_photon_cut[best_4g_idx3_m{m}]".format(m=mass))
        #df=df.Define("Photon_passEGMIso_neutral_plus_photon_gamma3_m{m}".format(m=mass), "Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx3_m{m}]".format(m=mass))
        #df=df.Define("Photon_valIso_chg_gamma3_m{m}".format(m=mass), "Photon_pfRelIso03_chg[best_4g_idx3_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutEGMIso_chg_gamma3_m{m}".format(m=mass), "Photon_looseEGMIso_chg_cut[best_4g_idx3_m{m}]".format(m=mass))
        #df=df.Define("Photon_passEGMIso_chg_gamma3_m{m}".format(m=mass), "Photon_passlooseEGMIso_charged[best_4g_idx3_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutcustomIso_neutral_plus_photon_gamma3_m{m}".format(m=mass), "Photon_loosecustomIso_neutral_plus_photon_cut[best_4g_idx3_m{m}]".format(m=mass))
        #df=df.Define("Photon_passcustomIso_neutral_plus_photon_gamma3_m{m}".format(m=mass), "Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx3_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutcustomIso_chg_gamma3_m{m}".format(m=mass), "Photon_loosecustomIso_chg_cut[best_4g_idx3_m{m}]".format(m=mass))
        #df=df.Define("Photon_passcustomIso_chg_gamma3_m{m}".format(m=mass), "Photon_passloosecustomIso_charged[best_4g_idx3_m{m}]".format(m=mass))

        #df=df.Define("Photon_valIso_neutral_plus_photon_gamma4_m{m}".format(m=mass), "Photon_pfRelIso03_all_quadratic[best_4g_idx4_m{m}]-Photon_pfRelIso03_chg[best_4g_idx4_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutEGMIso_neutral_plus_photon_gamma4_m{m}".format(m=mass), "Photon_looseEGMIso_neutral_plus_photon_cut[best_4g_idx4_m{m}]".format(m=mass))
        #df=df.Define("Photon_passEGMIso_neutral_plus_photon_gamma4_m{m}".format(m=mass), "Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx4_m{m}]".format(m=mass))
        #df=df.Define("Photon_valIso_chg_gamma4_m{m}".format(m=mass), "Photon_pfRelIso03_chg[best_4g_idx4_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutEGMIso_chg_gamma4_m{m}".format(m=mass), "Photon_looseEGMIso_chg_cut[best_4g_idx4_m{m}]".format(m=mass))
        #df=df.Define("Photon_passEGMIso_chg_gamma4_m{m}".format(m=mass), "Photon_passlooseEGMIso_charged[best_4g_idx4_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutcustomIso_neutral_plus_photon_gamma4_m{m}".format(m=mass), "Photon_loosecustomIso_neutral_plus_photon_cut[best_4g_idx4_m{m}]".format(m=mass))
        #df=df.Define("Photon_passcustomIso_neutral_plus_photon_gamma4_m{m}".format(m=mass), "Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx4_m{m}]".format(m=mass))
        #df=df.Define("Photon_cutcustomIso_chg_gamma4_m{m}".format(m=mass), "Photon_loosecustomIso_chg_cut[best_4g_idx4_m{m}]".format(m=mass))
        #df=df.Define("Photon_passcustomIso_chg_gamma4_m{m}".format(m=mass), "Photon_passloosecustomIso_charged[best_4g_idx4_m{m}]".format(m=mass))
       
        #here are event-level booleans to see if all of the best 4 photons pass the threshold cuts for min(neutral, photon) and charged iso cuts. see chelpers.h for functions
        #df=df.Define("best_4g_passLooseEGMIso_neutral_plus_photon_m{m}".format(m=mass), "(Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx1_m{m}]==1)&&(Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx2_m{m}]==1)&&(Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx3_m{m}]==1)&&(Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx4_m{m}]==1)".format(m=mass))
        #df=df.Define("best_4g_passLooseEGMIso_charged_m{m}".format(m=mass), "(Photon_passlooseEGMIso_charged[best_4g_idx1_m{m}]==1)&&(Photon_passlooseEGMIso_charged[best_4g_idx2_m{m}]==1)&&(Photon_passlooseEGMIso_charged[best_4g_idx3_m{m}]==1)&&(Photon_passlooseEGMIso_charged[best_4g_idx4_m{m}]==1)".format(m=mass))
        #df=df.Define("best_4g_passLooseEGMIso_m{m}".format(m=mass), "best_4g_passLooseEGMIso_neutral_plus_photon_m{m} && best_4g_passLooseEGMIso_charged_m{m}".format(m=mass))

        #df=df.Define("best_4g_passcustomIso_neutral_plus_photon_m{m}".format(m=mass), "(Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx1_m{m}]==1)&&(Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx2_m{m}]==1)&&(Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx3_m{m}]==1)&&(Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx4_m{m}]==1)".format(m=mass))
        #df=df.Define("best_4g_passcustomIso_charged_m{m}".format(m=mass), "(Photon_passloosecustomIso_charged[best_4g_idx1_m{m}]==1)&&(Photon_passloosecustomIso_charged[best_4g_idx2_m{m}]==1)&&(Photon_passloosecustomIso_charged[best_4g_idx3_m{m}]==1)&&(Photon_passloosecustomIso_charged[best_4g_idx4_m{m}]==1)".format(m=mass))
        #df=df.Define("best_4g_passLoosecustomIso_m{m}".format(m=mass), "best_4g_passcustomIso_neutral_plus_photon_m{m} && best_4g_passcustomIso_charged_m{m}".format(m=mass))
#==========================================================================================================================
#==========================================================================================================================

    def corrected_kinematic_vars(df, mass):
        #pt's
        df=df.Define('best_4g_phi1_gamma1_pt_m{}'.format(mass),'raw_best_4g_m{}[0]'.format(mass))
        df=df.Define('best_4g_phi1_gamma2_pt_m{}'.format(mass),'raw_best_4g_m{}[3]'.format(mass))
        df=df.Define('best_4g_phi2_gamma1_pt_m{}'.format(mass),'raw_best_4g_m{}[8]'.format(mass))
        df=df.Define('best_4g_phi2_gamma2_pt_m{}'.format(mass),'raw_best_4g_m{}[11]'.format(mass))

        #eta's
        df=df.Define('best_4g_phi1_gamma1_eta_m{}'.format(mass),'raw_best_4g_m{}[1]'.format(mass))
        df=df.Define('best_4g_phi1_gamma2_eta_m{}'.format(mass),'raw_best_4g_m{}[4]'.format(mass))
        df=df.Define('best_4g_phi2_gamma1_eta_m{}'.format(mass),'raw_best_4g_m{}[9]'.format(mass))
        df=df.Define('best_4g_phi2_gamma2_eta_m{}'.format(mass),'raw_best_4g_m{}[12]'.format(mass))

        #phi's
        df=df.Define('best_4g_phi1_gamma1_phi_m{}'.format(mass),'raw_best_4g_m{}[2]'.format(mass))
        df=df.Define('best_4g_phi1_gamma2_phi_m{}'.format(mass),'raw_best_4g_m{}[5]'.format(mass))
        df=df.Define('best_4g_phi2_gamma1_phi_m{}'.format(mass),'raw_best_4g_m{}[10]'.format(mass))
        df=df.Define('best_4g_phi2_gamma2_phi_m{}'.format(mass),'raw_best_4g_m{}[13]'.format(mass))

        #id's
        df=df.Define('best_4g_phi1_gamma1_id_m{}'.format(mass),'raw_best_4g_m{}[20]'.format(mass))
        df=df.Define('best_4g_phi1_gamma2_id_m{}'.format(mass),'raw_best_4g_m{}[21]'.format(mass))
        df=df.Define('best_4g_phi2_gamma1_id_m{}'.format(mass),'raw_best_4g_m{}[22]'.format(mass))
        df=df.Define('best_4g_phi2_gamma2_id_m{}'.format(mass),'raw_best_4g_m{}[23]'.format(mass))
        
        #lxy's
        df=df.Define('best_4g_phi1_dxy_m{}'.format(mass),'raw_best_4g_m{}[6]'.format(mass))
        df=df.Define('best_4g_phi2_dxy_m{}'.format(mass),'raw_best_4g_m{}[14]'.format(mass))
        return df
    
    def uncorrected_kinematic_vars(df, mass, scouting):
        if scouting==0:
            df=df.Define('Photon_pt_gamma1_m{m}'.format(m=mass),'Photon_pt[best_4g_idx1_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_pt_gamma2_m{m}'.format(m=mass),'Photon_pt[best_4g_idx2_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_pt_gamma3_m{m}'.format(m=mass),'Photon_pt[best_4g_idx3_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_pt_gamma4_m{m}'.format(m=mass),'Photon_pt[best_4g_idx4_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_hoe_gamma1_m{m}'.format(m=mass),'Photon_hoe[best_4g_idx1_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_hoe_gamma2_m{m}'.format(m=mass),'Photon_hoe[best_4g_idx2_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_hoe_gamma3_m{m}'.format(m=mass),'Photon_hoe[best_4g_idx3_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_hoe_gamma4_m{m}'.format(m=mass),'Photon_hoe[best_4g_idx4_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_sieie_gamma1_m{m}'.format(m=mass),'Photon_sieie[best_4g_idx1_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_sieie_gamma2_m{m}'.format(m=mass),'Photon_sieie[best_4g_idx2_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_sieie_gamma3_m{m}'.format(m=mass),'Photon_sieie[best_4g_idx3_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_sieie_gamma4_m{m}'.format(m=mass),'Photon_sieie[best_4g_idx4_m{m}]'.format(m=mass)) 
        else:
            df=df.Define('Photon_pt_gamma1_m{m}'.format(m=mass),'ScoutingPhoton_pt[best_4g_idx1_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_pt_gamma2_m{m}'.format(m=mass),'ScoutingPhoton_pt[best_4g_idx2_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_pt_gamma3_m{m}'.format(m=mass),'ScoutingPhoton_pt[best_4g_idx3_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_pt_gamma4_m{m}'.format(m=mass),'ScoutingPhoton_pt[best_4g_idx4_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_hoe_gamma1_m{m}'.format(m=mass),'ScoutingPhoton_hOverE[best_4g_idx1_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_hoe_gamma2_m{m}'.format(m=mass),'ScoutingPhoton_hOverE[best_4g_idx2_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_hoe_gamma3_m{m}'.format(m=mass),'ScoutingPhoton_hOverE[best_4g_idx3_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_hoe_gamma4_m{m}'.format(m=mass),'ScoutingPhoton_hOverE[best_4g_idx4_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_sieie_gamma1_m{m}'.format(m=mass),'ScoutingPhoton_sigmaIetaIeta[best_4g_idx1_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_sieie_gamma2_m{m}'.format(m=mass),'ScoutingPhoton_sigmaIetaIeta[best_4g_idx2_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_sieie_gamma3_m{m}'.format(m=mass),'ScoutingPhoton_sigmaIetaIeta[best_4g_idx3_m{m}]'.format(m=mass)) 
            df=df.Define('Photon_sieie_gamma4_m{m}'.format(m=mass),'ScoutingPhoton_sigmaIetaIeta[best_4g_idx4_m{m}]'.format(m=mass)) 
        return df

    def indices(df, mass):
        df=df.Define('best_4g_idx1_m{}'.format(mass),'raw_best_4g_m{}[24]'.format(mass)) 
        df=df.Define('best_4g_idx2_m{}'.format(mass),'raw_best_4g_m{}[25]'.format(mass)) 
        df=df.Define('best_4g_idx3_m{}'.format(mass),'raw_best_4g_m{}[26]'.format(mass)) 
        df=df.Define('best_4g_idx4_m{}'.format(mass),'raw_best_4g_m{}[27]'.format(mass)) 
        return df

    def masses(df, mass):
        df=df.Define('best_4g_phi1_mass_m{}'.format(mass),'raw_best_4g_m{}[16]'.format(mass))
        df=df.Define('best_4g_phi2_mass_m{}'.format(mass),'raw_best_4g_m{}[17]'.format(mass))
        df=df.Define('best_4g_uncorr_mass_m{}'.format(mass),'raw_best_4g_m{}[18]'.format(mass))
        df=df.Define('best_4g_corr_mass_m{}'.format(mass),'raw_best_4g_m{}[19]'.format(mass))
        return df
    
    def preselection(df, mass):
        df=df.Define('best_4g_preselected_m{m}'.format(m=mass),'Photon_preselection[best_4g_idx1_m{m}] && Photon_preselection[best_4g_idx2_m{m}] && Photon_preselection[best_4g_idx3_m{m}] && Photon_preselection[best_4g_idx4_m{m}]'.format(m=mass)) 
        return df
    
    def IDs(df, mass, scouting):
        df=df.Define('best_4g_sumID_m{}'.format(mass),'raw_best_4g_m{m}[20]+raw_best_4g_m{m}[21]+raw_best_4g_m{m}[22]+raw_best_4g_m{m}[23]'.format(m=mass))
        if scouting==0:
            preselection_str = " && ".join(f"(Photon_preselection[raw_best_4g_m{{m}}[{i}]]==1)" for i in range(24,28))
            idnoiso_str = " && ".join(f"(Photon_IdNoIso[raw_best_4g_m{{m}}[{i}]]==1)" for i in range(24,28))
            iso_str = " && ".join(f"(Photon_PassPhIso[raw_best_4g_m{{m}}[{i}]]==1)" for i in range(24,28))
            best_4g_ID_str_ggH4g="{} && {} && {}".format(preselection_str, idnoiso_str, iso_str)
            df=df.Define('best_4g_ID_m{}'.format(mass),best_4g_ID_str_ggH4g.format(m=mass))
            return df
        else:
            preselection_str = " && ".join(f"(Photon_preselection[raw_best_4g_m{{m}}[{i}]]==1)" for i in range(24,28))
            idnoiso_str = " && ".join(f"(Photon_IdNoIso[raw_best_4g_m{{m}}[{i}]]==1)" for i in range(24,28))
            best_4g_ID_str_ggH4g="{} && {}".format(preselection_str, idnoiso_str)
            df=df.Define('best_4g_ID_m{}'.format(mass),best_4g_ID_str_ggH4g.format(m=mass))
            return df
    

    for mass in phi_mass:
        #Photon_corrIso_m{} needs to be defined first since four_gamma relies on it
        if scouting==0:
            ggH4g=ggH4g.Define(f'Photon_corrIso_m{mass}','correct_gammaIso(Photon_pt,Photon_eta,Photon_phi,{},Photon_preselection)'.format(iso[data['era']]))
        else:
            ggH4g=ggH4g.Define(f'Photon_corrIso_m{mass}','correct_gammaIso(ScoutingPhoton_pt,ScoutingPhoton_eta,ScoutingPhoton_phi,{},Photon_preselection)'.format(iso_scouting[data['era']]))

        #organizing the branch definitions into functions that couple related branches
        ggH4g=four_gamma(ggH4g, mass, scouting)
        ggH4g=indices(ggH4g, mass)
        ggH4g=corrected_kinematic_vars(ggH4g, mass)
        ggH4g=uncorrected_kinematic_vars(ggH4g, mass, scouting)
        ggH4g=preselection(ggH4g, mass)
        ggH4g=isolation_vars(ggH4g, mass, scouting)
        ggH4g=masses(ggH4g, mass)
        ggH4g=IDs(ggH4g, mass, scouting)

        ggH4g=ggH4g.Define('non_MC_cut_m{}'.format(mass),'sample_isMC==0 && best_4g_uncorr_mass_m{m}<90|best_4g_uncorr_mass_m{m}>150'.format(m=mass))
        ggH4g=ggH4g.Define('best_4g_phi1_valid_m{}'.format(mass),'raw_best_4g_m{}[7]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_valid_m{}'.format(mass),'raw_best_4g_m{}[15]'.format(mass))


    #ggH4g=ggH4g.Filter('sample_isMC==1 | non_MC_cut_m{}==1','blinding_data_samples')

    actions.append(ggH4g.Snapshot('ggH4g', f"{sample}_ggH4g.root", cols, opts))
    save_report(ggH4g, "Report_ggH4g", f"{sample}_ggH4g", opts, actions)
    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, f"{sample}_ggH4g.root", "", opts))

    return actions

def Zee(data, sample):
    cols = "best_.*|sample_.*|^Photon_.*|^Electron_.*|Weight.*|^Gen.*|^weight.*|^TrigObj_.*|^event.*|^Pileup_.*|^run.*|gen.*|.*LHE.*|^PV.*|luminosity|Block|genWeight|HLT_passed|sorted_photon_pt|DST_PFScouting_DoubleEG|nScoutingPhoton|^ScoutingPhoton_.*|gg_scouting|^gg_scouting_.*|good_idx|pt_id|eta_id|phi_id|pho_pass_id|sorted_pt_id|nPhoID|pass_DST"

    ROOT.gInterpreter.Declare("""
    #include <ROOT/RVec.hxx>
    #include <Math/Vector4D.h>
    #include <algorithm>
    using ROOT::RVec;
    RVec<float> gg_info(const RVec<float>& pt,const RVec<float>& eta,const RVec<float>& phi) {
        RVec<float> out(3, -1.0);
        int n = pt.size();
        if (n < 2) return out;
        RVec<int> idx(n);
        
        for (int i=0; i<n; i++) idx[i]=i;
        std::sort(idx.begin(), idx.end(), [&](int a, int b){return pt[a]>pt[b];});
        int i1=idx[0];
        int i2=idx[1];
        ROOT::Math::PtEtaPhiMVector v1(pt[i1],eta[i1],phi[i1],0.0);
        ROOT::Math::PtEtaPhiMVector v2(pt[i2],eta[i2],phi[i2],0.0);
        float mass=(v1+v2).M();
        out[0]=mass;
        out[1]=float(i1);
        out[2]=float(i2);
        return out;
    }
    """)

    actions=[]

    dataframe =load_meta_data(data)
    Zee=dataframe["Events"]

    Zee=Zee.Define("pho_pass_id","""
    (ROOT::VecOps::RVec<char>)((
    (abs(ScoutingPhoton_eta)<1.3 &&
     ScoutingPhoton_sigmaIetaIeta<0.04 &&
     ScoutingPhoton_hOverE<0.2 &&
     ScoutingPhoton_ecalIso<40.0)
    ||
    (abs(ScoutingPhoton_eta)>1.3 && abs(ScoutingPhoton_eta)<2.5 &&
     ScoutingPhoton_sigmaIetaIeta<0.08 &&
     ScoutingPhoton_hOverE<0.14 &&
     ScoutingPhoton_ecalIso<40.0)))
    """)

    if data["isMC"]:
        Zee=Zee.Define("Pileup_weight", "getPUweight(Pileup_nPU, puWeight_UL{}, sample_isMC)".format(data["era"]))

    Zee=Zee.Define("nPhoID", "Sum(pho_pass_id)")
    Zee=Zee.Define("pass_DST", "DST_PFScouting_DoubleEG == 1")

    Zee=Zee.Define("good_idx", "Nonzero(pho_pass_id)")
    Zee=Zee.Define("pt_id",  "Take(ScoutingPhoton_pt,  good_idx)")
    Zee=Zee.Define("sorted_pt_id", "Reverse(Sort(pt_id))")
    Zee=Zee.Define("eta_id", "Take(ScoutingPhoton_eta, good_idx)")
    Zee=Zee.Define("phi_id", "Take(ScoutingPhoton_phi, good_idx)")
    Zee=Zee.Define("gg_scouting", "gg_info(pt_id, eta_id, phi_id)")
    Zee=Zee.Define("gg_scouting_mass", "gg_scouting[0]")

    actions.append(Zee.Snapshot('Zee', f"{sample}_Zee.root", cols, opts))
    save_report(Zee, "Report_Zee", f"{sample}_Zee", opts, actions)
    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, f"{sample}_Zee.root", "", opts))

    return actions 

def analysis(data,sample):
    actions=[]
    phi_mass=[7,15,20,30,40,50,55]
    actions.extend(ggH(data,phi_mass,sample))
    #actions.extend(Zee(data, sample))
    return actions
