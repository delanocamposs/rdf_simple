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
cols = "best_.*|sample_.*|^Photon_.*|Weight.*|^Gen.*|^weight.*|^TrigObj_.*|^event.*|^Pileup_.*|^run.*|gen.*|.*LHE.*|^PV.*|luminosity|Block|genWeight|HLT_passed"

iso = {'2016': 'Photon_pfRelIso03_all',
       '2017': 'Photon_pfRelIso03_all',
       '2018': 'Photon_pfRelIso03_all',
       '2022': 'Photon_pfRelIso03_all_quadratic',
       '2023': 'Photon_pfRelIso03_all_quadratic',
       '2024': 'Photon_pfRelIso03_all_quadratic'
       }

def get_ID_val(var):
    ID_type="Custom"
    cuts_ID={
            "EGM":{"endcap":{"iso":[0.1], 
                             "hoe":[0.0590], 
                             "sieie":[0.0272]}, 
                   "barrel":{"iso":[0.1], 
                             "hoe":[0.04596], 
                             "sieie":[0.0106]}},
            "Custom":{"endcap":{"iso":[0.3], 
                                "hoe":[0.2], 
                                "sieie":[0.045]}, 
                      "barrel":{"iso":[0.25], 
                                "hoe":[0.3], 
                                "sieie":[0.035]}}
            }

    return cuts_ID[ID_type]['barrel'][var][0], cuts_ID[ID_type]['endcap'][var][0]

ID_dict1_ggH4g={
 'Photon_preselection':'==1'}
ID_dict2_ggH4g={
 'Photon_IdNoIso':'==1'}

iso1, iso2 =get_ID_val("iso")
ID1_ggH4g = " && ".join(f"{branch}[raw_best_4g_m{{m}}[{i}]]{b}" for branch, b in ID_dict1_ggH4g.items() for i in range(24,28))
ID2_ggH4g = " && ".join(f"{branch}[raw_best_4g_m{{m}}[{i}]]{b}" for branch, b in ID_dict2_ggH4g.items() for i in range(24,28))
ID3_ggH4g = " && ".join(f"((Photon_isScEtaEB[raw_best_4g_m{{m}}[{i}]]==1 && Photon_corrIso_m{{m}}[raw_best_4g_m{{m}}[{i}]]<{iso1})||(Photon_isScEtaEE[raw_best_4g_m{{m}}[{i}]]==1 && Photon_corrIso_m{{m}}[raw_best_4g_m{{m}}[{i}]]<{iso2}))" for i in range(24,28))
#ID3_ggH4g = " && ".join(f"(Photon_PassPhIso[raw_best_4g_m{{m}}[{i}]]==1)" for i in range(24,28))

#only preselection and ID rules change for ggH4g_1bad. EE and EB logic identical
#ID1_1bad = "(((Photon_preselection[raw_best_4g_m{m}[24]])+(Photon_preselection[raw_best_4g_m{m}[25]])+(Photon_preselection[raw_best_4g_m{m}[26]])+(Photon_preselection[raw_best_4g_m{m}[27]]))==3)"
#ID2_1bad = "(((Photon_IdNoIso[raw_best_4g_m{m}[24]])+(Photon_IdNoIso[raw_best_4g_m{m}[25]])+(Photon_IdNoIso[raw_best_4g_m{m}[26]])+(Photon_IdNoIso[raw_best_4g_m{m}[27]]))==3)"

#strings to pass into Define in phi mass loop
#best_4g_ID_str_ggH4g_1bad = "{} && {} && {}".format(ID1_1bad, ID2_1bad, ID3_ggH4g)
best_4g_ID_str_ggH4g="{} && {}".format(ID1_ggH4g, ID2_ggH4g)
    
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
    
def photonAna(dataframe):

    # Overlap with loose leptons
    photons = dataframe.Define("Photon_muOverlap", "overlapClean(Photon_phi, Photon_eta, Muon_phi[loose_muon], Muon_eta[loose_muon])")
    photons = photons.Define("Photon_eleOverlap", "overlapClean(Photon_phi, Photon_eta, Electron_phi[loose_electron], Electron_eta[loose_electron])")
    photons = photons.Define("Photon_overlap", "Photon_muOverlap||Photon_eleOverlap")

    # Photon Preselection criteria
    photons = photons.Define("Photon_preselection", "Photon_pt>20&&!Photon_pixelSeed&&abs(Photon_eta)<2.5&&(abs(Photon_eta)>1.57||abs(Photon_eta)<1.44)&&(Photon_isScEtaEE||Photon_isScEtaEB)")
    
    photons = photons.Define("Photon_rho", "fixedGridRhoFastjetAll")
    photons=photons.Define("Photon_PassPhIso" , "passPhIso(Photon_vidNestedWPBitmap)")

    #etxra information about isolation. EGM recommends to only use for quick studies/not part of the official ID. use bitmap instead for the ID (above).
    photons=photons.Define("Photon_passlooseEGMIso_neutral_plus_photon", "Photon_corrRelIso_EGM_loose_np(Photon_pt, Photon_isScEtaEE, Photon_isScEtaEB, Photon_pfRelIso03_all-Photon_pfRelIso03_chg)")  #bools for neutral+photon
    photons=photons.Define("Photon_passlooseEGMIso_charged", "Photon_corrRelIso_EGM_loose_chg(Photon_isScEtaEE, Photon_isScEtaEB, Photon_pfRelIso03_chg)") #bools for charged
    photons=photons.Define("Photon_looseEGMIso_neutral_plus_photon_cut", "Photon_corrRelIso_EGM_loose_cut_np(Photon_pt, Photon_isScEtaEE, Photon_isScEtaEB)")
    photons=photons.Define("Photon_looseEGMIso_chg_cut", "Photon_corrRelIso_EGM_loose_cut_chg(Photon_isScEtaEE, Photon_isScEtaEB, Photon_pfRelIso03_chg)")

    photons=photons.Define("Photon_passloosecustomIso_neutral_plus_photon", "Photon_corrRelIso_custom_loose_np(Photon_pt, Photon_isScEtaEE, Photon_isScEtaEB, Photon_pfRelIso03_all-Photon_pfRelIso03_chg)") #returns vec of bools if photon passes custom neutral+photon iso custom
    photons=photons.Define("Photon_passloosecustomIso_charged", "Photon_corrRelIso_custom_loose_chg(Photon_isScEtaEE, Photon_isScEtaEB, Photon_pfRelIso03_chg)") #bools for charged custom
    photons=photons.Define("Photon_loosecustomIso_neutral_plus_photon_cut", "Photon_corrRelIso_custom_loose_cut_np(Photon_pt, Photon_isScEtaEE, Photon_isScEtaEB)")
    photons=photons.Define("Photon_loosecustomIso_chg_cut", "Photon_corrRelIso_custom_loose_cut_chg(Photon_isScEtaEE, Photon_isScEtaEB, Photon_pfRelIso03_chg)")

    # Common Photon ID definitions (No isolation)
    sieie1, sieie2 = get_ID_val("sieie")
    hoe1, hoe2 = get_ID_val("hoe")
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
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    ggH=dataframe["Events"]


    if data["isMC"]:
        ggH = ggH.Define("Pileup_weight", "getPUweight(Pileup_nPU, puWeight_UL{}, sample_isMC)".format(data["era"]))

    

    #Filter out muons above 10Gev and electrons above 15GeV
    ggH=electronAna(ggH)
    ggH=muonAna(ggH)
    
    ggH=ggH.Filter("Sum(loose_muon==1)==0",'muon_veto')
    ggH=ggH.Filter("Sum(loose_electron==1)==0",'electron_veto')

    #Add photon preselection and common photon ID definitions 
    ggH=photonAna(ggH)
    
    #filtering for all snapshots done here:
    ggH4g=ggH.Filter('nPhoton>3','at_least_4_photons')
    ggH4g=ggH4g.Filter('Sum(Photon_preselection==1)>3','at_least_3_preselected_photons')

    def four_gamma(df, mass):
        df=df.Define('raw_best_4g_m{}'.format(mass),"best_4gamma(Photon_pt,Photon_eta,Photon_phi,Photon_isScEtaEB,Photon_isScEtaEE,Photon_preselection,Photon_IdNoIso,Photon_corrIso_m{},{})".format(mass,float(mass)))
        return df

    def isolation_vars(df, mass):
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
 
        df=df.Define("best_4g_passBitMap_loose_iso_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx1_m{m}]==1 && Photon_PassPhIso[best_4g_idx2_m{m}]==1 && Photon_PassPhIso[best_4g_idx3_m{m}]==1 && Photon_PassPhIso[best_4g_idx4_m{m}]==1".format(m=mass))
        df=df.Define("best_4g_passBitMap_medium_iso_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx2_m{m}]==2 && Photon_PassPhIso[best_4g_idx2_m{m}]==2 && Photon_PassPhIso[best_4g_idx3_m{m}]==2 && Photon_PassPhIso[best_4g_idx4_m{m}]==2".format(m=mass))
        df=df.Define("best_4g_passBitMap_tight_iso_m{m}".format(m=mass), "Photon_PassPhIso[best_4g_idx3_m{m}]==3 && Photon_PassPhIso[best_4g_idx2_m{m}]==3 && Photon_PassPhIso[best_4g_idx3_m{m}]==3 && Photon_PassPhIso[best_4g_idx4_m{m}]==3".format(m=mass))

        #the relative isolation values for each of the 4 best photons. NOT RECOMMENDED TO USE BY EGM - ONLY FOR QUICK STUDIES/TESTS.
        df=df.Define("Photon_RelIso_gamma1_m{m}".format(m=mass), "Photon_pfRelIso03_all[best_4g_idx1_m{m}]".format(m=mass))
        df=df.Define("Photon_RelIso_gamma2_m{m}".format(m=mass), "Photon_pfRelIso03_all[best_4g_idx2_m{m}]".format(m=mass))
        df=df.Define("Photon_RelIso_gamma3_m{m}".format(m=mass), "Photon_pfRelIso03_all[best_4g_idx3_m{m}]".format(m=mass))
        df=df.Define("Photon_RelIso_gamma4_m{m}".format(m=mass), "Photon_pfRelIso03_all[best_4g_idx4_m{m}]".format(m=mass))

        #sanity checks to be sure the boolean logic on the isolation with the photons is correct. checking the best 4 photons individually. 
        #i check the n+p iso value, the n+p cut value, if it passes for all 4 then repeat for the charged. alot of work but crucial to make sure it works.
        df=df.Define("Photon_valIso_neutral_plus_photon_gamma1_m{m}".format(m=mass), "Photon_pfRelIso03_all[best_4g_idx1_m{m}]-Photon_pfRelIso03_chg[best_4g_idx1_m{m}]".format(m=mass))
        df=df.Define("Photon_cutEGMIso_neutral_plus_photon_gamma1_m{m}".format(m=mass), "Photon_looseEGMIso_neutral_plus_photon_cut[best_4g_idx1_m{m}]".format(m=mass))
        df=df.Define("Photon_passEGMIso_neutral_plus_photon_gamma1_m{m}".format(m=mass), "Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx1_m{m}]".format(m=mass))
        df=df.Define("Photon_valIso_chg_gamma1_m{m}".format(m=mass), "Photon_pfRelIso03_chg[best_4g_idx1_m{m}]".format(m=mass))
        df=df.Define("Photon_cutEGMIso_chg_gamma1_m{m}".format(m=mass), "Photon_looseEGMIso_chg_cut[best_4g_idx1_m{m}]".format(m=mass))
        df=df.Define("Photon_passEGMIso_chg_gamma1_m{m}".format(m=mass), "Photon_passlooseEGMIso_charged[best_4g_idx1_m{m}]".format(m=mass))
        df=df.Define("Photon_cutcustomIso_neutral_plus_photon_gamma1_m{m}".format(m=mass), "Photon_loosecustomIso_neutral_plus_photon_cut[best_4g_idx1_m{m}]".format(m=mass))
        df=df.Define("Photon_passcustomIso_neutral_plus_photon_gamma1_m{m}".format(m=mass), "Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx1_m{m}]".format(m=mass))
        df=df.Define("Photon_cutcustomIso_chg_gamma1_m{m}".format(m=mass), "Photon_loosecustomIso_chg_cut[best_4g_idx1_m{m}]".format(m=mass))
        df=df.Define("Photon_passcustomIso_chg_gamma1_m{m}".format(m=mass), "Photon_passloosecustomIso_charged[best_4g_idx1_m{m}]".format(m=mass))

        df=df.Define("Photon_valIso_neutral_plus_photon_gamma2_m{m}".format(m=mass), "Photon_pfRelIso03_all[best_4g_idx2_m{m}]-Photon_pfRelIso03_chg[best_4g_idx2_m{m}]".format(m=mass))
        df=df.Define("Photon_cutEGMIso_neutral_plus_photon_gamma2_m{m}".format(m=mass), "Photon_looseEGMIso_neutral_plus_photon_cut[best_4g_idx2_m{m}]".format(m=mass))
        df=df.Define("Photon_passEGMIso_neutral_plus_photon_gamma2_m{m}".format(m=mass), "Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx2_m{m}]".format(m=mass))
        df=df.Define("Photon_valIso_chg_gamma2_m{m}".format(m=mass), "Photon_pfRelIso03_chg[best_4g_idx2_m{m}]".format(m=mass))
        df=df.Define("Photon_cutEGMIso_chg_gamma2_m{m}".format(m=mass), "Photon_looseEGMIso_chg_cut[best_4g_idx2_m{m}]".format(m=mass))
        df=df.Define("Photon_passEGMIso_chg_gamma2_m{m}".format(m=mass), "Photon_passlooseEGMIso_charged[best_4g_idx2_m{m}]".format(m=mass))
        df=df.Define("Photon_cutcustomIso_neutral_plus_photon_gamma2_m{m}".format(m=mass), "Photon_loosecustomIso_neutral_plus_photon_cut[best_4g_idx2_m{m}]".format(m=mass))
        df=df.Define("Photon_passcustomIso_neutral_plus_photon_gamma2_m{m}".format(m=mass), "Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx2_m{m}]".format(m=mass))
        df=df.Define("Photon_cutcustomIso_chg_gamma2_m{m}".format(m=mass), "Photon_loosecustomIso_chg_cut[best_4g_idx2_m{m}]".format(m=mass))
        df=df.Define("Photon_passcustomIso_chg_gamma2_m{m}".format(m=mass), "Photon_passloosecustomIso_charged[best_4g_idx2_m{m}]".format(m=mass))

        df=df.Define("Photon_valIso_neutral_plus_photon_gamma3_m{m}".format(m=mass), "Photon_pfRelIso03_all[best_4g_idx3_m{m}]-Photon_pfRelIso03_chg[best_4g_idx3_m{m}]".format(m=mass))
        df=df.Define("Photon_cutEGMIso_neutral_plus_photon_gamma3_m{m}".format(m=mass), "Photon_looseEGMIso_neutral_plus_photon_cut[best_4g_idx3_m{m}]".format(m=mass))
        df=df.Define("Photon_passEGMIso_neutral_plus_photon_gamma3_m{m}".format(m=mass), "Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx3_m{m}]".format(m=mass))
        df=df.Define("Photon_valIso_chg_gamma3_m{m}".format(m=mass), "Photon_pfRelIso03_chg[best_4g_idx3_m{m}]".format(m=mass))
        df=df.Define("Photon_cutEGMIso_chg_gamma3_m{m}".format(m=mass), "Photon_looseEGMIso_chg_cut[best_4g_idx3_m{m}]".format(m=mass))
        df=df.Define("Photon_passEGMIso_chg_gamma3_m{m}".format(m=mass), "Photon_passlooseEGMIso_charged[best_4g_idx3_m{m}]".format(m=mass))
        df=df.Define("Photon_cutcustomIso_neutral_plus_photon_gamma3_m{m}".format(m=mass), "Photon_loosecustomIso_neutral_plus_photon_cut[best_4g_idx3_m{m}]".format(m=mass))
        df=df.Define("Photon_passcustomIso_neutral_plus_photon_gamma3_m{m}".format(m=mass), "Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx3_m{m}]".format(m=mass))
        df=df.Define("Photon_cutcustomIso_chg_gamma3_m{m}".format(m=mass), "Photon_loosecustomIso_chg_cut[best_4g_idx3_m{m}]".format(m=mass))
        df=df.Define("Photon_passcustomIso_chg_gamma3_m{m}".format(m=mass), "Photon_passloosecustomIso_charged[best_4g_idx3_m{m}]".format(m=mass))

        df=df.Define("Photon_valIso_neutral_plus_photon_gamma4_m{m}".format(m=mass), "Photon_pfRelIso03_all[best_4g_idx4_m{m}]-Photon_pfRelIso03_chg[best_4g_idx4_m{m}]".format(m=mass))
        df=df.Define("Photon_cutEGMIso_neutral_plus_photon_gamma4_m{m}".format(m=mass), "Photon_looseEGMIso_neutral_plus_photon_cut[best_4g_idx4_m{m}]".format(m=mass))
        df=df.Define("Photon_passEGMIso_neutral_plus_photon_gamma4_m{m}".format(m=mass), "Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx4_m{m}]".format(m=mass))
        df=df.Define("Photon_valIso_chg_gamma4_m{m}".format(m=mass), "Photon_pfRelIso03_chg[best_4g_idx4_m{m}]".format(m=mass))
        df=df.Define("Photon_cutEGMIso_chg_gamma4_m{m}".format(m=mass), "Photon_looseEGMIso_chg_cut[best_4g_idx4_m{m}]".format(m=mass))
        df=df.Define("Photon_passEGMIso_chg_gamma4_m{m}".format(m=mass), "Photon_passlooseEGMIso_charged[best_4g_idx4_m{m}]".format(m=mass))
        df=df.Define("Photon_cutcustomIso_neutral_plus_photon_gamma4_m{m}".format(m=mass), "Photon_loosecustomIso_neutral_plus_photon_cut[best_4g_idx4_m{m}]".format(m=mass))
        df=df.Define("Photon_passcustomIso_neutral_plus_photon_gamma4_m{m}".format(m=mass), "Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx4_m{m}]".format(m=mass))
        df=df.Define("Photon_cutcustomIso_chg_gamma4_m{m}".format(m=mass), "Photon_loosecustomIso_chg_cut[best_4g_idx4_m{m}]".format(m=mass))
        df=df.Define("Photon_passcustomIso_chg_gamma4_m{m}".format(m=mass), "Photon_passloosecustomIso_charged[best_4g_idx4_m{m}]".format(m=mass))
       
        #here are event-level booleans to see if all of the best 4 photons pass the threshold cuts for min(neutral, photon) and charged iso cuts. see chelpers.h for functions
        df=df.Define("best_4g_passLooseEGMIso_neutral_plus_photon_m{m}".format(m=mass), "(Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx1_m{m}]==1)&&(Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx2_m{m}]==1)&&(Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx3_m{m}]==1)&&(Photon_passlooseEGMIso_neutral_plus_photon[best_4g_idx4_m{m}]==1)".format(m=mass))
        df=df.Define("best_4g_passLooseEGMIso_charged_m{m}".format(m=mass), "(Photon_passlooseEGMIso_charged[best_4g_idx1_m{m}]==1)&&(Photon_passlooseEGMIso_charged[best_4g_idx2_m{m}]==1)&&(Photon_passlooseEGMIso_charged[best_4g_idx3_m{m}]==1)&&(Photon_passlooseEGMIso_charged[best_4g_idx4_m{m}]==1)".format(m=mass))
        df=df.Define("best_4g_passLooseEGMIso_m{m}".format(m=mass), "best_4g_passLooseEGMIso_neutral_plus_photon_m{m} && best_4g_passLooseEGMIso_charged_m{m}".format(m=mass))

        df=df.Define("best_4g_passcustomIso_neutral_plus_photon_m{m}".format(m=mass), "(Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx1_m{m}]==1)&&(Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx2_m{m}]==1)&&(Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx3_m{m}]==1)&&(Photon_passloosecustomIso_neutral_plus_photon[best_4g_idx4_m{m}]==1)".format(m=mass))
        df=df.Define("best_4g_passcustomIso_charged_m{m}".format(m=mass), "(Photon_passloosecustomIso_charged[best_4g_idx1_m{m}]==1)&&(Photon_passloosecustomIso_charged[best_4g_idx2_m{m}]==1)&&(Photon_passloosecustomIso_charged[best_4g_idx3_m{m}]==1)&&(Photon_passloosecustomIso_charged[best_4g_idx4_m{m}]==1)".format(m=mass))
        df=df.Define("best_4g_passLoosecustomIso_m{m}".format(m=mass), "best_4g_passcustomIso_neutral_plus_photon_m{m} && best_4g_passcustomIso_charged_m{m}".format(m=mass))
        return df

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
    
    def uncorrected_kinematic_vars(df, mass):
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
    
    def IDs(df, mass):
        df=df.Define('best_4g_sumID_m{}'.format(mass),'raw_best_4g_m{m}[20]+raw_best_4g_m{m}[21]+raw_best_4g_m{m}[22]+raw_best_4g_m{m}[23]'.format(m=mass))
        df=df.Define('best_4g_ID_m{}'.format(mass),best_4g_ID_str_ggH4g.format(m=mass))
        return df
        
        

    for mass in phi_mass:
        #Photon_corrIso_m{} needs to be defined first since four_gamma relies on it
        ggH4g=ggH4g.Define('Photon_corrIso_m{}'.format(mass),'correct_gammaIso(Photon_pt,Photon_eta,Photon_phi,{},Photon_preselection)'.format(iso[data['era']]))

        #organizing the branch definitions into functions that couple related branches
        ggH4g=four_gamma(ggH4g, mass)
        ggH4g=indices(ggH4g, mass)
        ggH4g=corrected_kinematic_vars(ggH4g, mass)
        ggH4g=uncorrected_kinematic_vars(ggH4g, mass)
        ggH4g=preselection(ggH4g, mass)
        ggH4g=isolation_vars(ggH4g, mass)
        ggH4g=masses(ggH4g, mass)
        ggH4g=IDs(ggH4g, mass)

        ggH4g=ggH4g.Define('non_MC_cut_m{}'.format(mass),'sample_isMC==0 && best_4g_uncorr_mass_m{m}<90|best_4g_uncorr_mass_m{m}>150'.format(m=mass))
        ggH4g=ggH4g.Define('best_4g_phi1_valid_m{}'.format(mass),'raw_best_4g_m{}[7]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_valid_m{}'.format(mass),'raw_best_4g_m{}[15]'.format(mass))


        
    ggH4g=ggH4g.Filter('sample_isMC==1 | non_MC_cut_m30==1','blinding_data_samples')


    # Create snapshots for ggH3g and ggH4g
    actions.append(ggH4g.Snapshot('ggH4g', f"{sample}_ggH4g.root", cols, opts))

    # Generate and save reports for ggH3g and ggH4g
    save_report(ggH4g, "Report_ggH4g", f"{sample}_ggH4g", opts, actions)

    # Create snapshots for the 'Runs' tree in both files
    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, f"{sample}_ggH4g.root", "", opts))

    return actions

def analysis(data,sample):
    actions=[]
    phi_mass=[7,15,20,30,40,50,55]
    actions.extend(ggH(data,phi_mass,sample))
    return actions
