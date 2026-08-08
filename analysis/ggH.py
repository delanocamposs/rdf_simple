import ROOT
import random
import numpy as np
ROOT.gInterpreter.Declare('#include "analysis/ddp_vertex.h"')
ROOT.gInterpreter.Declare('#include "common/scaleFactors.h"')
ROOT.gInterpreter.Declare('#include "common/signalEfficiency.h"')

opts = ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"
opts.fOverwriteIfExists = True

from common.pyhelpers import load_meta_data


#cols = "best_3g.*|best_4g.*|sample_.*|^Photon_.*|^Muon_.*|^Z_.*|Weight.*|^Gen.*|^weight.*|^TrigObj_.*|^event.*|^Electron_.*|^Pileup_.*|^run.*"

iso = {'2016preVFP': 'Photon_pfRelIso03_all',
       '2016postVFP': 'Photon_pfRelIso03_all',
       '2017': 'Photon_pfRelIso03_all',
       '2018': 'Photon_pfRelIso03_all',
       '2022preEE': 'Photon_pfRelIso03_all_quadratic',
       '2022postEE': 'Photon_pfRelIso03_all_quadratic',
       '2023preBPix': 'Photon_pfRelIso03_all_quadratic',
       '2023postBPix': 'Photon_pfRelIso03_all_quadratic',
       '2024': 'Photon_pfRelIso03_all_quadratic'}
iso_scouting={
        "2022":"ScoutingPhoton_ecalIso",
        "2023":"ScoutingPhoton_ecalIso",
        "2024":"ScoutingPhoton_ecalIso"}

def get_ID_val(var, ID_type):
    cuts_ID={
            "custom":{"endcap":{"iso":[0.3], 
                                "hoe":[0.2], 
                                "sieie":[0.045]}, 
                      "barrel":{"iso":[0.25], 
                                "hoe":[0.3], 
                                "sieie":[0.035]}}
            }

    return cuts_ID[ID_type]['barrel'][var][0], cuts_ID[ID_type]['endcap'][var][0]

egm_wp={"Loose":1,"Medium":2,"Tight":3}
run3_eras={"2022preEE","2022postEE","2023preBPix","2023postBPix","2024"}

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

def photonAna(dataframe, era):
    # Overlap with loose leptons
    #photons = dataframe.Define("Photon_muOverlap", "overlapClean(Photon_phi, Photon_eta, Muon_phi[loose_muon], Muon_eta[loose_muon])")
    #photons = photons.Define("Photon_eleOverlap", "overlapClean(Photon_phi, Photon_eta, Electron_phi[loose_electron], Electron_eta[loose_electron])")
    #photons = photons.Define("Photon_overlap", "Photon_muOverlap||Photon_eleOverlap")

    #this preselection should be the same for both custom ID and standard EGM ID
    photons = dataframe.Define("Photon_preselection", "Photon_pt>20&&!Photon_pixelSeed&&abs(Photon_eta)<2.5&&(abs(Photon_eta)>1.57||abs(Photon_eta)<1.44)&&(Photon_isScEtaEE||Photon_isScEtaEB)")
    #photons = photons.Define("Photon_rho", "fixedGridRhoFastjetAll")
    ph_iso_wp = "phIsoWP_Run3" if era in run3_eras else "phIsoWP_Run2"
    photons=photons.Define("Photon_PhIsoWP_EGM",f"{ph_iso_wp}(Photon_vidNestedWPBitmap)")
    photons=photons.Define("Photon_PassPhIso_LooseEGM","Photon_PhIsoWP_EGM>=1")
    sieie_customEB,sieie_customEE=get_ID_val("sieie","custom")
    HoE_customEB,HoE_customEE=get_ID_val("hoe","custom")
    hoe = "Photon_hoe_PUcorr" if era in run3_eras else "Photon_hoe"
    photons=photons.Define("Photon_IdNoIso_custom",f"(({hoe}<{HoE_customEB}&&Photon_isScEtaEB&&Photon_sieie<{sieie_customEB})||({hoe}<{HoE_customEE}&&Photon_isScEtaEE&&Photon_sieie<{sieie_customEE}))")
    # obviousyly built run 3 EGM Photon ID with: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun3
    for wp,level in egm_wp.items():
        photons=photons.Define(f"Photon_passEGM{wp}ID",f"Photon_cutBased>={level}")
        photons=photons.Define(f"Photon_passFullCutBasedID_{wp}EGM",f"Photon_preselection&&Photon_passEGM{wp}ID")
    return photons

def save_report(df, report_name, sample, opts, actions):
        report = ROOT.RDataFrame(1)
        r = df.Report()
        for cut in r:
            report = report.Define(f"report_{cut.GetName()}_all", f"{cut.GetAll()}")
            report = report.Define(f"report_{cut.GetName()}_pass", f"{cut.GetPass()}")
        actions.append(report.Snapshot(report_name, f"{sample}.root", "", opts))

def ggH(data,phi_mass,sample):

    def four_gamma(df, mass):
        df=df.Define(f'raw_best_4g_m{mass}',f"best_4gamma(Photon_pt,Photon_eta,Photon_phi,Photon_isScEtaEB,Photon_isScEtaEE,Photon_preselection,Photon_IdNoIso_custom,Photon_corrIso_m{mass},{float(mass)})")
        return df

    def isolation_vars(df, mass):
        for wp,level in egm_wp.items():
            for i in range(1,5):
                df=df.Define(f"Photon_passPhIso_{wp}EGM_gamma{i}_m{mass}",f"Photon_PhIsoWP_EGM[best_4g_idx{i}_m{mass}]>={level}")
 
        #EVENT-LEVEL BOOLEAN IF BEST 4 PHOTONS PASS LOOSE EGM ISOLATION. THIS IS THE RECOMMENDED METHOD TO ISO ID.
        for wp in egm_wp:
            decisions=" && ".join(f"Photon_passPhIso_{wp}EGM_gamma{i}_m{mass}" for i in range(1,5))
            df=df.Define(f"best_4g_passPhIso_{wp}EGM_m{mass}",decisions)
        return df

    def corrected_kinematic_vars(df, mass):
        #pt's
        df=df.Define(f'best_4g_phi1_gamma1_pt_m{mass}',f'raw_best_4g_m{mass}[0]')
        df=df.Define(f'best_4g_phi1_gamma2_pt_m{mass}',f'raw_best_4g_m{mass}[3]')
        df=df.Define(f'best_4g_phi2_gamma1_pt_m{mass}',f'raw_best_4g_m{mass}[8]')
        df=df.Define(f'best_4g_phi2_gamma2_pt_m{mass}',f'raw_best_4g_m{mass}[11]')

        #eta's
        df=df.Define(f'best_4g_phi1_gamma1_eta_m{mass}',f'raw_best_4g_m{mass}[1]')
        df=df.Define(f'best_4g_phi1_gamma2_eta_m{mass}',f'raw_best_4g_m{mass}[4]')
        df=df.Define(f'best_4g_phi2_gamma1_eta_m{mass}',f'raw_best_4g_m{mass}[9]')
        df=df.Define(f'best_4g_phi2_gamma2_eta_m{mass}',f'raw_best_4g_m{mass}[12]')

        #phi's
        df=df.Define(f'best_4g_phi1_gamma1_phi_m{mass}',f'raw_best_4g_m{mass}[2]')
        df=df.Define(f'best_4g_phi1_gamma2_phi_m{mass}',f'raw_best_4g_m{mass}[5]')
        df=df.Define(f'best_4g_phi2_gamma1_phi_m{mass}',f'raw_best_4g_m{mass}[10]')
        df=df.Define(f'best_4g_phi2_gamma2_phi_m{mass}',f'raw_best_4g_m{mass}[13]')

        #id's
        df=df.Define(f'best_4g_phi1_gamma1_id_customCorrIso_m{mass}',f'raw_best_4g_m{mass}[20]')
        df=df.Define(f'best_4g_phi1_gamma2_id_customCorrIso_m{mass}',f'raw_best_4g_m{mass}[21]')
        df=df.Define(f'best_4g_phi2_gamma1_id_customCorrIso_m{mass}',f'raw_best_4g_m{mass}[22]')
        df=df.Define(f'best_4g_phi2_gamma2_id_customCorrIso_m{mass}',f'raw_best_4g_m{mass}[23]')
        
        #lxy's
        df=df.Define(f'best_4g_phi1_dxy_m{mass}',f'raw_best_4g_m{mass}[6]')
        df=df.Define(f'best_4g_phi2_dxy_m{mass}',f'raw_best_4g_m{mass}[14]')
        return df
    
    def uncorrected_kinematic_vars(df, mass):
        df=df.Define(f'Photon_pt_gamma1_m{mass}',f'Photon_pt[best_4g_idx1_m{mass}]')
        df=df.Define(f'Photon_pt_gamma2_m{mass}',f'Photon_pt[best_4g_idx2_m{mass}]')
        df=df.Define(f'Photon_pt_gamma3_m{mass}',f'Photon_pt[best_4g_idx3_m{mass}]')
        df=df.Define(f'Photon_pt_gamma4_m{mass}',f'Photon_pt[best_4g_idx4_m{mass}]')
        df=df.Define(f'Photon_hoe_gamma1_m{mass}',f'Photon_hoe[best_4g_idx1_m{mass}]')
        df=df.Define(f'Photon_hoe_gamma2_m{mass}',f'Photon_hoe[best_4g_idx2_m{mass}]')
        df=df.Define(f'Photon_hoe_gamma3_m{mass}',f'Photon_hoe[best_4g_idx3_m{mass}]')
        df=df.Define(f'Photon_hoe_gamma4_m{mass}',f'Photon_hoe[best_4g_idx4_m{mass}]')
        df=df.Define(f'Photon_sieie_gamma1_m{mass}',f'Photon_sieie[best_4g_idx1_m{mass}]')
        df=df.Define(f'Photon_sieie_gamma2_m{mass}',f'Photon_sieie[best_4g_idx2_m{mass}]')
        df=df.Define(f'Photon_sieie_gamma3_m{mass}',f'Photon_sieie[best_4g_idx3_m{mass}]')
        df=df.Define(f'Photon_sieie_gamma4_m{mass}',f'Photon_sieie[best_4g_idx4_m{mass}]')
        return df

    def indices(df, mass):
        df=df.Define(f'best_4g_idx1_m{mass}',f'raw_best_4g_m{mass}[24]')
        df=df.Define(f'best_4g_idx2_m{mass}',f'raw_best_4g_m{mass}[25]')
        df=df.Define(f'best_4g_idx3_m{mass}',f'raw_best_4g_m{mass}[26]')
        df=df.Define(f'best_4g_idx4_m{mass}',f'raw_best_4g_m{mass}[27]')
        return df

    def masses(df, mass):
        df=df.Define(f'best_4g_phi1_mass_m{mass}',f'raw_best_4g_m{mass}[16]')
        df=df.Define(f'best_4g_phi2_mass_m{mass}',f'raw_best_4g_m{mass}[17]')
        df=df.Define(f'best_4g_uncorr_mass_m{mass}',f'raw_best_4g_m{mass}[18]')
        df=df.Define(f'best_4g_corr_mass_m{mass}',f'raw_best_4g_m{mass}[19]')
        df=df.Define(f'best_4g_pairing_score_m{mass}',f'raw_best_4g_m{mass}[28]')
        return df
    
    def preselection(df, mass):
        df=df.Define(f'best_4g_preselected_m{mass}',f'Photon_preselection[best_4g_idx1_m{mass}] && Photon_preselection[best_4g_idx2_m{mass}] && Photon_preselection[best_4g_idx3_m{mass}] && Photon_preselection[best_4g_idx4_m{mass}]')
        return df
    
    def IDs(df, mass):
        df=df.Define(f'best_4g_sumID_customCorrIso_m{mass}',f'raw_best_4g_m{mass}[20]+raw_best_4g_m{mass}[21]+raw_best_4g_m{mass}[22]+raw_best_4g_m{mass}[23]')
        gamma_labels = {1: 'phi1_gamma1', 2: 'phi1_gamma2', 3: 'phi2_gamma1', 4: 'phi2_gamma2'}
        preselection_str=" && ".join(f"(Photon_preselection[raw_best_4g_m{mass}[{i}]]==1)" for i in range(24,28))
        idnoiso_str=" && ".join(f"(Photon_IdNoIso_custom[raw_best_4g_m{mass}[{i}]]==1)" for i in range(24,28))
        iso_str=" && ".join(f"(Photon_PassPhIso_LooseEGM[raw_best_4g_m{mass}[{i}]]==1)" for i in range(24,28))
        df=df.Define(f'best_4g_ID_custom_m{mass}',f'{preselection_str} && {idnoiso_str} && {iso_str}')
        for wp in egm_wp:
            id_flags=[f'Photon_passFullCutBasedID_{wp}EGM[best_4g_idx{i}_m{mass}]' for i in range(1,5)]
            for i,flag in enumerate(id_flags,1):
                df=df.Define(f'best_4g_{gamma_labels[i]}_id_EGM_{wp}_m{mass}',flag)
            df=df.Define(f'best_4g_sumID_EGM_{wp}_m{mass}'," + ".join(id_flags))
            df=df.Define(f'best_4g_ID_EGM_{wp}_m{mass}'," && ".join(id_flags))
        return df

    def scale_factors(df, era):
        if era in ['2023preBPix','2023postBPix']:
            df=df.Define("pho_SFs_id_LooseEGM",f"scaleFactors_3d(Photon_phi,Photon_eta,Photon_pt,PHO_ID_{era}_sf,PHO_ID_{era}_binsX,PHO_ID_{era}_binsY,PHO_ID_{era}_binsZ,sample_isMC,Photon_passFullCutBasedID_LooseEGM)")
        else:
            df=df.Define("pho_SFs_id_LooseEGM",f"scaleFactors_2d(Photon_eta,Photon_pt,PHO_ID_{era}_sf,PHO_ID_{era}_binsX,PHO_ID_{era}_binsY,sample_isMC,Photon_passFullCutBasedID_LooseEGM)")
        df=df.Define("Photon_idSF_LooseEGM_val","pho_SFs_id_LooseEGM[0]")
        df=df.Define("Photon_idSF_LooseEGM_up","pho_SFs_id_LooseEGM[1]+pho_SFs_id_LooseEGM[0]")
        df=df.Define("Photon_idSF_LooseEGM_down","pho_SFs_id_LooseEGM[0]-pho_SFs_id_LooseEGM[1]")

        df=df.Define("pho_SFs_pix", f"getPixelSeedSF(Photon_isScEtaEB, Photon_isScEtaEE, hasPix_{era}_sf, sample_isMC, !Photon_pixelSeed)")
        df=df.Define("Photon_pixSF_val", "pho_SFs_pix[0]")
        df=df.Define("Photon_pixSF_up", "pho_SFs_pix[0]+pho_SFs_pix[1]")
        df=df.Define("Photon_pixSF_down", "pho_SFs_pix[0]-pho_SFs_pix[1]")
        return df

#analysis will actually begin now. all above was jsut organizing branch definitions
#+++++++++++++++++++++++++++++++++++++++++++++++


    era=data['era']
    cols = "best_.*|sample_.*|^Photon_.*|^Electron_.*|Weight.*|^Gen.*|^weight.*|^TrigObj_.*|^event.*|^Pileup_.*|^run.*|gen.*|.*LHE.*|^PV.*|luminosity|Block|genWeight|HLT_passed|sorted_photon_pt|Pass_L1_DoubleEG15_11|Pass_L1_DoubleEG16_11|Pass_L1_DoubleEG17_11|Pass_L1_DoubleEG_OR"
    actions=[]

    dataframe =load_meta_data(data)
    ggH=dataframe["Events"].Filter("isGoodLumi","passed_lumiFilter")

    if data["isMC"]:
        ggH = ggH.Define("Pileup_weight", f"getPUweight(Pileup_nPU, puWeight_{era}, sample_isMC)")

    #ggH=electronAna(ggH)
    #ggH=muonAna(ggH)
    
    #ggH=ggH.Filter("Sum(loose_muon==1)==0",'muon_veto')
    #ggH=ggH.Filter("Sum(loose_electron==1)==0",'electron_veto')

    ggH=photonAna(ggH,era)
    ggH4g=ggH.Filter('nPhoton>3','at_least_4_photons')
    ggH4g=ggH4g.Filter('Sum(Photon_preselection==1)>3','at_least_3_preselected_photons')

    for mass in phi_mass:
        #Photon_corrIso_m{} needs to be defined first since four_gamma relies on it
        ggH4g=ggH4g.Define(f'Photon_corrIso_m{mass}',f"correct_gammaIso(Photon_pt,Photon_eta,Photon_phi,{iso[data['era']]},Photon_preselection)")

        #organizing the branch definitions into functions that couple related branches
        ggH4g=four_gamma(ggH4g, mass)
        ggH4g=indices(ggH4g, mass)
        ggH4g=corrected_kinematic_vars(ggH4g, mass)
        ggH4g=uncorrected_kinematic_vars(ggH4g, mass)
        ggH4g=preselection(ggH4g, mass)
        ggH4g=isolation_vars(ggH4g, mass)
        ggH4g=masses(ggH4g, mass)
        ggH4g=IDs(ggH4g, mass)

        #ggH4g=ggH4g.Define(f'non_MC_cut_m{mass}',f'sample_isMC==0 && best_4g_uncorr_mass_m{mass}<90|best_4g_uncorr_mass_m{mass}>150')
        ggH4g=ggH4g.Define(f'best_4g_phi1_valid_m{mass}',f'raw_best_4g_m{mass}[7]')
        ggH4g=ggH4g.Define(f'best_4g_phi2_valid_m{mass}',f'raw_best_4g_m{mass}[15]')


    #ggH4g=ggH4g.Filter(f'sample_isMC==1 | non_MC_cut_m{m}==1','blinding_data_samples')

    ggH4g=ggH4g.Define("Photon_passFullCutBasedID_custom","Photon_preselection==1&&Photon_IdNoIso_custom==1&&Photon_PassPhIso_LooseEGM==1")
    ggH4g=scale_factors(ggH4g,era)
    actions.append(ggH4g.Snapshot('ggH4g', f"{sample}_ggH4g.root", cols, opts))

    save_report(ggH4g, "Report_ggH4g", f"{sample}_ggH4g", opts, actions)
    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, f"{sample}_ggH4g.root", "", opts))

    return actions

def Zee(data, sample):
    #analysis only used for tag and probe studies
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
    Zee=dataframe["Events"].Filter("isGoodLumi","passed_lumiFilter")

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
        Zee=Zee.Define("Pileup_weight",f"getPUweight(Pileup_nPU, puWeight_UL{data['era']}, sample_isMC)")

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
