import ROOT
ROOT.gInterpreter.Declare('#include "analysis/ddp_vertex.h"')
ROOT.gInterpreter.Declare('#include "common/scaleFactors.h"')
ROOT.gInterpreter.Declare('#include "common/signalEfficiency.h"')
opts = ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"
opts.fOverwriteIfExists = True

from common.pyhelpers import load_meta_data


cols = ".*best_2g.*|sample_.*|^Photon_.*|^Muon_.*|nMuon|nElectron|nPhoton|^Z.*|^W.*|Weight.*|^Gen.*|^weight.*|^TrigObj_.*|^event.*|^Electron_.*|^Pileup_.*|^MET_.*|^run.*|^luminosityBlock.*"


# Muon trigger[era][par], par = ['name', 'bits', 'pt']
# Name = branch name in tree
# Bits = trigger bits for HLT path
# pt = pt threshold for trigger
muTrig = {'2024': [{'name': 'HLT_IsoMu24', 'bits': 8, 'pt': 24}],
          '2023postBPix': [{'name': 'HLT_IsoMu24', 'bits': 8, 'pt': 24}],
          '2023preBPix': [{'name': 'HLT_IsoMu24', 'bits': 8, 'pt': 24}],
          '2022postEE': [{'name': 'HLT_IsoMu24', 'bits': 8, 'pt': 24}],
          '2022preEE': [{'name': 'HLT_IsoMu24', 'bits': 8, 'pt': 24}],
          '2018': [{'name': 'HLT_IsoMu24', 'bits': 8, 'pt': 24}],
          '2017': [{'name': 'HLT_IsoMu27', 'bits': 8, 'pt': 27}],
          '2016postVFP': [{'name': 'HLT_IsoMu24', 'bits': 8, 'pt': 24},
                          {'name': 'HLT_IsoTkMu24', 'bits': 1+8, 'pt': 24}], # Check if filter bits are correct
          '2016preVFP': [{'name': 'HLT_IsoMu24', 'bits': 8, 'pt': 24},
                          {'name': 'HLT_IsoTkMu24', 'bits': 1+8, 'pt': 24}]} # Check if filter bits are correct

eleTrig = {'2024': [{'name': 'HLT_Ele30_WPTight_Gsf', 'bits': 2, 'pt': 32}],
           '2023postBPix': [{'name': 'HLT_Ele30_WPTight_Gsf', 'bits': 2, 'pt': 32}],
           '2023preBPix': [{'name': 'HLT_Ele30_WPTight_Gsf', 'bits': 2, 'pt': 32}],
           '2022postEE': [{'name': 'HLT_Ele30_WPTight_Gsf', 'bits': 2, 'pt': 32}],
           '2022preEE': [{'name': 'HLT_Ele30_WPTight_Gsf', 'bits': 2, 'pt': 32}],
           '2018': [{'name': 'HLT_Ele32_WPTight_Gsf', 'bits': 2, 'pt': 32}],
           '2017': [{'name': 'HLT_Ele32_WPTight_Gsf', 'bits': 1024, 'pt': 32}],
           '2016postVFP': [{'name': 'HLT_Ele27_WPTight_Gsf', 'bits': 2, 'pt': 27}],
           '2016preVFP': [{'name': 'HLT_Ele27_WPTight_Gsf', 'bits': 2, 'pt': 27}]}
          
# Common Object ID:
def muonAna(dataframe, era = '2018'):
    # Common Muon ID definitions (No isolation)
    muons = dataframe.Define("loose_muon", "(Muon_pt>5&&abs(Muon_eta)<2.4&&abs(Muon_dxy)<0.2&&abs(Muon_dz)<0.5&&Muon_pfIsoId>1&&Muon_looseId>0)")
    muons = muons.Define("tight_muon", "(loose_muon&&Muon_tightId>0&&Muon_pfIsoId>3)")
    muons = muons.Define("overlap_muon", "(Muon_pt>3&&abs(Muon_eta)<2.4&&Muon_softId>0)")        

    
    muons = muons.Define("Muon_ntight", "Sum(tight_muon)")
    muons = muons.Define("Muon_nloose", "Sum(loose_muon)")
    #define a muon for photon overlap
    
    for trigger in muTrig[era]:
        muons = muons.Define("Muon_pass{}".format(trigger['name']), "matchTrigger(Muon_eta, Muon_phi, Muon_pdgId, TrigObj_eta, TrigObj_phi, TrigObj_pt, TrigObj_id, TrigObj_filterBits, {}, {})".format(trigger['bits'], trigger['pt']))
    muons = muons.Define("Muon_isTrigger", "||".join(["Muon_pass{}".format(trig['name']) for trig in muTrig[era]]))

    #For scale factors redefine them so the uncertainty is actually a variation up or down instead of just uncertainty
    
    muons = muons.Define("mu_SFs_reco", 'scaleFactors_2d(abs(Muon_eta), Muon_pt, MU_RECO_{era}_sf, MU_RECO_{era}_binsX, MU_RECO_{era}_binsY, sample_isMC, Muon_pt>10)'.format(era=era))
    muons = muons.Define("Muon_recoSF_val", "mu_SFs_reco[0]")
    muons = muons.Define("Muon_recoSF_up", "mu_SFs_reco[0]+mu_SFs_reco[1]")
    muons = muons.Define("Muon_recoSF_down", "mu_SFs_reco[0]-mu_SFs_reco[1]")
    
    muons = muons.Define("mu_SFs_id", 'scaleFactors_2d(abs(Muon_eta), Muon_pt, MU_ID_{era}_sf, MU_ID_{era}_binsX, MU_ID_{era}_binsY, sample_isMC, tight_muon)'.format(era=era))
    muons = muons.Define("Muon_idSF_val", "mu_SFs_id[0]")
    muons = muons.Define("Muon_idSF_up", "mu_SFs_id[0]+mu_SFs_id[1]")
    muons = muons.Define("Muon_idSF_down", "mu_SFs_id[0]-mu_SFs_id[1]")

    muons = muons.Define("mu_SFs_iso", 'scaleFactors_2d(abs(Muon_eta), Muon_pt, MU_ISO_{era}_sf, MU_ISO_{era}_binsX, MU_ISO_{era}_binsY, sample_isMC, tight_muon)'.format(era=era))
    muons = muons.Define("Muon_isoSF_val", "mu_SFs_iso[0]")
    muons = muons.Define("Muon_isoSF_up", "mu_SFs_iso[0]+mu_SFs_iso[1]")
    muons = muons.Define("Muon_isoSF_down", "mu_SFs_iso[0]-mu_SFs_iso[1]")

    muons = muons.Define("mu_SFs_trig", 'scaleFactors_3d(Muon_charge, Muon_eta, Muon_pt, MU_TRIG_{era}_sf, MU_TRIG_{era}_binsX, MU_TRIG_{era}_binsY, MU_TRIG_{era}_binsZ, sample_isMC, Muon_isTrigger)'.format(era=era))
    muons = muons.Define("Muon_trigSF_val", "mu_SFs_trig[0]")
    muons = muons.Define("Muon_trigSF_up", "mu_SFs_trig[0]+mu_SFs_trig[1]")
    muons = muons.Define("Muon_trigSF_down", "mu_SFs_trig[0]-mu_SFs_trig[1]")
    

    return muons

def electronAna(dataframe, era = '2018'):

    # Common Electron ID definitions
    electrons = dataframe.Define("Electron_scEta","Electron_eta + Electron_deltaEtaSC")
    electrons = electrons.Define("loose_electron", "Electron_pt>5&&abs(Electron_eta)<2.5&&(abs(Electron_scEta)>1.57||abs(Electron_scEta)<1.44)&&abs(Electron_dxy)<0.2&&abs(Electron_dz)<0.2&&Electron_lostHits<2&&Electron_convVeto&&Electron_cutBased>0")
    electrons = electrons.Define("tight_electron", "loose_electron&&Electron_cutBased>3")
    electrons = electrons.Define("Electron_ntight", "Sum(tight_electron)")
    electrons = electrons.Define("Electron_nloose", "Sum(loose_electron)")

    for trigger in eleTrig[era]:
        electrons = electrons.Define("Electron_pass{}".format(trigger['name']), "matchTrigger(Electron_eta, Electron_phi, Electron_pdgId, TrigObj_eta, TrigObj_phi, TrigObj_pt, TrigObj_id, TrigObj_filterBits, {}, {})".format(trigger['bits'], trigger['pt']))
    electrons = electrons.Define("Electron_isTrigger", "||".join(["Electron_pass{}".format(trig['name']) for trig in eleTrig[era]]))

    #For scale factors redefine them so the uncertainty is actually a variation up or down instead of just uncertainty

    
    electrons = electrons.Define("ele_SFs_id", 'scaleFactors_2d(Electron_eta, Electron_pt, ELE_ID_{era}_sf, ELE_ID_{era}_binsX, ELE_ID_{era}_binsY, sample_isMC, tight_electron)'.format(era=era))
    electrons = electrons.Define("Electron_idSF_val", "ele_SFs_id[0]")
    electrons = electrons.Define("Electron_idSF_up", "ele_SFs_id[0]+ele_SFs_id[1]")
    electrons = electrons.Define("Electron_idSF_down", "ele_SFs_id[0]-ele_SFs_id[1]")
    
    electrons = electrons.Define("ele_SFs_reco", 'scaleFactors_eleReco(Electron_eta, Electron_eta, ELE_RECO_ptBelow20_{era}_sf, ELE_RECO_ptBelow20_{era}_binsX, ELE_RECO_ptBelow20_{era}_binsY,  ELE_RECO_ptAbove20_{era}_sf, ELE_RECO_ptAbove20_{era}_binsX, ELE_RECO_ptAbove20_{era}_binsY, sample_isMC, tight_electron)'.format(era=era))
    electrons = electrons.Define("Electron_recoSF_val", "ele_SFs_reco[0]")
    electrons = electrons.Define("Electron_recoSF_up", "ele_SFs_reco[0]+ele_SFs_reco[1]")
    electrons = electrons.Define("Electron_recoSF_down", "ele_SFs_reco[0]-ele_SFs_reco[1]")
    
    electrons = electrons.Define("ele_SFs_trig", "scaleFactors_2d(Electron_eta, Electron_pt, ELE_TRIG_{era}_sf, ELE_TRIG_{era}_binsX, ELE_TRIG_{era}_binsY, sample_isMC, Electron_pt>35)".format(era=era))
    electrons = electrons.Define("Electron_trigSF_val", "ele_SFs_trig[0]")
    electrons = electrons.Define("Electron_trigSF_up", "ele_SFs_trig[0]+ele_SFs_trig[1]")
    electrons = electrons.Define("Electron_trigSF_down", "ele_SFs_trig[0]-ele_SFs_trig[1]")
    
    return electrons

# Must run muon + electron analyzer first to do overlap with leptons
def photonAna(dataframe, era = '2018'):
    # Overlap with loose leptons
    photons = dataframe.Define("Photon_muOverlap", "overlapClean(Photon_phi, Photon_eta, Muon_phi[overlap_muon], Muon_eta[overlap_muon],0.8,0.8)")
    photons = photons.Define("Photon_eleOverlap", "overlapClean(Photon_phi, Photon_eta, Electron_phi[loose_electron], Electron_eta[loose_electron],0.5,0.4)")
    photons = photons.Define("Photon_overlap", "Photon_muOverlap||Photon_eleOverlap")
    
    # Photon Preselection criteria
    photons = photons.Define("Photon_preselection", "Photon_pt>20&&!Photon_pixelSeed&&abs(Photon_eta)<2.5&&(abs(Photon_eta)>1.57||abs(Photon_eta)<1.44)&&(!Photon_overlap)&&(Photon_isScEtaEE||Photon_isScEtaEB)&&passPhIso(Photon_vidNestedWPBitmap)")

    # Common Photon ID definitions (No isolation)
    photons = photons.Define("Photon_passCutBasedID","Photon_preselection>0 && Photon_cutBased>0")
#    photons = photons.Define("Photon_passPhIso" , "passPhIso(Photon_vidNestedWPBitmap)")
#    photons = photons.Define("Photon_IDandIso" , "Photon_passCutBasedID&&Photon_passPhIso")


    #For scale factors redefine them so the uncertainty is actually a variation up or down instead of just uncertainty


    photons = photons.Define("pho_SFs_id", "scaleFactors_2d(Photon_eta, Photon_pt, PHO_ID_{era}_sf, PHO_ID_{era}_binsX, PHO_ID_{era}_binsY, sample_isMC, Photon_passCutBasedID)".format(era=era))
    photons = photons.Define("Photon_idSF_val", "pho_SFs_id[0]")
    photons = photons.Define("Photon_idSF_up", "pho_SFs_id[1]+pho_SFs_id[0]")
    photons = photons.Define("Photon_idSF_down", "pho_SFs_id[0]-pho_SFs_id[1]")

    photons = photons.Define("pho_SFs_pix", "getPixelSeedSF(Photon_isScEtaEB, Photon_isScEtaEE, hasPix_{}_sf, sample_isMC, !Photon_pixelSeed)".format(era))
    photons = photons.Define("Photon_pixSF_val", "pho_SFs_pix[0]")
    photons = photons.Define("Photon_pixSF_up", "pho_SFs_pix[0]+pho_SFs_pix[1]")
    photons = photons.Define("Photon_pixSF_down", "pho_SFs_pix[0]-pho_SFs_pix[1]")

    photons = photons.Define("Photon_energyScaleUp", "photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaledown_{}_val, PHO_scaledown_{}_bins, sample_isMC)".format(era, era))
    photons = photons.Define("Photon_energyScaleDown", "photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaledown_{}_val, PHO_scaleup_{}_bins, sample_isMC)".format(era, era))

    #Find the closest Jet and copy the pileUpId branch
    #photons = photons.Define("Photon_jetPUId", "photon_closest_jet_puID(Photon_jetIdx,Jet_pt,Jet_puId)")    
    #Find the closest Tau and copy the decay mode and the mass:
    photons=photons.Define("Photon_tauMass","photon_closest_one_prong_tau_mass(Photon_eta, Photon_phi,Tau_pt,Tau_eta, Tau_phi,Tau_decayMode,Tau_mass)")
    return photons    

def genAna(dataframe):
    gen = dataframe.Define("GenPart_ctau", "getctau(GenPart_dx, GenPart_dy, GenPart_dz, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
    gen = gen.Define("GenPart_ecalEta", "getGenScEta(GenPart_vx, GenPart_vy, GenPart_vz, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
    gen = gen.Define("GenPart_ecalPhi", "getGenScPhi(GenPart_vx, GenPart_vy, GenPart_vz, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
    gen = gen.Define("GenPart_isSignal", "isGenSignal(GenPart_pdgId, GenPart_genPartIdxMother)")
    gen = gen.Define("GenPart_nSignal", "Sum(GenPart_isSignal)")
    return gen

def makeZ(dataframe, lepton):
    # pT, eta, phi, mass, deltaR, deltaPhi
    Zs = dataframe.Define("Z_idx", "best_z({L}_pt, {L}_eta, {L}_phi, {L}_mass, {L}_charge, tight_{l})".format(L=lepton, l=lepton.lower()))
    Zs = Zs.Define("bestZ_info", "best_Z_info({L}_pt, {L}_eta, {L}_phi, {L}_mass, Z_idx)".format(L=lepton, l=lepton.lower()))
    Zs = Zs.Define("Z_pt", "bestZ_info[0]")
    Zs = Zs.Define("Z_eta", "bestZ_info[1]")
    Zs = Zs.Define("Z_phi", "bestZ_info[2]")
    Zs = Zs.Define("Z_mass", "bestZ_info[3]")
    Zs = Zs.Define("Z_deltaR", "bestZ_info[4]")
    Zs = Zs.Define("Z_deltaPhi", "bestZ_info[5]")

    return Zs

def makeZ_fsr(dataframe, lepton):
    Zgs = dataframe.Define("Zg_idx", "best_zg({L}_pt, {L}_eta, {L}_phi, {L}_mass, {L}_charge, tight_{l}, Photon_pt, Photon_eta, Photon_phi)".format(L=lepton, l=lepton.lower()))
    Zgs = Zgs.Define("bestZg_info", "best_Zg_info({L}_pt, {L}_eta, {L}_phi, {L}_mass, Photon_pt, Photon_eta, Photon_phi, Zg_idx)".format(L=lepton, l=lepton.lower()))
    Zgs = Zgs.Define("Zg_pt", "bestZg_info[0]")
    Zgs = Zgs.Define("Zg_eta", "bestZg_info[1]")
    Zgs = Zgs.Define("Zg_phi", "bestZg_info[2]")
    Zgs = Zgs.Define("Zg_mass", "bestZg_info[3]")
    Zgs = Zgs.Define("Zg_deltaR_l1g", "bestZg_info[4]")
    Zgs = Zgs.Define("Zg_deltaR_l2g", "bestZg_info[5]")
    Zgs = Zgs.Define("Zg_deltaPhi_l1g", "bestZg_info[6]")
    Zgs = Zgs.Define("Zg_deltaPhi_l2g", "bestZg_info[7]")
    Zgs = Zgs.Define("Zg_l1_idx", "Zg_idx[0]")
    Zgs = Zgs.Define("Zg_l2_idx", "Zg_idx[1]")
    Zgs = Zgs.Define("Zg_g_idx", "Zg_idx[2]")
    
    return Zgs

def makeW(dataframe, lepton):
    Ws = dataframe.Define("bestW_info", "best_W_info({L}_pt, {L}_eta, {L}_phi, {L}_mass, tight_{l}, MET_pt, MET_phi)".format(L=lepton, l = lepton.lower()))
    Ws = Ws.Define("W_pt", "bestW_info[0]")
    Ws = Ws.Define("W_eta", "bestW_info[1]")
    Ws = Ws.Define("W_phi", "bestW_info[2]")
    Ws = Ws.Define("W_mass", "bestW_info[3]")
    Ws = Ws.Define("W_deltaR", "bestW_info[4]")
    Ws = Ws.Define("W_deltaPhi", "bestW_info[5]")
    Ws = Ws.Define("W_mt", "bestW_info[6]")
    Ws = Ws.Define("W_l1_idx", "bestW_info[7]")
    return Ws
    
def zeeH(data,phi_mass,sample):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    ####################
    #ANALYSIS CODE HERE#        
    ####################
    #pass HLT
    zee = dataframe['Events'].Filter("isGoodLumi", "passed_lumiFilter")
    zee = zee.Filter('HLT_passed','passed_HLT')
    if data['isMC']:
        zee = zee.Define("Pileup_weight", "getPUweight(Pileup_nPU, puWeight_{},sample_isMC)".format(data['era']))
        if data['customNanoAOD']:
            zee = genAna(zee)
    #Apply Lepton ID (no ISO)
    zee = muonAna(zee, data['era'])
    zee = electronAna(zee, data['era'])
    #Look for + and - muons
    zee = zee.Filter("Sum(Electron_charge[tight_electron]==1)>0 && Sum(Electron_charge[tight_electron]==-1)>0","Tight_eplus_eminus")

    #create the best Zmumu candidate and filter
    zee = makeZ(zee, "Electron")  
    zee = zee.Filter("Z_mass>70&&Z_mass<110", "dielectron_mass_70to110")  
    #Apply Thresholds to the muon pts and cut on muon pf iso
    ptThresh = 35
    zee = zee.Filter("(Electron_pt[Z_idx[0]]>{p}||Electron_pt[Z_idx[1]]>{p})".format(p=ptThresh), "leading_electron_ptOver{}".format(int(ptThresh)))

    #require just a super cluster
    zee = zee.Filter("nPhoton>0","at_least_1_photon")
    
    #Apply photon ID (no ISO)
    zee = photonAna(zee, data['era'])

    # FSR Recovery with loose muons and preselection photons
    zee=zee.Define("Photon_isFSR","fsr_recovery(Z_idx,Electron_pt, Electron_eta, Electron_phi, Electron_mass,Photon_pt,Photon_eta,Photon_phi,Photon_preselection)")
    # Correct isolation for FSR
    #zee = zee.Define("Photon_pfRelIso03_fsrCorr", "correct_gammaIso_for_muons(Z_idx, Electron_pt, Electron_eta, Electron_phi, Photon_pt, Photon_eta, Photon_phi, Photon_pfRelIso03_all, Photon_isFSR)")
        
    ##### gamma gamma +X analysis ######
    ################################git ####

    #Al least Two good photons
    zee2g = zee.Filter("Sum(Photon_preselection==1)>1", "at_least_2_preselection_photons")
    #zee2g = zee2g.Filter("Sum(Photon_preselection==1 && Photon_isFSR==0)>1", "no_fsr_photons")
    for mass in phi_mass:
        zee2g=zee2g.Define('raw_best_2g_m{}'.format(mass),'best_2gamma(Photon_pt,Photon_eta,Photon_phi,Photon_isScEtaEB, Photon_isScEtaEE, Photon_preselection, Photon_passCutBasedID, {})'.format(float(mass)))
        zee2g=zee2g.Define('best_2g_dxy_m{}'.format(mass),'raw_best_2g_m{}[6]'.format(mass))
        zee2g=zee2g.Define('best_2g_valid_m{}'.format(mass),'raw_best_2g_m{}[7]'.format(mass))
        zee2g=zee2g.Define('best_2g_raw_mass_m{}'.format(mass),'raw_best_2g_m{}[8]'.format(mass))
        zee2g=zee2g.Define('best_2g_idx1_m{}'.format(mass),'raw_best_2g_m{}[9]'.format(mass))
        zee2g=zee2g.Define('best_2g_idx2_m{}'.format(mass),'raw_best_2g_m{}[10]'.format(mass))
        zee2g=zee2g.Define('best_2g_pt_m{}'.format(mass), 'raw_best_2g_m{}[15]'.format(mass))
        zee2g=zee2g.Define('fsr_best_2g_m{}_info'.format(mass), 'Zgg_fsr(Photon_pt, Photon_eta, Photon_phi, best_2g_idx1_m{}, best_2g_idx2_m{}, Electron_pt, Electron_eta, Electron_phi, Electron_mass, Z_idx)'.format(mass,mass))
        zee2g=zee2g.Define('best_2g_fsr1_m{}'.format(mass), "fsr_best_2g_m{}_info[0]".format(mass))
        zee2g=zee2g.Define('best_2g_fsr2_m{}'.format(mass), "fsr_best_2g_m{}_info[1]".format(mass))
        zee2g=zee2g.Define('best_2g_fsr3_m{}'.format(mass), "fsr_best_2g_m{}_info[2]".format(mass))        

    actions.append(zee2g.Snapshot('zee2g',sample+'.root',cols,opts))
    report = ROOT.RDataFrame(1)
    r = zee2g.Report()
    for cut in r:
        report = report.Define("report_{}_all".format(cut.GetName()), "{}".format(cut.GetAll()))
        report = report.Define("report_{}_pass".format(cut.GetName()), "{}".format(cut.GetPass()))
    actions.append(report.Snapshot("Report_zee2g", sample+'.root', "", opts))
    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, sample+'.root', "", opts))

    return actions

def wenuH(data,phi_mass,sample):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    ####################
    #ANALYSIS CODE HERE#        
    ####################
    wen = dataframe['Events'].Filter("isGoodLumi", "passed_lumiFilter")
    wen = wen.Filter('HLT_passed', 'passed_HLT')
    if data['isMC']:
        wen = wen.Define("Pileup_weight", "getPUweight(Pileup_nPU, puWeight_{}, sample_isMC)".format(data['era']))
        if data['customNanoAOD']:
            wen = genAna(wen)

    wen = muonAna(wen, data['era'])
    wen = electronAna(wen, data['era'])

    wen = wen.Filter("Electron_ntight==1", "exactly_1_tight_electron")
    ptThresh = 35
    wen = wen.Filter("Sum(Electron_pt[tight_electron]>{})>0".format(ptThresh), "electron_pt_over{}".format(ptThresh))
    wen = makeW(wen, "Electron")

    wen = wen.Filter("nPhoton>0", "at_least_1_photon")

    wen = photonAna(wen, data['era'])
    
    wen2g = wen.Filter('Sum(Photon_preselection==1)>1', "at_least_2_preselection_photons")
    
    for mass in phi_mass:
        wen2g=wen2g.Define('raw_best_2g_m{}'.format(mass),'best_2gamma(Photon_pt,Photon_eta,Photon_phi,Photon_isScEtaEB, Photon_isScEtaEE, Photon_preselection, Photon_passCutBasedID, {})'.format(float(mass)))
        wen2g=wen2g.Define('best_2g_dxy_m{}'.format(mass),'raw_best_2g_m{}[6]'.format(mass))
        wen2g=wen2g.Define('best_2g_valid_m{}'.format(mass),'raw_best_2g_m{}[7]'.format(mass))
        wen2g=wen2g.Define('best_2g_raw_mass_m{}'.format(mass),'raw_best_2g_m{}[8]'.format(mass))
        wen2g=wen2g.Define('best_2g_idx1_m{}'.format(mass),'raw_best_2g_m{}[9]'.format(mass))
        wen2g=wen2g.Define('best_2g_idx2_m{}'.format(mass),'raw_best_2g_m{}[10]'.format(mass))
        wen2g=wen2g.Define('best_2g_pt_m{}'.format(mass), 'raw_best_2g_m{}[15]'.format(mass))
        wen2g = wen2g.Define("misID_info_m{}".format(mass), "check_W_misID(Electron_pt[W_l1_idx], Electron_eta[W_l1_idx], Electron_phi[W_l1_idx], Electron_mass[W_l1_idx], Photon_pt, Photon_eta, Photon_phi, best_2g_idx1_m{m}, best_2g_idx2_m{m})".format(m=mass))
        wen2g = wen2g.Define("best_2g_misID1_m{}".format(mass), "misID_info_m{}[0]".format(mass))
        wen2g = wen2g.Define("best_2g_misID2_m{}".format(mass), "misID_info_m{}[1]".format(mass))
        wen2g = wen2g.Define("best_2g_misID3_m{}".format(mass), "misID_info_m{}[2]".format(mass))

    actions.append(wen2g.Snapshot('wen2g', sample+".root", cols, opts))
    report = ROOT.RDataFrame(1)
    r = wen2g.Report()
    for cut in r:
        report = report.Define("report_{}_all".format(cut.GetName()), "{}".format(cut.GetAll()))
        report = report.Define("report_{}_pass".format(cut.GetName()), "{}".format(cut.GetPass()))
    actions.append(report.Snapshot("Report_wen2g", sample+'.root', "", opts))
    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, sample+".root", "", opts))

    return actions

def wmunuH(data,phi_mass,sample):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    ####################
    #ANALYSIS CODE HERE#        
    ####################
    wmn = dataframe['Events'].Filter("isGoodLumi", "passed_lumiFilter")
    wmn = wmn.Filter('HLT_passed', 'passed_HLT')
    
    if data['isMC']:
        wmn = wmn.Define("Pileup_weight", "getPUweight(Pileup_nPU, puWeight_{}, sample_isMC)".format(data['era']))
        if data['customNanoAOD']:
            wmn = genAna(wmn)
    
    wmn = muonAna(wmn, data['era'])
    wmn = electronAna(wmn, data['era'])
    
    wmn = wmn.Filter("Muon_ntight==1", "exactly_1_tight_muon")
    ptThresh = 28 if data['era']=='2017' else 25
    wmn = wmn.Filter("Sum(Muon_pt[tight_muon]>{})>0".format(ptThresh), "muon_pt_over{}".format(ptThresh))
    wmn = makeW(wmn, "Muon")

    wmn = wmn.Filter("nPhoton>0", "at_least_1_photon")

    wmn = photonAna(wmn, data['era'])
    
    wmn2g = wmn.Filter('Sum(Photon_preselection==1)>1', "at_least_2_preselection_photons")
    
    for mass in phi_mass:
        wmn2g=wmn2g.Define('raw_best_2g_m{}'.format(mass),'best_2gamma(Photon_pt,Photon_eta,Photon_phi,Photon_isScEtaEB, Photon_isScEtaEE, Photon_preselection, Photon_passCutBasedID, {})'.format(float(mass)))
        wmn2g=wmn2g.Define('best_2g_dxy_m{}'.format(mass),'raw_best_2g_m{}[6]'.format(mass))
        wmn2g=wmn2g.Define('best_2g_valid_m{}'.format(mass),'raw_best_2g_m{}[7]'.format(mass))
        wmn2g=wmn2g.Define('best_2g_raw_mass_m{}'.format(mass),'raw_best_2g_m{}[8]'.format(mass))
        wmn2g=wmn2g.Define('best_2g_idx1_m{}'.format(mass),'raw_best_2g_m{}[9]'.format(mass))
        wmn2g=wmn2g.Define('best_2g_idx2_m{}'.format(mass),'raw_best_2g_m{}[10]'.format(mass))
        wmn2g=wmn2g.Define('best_2g_pt_m{}'.format(mass), 'raw_best_2g_m{}[15]'.format(mass))
    actions.append(wmn2g.Snapshot('wmn2g', sample+".root", cols,opts))
    report = ROOT.RDataFrame(1)
    r = wmn2g.Report()
    for cut in r:
        report = report.Define("report_{}_all".format(cut.GetName()), "{}".format(cut.GetAll()))
        report = report.Define("report_{}_pass".format(cut.GetName()), "{}".format(cut.GetPass()))
    actions.append(report.Snapshot("Report_wmn2g", sample+'.root', "", opts))

    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, sample+".root", "", opts))
        
    return actions



def wmugamma(data,sample):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    ####################
    #ANALYSIS CODE HERE#        
    ####################
    wmg = dataframe['Events'].Filter("isGoodLumi", "passed_lumiFilter")
    wmg = wmg.Filter('HLT_passed', 'passed_HLT')
    
    if data['isMC']:
        wmg = wmg.Define("Pileup_weight", "getPUweight(Pileup_nPU, puWeight_{}, sample_isMC)".format(data['era']))
        if data['customNanoAOD']:
            wmg = genAna(wmg)
    
    wmg = muonAna(wmg, data['era'])
    wmg = electronAna(wmg, data['era'])
    
    wmg = wmg.Filter("Muon_ntight==1", "exactly_1_tight_muon")
    ptThresh = 28 if data['era']=='2017' else 25
    wmg = wmg.Filter("Sum(Muon_pt[tight_muon]>{})>0".format(ptThresh), "muon_pt_over{}".format(ptThresh))
    wmg = makeW(wmg, "Muon")
    wmg = photonAna(wmg, data['era'])
    wmg = wmg.Filter("Sum(Photon_preselection)>0", "at_least_1_preselected_photon")

    actions.append(wmg.Snapshot('wmugamma', sample+".root", cols,opts))
    report = ROOT.RDataFrame(1)
    r = wmg.Report()
    for cut in r:
        report = report.Define("report_{}_all".format(cut.GetName()), "{}".format(cut.GetAll()))
        report = report.Define("report_{}_pass".format(cut.GetName()), "{}".format(cut.GetPass()))
    actions.append(report.Snapshot("Report_wmugamma", sample+'.root', "", opts))

    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, sample+".root", "", opts))        
    return actions


def wegamma(data,sample):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    ####################
    #ANALYSIS CODE HERE#        
    ####################
    wen = dataframe['Events'].Filter("isGoodLumi", "passed_lumiFilter")
    wen = wen.Filter('HLT_passed', 'passed_HLT')
    if data['isMC']:
        wen = wen.Define("Pileup_weight", "getPUweight(Pileup_nPU, puWeight_{}, sample_isMC)".format(data['era']))
        if data['customNanoAOD']:
            wen = genAna(wen)

    wen = muonAna(wen, data['era'])
    wen = electronAna(wen, data['era'])

    wen = wen.Filter("Electron_ntight==1", "exactly_1_tight_electron")
    ptThresh = 35
    wen = wen.Filter("Sum(Electron_pt[tight_electron]>{})>0".format(ptThresh), "electron_pt_over{}".format(ptThresh))
    wen = makeW(wen, "Electron")
    wen = photonAna(wen, data['era'])

    weg = wen.Filter('Sum(Photon_preselection==1)>0', "at_least_1_preselection_photons")
    
    actions.append(weg.Snapshot('wegamma', sample+".root", cols, opts))
    report = ROOT.RDataFrame(1)
    r = weg.Report()
    for cut in r:
        report = report.Define("report_{}_all".format(cut.GetName()), "{}".format(cut.GetAll()))
        report = report.Define("report_{}_pass".format(cut.GetName()), "{}".format(cut.GetPass()))
    actions.append(report.Snapshot("Report_wegamma", sample+'.root', "", opts))
    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, sample+".root", "", opts))

    return actions


def zmumuH(data,phi_mass,sample):
    actions=[]
    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)

    ####################
    #ANALYSIS CODE HERE#        
    ####################
    
    #pass HLT
    zmm = dataframe['Events'].Filter("isGoodLumi", "passed_lumiFilter")
    zmm = zmm.Filter('HLT_passed','passed_HLT')
    if data['isMC']:
        zmm = zmm.Define("Pileup_weight", "getPUweight(Pileup_nPU, puWeight_{},sample_isMC)".format(data['era']))
        if data['customNanoAOD']:
            zmm = genAna(zmm)
    #Apply Lepton ID (no ISO)
    zmm = muonAna(zmm, data['era'])
    zmm = electronAna(zmm, data['era'])
    #Look for + and - muons
    zmm = zmm.Filter("Sum(Muon_charge[tight_muon]==1)>0 && Sum(Muon_charge[tight_muon]==-1)>0","Tight_muplus_muminus")

    #create the best Zmumu candidate and filter
    zmm = makeZ(zmm, "Muon")  
    zmm = zmm.Filter("Z_mass>70&&Z_mass<110", "dimuon_mass_70to110")  

    #Apply Thresholds to the muon pts and cut on muon pf iso
    ptThresh = 28 if data['era']=='2017' else 25
    zmm = zmm.Filter("(Muon_pt[Z_idx[0]]>{p}||Muon_pt[Z_idx[1]]>{p})".format(p=ptThresh), "leading_mu_pt_over{}".format(ptThresh))

    #require just a super cluster
    zmm = zmm.Filter("nPhoton>0","at_least_1_photon")
    
    #Apply photon ID (no ISO)
    zmm = photonAna(zmm, data['era'])
    #actions.append(zmm.Snapshot("Events", "zmmg.root", cols, opts))

    # FSR Recovery with loose muons and preselection photons
    zmm=zmm.Define("Photon_isFSR","fsr_recovery(Z_idx,Muon_pt, Muon_eta, Muon_phi, Muon_mass,Photon_pt,Photon_eta,Photon_phi,Photon_preselection)");
    # Correct isolation for FSR
    #zmm = zmm.Define("Photon_pfRelIso03_fsrCorr", "correct_gammaIso_for_muons(Z_idx, Muon_pt, Muon_eta, Muon_phi, Photon_pt, Photon_eta, Photon_phi, Photon_pfRelIso03_all, Photon_isFSR)")
    
    ##### gamma gamma +X analysis ######
    ####################################

    #Al least Two good photons
    zmm2g=zmm.Filter('Sum(Photon_preselection==1)>1','at_least_2_preselection_photons')
    for mass in phi_mass:
        zmm2g=zmm2g.Define('raw_best_2g_m{}'.format(mass),'best_2gamma(Photon_pt,Photon_eta,Photon_phi,Photon_isScEtaEB, Photon_isScEtaEE, Photon_preselection, Photon_passCutBasedID, {})'.format(float(mass)))
        zmm2g=zmm2g.Define('best_2g_dxy_m{}'.format(mass),'raw_best_2g_m{}[6]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_valid_m{}'.format(mass),'raw_best_2g_m{}[7]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_raw_mass_m{}'.format(mass),'raw_best_2g_m{}[8]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_idx1_m{}'.format(mass),'raw_best_2g_m{}[9]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_idx2_m{}'.format(mass),'raw_best_2g_m{}[10]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_pt_m{}'.format(mass), 'raw_best_2g_m{}[15]'.format(mass))
        zmm2g=zmm2g.Define('fsr_best_2g_m{}_info'.format(mass), 'Zgg_fsr(Photon_pt, Photon_eta, Photon_phi, best_2g_idx1_m{}, best_2g_idx2_m{}, Muon_pt, Muon_eta, Muon_phi, Muon_mass, Z_idx)'.format(mass,mass))
        zmm2g=zmm2g.Define('best_2g_fsr1_m{}'.format(mass), "fsr_best_2g_m{}_info[0]".format(mass))
        zmm2g=zmm2g.Define('best_2g_fsr2_m{}'.format(mass), "fsr_best_2g_m{}_info[1]".format(mass))
        zmm2g=zmm2g.Define('best_2g_fsr3_m{}'.format(mass), "fsr_best_2g_m{}_info[2]".format(mass))
    actions.append(zmm2g.Snapshot('zmm2g',sample+'.root',cols, opts))
    #actions.append(zmm3g.Snapshot('Events',sample+'_zmm3g.root',cols))
    #actions.append(zmm4g.Snapshot('Events',sample+'_zmm4g.root',cols))
    report = ROOT.RDataFrame(1)
    r = zmm2g.Report()
    for cut in r:
        report = report.Define("report_{}_all".format(cut.GetName()), "{}".format(cut.GetAll()))
        report = report.Define("report_{}_pass".format(cut.GetName()), "{}".format(cut.GetPass()))
    actions.append(report.Snapshot("Report_zmm2g", sample+'.root', "", opts))
    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, sample+'.root', "", opts))
        #actions.append(dataframe[tree].Snapshot(tree, sample+'_zmm3g.root', "", opts))
        #actions.append(dataframe[tree].Snapshot(tree, sample+'_zmm4g.root', "", opts))
    #r=zmm2g.Report()
    #r.Print()

    return actions
    

def analysis(data,sample):
    phi_mass=[15,20,30,40,50,55]
    actions = []
    actions.extend(zmumuH(data,phi_mass,sample))
    actions.extend(zeeH(data,phi_mass,sample))
    actions.extend(wmunuH(data,phi_mass,sample))
    actions.extend(wenuH(data,phi_mass,sample))
    actions.extend(wmugamma(data,sample))
    actions.extend(wegamma(data,sample))
    return actions
