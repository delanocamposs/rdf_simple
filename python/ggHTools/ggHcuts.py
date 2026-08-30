from ggHparameters import (dxy_min, dxy_max, signal_xsec, BR,lower_sb, upper_sb, delta_r_cut)

#basic HLT passed cut
def trigger():
    return "HLT_passed==1"

#trigger_and_pT ensures photons in events pass speciofic pT requirements and HLT
#for Run3: pT1>32, pT2>20, pT3,pT4>10
#for Run2: pT1,pT2,pT3>22
def trigger_and_pT(mass):
    return f"best_4g_passTriggerProxy_m{mass}==1"

#valid lxy cut (-20<lxy<110)
def dxy_valid(mass, minimum=dxy_min, maximum=dxy_max):
    return (f"best_4g_phi1_dxy_m{mass}>{minimum} && best_4g_phi1_dxy_m{mass}<{maximum} && best_4g_phi2_dxy_m{mass}>{minimum} && best_4g_phi2_dxy_m{mass}<{maximum}")

#standard preselection cut for best 4 photons
def preselection(mass):
    photon_preselection = " && ".join(f"Photon_preselection[best_4g_idx{i}_m{mass}]==1" for i in range(1, 5))
    return photon_preselection

#4 best photons per event pass dR isolation
#not just chekcing the final pair, all combos
# we dont want any of the best 4 photons to be close
def deltaR(mass, minimum=delta_r_cut):
    indices = [f"best_4g_idx{i}_m{mass}" for i in range(1, 5)]
    pair_expressions = []
    for first, second in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
        pair_expressions.append(f"ROOT::VecOps::DeltaR(Photon_eta[{indices[first]}], Photon_eta[{indices[second]}], Photon_phi[{indices[first]}], Photon_phi[{indices[second]}])")
    minimum_expression = pair_expressions[0]
    for pair_expression in pair_expressions[1:]:
        minimum_expression = f"std::min({minimum_expression}, {pair_expression})"
    return f"({minimum_expression})>{minimum}"

#photon id definiton dynamic to working point based on recommended EGM ID
def photon_id(mass, workingpoint):
    working_points = {"LooseEGM": "Loose","MediumEGM": "Medium","TightEGM": "Tight"}
    wp = working_points[workingpoint]
    return f"best_4g_ID_EGM_{wp}_m{mass}==1"

def pileup():
    return "abs(Pileup_weight)<=10" ##this is genuinely a bandaid. found many MC events with large pileup weights 
                                    ## and we werent sure why. plans to iron that out, but for now this will do

#standard MC weight per event (not really a cut)
def mc_weight(sumw):
    return f"(genWeight / {sumw}) * {signal_xsec} * {BR} * Pileup_weight"

#restrict events to the sidebands defined by ggHparameters.py
def sidebands(mass, lower=lower_sb, upper=upper_sb):
    lo1,lo2 = lower
    up1,up2 = upper
    return (f"(best_4g_corr_mass_m{mass}>={lo1} && best_4g_corr_mass_m{mass}<{lo2}) || (best_4g_corr_mass_m{mass}>{up1} && best_4g_corr_mass_m{mass}<={up2})")

#helper function to easily combine cuts
def combine(*cuts):
    return " && ".join(f"({c})" for c in cuts if c)

#data cuts
def background_selection(mass, workingpoint):
    return combine(trigger_and_pT(mass),dxy_valid(mass),preselection(mass),deltaR(mass),photon_id(mass, workingpoint),sidebands(mass))

#signal MC cuts
def signal_selection(mass, workingpoint):
    return combine(trigger_and_pT(mass),dxy_valid(mass),preselection(mass),deltaR(mass),photon_id(mass, workingpoint),pileup())
