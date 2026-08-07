from ggHparameters import lxy1, lxy2, dxy_min, signal_window, signal_xsec, BR, lower_sb, upper_sb

triggers = {
    "current":"HLT_passed==1",
    "doublePhoton33":"HLT_DoublePhoton33_CaloIdL==1",
    "triplePhoton":"HLT_TriplePhoton_20_20_20_CaloIdLV2==1",
    "triplePhoton_R9IdVL":"HLT_TriplePhoton_20_20_20_CaloIdLV2_R9IdVL==1",
    "diphoton30_18":"HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId==1",
    "diphoton30_18_mass55":"HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_Mass55==1",
    "diphoton30PV_18PV":"HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55==1",
    "none":""
}

def trigger(name=None):
    if name is None:
        return "HLT_passed==1"
    if name not in triggers:
        raise KeyError(f"unknown trigger '{name}'.known: {sorted(triggers)}")
    return triggers[name]

def trigger_branch(name):
    cut = triggers[name]
    return cut.split("==")[0] if cut else None

def trigger_or(*names):
    parts = [triggers[n] for n in names if triggers.get(n)]
    return "(" + " || ".join(parts) + ")" if parts else ""

def dxy_valid(mass):
    return f"best_4g_phi1_dxy_m{mass}>{dxy_min} && best_4g_phi2_dxy_m{mass}>{dxy_min}"

def preselection(mass):
    return " && ".join(f"Photon_preselection[best_4g_idx{i}_m{mass}]==1" for i in range(1, 5))

def full_id(mass):
    return f"best_4g_ID_m{mass}==1 && best_4g_passBitMap_loose_iso_m{mass}==1"

def pileup():
    return "abs(Pileup_weight)<=10"

def mc_weight(sumw):
    return f"(genWeight / {sumw}) * {signal_xsec} * {BR} * Pileup_weight"

def blind(mass, window=signal_window):
    lo,hi = window
    return f"(best_4g_corr_mass_m{mass}<{lo})||(best_4g_corr_mass_m{mass}>{hi})"

def sidebands(mass, lower=lower_sb, upper=upper_sb):
    lo1,lo2 = lower
    up1,up2 = upper
    return (f"(best_4g_corr_mass_m{mass}>={lo1} && best_4g_corr_mass_m{mass}<{lo2}) || "
            f"(best_4g_corr_mass_m{mass}>{up1} && best_4g_corr_mass_m{mass}<={up2})")

def categories(mass, lxy1=lxy1, lxy2=lxy2):
    p1g=f"best_4g_phi1_dxy_m{mass}>{lxy1}"
    p1l=f"best_4g_phi1_dxy_m{mass}<{lxy1}"
    p2g=f"best_4g_phi2_dxy_m{mass}>{lxy2}"
    p2l=f"best_4g_phi2_dxy_m{mass}<{lxy2}"
    return {
        "displaced":f"({p1g})&&({p2g})",
        "asym": f"({p1g})&&({p2l})||({p1l})&&({p2g})",
        "prompt": f"({p1l})&&({p2l})",
        "none":f"({p1l})&&({p2l})||({p1l})&&({p2g})||({p1g})&&({p2l})||({p1g})&&({p2g})",
            }

def combine(*cuts):
    return " && ".join(f"({c})" for c in cuts if c)
