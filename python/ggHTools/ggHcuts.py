from ggHparameters import lxy1, lxy2, dxy_min, signal_window, signal_xsec, BR, lower_sb, upper_sb

def trigger():
    return "HLT_passed==1"

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
