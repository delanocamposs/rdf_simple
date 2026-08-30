import argparse
import os

from plotting import plot_upperLimits
from ggHparameters import signal_masses, signal_ctaus, order_fit, bins

finalstate = "4g"
physics = "ggH"
br_scale = 1e-4
yrange = (3e-6, 1e-2)
br_label = r"$\mathcal{B}(\Phi\rightarrow\gamma\gamma)=1$"
center_of_mass = {"Run2": 13, "Run3": 13.6, "Run2Run3": 13.6}


def run(era_token, masses=None, lifetimes=None, rescan=False):
    era, years = plot_upperLimits.resolve_era(era_token)
    masses = list(signal_masses if masses is None else masses)
    lifetimes = list(signal_ctaus if lifetimes is None else lifetimes)
    results_json = f"limits_UL_vs_mass_{era}.json"
    if os.path.exists(results_json) and not rescan:
        print(f"loading cached limits from {results_json} (pass -rescan to re-run)")
        results = plot_upperLimits.load_results(results_json)
    else:
        results = plot_upperLimits.scan_mass_lifetime(masses, lifetimes, era, years, bins=bins, finalstate=finalstate, physics=physics, order_fit=order_fit, results_json=results_json)
    return plot_upperLimits.plot_UL_vs_mass(results, era, lifetimes, total_lumi=plot_upperLimits.era_lumi(years), br_scale=br_scale, extra_labels=(br_label,), yrange=yrange, center_of_mass=center_of_mass[era])


if __name__ == "__main__":
    parser = argparse.ArgumentParser("UL vs mass limits", add_help=False)
    parser.add_argument("-h", action="help")
    parser.add_argument("-y", dest="year", type=str)
    parser.add_argument("-run2", dest="process_run2", action="store_true")
    parser.add_argument("-run3", dest="process_run3", action="store_true")
    parser.add_argument("-run23", dest="process_run2run3", action="store_true")
    parser.add_argument("-rescan", dest="rescan", action="store_true") ##to regenerate the UL values. otherwise no
    args = parser.parse_args()

    if args.process_run2:
        era_token = "Run2"
    elif args.process_run3:
        era_token = "Run3"
    elif args.process_run2run3:
        era_token = "Run2Run3"
    elif args.year:
        era_token = args.year
    else:
        parser.error("specify -y <year|Run2|Run3|Run2Run3>, or -run2/-run3/-run23")

    run(era_token, rescan=args.rescan)
