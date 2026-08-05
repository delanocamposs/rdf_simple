import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal")
from plotting.plot_postfit import plot
from ggHparameters import bins, signal_path, bkg_path
from run_datacard import run as build_datacard, combine_workflow
import argparse


def run(signal, bkg, cat, year, mass, lifetime, finalstate, physics, bins):
    #postfit plotting depends on the datacard + combine fits:
    #build the card, run combine (MultiDimFit + FitDiagnostics), then plot the postfit
    build_datacard(signal, bkg, cat, year, mass, lifetime, finalstate, physics, bins)
    combine_workflow(cat, year, mass, lifetime, finalstate, physics)
    plot(f"higgsCombineTest.MultiDimFit.mH125_m{mass}_ct{lifetime}_{cat}_{year}.root",
         f"fitDiagnosticsTest_m{mass}_ct{lifetime}_{cat}_{year}.root",
         f"{cat}", f"{year}", bins=bins, finalstate=finalstate, physics=physics,
         mass=mass, lifetime=lifetime)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Postfit plotting for year, category")
    parser.add_argument("-m", "--mass", type=str, help="mass of sample")
    parser.add_argument("-ct", "--ctau", type=str, help="lifetime of sample")
    parser.add_argument("-y", "--year", type=str, help="year of MC and data")
    parser.add_argument("-c1", "--cat1", dest="c1", type=str, help="choose one of: prompt, displaced, asym, none")
    parser.add_argument("-c2", "--cat2", dest="c2", type=str, help="choose one of: prompt, displaced, asym, none")
    parser.add_argument("-c3", "--cat3", dest="c3", type=str, help="choose one of: prompt, displaced, asym, none")
    parser.add_argument("-c4", "--cat4", dest="c4", type=str, help="choose one of: prompt, displaced, asym, none")

    parser.add_argument("-process_run2", "--process_run2", dest="process_run2", type=int, help="Postfit all Run 2. 1=yes, 0=no. Runs all mass/lifetime points for all categories")
    parser.add_argument("-process_run3", "--process_run3", dest="process_run3", type=int, help="Postfit all Run 3. 1=yes, 0=no. Runs all mass/lifetime points for all categories")

    args = parser.parse_args()
    mass = args.mass
    lifetime = args.ctau
    year = args.year
    process_run2 = args.process_run2
    process_run3 = args.process_run3

    if process_run2:
        categories = ["prompt", "asym", "displaced"]
        for m in [30]:
            for ct in [100]:
                for year in ["2018"]:  # no longer using 2016 because the trigger is too inefficient
                    sig = signal_path(m, ct, year)
                    bkg = bkg_path(year)
                    for cat in categories:
                        run(sig, bkg, cat, year, m, ct, "4g", "ggH", bins=bins)

    elif process_run3:
        categories = ["prompt", "asym", "displaced"]
        for m in [15, 20, 30, 40, 50, 55]:
            for ct in [0, 10, 20, 50, 100, 1000]:
                for year in ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]:
                    sig = signal_path(m, ct, year)
                    bkg = bkg_path(year)
                    for cat in categories:
                        run(sig, bkg, cat, year, m, ct, "4g", "ggH", bins=bins)

    else:
        sig = signal_path(mass, lifetime, year)
        bkg = bkg_path(year)
        categories = []
        for c in ["c1", "c2", "c3", "c4"]:
            argu = getattr(args, c)
            if argu is not None:
                categories.append(argu)
        for cat in categories:
            run(sig, bkg, cat, year, mass, lifetime, "4g", "ggH", bins=bins)
