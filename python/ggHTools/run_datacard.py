import argparse
import os
import shutil
import subprocess
import warnings

warnings.filterwarnings("ignore", message="The value of the smallest subnormal")

from datacard.ggHdatacardmaker import main
from ggHparameters import bins, bkg_path, signal_path


def combine_workflow(year, mass, lifetime, finalstate, physics):
    card_txt = f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{year}.txt"
    card_root = f"datacard_{physics}_{finalstate}_m{mass}_ct{lifetime}_{year}.root"

    subprocess.run(["text2workspace.py", card_txt, "-o", card_root], check=True)
    subprocess.run(["combine", card_root, "-M", "MultiDimFit", "--saveWorkspace", "--robustFit", "1", "--cminDefaultMinimizerStrategy", "2", "-m", "125"], check=True)
    subprocess.run(["combine", card_root, "-M", "FitDiagnostics", "--saveShapes", "--saveWorkspace", "--saveWithUncertainties", "--saveNormalizations", "--robustFit", "1", "--cminDefaultMinimizerStrategy", "2", "-m", "125"], check=True)
    subprocess.run(["mv", "higgsCombineTest.MultiDimFit.mH125.root", f"higgsCombineTest.MultiDimFit.mH125_m{mass}_ct{lifetime}_{year}.root"], check=True)
    subprocess.run(["mv", "fitDiagnosticsTest.root", f"fitDiagnosticsTest_m{mass}_ct{lifetime}_{year}.root"], check=True)


def run(signal, bkg, year, mass, lifetime, finalstate, physics, bins):
#colorful otuput to terminal
    mustard = "\033[38;5;136m"
    teal = "\033[38;5;44m"
    reset = "\033[0m"
    border = "============================================================"
    dirname = f"m{mass}_ct{lifetime}_{year}_{finalstate}_{physics}"
    print(f"{teal}{border}{reset}")
    print(f"{mustard}• searching for pre-existing directory: {dirname}{reset}")
    if os.path.isdir(dirname):
        shutil.rmtree(dirname)
        print(f"{mustard}  ↳ found. overriding previous datacard results.{reset}")
    else:
        print(f"{mustard}  ↳ not found, building directory: {dirname}{reset}")

    main(paths=[signal, bkg], isMC=[1, 0], trees=["ggH4g", "ggH4g"],var=f"best_4g_corr_mass_m{mass}", period=year, bins=bins,lifetime=lifetime, mass=mass, lumi_scaling=1)

    print(f"{mustard}• updated directory {dirname}.{reset}")


def process_points(years):
    for mass in [15, 20, 30, 40, 50, 55]:
        for lifetime in [0, 10, 20, 50, 100, 1000]:
            for year in years:
                run(signal_path(mass, lifetime, year), bkg_path(year), year,mass, lifetime, "4g", "ggH", bins=bins)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(" datacard processing by year")
    parser.add_argument("-m", "--mass", type=str, help="mass of sample")
    parser.add_argument("-ct", "--ctau", type=str, help="lifetime of sample")
    parser.add_argument("-y", "--year", type=str, help="year of MC and data")
    parser.add_argument("-process_run2", "--process_run2", dest="process_run2", type=int)
    parser.add_argument("-process_run3", "--process_run3", dest="process_run3", type=int)
    args = parser.parse_args()

    if args.process_run2:
        process_points(["2017", "2018"])
    elif args.process_run3:
        process_points(["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"])
    else:
        if not (args.mass and args.ctau and args.year):
            parser.error("mass, lifetime, and year must be specified")

        run(signal_path(args.mass, args.ctau, args.year), bkg_path(args.year),args.year, args.mass, args.ctau, "4g", "ggH", bins=bins)
