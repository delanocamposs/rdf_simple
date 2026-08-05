import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal")
from plotting import plot_summary
import argparse
import sys


def run(mass, lifetime, year, cat):
    plot_summary.run(mass, lifetime, year, cat)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Summary plots for year, mass, lifetime, category")
    parser.add_argument("-m", "--mass", type=str, help="mass of sample")
    parser.add_argument("-ct", "--ctau", type=str, help="lifetime of sample")
    parser.add_argument("-y", "--year", type=str, help="year of MC and data (single year or Run2/Run3/2022/2023 aggregates)")
    parser.add_argument("-c1", "--cat1", dest="c1", type=str, help="choose one of: prompt, asym, displaced, none")
    parser.add_argument("-c2", "--cat2", dest="c2", type=str, help="choose one of: prompt, asym, displaced, none")
    parser.add_argument("-c3", "--cat3", dest="c3", type=str, help="choose one of: prompt, asym, displaced, none")
    parser.add_argument("-c4", "--cat4", dest="c4", type=str, help="choose one of: prompt, asym, displaced, none")

    parser.add_argument("-process_run2", "--process_run2", dest="process_run2", type=int, help="Summary plots for all Run 2. 1=yes, 0=no. Runs all mass/lifetime/category points")
    parser.add_argument("-process_run3", "--process_run3", dest="process_run3", type=int, help="Summary plots for all Run 3. 1=yes, 0=no. Runs all mass/lifetime/category points")

    args = parser.parse_args()
    mass = args.mass
    lifetime = args.ctau
    year = args.year
    process_run2 = args.process_run2
    process_run3 = args.process_run3

    all_categories = ["prompt", "asym", "displaced"]

    if process_run2:
        for m in ["15", "20", "30", "40", "50", "55"]:
            for ct in ["0", "10", "20", "50", "100", "1000"]:
                for year in ["2017", "2018"]:  # no longer using 2016 because the trigger is too inefficient
                    for cat in all_categories:
                        run(m, ct, year, cat)

    elif process_run3:
        for m in ["15", "20", "30", "40", "50", "55"]:
            for ct in ["0", "10", "20", "50", "100", "1000"]:
                for year in ["2022", "2023", "2024"]:
                    for cat in all_categories:
                        run(m, ct, year, cat)

    else:
        categories = [getattr(args, c) for c in ["c1", "c2", "c3", "c4"] if getattr(args, c) is not None]
        if not (mass and lifetime and year and categories):
            print("\033[91mERROR: a single summary plot requires --mass, --ctau, --year and at least one category "
                  "(--cat1/--cat2/--cat3/--cat4). Provide all of these, or pass --process_run2 1 / --process_run3 1 "
                  "to run the full set.\033[0m")
            sys.exit(1)
        for cat in categories:
            run(mass, lifetime, year, cat)
