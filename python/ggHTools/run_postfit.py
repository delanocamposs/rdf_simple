import argparse
import os
import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal")
from datacard.ggHdatacardmaker import (main as make_datacard,recommended_photon_id)
from ggHparameters import (bins, bkg_path, lumi, order_fit, run2_data_years,run3_data_years, signal_ctaus, signal_masses,signal_path, year_config)
from plotting.plot_postfit import plot
from run_datacard import combine_workflow

finalstate = "4g"
physics = "ggH"
output_base = "postfit_plots"
run2_years = run2_data_years+["Run2"]
run3_years = run3_data_years+["Run3"]

def run(mass, ctau, year, output_dir=None, formats=("png", "pdf")):
    # Postfit plotting depends on the datacard + combine fits: build the card,
    # run combine (MultiDimFit + FitDiagnostics), then plot the result.
    config = year_config[year]
    signal_years = config["signal_years"]
    signals = [signal_path(mass, ctau, signal_year) for signal_year in signal_years]
    backgrounds = [bkg_path(data_year) for data_year in config["data_years"]]
    background_input = backgrounds[0] if len(backgrounds) == 1 else backgrounds
    make_datacard(paths=signals+[background_input],isMC=[1]*len(signals)+[0],trees=["ggH4g"]*(len(signals)+1),var=f"best_4g_corr_mass_m{mass}",period=year,bins=bins,lifetime=ctau,mass=mass,finalstate=finalstate,physics=physics,order_fit=order_fit,signal_lumis=[lumi[signal_year] for signal_year in signal_years])
    combine_workflow(year, mass, ctau, finalstate, physics)

    if output_dir is None:
        output_dir = os.path.join(output_base, year, recommended_photon_id(year))
    return plot(f"higgsCombineTest.MultiDimFit.mH125_m{mass}_ct{ctau}_{year}.root",f"fitDiagnosticsTest_m{mass}_ct{ctau}_{year}.root",year,bins=bins,finalstate=finalstate,physics=physics,mass=mass,lifetime=ctau,output_dir=output_dir,formats=formats)


def run_all_points(years, masses, ctaus, formats, output_dir=None):
    for year in years:
        for mass in masses:
            for ctau in ctaus:
                run(mass, ctau, year, output_dir=output_dir, formats=formats)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Inclusive postfit plotting by year", add_help=False)
    parser.add_argument("-h", action="help")
    parser.add_argument("-m", "--mass", dest="mass", type=str)
    parser.add_argument("-ct", "--ctau", dest="ctau", type=str)
    parser.add_argument("-y", "--year", dest="year", type=str,choices=sorted(year_config))
    parser.add_argument("-masses", "--masses", dest="masses", type=str,nargs="+")
    parser.add_argument("-ctaus", "--ctaus", dest="ctaus", type=str,nargs="+")
    parser.add_argument("-o", "--output-dir", dest="output_dir", default=None)
    parser.add_argument("-f", "--formats", dest="formats", nargs="+",choices=["png", "pdf"], default=["png", "pdf"])
    parser.add_argument("-run2", "--run2", "-process_run2", "--process_run2",dest="process_run2", nargs="?", const=1, default=0,type=int, choices=[0, 1])
    parser.add_argument("-run3", "--run3", "-process_run3", "--process_run3",dest="process_run3", nargs="?", const=1, default=0,type=int, choices=[0, 1])
    args = parser.parse_args()
    formats = tuple(args.formats)

    if args.process_run2 or args.process_run3:
        years = []
        if args.process_run2:
            years += run2_years
        if args.process_run3:
            years += run3_years
        masses = args.masses or ([args.mass] if args.mass else [str(mass) for mass in signal_masses])
        ctaus = args.ctaus or ([args.ctau] if args.ctau else [str(ctau) for ctau in signal_ctaus])
    else:
        years = [args.year] if args.year else []
        masses = args.masses or ([args.mass] if args.mass else [])
        ctaus = args.ctaus or ([args.ctau] if args.ctau else [])
        if not (years and masses and ctaus):
            parser.error("provide a year, at least one mass and at least one lifetime, or use -run2/-run3")

    run_all_points(years, masses, ctaus, formats, output_dir=args.output_dir)
