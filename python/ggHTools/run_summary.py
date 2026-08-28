import warnings
from plotting import plot_summary
from plotting.plot_summary import year_config
from datacard.ggHdatacardmaker import recommended_photon_id
from ggHparameters import signal_masses, signal_ctaus, run2_data_years, run3_data_years
import argparse
import os


run2_years=run2_data_years+["Run2"]
run3_years=run3_data_years+["Run3"]
output_base="summary_plots"

def run(mass,ctau,year,output_dir=None,formats=("png","pdf")):
    default_output_dir=os.path.join(output_base,year,recommended_photon_id(year))
    return plot_summary.run(mass,ctau,year,output_dir or default_output_dir,formats=formats)


def run_all_points(years,masses,ctaus,formats,output_dir=None):
    for year in years:
        for mass in masses:
            for ctau in ctaus:
                run(mass,ctau,year,output_dir=output_dir,formats=formats)

if __name__=="__main__":
    parser=argparse.ArgumentParser("Summary plots for year, mass and lifetime",add_help=False)
    parser.add_argument("-h",action="help")
    parser.add_argument("-m",dest="mass",type=str)
    parser.add_argument("-ct",dest="ctau",type=str)
    parser.add_argument("-y",dest="year",type=str,choices=sorted(year_config))
    parser.add_argument("-masses",dest="masses",type=str,nargs="+")
    parser.add_argument("-ctaus",dest="ctaus",type=str,nargs="+")
    parser.add_argument("-o",dest="output_dir",default=None)
    parser.add_argument("-f",dest="formats",nargs="+",choices=["png","pdf"],default=["png","pdf"])
    parser.add_argument("-run2",dest="process_run2",action="store_true")
    parser.add_argument("-run3",dest="process_run3",action="store_true")
    args=parser.parse_args()
    formats=tuple(args.formats)

    if args.process_run2 or args.process_run3:
        years=[]
        if args.process_run2:
            years+=run2_years
        if args.process_run3:
            years+=run3_years
        masses=args.masses or ([args.mass] if args.mass else [str(m) for m in signal_masses])
        ctaus=args.ctaus or ([args.ctau] if args.ctau else [str(c) for c in signal_ctaus])
    else:
        years=[args.year] if args.year else []
        masses=args.masses or ([args.mass] if args.mass else [])
        ctaus=args.ctaus or ([args.ctau] if args.ctau else [])
        if not (years and masses and ctaus):
            parser.error("provide a year, at least one mass and at least one lifetime, or use -run2/-run3")

    run_all_points(years,masses,ctaus,formats,output_dir=args.output_dir)
