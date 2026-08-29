# ggHTools:

ggHTools are tools for building the parametric datacards and producing all necessary plots related to the gluon fusion analysis

## setup
You need a combine environment set up to run the plotters. Set it up the same way as in VH or some other way. The instructions are at: https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/


## Shared cuts and parameters
To mitigate bugs in any results, all cuts are defined in a single place and should never be defined per plotter file. All plotters inherit the cuts applied to any dataframe (`ggHcuts.py`) and any parameters required for the analysis (`ggHparameters.py`). Both serve as the chief location to limit complications for plotters downstream.

## Plotters
The parent files that drive the plotters are prefaced with "run" and the children are prefaced with "plot" and are under plotting/. The role of the driver files is for an easy, top-level way to interact with the plotters adn to avoid complication of having to understand details under the hood in the children scripts. Each plot type has a driver file:

`run_summary.py` is the driver for `plotting/plot_summary.py`.

`run_upperLimits.py` is the driver for `plotting/plot_upperLimits.py`.

`run_postfit.py` is the driver for `plotting/plot_postfit.py`.

Each driver is easy to run and has configurability depending on what you want. You can run a single mass,lifetime,year point:

# Summary plots

```bash
python3 run_summary.py -m 40 -ct 100 -y 2023
```

Or run many masses or many lifetimes:

```bash
python3 run_summary.py -y 2023 -masses 20 40 -ctaus 0 100
```

Or run it for complete Run 2 or Run 3 across all lifetimes, masses, years (takes some time):

Generate a complete Run 2 or Run 3 scan:

```bash
python3 run_summary.py -run2
python3 run_summary.py -run3
```

Summary plots are written to `summary_plots/<year>/<photon_id>/` by default. You can override this with -o option. 

# Upper limits:

Run a single-year scan or a combined-era scan:

```bash
python3 run_upperLimits.py -y 2023
python3 run_upperLimits.py -run2
python3 run_upperLimits.py -run3
python3 run_upperLimits.py -run23
```

Cached results are loaded from `limits_UL_vs_mass_<era>.json` if you already ran the UL before. If you want to run UL again after some change, you can force a new scan with `-rescan`:

```bash
python3 run_upperLimits.py -y 2023 -rescan
```

# Postfit plots
  You get the point

# Datacards

One can easily create some datacards using the driver script `run_datacard.py` similar to the driver scripts for the plotters. `run_datacard.py` is the easy CLI way to interface with the datacard pipeline. Technically, the upper limit plotter will do the same thing, but if for some reason you wanted to build datacards without making a plot of the UL:

Build one datacard for a mass, lifetime, and data year:

```bash
python3 run_datacard.py -m 40 -ct 100 -y 2018
```

full Run 2:

```bash
python3 run_datacard.py -process_run2 1
```

full Run 3:

```bash
python3 run_datacard.py -process_run3 1
```

The single-point command creates the directory `m40_ct100_2018_4g_ggH` and writes the datacard-related ROOT and JSON files there. The driver automatically applies the recommended EGM photon ID for the selected period; it has no category or custom-ID option. It builds datacards only: converting them to workspaces or running Combine requires separate `text2workspace.py` and `combine` commands.

## More details
# Datacards
`datacard/` has files that go into the details on how the datacards are built. It deals alot with the fitting infrastructure and building RooWorkSpaces and putting the parametric model in the workspace and what not. Its the under-the-hood details on building the datacards.
# Plotting
`plotting/` of course contains the main plotters that do the dirty work to generate the canvases, but there are also some auxilliary scripts I added under `plotting/studies/.`These are plotters not used to necessarily show results but rather useful in studies we performed that are not necessarily the "killer plots" for the analysis
