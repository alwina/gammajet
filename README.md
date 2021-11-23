This is the code for the measurement of isolated photon-jet correlations in 5.02 TeV Pb-Pb collisions in ALICE.

The starting point is NTuples generated with [this repo](https://github.com/alwina/ntuple-gj) (which is forked from [this repo](https://github.com/yslai/ntuple-gj)) from ALICE data and simulations. The code to calculate the photon purity is a refined version of what's in [this repo](https://github.berkeley.edu/alwina/photon-correlations). The code to calculate the actual correlations is based on [this repo](https://github.com/ftoralesacosta/GH_Correlations).

## Configuration files
`yaml` is used for the configuration file. It contains nearly all of the information needed to run the analysis (after the ntuples are made). It's slightly recursive in the sense that the photon purity is part of the configuration, but the file is also used to calculate the photon purity. So after calculating the purity with a given configuration, one has to then go back and re-edit it to update the purity.

There are 3 levels of configuration files used in this analysis
- A global-level configuration which defines the observables being measured, and which is meant to be universal across all systems, datasets, etc that this analysis will run over
- A system-level configuration which defines things like the isolation and shower shape parameters or which kinds of jets to look at; this is meant to be the same across, say, multiple run periods or between small systems. It points to the global configuration
- A run-period-level (or equivalent) configuration, which defines things like filenames and the purity. It points to a system-level configuration

Anything in a more specific configuration file will overwrite anything in an earlier one. So, for example, the isolation can be defined in a run-period-level configuration, and it will overwrite whatever is in the system-level configuration. One thing to note is that, because of the way that python handles yaml and dictionaries and such, top-level nodes will be completely overwritten. So, for example, in order to change the isolation variable, one has to define the isolation variables *and also all of the cut values* in the run-period-level configuration. Just defining the `isovar` will cause the others to not be set.

## Purity calculation
To calculate the purity, the ntuples are first converted to CSVs with `make_cluster_csv.py`. Those CSVs are then read into pandas dataframes, from which template fits are done. The `DataframeCollection` and `PhotonPurity` classes handle most of this.

## Aux file
In the interest of not rerunning the entire ntuplizer and also not modifying 100+ GB files in-place, in order to use information that is derived from the ntuples but not saved in them, an auxiliary file gets created. At the moment, this includes `cluster_5x5all` and `jet_ak02tpc` information, but it can include anything that is better off calculated once to be used later. To create an aux file for an existing ntuple.root, within a python3 console (with ROOT), run the following:
```
from make_aux_file import createAuxFile
createAuxFile(ntuple.root)
```
This will generate a new file called `ntuple_AUX.root`. In principle, running `make_aux_file` with a configuration file passed as an argument will make an aux file for all files defined within a configuration, but that can take a very long time because the process of creating an aux file is quite slow (about 2000 events/minute). So it's probably better not to do it that way, but it is an option.

## C++ scripts
The same-event correlations, mixed-event correlations, and RooUnfold responses are done with C++ scripts. `cmake` (minimum version 3.10) is used to compile these scripts, which have the following requirements
- ROOT, for which `$ROOT_DIR` has to be defined and which should then be automatically found with `find_package(ROOT)`
- [yaml-cpp](https://github.com/jbeder/yaml-cpp), for which `$YAMLCPP` has to be defined
- HDF5 (mixed_event only), for which `$HDF5_DIR` has to be defined (for some reason, `find_package` fails here)
- [RooUnfold](https://gitlab.cern.ch/RooUnfold/RooUnfold) (response_matrix only), for which `$ROOUNFOLD_DIR` has to be defined

Generally speaking, here are the steps to compile the scripts (after building whatever necessary dependencies and/or loading whatever relevant modules):
1. `cd cpp`
2. `mkdir build`
3. `cd build`
4. `cmake ..` -- pay attention to any warnings/errors that come up
5. `make` (to compile all 3 scripts) or `make <filename>` (without the `.cc` extension to compile any particular script)

## Same-event correlations
Running the same-event correlation is maybe the easiest part. Once the purity has been defined in the configuration file, do
```
./cpp/build/same_event config-run-period.yaml
```

## Event mixing
The event mixing flow is unfortunately somewhat complicated. The [Event Mixing](https://github.com/ftoralesacosta/Event_Mixing) repo is included as a submodule for this directory. It is used to do the first two steps below. The rest of the steps happen in this repo.

The configuration file has a section within `filelists` for `mixing`. Within this section, each subsection has some `mixlabel`, for which the `triggered`, `mb`, and `pairing` files are defined. This way, as much or as little mixing can be done at once, and file names don't have to be changed all the time.

The steps are as follows, beginning with a `triggeredntuple.root` and a `mbntuple.root`:
1. Assuming the aux files exist, use `to_hdf5_with_aux` in `Event_Mixing` repo to convert to from ROOT files to HDF5
```
./to_hdf5_with_aux triggeredntuple.root triggeredntuple.hdf5
./to_hdf5_with_aux mbntuple.root mbntuple.root
```
The names of these HDF5 files go into the config file as the `triggered` and `mb` filenames within a given `mixlabel`.

2. Do the pairing. The name of this pairing file goes into the config file as the `pairing` filename within that same `mixlabel`. This only has to be done once per triggered-MB file pair, unless the pairing criteria change. So, for example, adding something to the aux file and remaking the HDF5 files in order to look at another variable does not require the pairing to be rerun, but creating new HDF5 files with a subset of the full events does require the pairing to be rerun (because there's a new triggered-MB file pair).

3. Open an interactive shell on cori in order to run things in parallel:
```
salloc -N 1 -C haswell -q interactive -t 04:00:00
```

4. Once in the shell, after loading the `root`, `cray-hdf5`, and `parallel` modules, there are a few options. The easiest thing to do is to modify `multi_parallel_run_mix.sh`, which will run various configurations one after the other. However, if there are too many, the 4-hour limit will abort whatever job is running at the time, which is probably less than ideal, so be careful with that. To run a single configuration in parallel, do
```
./parallel_run_mix.sh 4 config-run-period.yaml mixlabel
```
The output filenames of these correlations are generated as follows. First, the filename defined in `filelists.correlations.mixedevent` is used as the basis of creating the final filename; let's say it's called `mixedEvent.root`. Then, the `mixlabel` is appended, such that we get `mixedEvent_mixlabel.root`. Then, when running in parallel, each thread produces its own output, so we get something like `mixedEvent_mixlabel.rootmix_0.root`.
To run not in parallel with 10 mixed events/triggered event, you can do
```
./cpp/build/h5_mixed_event_parallel config-run-period.yaml mixlabel
```
and this will produce a file with whatever name is defined in `filelists.correlations.mixedevent`.

5. In order to combine the parallel outputs into one file and delete the individual outputs, do
```
hadd -f mixedEvent_mixlabel.root mixedEvent_mixlabel.rootmix_*
rm mixedEvent_mixlabel.rootmix_*
```

6. If the goal is to produce one final mixed event file (i.e. merging the various `mixlabels`), just do
```
hadd -f mixedEvent.root mixedEvent_*.root
```
and the individual outputs can again be removed if needed.

## Correlation calculation
The output of `same_event` and `mixed_event` is a pile of `THnSparse`s. The output of `response_matrix` is a pile of `RooUnfoldResponse`s. These are then parsed and manipulated in python using the `GammaJetCorrelation` class. Below is the workflow for using this; at some point in the future, there will be a notebook with examples of all of this added to this repo:
1. Read config files
2. Call `getAllCorr` with various parameters and filenames. This will return a dictionary of `GammaJetCorrelation` objects with the keys `[centrange][ptrange][observable]`.
3. Make whatever plots in whatever combinations. There is a `plotCorr` function defined in `GammaJetCorrelation` which takes the following keys and plots the following purity-weighted, per-event and per-trigger distributions:
  - `sesr`: same-event signal region
  - `sebr`: same-event background region
  - `mesr`: mixed-event signal region
  - `mebr`: mixed-event background region
  - `se`: same-event (`sesr` - `sebr`)
  - `me`: mixed-event (`mesr` - `mebr`)
  - `sr`: signal region (`sesr` - `mesr`)
  - `br`: background region (`sebr` - `mebr`)
  - `corr`: fully subtracted (`se` - `me` or `sr` - `br`)
