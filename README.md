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
In the interest of not rerunning the entire ntuplizer and also not modifying 100+ GB files in-place, in order to use information that is derived from the ntuples but not saved in them, an auxiliary file gets created. At the moment, this includes `cluster_5x5all`, `jet_ak02tpc`, `jet_ak02its`, and `jet_charged_truth_ak02` information, but it can include anything that is better off calculated once to be used later. To create an aux file for an existing ntuple.root, within a python3 console (with ROOT), run the following:
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
- [RooUnfold](https://gitlab.cern.ch/RooUnfold/RooUnfold) (process_gjmc only), for which `$ROOUNFOLD_DIR` has to be defined

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
The output of `same_event` and `mixed_event` is a pile of `THnSparse`s. These are then parsed and manipulated in python using the `GammaJetCorrelation2D` class. Below is the workflow for using this:
1. Read config files
2. Call `getAllCorr` with various parameters and filenames. This will return a dictionary of `GammaJetCorrelation2D` objects with the keys `[centrange][ptrange][observable]`.
3. Generally, you'll want to use the `corrth1` object, as this already has the bin widths divided out. Otherwise, you can take various projections (usually the X projection, which corresponds to the observable, rather than the Y projection, which is the jet pT projection), but keep in mind that the bin width has not been divided out if you do so.
4. `plotTH1` (in `utils.py`) is pretty cool. Check it out.

## Response matrices and (un)folding
`RooUnfoldResponse`s and their corresponding `TH4`s are created by running
```
./cpp/build/process_gjmc config/config-run-period.yaml
```

This will use the ntuples defined in `filelists.ntuples.gjmc` and will produce two files, one defined by `filelists.correlations.responsematrix` and the other defined by `filelists.correlations.gjmc`. The former contains two types of objects: `RooUnfoldResponse` objects and `THnSparse`s. These 4D histograms have the observable (`deltaphi` or `ptratio`) as the first two axes and the jet pT as the second two axes. By convention, reco is the first of the pair and truth is the second. So, for example, a TH4 might have the axes

1. deltaphi reco
2. deltaphi truth
3. jetpt reco
4. jetpt truth

Since there's no good way to get a 4D object out of a `RooUnfoldResponse`, we separately save the histograms for visualization purposes.

The `Response` and `Hist` objects are written out with a rough numbering scheme corresponding to which centrality and photon pT bin they correspond to. For example, `deltaphijetptResponse21` refers to the 2nd centrality bin and the 1st photon pT bin. The bins are defined in the configuration as `centralityranges` and `clusterptranges`, **both of which must be nonoverlapping**.

The `Unfolder2D` class defined in `unfolder.py` is essentially a wrapper for manipulating the `RooUnfoldResponse`. Because of how `SetRangeUser` and `Projection` interact in ROOT, you may have to create a new `TH2` that re-adds the jet pT bins that get cut and thrown out. A quick way to do this is something like

```
newth2 = ROOT.TH2F('', '', 6, 0, np.pi, 12, -10, 50)
for binx in range(1, 7):
    for biny in range(1, 9):
        newth2.SetBinContent(binx, biny + 4, oldth2.GetBinContent(binx, biny))
```

To use the `Unfolder2D` object, you can do something like the following:
```
uf = Unfolder2D()
uf.setMeasuredTH2(datath2)
uf.setResponseMatrix(roounfoldresponse)
uf.setTH4(responseth4)
```

Then, to do Bayesian unfolding with, say, 1, 3, and 5 iterations, just do
```
uf.unfoldBayes(1, 3, 5)
```
Note that this is not a list; if you have a list like `l = [1, 3, 5]`, then do `uf.unfoldBayes(*l)` to unpack it. Similarly, there is also `unfoldSvd` which works the same way. Once one (or both) of these is run, you can access the unfolded histograms with `uf.unfoldedBayesTH2[niterations]` or `uf.unfoldedSvdTH2[k]` as appropriate.

There are then some plots that can be automatically generated. First, `plotAllUnfolded` does what it sounds like. If you only want to plot a subset of the different numbers of iterations, you can do so by setting `bayesKeys=[1, 3]` (or, equivalently, `svdKeys=[3, 4]`).

Another plot is `plotBayesConvergence`, which divides the (n+1)th iteration by the nth iteration. Note that n+1 and n refer to the numbers passed into `unfoldBayes`, so in this example, it would do 1 / 0 (where 0 is the measured distribution), 3 / 1 and 5 / 3.

Then there's `shapeClosureTestMCDataRatio`, which applies the shape closure test and `statisticalClosureTest`. Ask someone who does jet substructure exactly what those test. Both plot the outcomes of these tests.

For all of these plots, a figure is not generated (so do `plt.figure` first) and legends, axis labels, etc are not set. The labels for the legend are automatically generated, but the legend itself is not automatically placed onto the plot.

## More about process_gjmc.cpp
One thing about generating response matrices is that you have to match reconstructed and truth objects; in this case, this means matching reconstructed and truth jets. In `process_gjmc.cpp` is a function called `matchJetsInEvent`; this is where the jet matching takes place, so any changes needed to the jet matching would go there.

In addition to the response matrices, another file is generated. This contains a few kinds of objects which are sometimes useful:
- `hTrig` and `hCorr` objects for dealing with the PYTHIA on either the truth or reco level
- Photon pT and phi resolution
- Jet pT and phi resolution
- Back-to-back (deltaphi > 7pi/8) jet pT and phi resolution

If you need to add anything that gets generated from processing the embedded gamma-jet MC that isn't some kind of response matrix, it's sensible to throw them into this file or to create a new file for whatever it is you're looking at.