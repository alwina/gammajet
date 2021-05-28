This is the code for the measurement of isolated photon-jet correlations in ALICE.

The starting point is NTuples generated with [this repo](https://github.com/alwina/ntuple-gj) (which is forked from [this repo](https://github.com/yslai/ntuple-gj)) from ALICE data and simulations. The code to calculate the photon purity is a refined version of what's in [this repo](https://github.berkeley.edu/alwina/photon-correlations). The code to calculate the actual correlations is based on [this repo](https://github.com/ftoralesacosta/GH_Correlations).

## Configuration files
`yaml` is used for the configuration file. It contains nearly all of the information needed to run the analysis (after the ntuples are made). It's slightly recursive in the sense that the photon purity is part of the configuration, but the file is also used to calculate the photon purity. So after calculating the purity with a given configuration, one has to then go back and re-edit it to update the purity.

## Purity calculation
To calculate the purity, the ntuples are first converted to CSVs. Those CSVs are then read into pandas dataframes, from which template fits are done. The `DataframeCollection` and `PhotonPurity` classes handle most of this.

## Correlation calculation
The same-event correlations, mixed-event correlations, and RooUnfold responses are done with C++ scripts. `cmake` (minimum version 3.10.2) is used to compile these scripts, which have the following requirements
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

Once built, run any of the scripts with the path to the config file as an argument. For example,

```
cpp/build/same_event config-pbpb15.yaml
```

The output of `same_event` and `mixed_event` is a pile of `THnSparse`s. The output of `response_matrix` is a pile of `RooUnfoldResponse`s. These are then parsed and manipulated in python using the `GammaJetCorrelation` class.