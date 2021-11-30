#!/bin/bash
config=config/config-pbpb18.yaml
./parallel_run_mix.sh 4 $config 18q_int7_1
./parallel_run_mix.sh 4 $config 18q_int7_2
./parallel_run_mix.sh 4 $config 18q_int7_3
./parallel_run_mix.sh 4 $config skimcent5090_18q_int7_1
./parallel_run_mix.sh 4 $config skimcent5090_18q_int7_2
./parallel_run_mix.sh 4 $config skimcent5090_18q_int7_3
