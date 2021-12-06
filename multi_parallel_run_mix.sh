#!/bin/bash
config=config/config-pbpb18.yaml
./parallel_run_mix.sh 4 $config 18q_int7_1
hadd -f root/mixedEventPbPb18_18q_int7_1.root root/mixedEventPbPb18_18q_int7_1.rootmix_*
rm root/mixedEventPbPb18_18q_int7_1.rootmix_*

./parallel_run_mix.sh 4 $config 18q_int7_2
hadd -f root/mixedEventPbPb18_18q_int7_2.root root/mixedEventPbPb18_18q_int7_2.rootmix*
rm root/mixedEventPbPb18_18q_int7_2.rootmix_*

./parallel_run_mix.sh 4 $config 18q_int7_3
hadd -f root/mixedEventPbPb18_18q_int7_3.root root/mixedEventPbPb18_18q_int7_3.rootmix*
rm root/mixedEventPbPb18_18q_int7_3.rootmix_*

hadd -f root/mixedEventPbPb18.root root/mixedEventPbPb18_18q_int7_*.root

./parallel_run_mix.sh 4 $config skimcent5090_18q_int7_1
./parallel_run_mix.sh 4 $config skimcent5090_18q_int7_2
./parallel_run_mix.sh 4 $config skimcent5090_18q_int7_3

hadd -f root/mixedEventPbPb18_skimcent5090.root root/mixedEventPbPb18_skimcent5090*rootmix*
rm root/mixedEventPbPb18_skimcent5090*rootmix*
