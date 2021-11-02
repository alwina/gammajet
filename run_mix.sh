#!/bin/bash
i=$1
nmix=$2
mix_start=$((i*nmix))
mix_end=$((mix_start+nmix))
echo $mix_start $mix_end
./cpp/build/h5_mixed_event config-pbpb18.yaml $mix_start $mix_end
