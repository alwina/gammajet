#!/bin/bash
i=$1
nmix=$2
config=$3
mixlabel=$4
mix_start=$((i*nmix))
mix_end=$((mix_start+nmix))
echo Starting $config $mixlabel $mix_start $mix_end
./cpp/build/h5_mixed_event_parallel $config $mixlabel $mix_start $mix_end
echo Finished $config $mixlabel $mix_start $mix_end
