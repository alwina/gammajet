#!/bin/bash
i=$1
nmix=$2
config=$3
mixlabel=$4
maxevents=$5
gjdir=$6

if [[ -z $gjdir ]]; then
	gjdir=$(pwd)
fi

mix_start=$((i*nmix))
mix_end=$((mix_start+nmix))
${gjdir}/cpp/build/mixed_event $config $mixlabel $mix_start $mix_end $maxevents
