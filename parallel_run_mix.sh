#!/bin/bash
nmix=$1
config=$2
mixlabel=$3
maxevents=$4
gjdir=$5

if [[ -z $gjdir ]]; then
	gjdir=$(pwd)
fi

seq 0 63 | parallel --ungroup ${gjdir}/run_mix.sh {%}-1 $nmix $config $mixlabel $maxevents $gjdir
