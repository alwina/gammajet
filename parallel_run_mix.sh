#!/bin/bash
nmix=$1
config=$2
mixlabel=$3
seq 0 63 | parallel --ungroup ./run_mix.sh {%}-1 $nmix $config $mixlabel
