#!/bin/bash
nmix=$1
config=$2
mixlabel=$3
seq 0 64 | parallel --ungroup ./run_mix.sh {%} $nmix $config $mixlabel
