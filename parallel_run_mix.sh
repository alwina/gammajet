#!/bin/bash
nmix=$1
seq 0 64 | parallel --ungroup ./run_mix.sh {%} $nmix
