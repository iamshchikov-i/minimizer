#!/bin/bash

SAMPLE_PATH=/common/home/yamshikov_i/build_minimizer/bin

export NCOUNT=1
export PRIMARY_ITERATIONS_TYPE=involve

#$SAMPLE_PATH/sample dims=1 epsPar=0.001 rPar=2.3 Nmax=100000000 epsErr=0.01
$SAMPLE_PATH/sample dims=1 epsPar=0.001 rPar=2.3 Nmax=100000000 epsErr=0.01 useThreads=1 threadsNum=4
#./sample dims=1 epsPar=0.001 rPar=2.0 Nmax=100000000 epsErr=0.01 useThreads=1 threadsNum=4 > th_4.txt
#./sample dims=1 epsPar=0.001 rPar=2.0 Nmax=100000000 epsErr=0.01 useThreads=1 threadsNum=8 > th_8.txt