#!/bin/bash

SAMPLE_PATH=/common/home/yamshikov_i/build_minimizer/bin

export NCOUNT=5000000

export PRIMARY_ITERATIONS_TYPE=none
$SAMPLE_PATH/sample dims=1 epsPar=0.001 rPar=2.3 Nmax=100000000 epsErr=0.01 > ./res/none_th_1.txt
$SAMPLE_PATH/sample dims=1 epsPar=0.001 rPar=2.3 Nmax=100000000 epsErr=0.01 useThreads=1 threadsNum=2 > ./res/none_th_2.txt 
$SAMPLE_PATH/sample dims=1 epsPar=0.001 rPar=2.3 Nmax=100000000 epsErr=0.01 useThreads=1 threadsNum=4 > ./res/none_th_4.txt 
$SAMPLE_PATH/sample dims=1 epsPar=0.001 rPar=2.3 Nmax=100000000 epsErr=0.01 useThreads=1 threadsNum=8 > ./res/none_th_8.txt

export PRIMARY_ITERATIONS_TYPE=separate
$SAMPLE_PATH/sample dims=1 epsPar=0.001 rPar=2.3 Nmax=100000000 epsErr=0.01 > ./res/sep_th_1.txt
$SAMPLE_PATH/sample dims=1 epsPar=0.001 rPar=2.3 Nmax=100000000 epsErr=0.01 useThreads=1 threadsNum=2 > ./res/sep_th_2.txt 
$SAMPLE_PATH/sample dims=1 epsPar=0.001 rPar=2.3 Nmax=100000000 epsErr=0.01 useThreads=1 threadsNum=4 > ./res/sep_th_4.txt 
$SAMPLE_PATH/sample dims=1 epsPar=0.001 rPar=2.3 Nmax=100000000 epsErr=0.01 useThreads=1 threadsNum=8 > ./res/sep_th_8.txt

export PRIMARY_ITERATIONS_TYPE=involve
$SAMPLE_PATH/sample dims=1 epsPar=0.001 rPar=2.3 Nmax=100000000 epsErr=0.01 > ./res/inv_th_1.txt
$SAMPLE_PATH/sample dims=1 epsPar=0.001 rPar=2.3 Nmax=100000000 epsErr=0.01 useThreads=1 threadsNum=2 > ./res/inv_th_2.txt 
$SAMPLE_PATH/sample dims=1 epsPar=0.001 rPar=2.3 Nmax=100000000 epsErr=0.01 useThreads=1 threadsNum=4 > ./res/inv_th_4.txt 
$SAMPLE_PATH/sample dims=1 epsPar=0.001 rPar=2.3 Nmax=100000000 epsErr=0.01 useThreads=1 threadsNum=8 > ./res/inv_th_8.txt