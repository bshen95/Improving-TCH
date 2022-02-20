#!/bin/bash
f="$(basename -- $1)";
outputdir=$1"/results/performance";
mkdir $1"/results";
mkdir $outputdir;

./bin/run_experiments tch $1 $outputdir
./bin/run_experiments tch_l_star $1 $outputdir
./bin/run_experiments tch_cpd $1 $outputdir
./bin/run_experiments rev_tch_cpd $1 $outputdir
./bin/run_experiments tch_cpd_l $1 $outputdir
./bin/run_experiments tch_ori_cpd $1 $outputdir

./bin/run_experiments htch $1 $outputdir
./bin/run_experiments htch_cpd $1 $outputdir
./bin/run_experiments ts_htch $1 $outputdir
./bin/run_experiments ts_htch_cpd $1 $outputdir
