#!/bin/bash
f="$(basename -- $1)";
outputdir=$1"/results/performance";
mkdir $1"/results";
mkdir $outputdir;

./bin/run_experiments btch $1 $outputdir
./bin/run_experiments btch_l $1 $outputdir
./bin/run_experiments ftch_tcpd $1 $outputdir
./bin/run_experiments ftch_rtpd $1 $outputdir
./bin/run_experiments ftch_l $1 $outputdir
./bin/run_experiments ftch_cpd $1 $outputdir

./bin/run_experiments bstch $1 $outputdir
./bin/run_experiments fstch_tcpd $1 $outputdir
./bin/run_experiments bmtch $1 $outputdir
./bin/run_experiments fmtch_tcpd $1 $outputdir
