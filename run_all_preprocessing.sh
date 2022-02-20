#!/bin/bash
f="$(basename -- $1)";
outputdir=$1"/results/preprocessing";
mkdir $1"/results";
mkdir $outputdir;

#./bin/run_ordering $1"/"$f".tpgr" $1"/"$f".btch" 32;
#
#hour=0
#while [ $hour -lt  24 ]
#do
#./bin/construct_hourly_tch $1"/"$f $hour 32
#((hour++))
#done
#
#
#./bin/construct_cpd $1"/"$f".tpgr" $outputdir;
#echo  ./bin/run_ordering $1"/"$f".tpgr" $1"/"$f".btch" 12;
#

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
#    hour=0
#    while [ $hour -lt  24 ]
#    do
#    ./bin/construct_hourly_tch $1"/"$f $hour 32
#    ((hour++))
#    done
#    #
#    ./bin/split_data  $1"/"$f".tpgr" 5
#    ./bin/construct_hourly_tch $1"/"$f 0 32
#    ./bin/run_ordering $1"/"$f".tpgr" $1"/"$f".btch" 32;
#    ./bin/generate_landmarks $1"/"$f".tpgr" 4
#    ./bin/generate_landmarks $1"/"$f".tpgr" 8
#    ./bin/generate_landmarks $1"/"$f".tpgr" 12
#    ./bin/generate_landmarks $1"/"$f".tpgr" 16

    ./bin/construct_fw_tch_cpd $1"/"$f".tpgr" $1"/"$f".btch" $outputdir;
#    ./bin/construct_bw_tch_cpd $1"/"$f".tpgr" $1"/"$f".btch" $outputdir;
#    ./bin/generate_queries $1"/"$f".tpgr" 10000
#    ./bin/construct_reverse_fw_tch_cpd $1"/"$f".tpgr" $1"/"$f".btch" $outputdir;
elif [[ "$OSTYPE" == "darwin"* ]]; then
# APPLE Machine
    # split the time dependent  data for building stch and mtch
    ./bin/split_data  $1"/"$f".tpgr"
    # build stch and mtch, input are number of core used
    ./bin/construct_stch_mtch $1"/"$f 12
    # build tch, input are number of core used
    ./bin/run_ordering $1"/"$f".tpgr" $1"/"$f".btch" 12;
    # build landmark heuristic
    ./bin/generate_landmarks $1"/"$f".tpgr" 4
    ./bin/generate_landmarks $1"/"$f".tpgr" 8
    ./bin/generate_landmarks $1"/"$f".tpgr" 12
    ./bin/generate_landmarks $1"/"$f".tpgr" 16

    # build forward tcpd
    ./bin/construct_fw_tch_cpd $1"/"$f".tpgr" $1"/"$f".btch" $outputdir;
    # build backward tcpd
    ./bin/construct_bw_tch_cpd $1"/"$f".tpgr" $1"/"$f".btch" $outputdir;
    # generate random queries
    ./bin/generate_queries $1"/"$f".tpgr" 10000

    # we only build RTPD for target row
    ./bin/construct_reverse_fw_tch_cpd $1"/"$f".tpgr" $1"/"$f".btch" $outputdir;
fi

###

