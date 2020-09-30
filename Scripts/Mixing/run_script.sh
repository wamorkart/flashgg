#!/bin/sh -e

JOBID=$1;
INPUTDIR=$2;
INPUTSTRING=$3;
OUTPUTDIR=$4;

cd $INPUTDIR/

echo -e "evaluate"
eval `scramv1 ru -sh`

echo -e "Compute SoB";
# python runSoBOptimization.py -s "${INPUTSTRING}"
python /afs/cern.ch/work/t/twamorka/flashgg_16aug2020/CMSSW_10_6_8/src/flashgg/Scripts/Mixing/H4GTreeMixing.py  "${INPUTSTRING}"

echo -e "DONE";
