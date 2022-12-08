#!/bin/sh

cd $TMPDIR
mkdir Job_$1_$2
ls
echo "============ Created the master run folder ============"

cd Job_$1_$2
cd /user/asahasra/CMSSW_12_4_4/src/
eval `scramv1 runtime -sh`
cd -

cp /user/asahasra/CMSSW_12_4_4/src/Run3ScoutingAnalysisTools/Analysis/test/ScoutingNanoAOD_cfg.py ./

ls
echo "============ Copied the python config file ============"
echo $3
#cp $3 ./
echo "============ Listing the directory contents pre cmsRun ============"
ls

echo "============ Running ============"

cmsRun ScoutingNanoAOD_cfg.py isMC=True GlobalTagMC=124X_mcRun3_2022_realistic_v10 output=QCDPt80To120_ntuple_$1_$2.root inputFile=file:$3

echo "============ Listing the directory contents post cmsRun ============"
ls

cp QCDPt80To120_ntuple_$1_$2.root /user/asahasra/CMSSW_12_4_4/src/Run3ScoutingAnalysisTools/Analysis/test
echo "============ All processes completed. Exiting now.  ============"
