#!/bin/sh

cd $TMPDIR
mkdir Job_$1_$2
echo "============ Created the master run folder ============"

cd Job_$1_$2
cd /user/asahasra/CMSSW_12_4_4/src/
eval `scramv1 runtime -sh`
cd $TMPDIR/Job_$1_$2
pwd

cp /user/asahasra/CMSSW_12_4_4/src/Run3ScoutingAnalysisTools/PlotAndStudyNtuples/robustanalyzer* ./

echo "============ Copied the cpp run files file. Listing directory contents ============"
ls
echo "======================= End listing the directory contents ========================"
echo ""

echo "Compiling the C++ files: "
g++ robustanalyzermain_forCondor.C robustanalyzer.C `root-config --cflags --glibs` -o robustanalyzer.out

echo "Compiling successful. Begin execution."

./robustanalyzer.out $2 49 "/pnfs/iihe/cms/store/user/asahasra/ScoutingPFRun3/Scouting2022B_crabRunSkim220823_asahasra/220823_185023/0001/scoutingNTuple_*.root" "hists_data_$1_$2.root"

echo ""
echo "Run over. Copy the output files."
ls

cp hists_data_$1_$2.root $CMSSW_BASE/src/Run3ScoutingAnalysisTools/PlotAndStudyNtuples/hists_data_$1_$2.root
