#!/bin/bash
#cd /user/asahasra/CMSSW_12_1_0_pre4/src/
#cd /afs/cern.ch/work/b/bgreenbe/public/CMSSW_12_3_0_pre6/src/
#cd /afs/cern.ch/work/b/bgreenbe/public/CMSSW_12_3_2/src/
cd /afs/cern.ch/work/b/bgreenbe/public/CMSSW_12_4_8/src/
eval `scram runtime -sh`

cd ./Run3ScoutingAnalysisTools/Analysis/test
#cd ../../CMSSW_12_3_0_pre6/src/Run3ScoutingAnalysisTools/Analysis/test

#ls
echo "=========== BEGINNING TO RUN THE SCRIPT ================"

export X509_USER_PROXY=$2
voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/b/bgreenbe/x509up_u104084

#testname="DoubleElectron_ScoutingSkim220410"
#testname="DoubleElectron_ScoutingSkim220411"
#testname="DoubleElectron_ScoutingSkim220517"
#testname="QCD_ScoutingSkim220510"
#testname="BuToKee_ScoutingSkim220909"
#testname="DYToLL_ScoutingSkim220510"
#testname="DYToLL_ScoutingSkim220927"
#testname="EtaTo2Mu2E_ScoutingSkim221006"
#testname="EtaTo2Mu2E_ScoutingSkimSeeded221208"
testname="EtaTo2Mu2E_ScoutingSkim230210"
#mkdir -p /eos/user/b/bgreenbe/scouting/ntuples/$testname
jobnum=$1
echo "script number $jobnum"
#how many files to do per job
nfiles=5 #8 #4821
##MAKE SURE TO CHANGE THIS BACK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
startnum=$(( jobnum * nfiles + 1 ))
endnum=$(( startnum + nfiles - 1 ))
echo "startnum $startnum , endnum $endnum, dirnum $dirnum"
#for num in $(ls /pnfs/iihe/cms/store/user/asahasra/DYToLL_M-50_TuneCP5_14TeV-pythia8/ScoutingSkim220201_DYToLLM50Run3Summer21_asahasra/220201_162520/0000/ | grep -o -E '[0-9]+')
#do
#    if [[ $num -eq 2022 ]]; then
#    continue
#    fi
#    if [[ $(($num%10)) -ne $1 ]]; then
#    continue
#    fi
#    cmsRun ScoutingNanoAOD_cfg.py isMC=True useWeights=False GlobalTagMC=112X_mcRun3_2021_realistic_v16 output=DYToLLM50_ScoutingSkim220127_$1_$num.root inputFile=file:/pnfs/iihe/cms/store/user/asahasra/DYToLL_M-50_TuneCP5_14TeV-pythia8/ScoutingSkim220201_DYToLLM50Run3Summer21_asahasra/220201_162520/0000/HLT2022_HLT_$num.root
#done
#for num in `seq 1 999`
for num in `seq $startnum $endnum`
do
    echo "starting num $num"
    dirnum=$(( num / 1000 ))
    #fulldir="/eos/user/b/bgreenbe/scouting/DoubleElectron_Pt-1To300-gun/ScoutingSkim220409_DoubleElectronGunRun3Summer21_bgreenbe/220408_150915/000$dirnum"
    #fulldir="/eos/user/b/bgreenbe/scouting/DoubleElectron_Pt-1To300-gun/ScoutingSkim220517_DoubleElectron_bgreenbe/220517_143258/000$dirnum"
    #fulldir="/eos/user/b/bgreenbe/scouting/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220510_QCDPt20-30_bgreenbe/220510_212445/000$dirnum"
    #fulldir="/eos/user/b/bgreenbe/scouting/BuToKee_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/ScoutingSkim220909_BuToKee_bgreenbe/220909_092238/000$dirnum"
    #fulldir="/eos/user/b/bgreenbe/scouting/CRAB_UserFiles/ScoutingSkim220919_EtaTo2Mu2E_bgreenbe/220919_101711/000$dirnum"
    #fulldir="/eos/user/b/bgreenbe/scouting/CRAB_UserFiles/ScoutingSkim220922_EtaTo2Mu2E_bgreenbe/220922_155408/000$dirnum"
    #fulldir="/eos/user/b/bgreenbe/scouting/CRAB_UserFiles/ScoutingSkim220923_EtaTo2Mu2E_bgreenbe/220923_074250/000$dirnum"
    #fulldir="/eos/user/b/bgreenbe/scouting/CRAB_UserFiles/ScoutingSkim221006_EtaTo2Mu2E_bgreenbe/221006_105951/000$dirnum"
    #fulldir="/eos/user/b/bgreenbe/scouting/CRAB_UserFiles/ScoutingSkim221208_EtaTo2Mu2E_bgreenbe/221208_170905/000$dirnum"
    fulldir="/eos/user/b/bgreenbe/scouting/CRAB_UserFiles/ScoutingSkim230210_EtaTo2Mu2E_bgreenbe/230210_105846/000$dirnum"
    #fulldir="/eos/user/b/bgreenbe/scouting/DYToLL_M-4To50_TuneCP5_14TeV-pythia8/ScoutingSkim220510_DYToLL_bgreenbe/220510_161328/000$dirnum"
    #fulldir="/eos/user/b/bgreenbe/scouting/DYToLL_M-4To50_TuneCP5_14TeV-pythia8/ScoutingSkim220927_DYToLL_bgreenbe/220927_113037/000$dirnum"
    #fulldir="/eos/user/b/bgreenbe/scouting/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220429_QCDPt20To30_bgreenbe/220429_000456/000$dirnum"
    #fulldir="/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/"
    #fulldir="/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220329_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220330_115058/0000/"
    #see if the file can be found
    fullpath="$fulldir/outputScoutingPF_$num.root"
    #xrdfs cms-xrd-global.cern.ch locate $fullpath
    ls $fullpath
#    xrdfs cms-xrd-global.cern.ch locate /store/user/asahasra/DoubleElectron_Pt-1To300-gun/ScoutingSkim220329_DoubleElectronGunRun3Summer21_asahasra/220330_092202/0000/outputScoutingPF_$num.root
    if [[ $? -eq 0 ]] 
    then
        comm="cmsRun ScoutingNanoAOD_cfg.py isMC=True useWeights=False GlobalTagMC=112X_mcRun3_2021_realistic_v16 output=/eos/user/b/bgreenbe/scouting/ntuples/$testname/${testname}_0_$num.root inputFile=file:$fulldir/outputScoutingPF_$num.root"
        #comm="cmsRun ScoutingNanoAOD_cfg.py isMC=True useWeights=False GlobalTagMC=112X_mcRun3_2021_realistic_v16 output=DoubleElectron_ScoutingSkim220506_0_$num.root inputFile=file:/afs/cern.ch/work/a/asahasra/public/tempfiles/outputScoutingPF_1.root"
        #comm="cmsRun ScoutingNanoAOD_cfg.py isMC=True useWeights=False GlobalTagMC=112X_mcRun3_2021_realistic_v16 output=/eos/user/b/bgreenbe/scouting/ntuples/QCD_ScoutingSkim220429/QCD_ScoutingSkim220429_0_$num.root inputFile=file:$fullpath"
        echo $comm
        comm2="echo trying file number $num"
        eval $comm2
        eval $comm
    else
        #echo "file number $num not found. Breaking!"
        echo "file number $num not found. Continuing."
        #break
    fi
done

echo "loop over!"

#for num in $(ls /pnfs/iihe/cms/store/user/asahasra/DYToLL_M-4To50_TuneCP5_14TeV-pythia8/ScoutingSkim220127_DYToLLM4To5Run3Summer21_asahasra/220128_130307/0000/ | grep -o -E '[0-9]+')
#do
#    if [[ $num -eq 2022 ]]; then
#    continue
#    fi
#    if [[ $(($num%10)) -ne $1 ]]; then
#    continue
#    fi
#    cmsRun ScoutingNanoAOD_cfg.py isMC=True useWeights=False GlobalTagMC=112X_mcRun3_2021_realistic_v16 output=DYToLLM4To50_ScoutingSkim220127_$1_$num.root inputFile=file:/pnfs/iihe/cms/store/user/asahasra/DYToLL_M-4To50_TuneCP5_14TeV-pythia8/ScoutingSkim220127_DYToLLM4To5Run3Summer21_asahasra/220128_130307/0000/HLT2022_HLT_$num.root
#done
#
#for num in $(ls /pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/ | grep -o -E '[0-9]+')
#do
#    if [[ $num -eq 2022 ]]; then
#    continue
#    fi
#    if [[ $(($num%10)) -ne $1 ]]; then
#    continue
#    fi
#    cmsRun ScoutingNanoAOD_cfg.py isMC=True useWeights=False GlobalTagMC=112X_mcRun3_2021_realistic_v16 output=QCDPt20To30EmEnriched_ScoutingSkim220127_$1_$num.root inputFile=file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-20To30_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt20To30EmEnrichedRun3Summer21_asahasra/220128_130326/0000/HLT2022_HLT_$num.root
#done
#
#for num in $(ls /pnfs/iihe/cms/store/user/asahasra/QCD_Pt-30To50_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt30To50EmEnrichedRun3Summer21_asahasra/220128_130341/0000/ | grep -o -E '[0-9]+')
#do
#    if [[ $num -eq 2022 ]]; then
#    continue
#    fi
#    if [[ $(($num%10)) -ne $1 ]]; then
#    continue
#    fi
#    cmsRun ScoutingNanoAOD_cfg.py isMC=True useWeights=False GlobalTagMC=112X_mcRun3_2021_realistic_v16 output=QCDPt30To50EmEnriched_ScoutingSkim220127_$1_$num.root inputFile=file:/pnfs/iihe/cms/store/user/asahasra/QCD_Pt-30To50_EMEnriched_TuneCP5_14TeV-pythia8/ScoutingSkim220127_QCDPt30To50EmEnrichedRun3Summer21_asahasra/220128_130341/0000/HLT2022_HLT_$num.root
#done
#
#UNCOMMENT BELOW LINE!!
#hadd -f DYToLLM50_ScoutingSkim220127_$1.root DYToLLM50_ScoutingSkim220127_$1_*.root
#hadd -f DYToLLM4To50_ScoutingSkim220127_$1.root DYToLLM4To50_ScoutingSkim220127_$1_*.root
#hadd -f QCDPt20To30EmEnriched_ScoutingSkim220127_$1.root QCDPt20To30EmEnriched_ScoutingSkim220127_$1_*.root
#hadd -f QCDPt30To50EmEnriched_ScoutingSkim220127_$1.root QCDPt30To50EmEnriched_ScoutingSkim220127_$1_*.root

#UNCOMMENT BELOW LINE!!!
#rm DYToLLM50_ScoutingSkim220127_$1_*.root
#rm DYToLLM4To50_ScoutingSkim220127_$1_*.root
#rm QCDPt20To30EmEnriched_ScoutingSkim220127_$1_*.root
#rm QCDPt30To50EmEnriched_ScoutingSkim220127_$1_*.root
#
echo "================= END OF SCRIPT RUN ===================="
#ls
