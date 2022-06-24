echo "Compiling the C++ files: "
g++ robustanalyzermain.C data_robustanalyzer.C `root-config --cflags --glibs` -o data_robustanalyzer.out

echo "Compiling successful. Begin execution."

#outname="DoubleElectronGun"
#DYToLLM50_ScoutingSkim220404_0_18.root
#outname="DYToLLM50_ScoutingSkim220404"
#outname="DoubleElectron_ScoutingSkim220411"
#outname="DoubleElectron_ScoutingSkim220517"
#outname="QCD_ScoutingSkim220510"
outname='DYToLL_ScoutingSkim220510'

#how many processing cores are available to us
nproc=6
eosdir="/eos/user/b/bgreenbe/scouting/ntuples/$outname"
#how many total files there are 
#ntot=366
#count the number of files stored in the eos directory of interest.
ntot=$( ls -l $eosdir | grep -v ^d | wc -l )
#ntot=100
nloop=$(( ntot / nproc + 1 ))
#what file to start the next hadd on
haddstart=0
for i in `seq 0 $nloop`
do
    for j in `seq 1 $nproc`
    do
        currnum=$(( $i * nproc + $j ))
        if [[ $currnum -gt $ntot ]]
        then
            break
        fi
        echo "analyzing job number $currnum"
        ./data_robustanalyzer.out $currnum $outname $(( j - 1 )) &
        eval "proc$j=$!"

        #./data_robustanalyzer.out 1 $outname &
        #proc1=$!
        #./data_robustanalyzer.out 2 $outname &
        #proc2=$!
        #./data_robustanalyzer.out 3 $outname &
        #proc3=$!
        #./data_robustanalyzer.out 4 $outname &
        #proc4=$!
        #./data_robustanalyzer.out 5 $outname &
        #proc5=$!
        #./data_robustanalyzer.out 6 $outname &
        #proc6=$!
    done

    #######while [ -d "/proc/${proc0}" -o -d "/proc/${proc1}" -o -d "/proc/${proc2}" -o -d "/proc/${proc3}" -o -d "/proc/${proc4}" -o -d "/proc/${proc5}" -o -d "/proc/${proc6}" ]
    while [ -d "/proc/${proc1}" -o -d "/proc/${proc2}" -o -d "/proc/${proc3}" -o -d "/proc/${proc4}" -o -d "/proc/${proc5}" -o -d "/proc/${proc6}" ]
    do
        echo "====Still executing===="
        sleep 1
    done
    #only hadd every 100 loops to save some time
    if [ $( expr $i % 100 ) == "99" ] 
    then
       nfiles=$(( i * nproc + nproc - haddstart ))
       python hadd_multi.py $outname $nfiles $haddstart
       rm -f hists_${outname}_*.root
       haddstart=$(( haddstart + nfiles - 1 )) 
       mv hists_${outname}.root hists_${outname}_${haddstart}.root
    fi
done

echo "Run over. Clean up and combine files."

rm data_robustanalyzer.out

#hadd -f hists_DYToLLM50.root hists_DYToLLM50_?.root
#hadd -f hists_${outname}.root hists_${outname}_*.root

#rm hists_${outname}_?.root
#rm hists_${outname}_*.root

#the final file will have the name of the last iteration of the loop. 
#hadd -f hists_${outname}.root hists_${outname}_*.root
python hadd_multi.py $outname $ntot $haddstart
#rm -f hists_${outname}_*.root
