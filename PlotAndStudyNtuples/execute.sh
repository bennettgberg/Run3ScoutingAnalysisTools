echo "Compiling the C++ files: "
g++ robustanalyzermain.C robustanalyzer.C `root-config --cflags --glibs` -o robustanalyzer.out

echo "Compiling successful. Begin execution."

./robustanalyzer.out 0 7 &
proc0=$!
./robustanalyzer.out 1 7 &
proc1=$!
./robustanalyzer.out 2 7 &
proc2=$!
./robustanalyzer.out 3 7 &
proc3=$!
./robustanalyzer.out 4 7 &
proc4=$!
./robustanalyzer.out 5 7 &
proc5=$!
./robustanalyzer.out 6 7 &
proc6=$!
./robustanalyzer.out 7 7 &
proc7=$!

while [ -d "/proc/${proc0}" -o -d "/proc/${proc1}" -o -d "/proc/${proc2}" -o -d "/proc/${proc3}" -o -d "/proc/${proc4}" -o -d "/proc/${proc5}" -o -d "/proc/${proc6}" -o -d "/proc/${proc7}" ]
#while [ -d "/proc/${proc0}" -o -d "/proc/${proc1}" -o -d "/proc/${proc2}" -o -d "/proc/${proc3}" ]
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

rm robustanalyzer.out

hadd -f hists_data.root hists_data_?.root
#hadd -f hists_DoubleElectronGunPt1To300.root hists_DoubleElectronGunPt1To300_?.root
#hadd -f hists_DoubleElectronGunPt1To300Old.root hists_DoubleElectronGunPt1To300Old_?.root
#hadd -f hists_QCDPt20To30EmEnriched.root hists_QCDPt20To30EmEnriched_?.root
#hadd -f hists_QCDPt30To50EmEnriched.root hists_QCDPt30To50EmEnriched_?.root

rm hists_data_?.root
#rm hists_DoubleElectronGunPt1To300_?.root
#rm hists_DoubleElectronGunPt1To300Old_?.root
#rm hists_QCDPt20To30EmEnriched_?.root
#rm hists_QCDPt30To50EmEnriched_?.root
