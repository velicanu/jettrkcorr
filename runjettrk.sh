if [ $# -ne 24 ]
then
  echo "Usage: ./runcorr.sh <condor_iter> <trackqual> <file-list> <tag> <nmin> <nmax> <pttrigmin> <pttrigmax> <ptassmin> <ptassmax> <leadingjetptlow> <leadingjetpthigh> <subleadingjetptlow> <subleadingjetpthigh> <centmin> <centmax> <ajmin> <ajmax> <dotrkcorr> <mccommand> <jetetamin> <jetetamax> <whichjet> <bklist>"
  exit 1
fi


# echo tar -xf tables.tar
# echo tar -xf tables.tar | bash

echo tar -xzf corrFilePbPb_20140429.tar.gz
echo tar -xzf corrFilePbPb_20140429.tar.gz | bash
echo tar -xzf mix.tgz
echo tar -xzf mix.tgz | bash


echo | awk -v i=$1 -v trkq=$2 -v flist=$3 -v tag=$4 -v nmin=$5 -v nmax=$6 -v pttrigmin=$7 -v pttrigmax=$8 -v ptassmin=$9 -v ptassmax=${10} -v leadingjetptlow=${11} -v leadingjetpthigh=${12} -v subleadingjetptlow=${13} -v subleadingjetpthigh=${14} -v centmin=${15} -v centmax=${16} -v ajmin=${17} -v ajmax=${18} -v dotrkcorr=${19} -v mccommand=${20} -v jetetamin=${21} -v jetetamax=${22} -v whichjet=${23} -v bklist=${24}  '{print "./runjettrk.exe "i" "trkq" "flist" "tag" "nmin" "nmax" "pttrigmin" "pttrigmax" "ptassmin" "ptassmax" "leadingjetptlow" "leadingjetpthigh" "subleadingjetptlow" "subleadingjetpthigh" "centmin" "centmax" "ajmin" "ajmax" "dotrkcorr" "mccommand" "jetetamin" "jetetamax" "whichjet" "bklist" "}' 
echo | awk -v i=$1 -v trkq=$2 -v flist=$3 -v tag=$4 -v nmin=$5 -v nmax=$6 -v pttrigmin=$7 -v pttrigmax=$8 -v ptassmin=$9 -v ptassmax=${10} -v leadingjetptlow=${11} -v leadingjetpthigh=${12} -v subleadingjetptlow=${13} -v subleadingjetpthigh=${14} -v centmin=${15} -v centmax=${16} -v ajmin=${17} -v ajmax=${18} -v dotrkcorr=${19} -v mccommand=${20} -v jetetamin=${21} -v jetetamax=${22}  -v whichjet=${23} -v bklist=${24}  '{print "./runjettrk.exe "i" "trkq" "flist" "tag" "nmin" "nmax" "pttrigmin" "pttrigmax" "ptassmin" "ptassmax" "leadingjetptlow" "leadingjetpthigh" "subleadingjetptlow" "subleadingjetpthigh" "centmin" "centmax" "ajmin" "ajmax" "dotrkcorr" "mccommand" "jetetamin" "jetetamax" "whichjet" "bklist" "}' | bash
# echo | awk -v tag=$4 -v user=$USER '{print "mkdir -p /net/hidsk0001/d00/scratch/"user"/corrhists/HIRun2013/merged"}'
# echo | awk -v tag=$4 -v user=$USER '{print "mkdir -p /net/hidsk0001/d00/scratch/"user"/corrhists/HIRun2013/merged"}' | bash

echo | awk -v tag=$4 -v user=$USER -v now=${24} '{print "mv "tag"*.root /mnt/hadoop/cms/store/user/"user"/unmerged/"}'
echo | awk -v tag=$4 -v user=$USER -v now=${24} '{print "mv "tag"*.root /mnt/hadoop/cms/store/user/"user"/unmerged/"}' | bash

echo rm *.root corrFilePbPb_20140429.tar.gz *.txt
echo rm -r tables | bash

# echo | awk -v tag=$4 '{print "rm "tag"*.root"}'
# echo | awk -v tag=$4 '{print "rm "tag"*.root"}' | bash

echo "job done successfully"
