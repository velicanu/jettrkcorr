if [ $# -ne 23 ]
then
  echo "Usage: ./psort.sh <trackqual> <file-list> <tag> <nmin> <nmax> <pttrigmin> <pttrigmax> <ptassmin> <ptassmax> <leadingjetptlow> <leadingjetpthigh> <subleadingjetptlow> <subleadingjetpthigh> <centmin> <centmax> <ajmin> <ajmax> <dotrkcorr> <mccommand> <jetetamin> <jetetamax> <whichjet> <bklist>"
  echo "akPu3Calo akPu3PF akVs3Calo akVs3PF"
  exit 1
fi

now="$3_$4_$5_$6_$7_$8_$9_${10}_${11}_${12}_${13}_${14}_${15}_${16}_${17}_${18}_${19}_${20}_${21}_$(date +"%Y_%m_%d__%H_%M_%S")"
len=`wc -l $2 | awk '{print $1}'`
njobs=$((len))

mkdir $now
cp $2 $now
cp ${23} $now
cp runjettrk.sh $now
# cp tables.tar $now
cp corrFilePbPb_20140429.tar.gz $now
cat runjettrk.condor | sed "s@log_flag@$now@g" | sed "s@dir_flag@$PWD/$now@g" | sed "s@user_flag@$USER@g" |  sed "s@arglist@$1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18} ${19} ${20} ${21} ${22} ${23}@g" | sed "s@transfer_filelist@$2,${23},runjettrk.exe,corrFilePbPb_20140429.tar.gz@g" | sed "s@njobs@$njobs@g" > $now/runjettrk.condor

NAME="runjettrk.C"
g++ $NAME $(root-config --cflags --libs) -Werror -Wall -O2 -o "${NAME/%.C/}.exe"
cp runjettrk.exe $now

fajmin=`bc <<< "${16}*10"`
fajmax=`bc <<< "${17}*10"`
iajmin=${fajmin%.*}
iajmax=${fajmax%.*}


echo -e "\E[32m"
echo "hadd -f merged/$3_${22}_trkqaul$1_nmin$4_nmax$5_cmin${14}_cmax${15}_aj${iajmin}_${iajmax}_tptmin$6_tptmax$7_aptmin$8_aptmax${9}.root unmerged/$3_${22}_trkqaul$1_nmin$4_nmax$5_cmin${14}_cmax${15}_aj${iajmin}_${iajmax}_tptmin$6_tptmax$7_aptmin$8_aptmax${9}_*.root"
echo "hadd -f merged/$3_${22}_trkqaul$1_nmin$4_nmax$5_cmin${14}_cmax${15}_aj${iajmin}_${iajmax}_tptmin$6_tptmax$7_aptmin$8_aptmax${9}.root unmerged/$3_${22}_trkqaul$1_nmin$4_nmax$5_cmin${14}_cmax${15}_aj${iajmin}_${iajmax}_tptmin$6_tptmax$7_aptmin$8_aptmax${9}_*.root" >> hadd.txt
# cat $now/runjettrk.condor

echo 
echo -e "\E[31mcondor_submit $now/runjettrk.condor"
tput sgr0 
echo "condor_submit $now/runjettrk.condor" >> condor.txt
echo
echo $now/runjettrk.sh 0 $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15} ${16} ${17} ${18} ${19} ${20} ${21} ${22} ${23}
echo -e "\E[33m"
echo
echo cd $now
echo
echo | awk -v i=0 -v trkq=$1 -v flist=$2 -v tag=$3 -v nmin=$4 -v nmax=$5 -v pttrigmin=$6 -v pttrigmax=$7 -v ptassmin=$8 -v ptassmax=${9} -v leadingjetptlow=${10} -v leadingjetpthigh=${11} -v subleadingjetptlow=${12} -v subleadingjetpthigh=${13} -v centmin=${14} -v centmax=${15} -v ajmin=${16} -v ajmax=${17} -v dotrkcorr=${18} -v mccommand=${19} -v jetetamin=${20} -v jetetamax=${21} -v whichjet=${22} -v bklist=${23} '{print "./runjettrk.exe "i" "trkq" "flist" "tag" "nmin" "nmax" "pttrigmin" "pttrigmax" "ptassmin" "ptassmax" "leadingjetptlow" "leadingjetpthigh" "subleadingjetptlow" "subleadingjetpthigh" "centmin" "centmax" "ajmin" "ajmax" "dotrkcorr" "mccommand" "jetetamin" "jetetamax" "whichjet" "bklist" "}' 
echo 
tput sgr0 
# condor_submit $now/runjettrk.condor

cd $now
# ln -s /net/hisrv0001/home/dav2105/git/jettrkcorr/final_hists_Vs3Calo
