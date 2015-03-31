if [ $# -lt 1 ]
then
  echo "Usage: ./balanced.sh <input-list>"
  exit 1
fi

nlines=`wc -l $1 | awk '{print $1}'`
for i in `seq 1 $nlines`
do
  subline=`head -n$i $1 | tail -n1`
  echo $subline | bash
  njobs=`condor_q $USER | grep jobs | awk '{print $1}'`
  while [ $njobs -gt 1000 ]
  do
    echo $njobs/1000
    sleep 60
    njobs=`condor_q $USER | grep jobs | awk '{print $1}'`
  done
  spaceleft=`quota | tail -n1 | awk '{print $2-$1}'`
  if [ $spaceleft -lt 100000 ]
  then
    echo "ERROR: Less than 100MB left in quota."
    echo "Last submitted job:"
    echo $subline
    exit 1
  fi
done

