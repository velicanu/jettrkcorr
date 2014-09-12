cp $1 $2.origin
root -l -q -b "skimTree.C+(\"$2.origin\",\"$2\")"
rm $2.origin
