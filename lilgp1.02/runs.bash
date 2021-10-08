#!/bin/bash

if [ $# -ne 4 ]
then
  echo "Invocation: $0 whichGP whichParam numRuns(random seed = 1..numRuns) outBaseName"
  exit
fi

whichGP=$1
whichParam=$2
numRuns=$3
outBaseName=$4

for (( i = 1; i <= ${numRuns}; i++ )) 
do
  out=${outBaseName}.r$i
  echo RUNNING ${whichGP} -f ${whichParam} -p output.basename=${out} -p random_seed=$i 
  "$whichGP" -f ${whichParam} -p output.basename=${out} -p random_seed=$i
  rm ${out}.bst ${out}.gen ${out}.his ${out}.prg ${out}.sys
done
echo ALL FINISHED -----------------------------------------------------------

