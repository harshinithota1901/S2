#!/bin/bash

if [ $# -ne 5 ]
then
  echo Invocation: $0 InputsDir OutputsDir kbConstraints RandomRepeats:1,2,etc  PopRepeats:500,1000,etc 
  exit
fi

inDir=$1
if [ ! -d ${inDir} ]
then
  echo ${inDir} input dir does not exist
  exit
fi

outDir=$2
if [ ! -d ${outDir} ]
then
  mkdir ${outDir}
else
  echo ${outDir} must be nonexistent
  exit
fi

kbInput=$3
if [ ! -e $kbInput ]
then
  echo $kbInput could not be open
  exit
fi

popRepeats=$4
randomRepeats=$5

for i in `ls ${inDir}/*`
do
  outBase=${outDir}/${i##*/}
  r=1
  while [ $r -le $randomRepeats ] 
  do
    p=1
    popSize=500
    while [ $p -le $popRepeats ]
    do
      out=${outBase}.P${popSize}.R${r}
      echo RUNNING acgp1.2 -f $i -p output.basename=${out} -p random_seed=$r -p pop_size=$popSize <${kbInput}
      acgp1.2 -f $i -p output.basename=${out} -p random_seed=$r -p pop_size=${popSize} <${kbInput}
      rm ${out}.bst ${out}.cnt ${out}.gen ${out}.his ${out}.prg ${out}.sys
      let p=${p}+1
      let popSize=${popSize}+500
    done
    let r=${r}+1
  done
  echo INPUT $i FINISHED
done
echo ALL FINISHED -----------------------------------------------------------
