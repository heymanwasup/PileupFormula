#! /bin/bash

source  /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh

# g++ `root-config --glibs --cflags` ./main.cxx -o main.exe

ctag=${1}
seed=$((${2}%100))
Nfills=${3}

outfile=/lustre/collider/chencheng/data/PU_formula/toyMC/raw_ctag${ctag}_seed${seed}.root

echo -e "running-${outfile}-${ctag}-${seed}-${Nfills}"

./main.exe ${outfile} ${ctag} ${seed} ${Nfills}

echo -e "finished-${outfile}-${ctag}-${seed}-${Nfills}"
