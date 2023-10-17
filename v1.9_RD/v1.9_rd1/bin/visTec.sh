#!/bin/bash
set -e

cd ../codes
ifort combine.f90 -o combine.x
echo "============================================"
echo "Code to combine component files initiated"
echo "--------------------------------------------"
./combine.x

echo "Visualize instantaneous solution(1) or fluctuations(2) "
read option

if [[ ${option} = 1 ]]; then
  cd ../input
  tec360 3D.dat
fi

if [[ ${option} = 2 ]]; then
  ifort fluc.f90 -o fluc.x
  echo "Code to find fluctuations initiated"
  echo "------------------------------------------"
  ./fluc.x
  cd ../output
  tec360 fluc.dat
fi

cd ../bin
