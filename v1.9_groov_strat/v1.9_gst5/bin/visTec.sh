#!/bin/bash
set -e

cd ../codes
ifort combine.f90 -o combine.x
echo "============================================"
echo "Code to combine component files initiated"
echo "--------------------------------------------"
./combine.x

echo "Visualize instantaneous solution(1) or fluctuations(2) or Q(3) "
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

if [[ ${option} = 3 ]]; then
  ifort qInvariant.f90 -o qInvariant.x
  echo "Code to find Q initiated"
  echo "------------------------------------------"
  ./qInvariant.x
  cd ../output
  tec360 qInvariant.dat
fi

cd ../bin
