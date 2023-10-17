#!/bin/bash
set -e

if [[ -e ../output/3D001.dat ]]; then

  source readParam.sh
  
  if [[ ! -d "../backup" ]]; then
  	mkdir ../backup
  fi
  
  read -n 17 -a tim < "../output/3D001.dat"
  
  var=$(date +%Hh%d%b)
  folder=date${var}_ndt${tim[2]}_${tag}
  
  if [[ ! -d ../backup/${folder} ]]; then
  	mkdir ../backup/${folder}
  fi

  for file in ../output/3D*; do cp "$file" ../backup/${folder} ;done
  echo "Solution backup created in <Working Dir>/backup"
  
else
  echo "Found no solution file to backup"
fi
