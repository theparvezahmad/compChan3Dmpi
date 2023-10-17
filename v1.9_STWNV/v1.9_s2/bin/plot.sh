#!/bin/bash
set -e

if [[ -e ../output/moment1st.dat ]]; then

  source readParam.sh

	echo "==========================================="
	echo "Case : ${tag}"
	echo "-------------------------------------------"
  echo "Parameters file param.dat read"
  
  fc=ifort
  mpifc=mpiifort
  avg=avg_${avgVer}
  mod=comdata
  
  cd ../codes
  
  rm -f ${avg}.o ${mod}.o ${avg}.x ${mod}.mod
  echo "Object files & executables from previous builds removed"
  
  ${fc} ${flags} -c ${mod}.f90
  ${fc} ${flags} -c ${avg}.f90
  ${fc} ${flags} ${mod}.o ${avg}.o -o ${avg}.x
  
  echo "Code for averaging over time initiated"
  echo "-------------------------------------------"
  
  ./${avg}.x
  
  tim=$(date +%Hh_%d%b)
  pl=plots${tim}_${tag}
  
  if [[ ! -d ../plots ]]; then
  mkdir ../plots
  fi
  
  cd ../plots
  
  if [[ ! -d $pl ]]; then
  mkdir $pl
  fi
  
  cd $pl
  
  gnuplot ../../codes/cmdGNUplot.txt
  
  if [[ -e plotInfo.dat ]]; then
    rm plotInfo.dat
  fi
  
  mv ../../codes/plotInfo.dat plotInfo.dat
  
  echo "${descp[11]}:${value[11]}" >> plotInfo.dat
  echo "${descp[12]}:${value[12]}" >> plotInfo.dat
  echo "${descp[13]}:${value[13]}" >> plotInfo.dat
  echo "${descp[14]}:${value[14]}" >> plotInfo.dat
  echo "${descp[16]}:${value[16]}" >> plotInfo.dat
  echo "${descp[17]}:${value[17]}" >> plotInfo.dat
  echo "${descp[18]}:${value[18]}" >> plotInfo.dat
  echo "${descp[21]}:${value[21]}" >> plotInfo.dat
  echo "${descp[23]}:${value[23]}" >> plotInfo.dat
  echo "${descp[24]}:${value[24]}" >> plotInfo.dat
  
  echo "-------------------------------------------"
  echo "Plot information saved in plotInfo.dat"

	else
		echo "Found no spatially averaged field"
	fi  

  cd ../..
  
  if [[ -d ../allPlots ]]; then
    mkdir ../allPlots/${tag}
    cp -R plots/${pl} ../allPlots/${tag}/
    echo "A copy of current plots saved in folder allPlots"
  fi
  
  cd bin
  
