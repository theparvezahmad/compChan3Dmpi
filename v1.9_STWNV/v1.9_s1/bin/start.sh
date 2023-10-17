#!/bin/bash
set -e
#-----------------------------------------------------------------------
echo "============================================================="
echo "Compressible Supersonic Wall-bounded Turbulent Channel Flow "
echo "-------------------------------------------------------------"
echo "Developer : Parvez Ahmad (M.Tech 2015-2017)"
echo "Supervisor: Prof. Mirza Faisal S. Baig"
echo "============================================================="
echo ""

source readParam.sh

echo "==========================================="
echo "Case : ${tag}"
echo "-------------------------------------------"
echo "Parameters file read"
#-----------------------------------------------------------------------
source backup.sh
#-----------------------------------------------------------------------
if [[ ! -d "../output" ]]; then
	mkdir ../output
	echo "output folder created"
fi

if [[ ! -d "../debug" ]]; then
	mkdir ../debug
	echo "debug folder created"
fi
#-----------------------------------------------------------------------
cd ../codes

if [[ -e main_${mainVer}.f90 ]]; then
	mv main_${mainVer}.f90 main_${mainVer}_${tag}.f90
	echo "main_${mainVer}.f90 renamed to main_${mainVer}_${tag}.f90"
fi

fc=ifort
mpifc=mpiifort
main=main_${mainVer}_${tag}
mod=comdata

rm -f ${main}.o ${mod}.o ${main}.x ${mod}.mod combine.x
echo "Object files & executables from previous builds removed"

echo "Code is being compiled..."
${fc} ${flags} -c ${mod}.f90
${mpifc} ${flags} -c ${main}.f90
${mpifc} ${flags} ${mod}.o ${main}.o -o ${main}.x

${fc} ${flags} combine.f90 -o combine.x

if [[ ${resOp} = 1 && -e ../output/3D001.dat &&  ${inputPath} = ../input/3D.dat ]]; then
	echo "Code to combine component files initiated"
	./combine.x
	echo "Combined solution file is created in <Working Dir>/input"
fi
#-----------------------------------------------------------------------

nP=$[xMPI*yMPI*zMPI]

export I_MPI_PIN_PROCESSOR_LIST=$list

nohup mpiexec -np ${nP} ./${main}.x >> ../terminal.dat &
echo "Main code execution initiated"
#-----------------------------------------------------------------------
