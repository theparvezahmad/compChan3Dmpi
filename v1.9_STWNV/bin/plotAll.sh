#!/bin/bash
set -e

if [[ -d ../allPlots ]]; then
	rm -rf ../allPlots/*
fi

cd ../v1.9_s1/bin
source plot.sh

cd ../../v1.9_s2/bin
source plot.sh

cd ../../v1.9_s3/bin
source plot.sh

cd ../../v1.9_s4/bin
source plot.sh

cd ../../v1.9_s5/bin
source plot.sh

cd ../../v1.9_s6/bin
source plot.sh

echo "========================================="
echo "Plots for all Cases s1 to s5 generated"
echo "========================================="
