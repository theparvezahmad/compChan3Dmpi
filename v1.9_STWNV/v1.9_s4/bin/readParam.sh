#!/bin/bash
set -e

i=1;
while IFS=':' read -a var
do
	descp[$i]=${var[0]};
	value[$i]=${var[1]};
	i=$((i+1));
  #echo "${var[0]}" "${var[1]}";
done < "../param.dat"

IFS=' '

tag=${value[3]};
resOp=${value[4]};

xMPI=${value[6]};
yMPI=${value[7]};
zMPI=${value[8]};
list=${value[9]};

Re=${value[11]};
Mac=${value[12]};
dt=${value[13]};
istat=${value[14]};

PexVal=${value[16]};
PeyVal=${value[17]};
PezVal=${value[18]};

debugOp=${value[20]};
inputPath=${value[21]};
flags=${value[22]};
mainVer=${value[23]};
avgVer=${value[24]};
