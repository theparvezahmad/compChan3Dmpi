# Author: Parvez Ahmad <parvez.cfd@gmail.com>| 25-02-17
.PHONY: cmd clobber clean
default: main

SHELL= /bin/bash
np=36
fc=ifort
mfc=mpiifort
dbg= -traceback -CB -warn all 
run= -O3#-CB#-O3# -ipo -xHost
mcode= v1.8_s2
fcode= avgV1.6
#-env I_MPI_PIN_PROCESSOR_LIST 16-63

main: cmd comData.o $(mcode).o
	$(mfc) $(run) comData.o $(mcode).o -o v1.8_s2.x
	nohup mpiexec -np 36 ./v1.8_s2.x &
	
avg: comData.o $(fcode).o
	$(fc) $(run) comData.o $(fcode).o -o $(fcode).x
	./$(fcode).x
	
comData.o: comData.f90
	$(fc) $(run) -c comData.f90

comData.mod: comData.f90
	$(fc) $(run) -c comData.f90

$(mcode).o: $(mcode).f90
	$(mfc) $(run) -c $(mcode).f90
	
$(fcode).o: $(fcode).f90
	$(fc) $(run) -c $(fcode).f90		

cmd:
	$(SHELL if [ ! -d "../output" ]; then mkdir ../output ; fi)
	$(SHELL if [ ! -d "../backup" ]; then mkdir ../backup ; fi)
  
clobber:
	$(SHELL if [ -d "../output" ]; then rm -r ../output ; fi)
	$(SHELL if [ -d "../plots" ]; then rm -r ../plots ; fi)
	rm -f *.x *.out *.o *.mod

clean:
	rm -f *.x *.o *.mod
