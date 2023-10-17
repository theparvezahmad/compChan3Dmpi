program chanCompr
!=======================================================================
!Developer : Parvez Ahmad (M.Tech 2015-2017)
!Supervisor: Prof. Mirza Faisal S. Baig
!-----------------------------------------------------------------------
!For best viewing, set Tab Width = 2
!For profuse comments, see the end
!=======================================================================
use module_comp
use module_const
use module_MPI
use module_read
use module_write

implicit none
integer:: nstep, iter
double precision:: time, tStart, tEnd

!========================Variable Declaration Ends======================
call readParam()
call calcMeshParam()
call calcCoeffs()
call setupMPI()
call writeParam()
call readInputFile()
call initSolV()

call share2P(solV)
call share2P(var)

DO WHILE (nstep < maxstep) ! Main time loop

	if(master .and. mod(nstep,10)==0) call cpu_time(tStart)
	nstep = nstep+1
	time  = time + dt
		
	do iter=1,maxiter ! RK2 Sub-steps
		call calcFlux()

		call calcSolV()
		call share2P(solV)

		call calcPrim()
		call share2P(var)

		call calcBC()
	enddo

	call checkSoln()
	call writeSoln()	
	call writeBulk()
	call calcTauWall()
	call avgSpatial()

	if(master .and. mod(nstep,10)==0) then
		call cpu_time(tEnd)
		write(6,'(I10,F10.5,A,F6.3,2(1X,A))') nstep,time," | CPU Time:",tEnd-tStart,"sec for 10 iters | ", dateTime()
		write(6,'(A)')"-----------------------------------------------------------------------------------"
	endif
		
END DO

call cleanupMPI()

!=======================================================================
if(master) write(*,*) "Compressible Channel Flow Solver exited at :",dateTime()
!=======================================================================
end program chanCompr
