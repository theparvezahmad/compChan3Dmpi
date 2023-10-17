module module_write
use mpi
use module_read
implicit none

double precision,protected			::tmpA1(5)
integer,dimension(5)  :: plane =[5,10,70,100,125]

contains
subroutine writeBulk()
implicit none
if (mod(nstep,10).eq.0) then
		!------------------------Write Kinetic Energy-----------------------
		tmpR1=0.0d0
		tmpR2	=0.0d0
		do k=bs(3),es(3)
			do j=bs(2),es(2)
				do i=bs(1),es(1)
					tmpR1 = tmpR1 + sum(var(i,j,k,2:4)**two)  
				enddo
			enddo
		enddo

		CALL MPI_Reduce(tmpR1,tmpR2,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
		
		if(master) then
			tmpR2 = tmpR2/product(nP)
			open(10,file="../output/KE.dat",access="append")
			write(10,108) time,tmpR2
			close(10)
		endif
		108	format(2(1X,F15.7))
		!------------------------Write Bulk Density-------------------------
		Ayz = x3(nP(3))*x2(nP(2))
		
		tmpA1	=0.0d0
		do i=1,5
			tmpR1	=0.0d0
			
			if(plane(i) > bn(1) .and. plane(i) < en(1)) then
				do k=bs(3),en(3)
					do j=bs(2),en(2)
						dAyz = dy*(x3(k+1)-x3(k))
						tmpR1 =tmpR1 + quar*(	var(plane(i),j  ,k  ,1) +	&
																	var(plane(i),j+1,k  ,1) +	&
																	var(plane(i),j  ,k+1,1) +	&
																	var(plane(i),j+1,k+1,1)		&
																)*dAyz
					enddo
				enddo
			endif
			
			CALL MPI_Reduce(tmpR1,tmpA1(i),1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
		enddo  

		if(master) then
			tmpA1=tmpA1/Ayz
			open (unit=10,file="../output/rhobulk.dat",access="append")
			write(10,109)time,tmpA1(1),tmpA1(2),tmpA1(3),tmpA1(4),tmpA1(5)
			close(10)
		endif
		109	format(6(1X,F15.7))
		!------------------------Write Bulk Velocity------------------------ 
		tmpA1	=0.0d0
		do i=1,5
			tmpR1	=0.0d0
			
			if(plane(i) > bn(1) .and. plane(i) < en(1)) then
				do k=bs(3),en(3)
					do j=bs(2),en(2)
						dAyz = dy*(x3(k+1)-x3(k))
						tmpR1 =tmpR1 + quar*(	var(plane(i),j  ,k  ,2) +	&
																	var(plane(i),j+1,k  ,2) +	&
																	var(plane(i),j  ,k+1,2) +	&
																	var(plane(i),j+1,k+1,2)		&
																)*dAyz
					enddo  
				enddo  
			endif
			
			CALL MPI_Reduce(tmpR1,tmpA1(i),1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
		enddo  
		
		if(master) then
			tmpA1=tmpA1/Ayz
			open (unit=10,file="../output/Ubulk.dat",access="append")
			write(10,109)time,tmpA1(1),tmpA1(2),tmpA1(3),tmpA1(4),tmpA1(5)
			close(10)
		endif
		
		!--------------------------Write Massflow--------------------------- 
		tmpA1	=0.0d0
		do i=1,5
			tmpR1	=0.0d0
			
			if(plane(i) > bn(1) .and. plane(i) < en(1)) then
				do k=bs(3),en(3)
					do j=bs(2),en(2)
						dAyz = dy*(x3(k+1)-x3(k))
						tmpR1 =tmpR1 + quar*(	var(plane(i),j  ,k  ,1)*var(plane(i),j  ,k  ,2) +	&
																	var(plane(i),j+1,k  ,1)*var(plane(i),j+1,k  ,2) + &
																	var(plane(i),j  ,k+1,1)*var(plane(i),j  ,k+1,2) +	&
																	var(plane(i),j+1,k+1,1)*var(plane(i),j+1,k+1,2)		&
																)*dAyz
					enddo  
				enddo  
			endif
			
			CALL MPI_Reduce(tmpR1,tmpA1(i),1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
		enddo  

		if(master) then
			tmpA1=tmpA1/Ayz
			open (unit=10,file="../output/massflow.dat",access="append")
			write(10,109)time,tmpA1(1),tmpA1(2),tmpA1(3),tmpA1(4),tmpA1(5)
			close(10)
		endif
		
	endif	
end subroutine writeBulk

subroutine writeSoln()
implicit none
character(len=30),private			:: filename

if(mod(nstep,100).eq.0) then
	
	if(mod(nstep/100,2).eq.0) then
		WRITE(filename,'(a,i3.3,a)')"../output/3D",id+1,".dat"
	else	
		WRITE(filename,'(a,i3.3,a)')"../output/3D",id+1,"_.dat"
	endif
		
	OPEN (UNIT=10,FILE=filename)

	write(10,'(A,F20.10,I9,A)') 'TITLE ="',time,nstep,'"'
	write(10,*) 'Variables ="x","y","z","Rho","U","V","W","E","T","P","Pattern"'
	write(10,'(A,I3,A,I3,A,I3,A)') 'ZONE k=',es(3)-bs(3)+1,',j=',es(2)-bs(2)+1,',i=',es(1)-bs(1)+1,',DATAPACKING="POINT"'

	do k=bs(3),es(3)
		write(10,*)
		do j=bs(2),es(2)
			write(10,*)
			do i=bs(1),es(1)
				write(10,'(3(1X,F10.6),7(1X,E22.15),1X,I1)') x1(i),x2(j),x3(k),var(i,j,k,1),var(i,j,k,2),var(i,j,k,3),var(i,j,k,4)	&
																																		,var(i,j,k,5),var(i,j,k,6),var(i,j,k,7),pattern(i,j,k)
			end do
		end do
	enddo
	close(10)

	!-------------------------------------------------------------------
	call mpi_gather(bs(1),3,mpi_integer,bb(1,0),3,mpi_integer,0,mpi_comm_world,ierr)
	call mpi_gather(es(1),3,mpi_integer,ee(1,0),3,mpi_integer,0,mpi_comm_world,ierr)
	!-----------------------Write Domain Info---------------------------
	if(master) then
		open(unit=10,file="../output/domInfo.txt")
		write(10,'(3(i4))') nP(1),nP(2),nP(3)
		write(10,'(i3)') nproc
		do i=0,nproc-1
			write(10,'(6(i4,1x))') bb(1,i),ee(1,i),bb(2,i),ee(2,i),bb(3,i),ee(3,i)
		enddo
		close(10)
	endif
	
endif

end subroutine writeSoln

subroutine writeParam()
implicit none
if(master) then
	write(6,'(A)')"==================================================================================="
	write(6,'(2A)')"Compressible Channel Flow Solver started at :",dateTime()
	write(6,'(A)')"==================================================================================="
	write(6,'(A)')
	write(6,'(A)')"==================================================================================="
	write(6,'(A)')"Parameters read from param.dat"
	write(6,'(A)')"-----------------------------------------------------------------------------------"
	write(6,'(2A)')"Case Description        : ",descp
	write(6,'(2A)')"Path to input file      : ",inputPath
	write(6,'(A,3(I3,1X))')"Process Topology        : ",topo(1),topo(2),topo(3)
	write(6,'(A,F10.4)'   )"Reynolds Number         : ",Re
	write(6,'(A,F10.2)'   )"Mach Number             : ",Mac
	write(6,'(A,F10.8)'   )"Timestep                : ",dt
	write(6,'(A,I10)'     )"Freq of stat calculation: ",istat
	write(6,'(A,F10.4)'   )"Peclet in x-direction   : ",PexVal
	write(6,'(A,F10.4)'   )"Peclet in y-direction   : ",PeyVal
	write(6,'(A,F10.4)'   )"Peclet in z-direction   : ",PezVal
	write(6,'(A,I10)'     )"Restart Option  [1/2]   : ",resOp
	write(6,'(A,I10)'     )"Debug Option    [0/1]   : ",debugOp
	write(6,'(A)')"==================================================================================="
endif
end subroutine writeParam

end module module_write
