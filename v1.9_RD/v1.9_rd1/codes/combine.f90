program combine
!Author : Parvez Ahmad <pahmed333@gmail.com>
!It combines separate files into one
!Last modified on 01/06/2017
!=======================================================================
implicit none
integer,parameter																::nDim=3
DOUBLE PRECISION																::time,nstep
DOUBLE PRECISION,dimension(:,:,:,:),allocatable	::var
integer,dimension(:,:,:),allocatable						::Pattern
DOUBLE PRECISION,dimension(:)    ,allocatable		::x1,x2,x3
integer,allocatable,dimension(:,:) 							::bb,ee
INTEGER																					::id,nproc,i,j,k,st,nX,nY,nZ
CHARACTER(LEN=50)																::filename
!====================Read Domain Info===================================
open(unit=15,file="../output/domInfo.txt")
read(15,105) nX,nY,nZ
read(15,106) nproc
allocate(bb(nDim,0:nproc-1),ee(nDim,0:nproc-1))
do i=0,nproc-1
	read(15,107) bb(1,i),ee(1,i),bb(2,i),ee(2,i),bb(3,i),ee(3,i)
enddo
close(15)
!-----------------------------------------------------------------------
105 format(3(i4))
106 format(i3)
107 format(6(i4,1x))

allocate(var(nX,nY,nZ,7))
allocate(x1(nX),x2(nY),x3(nZ))
allocate(Pattern(nX,nY,nZ))

!====================Read component files===============================
do id=0,nproc-1
	write(filename,'(a,i3.3,a)')'../output/3D',id+1,'.dat'
	open(unit=12,file=filename)
	
	read(12,100)time,nstep
	read(12,*)
	read(12,*)

	do k=bb(3,id),ee(3,id)
		read(12,*)
		do j=bb(2,id),ee(2,id)
		read(12,*)
			do i=bb(1,id),ee(1,id)
				read(12,'(3(1X,F10.6),7(1X,E22.15),1X,I1)',iostat=st) x1(i),x2(j),x3(k),var(i,j,k,1),var(i,j,k,2),var(i,j,k,3),var(i,j,k,4)	&
																														 ,var(i,j,k,5),var(i,j,k,6),var(i,j,k,7),pattern(i,j,k)
																														 
					if(st/=0)	then
						write(6,'(A,5(1X,I4))') "Error in reading component input file at",id,i,j,k,st
						stop
					endif	
																																		 
			enddo
		enddo
	enddo
	
	close(12)
enddo

write(6,*)"Component files read successfully"
100	format(8X,F20.10,I9,1X)
101	format(3(1X,F10.6),7(1X,E22.15))
!-----------------------------------------------------------------------
!=====================Write combined file===============================
open(unit=12,file='../input/3D.dat')
write(12,103) 'TITLE ="',time,nstep,'"'
write(12,*)'Variables = "x","y","z","Rho","u","v","w","e","T","P","Pattern"'
write(12,104)'ZONE k=',nZ,',j=',nY,',i=',nX,',DATAPACKING=POINT'
do k=1,nZ
	write(12,*)
	do j=1,nY
		write(12,*)
		do i=1,nX
			write(12,'(3(1X,F10.6),7(1X,E22.15),1X,I1)')x1(i),x2(j),x3(k),var(i,j,k,1),var(i,j,k,2),var(i,j,k,3),var(i,j,k,4)	&
																								 ,var(i,j,k,5),var(i,j,k,6),var(i,j,k,7),pattern(i,j,k)
		enddo
	enddo
enddo
close(12)
!-----------------------------------------------------------------------
103	format(A,F20.10,I9,A)
104	format(A,I3,A,I3,A,I3,A)

write(6,*)"Component files combined successfully"
end program combine
