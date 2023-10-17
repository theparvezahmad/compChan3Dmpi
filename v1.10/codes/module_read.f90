module module_read
use module_comp, only:var
use module_const, only:i,j,k,st
implicit none

integer,protected::resOp,topo(3), istat, debugOp
double precision,protected:: Re,invRe,Mac,dt,PexVal,PeyVal,PezVal,ggM

character(len=:),allocatable,protected::descp,inputPath

contains
subroutine readParam()
implicit none
character(len=10),private::descp_
character(len=100),private::inputPath_

open(unit=10,file="../param.dat",status="old")

read(10,*)
read(10,*)
read(10,'(48X,A)') descp_
read(10,'(48X,i1)') resOp
read(10,*)
read(10,'(48X,i2)') topo(1)
read(10,'(48X,i2)') topo(2)
read(10,'(48X,i2)') topo(3)
read(10,*)
read(10,*)
read(10,'(48X,F4.0)') Re
read(10,'(48X,f3.1)') Mac
read(10,'(48X,f6.4)') dt
read(10,'(48X,i2)') istat
read(10,*)
read(10,'(48X,f3.1)') PexVal
read(10,'(48X,f3.1)') PeyVal
read(10,'(48X,f3.1)') PezVal
read(10,*)
read(10,'(48X,i1)') debugOp
read(10,'(48X,A)') inputPath_

close(10)

invRe = 1.0d0/Re
ggM = gama*(gama-1)*Mac**two

allocate(character(len=len(trim(descp_))) :: descp)
allocate(character(len=len(trim(inputPath_ ))) :: inputPath )

descp=trim(descp_)
inputPath =trim(inputPath_)

end subroutine readParam

subroutine readInputFile()
implicit none
select case(resOp)
	case(1)
		!======================Read Combined Input File=====================
		open(unit=10,file=inputPath)

		read (10,100)time,nstep
		read (10,*)
		read (10,*)

		100	format(8X,F20.10,I9,1X)

		do k=1,nP(3)
			read(10,*)
			do j=1,nP(2)
				read(10,*)
				do i=1,nP(1)
					read(10,'(3(1X,F10.6),7(1X,E22.15))',iostat=st) x1(i),x2(j),x3(k),var(i,j,k,1),var(i,j,k,2),var(i,j,k,3),var(i,j,k,4)	&
																															 ,var(i,j,k,5),var(i,j,k,6),var(i,j,k,7)
																										 
					if(st/=0)	then
						write(6,'(A,4(1X,I4))') "Error in reading combined input file at",i,j,k,st
						stop
					endif
						
				end do
			end do
		end do

		close(10)
		!101	format(3(1X,F10.6),7(1X,E22.15))
	case(2)
		!======================Read Component Input File====================
		WRITE(filename,'(a,i3.3,a)')"../output/3D",id+1,".dat"
		OPEN (UNIT=10,FILE=filename)

		read(10,100) time,nstep
		read(10,*)
		read(10,*)

		do k=bs(3),es(3)
			read(10,*)
			do j=bs(2),es(2)
				read(10,*)
				do i=bs(1),es(1)
					read(10,'(3(1X,F10.6),7(1X,E22.15))',iostat=st) x1(i),x2(j),x3(k),var(i,j,k,1),var(i,j,k,2),var(i,j,k,3),var(i,j,k,4)	&
																															 ,var(i,j,k,5),var(i,j,k,6),var(i,j,k,7)
					if(st/=0)	then
						write(6,'(A,5(1X,I4))') "Error in reading component input file at",id,i,j,k,st
						write(6,'(A)') "Change in MPI process topology since last run could possibly be a reason"
						stop
					endif
																										 
				end do
			end do
		enddo
		close(10)
		!-----------------------------------------------------------------------
end select

if(master) write(6,'(A)') "Input file successfully read"
end subroutine readInputFile

end module_read
