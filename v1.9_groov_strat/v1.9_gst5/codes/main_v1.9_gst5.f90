program chanCompr
!=======================================================================
use mpi
use comData
implicit none
!=======================Variable Declaration Starts=====================
character(len=30)			:: filename
logical								:: isBlown=.false.
integer								:: i,j,k,m,p,nstep,iter,st,nXY,pattern(nP(1),nP(2),nP(3))
integer,dimension(5)  :: plane =[5,10,70,100,125]
double precision 			:: tStart, tEnd, tauxzL,tauxzU
double precision 			:: time, Pe, flux,dAyz,Ayz
double precision 			:: f1d2c0, f1d2c1, f1d2c2, b1d2c0, b1d2c1, b1d2c2
double precision 			:: f2d1c0, f2d1c1, f2d1c2, b2d1c0, b2d1c1, b2d1c2
double precision 			:: dz1,dzN,alpha

!--------------------------Temporary Variables--------------------------
integer								::tmpI1,tmpI2,tmpI3,tmpI4,tmpI5
double precision			::tmpR1,tmpR2,tmpR3,tmpR4,tmpA1(5),loc(nP(1),nP(2),7)
double precision			::Ut,Vt,Wt,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
double precision			::dTdx,dTdy,dTdz,divg,muT,ktT,PressT,dissF,dissG,dissH
double precision			::dudxdz,dvdydz,d2wdz2,drhow2dz,lhs
integer,allocatable,dimension(:,:) 	::tArr,bb,ee

!------------------------MPI Variables----------------------------------
integer,dimension(3)		:: bs,es,bn,en,siz,des,src,myDim,myCoord
integer		:: liney,yz1p, xz1p, xy1p, yz1p5v, xz1p5v, xy1p5v, yz1p7v, xz1p7v, xy1p7v
integer		:: lineyN, yzN1p, xzN1p, xyN1p, yzN2p, xzN2p, xyN2p, yz2p5v, xz2p5v, xy2p5v
integer		:: comm3d,nproc, id, ierr, sizDP, STATUS(mpi_status_size)
logical 	:: myperiod(3), bottom, top, east, west, north, south, master

!-----------------------------------------------------------------------
double precision,dimension(nP(3),7) :: stat1Spc	! In order of rho,u,v,w,e,t,p 
double precision,dimension(nP(3),7) :: stat2Spc ! In order of rho,u,v,w,e,t,p
double precision,dimension(nP(3),9) :: corrSpc 	! In order of rhow,rhot,uv,uw,vw,vt,tw,pw,pt

double precision,dimension(0:nP(1)+1,0:nP(2)+1,0:nP(3)+1,5,2)	:: solV
double precision,dimension(nP(1),nP(2),nP(3),7)								:: var ! In order of rho,U,V,W,E,T,P
double precision,dimension(nP(1),nP(2),nP(3),5)								:: fD,gD,hD
double precision,dimension(nP(1),nP(2),nP(3),5)								:: srcV
double precision,dimension(nP(1),nP(2),nP(3),5,2)							:: F,G,H,R
                                    
double precision,dimension(nP(1),nP(2),nP(3))									:: mu,kt

!-------------------------Debug Variables-------------------------------
double precision,dimension(nP(3))			:: pl

!========================Variable Declaration Ends======================
call readParam()
call setupMesh()
!-----------------------------------------------------------------------
dz1		= x3(2)-x3(1)
alpha	= (x3(3)-x3(2))/dz1

tmpR1	= dz1*alpha*(alpha + one)
f1d2c0	= (one - (alpha + one)**two)/tmpR1
f1d2c1	= ((alpha + one)**two)/tmpR1
f1d2c2	= -one / tmpR1

tmpR1	= haf*dz1**two*alpha*(alpha + one)
f2d1c0	= alpha/tmpR1
f2d1c1	= -(alpha + one)/tmpR1
f2d1c2	= one / tmpR1
!-----------------------------------------------------------------------
dzN		= x3(nP(3))-x3(nP(3)-1)
alpha	= (x3(nP(3)-1)-x3(nP(3)-2))/dzN

tmpR1	= dzN*alpha*(alpha + one)
b1d2c0	= ((alpha + one)**two - one)/tmpR1
b1d2c1	= -((alpha + one)**two)/tmpR1
b1d2c2	= one / tmpR1

tmpR1	= haf*dzN**two*alpha*(alpha + one)
b2d1c0	= alpha/tmpR1
b2d1c1	= -(alpha + one)/tmpR1
b2d1c2	= one / tmpR1
!-----------------------------------------------------------------------

!==============================MPI======================================
CALL mpi_init(ierr)
CALL mpi_comm_rank(mpi_comm_world,id,ierr)
CALL mpi_comm_size(mpi_comm_world,nproc,ierr)

call mpi_cart_create(mpi_comm_world,3,topo,period,.false.,comm3d,ierr)
call mpi_cart_get(comm3d,3,myDim,myperiod,mycoord,ierr)

call mpi_cart_shift(comm3d,0,1,src(1),des(1),ierr)
call mpi_cart_shift(comm3d,1,1,src(2),des(2),ierr)
call mpi_cart_shift(comm3d,2,1,src(3),des(3),ierr)

allocate(tArr(0:maxval(myDim)-1,3))
allocate(bb(3,0:nproc-1),ee(3,0:nproc-1))

master=(id==0)
east	=(mycoord(1)==myDim(1)-1)
west	=(mycoord(1)==0)
north	=(mycoord(2)==myDim(2)-1)
south	=(mycoord(2)==0)
top		=(mycoord(3)==myDim(3)-1)
bottom=(mycoord(3)==0)

do k=1,3

	tmpI1=nP(k)/myDim(k)
	tmpI2=myDim(k)-mod(nP(k),myDim(k))

	tArr(0,1)=1
	do i=0,myDim(k)-1
		if(i==tmpI2) tmpI1=tmpI1+1
		tArr(i,2)=tArr(i,1)+tmpI1-1
		tArr(i,3)=tArr(i,2)-tArr(i,1)+1
		if(i==myDim(k)-1) exit
		tArr(i+1,1)=tArr(i,2)+1
	enddo

	do i=0,myDim(k)-1
		if(i==mycoord(k)) then
			bs(k)=tArr(i,1)
			es(k)=tArr(i,2)
			siz(k)=tArr(i,3)
		endif
	enddo

	bn(k)=bs(k)
	en(k)=es(k)
	
	if((west .and. k==1) .or. (south 	.and. k==2) .or. (bottom .and. k==3)) then
		bn(k)=bs(k)+1
		siz(k)=siz(k)-1
	endif
		
	if((east .and. k==1) .or. (north 	.and. k==2) .or. (top 	 .and. k==3))	then
		en(k)=es(k)-1
		siz(k)=siz(k)-1
	endif
		
enddo
!-----------------------------------------------------------------------
CALL mpi_type_extent (mpi_double_precision,sizDP,ierr)
CALL mpi_type_vector (siz(2), 1, nP(1), mpi_double_precision,liney ,ierr)

CALL mpi_type_hvector(siz(3), 1			, nP(1)*nP(2)*sizDP	,liney								,yz1p,ierr)
CALL mpi_type_vector (siz(3), siz(1), nP(1)*nP(2)				, mpi_double_precision,xz1p,ierr)
CALL mpi_type_vector (siz(2), siz(1), nP(1)							, mpi_double_precision,xy1p,ierr)

CALL mpi_type_commit(yz1p,ierr)
CALL mpi_type_commit(xz1p,ierr)
CALL mpi_type_commit(xy1p,ierr)

CALL mpi_type_hvector(5, 1 , product(nP)*sizDP , yz1p , yz1p5v ,ierr)
CALL mpi_type_hvector(5, 1 , product(nP)*sizDP , xz1p , xz1p5v ,ierr)
CALL mpi_type_hvector(5, 1 , product(nP)*sizDP , xy1p , xy1p5v ,ierr)

CALL mpi_type_commit(yz1p5v,ierr)
CALL mpi_type_commit(xz1p5v,ierr)
CALL mpi_type_commit(xy1p5v,ierr)

CALL mpi_type_hvector(7, 1 , product(nP)*sizDP , yz1p , yz1p7v ,ierr)
CALL mpi_type_hvector(7, 1 , product(nP)*sizDP , xz1p , xz1p7v ,ierr)
CALL mpi_type_hvector(7, 1 , product(nP)*sizDP , xy1p , xy1p7v ,ierr)

CALL mpi_type_commit(yz1p7v,ierr)
CALL mpi_type_commit(xz1p7v,ierr)
CALL mpi_type_commit(xy1p7v,ierr)
!-----------------------------------------------------------------------
CALL mpi_type_vector (siz(2), 1, nP(1)+2, mpi_double_precision,lineyN ,ierr)

CALL mpi_type_hvector(siz(3), 1, (nP(1)+2)*(nP(2)+2)*sizDP,lineyN,yzN1p,ierr)
CALL mpi_type_vector (siz(3), siz(1), (nP(1)+2)*(nP(2)+2), mpi_double_precision,xzN1p ,ierr)
CALL mpi_type_vector (siz(2), siz(1), nP(1)+2, mpi_double_precision,xyN1p ,ierr)

CALL mpi_type_hvector(2, 1, sizDP      								, yzN1p,yzN2p,ierr)
CALL mpi_type_hvector(2, 1, (nP(1)+2)*sizDP   				, xzN1p,xzN2p,ierr)
CALL mpi_type_hvector(2, 1, (nP(1)+2)*(nP(2)+2)*sizDP	, xyN1p,xyN2p,ierr)

CALL mpi_type_hvector(5, 1 , product(nP+2)*sizDP , yzN2p , yz2p5v ,ierr)
CALL mpi_type_hvector(5, 1 , product(nP+2)*sizDP , xzN2p , xz2p5v ,ierr)
CALL mpi_type_hvector(5, 1 , product(nP+2)*sizDP , xyN2p , xy2p5v ,ierr)

CALL mpi_type_commit(yz2p5v,ierr)
CALL mpi_type_commit(xz2p5v,ierr)
CALL mpi_type_commit(xy2p5v,ierr)

CALL mpi_barrier(mpi_comm_world,ierr)
!-----------------------------MPI---------------------------------------

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
					read(10,101,iostat=st) x1(i),x2(j),x3(k),var(i,j,k,1),var(i,j,k,2),var(i,j,k,3),var(i,j,k,4)	&
																															 ,var(i,j,k,5),var(i,j,k,6),var(i,j,k,7)
																										 
					if(st/=0)	then
						write(6,'(A,4(1X,I4))') "Error in reading combined input file at",i,j,k,st
						stop
					endif
						
				end do
			end do
		end do

		close(10)
		101	format(3(1X,F10.6),7(1X,E22.15))
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
					read(10,101,iostat=st) x1(i),x2(j),x3(k),var(i,j,k,1),var(i,j,k,2),var(i,j,k,3),var(i,j,k,4)	&
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

!--------------------Initial Solution Vector----------------------------
do k=bs(3),es(3)
	do j=bs(2),es(2)
		do i=bs(1),es(1)
			solV(i,j,k,1,1)		=	var(i,j,k,1)
			solV(i,j,k,2:5,1)	=	var(i,j,k,1) * var(i,j,k,2:5)
			
			mu(i,j,k)					= var(i,j,k,6)**0.7d0 
			kt(i,j,k)					= var(i,j,k,6)**0.7d0
		enddo
	enddo
enddo

!=====================================================================	
call mpi_sendrecv(solV(en(1)-1,bn(2),bn(3),1,1),1,yz2p5v,des(1),50,solV(bn(1)-2,bn(2),bn(3),1,1),1,yz2p5v,src(1),50,mpi_comm_world,STATUS,ierr)
call mpi_sendrecv(solV(bn(1)  ,bn(2),bn(3),1,1),1,yz2p5v,src(1),50,solV(en(1)+1,bn(2),bn(3),1,1),1,yz2p5v,des(1),50,mpi_comm_world,STATUS,ierr)

call mpi_sendrecv(solV(bn(1),en(2)-1,bn(3),1,1),1,xz2p5v,des(2),50,solV(bn(1),bn(2)-2,bn(3),1,1),1,xz2p5v,src(2),50,mpi_comm_world,STATUS,ierr)
call mpi_sendrecv(solV(bn(1),bn(2)  ,bn(3),1,1),1,xz2p5v,src(2),50,solV(bn(1),en(2)+1,bn(3),1,1),1,xz2p5v,des(2),50,mpi_comm_world,STATUS,ierr)

call mpi_sendrecv(solV(bn(1),bn(2),en(3)-1,1,1),1,xy2p5v,des(3),50,solV(bn(1),bn(2),bn(3)-2,1,1),1,xy2p5v,src(3),50,mpi_comm_world,STATUS,ierr)
call mpi_sendrecv(solV(bn(1),bn(2),bn(3)  ,1,1),1,xy2p5v,src(3),50,solV(bn(1),bn(2),en(3)+1,1,1),1,xy2p5v,des(3),50,mpi_comm_world,STATUS,ierr)
!---------------------------------------------------------------------
call mpi_sendrecv(var(en(1),bn(2),bn(3),1),1,yz1p7v,des(1),50,var(bn(1)-1,bn(2),bn(3),1),1,yz1p7v,src(1),50,mpi_comm_world,STATUS,ierr)
call mpi_sendrecv(var(bn(1),bn(2),bn(3),1),1,yz1p7v,src(1),50,var(en(1)+1,bn(2),bn(3),1),1,yz1p7v,des(1),50,mpi_comm_world,STATUS,ierr)

call mpi_sendrecv(var(bn(1),en(2),bn(3),1),1,xz1p7v,des(2),50,var(bn(1),bn(2)-1,bn(3),1),1,xz1p7v,src(2),50,mpi_comm_world,STATUS,ierr)
call mpi_sendrecv(var(bn(1),bn(2),bn(3),1),1,xz1p7v,src(2),50,var(bn(1),en(2)+1,bn(3),1),1,xz1p7v,des(2),50,mpi_comm_world,STATUS,ierr)

call mpi_sendrecv(var(bn(1),bn(2),en(3),1),1,xy1p7v,des(3),50,var(bn(1),bn(2),bn(3)-1,1),1,xy1p7v,src(3),50,mpi_comm_world,STATUS,ierr)
call mpi_sendrecv(var(bn(1),bn(2),bn(3),1),1,xy1p7v,src(3),50,var(bn(1),bn(2),en(3)+1,1),1,xy1p7v,des(3),50,mpi_comm_world,STATUS,ierr)

!========================Main time loop (starts)========================
DO WHILE (nstep < maxstep)

	if(master .and. mod(nstep,10)==0) call cpu_time(tStart)
	nstep = nstep+1
	time  = time + dt
	
	!========================RK2 Sub-steps (starts)=======================
	do iter=1,maxiter

		do k=bn(3),en(3)
			do j=bn(2),en(2)
				do i=bn(1),en(1)
					
					Pe = Re*abs(solV(i,j,k,2,1))*dx/mu(i,j,k) 

					if(Pe.gt.PexVal) then
						tmpR1  = ( var(i,j,k,2)+var(i+1,j,k,2) )*haf
						tmpR2  = ( var(i,j,k,2)+var(i-1,j,k,2) )*haf

						if(tmpR1 .ge. 0.0d0 .and. tmpR2 .ge. 0.0d0) then
							F(i,j,k,1:5,1) =(	tmpR1*( 6.0d0*solV (i  ,j,k,1:5,1) + 3.0d0*solV (i+1,j,k,1:5,1) - 1.0d0*solV (i-1,j,k,1:5,1) ) -	& 
																tmpR2*( 6.0d0*solV (i-1,j,k,1:5,1) + 3.0d0*solV (i  ,j,k,1:5,1) - 1.0d0*solV (i-2,j,k,1:5,1) )		&
															)*invdx/8.0d0 
						endif

						if(tmpR1 .lt. 0.0d0 .and. tmpR2 .lt. 0.0d0) then
							F(i,j,k,1:5,1) =(	tmpR1*( 6.0d0*solV (i+1,j,k,1:5,1) + 3.0d0*solV (i  ,j,k,1:5,1) - 1.0d0*solV (i+2,j,k,1:5,1) ) -	&
																tmpR2*( 6.0d0*solV (i  ,j,k,1:5,1) + 3.0d0*solV (i-1,j,k,1:5,1) - 1.0d0*solV (i+1,j,k,1:5,1) ) 		&
															)*invdx/8.0d0
						endif

						if(tmpR1 .lt. 0.0d0 .and. tmpR2 .ge. 0.0d0) then
							F(i,j,k,1:5,1) =(	tmpR1*( 6.0d0*solV (i+1,j,k,1:5,1) + 3.0d0*solV (i,j,k,1:5,1) - 1.0d0*solV (i+2,j,k,1:5,1) ) -	&
																tmpR2*( 6.0d0*solV (i-1,j,k,1:5,1) + 3.0d0*solV (i,j,k,1:5,1) - 1.0d0*solV (i-2,j,k,1:5,1) ) 		&
															)*invdx/8.0d0 
						endif

						if(tmpR1 .ge. 0.0d0 .and. tmpR2 .lt. 0.0d0) then
							F(i,j,k,1:5,1) =(	tmpR1*( 6.0d0*solV (i,j,k,1:5,1) + 3.0d0*solV (i+1,j,k,1:5,1) - 1.0d0*solV (i-1,j,k,1:5,1) ) -	&
																tmpR2*( 6.0d0*solV (i,j,k,1:5,1) + 3.0d0*solV (i-1,j,k,1:5,1) - 1.0d0*solV (i+1,j,k,1:5,1) ) 		&
															)*invdx/8.0d0 
						endif

					else   ! Fourth order Central Differencing																					
						F(i,j,k,1:5,1) = 	( -				var(i+2,j,k,2)*solV (i+2,j,k,1:5,1)	&
																+	8.0d0*var(i+1,j,k,2)*solV (i+1,j,k,1:5,1)	&
																-	8.0d0*var(i-1,j,k,2)*solV (i-1,j,k,1:5,1)	&
																+				var(i-2,j,k,2)*solV (i-2,j,k,1:5,1)	&
															)*invdx*one12th																					
					endif
					!-------------------------------------------------------------
					Pe = Re*abs(solV(i,j,k,3,1))*dy/mu(i,j,k)
					
					if(Pe .gt. PeyVal) then
						tmpR1  = ( var(i,j,k,3)+var(i,j+1,k,3) )*haf
						tmpR2  = ( var(i,j,k,3)+var(i,j-1,k,3) )*haf

						if(tmpR1.ge.0.0d0.and.tmpR2.ge.0.0d0) then
							G(i,j,k,1:5,1) =(	tmpR1*( 6.0d0*solV (i,j  ,k,1:5,1) + 3.0d0*solV (i,j+1,k,1:5,1) - 1.0d0*solV (i,j-1,k,1:5,1) ) -	& 
																tmpR2*( 6.0d0*solV (i,j-1,k,1:5,1) + 3.0d0*solV (i,j  ,k,1:5,1) - 1.0d0*solV (i,j-2,k,1:5,1) ) 		&
															)*invdy/8.0d0 
						endif

						if(tmpR1.lt.0.0d0.and.tmpR2.lt.0.0d0) then
							G(i,j,k,1:5,1) =(	tmpR1*( 6.0d0*solV (i,j+1,k,1:5,1) + 3.0d0*solV (i,j  ,k,1:5,1) - 1.0d0*solV (i,j+2,k,1:5,1) ) -	&
																tmpR2*( 6.0d0*solV (i,j  ,k,1:5,1) + 3.0d0*solV (i,j-1,k,1:5,1) - 1.0d0*solV (i,j+1,k,1:5,1) ) 		&
															)*invdy/8.0d0 
						endif

						if(tmpR1.lt.0.0d0.and.tmpR2.ge.0.0d0) then
							G(i,j,k,1:5,1) =(	tmpR1*( 6.0d0*solV (i,j+1,k,1:5,1) + 3.0d0*solV (i,j,k,1:5,1) - 1.0d0*solV (i,j+2,k,1:5,1) ) -	&
																tmpR2*( 6.0d0*solV (i,j-1,k,1:5,1) + 3.0d0*solV (i,j,k,1:5,1) - 1.0d0*solV (i,j-2,k,1:5,1) ) 		&
															)*invdy/8.0d0 
						endif

						if(tmpR1.ge.0.0d0.and.tmpR2.lt.0.0d0) then
							G(i,j,k,1:5,1) =(	tmpR1*( 6.0d0*solV (i,j,k,1:5,1) + 3.0d0*solV (i,j+1,k,1:5,1) - 1.0d0*solV (i,j-1,k,1:5,1) ) -	&
																tmpR2*( 6.0d0*solV (i,j,k,1:5,1) + 3.0d0*solV (i,j-1,k,1:5,1) - 1.0d0*solV (i,j+1,k,1:5,1) ) 		&
															)*invdy/8.0d0 
						endif

					else   ! Fourth order Central Differencing
						G(i,j,k,1:5,1) = 	( -				var(i,j+2,k,3)*solV (i,j+2,k,1:5,1)	&
																+	8.0d0*var(i,j+1,k,3)*solV (i,j+1,k,1:5,1)	&
																-	8.0d0*var(i,j-1,k,3)*solV (i,j-1,k,1:5,1)	&
																+				var(i,j-2,k,3)*solV (i,j-2,k,1:5,1)	&
															)*invdy*one12th	

					endif  
					!-------------------------------------------------------------
					Pe = Re*abs(solV(i,j,k,4,1))*dXi*Jac(k)/mu(i,j,k) 

					if (k .lt. 3 .or. k .gt. (nP(3)-2)) then
						H(i,j,k,1:5,1) = ( var(i,j,k+1,4)*solV (i,j,k+1,1:5,1)-var(i,j,k-1,4)*solV (i,j,k-1,1:5,1) )*invdXi*invJac(k)*haf
						
					elseif(Pe .gt. PezVal) then
					
						tmpR1  = ( var(i,j,k,4) + var(i,j,k+1,4) )*haf
						tmpR2  = ( var(i,j,k,4) + var(i,j,k-1,4) )*haf

						if(tmpR1.ge.0.0d0.and.tmpR2.ge.0.0d0) then
							H(i,j,k,1:5,1) =(	tmpR1*( 6.0d0*solV (i,j,k,1:5,1  ) + 3.0d0*solV (i,j,k+1,1:5,1) - 1.0d0*solV (i,j,k-1,1:5,1) ) -	& 
																tmpR2*( 6.0d0*solV (i,j,k-1,1:5,1) + 3.0d0*solV (i,j,k,1:5,1  ) - 1.0d0*solV (i,j,k-2,1:5,1) )		& 
															)*invdXi*invJac(k)/8.0d0
						endif

						if(tmpR1 .lt. 0.0d0 .and. tmpR2 .lt. 0.0d0) then
							H(i,j,k,1:5,1) =(	tmpR1*( 6.0d0*solV (i,j,k+1,1:5,1) + 3.0d0*solV (i,j,k  ,1:5,1) - 1.0d0*solV (i,j,k+2,1:5,1) ) -	&
																tmpR2*( 6.0d0*solV (i,j,k  ,1:5,1) + 3.0d0*solV (i,j,k-1,1:5,1) - 1.0d0*solV (i,j,k+1,1:5,1) ) 		&
															)*invdXi*invJac(k)/8.0d0
						endif

						if(tmpR1.lt.0.0d0.and.tmpR2.ge.0.0d0) then
							H(i,j,k,1:5,1) =(	tmpR1*( 6.0d0*solV (i,j,k+1,1:5,1) + 3.0d0*solV (i,j,k,1:5,1) - 1.0d0*solV (i,j,k+2,1:5,1) ) -	&
																tmpR2*( 6.0d0*solV (i,j,k-1,1:5,1) + 3.0d0*solV (i,j,k,1:5,1) - 1.0d0*solV (i,j,k-2,1:5,1) ) 		&
															)*invdXi*invJac(k)/8.0d0
						endif

						if(tmpR1.ge.0.0d0.and.tmpR2.lt.0.0d0) then
							H(i,j,k,1:5,1) =(	tmpR1*( 6.0d0*solV (i,j,k,1:5,1) + 3.0d0*solV (i,j,k+1,1:5,1) - 1.0d0*solV (i,j,k-1,1:5,1) ) -	&
																tmpR2*( 6.0d0*solV (i,j,k,1:5,1) + 3.0d0*solV (i,j,k-1,1:5,1) - 1.0d0*solV (i,j,k+1,1:5,1) ) 	&
															)*invdXi*invJac(k)/8.0d0 
						endif

					else   ! Fourth order Central Differencing
						H(i,j,k,1:5,1) = 	( -				var(i,j,k+2,4)*solV (i,j,k+2,1:5,1)	&
																+	8.0d0*var(i,j,k+1,4)*solV (i,j,k+1,1:5,1)	&
																-	8.0d0*var(i,j,k-1,4)*solV (i,j,k-1,1:5,1)	&
																+				var(i,j,k-2,4)*solV (i,j,k-2,1:5,1)	&
															)*invdXi*invJac(k)*one12th																						
					endif

				enddo
			enddo
		enddo
		!-------------------------------------------------------------------
		!===========Discretization of Diffusive terms in Flux Form==========
		do k=bs(3),es(3)
			do j=bn(2),en(2)
				do i=bn(1),en(1)
					
					Ut	=var(i,j,k,2)
					Vt	=var(i,j,k,3)
					Wt	=var(i,j,k,4)
					!-------------------------------------------------------------
					dudx	=(var(i+1,j,k,2)-var(i-1,j,k,2))*haf*invdx
					dudy	=(var(i,j+1,k,2)-var(i,j-1,k,2))*haf*invdy
					!-------------------------------------------------------------
					dvdx	=(var(i+1,j,k,3)-var(i-1,j,k,3))*haf*invdx
					dvdy	=(var(i,j+1,k,3)-var(i,j-1,k,3))*haf*invdy
					!-------------------------------------------------------------
					dwdx	=(var(i+1,j,k,4)-var(i-1,j,k,4))*haf*invdx
					dwdy	=(var(i,j+1,k,4)-var(i,j-1,k,4))*haf*invdy
					!-------------------------------------------------------------
					dTdx	=(var(i+1,j,k,6)-var(i-1,j,k,6))*haf*invdx
					dTdy	=(var(i,j+1,k,6)-var(i,j-1,k,6))*haf*invdy
					!-------------------------------------------------------------
					if(k==1) then
						dudz	=(var(i,j,2,2)-var(i,j,1,2))/dz1
						dvdz	=(var(i,j,2,3)-var(i,j,1,3))/dz1
						dwdz	=(var(i,j,2,4)-var(i,j,1,4))/dz1
						dTdz	=(var(i,j,2,6)-var(i,j,1,6))/dz1
					elseif(k==nP(3)) then
						dudz	=(var(i,j,nP(3),2)-var(i,j,nP(3)-1,2))/dzN
						dvdz	=(var(i,j,nP(3),3)-var(i,j,nP(3)-1,3))/dzN
						dwdz	=(var(i,j,nP(3),4)-var(i,j,nP(3)-1,4))/dzN
						dTdz	=(var(i,j,nP(3),6)-var(i,j,nP(3)-1,6))/dzN
					else
						dudz	=(var(i,j,k+1,2)-var(i,j,k-1,2))*haf*invdXi*invJac(k)
						dvdz	=(var(i,j,k+1,3)-var(i,j,k-1,3))*haf*invdXi*invJac(k)
						dwdz	=(var(i,j,k+1,4)-var(i,j,k-1,4))*haf*invdXi*invJac(k)
						dTdz	=(var(i,j,k+1,6)-var(i,j,k-1,6))*haf*invdXi*invJac(k)
					endif
					!-------------------------------------------------------------
					divg	=dudx + dvdy + dwdz
					muT		=mu(i,j,k)
					ktT		=kt(i,j,k)
					PressT=var(i,j,k,7)
					
					dissF	=ggM*invRe*(  two3rd*muT*Ut*divg - two*muT*Ut*dudx - muT*Vt*(dvdx + dudy) - muT*Wt*(dwdx + dudz) )
					dissG	=ggM*invRe*( -muT*Ut*(dvdx + dudy) + two3rd*muT*Vt*divg - two*muT*Vt*dvdy - muT*Wt*(dwdy + dvdz) )
					dissH	=ggM*invRe*( -muT*Ut*(dwdx + dudz) - muT*Vt*(dwdy + dvdz) + two3rd*muT*Wt*divg - two*muT*Wt*dwdz )
					
					!-------------------------------------------------------------
					fD(i,j,k,1)	=	zer
					
					fD(i,j,k,2)	=	PressT + two3rd*muT*invRe*divg - two*muT*invRe*dudx 

					fD(i,j,k,3)	=	-muT*invRe*( dvdx + dudy )

					fD(i,j,k,4)	=	-muT*invRe*( dwdx + dudz )
					
					fD(i,j,k,5)	=	-gama*ktT*invRe*invPr*dTdx + ggM*Ut*PressT + dissF

					!-------------------------------------------------------------
					gD(i,j,k,1)	=	zer
					
					gD(i,j,k,2)	=	-muT*invRe*(dvdx + dudy) 

					gD(i,j,k,3)	=	PressT + two3rd*muT*invRe*divg - two*muT*invRe*dvdy

					gD(i,j,k,4)	=	-muT*invRe*(dwdy + dvdz)

					gD(i,j,k,5)	=	-gama*ktT*invRe*invPr*dTdy + ggM*Vt*PressT + dissG

					!-------------------------------------------------------------
					hD(i,j,k,1)	=	zer
					
					hD(i,j,k,2)	=	-muT*invRe*(dwdx + dudz)

					hD(i,j,k,3)	=	-muT*invRe*(dwdy + dvdz)

					hD(i,j,k,4)	=	PressT + two3rd*muT*invRe*divg - two*muT*invRe*dwdz 

					hD(i,j,k,5)	=	-gama*ktT*invRe*invPr*dTdz + ggM*Wt*PressT + dissH
					!-------------------------------------------------------------
				enddo 
			enddo
		enddo
		
		!-------------------------------------------------------------------
		call mpi_sendrecv(fD(en(1),bn(2),bn(3),1),1,yz1p5v,des(1),50,fD(bn(1)-1,bn(2),bn(3),1),1,yz1p5v,src(1),50,mpi_comm_world,STATUS,ierr)
		call mpi_sendrecv(fD(bn(1),bn(2),bn(3),1),1,yz1p5v,src(1),50,fD(en(1)+1,bn(2),bn(3),1),1,yz1p5v,des(1),50,mpi_comm_world,STATUS,ierr)

		call mpi_sendrecv(gD(bn(1),en(2),bn(3),1),1,xz1p5v,des(2),50,gD(bn(1),bn(2)-1,bn(3),1),1,xz1p5v,src(2),50,mpi_comm_world,STATUS,ierr)
		call mpi_sendrecv(gD(bn(1),bn(2),bn(3),1),1,xz1p5v,src(2),50,gD(bn(1),en(2)+1,bn(3),1),1,xz1p5v,des(2),50,mpi_comm_world,STATUS,ierr)

		call mpi_sendrecv(hD(bn(1),bn(2),en(3),1),1,xy1p5v,des(3),50,hD(bn(1),bn(2),bn(3)-1,1),1,xy1p5v,src(3),50,mpi_comm_world,STATUS,ierr)
		call mpi_sendrecv(hD(bn(1),bn(2),bn(3),1),1,xy1p5v,src(3),50,hD(bn(1),bn(2),en(3)+1,1),1,xy1p5v,des(3),50,mpi_comm_world,STATUS,ierr)
		!-------------------------------------------------------------------

		do k=bn(3),en(3)
			do j=bn(2),en(2)
				do i=bn(1),en(1)
					F(i,j,k,:,2) = ( fD(i+1,j,k,:) - fD(i-1,j,k,:) )*haf*invdx

					G(i,j,k,:,2) = ( gD(i,j+1,k,:) - gD(i,j-1,k,:) )*haf*invdy

					H(i,j,k,:,2) = ( hD(i,j,k+1,:) - hD(i,j,k-1,:) )*haf*invdXi*invJac(k)
				enddo
			enddo
		enddo
		!---------------------------------------------------------------
		do k=bn(3),en(3)
			do j=bn(2),en(2)
				do i=bn(1),en(1)
					srcV(i,j,k,1) 	= zer
					srcV(i,j,k,2) 	= haf*(tauxzL-tauxzU)*stat1Spc(k,1)
					srcV(i,j,k,3:4) = zer 
					srcV(i,j,k,5) 	= ggM*var(i,j,k,2)*haf*(tauxzL-tauxzU)*stat1Spc(k,1)
					
					R(i,j,k,:,1) = srcV(i,j,k,:) - (F(i,j,k,:,1) + F(i,j,k,:,2) + G(i,j,k,:,1) + G(i,j,k,:,2) + H(i,j,k,:,1) + H(i,j,k,:,2))

					if (iter.eq.1) then
						solV(i,j,k,:,2)	= solV(i,j,k,:,1)
						solV(i,j,k,:,1)	= solV(i,j,k,:,1) + haf*dt*R(i,j,k,:,1)
					endif

					if (iter.eq.2) then
						solV(i,j,k,:,1) = solV(i,j,k,:,2) + dt*( R(i,j,k,:,1) )
					endif
					
				enddo
			enddo
		enddo
		
		!-------------------------------------------------------------------	
		call mpi_sendrecv(solV(en(1)-1,bn(2),bn(3),1,1),1,yz2p5v,des(1),50,solV(bn(1)-2,bn(2),bn(3),1,1),1,yz2p5v,src(1),50,mpi_comm_world,STATUS,ierr)
		call mpi_sendrecv(solV(bn(1)  ,bn(2),bn(3),1,1),1,yz2p5v,src(1),50,solV(en(1)+1,bn(2),bn(3),1,1),1,yz2p5v,des(1),50,mpi_comm_world,STATUS,ierr)

		call mpi_sendrecv(solV(bn(1),en(2)-1,bn(3),1,1),1,xz2p5v,des(2),50,solV(bn(1),bn(2)-2,bn(3),1,1),1,xz2p5v,src(2),50,mpi_comm_world,STATUS,ierr)
		call mpi_sendrecv(solV(bn(1),bn(2)  ,bn(3),1,1),1,xz2p5v,src(2),50,solV(bn(1),en(2)+1,bn(3),1,1),1,xz2p5v,des(2),50,mpi_comm_world,STATUS,ierr)

		call mpi_sendrecv(solV(bn(1),bn(2),en(3)-1,1,1),1,xy2p5v,des(3),50,solV(bn(1),bn(2),bn(3)-2,1,1),1,xy2p5v,src(3),50,mpi_comm_world,STATUS,ierr)
		call mpi_sendrecv(solV(bn(1),bn(2),bn(3)  ,1,1),1,xy2p5v,src(3),50,solV(bn(1),bn(2),en(3)+1,1,1),1,xy2p5v,des(3),50,mpi_comm_world,STATUS,ierr)
		!-------------------------------------------------------------------
		
		! Finding Primitive Variables
		do k=bn(3),en(3)
			do j=bn(2),en(2)
				do i=bn(1),en(1)
					var(i,j,k,1  )	= solV(i,j,k,1,1)
					var(i,j,k,2:5)	= solV(i,j,k,2:5,1) / solV(i,j,k,1,1)
					var(i,j,k,6)		= var(i,j,k,5) - ggM*haf*sum(var(i,j,k,2:4)**two)
					var(i,j,k,7)		= var(i,j,k,1)*var(i,j,k,6) / (gama*Mac**two)
				enddo
			enddo
		enddo
		
		!-------------------------------------------------------------------
		call mpi_sendrecv(var(en(1),bn(2),bn(3),1),1,yz1p7v,des(1),50,var(bn(1)-1,bn(2),bn(3),1),1,yz1p7v,src(1),50,mpi_comm_world,STATUS,ierr)
		call mpi_sendrecv(var(bn(1),bn(2),bn(3),1),1,yz1p7v,src(1),50,var(en(1)+1,bn(2),bn(3),1),1,yz1p7v,des(1),50,mpi_comm_world,STATUS,ierr)

		call mpi_sendrecv(var(bn(1),en(2),bn(3),1),1,xz1p7v,des(2),50,var(bn(1),bn(2)-1,bn(3),1),1,xz1p7v,src(2),50,mpi_comm_world,STATUS,ierr)
		call mpi_sendrecv(var(bn(1),bn(2),bn(3),1),1,xz1p7v,src(2),50,var(bn(1),en(2)+1,bn(3),1),1,xz1p7v,des(2),50,mpi_comm_world,STATUS,ierr)

		call mpi_sendrecv(var(bn(1),bn(2),en(3),1),1,xy1p7v,des(3),50,var(bn(1),bn(2),bn(3)-1,1),1,xy1p7v,src(3),50,mpi_comm_world,STATUS,ierr)
		call mpi_sendrecv(var(bn(1),bn(2),bn(3),1),1,xy1p7v,src(3),50,var(bn(1),bn(2),en(3)+1,1),1,xy1p7v,des(3),50,mpi_comm_world,STATUS,ierr)
		!-------------------------------------------------------------------
				
		! Boundary Conditions
		do j=bs(2),es(2)
			do i=bs(1),es(1)
			
				if(bottom) then
					!tmpR1=(x1(nP(1))+x2(nP(2)))/8.0
					!tmpR1=(x1(nP(1))+x2(nP(2)))/16.0
					tmpR1=x2(nP(2))/8.0
					!if(mod(floor((x1(i)-x2(j)+x2(nP(2)))/tmpR1),2)==0) then
					!if(mod(floor((0.364*x1(i)-x2(j)+x2(nP(2)))/tmpR1),2)==0) then
					if(mod(floor((x2(j)-tmpR1*sin(x1(i)))/tmpR1)+1,3)==0) then
						pattern(i,j,1) = 1
						var(i,j,1,2) = (0.0 - f1d2c1*var(i,j,2,2) - f1d2c2*var(i,j,3,2) ) / f1d2c0
						var(i,j,1,3) = (0.0 - f1d2c1*var(i,j,2,3) - f1d2c2*var(i,j,3,3) ) / f1d2c0
					else
						pattern(i,j,1) = 0
						var(i,j,1,2) = 0.0d0
						var(i,j,1,3) = 0.0d0
					endif
					
					var(i,j,1,4) = 0.0d0
					var(i,j,1,6) = 0.3d0
					var(i,j,1,5) = var(i,j,1,6) +	ggM*haf*sum(var(i,j,1,2:4)**two)
					var(i,j,1,1) = gama*Mac**two*var(i,j,1,7) / var(i,j,1,6)
				endif
				
				if(top) then
					var(i,j,nP(3),2:4)= 0.0d0
					var(i,j,nP(3),6) 	= 1.0d0
					var(i,j,nP(3),5) 	= var(i,j,nP(3),6) +	ggM*haf*sum(var(i,j,nP(3),2:4)**2)
					var(i,j,nP(3),1) 	= gama*Mac**two*var(i,j,nP(3),7) / var(i,j,nP(3),6)
				endif	
								
			enddo
		enddo
				
		mu(bn(1)-1:en(1)+1,bn(2)-1:en(2)+1,bs(3):es(3)) = var(bn(1)-1:en(1)+1,bn(2)-1:en(2)+1,bs(3):es(3),6)**0.7d0
				
		do j=bn(2),en(2)
			do i=bn(1),en(1)	
				!-------------------------Pressure BC---------------------------
				!Simplifying z-momentum equation assuming impervious wall
				if(bottom) then
					dudxdz	= (	(		f1d2c0*mu(i+1,j,1)*var(i+1,j,1,2)	&
												+ f1d2c1*mu(i+1,j,2)*var(i+1,j,2,2)	&
												+ f1d2c2*mu(i+1,j,3)*var(i+1,j,3,2)	&
											) -																		&
											(		f1d2c0*mu(i-1,j,1)*var(i-1,j,1,2)	&
												+ f1d2c1*mu(i-1,j,2)*var(i-1,j,2,2)	&
												+ f1d2c2*mu(i-1,j,3)*var(i-1,j,3,2)	&
											)																			&
										) *haf*invdx
										
					dvdydz	= (	(		f1d2c0*mu(i,j+1,1)*var(i,j+1,1,3)	&
												+ f1d2c1*mu(i,j+1,2)*var(i,j+1,2,3)	&
												+ f1d2c2*mu(i,j+1,3)*var(i,j+1,3,3)	&
											) -																		&
											(		f1d2c0*mu(i,j-1,1)*var(i,j-1,1,3)	&
												+ f1d2c1*mu(i,j-1,2)*var(i,j-1,2,3)	&
												+ f1d2c2*mu(i,j-1,3)*var(i,j-1,3,3)	&
											)																			&
										) *haf*invdy
										
					d2wdz2	= 	f2d1c0*mu(i,j,1)*var(i,j,1,4)	&
										+ f2d1c1*mu(i,j,2)*var(i,j,2,4)	&
										+ f2d1c2*mu(i,j,3)*var(i,j,3,4)
					
					drhow2dz= 	f1d2c0*var(i,j,1,1)*var(i,j,1,4)**two	&
										+	f1d2c1*var(i,j,2,1)*var(i,j,2,4)**two	&
										+	f1d2c2*var(i,j,3,1)*var(i,j,3,4)**two
					
					lhs			= invRe*(-dudxdz -dvdydz + two3rd*(dudxdz + dvdydz - two*d2wdz2)) + drhow2dz
					
					var(i,j,1,7) = -(lhs + f1d2c1*var(i,j,2,7) + f1d2c2*var(i,j,3,7) ) / f1d2c0
				endif
				
				if(top) then								
					dudxdz	= (	(		b1d2c0*mu(i+1,j,nP(3)  )*var(i+1,j,nP(3)  ,2)	&
												+ b1d2c1*mu(i+1,j,nP(3)-1)*var(i+1,j,nP(3)-1,2)	&
												+ b1d2c2*mu(i+1,j,nP(3)-2)*var(i+1,j,nP(3)-2,2)	&
											)	-	&
											(		b1d2c0*mu(i-1,j,nP(3)  )*var(i-1,j,nP(3)  ,2)	&
												+ b1d2c1*mu(i-1,j,nP(3)-1)*var(i-1,j,nP(3)-1,2)	&
												+ b1d2c2*mu(i-1,j,nP(3)-2)*var(i-1,j,nP(3)-2,2)	&
											)		&
										) *haf*invdx
										
					dvdydz	= (	(		b1d2c0*mu(i,j+1,nP(3)  )*var(i,j+1,nP(3)  ,3)	&
												+	b1d2c1*mu(i,j+1,nP(3)-1)*var(i,j+1,nP(3)-1,3)	&
												+ b1d2c2*mu(i,j+1,nP(3)-2)*var(i,j+1,nP(3)-2,3)	&
											)	-																								&
											(		b1d2c0*mu(i,j-1,nP(3)  )*var(i,j-1,nP(3)  ,3)	&
												+ b1d2c1*mu(i,j-1,nP(3)-1)*var(i,j-1,nP(3)-1,3)	&
												+ b1d2c2*mu(i,j-1,nP(3)-2)*var(i,j-1,nP(3)-2,3)	&
											)																									&
										)	*haf*invdy
										
					d2wdz2	= 	b2d1c0*mu(i,j,nP(3)  )*var(i,j,nP(3)  ,4)	&
										+	b2d1c1*mu(i,j,nP(3)-1)*var(i,j,nP(3)-1,4)	&
										+	b2d1c2*mu(i,j,nP(3)-2)*var(i,j,nP(3)-2,4)
					
					drhow2dz= 	b1d2c0*var(i,j,nP(3)  ,1)*var(i,j,nP(3)  ,4)**two	&
										+	b1d2c1*var(i,j,nP(3)-1,1)*var(i,j,nP(3)-1,4)**two	&
										+	b1d2c2*var(i,j,nP(3)-2,1)*var(i,j,nP(3)-2,4)**two
					
					lhs			= invRe*(-dudxdz -dvdydz + two3rd*(dudxdz + dvdydz - two*d2wdz2)) + drhow2dz
					
					var(i,j,nP(3),7) = -(lhs + b1d2c1*var(i,j,nP(3)-1,7) + b1d2c2*var(i,j,nP(3)-2,7) ) / b1d2c0
					!-------------------------------------------------------------
				endif
				
			enddo
		enddo
		kt(bs(1):es(1),bs(2):es(2),bs(3):es(3)) = var(bs(1):es(1),bs(2):es(2),bs(3):es(3),6)**0.7d0
		
	enddo
	!=========================RK2 Sub-steps (ends)========================
	!=====================Check for unphysical values=====================
	do k=bs(3),es(3)
		do j=bs(2),es(2)
			do i=bs(1),es(1)
			
				if (var(i,j,k,1) < 0.0d0 .or. var(i,j,k,6) < 0.0d0) then
					write(6,102) nstep,id,i,j,k,var(i,j,k,1),var(i,j,k,6)
					isBlown=.true.
				endif
				
			enddo
		enddo
	enddo
	
	if(isBlown) then
	
		WRITE(filename,'(a,i3.3,a)')"../debug/3D",id+1,".dat"
		OPEN (UNIT=10,FILE=filename)

		write(10,103) 'TITLE ="',time,nstep,'"'
		write(10,*) 'Variables ="x","y","z","Rho","U","V","W","E","T","P"'
		write(10,104) 'ZONE k=',es(3)-bs(3)+1,',j=',es(2)-bs(2)+1,',i=',es(1)-bs(1)+1,',DATAPACKING="POINT"'

		do k=bs(3),es(3)
			write(10,*)
			do j=bs(2),es(2)
				write(10,*)
				do i=bs(1),es(1)
					write(10,116) x1(i),x2(j),x3(k),var(i,j,k,1),var(i,j,k,2),var(i,j,k,3),var(i,j,k,4)	&
																											,var(i,j,k,5),var(i,j,k,6),var(i,j,k,7)
				end do
			end do
		enddo
		close(10)	
		stop
	endif
	
	102 format(I10,4(I4),2(F8.4)) 
	116 format(3(1X,F7.4),7(1X,E12.5))
	
	!==================Write output files (starts)========================
	if(mod(nstep,100).eq.0) then
		
		if(mod(nstep/100,2).eq.0) then
			WRITE(filename,'(a,i3.3,a)')"../output/3D",id+1,".dat"
		else	
			WRITE(filename,'(a,i3.3,a)')"../output/3D",id+1,"_.dat"
		endif
			
		OPEN (UNIT=10,FILE=filename)

		write(10,103) 'TITLE ="',time,nstep,'"'
		write(10,*) 'Variables ="x","y","z","Rho","U","V","W","E","T","P","Pattern"'
		write(10,104) 'ZONE k=',es(3)-bs(3)+1,',j=',es(2)-bs(2)+1,',i=',es(1)-bs(1)+1,',DATAPACKING="POINT"'

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
		103	format(A,F20.10,I9,A)
		104	format(A,I3,A,I3,A,I3,A)

		!-------------------------------------------------------------------
		call mpi_gather(bs(1),3,mpi_integer,bb(1,0),3,mpi_integer,0,mpi_comm_world,ierr)
		call mpi_gather(es(1),3,mpi_integer,ee(1,0),3,mpi_integer,0,mpi_comm_world,ierr)
		!-----------------------Write Domain Info---------------------------
		if(master) then
			open(unit=10,file="../output/domInfo.txt")
			write(10,105) nP(1),nP(2),nP(3)
			write(10,106) nproc
			do i=0,nproc-1
				write(10,107) bb(1,i),ee(1,i),bb(2,i),ee(2,i),bb(3,i),ee(3,i)
			enddo
			close(10)
			105 format(3(i4))
			106 format(i3)
			107 format(6(i4,1x))	
		endif
		
	endif
	!==================Write output files (ends)==========================
	
	!==================Write Diagnostic Files (starts)====================
	
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
	!====================Write Diagnostic Files (ends)====================

	!==================Wall shear stress (starts)=========================
	if(mod(nstep,istat)==0) then
		!-------------------------------------------------------------------
		tmpR1=0.0d0
		if(bottom) then
		
			do j=bn(2),en(2)
				do i=bn(1),en(1)
					dwdx = mu(i,j,1) * (var(i+1,j,1,4) - var(i-1,j,1,4)) *haf*invdx !zero for impervious bottom wall

					dudz = mu(i,j,1)*	(		f1d2c0*var(i,j,1,2)	&
															+ f1d2c1*var(i,j,2,2) &
															+ f1d2c2*var(i,j,3,2)	&
														)

					tmpR1 = tmpR1 + (dwdx + dudz)*invRe
				enddo
			enddo
			
		endif
		
		CALL MPI_Reduce(tmpR1,tauxzL,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
		tauxzL= tauxzL/dble((nP(1)-2)*(nP(2)-2))
		
		call mpi_bcast(tauxzL,1,mpi_double_precision,0,mpi_comm_world,ierr)
		!-------------------------------------------------------------------
		
		tmpR1=0.0d0
		if(top) then

			do j=bn(2),en(2)
				do i=bn(1),en(1)
					dwdx = mu(i,j,nP(3)) * (var(i+1,j,nP(3),4) - var(i-1,j,nP(3),4)) *haf*invdx !zero for impervious top wall
					
					dudz = mu(i,j,nP(3)) * (		b1d2c0*var(i,j,nP(3)  ,2)	&
																		+ b1d2c1*var(i,j,nP(3)-1,2) &
																		+	b1d2c2*var(i,j,nP(3)-2,2)	&
																	)

					tmpR1 = tmpR1 + (dwdx + dudz)*invRe
				enddo
			enddo
			
		endif
		
		CALL MPI_Reduce(tmpR1,tauxzU,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
		tauxzU= tauxzU/dble(((nP(1)-2)*(nP(2)-2)))
		
		call mpi_bcast(tauxzU,1,mpi_double_precision,0,mpi_comm_world,ierr)
		!-------------------------------------------------------------------

		if(master) then
			open(10,file="../output/tauWall.dat",access="append")
			write(10,110)time, tauxzL, tauxzU
			close(10) 
		endif
		110	format(3(1X,F15.7))
				
	end if
	!==================Wall shear stress (ends)===========================
	
	!==================Spatial Statistics (starts)========================
	
	if(mod(nstep,istat)==0) then
		nXY=nP(1)*nP(2)
		
		do k=1,nP(3)
			loc=0.0d0
			
			if(k >= bs(3) .and. k <= es(3)) then
				do j=bs(2),es(2)
					do i=bs(1),es(1)
						loc(i,j,1:7)	= var(i,j,k,1:7)
					enddo
				enddo				
			endif
			
			do p=1,7
				call mpi_reduce(sum(loc(:,:,p)   )/nXY,stat1Spc(k,p),1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
				call mpi_reduce(sum(loc(:,:,p)**2)/nXY,stat2Spc(k,p),1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
			enddo
			
			call mpi_reduce(sum(loc(:,:,1)*loc(:,:,4))/nXY,corrSpc(k,1),1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
			call mpi_reduce(sum(loc(:,:,1)*loc(:,:,6))/nXY,corrSpc(k,2),1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
			call mpi_reduce(sum(loc(:,:,2)*loc(:,:,3))/nXY,corrSpc(k,3),1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
			call mpi_reduce(sum(loc(:,:,2)*loc(:,:,4))/nXY,corrSpc(k,4),1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
			call mpi_reduce(sum(loc(:,:,3)*loc(:,:,4))/nXY,corrSpc(k,5),1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
			call mpi_reduce(sum(loc(:,:,3)*loc(:,:,6))/nXY,corrSpc(k,6),1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
			call mpi_reduce(sum(loc(:,:,6)*loc(:,:,4))/nXY,corrSpc(k,7),1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
			call mpi_reduce(sum(loc(:,:,7)*loc(:,:,4))/nXY,corrSpc(k,8),1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
			call mpi_reduce(sum(loc(:,:,7)*loc(:,:,6))/nXY,corrSpc(k,9),1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
		enddo
		
		call mpi_bcast(stat1Spc(1,1),nP(3),mpi_double_precision,0,mpi_comm_world,ierr)

		if(master) then
			!-----------------------------------------------------------------
			open(unit=10,file="../output/moment1st.dat",access="append")
			write(10,111,advance="no") time
			do k=1,nP(3)
				write(10,112,advance="no") stat1Spc(k,1),stat1Spc(k,2),stat1Spc(k,3),stat1Spc(k,4)	&
																								,stat1Spc(k,5),stat1Spc(k,6),stat1Spc(k,7)
			enddo
			write (10, *)
			close(10)
			111	format(1(1X,F15.7))
			112	format(7(1X,F15.7))
			!-----------------------------------------------------------------
			open(unit=10,file="../output/moment2nd.dat",access="append")
			write(10,111,advance="no") time
			do k=1,nP(3)
				write(10,112,advance="no") stat2Spc(k,1),stat2Spc(k,2),stat2Spc(k,3),stat2Spc(k,4)	&
																								,stat2Spc(k,5),stat2Spc(k,6),stat2Spc(k,7)
			enddo
			write (10, *)
			close(10)
			!-----------------------------------------------------------------
			open(unit=10,file="../output/momentJoint.dat",access="append")
			write(10,111,advance="no")time
			do k=1,nP(3)
				write(10,113,advance="no") corrSpc(k,1),corrSpc(k,2),corrSpc(k,3)	&
																	,corrSpc(k,4),corrSpc(k,5),corrSpc(k,6)	&
																	,corrSpc(k,7),corrSpc(k,8),corrSpc(k,9)
			enddo
			write (10,*)
			close(10)
			113	format(9(1X,F15.7))
			!-----------------------------------------------------------------
		endif
		
	endif
	!====================Spatial Statistics (ends)========================
	
	if(master .and. mod(nstep,10)==0) then
		call cpu_time(tEnd)
		write(6,114) nstep,time," | CPU Time:",tEnd-tStart,"sec for 10 iters | ", dateTime()
		write(6,'(A)')"-----------------------------------------------------------------------------------"
	endif
	114  format(I10,F10.5,A,F6.3,2(1X,A))
	
END DO
!========================Main time loop (ends)==========================

CALL mpi_type_free(yz1p,ierr)
CALL mpi_type_free(xz1p,ierr)
CALL mpi_type_free(xy1p,ierr)

CALL mpi_type_free(yz1p5v,ierr)
CALL mpi_type_free(xz1p5v,ierr)
CALL mpi_type_free(xy1p5v,ierr)

CALL mpi_type_free(yz2p5v,ierr)
CALL mpi_type_free(xz2p5v,ierr)
CALL mpi_type_free(xy2p5v,ierr)

CALL mpi_type_free(yz1p7v,ierr)
CALL mpi_type_free(xz1p7v,ierr)
CALL mpi_type_free(xy1p7v,ierr)

CALL mpi_finalize(ierr)
!=======================================================================
if(master) write(*,*) "Compressible Channel Flow Solver exited at :",dateTime()
!=======================================================================
end program chanCompr
