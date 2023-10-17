module module_MPI
use module_read
use module_comp
use module_const, only: i,k,nP
implicit none

integer,protected		:: id
integer,dimension(3),protected		:: bs,es,bn,en
integer,dimension(3),private			:: siz,des,src
integer,private		:: liney,yz1p, xz1p, xy1p, yz1p5v, xz1p5v, xy1p5v, yz2p7v, xz2p7v, xy2p7v
integer,private		:: lineyN, yzN1p, xzN1p, xyN1p, yzN2p, xzN2p, xyN2p, yz2p5v, xz2p5v, xy2p5v
integer,private		:: ierr,STATUS(mpi_status_size)
integer,allocatable,dimension(:,:),public 	::bb,ee
logical,protected 	:: bottom, top, master

contains
subroutine setupMPI()
implicit none

integer,private		:: comm3d,nproc,sizDP
integer,dimension(3),private		:: myDim,myCoord
integer,allocatable,dimension(:,:),private 	::tArr
logical,private 	:: myperiod(3),east, west, north, south

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

CALL mpi_type_hvector(7, 1 , product(nP+2)*sizDP , yzN2p , yz2p7v ,ierr)
CALL mpi_type_hvector(7, 1 , product(nP+2)*sizDP , xzN2p , xz2p7v ,ierr)
CALL mpi_type_hvector(7, 1 , product(nP+2)*sizDP , xyN2p , xy2p7v ,ierr)

CALL mpi_type_commit(yz2p7v,ierr)
CALL mpi_type_commit(xz2p7v,ierr)
CALL mpi_type_commit(xy2p7v,ierr)

CALL mpi_barrier(mpi_comm_world,ierr)

end subroutine setupMPI

subroutine share1P()
implicit none

call mpi_sendrecv(fD(en(1),bn(2),bn(3),1),1,yz1p5v,des(1),50,fD(bn(1)-1,bn(2),bn(3),1),1,yz1p5v,src(1),50,mpi_comm_world,STATUS,ierr)
call mpi_sendrecv(fD(bn(1),bn(2),bn(3),1),1,yz1p5v,src(1),50,fD(en(1)+1,bn(2),bn(3),1),1,yz1p5v,des(1),50,mpi_comm_world,STATUS,ierr)

call mpi_sendrecv(gD(bn(1),en(2),bn(3),1),1,xz1p5v,des(2),50,gD(bn(1),bn(2)-1,bn(3),1),1,xz1p5v,src(2),50,mpi_comm_world,STATUS,ierr)
call mpi_sendrecv(gD(bn(1),bn(2),bn(3),1),1,xz1p5v,src(2),50,gD(bn(1),en(2)+1,bn(3),1),1,xz1p5v,des(2),50,mpi_comm_world,STATUS,ierr)

call mpi_sendrecv(hD(bn(1),bn(2),en(3),1),1,xy1p5v,des(3),50,hD(bn(1),bn(2),bn(3)-1,1),1,xy1p5v,src(3),50,mpi_comm_world,STATUS,ierr)
call mpi_sendrecv(hD(bn(1),bn(2),bn(3),1),1,xy1p5v,src(3),50,hD(bn(1),bn(2),en(3)+1,1),1,xy1p5v,des(3),50,mpi_comm_world,STATUS,ierr)

end subroutine share1P

subroutine share2P(var,noOfVar)
implicit none
double precision,dimension(:,:,:,:),intent(in):: var
integer,intent(in)::noOfVar

if(noOfVar==5) then
	yz2p_v = yz2p5v
	xz2p_v = xz2p5v
	xy2p_v = xy2p5v
elseif(noOfVar==7) then
	yz2p_v = yz2p7v
	xz2p_v = xz2p7v
	xy2p_v = xy2p7v
else
		write(6,'(A,1X,I2,1X,A)')"No in-built mechanism to transfer",noOfVar,"planes"
endif		
	
call mpi_sendrecv(var(en(1)-1,bn(2),bn(3),1),1,yz2p_v,des(1),50,var(bn(1)-2,bn(2),bn(3),1),1,yz2p_v,src(1),50,mpi_comm_world,STATUS,ierr)
call mpi_sendrecv(var(bn(1)  ,bn(2),bn(3),1),1,yz2p_v,src(1),50,var(en(1)+1,bn(2),bn(3),1),1,yz2p_v,des(1),50,mpi_comm_world,STATUS,ierr)

call mpi_sendrecv(var(bn(1),en(2)-1,bn(3),1),1,xz2p_v,des(2),50,var(bn(1),bn(2)-2,bn(3),1),1,xz2p_v,src(2),50,mpi_comm_world,STATUS,ierr)
call mpi_sendrecv(var(bn(1),bn(2)  ,bn(3),1),1,xz2p_v,src(2),50,var(bn(1),en(2)+1,bn(3),1),1,xz2p_v,des(2),50,mpi_comm_world,STATUS,ierr)

call mpi_sendrecv(var(bn(1),bn(2),en(3)-1,1),1,xy2p_v,des(3),50,var(bn(1),bn(2),bn(3)-2,1),1,xy2p_v,src(3),50,mpi_comm_world,STATUS,ierr)
call mpi_sendrecv(var(bn(1),bn(2),bn(3)  ,1),1,xy2p_v,src(3),50,var(bn(1),bn(2),en(3)+1,1),1,xy2p_v,des(3),50,mpi_comm_world,STATUS,ierr)

end subroutine share2P

subroutine cleanupMPI()
implicit none
CALL mpi_type_free(yz1p,ierr)
CALL mpi_type_free(xz1p,ierr)
CALL mpi_type_free(xy1p,ierr)

CALL mpi_type_free(yz1p5v,ierr)
CALL mpi_type_free(xz1p5v,ierr)
CALL mpi_type_free(xy1p5v,ierr)

CALL mpi_type_free(yz2p5v,ierr)
CALL mpi_type_free(xz2p5v,ierr)
CALL mpi_type_free(xy2p5v,ierr)

CALL mpi_type_free(yz2p7v,ierr)
CALL mpi_type_free(xz2p7v,ierr)
CALL mpi_type_free(xy2p7v,ierr)

CALL mpi_finalize(ierr)
end subroutine cleanupMPI

end module module_MPI
