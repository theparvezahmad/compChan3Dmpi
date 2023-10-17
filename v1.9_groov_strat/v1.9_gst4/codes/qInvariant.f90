program qInvariant
use comData, only : nP,haf
implicit none
!=======================================================================
integer						:: i,j,k,jcnt
double precision	:: Jt(3,3),As(2),S(3,3),R(3,3),dx(3,2),x(nP(1)),y(nP(2)),z(nP(3)),P
double precision	:: Ut(nP(1),nP(2),nP(3),3),Q(nP(1),nP(2),nP(3))
!-----------------------------------------------------------------------
open (unit=10,file="../output/fluc.dat",status="old")
read (10,*)
read (10,*)
read (10,*)

do k=1,nP(3)
	read(10,*)
	do j=1,nP(2)
		read(10,*)
		do i=1,nP(1)
			read(10,'(3(1X,F10.6),23X,3(1X,E22.15))') x(i),y(j),z(k),Ut(i,j,k,1),Ut(i,j,k,2),Ut(i,j,k,3)
		end do
	end do
end do
close(10)
write(6,*) "File fluc.dat read successfully"
!-----------------------------------------------------------------------
do k=2,nP(3)-1
	do j=2,nP(2)-1
		do i=2,nP(1)-1
		
			dx(1,1)=x(i+1)-x(i)
			dx(2,1)=y(j+1)-y(j)
			dx(3,1)=z(k+1)-z(k)

			dx(1,2)=x(i)-x(i-1)
			dx(2,2)=y(j)-y(j-1)
			dx(3,2)=z(k)-z(k-1)

			do jcnt=1,3
				! x-derivatives        
				Jt(jcnt,1)=((dx(1,2)/dx(1,1))*(Ut(i+1,j,k,jcnt)-Ut(i,j,k,jcnt))+(dx(1,1)/dx(1,2))*(Ut(i,j,k,jcnt)-Ut(i-1,j,k,jcnt)))/(dx(1,1)+dx(1,2))
				! y-derivatives
				Jt(jcnt,2)=((dx(2,2)/dx(2,1))*(Ut(i,j+1,k,jcnt)-Ut(i,j,k,jcnt))+(dx(2,1)/dx(2,2))*(Ut(i,j,k,jcnt)-Ut(i,j-1,k,jcnt)))/(dx(2,1)+dx(2,2))
				! z-derivatives
				Jt(jcnt,3)=((dx(3,2)/dx(3,1))*(Ut(i,j,k+1,jcnt)-Ut(i,j,k,jcnt))+(dx(3,1)/dx(3,2))*(Ut(i,j,k,jcnt)-Ut(i,j,k-1,jcnt)))/(dx(3,1)+dx(3,2))
			enddo

			S(1,1)=Jt(1,1)
			S(1,2)=haf*(Jt(1,2)+Jt(2,1))
			S(1,3)=haf*(Jt(1,3)+Jt(3,1))

			S(2,1)=haf*(Jt(2,1)+Jt(1,2))
			S(2,2)=Jt(2,2)
			S(2,3)=haf*(Jt(3,2)+Jt(2,3))

			S(3,1)=haf*(Jt(3,1)+Jt(1,3))
			S(3,2)=haf*(Jt(3,2)+Jt(2,3))
			S(3,3)=Jt(3,3)

			! Note s(ii,jj) = s(jj,ii)

			R(1,1)=0.0
			R(1,2)=haf*(Jt(1,2)-Jt(2,1))
			R(1,3)=haf*(Jt(1,3)-Jt(3,1))

			R(2,1)=haf*(Jt(2,1)-Jt(1,2))
			R(2,2)=0.0
			R(2,3)=haf*(Jt(2,3)-Jt(3,2))

			R(3,1)=haf*(Jt(3,1)-Jt(1,3))
			R(3,2)=haf*(Jt(3,2)-Jt(2,3))
			R(3,3)=0.0
			
			As(1)=sum(S**2.0)
			As(2)=sum(-R**2.0)

			P = -(S(1,1)+S(2,2)+S(3,3))
			Q(i,j,k) =haf*( P*P-As(1)-As(2))
			
		end do
	end do
end do
!-----------------------------------------------------------------------
open (unit=10,file="../output/qInvariant.dat",status="unknown")
write (10,*)'TITLE ="qInvariant"'
write (10,*)'Variables = "x","y","z","Q"'
write (10,'(3(A,I3),A)')'ZONE i=',nP(1),',j=',nP(2),',k=',nP(3),',DATAPACKING="POINT"'

do k=1,nP(3)
	do j=1,nP(2)
		do i=1,nP(1)
			write(10,'(3(1X,F10.6),1X,E22.15)')x(i),y(j),z(k),Q(i,j,k)
		end do
	end do
end do
close(10)

print *, "File qInvariant.dat created"
!-----------------------------------------------------------------------
end program qInvariant
