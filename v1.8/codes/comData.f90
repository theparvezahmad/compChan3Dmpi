module comData
implicit none
!=======================================================================
integer,parameter,dimension(3)					::nP			=[130,130,160]

INTEGER, PARAMETER											::maxiter	=2,&
																					maxstep	=int(5E7),&
																					istat		=10,&
																					lowWall	=1,&		!1)No Slip 2)Free Shear
																					inputMode =2		!1)Combined Input 2)Component Input
																					
DOUBLE PRECISION, PARAMETER							::PexVal	=2.0d0,&
																					PeyVal	=2.0d0,&
																					PezVal	=0.2d0,&
																					dt			=1.0d-4,&
																					Re			=3000.0d0,&
																					invRe		=1.0d0/3000.0d0,&
																					invPr		=1.0d0/0.71d0,&
																					gama		=1.4d0,&
																					Mac			=1.5d0,&	
																					one			=1.0d0,&
																					two			=2.0d0,&
																					zer			=0.0d0,&
																					quar		=0.25d0,&
																					haf			=0.5d0,&
																					two3rd	=2.0d0/3.0d0,&
																					one12th	=1.0d0/12.0d0,&
																					ggM			=gama*(gama-1)*Mac**2.0d0
																					
CHARACTER(LEN=30), PARAMETER						::inSoln	="../../Exp1/3Dcmprs.tp",&
																					outSoln	="../3Dm15uc.tp"
																					
LOGICAL,PARAMETER												::isDebug	=.false.,	&	
																					isNow		=.false.
!-----------------------------------------------------------------------
double precision 	:: dXi,dx,dy,Pi,x1(nP(1)),x2(nP(2)),invdx,invdy,invdXi
double precision,dimension(nP(3)) :: x3,jac,jac2,Xi,invJac

contains
subroutine setupMesh()
implicit none
integer					:: i,j,k

Pi=4.d0*datan(1.0d0)

x1(1) = 0.0d0
x2(1) = 0.0d0
Xi(1) =0.0d0
x3(1) =0.0d0

x1(nP(1)) = 4.0d0*Pi
x2(nP(2)) = 4.0d0*Pi/3.0d0
x3(nP(3)) = 2.0d0

dx = x1(nP(1))/dble(nP(1)-1)
dy = x2(nP(2))/dble(nP(2)-1)
dXi= x3(nP(3))/dble(nP(3)-1)

invdx = 1.0d0/dx
invdy = 1.0d0/dy
invdXi= 1.0d0/dXi

do i=1,nP(1)-1
	x1(i+1) =x1(i) + dx    
enddo

do j=1,nP(2)-1
	x2(j+1) =x2(j) + dy    
enddo

do k=1,nP(3)-1
	Xi(k+1) =Xi(k) + dXi
enddo

do k=1,nP(3)
	x3(k)=         1.0d0-dCos(0.50d0*Xi(k)*Pi)
	jac(k)=    0.50d0*Pi* dSin(0.50d0*Xi(k)*Pi)
	jac2(k)= 0.25d0*Pi*Pi*dCos(0.50d0*Xi(k)*Pi)
enddo
invJac(2:nP(3))	= 1.0d0/jac(2:nP(3))

if(isDebug) then
	OPEN (UNIT=10,FILE='Mesh.dat')

		write(10,*) 'TITLE ="Mesh"'
		write(10,*) 'Variables ="x","y","z"'
		write(10,104) 'ZONE k=',nP(3),',j=',nP(2),',i=',nP(1),',DATAPACKING="POINT"'

		do k=1,nP(3)
			write(10,*)
			do j=1,nP(2)
				write(10,*)
				do i=1,nP(1)
					write(10,117) x1(i),x2(j),x3(k)
				end do
			end do
		enddo
		close(10)	

endif
104	format(3(A,I3),A)
117 format(3(F6.4))

end subroutine setupMesh

function dateTime()

implicit none
character(len=30)::dateTime
character(len=3):: ampm
integer:: d,h,m,n,s,y,mm,values(8)
character(len=3), parameter, dimension(12) :: &
month=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

call date_and_time(values=values)

y=values(1)
m=values(2)
d=values(3)
h=values(5)
n=values(6)
s=values(7)
mm=values(8)

if(h<12) then
	ampm='AM'
elseif(h==12) then
	if(n==0 .and. s==0) then
		ampm='Noon'
	else
		ampm='PM'
	endif
else
	h=h-12
	if(h<12) then
		ampm='PM'
	elseif(h==12) then
		if(n==0 .and. s==0) then
			ampm='Midnight'
		else
			ampm='AM'
		endif
	endif
endif

write(dateTime,'(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)')&
d,trim(month(m)),y,h,':',n,':',s,'.',mm,trim(ampm)
end function dateTime

end module comData
