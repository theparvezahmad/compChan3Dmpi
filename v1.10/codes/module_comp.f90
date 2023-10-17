module module_comp
use module_const
use module_read
use module_MPI

implicit none

double precision,protected 									:: dXi,dx,dy,invdx,invdy,invdXi
double precision,dimension(nP(1)),public 	:: x1
double precision,dimension(nP(2)),public 	:: x2
double precision,dimension(nP(3)),public 	:: x3
double precision,dimension(nP(3)),protected 	:: jac,jac2,Xi,invJac
!=======================Variable Declaration Starts=====================
integer								:: m
double precision 			:: tauxzL,tauxzU
double precision 			:: dAyz,Ayz

!--------------------------Temporary Variables--------------------------
integer,allocatable,dimension(:,:) 	::tArr,bb,ee
!-----------------------------------------------------------------------
double precision,dimension(nP(3),7),protected :: stat1Spc	! In order of rho,u,v,w,e,t,p 
double precision,dimension(0:nP(1)+1,0:nP(2)+1,0:nP(3)+1,5,2),protected	:: solV
double precision,dimension(0:nP(1)+1,0:nP(2)+1,0:nP(3)+1,7),protected		:: var ! In order of rho,U,V,W,E,T,P
double precision,dimension(nP(1),nP(2),nP(3),5,2)	,protected						:: F,G,H,R
double precision,dimension(nP(1),nP(2),nP(3)),protected								:: mu,kt


contains

subroutine createMesh()
implicit none

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

if(debugOp==1) then
	OPEN (UNIT=10,FILE='Mesh.dat')

	write(10,*) 'TITLE ="Mesh"'
	write(10,*) 'Variables ="x","y","z"'
	write(10,'(3(A,I3),A)') 'ZONE k=',nP(3),',j=',nP(2),',i=',nP(1),',DATAPACKING="POINT"'

	do k=1,nP(3)
		write(10,*)
		do j=1,nP(2)
			write(10,*)
			do i=1,nP(1)
				write(10,'(3(F6.4))') x1(i),x2(j),x3(k)
			end do
		end do
	enddo
	close(10)	

endif

end subroutine createMesh

subroutine calcMeshParam()
implicit none

Xi(1) =0.0d0
dx = x1(nP(1))/dble(nP(1)-1)
dy = x2(nP(2))/dble(nP(2)-1)
dXi= x3(nP(3))/dble(nP(3)-1)

invdx = 1.0d0/dx
invdy = 1.0d0/dy
invdXi= 1.0d0/dXi

do k=1,nP(3)-1
	Xi(k+1) =Xi(k) + dXi
enddo

do k=1,nP(3)
	jac(k)=    0.50d0*Pi* dsin(0.50d0*Xi(k)*Pi)
	jac2(k)= 0.25d0*Pi*Pi*dcos(0.50d0*Xi(k)*Pi)
enddo
invJac(2:nP(3))	= 1.0d0/jac(2:nP(3))

end subroutine calcMeshParam

subroutine initSolV()
implicit none

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

end subroutine initSolV

subroutine calcFlux()
implicit none

double precision,private			::Ut,Vt,Wt,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,Pe
double precision,private			::dTdx,dTdy,dTdz,divg,muT,ktT,PressT,dissF,dissG,dissH
double precision,dimension(nP(1),nP(2),nP(3),5),private							:: fD,gD,hD

do k=bn(3),en(3)
	do j=bn(2),en(2)
		do i=bn(1),en(1)
			
			Pe = abs(solV(i,j,k,2,1))*dx/(mu(i,j,k)*invRe)

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
			Pe = abs(solV(i,j,k,3,1))*dy/(mu(i,j,k)*invRe)
			
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
			Pe = abs(solV(i,j,k,4,1))*dXi*Jac(k)/(mu(i,j,k)*invRe) 

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

call share1P()

do k=bn(3),en(3)
	do j=bn(2),en(2)
		do i=bn(1),en(1)
			F(i,j,k,:,2) = ( fD(i+1,j,k,:) - fD(i-1,j,k,:) )*haf*invdx

			G(i,j,k,:,2) = ( gD(i,j+1,k,:) - gD(i,j-1,k,:) )*haf*invdy

			H(i,j,k,:,2) = ( hD(i,j,k+1,:) - hD(i,j,k-1,:) )*haf*invdXi*invJac(k)
		enddo
	enddo
enddo

end subroutine calcFlux

subroutine calcSolV()
implicit none
double precision,dimension(nP(1),nP(2),nP(3),5),private							:: srcV
double precision,dimension(nP(1),nP(2),nP(3),5)	,private						:: R

do k=bn(3),en(3)
	do j=bn(2),en(2)
		do i=bn(1),en(1)
			srcV(i,j,k,1) 	= zer
			srcV(i,j,k,2) 	= haf*(tauxzL-tauxzU)*stat1Spc(k,1)
			srcV(i,j,k,3:4) = zer 
			srcV(i,j,k,5) 	= ggM*var(i,j,k,2)*haf*(tauxzL-tauxzU)*stat1Spc(k,1)
			
			R(i,j,k,:) = srcV(i,j,k,:) - (F(i,j,k,:,1) + F(i,j,k,:,2) + G(i,j,k,:,1) + G(i,j,k,:,2) + H(i,j,k,:,1) + H(i,j,k,:,2))

			if (iter.eq.1) then
				solV(i,j,k,:,2)	= solV(i,j,k,:,1)
				solV(i,j,k,:,1)	= solV(i,j,k,:,1) + haf*dt*R(i,j,k,:)
			endif

			if (iter.eq.2) then
				solV(i,j,k,:,1) = solV(i,j,k,:,2) + dt*( R(i,j,k,:) )
			endif
			
		enddo
	enddo
enddo

end subroutine calcSolV

subroutine checkSoln()
implicit none
character(len=30),private			:: filename
logical,private						:: isBlown=.false.

do k=bs(3),es(3)
	do j=bs(2),es(2)
		do i=bs(1),es(1)
		
			if (var(i,j,k,1) < 0.0d0 .or. var(i,j,k,6) < 0.0d0) then
				write(6,'(I10,4(I4),2(F8.4))') nstep,id,i,j,k,var(i,j,k,1),var(i,j,k,6)
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
				write(10,'(3(1X,F7.4),7(1X,E12.5))') x1(i),x2(j),x3(k),var(i,j,k,1),var(i,j,k,2),var(i,j,k,3),var(i,j,k,4)	&
																										,var(i,j,k,5),var(i,j,k,6),var(i,j,k,7)
			end do
		end do
	enddo
	close(10)	
	stop
endif
	
end subroutine checkSoln

subroutine calcPrim()
implicit none

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

end subroutine calcPrim

subroutine calcBC
implicit none

double precision,private			::dudxdz,dvdydz,d2wdz2,drhow2dz,lhs

do j=bs(2),es(2)
	do i=bs(1),es(1)
	
		if(bottom) then
			var(i,j,1,2:4)	= 0.0d0
			var(i,j,1,6)		= 1.0d0
			var(i,j,1,5)		= var(i,j,1,6) +	ggM*haf*sum(var(i,j,1,2:4)**two)
			var(i,j,1,1)		= gama*Mac**two*var(i,j,1,7) / var(i,j,1,6)
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

end subroutine calcBC

subroutine calcTauWall()
implicit none

double precision,private :: dwdx,dudz

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
			write(10,'(3(1X,F15.7))')time, tauxzL, tauxzU
			close(10) 
		endif
				
	end if
end subroutine calcTauWall

subroutine avgSpatial() !call it at step 10 & remove cond inside
implicit none
integer,private :: p,nXY
double precision,private			::loc(nP(1),nP(2),7)
double precision,dimension(nP(3),7),private :: stat2Spc ! In order of rho,u,v,w,e,t,p
double precision,dimension(nP(3),9),private :: corrSpc 	! In order of rhow,rhot,uv,uw,vw,vt,tw,pw,pt

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
end subroutine avgSpatial

end module module_comp
