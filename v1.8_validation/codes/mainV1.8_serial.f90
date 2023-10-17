program chanCompr
!=======================================================================
use comData
implicit none
!=======================Variable Declaration Starts=====================
character(len=30)			:: filename
logical								:: isBlown=.false.
integer								:: i,j,k,p,nstep,iter,st,nXY
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
double precision			::dudxdz,dvdydz,d2wdz2,drhow2dz,rhs
!integer,allocatable,dimension(:,:) 	::tArr,bb,ee

!-----------------------------------------------------------------------
double precision,dimension(nP(3),7) :: stat1Spc	! In order of rho,u,v,w,e,t,p 
double precision,dimension(nP(3),7) :: stat2Spc ! In order of rho,u,v,w,e,t,p
double precision,dimension(nP(3),9) :: corrSpc 	! In order of rhow,rhot,uv,uw,vw,vt,tw,pw,pt

double precision,dimension(0:nP(1)+1,0:nP(2)+1,0:nP(3)+1,5,2)	:: solV
double precision,dimension(nP(1),nP(2),nP(3),7)								:: var ! In order of rho,U,V,W,E,T,P
double precision,dimension(nP(1),nP(2),nP(3),2:5)							:: fD,gD,hD
double precision,dimension(nP(1),nP(2),nP(3),5)								:: srcV
double precision,dimension(nP(1),nP(2),nP(3),5,2)							:: F,G,H,R
                                    
double precision,dimension(nP(1),nP(2),nP(3))									:: mu,kt

!-------------------------Debug Variables-------------------------------
double precision,dimension(nP(3))			:: pl
!integer,allocatable,dimension(:,:) 		::db(r)

!========================Variable Declaration Ends======================
F=zer; G=zer; H=zer; R=zer; srcV=zer; solV=zer
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

write(*,*) "=================================================================="
write(*,*) "Compressible Channel Flow Solver started at :",dateTime()
write(*,*) "=================================================================="

if(isDebug) then
	do k=1,nP(3)
		pl(k)=poly(x3(k))
	enddo
	
	tmpR1=(		f1d2c0*pl(1)	&
					+ f1d2c1*pl(2) 	&
					+ f1d2c2*pl(3)	&
				)	
	
	tmpR2=(		b1d2c0*pl(nP(3)  )	&
					+ b1d2c1*pl(nP(3)-1) 	&
					+ b1d2c2*pl(nP(3)-2)	&
				)
				
	tmpR3=(		f2d1c0*pl(1)	&
					+ f2d1c1*pl(2) 	&
					+ f2d1c2*pl(3)	&
				)
	
	tmpR4=(		b2d1c0*pl(nP(3)  )	&
					+ b2d1c1*pl(nP(3)-1) 	&
					+ b2d1c2*pl(nP(3)-2)	&
				)	
				
	write(6,*) "Bottom 1st Derivative: ",tmpR1,dpoly (x3(1))
	write(6,*) "Top    1st Derivative: ",tmpR2,dpoly (x3(nP(3)))
	write(6,*) "Bottom 2nd Derivative: ",tmpR3,d2poly(x3(1))
	write(6,*) "Top    2nd Derivative: ",tmpR4,d2poly(x3(nP(3)))

endif
!======================Read Combined Input File========================
open(unit=10,file=inSoln)

read (10,100)time,nstep
read (10,*)
read (10,*)

100	format(8X,F20.10,I9,1X)

do i=1,nP(1)
  read(10,*)
  do j=1,nP(2)
    read(10,*)
		do k=1,nP(3)
			read(10,101,iostat=st) x1(i),x2(j),x3(k),var(i,j,k,1),var(i,j,k,2),var(i,j,k,3),var(i,j,k,4)	&
																								 ,var(i,j,k,5),var(i,j,k,6),var(i,j,k,7)
																								 
			if(st/=0)	then
				write(6,115) "Error in reading input file at",i,j,k,st
				stop
			endif
				
		end do
	end do
end do

close(10)
115 format (A,4(1X,I4))
101	format(3(1X,F10.6),7(1X,E22.15))

write(6,*) "Input file successfully read"

!--------------------Initial Solution Vector----------------------------
do k=1,nP(3)
	do j=1,nP(2)
		do i=1,nP(1)
			solV(i,j,k,1,1)		=	var(i,j,k,1)
			solV(i,j,k,2:5,1)	=	var(i,j,k,1) * var(i,j,k,2:5)
			
			mu(i,j,k)					= var(i,j,k,6)**0.7d0 
			kt(i,j,k)					= var(i,j,k,6)**0.7d0
		enddo
	enddo
enddo

!========================Main time loop (starts)========================
DO WHILE (nstep < maxstep)

	if(mod(nstep,10)==0) call cpu_time(tStart)
	nstep = nstep+1
	time  = time + dt
	
	!========================RK2 Sub-steps (starts)=======================
	do iter=1,maxiter
	
	solV(0      ,:,:,:,1) = solV(nP(1)-2,:,:,:,1)
	solV(1      ,:,:,:,1) = solV(nP(1)-1,:,:,:,1)
	solV(nP(1)  ,:,:,:,1) = solV(2      ,:,:,:,1)
	solV(nP(1)+1,:,:,:,1) = solV(3      ,:,:,:,1)

	solV(:,0      ,:,:,1) = solV(:,nP(2)-2,:,:,1)
	solV(:,1      ,:,:,1) = solV(:,nP(2)-1,:,:,1)
	solV(:,nP(2)  ,:,:,1) = solV(:,2      ,:,:,1)
	solV(:,nP(2)+1,:,:,1) = solV(:,3      ,:,:,1)
		
		do k=2,nP(3)-1
			do j=2,nP(2)-1
				do i=2,nP(1)-1
					
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
						F(i,j,k,1:5,1) = var(i,j,k,2)*(      -solV (i+2,j,k,1:5,1) +	&
																						8.0d0*solV (i+1,j,k,1:5,1) -	&
																						8.0d0*solV (i-1,j,k,1:5,1) +	&
																					        solV (i-2,j,k,1:5,1) 		&
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
						G(i,j,k,1:5,1) = var(i,j,k,3)*(      -solV (i,j+2,k,1:5,1) +	&
																						8.0d0*solV (i,j+1,k,1:5,1) -	&
																						8.0d0*solV (i,j-1,k,1:5,1) +	&
																						      solV (i,j-2,k,1:5,1) 		&
																					)*invdy*one12th
					endif  
					!-------------------------------------------------------------
					Pe = Re*abs(solV(i,j,k,4,1))*dXi*invJac(k)/mu(i,j,k) 

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
						H(i,j,k,1:5,1) = var(i,j,k,4)*(      -solV (i,j,k+2,1:5,1) +	&
																						8.0d0*solV (i,j,k+1,1:5,1) -	&
																						8.0d0*solV (i,j,k-1,1:5,1) +	&
																						      solV (i,j,k-2,1:5,1) 		&
																					)*invdXi*invJac(k)*one12th
					endif

				enddo
			enddo
		enddo
		!-------------------------------------------------------------------
		!===========Discretization of Diffusive terms in Flux Form==========
		do k=1,nP(3)
			do j=2,nP(2)-1
				do i=2,nP(1)-1
					
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
					dissH	=ggM*invRe*( -muT*Ut*(dwdx + dudz) - muT*Vt*(dwdy + dvdz) + two3rd*muT*Wt*divg - two*muT*Wt*dwdz  )
					
					!-------------------------------------------------------------
					fD(i,j,k,2)	=	PressT + two3rd*muT*invRe*divg - two*muT*invRe*dudx 

					fD(i,j,k,3)	=	-muT*invRe*( dvdx + dudy )

					fD(i,j,k,4)	=	-muT*invRe*( dwdx + dudz )
					
					fD(i,j,k,5)	=	-gama*ktT*invRe*invPr*dTdx + ggM*Ut*PressT + dissF

					!-------------------------------------------------------------
					gD(i,j,k,2)	=	-muT*invRe*(dvdx + dudy) 

					gD(i,j,k,3)	=	PressT + two3rd*muT*invRe*divg - two*muT*invRe*dvdy

					gD(i,j,k,4)	=	-muT*invRe*(dwdy + dvdz)

					gD(i,j,k,5)	=	-gama*ktT*invRe*invPr*dTdy + ggM*Vt*PressT + dissG

					!-------------------------------------------------------------
					hD(i,j,k,2)	=	-muT*invRe*(dwdx + dudz)

					hD(i,j,k,3)	=	-muT*invRe*(dwdy + dvdz)

					hD(i,j,k,4)	=	PressT + two3rd*muT*invRe*divg - two*muT*invRe*dwdz 

					hD(i,j,k,5)	=	-gama*ktT*invRe*invPr*dTdz + ggM*Wt*PressT + dissH
					!-------------------------------------------------------------
				enddo 
			enddo
		enddo
		
		fD(1      ,:,:,:) = fD(nP(1)-1,:,:,:)
		fD(nP(1)  ,:,:,:) = fD(2      ,:,:,:)

		gD(:,1      ,:,:) = gD(:,nP(2)-1,:,:)
		gD(:,nP(2)  ,:,:) = gD(:,2      ,:,:)
		
		do k=2,nP(3)-1
			do j=2,nP(2)-1
				do i=2,nP(1)-1
					F(i,j,k,2:5,2) = ( fD(i+1,j,k,2:5) - fD(i-1,j,k,2:5) )*haf*invdx

					G(i,j,k,2:5,2) = ( gD(i,j+1,k,2:5) - gD(i,j-1,k,2:5) )*haf*invdy

					H(i,j,k,2:5,2) = ( hD(i,j,k+1,2:5) - hD(i,j,k-1,2:5) )*haf*invdXi*invJac(k)
				enddo
			enddo
		enddo
		!---------------------------------------------------------------
		do k=2,nP(3)-1
			do j=2,nP(2)-1
				do i=2,nP(1)-1
					srcV(i,j,k,2) = haf*(tauxzL-tauxzU)*stat1Spc(k,1)
					srcV(i,j,k,5) = ggM*var(i,j,k,2)*haf*(tauxzL-tauxzU)*stat1Spc(k,1)
					
					R(i,j,k,:,1) = srcV(i,j,k,:) - ( sum(F(i,j,k,:,:),2) + sum(G(i,j,k,:,:),2) + sum(H(i,j,k,:,:),2) )

					if (iter.eq.1) then
						R(i,j,k,:,2)		= R(i,j,k,:,1)
						solV(i,j,k,:,2)	= solV(i,j,k,:,1)

						solV(i,j,k,:,1)	= solV(i,j,k,:,1) + dt*R(i,j,k,:,1)
					endif

					if (iter.eq.2) then
						solV(i,j,k,:,1) = solV(i,j,k,:,2) + haf*dt*( R(i,j,k,:,2) + R(i,j,k,:,1) )
					endif
					
				enddo
			enddo
		enddo
		
		! Finding Primitive Variables
		do k=2,nP(3)-1
			do j=2,nP(2)-1
				do i=2,nP(1)-1
					var(i,j,k,1  )	= solV(i,j,k,1,1)
					var(i,j,k,2:5)	= solV(i,j,k,2:5,1) / solV(i,j,k,1,1)
					var(i,j,k,6)		= var(i,j,k,5) - ggM*haf*sum(var(i,j,k,2:4)**two)
					var(i,j,k,7)		= var(i,j,k,1)*var(i,j,k,6) / (gama*Mac**two)
				enddo
			enddo
		enddo
				
		! Boundary Conditions
		var(1      ,:,:,:) = var(nP(1)-1,:,:,:)
		var(nP(1)  ,:,:,:) = var(2      ,:,:,:)

		var(:,1      ,:,:) = var(:,nP(2)-1,:,:)
		var(:,nP(2)  ,:,:) = var(:,2      ,:,:)
		
		
		do j=2,nP(2)-1
			do i=2,nP(1)-1

				var(i,j,1,2:4) = 0.0d0
				var(i,j,1,6) = 1.0d0
				var(i,j,1,5) = var(i,j,1,6) +	ggM*haf*sum(var(i,j,1,2:4)**two)
				var(i,j,1,1) = gama*Mac**two*var(i,j,1,7) / var(i,j,1,6) 
				!-------------------------Pressure BC-------------------------
				!Simplifying z-momentum equation assuming impervious wall		
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
				
				rhs			= invRe*(-dudxdz -dvdydz + two3rd*(dudxdz + dvdydz - two*d2wdz2)) + drhow2dz
				
				var(i,j,1,7) = -(rhs + f1d2c1*var(i,j,2,7) + f1d2c2*var(i,j,3,7) ) / f1d2c0
				!-------------------------------------------------------------

				var(i,j,nP(3),2:4)= 0.0d0
				var(i,j,nP(3),6) 	= 1.0d0
				var(i,j,nP(3),5) 	= var(i,j,nP(3),6) +	ggM*haf*sum(var(i,j,nP(3),2:4)**2)
				var(i,j,nP(3),1) 	= gama*Mac**two*var(i,j,nP(3),7) / var(i,j,nP(3),6) 
				!-------------------------Pressure BC-------------------------
				!Simplifying z-momentum equation assuming impervious wall				
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
				
				rhs			= invRe*(-dudxdz -dvdydz + two3rd*(dudxdz + dvdydz - two*d2wdz2)) + drhow2dz
				
				var(i,j,nP(3),7) = -(rhs + b1d2c1*var(i,j,nP(3)-1,7) + b1d2c2*var(i,j,nP(3)-2,7) ) / b1d2c0
				!---------------------------------------------------------------
			
			enddo
		enddo
		
		! Updating Thermophysical Properties
		mu(:,:,:) = var(:,:,:,6)**0.7d0
		kt(:,:,:) = var(:,:,:,6)**0.7d0
		
	enddo
	!=========================RK2 Sub-steps (ends)========================
	!=====================Check for unphysical values=====================
	do k=1,nP(3)
		do j=1,nP(2)
			do i=1,nP(1)
			
				if (var(i,j,k,1) < 0.0d0 .or. var(i,j,k,6) < 0.0d0) then
					write(6,102) nstep,i,j,k,var(i,j,k,1),var(i,j,k,6)
					isBlown=.true.
				endif
				
			enddo
		enddo
	enddo
	
	if(isBlown) then
	
		OPEN (UNIT=10,FILE="../debug/3D.dat")

		write(10,103) 'TITLE ="',time,nstep,'"'
		write(10,*) 'Variables ="x","y","z","Rho","U","V","W","E","T","P"'
		write(10,104) 'ZONE k=',nP(3),',j=',nP(2),',i=',nP(1),',DATAPACKING="POINT"'

		do k=1,nP(3)
			write(10,*)
			do j=1,nP(2)
				write(10,*)
				do i=1,nP(1)
					write(10,116) x1(i),x2(j),x3(k),var(i,j,k,1),var(i,j,k,2),var(i,j,k,3),var(i,j,k,4)	&
																											,var(i,j,k,5),var(i,j,k,6),var(i,j,k,7)
				end do
			end do
		enddo
		close(10)	
		stop
	endif
	
	102 format(I10,4(I4),2(F8.4)) 
	116 format(3(1X,F6.4),7(1X,E10.5))
	
	!==================Write output files (starts)========================
	if(mod(nstep,1000).eq.0) then
		
		OPEN (UNIT=10,FILE='../debug/3D.dat')

		write(10,103) 'TITLE ="',time,nstep,'"'
		write(10,*) 'Variables ="x","y","z","Rho","U","V","W","E","T","P"'
		write(10,104) 'ZONE k=',nP(3),',j=',nP(2),',i=',nP(1),',DATAPACKING="POINT"'

		do k=1,nP(3)
			write(10,*)
			do j=1,nP(2)
				write(10,*)
				do i=1,nP(1)
					write(10,101) x1(i),x2(j),x3(k),var(i,j,k,1),var(i,j,k,2),var(i,j,k,3),var(i,j,k,4)	&
																										,var(i,j,k,5),var(i,j,k,6),var(i,j,k,7)
				end do
			end do
		enddo
		close(10)	
		103	format(A,F20.10,I9,A)
		104	format(A,I3,A,I3,A,I3,A)

	endif
	!==================Write output files (ends)==========================
	
	!==================Write Diagnostic Files (starts)====================
	
	if (mod(nstep,10).eq.0) then
		!------------------------Write Kinetic Energy-----------------------
		tmpR1=0.0d0
		do k=1,nP(3)
			do j=1,nP(2)
				do i=1,nP(1)
					tmpR1 = tmpR1 + sum(var(i,j,k,2:4)**two)  
				enddo
			enddo
		enddo
			tmpR1 = tmpR1/product(nP)
			open(10,file='../output/KE.dat',access='append')
			write(10,108) time,tmpR1
			close(10)
		108	format(2(1X,F15.7))
		!------------------------Write Bulk Density-------------------------
		Ayz = (x3(nP(3))-x3(2))*(x2(nP(2))-x2(2))
		
		tmpA1	=0.0d0
		do i=1,5
			
				do k=2,nP(3)-1
					do j=2,nP(2)-1
						dAyz = dy*(x3(k+1)-x3(k))
						tmpA1(i) =tmpA1(i) + quar*(	var(plane(i),j  ,k  ,1) +	&
																	var(plane(i),j+1,k  ,1) +	&
																	var(plane(i),j  ,k+1,1) +	&
																	var(plane(i),j+1,k+1,1)		&
																)*dAyz
					enddo
				enddo
		enddo  

			tmpA1=tmpA1/Ayz
			open (unit=10,file='../output/rhobulk.dat',access='append')
			write(10,109)time,tmpA1(1),tmpA1(2),tmpA1(3),tmpA1(4),tmpA1(5)
			close(10)
		109	format(6(1X,F15.7))
		!------------------------Write Bulk Velocity------------------------ 
		tmpA1	=0.0d0
		do i=1,5
			
				do k=2,nP(3)-1
					do j=2,nP(2)-1
						dAyz = dy*(x3(k+1)-x3(k))
						tmpA1(i) =tmpA1(i) + quar*(	var(plane(i),j  ,k  ,2) +	&
																	var(plane(i),j+1,k  ,2) +	&
																	var(plane(i),j  ,k+1,2) +	&
																	var(plane(i),j+1,k+1,2)		&
																)*dAyz
					enddo
				enddo
		enddo
		
			tmpA1=tmpA1/Ayz
			open (unit=10,file='../output/Ubulk.dat',access='append')
			write(10,109)time,tmpA1(1),tmpA1(2),tmpA1(3),tmpA1(4),tmpA1(5)
			close(10)
		
		!--------------------------Write Massflow--------------------------- 
		tmpA1	=0.0d0
		do i=1,5
				do k=2,nP(3)-1
					do j=2,nP(2)-1
						dAyz = dy*(x3(k+1)-x3(k))
						tmpA1(i) =tmpA1(i) + quar*(	var(plane(i),j  ,k  ,1)*var(plane(i),j  ,k  ,2) +	&
																	var(plane(i),j+1,k  ,1)*var(plane(i),j+1,k  ,2) + &
																	var(plane(i),j  ,k+1,1)*var(plane(i),j  ,k+1,2) +	&
																	var(plane(i),j+1,k+1,1)*var(plane(i),j+1,k+1,2)		&
																)*dAyz
					enddo  
				enddo  
		enddo  

			tmpA1=tmpA1/Ayz
			open (unit=10,file='../output/massflow.dat',access='append')
			write(10,109)time,tmpA1(1),tmpA1(2),tmpA1(3),tmpA1(4),tmpA1(5)
			close(10)
		
	endif	
	!====================Write Diagnostic Files (ends)====================
	
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<Check Residual>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	!==================Wall shear stress (starts)=========================
	if(mod(nstep,istat)==0) then
		!-------------------------------------------------------------------
		tauxzL=0.0d0
		
			do j=2,nP(2)-1
				do i=2,nP(1)-1
					dwdx = mu(i,j,1) * (var(i+1,j,1,4) - var(i-1,j,1,4)) *haf*invdx !zero for impervious bottom wall

					dudz = mu(i,j,1)*	(		f1d2c0*var(i,j,1,2)	&
															+ f1d2c1*var(i,j,2,2) &
															+ f1d2c2*var(i,j,3,2)	&
														)

					tauxzL = tauxzL + (dwdx + dudz)/Re
				enddo
			enddo
		
		tauxzL= tauxzL/dble((nP(1)-2)*(nP(2)-2))
		!-------------------------------------------------------------------
		
		tauxzU=0.0d0

			do j=2,nP(2)-1
				do i=2,nP(1)-1
					dwdx = mu(i,j,nP(3)) * (var(i+1,j,nP(3),4) - var(i-1,j,nP(3),4)) *haf*invdx !zero for impervious top wall
					
					dudz = mu(i,j,nP(3)) * (		b1d2c0*var(i,j,nP(3)  ,2)	&
																		+ b1d2c1*var(i,j,nP(3)-1,2) &
																		+	b1d2c2*var(i,j,nP(3)-2,2)	&
																	)

					tauxzU = tauxzU + (dwdx + dudz)/Re
				enddo
			enddo
		
		tauxzU= tauxzU/dble(((nP(1)-2)*(nP(2)-2)))
		
		!-------------------------------------------------------------------
			open(10,file='../output/tauWall.dat',access='append')
			write(10,110)time, tauxzL, tauxzU
			close(10)
		110	format(3(1X,F15.7))
				
	end if
	!==================Wall shear stress (ends)===========================
	
	!==================Velocity Correlations (starts)=====================
	! To be included later
	!==================Velocity Correlations (ends)=======================
	
	!==================Spatial Statistics (starts)========================
	
	if(mod(nstep,istat)==0) then
		nXY=(nP(1)-1)*(nP(2)-1)
		
		do k=1,nP(3)
			loc=0.0d0
				do j=2,nP(2)-1
					do i=2,nP(1)-1
						loc(i,j,1:7)	= quar*(var(i,j,k,1:7) + var(i+1,j,k,1:7) + var(i,j+1,k,1:7) + var(i+1,j+1,k,1:7) )
					enddo
				enddo		
			
			do p=1,7
				stat1Spc(k,p) = sum(loc(:,:,p)   )/nXY
				stat2Spc(k,p) = sum(loc(:,:,p)**2)/nXY
			enddo

			corrSpc(k,1) = sum(loc(:,:,1)*loc(:,:,4))/nXY
			corrSpc(k,2) = sum(loc(:,:,1)*loc(:,:,6))/nXY
			corrSpc(k,3) = sum(loc(:,:,2)*loc(:,:,3))/nXY
			corrSpc(k,4) = sum(loc(:,:,2)*loc(:,:,4))/nXY
			corrSpc(k,5) = sum(loc(:,:,3)*loc(:,:,4))/nXY
			corrSpc(k,6) = sum(loc(:,:,3)*loc(:,:,6))/nXY
			corrSpc(k,7) = sum(loc(:,:,6)*loc(:,:,4))/nXY
			corrSpc(k,8) = sum(loc(:,:,7)*loc(:,:,4))/nXY
			corrSpc(k,9) = sum(loc(:,:,7)*loc(:,:,6))/nXY
		enddo
		
			!-----------------------------------------------------------------
			open(unit=10,file='../output/moment1st.dat',access='append')
			write(10,111,advance='no') time
			do k=1,nP(3)
				write(10,112,advance='no') stat1Spc(k,1),stat1Spc(k,2),stat1Spc(k,3),stat1Spc(k,4)	&
																								,stat1Spc(k,5),stat1Spc(k,6),stat1Spc(k,7)
			enddo
			write (10, *)
			close(10)
			111	format(1(1X,F15.7))
			112	format(7(1X,F15.7))
			!-----------------------------------------------------------------
			open(unit=10,file='../output/moment2nd.dat',access='append')
			write(10,111,advance='no') time
			do k=1,nP(3)
				write(10,112,advance='no') stat2Spc(k,1),stat2Spc(k,2),stat2Spc(k,3),stat2Spc(k,4)	&
																								,stat2Spc(k,5),stat2Spc(k,6),stat2Spc(k,7)
			enddo
			write (10, *)
			close(10)
			!-----------------------------------------------------------------
			open(unit=10,file='../output/momentJoint.dat',access='append')
			write(10,104,advance='no')time
			do k=1,nP(3)
				write(10,113,advance='no') corrSpc(k,1),corrSpc(k,2),corrSpc(k,3)	&
																	,corrSpc(k,4),corrSpc(k,5),corrSpc(k,6)	&
																	,corrSpc(k,7),corrSpc(k,8),corrSpc(k,9)
			enddo
			write (10,*)
			close(10)
			113	format(9(1X,F15.7))
			!-----------------------------------------------------------------
		
	endif
	!====================Spatial Statistics (ends)========================
	
	if(mod(nstep,10)==0) then
		call cpu_time(tEnd)
		write(*,114) nstep,time,' | CPU Time:',tEnd-tStart,'sec for 10 iters | ', dateTime()
		write(*,*) '--------------------------------------------------------'
	endif
	114  format(I10,F10.5,A,F6.3,2(1X,A))
	
END DO
!========================Main time loop (ends)==========================

write(*,*) "Compressible Channel Flow Solver exited at :",dateTime()
!=======================================================================
contains

function poly(x) result(f)
double precision,intent(in) :: x
double precision 						:: f
f		= x**3.0d0 - 5.0d0*x**2.0d0 + 3.0d0*x - 10.0d0
end function poly

function dpoly(x) result(f)
double precision,intent(in) :: x
double precision 						:: f
f	= 3.0d0*x**2.0d0 - 10.0d0*x + 3.0d0
end function dpoly

function d2poly(x) result(f)
double precision,intent(in) :: x
double precision 						:: f
f		= 6.0d0*x - 10.0d0
end function d2poly

end program chanCompr
