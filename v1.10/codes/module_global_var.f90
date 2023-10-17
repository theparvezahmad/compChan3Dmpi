module module_global_var

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

integer,parameter,dimension(3)					::nP			=[130,130,160]

logical,parameter,dimension(3)				::period	=[.true.,.true.,.false.]

INTEGER, PARAMETER											::maxiter	=2,&
																					maxstep	=int(5E7)
																					
DOUBLE PRECISION, PARAMETER		::gama		=1.4d0,&
																					Pi=4.d0*datan(1.0d0),&
																					invPr		=1.0d0/0.71d0,&
																					one			=1.0d0,&
																					two			=2.0d0,&
																					zer			=0.0d0,&
																					quar		=0.25d0,&
																					haf			=0.5d0,&
																					two3rd	=2.0d0/3.0d0,&
																					one12th	=1.0d0/12.0d0
																					
integer,public :: i,j,k,st,tmpI1,tmpI2,tmpI3,tmpI4,tmpI5
double precision,public	:: tmpR1,tmpR2,tmpR3,tmpR4,tmpR5

double precision,protected 			:: dz1,dzN
double precision,protected 			:: f1d2c0, f1d2c1, f1d2c2, b1d2c0, b1d2c1, b1d2c2
double precision,protected 			:: f2d1c0, f2d1c1, f2d1c2, b2d1c0, b2d1c1, b2d1c2	

integer,dimension(3),protected		:: bs,es,bn,en
logical,protected 	:: bottom, top, master


integer,protected::resOp,topo(3), istat, debugOp
double precision,protected:: Re,invRe,Mac,dt,PexVal,PeyVal,PezVal,ggM

character(len=:),allocatable,protected::descp,inputPath
