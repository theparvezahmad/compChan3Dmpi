module module_const
use module_read, only:x3
implicit none
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

contains

subroutine calcCoeffs()
implicit none
double precision,private 			:: alpha

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
end subroutine setupCoeffs

end module module_const
