program avgTauWall
implicit none
integer::st,k,cnt
double precision,allocatable,dimension(:):: time,tauxzL,tauxzU
!-----------------------------------------------------------------------
open(unit=10,file="../output/tauWall.dat",status="old")
cnt=0
do
	read(10,*,iostat=st)
	if(st/=0) exit
	cnt=cnt+1
enddo
cnt=cnt-1
close(10)

write(6,*) "No of lines read :",cnt
!-----------------------------------------------------------------------
allocate(time(cnt))
allocate(tauxzL(cnt))
allocate(tauxzU(cnt))
!-----------------------------------------------------------------------
open(unit=10,file="../output/tauWall.dat",status="old")
do k=1,cnt
	read(10,'(3(1X,F15.7))',iostat=st) time(k), tauxzL(k), tauxzU(k)
enddo
if(st/=0) stop
close(10)

write(6,*) "File tauWall.dat successfully read"
!-----------------------------------------------------------------------
open(unit=10,file="../output/tauWallAvg.dat",status="unknown")
do k=1,cnt
	write(10,'(4(1X,F15.7))') time(k), tauxzL(k), tauxzU(k), (tauxzL(k)-tauxzU(k))/2.0
enddo
close(10)

write(6,*) "File tauWallAvg.dat created"

end program avgTauWall
