PROGRAM manipulate
IMPLICIT NONE
Integer :: n,ierr,offset=500
Integer,PARAMETER :: nhmx=50000 ! The max number of thermo points
INTEGER :: nh
Real(8) :: tstart,tstop,th(nhmx),t9h(nhmx),rhoh(nhmx)
Real(8) :: yeh(nhmx),tdelstart       !NSE
Character (LEN=80) :: thermo_desc
 Open(10,file='th_nuwind_mb_s200t010')
 Read(10,"(a)") thermo_desc
 Read(10,*) tstart
 Read(10,*) tstop
 Read(10,*) tdelstart
 yeh = 0.0 
 Do n=1,nhmx
   Read(10,*,IOSTAT=ierr) th(n),t9h(n),rhoh(n),yeh(n)
   If(ierr==-1) Then
     Write(*,*) 'End of Thermo File Reached after',n,' records'
     nh=n
     Exit
   EndIf
 EndDo
CLOSE(10)
yeh=0.4
OPEN(20,file='th_nuwind_new',FORM='FORMATTED',STATUS='NEW')
 WRITE(20,"(a)") thermo_desc
 WRITE(20,*) tstart
 WRITE(20,*) tstop
 WRITE(20,*) tdelstart
 Do n=1,nh-1-offset
   WRITE(20,*) th(n),t9h(n+offset),rhoh(n+offset),yeh(n+offset)
 EndDo
CLOSE(20)
END PROGRAM