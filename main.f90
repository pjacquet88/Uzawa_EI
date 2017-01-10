PROGRAM main
  USE Solver
  USE Tools
  IMPLICIT NONE
  REAL*8,DIMENSION(:),ALLOCATABLE:: U1,U2,F1,F2,P,P0,P_int
  INTEGER :: i,j,time
  REAL*8 :: xi,yi,t1,t2,normDiff,normExact,x,y

  CALL CPU_TIME(t1)
  CALL read_param("param.txt")
  ALLOCATE(U1((nx-1)*(ny-1)),U2((nx-1)*(ny-1)),F1((nx-1)*(ny-1)),F2((nx-1)*(ny-1)),P(nx*ny),P0(nx*ny))

  U1 = 0
  U2 = 0
  P0 = 0
  P_int=P0

  DO time=1,10
     F1=0
     F2=0
     DO i=1,ny-1
        DO j=1,nx-1
           F1((i-1)*(nx-1)+j) = fx(j*dx,i*dy,(time-1)*dt,testcase)
           F2((i-1)*(nx-1)+j) = fy(j*dx,i*dy,(time-1)*dt,testcase)
        END DO
        F1((i-1)*(nx-1)+1) = F1((i-1)*(nx-1)+1) - coeffX*gx(0.d1,i*dy,(time-1)*dt,testcase)
        F2((i-1)*(nx-1)+1) = F2((i-1)*(nx-1)+1) - coeffX*gy(0.d1,i*dy,(time-1)*dt,testcase)
        F1(i*(nx-1)) = F1(i*(nx-1)) - coeffX*gx(lx,i*dy,(time-1)*dt,testcase)
        F2(i*(nx-1)) = F2(i*(nx-1)) - coeffX*gy(lx,i*dy,(time-1)*dt,testcase)
     END DO
     DO j=1,nx-1
        F1(j) = F1(j) - coeffY*hx(j*dx,0.d1,(time-1)*dt,testcase)
        F2(j) = F2(j) - coeffY*hy(j*dx,0.d1,(time-1)*dt,testcase)
        F1((ny-2)*(nx-1) + j) = F1((ny-2)*(nx-1) + j) - coeffY*hx(j*dx,ly,(time-1)*dt,testcase)
        F2((ny-2)*(nx-1) + j) = F2((ny-2)*(nx-1) + j) - coeffY*hy(j*dx,ly,(time-1)*dt,testcase)
     END DO

     F1 = F1 + U1/dt
     F2 = F2 + U2/dt

     CALL Uzawa(P_int,F1,F2,U1,U2,P,(time-1)*dt)
     P_int=P
  END DO
  !
  ! OPEN(unit=12,file='sol.dat',form='formatted',status='unknown',action='write')
  ! OPEN(unit=13,file='solExact.dat',form='formatted',status='unknown',action='write')
  !
  ! normDiff = 0.
  ! normExact = 0.
  ! DO i=1,ny-1
  !    DO j=1,nx-1
  !       xi = dx*j
  !       yi = dy*i
  !       WRITE(12,*) xi,yi,U((i-1)*(nx-1)+j)
  !       IF (p==2) THEN
  !          WRITE(13,*) xi,yi,sin(xi)+cos(yi)
  !          normDiff = max(abs(U((i-1)*(nx-1)+j)-(sin(xi)+cos(yi))),normDiff)
  !          normExact =  max(abs(sin(xi)+cos(yi)),normExact)
  !       ELSE IF (p==1) THEN
  !          WRITE(13,*) xi,yi,xi*(1-xi)*yi*(1-yi)
  !          normDiff = max(abs(U((i-1)*(nx-1)+j)-(xi*(1-xi)*yi*(1-yi))),normDiff)
  !          normExact = max(xi*(1-xi)*yi*(1-yi),normExact)
  !       END IF
  !    END DO
  !    WRITE(12,*)
  !    WRITE(13,*)
  ! END DO
  ! IF ((p==1).OR.(p==2)) THEN
  !    WRITE (*,*) "Erreur relative en norme infinie : ",normDiff/normExact
  ! END IF

  !!stop
  CLOSE(12)
  CLOSE(13)
  CLOSE(14)

open(unit=4,file='pression_e.dat',action='write')
open(unit=3,file='pression.dat',action='write')
do j=1,ny
  do i=1,nx
    x=(i-0.5)*dx
    y=(j-0.5)*dy
write(4,*),x,y,x+y-1
write(3,*),x,y,P(bij(i,j,nx))
  end do
  write(3,*),' '
  write(4,*),' '
end do
close(3)
close(4)

open(unit=4,file='vitesse_x_e.dat',action='write')
open(unit=3,file='vitesse_x.dat',action='write')
do j=1,ny-1
  do i=1,nx-1
    x=(i)*dx
    y=(j)*dy
write(4,*),x,y,x
write(3,*),x,y,U1(bij(i,j,nx-1))
  end do
  write(3,*),' '
  write(4,*),' '
end do
close(3)
close(4)


open(unit=4,file='vitesse_y_e.dat',action='write')
open(unit=3,file='vitesse_y.dat',action='write')
do j=1,ny-1
  do i=1,nx-1
    x=(i)*dx
    y=(j)*dy
write(4,*),x,y,-y
write(3,*),x,y,U2(bij(i,j,nx-1))
  end do
  write(3,*),' '
  write(4,*),' '
end do
close(3)
close(4)



  DEALLOCATE(U1,U2,F1,F2,P,P0)

  CALL CPU_TIME(t2)
  WRITE(*,*) "Temps total (s) : ",t2-t1

END PROGRAM main
