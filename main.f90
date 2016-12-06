PROGRAM main
  USE Functions
  USE Solver
  USE Tools
  IMPLICIT NONE
  REAL*8,DIMENSION(:),ALLOCATABLE:: U1,U2,F1,F2,P,P0
  INTEGER :: i,j,time
  REAL*8 :: xi,yi,t1,t2,normDiff,normExact

  CALL CPU_TIME(t1)
  CALL read_param("param.txt")
  ALLOCATE(U1((nx-1)*(ny-1)),U2((nx-1)*(ny-1)),F1((nx-1)*(ny-1)),F2((nx-1)*(ny-1)),P(nx*ny),P0(nx*ny))

  U1 = 0
  U2 = 0
  P0 = 0

  DO time=1,10
     F1=0
     F2=0
     DO i=1,ny-1
        DO j=1,nx-1
           F1((i-1)*(nx-1)+j) = f(j*dx,i*dy,(time-1)*dt,testcase)
           F2((i-1)*(nx-1)+j) = f(j*dx,i*dy,(time-1)*dt,testcase)
        END DO
        F1((i-1)*(nx-1)+1) = F1((i-1)*(nx-1)+1) - coeffX*g(0.d1,i*dy,(time-1)*dt,testcase)
        F2((i-1)*(nx-1)+1) = F2((i-1)*(nx-1)+1) - coeffX*g(0.d1,i*dy,(time-1)*dt,testcase)
        F1(i*(nx-1)) = F1(i*(nx-1)) - coeffX*g(lx,i*dy,(time-1)*dt,testcase)
        F2(i*(nx-1)) = F2(i*(nx-1)) - coeffX*g(lx,i*dy,(time-1)*dt,testcase)
     END DO
     DO j=1,nx-1
        F1(j) = F1(j) - coeffY*g(j*dx,0.d1,(time-1)*dt,testcase)
        F2(j) = F2(j) - coeffY*g(j*dx,0.d1,(time-1)*dt,testcase)
        F1((ny-2)*(nx-1) + j) = F1((ny-2)*(nx-1) + j) - coeffY*g(j*dx,ly,(time-1)*dt,testcase)
        F2((ny-2)*(nx-1) + j) = F2((ny-2)*(nx-1) + j) - coeffY*g(j*dx,ly,(time-1)*dt,testcase)
     END DO

     F1 = F1 + U1/dt
     F2 = F2 + U2/dt
     CALL Uzawa(P0,F1,F2,U1,U2,P)

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
  DEALLOCATE(U1,U2,F1,F2,P,P0)

  CALL CPU_TIME(t2)
  WRITE(*,*) "Temps total (s) : ",t2-t1

END PROGRAM main
