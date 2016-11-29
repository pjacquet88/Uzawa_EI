PROGRAM main_poisson
  USE Functions
  USE Solver
  use Tools
  IMPLICIT NONE
  !REAL*8,DIMENSION(:,:),ALLOCATABLE :: A,B,C
  REAL*8,DIMENSION(:),ALLOCATABLE :: U,X,RHS
  INTEGER :: i,j,time
  REAL*8 :: xi,yi,t1,t2,normDiff,normExact

  CALL CPU_TIME(t1)
  CALL read_param("param.txt")
  ALLOCATE(U(nl),X(nl),RHS(nl))



  U = 0
  WRITE(*,*) "nl ",size(RHS)
  DO time=1,10

     RHS=0
     DO i=1,ny-1
        DO j=1,nx-1
           RHS((i-1)*(nx-1)+j) = f(j*dx,i*dy,(time-1)*dt,p)
        END DO
        RHS((i-1)*(nx-1)+1) = RHS((i-1)*(nx-1)+1) - coeffX*g(0.d1,i*dy,(time-1)*dt,p)
        RHS(i*(nx-1)) = RHS(i*(nx-1)) - coeffX*g(lx,i*dy,(time-1)*dt,p)
     END DO
     DO j=1,nx-1
        RHS(j) = RHS(j) - coeffY*g(j*dx,0.d1,(time-1)*dt,p)
        RHS((ny-2)*(nx-1) + j) = RHS((ny-2)*(nx-1) + j) - coeffY*g(j*dx,ly,(time-1)*dt,p)
     END DO

     RHS = RHS + U/dt
     CALL gradConjA(U,RHS)

  END DO

  OPEN(unit=12,file='sol.dat',form='formatted',status='unknown',action='write')
  OPEN(unit=13,file='solExact.dat',form='formatted',status='unknown',action='write')

  normDiff = 0.
  normExact = 0.
  DO i=1,ny-1
     DO j=1,nx-1
        xi = dx*j
        yi = dy*i
        WRITE(12,*) xi,yi,U((i-1)*(nx-1)+j)
        IF (p==2) THEN
           WRITE(13,*) xi,yi,sin(xi)+cos(yi)
           normDiff = max(abs(U((i-1)*(nx-1)+j)-(sin(xi)+cos(yi))),normDiff)
           normExact =  max(abs(sin(xi)+cos(yi)),normExact)
        ELSE IF (p==1) THEN
           WRITE(13,*) xi,yi,xi*(1-xi)*yi*(1-yi)
           normDiff = max(abs(U((i-1)*(nx-1)+j)-(xi*(1-xi)*yi*(1-yi))),normDiff)
           normExact = max(xi*(1-xi)*yi*(1-yi),normExact)
        END IF
     END DO
     WRITE(12,*)
     WRITE(13,*)
  END DO
  IF ((p==1).OR.(p==2)) THEN
     WRITE (*,*) "Erreur relative en norme infinie : ",normDiff/normExact
  END IF

  !!stop
  CLOSE(12)
  CLOSE(13)
  CLOSE(14)
  DEALLOCATE(U,X,RHS)

  CALL CPU_TIME(t2)
  WRITE(*,*) "Temps total (s) : ",t2-t1

END PROGRAM main_poisson
