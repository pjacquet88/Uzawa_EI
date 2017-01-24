!MAIN POUR LE PROBLEME DE POISSON

PROGRAM main_poisson
  USE Functions
  USE Solver
  use Tools
  IMPLICIT NONE
  !REAL*8,DIMENSION(:,:),ALLOCATABLE :: A,B,C
  REAL*8,DIMENSION(:),ALLOCATABLE :: U,X,RHS,Ue
  INTEGER :: i,j,time,k
  REAL*8 :: xi,yi,t1,t2,normDiff,normExact
  CALL CPU_TIME(t1)


OPEN(unit=56,file='errorLap.dat',action='write')
  CALL read_param("param.txt")
DO nx=8,100

  print*,'nx= ',nx
  ny=nx
  nl=(nx-1)*(ny-1)
  dx=lx/nx
  dy=ly/ny
  coeffDiag=2*visc*(1/(dx**2)+1/(dy**2))
  coeffX=-visc*(1/(dx**2))
  coeffY=-visc*(1/(dy**2))


  ALLOCATE(U(nl),X(nl),RHS(nl),Ue(nl))



  U = 0
  WRITE(*,*) "nl ",size(RHS)
  DO time=1,1

     RHS=0
     DO i=1,ny-1
        DO j=1,nx-1
           RHS((i-1)*(nx-1)+j) = f(j*dx,i*dy,(time-1)*dt,testcase)
        END DO
        RHS((i-1)*(nx-1)+1) = RHS((i-1)*(nx-1)+1) - coeffX*g(0.d1,i*dy,(time-1)*dt,testcase)
        RHS(i*(nx-1)) = RHS(i*(nx-1)) - coeffX*g(lx,i*dy,(time-1)*dt,testcase)
     END DO
     DO j=1,nx-1
        RHS(j) = RHS(j) - coeffY*g(j*dx,0.d1,(time-1)*dt,testcase)
        RHS((ny-2)*(nx-1) + j) = RHS((ny-2)*(nx-1) + j) - coeffY*g(j*dx,ly,(time-1)*dt,testcase)
     END DO

     RHS = RHS! + U/dt
     CALL gradConjA(U,RHS)

  END DO

  !OPEN(unit=12,file='sol.dat',form='formatted',status='unknown',action='write')
  !OPEN(unit=13,file='solExact.dat',form='formatted',status='unknown',action='write')

  normDiff = 0.
  normExact = 0.
  DO i=1,ny-1
     DO j=1,nx-1
        xi = dx*j
        yi = dy*i
        !WRITE(12,*) xi,yi,U((i-1)*(nx-1)+j)
        IF (testcase==2) THEN
           !WRITE(13,*) xi,yi,sin(xi)+cos(yi)
           !normDiff = max(abs(U((i-1)*(nx-1)+j)-(sin(xi)+cos(yi))),normDiff)
           !normExact =  max(abs(sin(xi)+cos(yi)),normExact)
           Ue(bij(j,i,nx-1))=sin(xi)+cos(yi)
        ELSE IF (testcase==1) THEN
           !WRITE(13,*) xi,yi,xi*(1-xi)*yi*(1-yi)
           !normDiff = max(abs(U((i-1)*(nx-1)+j)-(xi*(1-xi)*yi*(1-yi))),normDiff)
           !normExact = max(xi*(1-xi)*yi*(1-yi),normExact)
            Ue(bij(j,i,nx-1))=xi*(1-xi)*yi*(1-yi)
        END IF
     END DO
    ! WRITE(12,*)
     !WRITE(13,*)
  END DO
  IF ((testcase==1).OR.(testcase==2)) THEN
     WRITE (56,*),dx, (norm(Ue-U)/norm(Ue))
  END IF

  !!stop
  !CLOSE(12)
!  CLOSE(13)
!  CLOSE(14)
  DEALLOCATE(U,X,RHS,Ue)
END DO

close(56)

  CALL CPU_TIME(t2)
  WRITE(*,*) "Temps total (s) : ",t2-t1

END PROGRAM main_poisson
