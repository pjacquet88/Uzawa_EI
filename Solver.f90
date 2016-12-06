MODULE Solver
USE Tools
USE Functions
IMPLICIT NONE

CONTAINS
  SUBROUTINE gradConjA(X,RHS)
    IMPLICIT NONE
    REAL*8,DIMENSION(:),INTENT(OUT) :: X
    REAL*8,DIMENSION(:),INTENT(IN) :: RHS
    REAL*8,DIMENSION(:),ALLOCATABLE :: R,D,Omega,Rnext
    REAL*8 :: alpha,beta,normR
    INTEGER :: n,iter

    n = SIZE(X)
    ALLOCATE(R(n),D(n),Omega(n),Rnext(n))

    X = 0
    iter=0
    R = matmulA(X) - RHS
    D = R
    normR = norm(R)

    DO WHILE ((normR > epsilonGC).AND.(iter<nIterMaxGC))
       Omega = matmulA(D)
       alpha = dot_product(D,R)/dot_product(D,Omega)
       X = X - alpha*D
       Rnext = R - alpha*Omega
       beta = dot_product(Rnext,Rnext)/dot_product(R,R)
       D = Rnext + beta*D
       R = Rnext
       normR = norm(R)
       iter = iter+1
    END DO

    WRITE(*,*) "RESIDU GC : ",normR
    WRITE(*,*) "ITERATION GC : ",iter

    DEALLOCATE(R,D,Omega,Rnext)
  END SUBROUTINE GradConjA

  SUBROUTINE Uzawa(P0,F1,F2,U1,U2,P,t)
    IMPLICIT NONE
    REAL*8,DIMENSION(:),INTENT(IN)::P0,F1,F2
    REAL*8,INTENT(IN)::t
    REAL*8,DIMENSION((nx-1)*(ny-1)),INTENT(OUT)::U1,U2
    REAL*8,DIMENSION(nx*ny),INTENT(OUT)::P
    REAL*8,DIMENSION((nx-1)*(ny-1))::RHS1,RHS2
    INTEGER::iter

    iter=1
    P = P0
    DO WHILE ((mu*norm((B1(U1,t) + B2(U2,t))) > epsilonUZ).AND.(iter<nIterMaxUZ))
      RHS1 = F1 - B1transpose(P)
      WRITE(*,*) "RHS1 : max ",maxval(F1)," min ",minval(F1)
      CALL gradConjA(U1,RHS1)
      RHS2 = F2 - B2transpose(P)
      WRITE(*,*) "RHS2 : max ",maxval(RHS2)," min ",minval(RHS2)
      CALL gradConjA(U2,RHS2)

      WRITE(*,*) "UX : max ",maxval(U1)," min ",minval(U1)
      WRITE(*,*) "UY : max ",maxval(U2)," min ",minval(U2)

      P = P + mu*(B1(U1,t) + B2(U2,t))
      WRITE(*,*) "P : max ",maxval(P)," min ",minval(P)

      iter = iter + 1
    END DO

    WRITE(*,*) "RESIDU UZAWA : ",mu*norm((B1(U1,t) + B2(U2,t)))
    WRITE(*,*) "ITERATION UZAWA : ",iter

  END SUBROUTINE


  FUNCTION B1transpose(P)
    IMPLICIT NONE
    REAL*8,DIMENSION(:),INTENT(IN)               ::   P
    REAL*8,DIMENSION((nx-1)*(ny-1)) :: B1transpose
    INTEGER :: i,j

    DO i=1,nx-1
      DO j=1,ny-1
        B1transpose(bij(i,j,nx-1))=P(bij(i+1,j,nx))+P(bij(i+1,j+1,nx))-P(bij(i,j+1,nx))-P(bij(i,j,nx))
        B1transpose=B1transpose/(2*dx)
      END DO
    END DO

  END FUNCTION B1transpose

  FUNCTION B2transpose(P)
    IMPLICIT NONE
    REAL*8,DIMENSION(:),INTENT(IN)               ::   P
    REAL*8,DIMENSION((nx-1)*(ny-1)) :: B2transpose
    INTEGER :: i,j

    DO i=1,nx-1
      DO j=1,ny-1
        B2transpose(bij(i,j,nx-1))=-P(bij(i+1,j,nx))+P(bij(i+1,j+1,nx))+P(bij(i,j+1,nx))-P(bij(i,j,nx))
        B2transpose=B2transpose/(2*dy)
      END DO
    END DO

  END FUNCTION B2transpose

  FUNCTION B1(U1,t)
    IMPLICIT NONE
    REAL*8,DIMENSION(:),INTENT(IN)::U1
    REAL*8,INTENT(IN)::t
    REAL*8,DIMENSION(nx*ny)::B1
    INTEGER::i,j

    B1(bij(1,1,nx)) = -(U1(bij(1,1,nx-1)) + gx(dx,0.d0,t,testcase) - gx(0.d0,dy,t,testcase) - gx(0.d0,0.d0,t,testcase))/(2*dx)
    B1(bij(1,ny,nx)) = -(gx(dx,ly,t,testcase) + U1(bij(1,ny-1,nx-1)) - gx(0.d0,ly,dt,testcase) - gx(0.d0,ly-dy,t,testcase))/(2*dx)
    B1(bij(nx,1,nx)) = -(gx(lx,0.d0,t,testcase) + gx(lx,dy,t,testcase) - U1(bij(nx-1,1,nx-1)) - gx(lx-dx,0.d0,t,testcase))/(2*dx)
    B1(bij(nx,ny,nx)) = -(gx(lx,ly,t,testcase) + gx(lx,ly-dy,t,testcase) - gx(lx-dx,ly,t,testcase) - U1(bij(nx-1,ny-1,nx-1)))/(2*dx)

    DO i=2,nx-1
      B1(bij(i,1,nx)) = -(U1(bij(i,1,nx-1)) + gx((i+1)*dx,0.d0,t,testcase) - U1(bij(i-1,1,nx-1)) - gx(i*dx,0.d0,t,testcase))/(2*dx)
      DO j=2,ny-1
        B1(bij(i,j,nx)) = -(U1(bij(i,j,nx-1)) + U1(bij(i,j-1,nx-1)) - U1(bij(i-1,j,nx-1)) - U1(bij(i-1,j-1,nx-1)))/(2*dx)
      END DO
      B1(bij(i,ny,nx)) = -(gx((i+1)*dx,ly,t,testcase) + U1(bij(i,ny-1,nx-1)) - gx(i*dx,ly,t,testcase) - &
      U1(bij(i-1,ny-1,nx-1)))/(2*dx)
    END DO

    DO j=2,ny-1
      B1(bij(1,j,nx)) = -(U1(bij(1,j,nx-1)) + U1(bij(1,j-1,nx-1)) - gx(0.d0,(j+1)*dy,t,testcase) - gx(0.d0,j*dy,t,testcase))/(2*dx)
      B1(bij(nx,j,nx)) = -(gx(lx,(j+1)*dy,t,testcase) + gx(lx,j*dy,t,testcase) - U1(bij(nx-1,j,nx-1)) - &
      U1(bij(nx-1,j-1,nx-1)))/(2*dx)
    END DO

  END FUNCTION B1

  FUNCTION B2(U2,t)
    IMPLICIT NONE
    REAL*8,DIMENSION(:),INTENT(IN)::U2
    REAL*8,INTENT(IN)::t
    REAL*8,DIMENSION(nx*ny)::B2
    INTEGER::i,j

    B2(bij(1,1,nx)) = -(U2(bij(1,1,nx-1)) - gy(dx,0.d0,t,testcase) + gy(0.d0,dy,t,testcase) - gy(0.d0,0.d0,t,testcase))/(2*dy)
    B2(bij(1,ny,nx)) = -(gy(dx,ly,t,testcase) - U2(bij(1,ny-1,nx-1)) + gy(0.d0,ly,dt,testcase) - gy(0.d0,ly-dy,t,testcase))/(2*dy)
    B2(bij(nx,1,nx)) = -(-gy(lx,0.d0,t,testcase) + gy(lx,dy,t,testcase) + U2(bij(nx-1,1,nx-1)) - gy(lx-dx,0.d0,t,testcase))/(2*dy)
    B2(bij(nx,ny,nx)) = -(gy(lx,ly,t,testcase) - gy(lx,ly-dy,t,testcase) + gy(lx-dx,ly,t,testcase) - U2(bij(nx-1,ny-1,nx-1)))/(2*dy)

    DO i=2,nx-1
      B2(bij(i,1,nx)) = -(U2(bij(i,1,nx-1)) - gy((i+1)*dx,0.d0,t,testcase) + U2(bij(i-1,1,nx-1)) - gy(i*dx,0.d0,t,testcase))/(2*dy)
      DO j=2,ny-1
        B2(bij(i,j,nx)) = -(U2(bij(i-1,j,nx-1)) + U2(bij(i,j,nx-1)) - U2(bij(i-1,j-1,nx-1)) - U2(bij(i,j-1,nx-1)))/(2*dy)
      END DO
      B2(bij(i,ny,nx)) = -(gy((i+1)*dx,ly,t,testcase) - U2(bij(i,ny-1,nx-1)) + gy(i*dx,ly,t,testcase) - &
      U2(bij(i-1,ny-1,nx-1)))/(2*dy)
    END DO

    DO j=2,ny-1
      B2(bij(1,j,nx)) = -(U2(bij(1,j,nx-1)) - U2(bij(1,j-1,nx-1)) + gy(0.d0,(j+1)*dy,t,testcase) - gy(0.d0,j*dy,t,testcase))/(2*dy)
      B2(bij(nx,j,nx)) = -(gy(lx,(j+1)*dy,t,testcase) - gy(lx,j*dy,t,testcase) + U2(bij(nx-1,j,nx-1)) - &
      U2(bij(nx-1,j-1,nx-1)))/(2*dy)
    END DO
  END FUNCTION B2




END MODULE Solver
