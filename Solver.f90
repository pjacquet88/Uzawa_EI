MODULE Solver
USE Tools
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

  SUBROUTINE Uzawa(P0,F1,F2,U1,U2,P)
    IMPLICIT NONE
    REAL*8,DIMENSION(:),INTENT(IN)::P0,F1,F2
    REAL*8,DIMENSION((nx-1)*(ny-1)),INTENT(OUT)::U1,U2
    REAL*8,DIMENSION(nx*ny),INTENT(OUT)::P
    REAL*8,DIMENSION((nx-1)*(ny-1))::RHS1,RHS2
    INTEGER::iter

    iter=1
    P = P0
    DO WHILE ((mu*norm((B1(U1) + B2(U2))) > epsilonUZ).AND.(iter<nIterMaxUZ))
      RHS1 = F1 - B1transpose(P)
      CALL gradConjA(U1,RHS1)
      RHS2 = F2 - B2transpose(P)
      CALL gradConjA(U2,RHS2)

      P = P + mu*(B1(U1) + B2(U2))
      iter = iter + 1
    END DO

    WRITE(*,*) "RESIDU UZAWA : ",mu*norm((B1(U1) + B2(U2)))
    WRITE(*,*) "ITERATION UZAWA : ",iter

  END SUBROUTINE


END MODULE Solver
