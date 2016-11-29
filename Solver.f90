MODULE Solver
USE Tools
IMPLICIT NONE

SUBROUTINE gradConjA(X,RHS)
  IMPLICIT NONE
  REAL*8,DIMENSION(:),INTENT(OUT) :: X
  REAL*8,DIMENSION(:),INTENT(IN) :: RHS
  REAL*8,DIMENSION(:),ALLOCATABLE :: R,D,Omega,Rnext
  REAL*8 :: alpha,beta,normR
  INTEGER :: n

  n = SIZE(X)
  ALLOCATE(R(n),D(n),Omega(n),Rnext(n))

  X = 0    
  R = matmulA(coeffDiag,coeffX,coeffY,X) - RHS
  D = R
  normR = norm(R)
  
  DO WHILE (normR > epsilonGC)
     Omega = matmulA(coeffDiag,coeffX,coeffY,D)
     alpha = dot_product(D,R)/dot_product(D,Omega)
     X = X - alpha*D
     Rnext = R - alpha*Omega
     beta = dot_product(Rnext,Rnext)/dot_product(R,R)
     D = Rnext + beta*D
     R = Rnext
     normR = norm(R)
  END DO
  WRITE(*,*) normR
  
  DEALLOCATE(R,D,Omega,Rnext)
END SUBROUTINE GradConjA
  
END MODULE Solver