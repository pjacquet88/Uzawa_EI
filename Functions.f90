MODULE Functions
  USE Tools
  IMPLICIT NONE

CONTAINS
  FUNCTION f(x,y,t,n)
    IMPLICIT NONE
    REAL*8,INTENT(IN) :: x,y,t
    INTEGER,INTENT(IN) :: n
    real*8 :: f

    SELECT CASE(n)
    CASE(1)
       f = 2*(y-y*y+x-x*x)
    CASE(2)
       f = sin(x) + cos(y)
    CASE(3)
       !f = exp(-(x-Lx/2)*(x-Lx/2))*exp(-(y-Ly/2)*(y-Ly/2))*cos(pi/2*t)
    CASE default
       f = -1
       WRITE(*,*) "bad case"
    END SELECT
  END FUNCTION f

  FUNCTION g(x,y,t,n)
    IMPLICIT NONE
    REAL*8,INTENT(IN) :: x,y,t
    INTEGER,INTENT(IN) :: n
    REAL*8 :: g

    SELECT CASE(n)
    CASE(1)
       g = 0
    CASE(2)
       g = sin(x) + cos(y)
    CASE(3)
       g = 0
    CASE default
       g = -1
       WRITE(*,*) "bad case"
    END SELECT
  END FUNCTION g

END MODULE Functions
