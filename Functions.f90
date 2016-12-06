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
    CASE default
       f = -1
       WRITE(*,*) "bad case"
    END SELECT
  END FUNCTION f

  FUNCTION fx(x,y,t,n)
    IMPLICIT NONE
    REAL*8,INTENT(IN) :: x,y,t
    INTEGER,INTENT(IN) :: n
    real*8 :: fx

    SELECT CASE(n)
    CASE(1)
       fx = 1
    CASE default
       fx = -1
       WRITE(*,*) "bad case"
    END SELECT
  END FUNCTION fx

  FUNCTION fy(x,y,t,n)
    IMPLICIT NONE
    REAL*8,INTENT(IN) :: x,y,t
    INTEGER,INTENT(IN) :: n
    real*8 :: fy

    SELECT CASE(n)
    CASE(1)
       fy = 1
    CASE default
       fy = -1
       WRITE(*,*) "bad case"
    END SELECT
  END FUNCTION fy

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

  FUNCTION gx(x,y,t,n)
    IMPLICIT NONE
    REAL*8,INTENT(IN) :: x,y,t
    INTEGER,INTENT(IN) :: n
    REAL*8 :: gx

    SELECT CASE(n)
    CASE(1)
       gx = x
    CASE default
       gx = -1
       WRITE(*,*) "bad case"
    END SELECT
  END FUNCTION gx

  FUNCTION gy(x,y,t,n)
    IMPLICIT NONE
    REAL*8,INTENT(IN) :: x,y,t
    INTEGER,INTENT(IN) :: n
    REAL*8 :: gy

    SELECT CASE(n)
    CASE(1)
       gy = -y
    CASE default
       gy = -1
       WRITE(*,*) "bad case"
    END SELECT
  END FUNCTION gy



END MODULE Functions
