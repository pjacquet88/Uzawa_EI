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
    CASE(0)
       fx=0.0
    CASE(1)
       fx = 1
     CASE(2)
       fx = 1
     CASE(3)
       fx = 0
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
    CASE(0)
       fy = 0.0 ! -gravity
    CASE(1)
       fy = 1
     CASE(2)
       fy = 1
     CASE(3)
       fy = 0
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
    CASE(0)
       gx = 23.12
    CASE(1)
       gx = x
     CASE(2)
       gx = x
     CASE(3)
       gx =Ly*Ly*(1-y/Ly)*y/Ly
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
    CASE(0)
       gy = 21.61
    CASE(1)
       gy = -y
     CASE(2)
       gy = y
     CASE(3)
       gy = 0.0
    CASE default
       gy = -1
       WRITE(*,*) "bad case"
    END SELECT
  END FUNCTION gy

  FUNCTION hx(x,y,t,n)
    IMPLICIT NONE
    REAL*8,INTENT(IN) :: x,y,t
    INTEGER,INTENT(IN) :: n
    REAL*8 :: hx

    SELECT CASE(n)
    CASE(0)
       hx = 23.12
    CASE(1)
       hx = x
     CASE(2)
       hx = x
     CASE(3)
       hx = 0
    CASE default
       hx = -1
       WRITE(*,*) "bad case"
    END SELECT
  END FUNCTION hx

  FUNCTION hy(x,y,t,n)
    IMPLICIT NONE
    REAL*8,INTENT(IN) :: x,y,t
    INTEGER,INTENT(IN) :: n
    REAL*8 :: hy

    SELECT CASE(n)
    CASE(0)
       hy = 21.61
    CASE(1)
       hy = -y
     CASE(2)
       hy = y
     CASE(3)
       hy = 0
    CASE default
       hy = -1
       WRITE(*,*) "bad case"
    END SELECT
  END FUNCTION hy



END MODULE Functions
