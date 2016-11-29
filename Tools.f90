MODULE Tools
IMPLICIT NONE

INTEGER,PUBLIC          ::    nx,ny
REAL*8,PUBLIC           ::    lx,ly,dx,dy
REAL*8,PUBLIC           ::    t_end,dt, visc, mu
REAL*8,PUBLIC,PARAMETER ::    g = 9.80665
REAL*8,PUBLIC,PARAMETER ::    pi = 3.14159265359
REAL*8,PUBLIC           ::    coeffX,coeffY,coeffDiag
REAL*8,PUBLIC           ::    epsilonGC,epsilonUZ
INTEGER,PUBLIC          ::    nIterMaxGC, nIterMaxUZ

PUBLIC :: bij ,read_param,  norme, TEST_read_param
CONTAINS

FUNCTION bij(i,j,N)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::i,j,N
  INTEGER           :: bij

  bij=i+(j-1)*N
END FUNCTION bij

SUBROUTINE read_param(filename)
  IMPLICIT NONE
	CHARACTER(len=*), INTENT(IN)::filename
  CHARACTER(len=1)            :: dummy_string

	open(unit=11, file=filename,form='formatted', action="read", status="unknown")
  READ(11,*) dummy_string
  READ(11,*) dummy_string
  READ(11,*) dummy_string
  READ(11,*) lx
  READ(11,*) dummy_string
  READ(11,*) ly
  READ(11,*) dummy_string
  READ(11,*) nx
  READ(11,*) dummy_string
  READ(11,*) ny
  READ(11,*) dummy_string
  READ(11,*) t_end
  READ(11,*) dummy_string
  READ(11,*) dt
  READ(11,*) dummy_string
  READ(11,*) visc
  READ(11,*) dummy_string
  READ(11,*) mu
  READ(11,*) dummy_string
  READ(11,*) epsilonGC
  READ(11,*) dummy_string
  READ(11,*) epsilonUZ
  READ(11,*) dummy_string
  READ(11,*) nIterMaxGC
  READ(11,*) dummy_string
  READ(11,*) nIterMaxUZ
  close(11)

  dx=lx/nx
  dy=ly/ny
  coeffDiag=1/dt+2*visc*(1/(dx**2)+1/(dy**2))
  coeffX=-visc*(1/(dx**2))
  coeffY=-visc*(1/(dy**2))
END SUBROUTINE read_param

SUBROUTINE TEST_read_param()
  IMPLICIT NONE

  PRINT*,'lx',lx
  PRINT*,'ly',ly
  PRINT*,'nx',nx
  PRINT*,'ny',ny
  PRINT*,'t_end',t_end
  PRINT*,'dt',dt
  PRINT*,'visc',visc
  PRINT*,'mu',mu
  PRINT*,'epsilonGC',epsilonGC
  PRINT*,'epsilonUZ',epsilonUZ
  PRINT*,'nIterMaxGC',nIterMaxGC
  PRINT*,'nIterMaxUZ',nIterMaxUZ

END SUBROUTINE TEST_read_param

FUNCTION norme(V)
  IMPLICIT NONE
  REAL*8,DIMENSION(:),INTENT(IN) :: V
  REAL*8                         :: norme

  norme=dot_product(V,V)
  norme=sqrt(norme)

END FUNCTION








END MODULE
