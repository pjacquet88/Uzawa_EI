MODULE Tools
  IMPLICIT NONE

  INTEGER,PUBLIC          ::    nx,ny,nl,p
  REAL*8,PUBLIC           ::    lx,ly,dx,dy
  REAL*8,PUBLIC           ::    t_end,dt, visc, mu
  REAL*8,PUBLIC,PARAMETER ::    gravity = 9.80665
  REAL*8,PUBLIC,PARAMETER ::    pi = 3.14159265359
  REAL*8,PUBLIC           ::    coeffX,coeffY,coeffDiag
  REAL*8,PUBLIC           ::    epsilonGC,epsilonUZ
  INTEGER,PUBLIC          ::    nIterMaxGC, nIterMaxUZ

  PUBLIC :: bij ,read_param,  norm, TEST_read_param , printvector
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

    OPEN(unit=11, file=filename,form='formatted', action="read", status="unknown")
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
    READ(11,*) dummy_string
    READ(11,*) p
    CLOSE(11)

    nl=(nx-1)*(ny-1)
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

  FUNCTION norm(V)
    IMPLICIT NONE
    REAL*8,DIMENSION(:),INTENT(IN) :: V
    REAL*8                         :: norm

    norm=DOT_PRODUCT(V,V)
    norm=SQRT(norm)

  END FUNCTION norm

  SUBROUTINE printvector(U,Nx,Ny,dx,dy,N)
    REAL*8,DIMENSION(Nx*Ny),INTENT(in)::U
    INTEGER,INTENT(in)::Nx,Ny,N
    REAL*8, INTENT(in)::dx,dy
    INTEGER::i,j
    CHARACTER(len=20)::F_NAME

    IF (N<10) THEN
       F_NAME='fichier/T'
       WRITE(F_NAME (10:10),'(I1)') N
       F_NAME(11:14)= '.dat'

    ELSEIF ((N>=10).AND.(N<100)) THEN
       F_NAME='fichier/T'
       WRITE(F_NAME (10:11),'(I2)') N
       F_NAME(12:15)= '.dat'

    ELSEIF ((N>=100).AND.(N<1000)) THEN
       F_NAME='fichier/T'
       WRITE(F_NAME (10:12),'(I3)') N
       F_NAME(13:16)= '.dat'

    ELSE IF ((N>=1000).AND.(N<10000)) THEN
       F_NAME='fichier/T'
       WRITE(F_NAME (10:13),'(I4)') N
       F_NAME(14:17)= '.dat'

    ELSE IF ((N>=10000).AND.(N<100000)) THEN
       F_NAME='fichier/T'
       WRITE(F_NAME (10:14),'(I5)') N
       F_NAME(15:18)= '.dat'
    END IF

    OPEN(unit=2, file=F_NAME, action="write")


    DO i=1,Nx
       DO j=1,Ny


          WRITE(2,*),i*dx,j*dy,U(bij(i,j,Ny))

       END DO
       WRITE(2,*)
    END DO
    CLOSE(2)
  END SUBROUTINE printvector

  FUNCTION matmulA(X)
    IMPLICIT NONE
    REAL*8,DIMENSION(:),INTENT(IN) :: X
    REAL*8,DIMENSION(SIZE(X)) :: matmulA
    INTEGER :: i,j

    DO j=1,ny-1
       IF (j==1) THEN
          DO i=1,nx-1
             IF (i==1) THEN
                matmulA((j-1)*(nx-1)+i) = coeffDiag*X(1) + coeffX*X(2) + coeffY*X(nx)
             ELSE IF (i==nx) THEN
                matmulA((j-1)*(nx-1)+i) = coeffX*X(nx-2) + coeffDiag*X(nx-1) + coeffY*X(2*(nx-1))
             ELSE
                matmulA((j-1)*(nx-1)+i) = coeffX*X(i-1) + coeffDiag*X(i) + coeffX*X(i+1) + coeffY*X(nx-1+i)
             END IF
          END DO
       ELSE IF (j==ny-1) THEN
          DO i=1,nx-1
             IF (i==1) THEN
                matmulA((j-1)*(nx-1)+i) = coeffY*X(((ny-1)-2)*(nx-1)+1) + coeffDiag*X(((ny-1)-1)*(nx-1)+1) + &
                coeffX*X(((ny-1)-1)*(nx-1)+2)
             ELSE IF (i==nx-1) THEN
                matmulA((j-1)*(nx-1)+i) = coeffY*X(((ny-1)-1)*(nx-1)) + coeffX*X((ny-1)*(nx-1)-1) + &
                coeffDiag*X((ny-1)*(nx-1))
             ELSE
                matmulA((j-1)*(nx-1)+i) = coeffY*X(((ny-1)-2)*(nx-1)+i) + coeffX*X(((ny-1)-1)*(nx-1) + i-1) + &
                coeffDiag*X(((ny-1)-1)*(nx-1) + i) + coeffX*X(((ny-1)-1)*(nx-1) + i+1)
             END IF
          END DO
       ELSE
          DO i=1,nx-1
             IF (i==1) THEN
                matmulA((j-1)*(nx-1)+i) = coeffY*X((j-2)*(nx-1)+1) + coeffDiag*X((j-1)*(nx-1)+1) + &
                coeffX*X((j-1)*(nx-1)+2) + coeffY*X(j*(nx-1)+1)
             ELSE IF (i==nx-1) THEN
                matmulA((j-1)*(nx-1)+i) = coeffY*X((j-2)*(nx-1)+(nx-1)) + coeffX*X((j-1)*(nx-1)+(nx-1)-1) + &
                coeffDiag*X((j-1)*(nx-1)+(nx-1)) + coeffY*X(j*(nx-1)+(nx-1))
             ELSE
                matmulA((j-1)*(nx-1)+i) = coeffY*X((j-2)*(nx-1)+i) + coeffX*X((j-1)*(nx-1) + i-1) + &
                coeffDiag*X((j-1)*(nx-1) + i) + coeffX*X((j-1)*(nx-1) + i+1) + coeffY*X(j*(nx-1)+i)
             END IF
          END DO
       END IF
    END DO
  END FUNCTION matmulA



END MODULE Tools
