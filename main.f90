PROGRAM main
  USE Solver
  USE Tools
  IMPLICIT NONE
  REAL*8,DIMENSION(:),ALLOCATABLE:: U1,U2,F1,F2,P,P0,P_int,TEST,KER,Uex,Uey
  INTEGER :: i,j,time
  REAL*8 :: xi,yi,t1,t2,normDiff,normExact,x,y,integrale,erreurvx,erreurvy

  CALL CPU_TIME(t1)
  CALL read_param("param.txt")
  ALLOCATE(U1((nx-1)*(ny-1)),U2((nx-1)*(ny-1)),F1((nx-1)*(ny-1)),F2((nx-1)*(ny-1)),P(nx*ny),P0(nx*ny),TEST((nx)*(ny))&
  ,KER((nx-1)*(ny-1)),Uex((nx-1)*(ny-1)),Uey((nx-1)*(ny-1)))

  U1 = 0
  U2 = 0
  P0 = 0
  P_int=P0

  DO time=1,1    ! A CHANGER POUR PASSER EN STATIONNAIRE
     F1=0
     F2=0
     DO i=1,ny-1
        DO j=1,nx-1
           F1((i-1)*(nx-1)+j) = fx(j*dx,i*dy,(time-1)*dt,testcase)
           F2((i-1)*(nx-1)+j) = fy(j*dx,i*dy,(time-1)*dt,testcase)
        END DO
        F1((i-1)*(nx-1)+1) = F1((i-1)*(nx-1)+1) - coeffX*gx(0.d1,i*dy,(time-1)*dt,testcase)
        F2((i-1)*(nx-1)+1) = F2((i-1)*(nx-1)+1) - coeffX*gy(0.d1,i*dy,(time-1)*dt,testcase)
        F1(i*(nx-1)) = F1(i*(nx-1)) - coeffX*gx(lx,i*dy,(time-1)*dt,testcase)
        F2(i*(nx-1)) = F2(i*(nx-1)) - coeffX*gy(lx,i*dy,(time-1)*dt,testcase)
     END DO
     DO j=1,nx-1
        F1(j) = F1(j) - coeffY*hx(j*dx,0.d1,(time-1)*dt,testcase)
        F2(j) = F2(j) - coeffY*hy(j*dx,0.d1,(time-1)*dt,testcase)
        F1((ny-2)*(nx-1) + j) = F1((ny-2)*(nx-1) + j) - coeffY*hx(j*dx,ly,(time-1)*dt,testcase)
        F2((ny-2)*(nx-1) + j) = F2((ny-2)*(nx-1) + j) - coeffY*hy(j*dx,ly,(time-1)*dt,testcase)
     END DO

  !   F1 = F1 + U1/dt             ! A DECOMMENTER POUR PASSER EN STATIONNAIRE

  !   F2 = F2 + U2/dt             ! A CHANGER POUR PASSER EN STATIONNAIRE

     CALL Uzawa(P_int,F1,F2,U1,U2,P,(time-1)*dt)
     P_int=P
  END DO





!*****************************   SORTIE DES RESULTATS FINAUX *************************************
open(unit=3,file='pression.dat',action='write')
do j=1,ny
  do i=1,nx
    x=(i-0.5)*dx
    y=(j-0.5)*dy
write(3,*),x,y,P(bij(i,j,nx))
  end do
  write(3,*),' '
end do
close(3)

open(unit=3,file='vitesse_x.dat',action='write')
do j=1,ny-1
  do i=1,nx-1
    x=(i)*dx
    y=(j)*dy
write(3,*),x,y,U1(bij(i,j,nx-1))
  end do
  write(3,*),' '
end do
close(3)

open(unit=3,file='vitesse_y.dat',action='write')
do j=1,ny-1
  do i=1,nx-1
    x=(i)*dx
    y=(j)*dy
write(3,*),x,y,U2(bij(i,j,nx-1))
  end do
  write(3,*),' '
end do
close(3)
!*******************************************************************************

!********************** CALCULS FINAUX******************************************
! PARTIE NON PARAMETRéE (Calcul d'erreur à adapter avec le cas testé puis à recompiler)

  CALL CPU_TIME(t2)
  WRITE(*,*) "Temps total (s) : ",t2-t1

integrale=0.0d0
do i=1,size(P)
  integrale=integrale+P(i)*dx*dy
end do
integrale=integrale/(Lx*Ly)

print*,'Integrale de P = ', integrale


  do i=1,nx-1
      x=(i-0.5)*dx
    do j=1,ny-1
      y=(j-0.5)*dy
      Uex(bij(i,j,nx-1))=(Ly-y)*y
      Uey(bij(i,j,nx-1))=0.0
        end do
  end do

erreurvx=0.0
erreurvy=0.0

erreurvx=norm(Uex-U1)/norm(Uex)
erreurvy=norm(Uey-U2)/norm(Uey)
PRINT*,'ATTENTION LES RESULTATS CI DESSOUS DEPENDENT DE L IMPLEMENTATION MANUELLE DES SOLUTIONS EXACTES, A LA FIN DE MAIN'
print*,'erreur vitesse selon x',erreurvx
print*,'erreur vitesse selon y',erreurvy
DEALLOCATE(U1,U2,F1,F2,P,P0)
END PROGRAM main
