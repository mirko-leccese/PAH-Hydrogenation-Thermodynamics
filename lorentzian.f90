SUBROUTINE Lorentzian(s, W, npoint, seq, n)

  USE thermo_constants
  USE hydrogen_constants

INTEGER(4) :: n
LOGICAL :: seq
REAL(8) :: W                    !Width of the Lorentzian
REAL(8) :: lgamma
REAL(8), DIMENSION(0:n,1:2) :: A      !Abundancies
REAL(8), ALLOCATABLE, DIMENSION(:) :: grid 
REAL(8) :: dx, den, Lor
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: fx
REAL(8), ALLOCATABLE, DIMENSION(:) :: Tot
INTEGER(4) :: i, npoint, j, s
CHARACTER(len=20) :: filenameopt
CHARACTER(len=20) :: lorentzname


dx=real(n)/npoint

write(*,*) dx

ALLOCATE(grid(0:npoint))
ALLOCATE(fx(0:n,0:npoint))
ALLOCATE(Tot(0:npoint))

DO i=0,npoint
        grid(i)=i*dx
        !WRITE(*,*) grid(i)
END DO


Tot=0.0d0
fx=0.0d0
A=0.0d0
WRITE(filenameopt,'(I2.2,"-seqopt.out")') s
WRITE(lorentzname, '(I2.2,"-lorentz.out")' ) s  
OPEN(Unit=2000, File=filenameopt)    
OPEN(Unit=2001, File=lorentzname)
lgamma=W/2.0d0
DO i=0,n
	READ(2000,*) (A(i,j), j=1,2)
END DO
DO j=0,n
	DO i=0,npoint
		den=(grid(i)-j)**2+lgamma**2
		Lor=(1.0d0/pi)*(lgamma/den)
		fx(j,i)=A(j,2)*Lor
		!WRITE(*,*) grid(i), fx(j,i) 
	END DO
END DO

DO k=0, npoint
	DO j=0,n
		Tot(k)=Tot(k)+fx(j,k)
	END DO
	write(2001,*) grid(k), Tot(k)
END DO

CLOSE(2000)
CLOSE(2001)

DEALLOCATE(grid)
DEALLOCATE(Tot)
DEALLOCATE(fx)


END SUBROUTINE
                    



