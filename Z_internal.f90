SUBROUTINE partition(T, nu, n, label, F, Fmax, ele, symm, Trot, Zvib, Zrot, Zint)

!Compute the internal partition function at the kth temperatures
!"label" refers to the specific hydrogenated molecule in exam and F is the number
!of its vibrational degrees of freedom
USE thermo_constants

IMPLICIT NONE

!Computing internal partition function

!Declaration of input variables
INTEGER(4) :: n, label, F, Fmax
REAL(8) :: T              
REAL(8), DIMENSION(0:n,Fmax) :: nu 
REAL(8), DIMENSION(0:n,2) :: ele
REAL(8), DIMENSION(0:n) :: symm
REAL(8), DIMENSION(0:n,3) :: Trot

!Declaration of side and output variables
INTEGER(4) :: Nele
REAL(8) :: theta, expo
REAL(8) :: Zvib, Zele, Zrot
REAL(8) :: Zint

!Declaration of counters
INTEGER :: i

!Computing Vibrational Partition function
Zvib=1.0d0
DO i=1,F
	expo = (-h*c*nu(label,i)*100.0d0)/(kb*T )
        Zvib = Zvib*( 1.0d0/(1.0d0-exp(expo) ) )
END DO
!Computing Rotational Partition Function
 theta = ( Trot(label,1)*Trot(label,2)*Trot(label,3) )**onethird
 Zrot = (1.0d0/symm(label))*((T/theta)**threesecond)*sqrt(pi)
!Computing Electronical Partition Function
Nele = ele(label,1)+ele(label,2)
Zele = 1 + mod(Nele,2)
!Computing Total Partition Function
Zint=Zvib*Zrot*Zele

END SUBROUTINE
          
