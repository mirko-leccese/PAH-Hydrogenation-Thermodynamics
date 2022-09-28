SUBROUTINE dbpressure(T, molar, n, label, dbp)

! This subroutine computes the De Broglie Pressure of a given PAH specie (label)
! at given thermodynamical conditions:
!
! Args:
!   T (real): temperature
!   molar (array): array of PAH molar masses
!   n (integer): number of hydrogenation levels (which sets the dimension of the array, i.e. the number of possible 
!                hydrogenated species, 24 for coronene)
!   label (integer): number identyfing the specific hydrogenated specie for which we want to compute "dbp"
!
! Out:
!   dbp (real): DeBroglie pressure of the target hydrogenated specie
!   

USE thermo_constants

IMPLICIT NONE

! Declaring input variables:
INTEGER(4) :: n, label
REAL(8), DIMENSION(0:n) :: molar
REAL(8) :: T

! Declaring side and output variables:
REAL(8) :: m, rad 
REAL(8) :: vb, lambda, dbp

! Initialization
lambda=0.0d0
vb=0.0d0
dbp=0.0d0

! Converting molar mass(amu) to Kg
m = molar(label)*amutog*gtokg   
rad = sqrt(2.0d0*pi*m*kb*T)

lambda = h/rad
vb = lambda**3.0d0
dbp = ( (kb*T)/vb )*Patobar

END SUBROUTINE
