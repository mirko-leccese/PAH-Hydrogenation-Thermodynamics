SUBROUTINE dbpressure(T, molar, n, label, dbp)

!Computing De Broglie Pressure (bar) and Thermal Wavelenght (m) at the kth temperature
!for the ith=label hydrogenated molecule

USE thermo_constants

IMPLICIT NONE

!Declaring input variables
INTEGER(4) :: n, label
REAL(8), DIMENSION(0:n) :: molar
REAL(8) :: T

!Declaring side and output variables
REAL(8) :: m, rad 
REAL(8) :: vb, lambda, dbp

!Initialization
lambda=0.0d0
vb=0.0d0
dbp=0.0d0

m = molar(label)*amutog*gtokg   !Converting molar mass(amu) to Kg
rad = sqrt(2.0d0*pi*m*kb*T)

lambda = h/rad
vb = lambda**3.0d0
dbp = ( (kb*T)/vb )*Patobar

END SUBROUTINE
