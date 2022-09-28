MODULE  hydrogen_constants

! Thermodynamic parameters of H2 from Nist Webbook
! https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Units=SI&Mask=1#Thermo-Gas
REAL(8), PARAMETER ::  A1 =  33.066178d0,  A2 =  18.563083d0,  A3 = 43.413560d0
REAL(8), PARAMETER ::  B1 = -11.363417d0,  B2 =  12.257357d0,  B3 = -4.293079d0
REAL(8), PARAMETER ::  C1 =  11.432816d0,  C2 =  -2.859786d0,  C3 = 1.272428d0
REAL(8), PARAMETER ::  D1 =  -2.772874d0,  D2 =   0.268238d0,  D3 = -0.096876d0
REAL(8), PARAMETER ::  E1 =  -0.158558d0,  E2 =   1.977990d0,  E3 = -20.533862d0
REAL(8), PARAMETER ::  F1 =  -9.980797d0,  F2 =  -1.147438d0,  F3 = -38.515158d0
REAL(8), PARAMETER ::  G1 = 172.707974d0,  G2 = 156.288133d0,  G3 = 162.081354d0

!Relevant thermodynamical and chemical values:
REAL(8), PARAMETER :: nuh2 = 4543.5913d+2               ! Frequency (m-1) of the vibrational mode
REAL(8), PARAMETER :: molarh2 = 2.01565d0               ! Molar mass H2 (amu)
REAL(8), PARAMETER :: molarh=1.00783d0                  ! Molar mass H (amu) 
REAL(8), PARAMETER :: B=60.9d+2							! Rotational constant of H2 in m^-1
!REAL(8), PARAMETER :: B=159.29239792d+2				! Rotational constant adjucted to ensure continuity with Shomate, m^-1
REAL(8) :: Edft=-1.157311d0                             ! DFT Energy (ZPE Corrected) at M06-2X, 6-311G(d,p) level (hartree)
REAL(8) :: EdftH=-0.498302d0                            ! DFT Energy (ZPE Corrected) at M06-2X, 6-311G(d,p) level (hartree) for Atomic Hydrogen

END MODULE


