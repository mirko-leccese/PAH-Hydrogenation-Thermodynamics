MODULE thermo_constants

!Declaration of parameters
REAL(8), PARAMETER :: h=6.62607015d-34          !Plank constant (J/s)
REAL(8), PARAMETER :: c=2.99792458d+8           !Speed of light (m/s)
REAL(8), PARAMETER :: kb=1.3806485279d-23       !Boltzmann constant (J/k)
REAL(8), PARAMETER :: amutog=1.6603145d-24      !Amu to g
REAL(8), PARAMETER :: pi=3.14159d0              !pi greek
REAL(8), PARAMETER :: onethird=1.0d0/3.0d0
REAL(8), PARAMETER :: threesecond=3.0d0/2.0d0
REAL(8), PARAMETER :: gtokg=1.0d-3              !From g to Kg
REAL(8), PARAMETER :: Patobar=1.0d-5            !From Pascal to bar
REAL(8), PARAMETER :: hartreetoeV=27.211399d0   !Hartree to eV
!Useful constants
REAL(8), PARAMETER :: Vstan=0.0224d0            !Standard volume (m^3)
REAL(8) :: pstan=1.0d0                   	!Standard pressure (bar)
REAL(8) :: Tstan=298.15d0                	!Standard temperature (k)
REAL(8) :: R=8.314459848d0               	!Gas constant (J/mol*K)
REAL(8) :: hartreetoJmol=2625.5002d+3		!Conversion from hartree to J/mol
REAL(8) :: JmoltoeV=1.03642723013d-5		!Conversion from J/mol to eV
REAL(8), PARAMETER :: Nav=6.02214085774d+23	!Avogadro number 


END MODULE






