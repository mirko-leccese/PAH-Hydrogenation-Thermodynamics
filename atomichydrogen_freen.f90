SUBROUTINE hatomfree(T, pH, lambdaH, Gah)

USE hydrogen_constants
USE thermo_constants

!Computing free energy of atomic hydrogen

REAL(8) :: T                    !Temperature
REAL(8) :: pH                   !Pressure
REAL(8) :: lambdaH, dbpH, Gah   !Thermal Wavelenght, DeBroglie pressure, Free energy
REAL(8) :: mH                   !Atomic mass of H

mH=molarh*amutog*gtokg
lambdaH = h/sqrt(2.0d0*pi*mH*kb*T)
dbpH = (kb*T/lambdaH**3.0d0)*Patobar
Gah = (EdftH*hartreetoJmol + R*T*log(pH/dbpH))*JmoltoeV

END SUBROUTINE

