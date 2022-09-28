SUBROUTINE hydrogenfree(T, pH2, lambdaH2, ZvH2, ZrH2, ZiH2, GH2)

USE hydrogen_constants
USE thermo_constants 

IMPLICIT NONE

!Declaring input variables
REAL(8) :: T
REAL(8) :: pH2

!Declaring side and output variables
REAL(8) :: Tred, Gt
REAL(8) :: ex
REAL(8) :: ZtH2, ZvH2, ZrH2, ZiH2
REAL(8) :: lambdaH2, mH2, dbpH2
REAL(8) :: HH2, SH2, GH2

!Initialization
ZvH2=0.0d0
ZrH2=0.0d0
ZiH2=0.0d0
lambdaH2=0.0d0

IF (T < 298.15d0) THEN
        !Computing G from partition function
         ex= -h*c*nuh2/(kb*T)
         ZvH2=1.0d0/(1.0d0-exp(ex))
         ZrH2=(1.0d0/2.0d0)*(kb*T)/(h*c*B)
         mH2 = molarh2*amutog*gtokg
         lambdaH2 = h/sqrt(2.0d0*pi*mH2*kb*T)
         dbpH2 = ( kb*T/lambdaH2**3.0d0 )*Patobar
         ZiH2=ZrH2*ZvH2
         GH2 = (Edft*hartreetoJmol + R*T*log(pH2/dbpH2) - R*T*log(ZiH2))*JmoltoeV
ELSE IF ( T >= 298.15d0 .AND. T <= 1000.0d0 ) THEN
	!For high temperatures, we do not use harmonic oscillator and rigid rotator 
	!approximation, but the Shomate equation. Hence we fix
	lambdaH2=0.0d0
	ZvH2=0.0d0
	ZrH2=0.0d0
	ZiH2=0.0d0
        !Computing enthalpy
        HH2 = (Edft*hartreetoJmol) + (7.0d0/2.0d0)*R*T
        !Computing entrophy
        Tred=T/1.0d+3
        SH2 = A1*log(Tred) + B1*Tred + C1*Tred**2/2.0d0 + D1*Tred**3/3.0d0 - E1/(2.0d0*Tred**2) + G1
       !Computing pressure-independent G
        Gt = HH2 - T*SH2

        !Computing Free Energy at T,pH2

        GH2 = ( Gt + R*T*log(pH2/pstan) )*JmoltoeV
ELSE IF ( T > 1000.0d0 .AND. T <= 2500.0d0 ) THEN
        !For high temperatures, we do not use harmonic oscillator and rigid rotator 
        !approximation, but the Shomate equation. Hence we fix
        lambdaH2=0.0d0
        ZvH2=0.0d0
        ZrH2=0.0d0
            !Computing enthalpy
        HH2 = (Edft*hartreetoJmol) + (7.0d0/2.0d0)*R*T
        !Computing entrophy
        Tred=T/1.0d+3
        SH2 = A2*log(Tred) + B2*Tred + C2*Tred**2/2.0d0 + D2*Tred**3/3.0d0 - E2/(2.0d0*Tred**2) + G2
       !Computing pressure-independent G
        Gt = HH2 - T*SH2

        !Computing Free Energy at T,pH2

        GH2 = ( Gt + R*T*log(pH2/pstan) )*JmoltoeV
ELSE
	STOP
END IF

END SUBROUTINE
