PROGRAM HydrogenationPAHs

USE thermo_constants
USE hydrogen_constants

IMPLICIT NONE

! Declaration of thermodynamical input data:
REAL(8), DIMENSION(:,:), ALLOCATABLE :: nu				! Wavenumber (cm^-1)
REAL(8), DIMENSION(:,:), ALLOCATABLE :: Ele				! alpha and beta electrons
REAL(8), DIMENSION(:), ALLOCATABLE :: molar				! molar mass (amu)
REAL(8), DIMENSION(:), ALLOCATABLE :: symm				! symmetry number
REAL(8), DIMENSION(:,:), ALLOCATABLE :: Trot			! rotational temperature (K)

! Declaration of general input parameters:
INTEGER(4) :: Ni, n										! Number of atoms of pristine molecule and number of hydrogenations
INTEGER(4) :: tmod, pmod								! Parameters defining whether a T or P analysis must be performed (1/0)
integer(4) :: walkmode									! Parameter defining whether a Random Walk must be performed (1/0)
INTEGER(4) :: nT, nP									! Number of temperatures, pressures and fractional hydrogenation level to sample                           
INTEGER(4) ::  Fmin, Fmax, F							! Min, Max and intermediate vibrational degrees of freedom
REAL(8) :: Tin, Tfin, dT								! Initial and final temperatures and dT (K)
REAL(8) :: Pin, Pfin, dP								! Initial and final pressures and dP (K)
REAL(8) :: Th											! Temperature of H and H2 in random walk mode
REAL(8) :: pH2, pH										! Partial pressure of molecular and atomic hydrogen (bar).
														! NB: in T analysis mode (tmod=1) pH,pH2 are set equal to Pin
REAL(8) :: pCn											! Partial pressure of PAHs. 
														! NB: in T analysis mode (tmod=1), pCn is set equal to Pin for each PAH
REAL(8), DIMENSION(:), ALLOCATABLE :: T					! Temperature array
REAL(8), DIMENSION(:), ALLOCATABLE :: P					! Pressure array
REAL(8), DIMENSION(:), ALLOCATABLE :: E					! DFT Energy ZPE Corrected (eV) array

! Declaration of required input variables for Random Walk (Mixture Free Energy Optimization) mode:
INTEGER(4) :: seq										! Parameter defining whether you want to perform a sequential Opt or not (1/0)
INTEGER(4) :: npoint									! Number of points used in Lorentzian fitting
INTEGER(4) :: maxcounter, maxwrong						! Maxcounter: maximum number of opt steps at a given hydrogenation level before variation interval redefintion
														! Maxwrong: maximum number of opt attempts to fail before moving to a next hydrogenation level
REAL(8) :: incr											! Fraction of the most abundant specie that is added to the "forming" one 
REAL(8) :: fact1, threshdelta							! Fact1: factor that allows to re-size the variation interval during the opt
														! Threshdelta: threshold for re-sizing the variation interval
REAL(8) :: Y1,  W1										! Y1: factor that defines the extent of the variation interval at the outset
														! W1 is an internal parameter, not read in the input
REAL(8) :: threshold									! Threshold for final convergency 
REAL(8) :: W											! FWHM of the Lorentzian function 
REAL(8), DIMENSION(:), ALLOCATABLE :: Ncinp, Nhinp		! Number of bare PAHs and H2 molecules in the input configuration
REAL(8), DIMENSION(1:2) :: Ghyd							! Gibbs Free Energy array (H/H2 molecules)
REAL(8) :: ihydro										! Population of intermediate H species
REAL(8) :: Ntot											! Total number of molecules in the mixture
REAL(8), DIMENSION(:), ALLOCATABLE :: xctrial, xhtrial  ! Trial molar fractions of PAH and hydrogen molecules
REAL(8), DIMENSION(:), ALLOCATABLE :: Nctrial, Nhtrial  ! Trial number of PAH and hydrogen molecules
REAL(8) :: G0, Gold, Gmixn, Gmixtot                     ! Initial Free Energy, Free Energy to check, 
														! Mixture Free Energy of PAH component, Total Mixture Free Energy (eV)
REAL(8), DIMENSION(:), ALLOCATABLE :: Gn                ! PAH Free Energy array 
REAL(8) :: deltaGmix									! Computed Free energy difference: Gmixtot-G0
INTEGER(4) :: mostabun									! "n" (hydrogenation level) of Most abundant specie
REAL(8) :: Nctot0, Nhtot0, sumH							! Internal variables of the Random Walk 

! Declaration of output variables of Temperature and Pressure Analysis:
REAL(8) :: Zvib, Zrot, Zint, dbp						! Vib, Rot, Int partition functions of PAHs, DeBroglie Pressure (bar)
REAL(8) :: alpha, beta, gamma, G						! Contributing terms to the Free Energy of PAH and total Free Energy (eV)
REAL(8) :: ZrH2, ZvH2, ZiH2, lambdaH2, GH2				! Rot, Vib, Int partition functions, Thermal Wavelenght and Free energy (eV) of H2
REAL(8) :: lambdaH, Gah									! Thermal Wavelenght (m) and Free Energy of H (eV)

! Declaration of phase diagram variables:
INTEGER(4) :: phasemod									! Parameter defining whether a phase diagram must be computed or not (1/0)
REAL(8), DIMENSION(:), ALLOCATABLE :: fr				! Fractional hydrogenation level array
REAL(8), PARAMETER :: fi=0.0d0, ff=1.0d0				! Internal parameters
REAL(8) :: df											! Delta for fractional hydrogenation levels 
INTEGER(4) :: nf										! Number of fractional hydrogenation levels to check 
REAL(8), dimension(:), allocatable :: DGf  				! Formation Free energy difference 
REAL(8) :: dbp0, Z0										! De Broglie pressure and Internal Partition Functions of the bare PAH
REAL(8) :: p0											! Partial pressure of the bare PAH
REAL(8) :: minen										! Lowest Free Energy
INTEGER(4) :: minlabel									! Label (n) of the most stable specie

! Declaration of counters:
INTEGER(4) :: counter, wrong     						! Internal counter used in sequential Random Walk mode
INTEGER(4) :: i, m, s, l, k, j, label, iter
INTEGER(4) :: nright, nwrong, nneg, ntotal

! Declaration of filename:
CHARACTER(len=20) :: filename							! File containing relevant thermodynamical data
CHARACTER(len=20) :: zvibname, zrotname, zintname		! Output filename for partition functions
CHARACTER(len=20) :: dbpname							! Output filename for De Broglie Pressure
CHARACTER(len=20) :: alphaname, betaname, gammaname		! Output filename where contributing terms to the Free Energy are written
CHARACTER(len=20) :: freename							! Output filename for Free Energy
CHARACTER(len=20) :: optname							! Output filename where optimized compositions of sequential sub-random walk are written 

OPEN(unit=3333, File='history-random.log', Position='Append')

! Initialization of counter:
m=0
s=0
l=0
k=0
j=0
label=0
iter=0

READ(*,*) Ni, n
READ(*,*) tmod
READ(*,*) pmod
READ(*,*) walkmode
READ(*,*) phasemod
READ(*,*) Tin, Tfin, nT
READ(*,*) Pin, Pfin, nP
! Read H/H2 temperature:
READ(*,*) Th
! Read number of fractional hydrogenation level to sample in Phase Diagram mode:		
READ(*,*) nf 			


ALLOCATE(Ncinp(0:n))
ALLOCATE(Nhinp(1:2))

READ(*,*) Ncinp(0)
READ(*,*) Ncinp(1)
READ(*,*) Ncinp(n)
! Read Population of intermediate hydrogenated species:
READ(*,*) ihydro		
DO l=2,n-1
	Ncinp(l)= ihydro
END DO
! Read number of H atoms:
READ(*,*) Nhinp(1)		
! Read number of H2 molecules:
READ(*,*)  Nhinp(2)		
READ(*,*) incr
READ(*,*) Y1			
READ(*,*) fact1
READ(*,*) threshdelta, maxcounter
READ(*,*) maxwrong
READ(*,*) threshold
READ(*,*) seq			
READ(*,*) W, npoint		

! Computing min and max number of vibrational degrees of freedom:
Fmin=3*Ni-6                           
Fmax=3*(Ni+n)-6                       
WRITE(*,*) ' '
WRITE(*,*) 'Min Vibrational Degrees of Freedom =',  Fmin, 'Max Vibrational Degrees of Freedom = ', Fmax

ALLOCATE(nu(0:n,Fmax))
ALLOCATE(E(0:n))
ALLOCATE(Ele(0:n,2))
ALLOCATE(molar(0:n))
ALLOCATE(Trot(0:n,3))
ALLOCATE(symm(0:n))

! Reading data from thermo.dat files:
WRITE(*,*) ''
WRITE(*,*) 'Reading thermodynamical data from thermo.dat... '
DO s=0,n
	WRITE(filename,'(I2.2,"-thermo.dat")') s
	WRITE(*,*) filename
	OPEN(Unit=1, File=filename)
	READ(1,*) E(s)
	READ(1,*) (ele(s,j), j=1,2)
	READ(1,*) molar(s)
	READ(1,*) symm(s)
	READ(1,*) (Trot(s,j), j=1,3)
	F=3*(Ni+s)-6
	READ(1,*) (nu(s,j), j=1,F)
	CLOSE(1)
END DO

IF ( tmod.eq.1) then
	! Perform a Temperature Analysis: This mod provides Rotational, Vibrational 
 	! and Internal Partition Function as functions of temperature. The analysis 
 	! is performed also on hydrogen atom and hydrogen molecule. 

	! Generating Temperature array:
	ALLOCATE(T(0:nT-1))
	DO l=0,nT-1
	dT = (Tfin-Tin)/(nT-1)
		T(l)=dT*l+Tin
		WRITE(*,*) 'Temperature = ', T(l)
	END DO

	! Starting analysis on PAHs species:

	! In this mod, pressure is held fixed, hence every PAHs has the same pressure, set equal to Pin
	pCn=Pin							
	DO label=0,n
		F=3*(Ni+label)-6
		WRITE(zvibname,'(I2.2,"-Zvib.t")') label
		WRITE(zrotname, '(I2.2, "-Zrot.t")') label
		WRITE(zintname, '(I2.2,"-Zint.t")') label
		WRITE(dbpname, '(I2.2,"-dbp.t")') label
		WRITE(alphaname, '(I2.2,"-RTlnP.t")') label
		WRITE(betaname, '(I2.2,"-RTlnZ.t")' ) label
		WRITE(gammaname, '(I2.2,"-diffRT.t")' ) label
		WRITE(freename, '(I2.2, "-G.t")' ) label
		OPEN(Unit=2, File=zvibname, Position='Append')
		OPEN(Unit=3, File=zrotname, Position='Append')
		OPEN(Unit=4, File=zintname, Position='Append')
		OPEN(Unit=5, File=dbpname, Position='Append')
		OPEN(Unit=6, File=alphaname, Position='Append')
		OPEN(Unit=7, File=betaname, Position='Append')
		OPEN(Unit=8, File=gammaname, Position='Append')
		OPEN(Unit=9, File=freename, Position='Append')
		DO k=0,nT-1
			! Initialization 
			G=0.0d0
			CALL partition( T(k), nu, n, label, F, Fmax, ele, symm, Trot, Zvib, Zrot, Zint)
			CALL dbpressure( T(k), molar, n, label, dbp)
			WRITE(2,*) T(k), Zvib
			WRITE(3,*) T(k), Zrot
			WRITE(4,*) T(k), Zint
			WRITE(5,*) T(k), dbp			!(bar)
			! Computing contributions to free energy and free energy at fixed pressure
			! alpha = RTln(pCn/dbp), beta=RTln(Zn) , gamma=alpha-beta   (eV)
			alpha= ( R*T(k)*log(pCn/dbp) )*JmoltoeV
			beta=( R*T(k)*log(Zint) )*JmoltoeV
			gamma=alpha-beta
			G=E(label)*hartreetoJmol*JmoltoeV+gamma
			WRITE(6,*) T(k), alpha
			WRITE(7,*) T(k), beta
			WRITE(8,*) T(k), gamma
			WRITE(9,*) T(k), G
		END DO
		CLOSE(2)
		CLOSE(3)
		CLOSE(4)
		CLOSE(5)
		CLOSE(6)
		CLOSE(7)
		CLOSE(8)
		CLOSE(9)
	END DO
	! Starting analysis on H/H2:
	pH2=Pin
	pH=Pin
	OPEN(Unit=100, File='h2-G.t', Position='Append')
	OPEN(Unit=101, File='h2-lambda.t', Position='Append')
	OPEN(Unit=102, File='h2-Zvib.t', Position='Append')
	OPEN(Unit=103, File='h2-Zrot.t', Position='Append')
	OPEN(Unit=104, File='h2-Zint.t', Position='Append')
	OPEN(Unit=105, File='h-lambda.t', Position='Append')
	OPEN(Unit=106, File='h-G.t', Position='Append')
	DO k=0,nT-1
		CALL hydrogenfree( T(k), pH2, lambdaH2, ZvH2, ZrH2, ZiH2, GH2)
		CALL hatomfree(T(k), pH, lambdaH, Gah)
		WRITE(100,*) T(k), GH2
		WRITE(101,*) T(k), lambdaH2
		WRITE(102,*) T(k), ZvH2
		WRITE(103,*) T(k), ZrH2
		WRITE(104,*) T(k), ZiH2
		WRITE(105,*) T(k), lambdaH
		WRITE(106,*) T(k), Gah
	END DO
	CLOSE(100)
	CLOSE(101)
	CLOSE(102)
	CLOSE(103)
	CLOSE(104)
	CLOSE(105)
	CLOSE(106)
	DEALLOCATE(T)
ELSE IF (tmod.eq.0) then
	WRITE(*,*) ' No Temperature Analysis requested '
ELSE 
	WRITE(*,*) ' Warning: Invalid value for tmod! '
END IF

IF ( pmod.eq.1) then
	! Perform a Pressure Analysis: this mod provides dependence of 
	! free energies with respect to pressure at fixed temperature.
	! In this mod, T is set equal to Ti
	ALLOCATE(P(0:nP-1))
        DO l=0,nP-1
        dP = (Pfin-Pin)/(nP-1)
                P(l)=dP*l+Pin
                WRITE(*,*) 'Pressure = ', P(l)
        END DO
	! Starting analysis on PAHs species:
	DO label=0,n
		F=3*(Ni+label)-6
		CALL partition( Tin, nu, n, label, F, Fmax, ele, symm, Trot, Zvib, Zrot, Zint)
        CALL dbpressure( Tin, molar, n, label, dbp)
		! Computing pressure-independent terms:
		beta = ( R*Tin*log(Zint) )*JmoltoeV
		WRITE(alphaname, '(I2.2,"-RTlnP.p")') label
		WRITE(gammaname, '(I2.2, "-diffP.p")') label
		WRITE(freename, '(I2.2, "-G.p")' ) label
		OPEN(Unit=10, File=alphaname, Position='Append')
		OPEN(Unit=11, File=gammaname, Position='Append')
		OPEN(Unit=12, File=freename, Position='Append')
		DO k=0,nP-1
			! Initialization
			G=0.0d0
			!Computing contributions to free energy and free energy at fixed temperature
            !alpha = RTln(pCn/dbp), beta=RTln(Zn) , gamma=alpha-beta   (eV)
            alpha= ( R*Tin*log(P(k)/dbp) )*JmoltoeV
            gamma=alpha-beta
            G=E(label)*hartreetoJmol*JmoltoeV+gamma
			WRITE(10,*) P(k), alpha
			WRITE(11,*) P(k), gamma
            WRITE(12,*) P(k), G
        END DO
		CLOSE(10)
		CLOSE(11)
		CLOSE(12)
	END DO
	! Starting analysis on H/H2:
	OPEN(Unit=107, File='h2-G.p', Position='Append')
    OPEN(Unit=108, File='h-G.p', Position='Append')
    DO k=0,nP-1
        CALL hydrogenfree( Tin, P(k), lambdaH2, ZvH2, ZrH2, ZiH2, GH2)
        CALL hatomfree(Tin, P(k), lambdaH, Gah)
		WRITE(107,*) P(k), GH2
        WRITE(108,*) P(k), Gah
    END DO
	CLOSE(107)
	CLOSE(108)
	DEALLOCATE(P)
ELSE IF (pmod.eq.0) THEN
    WRITE(*,*) ' No Pressure Analysis requested '
ELSE 
    WRITE(*,*) ' Warning: Invalid value for pmod! '
END IF

! Random Walk at initial Tin And Pin:
IF (walkmode.eq.1) THEN
    WRITE(*,*) ' '
    WRITE(*,*) ' Input Mixture Composition:'
    WRITE(*,*) ' N 0H     = ', Ncinp(0),  'N iH     = ', Ncinp(1), '      N', n,'H =', Ncinp(n)
    WRITE(*,*) ' N H atom = ', Nhinp(1),  'N H2 mol = ', Nhinp(2)
    WRITE(*,*) ' '
	ALLOCATE(Nctrial(0:n))
	ALLOCATE(Nhtrial(1:2))
 	ALLOCATE(xctrial(0:n))
	ALLOCATE(xhtrial(1:2))
	ALLOCATE(Gn(0:n))
	! Computing total number of PAH molecules:
	Nctot0=0.0d0
	sumH=0.0d0
	DO label=0,n
		Nctot0=Nctot0+Ncinp(label)
		sumH=sumH+label*Ncinp(label)
	END DO
	! Computing total number of molecules in the mixture:
	Nhtot0=Nhinp(1)+2.0d0*Nhinp(2)+sumH
	
	Gold = 0.0d0        !Check on G
    ! ---------------- guess and initial G of the mixture ----------------!
    WRITE(*,*) 'Calculating Free Energy of the input mixture at given conditions...'
    Nhtrial = Nhinp
    Nctrial = Ncinp
    Ntot = sum(Ncinp)
    Ntot = Ntot + sum(Nhinp)
	WRITE(*,*) ' Total Number of Molecules in the mixture:', Ntot
	! Computing molar fractions:
    xhtrial = Nhinp / Ntot
    xctrial = Ncinp / Ntot
    pH2 = xhtrial(2)*Pin
    pH = xhtrial(1)*Pin
	! Computing H/H2 Gibbs Free Energy calling external subroutines:
    CALL hydrogenfree(Th, pH2, lambdaH2, ZvH2, ZrH2, ZiH2, Ghyd(2))	! H2 Gibbs free energy
    CALL hatomfree(Th, pH, lambdaH, Ghyd(1))          			! H  Gibbs free energy

    Gmixn = 0.0d0
    ! Cycle over hydrogenated species:
    DO label=0,n
    	F=3*(Ni+label)-6                                                                      
        CALL partition(Tin, nu, n, label, F, Fmax, ele, symm, Trot, Zvib, Zrot, Zint)                 
        CALL dbpressure(Tin, molar, n, label, dbp)
		! Computing partial Pressures in the mixture for each PAH specie:
        pCn=xctrial(label)*Pin                                                                  
		IF ( pCn < 1d-50 ) THEN
			Gn(label) = 0.0d0 
	   	ELSE
			Gn(label) = (E(label)*hartreetoJmol + R*Tin*log(pCn/dbp) - R*Tin*log(Zint))*JmoltoeV  
        END IF  
		! Computing of Free Energy for the "carbon" component:
		Gmixn=Gmixn+Nctrial(label)*Gn(label)                                                     
    END DO
	! Adding H/H2 contributions:
    G0=Gmixn+Nhtrial(1)*Ghyd(1)+Nhtrial(2)*Ghyd(2)                                             
    Gold = G0
    WRITE(*,*) '-------------'
    WRITE(*,*) ' Temperature (K) = ', Tin,' Pressure (bar) = ', Pin,' Gmixture =', G0
    ! ----------------------------------------------------! 
	
	! Sequential Random Walk:
	s=1
	iter=1
	counter=1
	W1=Y1
	nright=0
	nneg=0
	nwrong=0
	WRITE(3333,*) ' RANDOM WALK: '
	WRITE(3333,*) ' '
	WRITE(3333,*) ' Variation interval defined by:'
	WRITE(3333,*) 'Y1 = ', W1
	WRITE(3333,*) ' '
	WRITE(3333,*) ' RANDOM WALK - START'
	walk: DO WHILE (s <= n)
		ntotal=ntotal+1
		OPEN(unit=6666, File='history-steps.log', Action='write', Status='Replace')
		WRITE(6666,*) ' ntotal =', ntotal
		WRITE(6666,*) ' nright = ', nright, ' nneg = ', nneg, ' nwrong = ', nwrong
		WRITE(6666,*) ' %neg = ', (real(nneg)/real(ntotal))*100.0d0 , '%wrong = ', (real(nwrong)/real(ntotal))*100.0d0
		IF (counter > maxcounter ) then
			WRITE(*,*) ' Warning: VARIATION INTERVAL set to half the original '
			WRITE(3333,*) ' Counter exceed maxcounter, counter = ', counter, ' at iter = ', iter
			WRITE(3333,*) ' VARIATION INTERVAL modified by a factor = ', fact1
			WRITE(3333,*) ' Y1 = ', W1
			WRITE(3333,*) ' '
			W1=W1*fact1
			counter=1
		ELSE IF ( counter == 0 ) then
			W1=Y1
		END IF
		 
		IF (wrong > maxwrong ) then
			WRITE(*,*) ' Warning: Not converged after ', wrong, ' iterations! '
			wrong=1
			WRITE(optname, '(I2.2,"-seqopt.out")') s
			OPEN(Unit=1100, File=optname)
			DO m=0,n
				WRITE(1100,*) m, Ncinp(m), Tin, Pin
			END DO
			WRITE(1100,*) "-1", Nhinp(1), Tin, Pin
			WRITE(1100,*) "-2", Nhinp(2), Tin, Pin
			CLOSE(1100)
			CALL Lorentzian(s, W, npoint, seq, n)
			IF ( s == n ) then 
				STOP
			ELSE 
				s=s+1  
				! When a new hydrogenation starts, the number of molecules of the "forming" specie is 
				! redefined as the 10% of the most abundant specie currently in the mixture. 
				! Doing so, the "variation space" for the forming species is not vanishing. 
				! This choice should not have a great impact on the currently reached equilibrium. 
				Ncinp(s)=Ncinp(s)+MAXVAL(Ncinp)*0.01d0
				mostabun=MAXLOC(Ncinp, DIM=1) - 1 
				Ncinp(mostabun)=Ncinp(mostabun)-MAXVAL(Ncinp)*0.01d0 
			END IF                         
			counter=0
			WRITE(*,*) s
		END IF
		
		CALL RandomWalk(W1, s, n, Nctot0, Nhtot0, Ncinp, Nhinp, Nctrial, Nhtrial, xctrial, xhtrial)
		pH2=xhtrial(2)*Pin
		pH=xhtrial(1)*Pin
		CALL hydrogenfree(Th, pH2, lambdaH2, ZvH2, ZrH2, ZiH2, Ghyd(2))        ! H2 Gibbs free energy
		CALL hatomfree(Th, pH, lambdaH, Ghyd(1))                               ! H  Gibbs free energy
		! Cycle over PAH species: 
		Gmixn = 0.0d0
		DO label=0,n
			F=3*(Ni+label)-6 
			CALL partition(Tin, nu, n, label, F, Fmax, ele, symm, Trot, Zvib, Zrot, Zint)
			CALL dbpressure(Tin, molar, n, label, dbp)
			pCn=xctrial(label)*Pin
			IF (pCn == 0.0d0 ) THEN
				Gn(label) = 0.0d0
			ELSE
				Gn(label) = (E(label)*hartreetoJmol + R*Tin*log(pCn/dbp) - R*Tin*log(Zint))*JmoltoeV
			END IF
				Gmixn=Gmixn+Nctrial(label)*Gn(label)                                                     
			END DO
			Gmixtot=Gmixn+Nhtrial(1)*Ghyd(1)+Nhtrial(2)*Ghyd(2)
			IF (Gmixtot < Gold) THEN
				wrong=1
				nright=nright+1
				OPEN(Unit=1000, File='G-RandomWalk.out', Position='Append')
				WRITE(1000,*) iter, Gmixtot-G0
				CLOSE(1000)
				! Counter of Total Iterations in RW:
				iter=iter+1					
				Ncinp=Nctrial
                Nhinp=Nhtrial
				deltaGmix=abs(Gmixtot-Gold)
				IF (deltaGmix <= threshdelta ) THEN
					counter=counter+1
				ELSE 
					counter=1
				END IF
				WRITE(*,*) iter-1, deltaGmix
                Gold=Gmixtot
				IF (deltaGmix < threshold ) THEN
					WRITE(3333,*) ' Optimized mixture found at iter =', iter-1
					WRITE(3333,*) ' '
					WRITE(3333,*) ' Adding specie with n = ', s+1
					WRITE(optname, '(I2.2,"-seqopt.out")') s
					OPEN(Unit=1100, File=optname)
					DO m=0,n
						WRITE(1100,*) m, Nctrial(m), xctrial(m), Tin, Pin
					END DO
					WRITE(1100,*) "-1", Nhtrial(1), xhtrial(1), Tin, Pin
					WRITE(1100,*) "-2", Nhtrial(2), xhtrial(2), Tin, Pin
					CLOSE(1100)
					CALL Lorentzian(s, W, npoint, seq, n)
					s=s+1
					! When a new hydrogenation starts, the number of molecules of the "forming" specie is 
					! redefined as the 10% of the most abundant specie currently in the mixture. 
					! Doing so, the "variation space" for the forming species is not vanishing. 
					! This choice should not have a great impact on the currently reached equilibrium.
					IF (s<= n) THEN
                        Ncinp(s)=Ncinp(s)+MAXVAL(Ncinp)*incr
                        mostabun=MAXLOC(Ncinp, DIM=1) -1 
                        Ncinp(mostabun)=Ncinp(mostabun)-MAXVAL(Ncinp)*incr	
					END IF 			
					counter=0
					WRITE(*,*) s
				END IF
			ELSE
				nwrong=nwrong+1
				wrong=wrong+1
			END IF
		CLOSE(6666)
        END DO walk
	WRITE(3333,*) ' Random Walk finished at iter =', iter
	CLOSE(3333)
	DEALLOCATE(Nctrial, Nhtrial, xctrial, xhtrial, Gn)
END IF

IF (phasemod.eq.1) THEN
! PHASE DIAGRAM Mode:
! Computing T, P and f grid:
	ALLOCATE(fr(0:nf-1))
    DO l=0,nf-1
    	df = (ff-fi)/(nf-1)
        fr(l)=df*l+fi
        WRITE(*,*) 'Fractional hydrogenation = ', fr(l)
    END DO

	ALLOCATE(T(0:nT-1))
    DO l=0,nT-1
        dT = (Tfin-Tin)/(nT-1)
        T(l)=dT*l+Tin
        WRITE(*,*) 'Temperature = ', T(l)
    END DO

	ALLOCATE(P(0:nP-1))
    DO l=0,nP-1
        dP = (Pfin-Pin)/(nP-1)
        P(l)=dP*l+Pin
        WRITE(*,*) 'Pressure = ', P(l)
    END DO


	WRITE(*,*) ''
    WRITE(*,*) '-------------'
    WRITE(*,*) 'PHASE DIAGRAM:'
    WRITE(*,*) ' '
    ! Starting simulation:
    ! Initiliazing Gn:
    ALLOCATE(Gn(0:n))
    Gn=0.0d0
    Ghyd=0.0d0
    ALLOCATE(DGf(1:n))
    OPEN(UNIT=500, FILE='Phase-diagram.dat', ACTION='Write', POSITION='Append')
    WRITE(500,*) 'T (kelvin)            P (bar)                 f               Label'
    WRITE(*,*) 'Calculating Free energy differences...'

	! Cycle over T:
	DO k=0,nT-1
		F=3*(Ni+0)-6
		! Computing De Broglie pressure and Partition function of bare molecules:
        CALL partition(T(k), nu, n, 0, F, Fmax, ele, symm, Trot, Zvib, Zrot, Z0)
		CALL dbpressure(T(k), molar, n, 0, dbp0)    							
		!Cycle over P
        DO s=0,nP-1
        	! Cycle over f:
            DO i=1, nf-2
                pH2=fr(i)*P(s)  			
				! Computing partial pressure of bare molecules:						
                p0=P(s)-pH2     		
				! Computing G at given thermodynamical conditions							
                Gn(0) = ( E(0)*hartreetoJmol + R*T(k)*log(p0/dbp0) - R*T(k)*log(Z0))*JmoltoeV           
				CALL hydrogenfree(T(k), pH2, lambdaH2, ZvH2, ZrH2, ZiH2, Ghyd(2))
                ! Cycle over hydrogenated molecules:
                DO label=1,n
                    F=3*(Ni+label)-6        						
                    CALL partition(T(k), nu, n, label, F, Fmax, ele, symm, Trot, Zvib, Zrot, Zint)    
                    CALL dbpressure(T(k),  molar, n, label, dbp)
                    Gn(label) = (E(label)*hartreetoJmol + R*T(k)*log(P(s)/dbp) - R*T(k)*log(Zint))*JmoltoeV   
					! Computing DeltaG:
                    DGf(label) = Gn(label)-Gn(0)-(label/2.0d0)*Ghyd(2) 
                END DO
					! Looking for smallest values at given conditions:
                    minlabel=MINLOC(DGf, DIM=1) 
                    minen=MINVAL(DGf)
				WRITE(500,*) T(k), P(s), fr(i), "n = ", minlabel, minen
	        END DO
        END DO
	END DO

	WRITE(*,*) 'Phase diagram done!!'

	DEALLOCATE(DGf)
	DEALLOCATE(Gn)
	DEALLOCATE(T,P)
END IF

END PROGRAM
