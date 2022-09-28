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


!Initialization of counter
m=0
s=0
l=0
k=0
j=0
label=0
iter=0

read(*,*) Ni, n
read(*,*) tmod
read(*,*) pmod
read(*,*) walkmode
read(*,*) phasemod
read(*,*) Tin, Tfin, nT
read(*,*) Pin, Pfin, nP
read(*,*) Th			!Hydrogen and Molecular Hydrogen temperature
read(*,*) nf 			!Number of fractional hydrogenation level to sample in Phase Diagram mode


allocate(Ncinp(0:n))
allocate(Nhinp(1:2))

read(*,*) Ncinp(0)
read(*,*) Ncinp(1)
read(*,*) Ncinp(n)
read(*,*) ihydro		!Population of intermediate hydrogenated species
do l=2,n-1
	Ncinp(l)= ihydro
end do
read(*,*) Nhinp(1)		!Hydrogen Atoms
read(*,*)  Nhinp(2)		!Hydrogen Molecules
read(*,*) incr
read(*,*) Y1			!Factor that controls the interval of random variations
read(*,*) fact1
read(*,*) threshdelta, maxcounter
read(*,*) maxwrong
read(*,*) threshold
read(*,*) seq			!Check sequential random walk mode
read(*,*) W, npoint		!Width of lorenzian smearing and number of points

Fmin=3*Ni-6                           !Min vibrational degrees of freedom
Fmax=3*(Ni+n)-6                       !Max vibrational degrees of freedom
WRITE(*,*) ' '
WRITE(*,*) 'Min Vibrational Degrees of Freedom =',  Fmin, 'Max Vibrational Degrees of Freedom = ', Fmax

ALLOCATE(nu(0:n,Fmax))
ALLOCATE(E(0:n))
ALLOCATE(Ele(0:n,2))
ALLOCATE(molar(0:n))
ALLOCATE(Trot(0:n,3))
ALLOCATE(symm(0:n))

  !Reading data from thermo.dat
write(*,*) ''
write(*,*) 'Reading thermodynamical data from thermo.dat... '
do s=0,n
	write(filename,'(I2.2,"-thermo.dat")') s
	write(*,*) filename
	open(Unit=1, File=filename)
	read(1,*) E(s)
	read(1,*) (ele(s,j), j=1,2)
	read(1,*) molar(s)
	read(1,*) symm(s)
	read(1,*) (Trot(s,j), j=1,3)
	F=3*(Ni+s)-6
	read(1,*) (nu(s,j), j=1,F)
	close(1)
end do

if ( tmod.eq.1) then
 !Perform a Temperature Analysis: This mod provides Rotational, Vibrational 
 !and Internal Partition Function as functions of temperature. The analysis 
 !is performed also on hydrogen atom and hydrogen molecule. 
	allocate(T(0:nT-1))
	do l=0,nT-1
	dT = (Tfin-Tin)/(nT-1)
		T(l)=dT*l+Tin
		write(*,*) 'Temperature = ', T(l)
	end do
	!Analysis on PAHs species
	pCn=Pin							!In this mod, pressure is held fixed, hence every PAHs
								!has the same pressure, set equal to Pin
	do label=0,n
		F=3*(Ni+label)-6
		write(zvibname,'(I2.2,"-Zvib.t")') label
		write(zrotname, '(I2.2, "-Zrot.t")') label
		write(zintname, '(I2.2,"-Zint.t")') label
		write(dbpname, '(I2.2,"-dbp.t")') label
		write(alphaname, '(I2.2,"-RTlnP.t")') label
		write(betaname, '(I2.2,"-RTlnZ.t")' ) label
		write(gammaname, '(I2.2,"-diffRT.t")' ) label
		write(freename, '(I2.2, "-G.t")' ) label
		open(Unit=2, File=zvibname, Position='Append')
		open(Unit=3, File=zrotname, Position='Append')
		open(Unit=4, File=zintname, Position='Append')
		open(Unit=5, File=dbpname, Position='Append')
		open(Unit=6, File=alphaname, Position='Append')
		open(Unit=7, File=betaname, Position='Append')
		open(Unit=8, File=gammaname, Position='Append')
		open(Unit=9, File=freename, Position='Append')
		do k=0,nT-1
			!Initialization 
			G=0.0d0
			CALL partition( T(k), nu, n, label, F, Fmax, ele, symm, Trot, Zvib, Zrot, Zint)
			CALL dbpressure( T(k), molar, n, label, dbp)
			write(2,*) T(k), Zvib
			write(3,*) T(k), Zrot
			write(4,*) T(k), Zint
			write(5,*) T(k), dbp			!(bar)
			!Computing contributions to free energy and free energy at fixed Pressure
			!alpha = RTln(pCn/dbp), beta=RTln(Zn) , gamma=alpha-beta   in eV
			alpha= ( R*T(k)*log(pCn/dbp) )*JmoltoeV
			beta=( R*T(k)*log(Zint) )*JmoltoeV
			gamma=alpha-beta
			G=E(label)*hartreetoJmol*JmoltoeV+gamma
			write(6,*) T(k), alpha
			write(7,*) T(k), beta
			write(8,*) T(k), gamma
			write(9,*) T(k), G
		end do
		close(2)
		close(3)
		close(4)
		close(5)
		close(6)
		close(7)
		close(8)
		close(9)
	end do
	!Analysis on Hydrogen Atom and Molecule
	pH2=Pin
	pH=Pin
	open(Unit=100, File='h2-G.t', Position='Append')
	open(Unit=101, File='h2-lambda.t', Position='Append')
	open(Unit=102, File='h2-Zvib.t', Position='Append')
	open(Unit=103, File='h2-Zrot.t', Position='Append')
	open(Unit=104, File='h2-Zint.t', Position='Append')
	open(Unit=105, File='h-lambda.t', Position='Append')
	open(Unit=106, File='h-G.t', Position='Append')
	do k=0,nT-1
		CALL hydrogenfree( T(k), pH2, lambdaH2, ZvH2, ZrH2, ZiH2, GH2)
		CALL hatomfree(T(k), pH, lambdaH, Gah)
		write(100,*) T(k), GH2
		write(101,*) T(k), lambdaH2
		write(102,*) T(k), ZvH2
		write(103,*) T(k), ZrH2
		write(104,*) T(k), ZiH2
		write(105,*) T(k), lambdaH
		write(106,*) T(k), Gah
	end do
	close(100)
	close(101)
	close(102)
	close(103)
	close(104)
	close(105)
	close(106)
	deallocate(T)
else if (tmod.eq.0) then
	write(*,*) ' No Temperature Analysis requested '
else 
	write(*,*) ' Warning: Invalid value for tmod! '
end if

if ( pmod.eq.1) then
	 !Perform a Pressure Analysis: This mod provides dependence of 
	 !free energies with respect to pressure at fixed temperature.
	 !In this mod, T is set equal to Ti
	allocate(P(0:nP-1))
        do l=0,nP-1
        dP = (Pfin-Pin)/(nP-1)
                P(l)=dP*l+Pin
                write(*,*) 'Pressure = ', P(l)
        end do
	!Analysis on PAHs species
	do label=0,n
		F=3*(Ni+label)-6
		CALL partition( Tin, nu, n, label, F, Fmax, ele, symm, Trot, Zvib, Zrot, Zint)
        	CALL dbpressure( Tin, molar, n, label, dbp)
		!Computing term of Free energy independent from pressure
		beta = ( R*Tin*log(Zint) )*JmoltoeV
		write(alphaname, '(I2.2,"-RTlnP.p")') label
		write(gammaname, '(I2.2, "-diffP.p")') label
		write(freename, '(I2.2, "-G.p")' ) label
		open(Unit=10, File=alphaname, Position='Append')
		open(Unit=11, File=gammaname, Position='Append')
		open(Unit=12, File=freename, Position='Append')
		do k=0,nP-1
			!Initialization
			G=0.0d0
			!Computing contributions to free energy and free energy at fixed temperature
                        !alpha = RTln(pCn/dbp), beta=RTln(Zn) , gamma=alpha-beta   in eV
                        alpha= ( R*Tin*log(P(k)/dbp) )*JmoltoeV
                        gamma=alpha-beta
                        G=E(label)*hartreetoJmol*JmoltoeV+gamma
			write(10,*) P(k), alpha
			write(11,*) P(k), gamma
                        write(12,*) P(k), G
                end do
		close(10)
		close(11)
		close(12)
	end do
	!Analysis on Hydrogen Atom And Hydrogen Molecule
	open(Unit=107, File='h2-G.p', Position='Append')
        open(Unit=108, File='h-G.p', Position='Append')
        do k=0,nP-1
                CALL hydrogenfree( Tin, P(k), lambdaH2, ZvH2, ZrH2, ZiH2, GH2)
                CALL hatomfree(Tin, P(k), lambdaH, Gah)
		write(107,*) P(k), GH2
                write(108,*) P(k), Gah
        end do
	close(107)
	close(108)
	deallocate(P)
else if (pmod.eq.0) then
        write(*,*) ' No Pressure Analysis requested '
else 
        write(*,*) ' Warning: Invalid value for pmod! '
end if

!Random Walk at Tin And Pin
if (walkmode.eq.1) then
        write(*,*) ' '
        write(*,*) ' Input Mixture Composition:'
        write(*,*) ' N 0H     = ', Ncinp(0),  'N iH     = ', Ncinp(1), '      N', n,'H =', Ncinp(n)
        write(*,*) ' N H atom = ', Nhinp(1),  'N H2 mol = ', Nhinp(2)
        write(*,*) ' '
	allocate(Nctrial(0:n))
	allocate(Nhtrial(1:2))
 	allocate(xctrial(0:n))
	allocate(xhtrial(1:2))
	allocate(Gn(0:n))
	!Total number of pahs
	Nctot0=0.0d0
	sumH=0.0d0
	do label=0,n
		Nctot0=Nctot0+Ncinp(label)
		sumH=sumH+label*Ncinp(label)
	end do
	Nhtot0=Nhinp(1)+2.0d0*Nhinp(2)+sumH
	
	Gold = 0.0d0        !Check on G
        ! ---------------- guess and initial G of the mixture ----------------!
        write(*,*) 'Calculating Free Energy of the input mixture at given conditions...'
        Nhtrial = Nhinp
        Nctrial = Ncinp
        Ntot = sum(Ncinp)
        Ntot = Ntot + sum(Nhinp)
	write(*,*) ' Total Number of Molecules in the mixture:', Ntot
        xhtrial = Nhinp / Ntot
        xctrial = Ncinp / Ntot
        pH2 = xhtrial(2)*Pin
        pH = xhtrial(1)*Pin
        CALL hydrogenfree(Th, pH2, lambdaH2, ZvH2, ZrH2, ZiH2, Ghyd(2))	! H2 Gibbs free energy
        CALL hatomfree(Th, pH, lambdaH, Ghyd(1))          			! H  Gibbs free energy

        Gmixn = 0.0d0
        !Cycle over hydrogenated species
        do label=0,n
           F=3*(Ni+label)-6                                                                      
           CALL partition(Tin, nu, n, label, F, Fmax, ele, symm, Trot, Zvib, Zrot, Zint)                 
           CALL dbpressure(Tin, molar, n, label, dbp)
           pCn=xctrial(label)*Pin                                                                  ! Partial Pressure in the mixture of each specie
	   if ( pCn < 1d-50 ) then
		Gn(label) = 0.0d0 
	   else
		Gn(label) = (E(label)*hartreetoJmol + R*Tin*log(pCn/dbp) - R*Tin*log(Zint))*JmoltoeV  
           end if   
	Gmixn=Gmixn+Nctrial(label)*Gn(label)                                                    ! Computing of PAH mixture only 
        end do
        G0=Gmixn+Nhtrial(1)*Ghyd(1)+Nhtrial(2)*Ghyd(2)                                             !Adding hydrogen contributions
        Gold = G0
        write(*,*) '-------------'
        write(*,*) ' Temperature (K) = ', Tin,' Pressure (bar) = ', Pin,' Gmixture =', G0
        ! ----------------------------------------------------! 
	!Random Walk Mode

	!Sequential Random Walk
	s=1
	iter=1
	counter=1
	W1=Y1
	nright=0
	nneg=0
	nwrong=0
	write(3333,*) ' RANDOM WALK: '
	write(3333,*) ' '
	write(3333,*) ' Variation interval defined by:'
	write(3333,*) 'Y1 = ', W1
	write(3333,*) ' '
	write(3333,*) ' RANDOM WALK - START'
	walk: do while (s <= n)
		ntotal=ntotal+1
		open(unit=6666, File='history-steps.log', Action='Write', Status='Replace')
		write(6666,*) ' ntotal =', ntotal
		write(6666,*) ' nright = ', nright, ' nneg = ', nneg, ' nwrong = ', nwrong
		write(6666,*) ' %neg = ', (real(nneg)/real(ntotal))*100.0d0 , '%wrong = ', (real(nwrong)/real(ntotal))*100.0d0
		if (counter > maxcounter ) then
			write(*,*) ' Warning: VARIATION INTERVAL set to half the original '
			write(3333,*) ' Counter exceed maxcounter, counter = ', counter, ' at iter = ', iter
			write(3333,*) ' VARIATION INTERVAL modified by a factor = ', fact1
			write(3333,*) ' Y1 = ', W1
			write(3333,*) ' '
			W1=W1*fact1
			!W2=W2*fact2
			counter=1
		else if ( counter == 0 ) then
			W1=Y1
		end if
		 
		if (wrong > maxwrong ) then
			write(*,*) ' Warning: Not converged after ', wrong, ' iterations! '
			wrong=1
			write(optname, '(I2.2,"-seqopt.out")') s
			open(Unit=1100, File=optname)
			do m=0,n
				write(1100,*) m, Ncinp(m), Tin, Pin
			end do
			write(1100,*) "-1", Nhinp(1), Tin, Pin
			write(1100,*) "-2", Nhinp(2), Tin, Pin
			close(1100)
			CALL Lorentzian(s, W, npoint, seq, n)
			if ( s == n ) then 
				STOP
			else 
				s=s+1  
!When a new hydrogenation starts, the number of molecules of the incoming specie is 
!redefined as the 10% of the most abundant species in the mixture. In this way, the 
!variation space of this new specie is increased.
				Ncinp(s)=Ncinp(s)+MAXVAL(Ncinp)*0.01d0
				mostabun=MAXLOC(Ncinp, DIM=1) - 1 
				Ncinp(mostabun)=Ncinp(mostabun)-MAXVAL(Ncinp)*0.01d0 
			end if                         
			counter=0
			write(*,*) s
		end if
		
		 CALL RandomWalk(W1, s, n, Nctot0, Nhtot0, Ncinp, Nhinp, Nctrial, Nhtrial, xctrial, xhtrial)
		!if (Nctrial(s)<0.0d0.or.Nhtrial(1)<0.0d0) then 
		!	nneg=nneg+1
		!	close(6666)
		!	cycle walk
		!else
			pH2=xhtrial(2)*Pin
			pH=xhtrial(1)*Pin
			CALL hydrogenfree(Th, pH2, lambdaH2, ZvH2, ZrH2, ZiH2, Ghyd(2))        ! H2 Gibbs free energy
			CALL hatomfree(Th, pH, lambdaH, Ghyd(1))                               ! H  Gibbs free energy
			!Cycle over molecular species 
			Gmixn = 0.0d0
			do label=0,n
				F=3*(Ni+label)-6 
				CALL partition(Tin, nu, n, label, F, Fmax, ele, symm, Trot, Zvib, Zrot, Zint)
				CALL dbpressure(Tin, molar, n, label, dbp)
				pCn=xctrial(label)*Pin
				if (pCn == 0.0d0 ) then
					Gn(label) = 0.0d0
				else
					Gn(label) = (E(label)*hartreetoJmol + R*Tin*log(pCn/dbp) - R*Tin*log(Zint))*JmoltoeV
				end if
				Gmixn=Gmixn+Nctrial(label)*Gn(label)                                                    ! Computing of PAH mixture only 
			end do
			Gmixtot=Gmixn+Nhtrial(1)*Ghyd(1)+Nhtrial(2)*Ghyd(2)
			if (Gmixtot < Gold) then
				wrong=1
				nright=nright+1
				open(Unit=1000, File='G-RandomWalk.out', Position='Append')
				write(1000,*) iter, Gmixtot-G0
				close(1000)
				iter=iter+1					!Counter of Total Iterations in RW
				Ncinp=Nctrial
                 		Nhinp=Nhtrial
				deltaGmix=abs(Gmixtot-Gold)
				if (deltaGmix <= threshdelta ) then
					counter=counter+1
				else 
					counter=1
				end if
				write(*,*) iter-1, deltaGmix
                 		Gold=Gmixtot
				if (deltaGmix < threshold ) then
					write(3333,*) ' Optimized mixture found at iter =', iter-1
					write(3333,*) ' '
					write(3333,*) ' Adding specie with n = ', s+1
					write(optname, '(I2.2,"-seqopt.out")') s
					open(Unit=1100, File=optname)
					do m=0,n
						write(1100,*) m, Nctrial(m), xctrial(m), Tin, Pin
					end do
					write(1100,*) "-1", Nhtrial(1), xhtrial(1), Tin, Pin
					write(1100,*) "-2", Nhtrial(2), xhtrial(2), Tin, Pin
					close(1100)
					CALL Lorentzian(s, W, npoint, seq, n)
					s=s+1
!!When a new hydrogenation starts, the number of molecules of the incoming specie is 
!redefined as the 10% of the most abundant species in the mixture. In this way, the 
!variation space of this new specie is increased.
				if (s<= n) then
                                	Ncinp(s)=Ncinp(s)+MAXVAL(Ncinp)*incr
                                	mostabun=MAXLOC(Ncinp, DIM=1) -1 
                                	Ncinp(mostabun)=Ncinp(mostabun)-MAXVAL(Ncinp)*incr	
				end if 			
					counter=0
					write(*,*) s
				end if
			else
				nwrong=nwrong+1
				wrong=wrong+1
			end if
		close(6666)
        end do walk
	write(3333,*) ' Random Walk finished at iter =', iter
	close(3333)
	deallocate(Nctrial, Nhtrial, xctrial, xhtrial, Gn)
end if

if (phasemod.eq.1) then
        !Perfrom phase diagram 
	allocate(fr(0:nf-1))
        do l=0,nf-1
        df = (ff-fi)/(nf-1)
                fr(l)=df*l+fi
                write(*,*) 'Fractional hydrogenation = ', fr(l)
        end do

	allocate(T(0:nT-1))
        do l=0,nT-1
        	dT = (Tfin-Tin)/(nT-1)
                T(l)=dT*l+Tin
                write(*,*) 'Temperature = ', T(l)
        end do

	allocate(P(0:nP-1))
        do l=0,nP-1
        	dP = (Pfin-Pin)/(nP-1)
                P(l)=dP*l+Pin
                write(*,*) 'Pressure = ', P(l)
        end do


	    WRITE(*,*) ''
    write(*,*) '-------------'
    write(*,*) 'PHASE DIAGRAM:'
    write(*,*) ' '
    !Simulating phase diagram
    !Initiliazing Gn
    allocate(Gn(0:n))
    Gn=0.0d0
    Ghyd=0.0d0
    allocate(DGf(1:n))
    OPEN(UNIT=500, FILE='Phase-diagram.dat', ACTION='Write', POSITION='Append')
    write(500,*) 'T (kelvin)            P (bar)                 f               Label'
    write(*,*) 'Calculating Free energy differences...'

	do k=0,nT-1
		F=3*(Ni+0)-6
        	CALL partition(T(k), nu, n, 0, F, Fmax, ele, symm, Trot, Zvib, Zrot, Z0)
		CALL dbpressure(T(k), molar, n, 0, dbp0)    							!Computing debroglie pressure of pristine
		!Cycle over P
        	do s=0,nP-1
                	!Cycle over f
                	do i=1, nf-2
                		pH2=fr(i)*P(s)  									!Computing partial pressure of H2
                		p0=P(s)-pH2     									!Computing partial pressure of pristine
                		Gn(0) = ( E(0)*hartreetoJmol + R*T(k)*log(p0/dbp0) - R*T(k)*log(Z0))*JmoltoeV           !Computing G of pristine at given conditions
				CALL hydrogenfree(T(k), pH2, lambdaH2, ZvH2, ZrH2, ZiH2, Ghyd(2))
                        	!Cycle over molecules 
                        	do label=1,n
                                	F=3*(Ni+label)-6        						!Computing vib degrees of the 'label'th specie
                                	CALL partition(T(k), nu, n, label, F, Fmax, ele, symm, Trot, Zvib, Zrot, Zint)    !Computing Z of the 'label'th specie
                                	CALL dbpressure(T(k),  molar, n, label, dbp)
                                	Gn(label) = (E(label)*hartreetoJmol + R*T(k)*log(P(s)/dbp) - R*T(k)*log(Zint))*JmoltoeV   !Computing G
                                	DGf(label) = Gn(label)-Gn(0)-(label/2.0d0)*Ghyd(2) !Computing DeltaG
                        	end do
                        	minlabel=MINLOC(DGf, DIM=1) !Looking for smallest values at given conditions
                        	minen=MINVAL(DGf)
				write(500,*) T(k), P(s), fr(i), "n = ", minlabel, minen
	               	end do
        	end do
	end do

	write(*,*) 'Phase diagram done!!'

	deallocate(DGf)
	deallocate(Gn)
	deallocate(T,P)
end if

end program
