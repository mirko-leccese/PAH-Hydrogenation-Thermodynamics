SUBROUTINE RandomWalk(W1, fix, n, Nctot0, Nhtot0, Ncinp, Nhinp, Nctrial, Nhtrial, xctrial, xhtrial)

! This subroutine generates a random variations on each mixture components:
!
! Args:
!   W1 (real)     : parameter that allows to tune the variation interval
!   fix (integer) : the last specie whose abundance needs to be modified at the given step
!	  n (integer)   : total number of hydrogenations 
!	  Nctot0 (real) : initial total number of PAH molecules 
!	  Nhtot0 (real) : initial total number of H/H2 molecules
!   Ncinp (real, array)  : input total number of PAH molecules at the given optimization step
!   Nhinp (real, array)  : input total numbers of H/H2 molecules at the given optimization step
!
! Out:
!	  Nctrial (real, array): trial number of PAH molecules 
!   Nhtrial (real, array): trial number of H/H2 molecules 
!   xctrial (real, array): trial molar fractions of PAH molecules 
!   xhtrial (real, array): trial molar fractions of H/H2 molecules

REAL(8), DIMENSION(0:n) :: Ncinp
REAL(8), DIMENSION(1:2) :: Nhinp
REAL(8), DIMENSION(0:n) :: Nctrial
REAL(8), DIMENSION(1:2) :: Nhtrial   ! 1 is for H and 2 for H2 
REAL(8), DIMENSION(0:n) :: dNc
REAL(8), DIMENSION(1:2) :: dNh
REAL(8), DIMENSION(0:n) :: xctrial
REAL(8), DIMENSION(1:2) :: xhtrial
REAL(8) :: Nctot0, Nctot, Nhtot0, Nhtot
REAL(8) :: highv, lowv, mean, Th, checkiter
REAL(8) :: W1  
REAL(8) :: ran, i0, i1, diff
REAL(8) :: dsum0, dsum1
REAL(8) :: dN1
INTEGER(4) :: l 
INTEGER(4) :: k, m, q, fix 
  
dsum0 = 0.0d0
dsum1 = 0.0d0

! Finding most and least abundant specie in the current mixture   
highv = maxval(Ncinp)
lowvv = minval(Ncinp)

DO l=0,fix
  CALL random_number(ran)
	dNc(l) = (2.0d0*(ran-0.50d0)*Ncinp(l))/W1
	Nctrial(l) = Ncinp(l)+dNc(l)
	Nctot = Nctot+Nctrial(l)
END DO

Nctrial = Nctrial*(Nctot0/Nctot)

! Other gradient components
! IF ( fix < n ) THEN 
!	  DO l=fix+1, n
!		  dNc(l) = 0.0d0
!		  Nctrial(l)=Ncinp(l)+dNc(l)
!		  Ntot=Ntot+Nctrial(l)
!	  END DO
! END IF 
       
! Let's compute the true variation of H due to PAhs
DO l = 0,fix
	dN1 = Nctrial(l)-Ncinp(l)
	dsum1 = dsum1+l*dN1
END DO

CALL random_number(ran)
dNh(1) = (2.0*(ran-0.5d0)*Nhinp(1))/W1
dNh(2) = (1.0d0/2.0d0)*(-dNh(1)-dsum1) 
Nhtrial(1) = Nhinp(1)+dNh(1)
Nhtrial(2) = Nhinp(2)+dNh(2)

! New version where both H and H2 are randomly chosed and then the hydrogen population is renormalized
! CALL random_number(ran)
! dNh(1)=(2.0*(ran-0.5d0)*Nhinp(1))/W1
! CALL random_number(ran)
! dNh(2)=(2.0*(ran-0.5d0)*Nhinp(2))/W1
!
! Nhtrial(1)=Nhinp(1)+dNh(1)
! Nhtrial(2)=Nhinp(2)+dNh(2)

! Unrenormalized H pop is
! Nhtot=Nhtrial(1)+2.0d0*Nhtrial(2)+dsum1
! Renormalizing H pop
! hratio=(Nhtot0-dsum1)/(Nhtot-dsum1)
! Nhtrial(1)=Nhtrial(1)*hratio
! Nhtrial(2)=Nhtrial(2)*hratio

DO l = 0,n
	xctrial(l) = Nctrial(l)/Ntot
END DO

xhtrial(1) = Nhtrial(1)/Ntot
xhtrial(2) = Nhtrial(2)/Ntot

END SUBROUTINE RandomWalk
