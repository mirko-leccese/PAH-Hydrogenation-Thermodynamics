SUBROUTINE RandomWalk(W1, fix, n, Nctot0, Nhtot0, Ncinp, Nhinp, Nctrial, Nhtrial, xctrial, xhtrial)
  
  REAL(8), DIMENSION(0:n) :: Ncinp
  REAL(8), DIMENSION(1:2) :: Nhinp
  REAL(8), DIMENSION(0:n) :: Nctrial
  REAL(8), DIMENSION(1:2) :: Nhtrial   !1 is for H and 2 for H2 
  REAL(8), DIMENSION(0:n) :: dNc
  REAL(8), DIMENSION(1:2) :: dNh
  REAL(8), DIMENSION(0:n) :: xctrial
  REAL(8), DIMENSION(1:2) :: xhtrial
  real(8) :: Nctot0, Nctot, Nhtot0, Nhtot
  real(8) :: checkC, checkH
  real(8) :: highv, lowv, mean, Th, checkiter
  INTEGER(4) :: k, m, q, fix 
  LOGICAL :: seq
    REAL(8) :: W1  !In the original version, C=1.d-2
  INTEGER(4) :: l 
  
  REAL(8) :: ran, i0, i1, diff
  REAL(8) :: dsum0, dsum1
  real(8) :: dN1
  !REAL(8), PARAMETER :: Nav=6.02214085774d+23
  
  dsum0=0.0d0
  dsum1=0.0d0
  
  highv=maxval(Ncinp)
  lowvv=minval(Ncinp)
  

!open(Unit=2000,File='Variations.dat', Position='Append')

DO l=0,fix
	
	CALL random_number(ran)
	!dNc(l)=diff*ran+i0
	dNc(l)=(2.0d0*(ran-0.50d0)*Ncinp(l))/W1
	!dsum0=dsum0+dNc(l)
	!dsum1=dsum1+l*dNc(l)

	Nctrial(l)=Ncinp(l)+dNc(l)
	!write(2000,*) l, Ncinp(l), dNc(l), Nctrial(l)
	Nctot=Nctot+Nctrial(l)
END DO

Nctrial=Nctrial*(Nctot0/Nctot)

DO l=0,fix
	checkC=checkC+Nctrial(l)
END DO

!write(2000,*) Nctrial(0), Nctrial(1), Nctrial(2), checksum 

!dNc(fix)=-dsum0
!dsum1=dsum1+fix*dNc(fix)
!Nctrial(fix)=Ncinp(fix)+dNc(fix)
!write(2000,*) fix, Ncinp(fix), dNc(fix), Nctrial(fix)
!Ntot=Ntot+Nctrial(fix)

!write(2000,*) Nctrial(1)
!Other gradient components

!IF ( fix < n ) THEN 
!	DO l=fix+1, n
!		dNc(l) = 0.0d0
!		Nctrial(l)=Ncinp(l)+dNc(l)
!		Ntot=Ntot+Nctrial(l)
!	END DO
!END IF 
       
!Let's compute the true variation of H due to PAhs
do l=0,fix
	dN1=Nctrial(l)-Ncinp(l)
	dsum1=dsum1+l*dN1
end do
!Old version where only H was selected randmoly

CALL random_number(ran)
dNh(1)=(2.0*(ran-0.5d0)*Nhinp(1))/W1
!write(*,*) diff, dNh(1)
dNh(2)=(1.0d0/2.0d0)*(-dNh(1)-dsum1) !Old version
Nhtrial(1)=Nhinp(1)+dNh(1)
Nhtrial(2)=Nhinp(2)+dNh(2)
!write(2000,*) 'H', Nhinp(1), dNh(1), Nhtrial(1)
!write(2000,*) 'H2', Nhinp(2), dNh(2), Nhtrial(2)

!New version where both H and H2 are randomly chosed and then the hydrogen population is renormalized
!CALL random_number(ran)
!dNh(1)=(2.0*(ran-0.5d0)*Nhinp(1))/W1
!CALL random_number(ran)
!dNh(2)=(2.0*(ran-0.5d0)*Nhinp(2))/W1
!
!Nhtrial(1)=Nhinp(1)+dNh(1)
!Nhtrial(2)=Nhinp(2)+dNh(2)

!Unrenormalized H pop is
!Nhtot=Nhtrial(1)+2.0d0*Nhtrial(2)+dsum1
!Renormalizing H pop
!hratio=(Nhtot0-dsum1)/(Nhtot-dsum1)
!Nhtrial(1)=Nhtrial(1)*hratio
!Nhtrial(2)=Nhtrial(2)*hratio

!Checking conservation of H
checkH = Nhtrial(1)+2.0d0*Nhtrial(2)+dsum1 

!write(2000,*) Nctrial(1), Nhtrial(1), Nhtrial(2), Nctot0, Nhtot0, checkC, checkH

Ntot=Nctot0+Nhtrial(1)+Nhtrial(2)
!close(2000)
        
DO l=0,n
	xctrial(l)=Nctrial(l)/Ntot
END DO

xhtrial(1)=Nhtrial(1)/Ntot
xhtrial(2)=Nhtrial(2)/Ntot

END SUBROUTINE RandomWalk
