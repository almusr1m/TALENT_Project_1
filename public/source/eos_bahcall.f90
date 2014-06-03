!***********************************************************************
! EoS Replacement based on Bahcall () for XNet 6.0 7/6/2010
! This file contains routines which calculate EoS quantites needed to
! calculate screening corrections for reaction rates.
!***********************************************************************

Subroutine eos_initialize
!-----------------------------------------------------------------------
! This routine initializes the EoS
!-----------------------------------------------------------------------
End Subroutine eos_initialize


Subroutine eos_interface(t9,rho,y,ztilde,zinter,lambda0,gammae)
!-----------------------------------------------------------------------
! This routine Bahcall's approach to calculate the factors needed for 
! screening from the input temperature, density and composition.
!-----------------------------------------------------------------------
  Use constants
  Use controls
  Real(8), Intent(in) :: t9, rho, y(*)
  Real(8), Intent(out) :: ztilde, zinter, lambda0, gammae
  Real(8) :: ye,ytot,bkt,abar,zbar,z2bar,zibar
  Real(8) :: etae_mb,sratio,efermkt,rel_ef,emass,efc,ae
  Real(8) :: onethird=1./3.,twothird=2./3.

! Calculate Ye and other needed moments of the abundance distribution
  call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)  

! Calculate electon distribution
  bkt=bok*T9
  emass=ele_en/clt**2
  etae_mb=log(rho*avn*ye*(2.0*pi*hbar**2/(emass*bkt))**1.5)
  rel_ef=hbar*(3.0*pi**2*rho*avn*ye)**onethird/(emass*clt)
  efermkt=ele_en*(sqrt(1+rel_ef**2)-1)/bkt
  efc=.5*hbar**2*(3.0*pi**2*rho*avn*ye)**twothird/(emass*bkt)
 !Write(lun_diag,'(a4,6es12.5)') 'MUb',bkt,etae_mb,efermkt,efc,rel_ef

! Calculate ratio f'/f for electrons (Salpeter, Eq. 24)
  call salpeter_ratio(efermkt,sratio)
  ztilde=sqrt(z2bar+zbar*sratio)
  
! Calculate plasma quantities
  lambda0=sqrt(4*pi*rho*avn*ytot)*(e2/bkt)**1.5 ! DGC, Eq. 3
  ae=(3./(4.*pi*avn*rho*ye))**onethird ! electron-sphere radius
  gammae=e2/(ae*bkt) ! electron Coulomb coupling parameter 
  zinter=zibar/(ztilde**.58*zbar**.28)
  If(idiag>=0) Write(lun_diag,'(a14,9es12.5)') 'Bahcall SCRN', t9, rho, ye, z2bar, zbar, sratio, ztilde,ztilde*lambda0, gammae
  
  Return
End Subroutine eos_interface

Subroutine salpeter_ratio(efmkt,ratio)
!-----------------------------------------------------------------------------
! This routine calculates the Salpeter (1954) ratio f'/f(eta) needed for
! electron screening, using a fit to Figure 24 in that paper.  
! efmkt is the ratio of electron chemical potential to kT.
!-----------------------------------------------------------------------------
  Real(8) :: efmkt,lefmkt,ratio
  lefmkt=log10(efmkt)
  If(lefmkt<=-1.578347) Then ! Bahcall uses limit of -2, but this gives ratio slightly above 1.
    ratio=1.0
  Elseif(lefmkt>=1.5) Then
    ratio=0.0
  Else
    ratio=0.75793-(0.54621*lefmkt)-(0.30964*lefmkt**2)+(0.12535*lefmkt**3) &
&    +(0.1203*lefmkt**4)-(0.012857*lefmkt**5)-(0.014768*lefmkt**6)
Endif

End Subroutine salpeter_ratio









