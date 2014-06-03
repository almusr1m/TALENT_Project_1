!***********************************************************************
! Screening for XNet 6.0 7/6/2010
! This file contains the routines needed to calculate screening 
! corrections for reaction rates.
!***********************************************************************

Subroutine screening()
!-----------------------------------------------------------------------------
! This routine calculates the screening factors necessary for Xnet.
! The HELMHOLTZ Equation of State (Timmes & Swesty 1999) is Used to 
! determine the electron distribution and chemical potential.
!
! References:
! Weak Screening:         
!    Salpeter (1954) Aust J Phys 7 353.
! Intermediate Screening: 
!    DeWitt, Graboske & Cooper (1973) ApJ 181 439.
!    Graboske, DeWitt, Grossman & Cooper (1973) ApJ 181 457.
! Strong Screening:       
!    DeWitt & Slattery (1999) Contrib Plasma Phys 39 97.
!-----------------------------------------------------------------------------
  Use controls
  Use constants
  Use abundances
  Use conditions
  Use nuclear_data              
  Use cross_sect_data
  
  Integer :: iz1,iz2,iz3,j,mu
  Real(8) :: z1,z2,z3,z12,z123,a1,a2,a3,a12,a123
  Real(8) :: ztilde,zinter,lambda0,gammae,gamma12,gamma123
  Real(8) :: h0,hw,hi,hs,lambda12,gamma,lambda123,z,hi2
  Real(8) :: onethird=1./3.,twothird=2./3.,fivethird=5./3.,bkt
  Real(8) :: fhs(0:izmax+2),fhi(0:izmax+2)
  Real(8) :: cds(5)=(/-.899172,0.602249,-.274823,-1.401915,0.3230064/)

! Call EOS to get plasma quantities
  call eos_interface(t9t,rhot,yt,ztilde,zinter,lambda0,gammae)

!-----------------------------------------------------------------------------
! Calculate screening energies as a function of Z, for prescriptions that 
! follow this approach 
!-----------------------------------------------------------------------------
  fhi(0)=0.0
  fhs(0)=0.0
  Do j=1,izmax+2
    z=Real(j)
    fhi(j)=0.38*zinter*lambda0**.86*(z)**1.86
    gamma=(z)**fivethird*gammae                            
    fhs(j)=cds(1)*gamma+cds(2)*gamma**cds(5)/cds(5)+cds(3)*log(gamma)+cds(4)
  Enddo

! Loop over 1 reactanct reactions to build screening factors.      
  h1=0.0

! Loop over 2 reactanct reactions to build screening factors.      
  Do mu=1,nreac(2)                                                      
    z1=zz(n2i(1,mu)) ; z2=zz(n2i(2,mu))
    iz1=int(z1)      ; iz2=int(z2)                                                  
    a1=aa(n2i(1,mu)) ; a2=aa(n2i(2,mu))                                                  
    z12=z1*z2
    If(z12==0.0.or.iscrn==0) Then
      h2(mu)=0.0                                                        
    Else

! Weak and intermediate screening factors, Table 4 of Graboske et al.         
      lambda12=z1*z2*ztilde*lambda0    
      hw=lambda12
!     hw=lambda12*(1.+lambda12*(log(lambda12)+0.8364))
      hi=0.38*zinter*lambda0**.86*((z1+z2)**1.86-z1**1.86-z2**1.86)
      hi2=fhi(iz1+iz2)-fhi(iz1)-fhi(iz2)                                                     

! Strong screening from Dewitt & Slattery using linear mixing.
      hs=fhs(iz1)+fhs(iz2)-fhs(iz1+iz2)

! Strong Screening from Itoh, non-radial component
      gamma12=2.0*z1*z2*gammae/(z1**onethird+z2**onethird)
      hs=1.25*gamma12

! Select Screening factor
!     h2(mu)=min(hw,hi,hs)                                                         
      If(lambda12<0.1) Then
        h2(mu)=hw
      Elseif(lambda12<2.0) Then
        h2(mu)=hi
      Elseif(lambda12>5.0) Then
        h2(mu)=hs
      Else
        h2(mu)=min(hi,hs)
      Endif
      gamma=(z12)**fivethird*gammae
      If(idiag>3) Write(lun_diag,'(3a5,i6,8es12.5)') 'H2', &
     &  nname(n2i(1:2,mu)),mu, lambda12,gamma12,h2(mu),hw,hi,hi2,hs
    EndIf
  Enddo

! Loop over 3 reactant reactions to build screening factors.
  Do mu=1,nreac(3)                                                      
    z1=zz(n3i(1,mu)) ; z2=zz(n3i(2,mu)) ; z3=zz(n3i(3,mu))
    iz1=int(z1)      ; iz2=int(z2)      ; iz3=int(z3)                                                  
    a1=aa(n3i(1,mu)) ; a2=aa(n3i(2,mu)) ; a3=aa(n3i(3,mu))
    z12 =z1*z2
    z123=z1*z2*z3
    If(z123==0.0.or.iscrn==0) Then
      h3(mu)=0.0
    Else

! Weak and intermediate screening factors, Table 4 of Graboske et al.         
      lambda123=(z1*z2+z1*z3+z2*z3)*ztilde*lambda0    
      hw=lambda123
      hi=fhi(iz1+iz2+iz3)-fhi(iz1)-fhi(iz2)-fhi(iz3)

! Strong screening from Dewitt & Slattery using linear mixing.
      hs=fhs(iz1)+fhs(iz2)+fhs(iz3)-fhs(iz1+iz2+iz3)

! Strong Screening from Itoh, non-radial component
      gamma123=2.0*gammae* &
    & ((z12/(z1**onethird+z2**onethird))+(z123/(z12**onethird+z3**onethird)))
!     hs=1.25*gamma123

! Select Screening factor
!     h3(mu)=min(hw,hi,hs)
      If(lambda12<0.1) Then
        h3(mu)=hw
      Elseif(lambda12<2.0) Then
        h3(mu)=hi
      Elseif(lambda12>5.0) Then
        h3(mu)=hs
      Else
        h3(mu)=min(hi,hs)
      Endif
      If(idiag>3) Write(lun_diag,'(4a5,i6,8es12.5)') 'H3',nname(n3i(1:3,mu)),mu, lambda123,gamma123,h3(mu),hw,hi,hs
    EndIf
  Enddo

  Return                                                            
End Subroutine screening

Subroutine y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)  
!------------------------------------------------------------------------------  
! This routine calculates moments of the abundance distribution for the EOS.
!------------------------------------------------------------------------------  
  Use controls
  Use nuclear_data
  Real(8), Intent(in)  :: y(ny)
  Real(8), Intent(out) :: ye,abar,zbar,z2bar,zibar
  Real(8)              :: ytot,atot,ztot

! Calculate abundance moments
  ytot =sum(y(:))
  atot =sum(aa*y)
  ztot =sum(zz*y)
  abar =atot/ytot
  zbar =ztot/ytot
  z2bar=sum(zz*zz*y)/ytot
  zibar=sum(zz**1.58*y)/ytot
  ye=ztot
  Write(lun_diag,'(a4,6es12.5)') 'YMom',ytot,abar,zbar,z2bar,zibar,ye
  
  Return                                                                    
End Subroutine y_moment                                                                       






                                                                                                                                                                    
