!***********************************************************************
! Interface to HELMHOLTZ EoS for XNet 6.0 7/6/2010
! This file contains routines which calculate EoS quantites needed to
! calculate screening corrections for reaction rates.
!***********************************************************************

Subroutine eos_initialize
!-----------------------------------------------------------------------
! This routine initializes the Helmholtz EoS
!-----------------------------------------------------------------------
call read_helm_table
End Subroutine eos_initialize


Subroutine eos_interface(t9,rho,y,ztilde,zinter,lambda0,gammae)
!-----------------------------------------------------------------------
! This routine calls the Helmholtz EOS with the input temperature, 
! density and composition.  It returns the factors needed for screening.
!-----------------------------------------------------------------------
  Use constants
  Use controls
  include 'vector_eos.dek'
  Real(8), Intent(in) :: t9, rho, y(*)
  Real(8), Intent(out) :: ztilde, zinter, lambda0, gammae
  Real(8) :: ye,ytot,bkt,abar,zbar,z2bar,zibar
  Real(8) :: etae,sratio,efermkt,rel_ef,emass,efc,ae
  Real(8) :: onethird=1./3.,twothird=2./3.

! Calculate Ye and other needed moments of the abundance distribution
  call y_moment(y,ye,ytot,abar,zbar,z2bar,zibar)  

! Load input variables for the eos 
  jlo_eos=1
  jhi_eos=1
  den_row(1)=rho
  temp_row(1)=t9*1e9
  abar_row(1)=abar
  zbar_row(1)=ye*abar

! Call the eos
  call helmeos
  etae=etaele_row(1) 
  
! Calculate electon distribution
  bkt=bok*T9 
  emass=ele_en/clt**2
  rel_ef=hbar*(3.0*pi**2*rho*avn*ye)**onethird/(emass*clt)
  efermkt=ele_en*(sqrt(1+rel_ef**2)-1)/bkt
  efc=.5*hbar**2*(3.0*pi**2*rho*avn*ye)**twothird/(emass*bkt)
 !Write(lun_diag,'(a4,6es12.5)') 'MUh',bkt,etae,efermkt,efc,rel_ef

! Calculate ratio f'/f for electrons (Salpeter, Eq. 24)
  call salpeter_ratio(etae,sratio)
  ztilde=sqrt(z2bar+zbar*sratio)
  
! Calculate plasma quantities
  lambda0=sqrt(4*pi*rho*avn*ytot)*(e2/bkt)**1.5 ! DGC, Eq. 3
  ae=(3./(4.*pi*avn*rho*ye))**onethird ! electron-sphere radius
  gammae=e2/(ae*bkt) ! electron Coulomb coupling parameter 
  zinter=zibar/(ztilde**.58*zbar**.28)
  If(idiag>=0) Write(lun_diag,'(a14,9es12.5)') 'Helmholtz EOS', t9, rho, ye,z2bar,zbar,sratio, ztilde, ztilde*lambda0, gammae
  
  Return
End Subroutine eos_interface

Subroutine salpeter_ratio(eta,ratio)
!-----------------------------------------------------------------------------
! This routine calculates the salpeter (1954) ratio f'/f(eta) needed for 
! electron screening.  eta is the ratio of electron chemical potential 
! to kT.  
!
! Calculation Uses Fermi function relation d/dx f_(k+1) = (k+1) f_k and 
! the rational function expansions of Antia (1993; ApJS 84 101) for the
! Fermi-Dirac integrals of order 1/2 and -1/2. Thanks to Frank Timmes.
!-----------------------------------------------------------------------------
  
  Integer :: i
  Real(8) :: eta,num,den,x,xn(12),fermip,fermim,ratio
! Coefficients for F_+1/2 
  Real(8) :: ap(8) = (/ 5.75834152995465d+06, 1.30964880355883d+07, &
 &                      1.07608632249013d+07, 3.93536421893014d+06, &
 &                      6.42493233715640d+05, 4.16031909245777d+04, &
 &                      7.77238678539648d+02, 1.0d+00/)
  Real(8) :: bp(8) = (/ 6.49759261942269d+06, 1.70750501625775d+07, &
 &                      1.69288134856160d+07, 7.95192647756086d+06, &
 &                      1.83167424554505d+06, 1.95155948326832d+05, &
 &                      8.17922106644547d+03, 9.02129136642157d+01/)
  Real(8) :: cp(11)= (/ 4.85378381173415d-14, 1.64429113030738d-11, &
 &                      3.76794942277806d-09, 4.69233883900644d-07, &
 &                      3.40679845803144d-05, 1.32212995937796d-03, &
 &                      2.60768398973913d-02, 2.48653216266227d-01, &
 &                      1.08037861921488d+00, 1.91247528779676d+00, &
 &                      1.0d+00/)
  Real(8) :: dp(12)= (/ 7.28067571760518d-14, 2.45745452167585d-11, &
 &                      5.62152894375277d-09, 6.96888634549649d-07, &
 &                      5.02360015186394d-05, 1.92040136756592d-03, &
 &                      3.66887808002874d-02, 3.24095226486468d-01, &
 &                      1.16434871200131d+00, 1.34981244060549d+00, &
 &                      2.01311836975930d-01,-2.14562434782759d-02/)
! Coefficients for F_-1/2 
  Real(8) :: am(8) = (/ 1.71446374704454d+07, 3.88148302324068d+07, &
 &                      3.16743385304962d+07, 1.14587609192151d+07, &
 &                      1.83696370756153d+06, 1.14980998186874d+05, &
 &                      1.98276889924768d+03, 1.0d+00/)
  Real(8) :: bm(8) = (/ 9.67282587452899d+06, 2.87386436731785d+07, &
 &                      3.26070130734158d+07, 1.77657027846367d+07, &
 &                      4.81648022267831d+06, 6.13709569333207d+05, &
 &                      3.13595854332114d+04, 4.35061725080755d+02/)
  Real(8) :: cm(12)= (/-4.46620341924942d-15,-1.58654991146236d-12, &
 &                     -4.44467627042232d-10,-6.84738791621745d-08, &
 &                     -6.64932238528105d-06,-3.69976170193942d-04, &
 &                     -1.12295393687006d-02,-1.60926102124442d-01, &
 &                     -8.52408612877447d-01,-7.45519953763928d-01, &
 &                      2.98435207466372d+00, 1.0d+00/)
  Real(8) :: dm(12)= (/-2.23310170962369d-15,-7.94193282071464d-13, &
 &                     -2.22564376956228d-10,-3.43299431079845d-08, &
 &                     -3.33919612678907d-06,-1.86432212187088d-04, &
 &                     -5.69764436880529d-03,-8.34904593067194d-02, &
 &                     -4.78770844009440d-01,-4.99759250374148d-01, &
 &                      1.86795964993052d+00, 4.16485970495288d-01/)

! Rational functions have two domains of applicability.
! Lower Branch of rational functions
  If (eta < 2.0) Then
    x = exp(eta)
    xn(1) = 1.0
    Do i=2,8
      xn(i) = x*xn(i-1)
    Enddo

! Calculate the +1/2 fermin integral
    num = sum(ap*xn(1:8))
    den = sum(bp*xn(1:8))
    fermip = x*num/den

! Calculate the -1/2 fermin integral
    num = sum(am*xn(1:8))
    den = sum(bm*xn(1:8))
    fermim = x*num/den

! Upper Branch of rational functions
  Else
    x = 1.0/eta**2
    xn(1) = 1.0
    Do i=2,12
      xn(i) = x*xn(i-1)
    Enddo

! Calculate the +1/2 fermin integral       
    num = sum(cp*xn(1:11))
    den = sum(dp*xn(1:12))
    fermip = eta*sqrt(eta)*num/den

! Calculate the -1/2 fermin integral
    num = sum(cm*xn(1:12))
    den = sum(dm*xn(1:12))
    fermim = sqrt(eta)*num/den
  Endif

! Evalutate the salpeter ratio
  ratio = 0.5 * fermim/fermip
! write(lun_diag,"(1x,4es12.4)") eta,ratio,fermim,fermip
  Return
End Subroutine salpeter_ratio







                                                                                                                                                                    
