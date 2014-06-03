!*******************************************************************************
! Jacobian Sparce for MA48, part of XNet 7, 6/2/10
!
! The routines in this file are used to replace the standard dense Jacobian 
! and associated solver with the Harwell MA48 sparse solver package.  
!
! The bulk of the computational cost of the network (60-95%) is the solving 
! of the matrix equation.  Careful selection of the matrix solver is therefore 
! very important to fast computation.  For networks from a few dozen up to a 
! couple hundred species, hand tuned dense solvers such as those supplied by the 
! hardware manufacturer (often LAPACK) or third-parties like NAG, IMSL, etc. are 
! the fastest. However for larger matrices, sparse solvers are faster.  MA48 is 
! an obsolescent solver from the Harwell Subroutine Library, which is therefore 
! publically available.  We find it to be less bulletproof and slower at large 
! network sizes than PARDISO, but the availablilty of the MA48 source makes it 
! valuable for some applications. 
!*******************************************************************************

Module jacobian_data
!===============================================================================
! Contains data for use in the sparse solver.
!===============================================================================
  Use nuclear_data
  Real(8), Dimension(:), Allocatable :: vals,tvals,sident,wB,wC
  Integer, Dimension(:), Allocatable :: ridx,cidx
  Integer, Dimension(:), Allocatable :: irn,jcn
  Integer, Dimension(:), Allocatable :: ns11,ns21,ns22
  Integer, Dimension(:), Allocatable :: ns31,ns32,ns33
  Integer, Dimension(:), Allocatable :: keep,iwA,iwB,iwC
  Integer :: lval,lia,lia_min,jobA,jobB,jobC,icntl(20),info(20)
  Real(8) :: cntl(10),rinfo(10),maxerr
  Logical :: trans = .false.
!$OMP THREADPRIVATE(tvals,vals,lia,lia_min,jcn,irn,wB,wC,iwA,iwB,iwC, &
!$OMP   keep,info,rinfo,cntl,icntl,jobA,jobB,jobC)

  Real(8), Dimension(:), Allocatable :: mc29r,mc29c,mc29w
  Integer :: ifail
  Logical :: scale_jac
!$OMP THREADPRIVATE(mc29r,mc29c,mc29w,ifail,scale_jac)
End Module jacobian_data
  
Subroutine read_jacobian_data(data_dir)
!===============================================================================
! Reads in data necessary to use sparse solver and initializes the Jacobian data.
!===============================================================================
  Use controls
  Use reac_rate_data
  Use jacobian_data
       
  Character (LEN=*),  Intent(in)  :: data_dir
  Integer :: i,pb(ny+1)
  
!$OMP PARALLEL DEFAULT(SHARED)
! Set default values for MA48 control parameters
  Call MA48ID(cntl,icntl)

! Set error/diagnostic output to XNet diagnostic file (default=6)
  icntl(1) = lun_diag
  icntl(2) = lun_diag

! Set level of MA48 error/diagnostic output (default=2)
  If(idiag>=4) Then
    icntl(3) = 3
  EndIf

! Limit pivot search to maximum icntl(4) columns or use special search technique for 0 (default=3)
! icntl(4) = 0

! Set block size (default=32)
  If (ny<32) Then
    icntl(5) = 16
  ElseIf (ny<128) Then
    icntl(5) = 32
  ElseIf (ny<512) Then
    icntl(5) = 64
  ElseIf (ny<2048) Then
    icntl(5) = 128
  Else
    icntl(5) = 256  
  EndIf

! Set the minimum size of a block of the block triangular form other than the final block (defualt=1)
! icntl(6) = ny

! Set option to move each column for which iwA(j)=0 to the end of the pivot sequence within its block (default=0)
! icntl(8) = 1

! Set option to automatically reset to job=1 if MA48BD fails for job=2
! icntl(11) = 1

! Set the pivoting threshold (near 0.0 emphasizes sparsity; near 1.0 emphasizes stability) (default=0.1)
! cntl(2) = 0.5
!$OMP END PARALLEL
  
  Open(600,file=trim(data_dir)//"/sparse_ind",status='old',form='unformatted')
  
  Read(600) lval

  Allocate(ridx(lval),cidx(lval),sident(lval))
!$OMP PARALLEL DEFAULT(SHARED)
  Allocate(tvals(lval))
!$OMP END PARALLEL
  
!$OMP PARALLEL DEFAULT(SHARED)
  lia = 8*lval
! These are variables to be used by the MA48 solver
  Allocate(vals(lia),jcn(lia),irn(lia))
  Allocate(wB(ny),wC(4*ny))
  Allocate(iwA(9*ny),iwB(4*ny),iwC(ny))
  If (icntl(8) == 0) Then
    Allocate(keep(max(ny/icntl(6),1)))
  Else
    Allocate(keep(6*ny+4*ny/icntl(6)+7))
  EndIf

! These are variables to be used by MC29 if scaling is deemed necessary
  Allocate(mc29r(ny),mc29c(ny),mc29w(5*ny))
  mc29r = 0.0
  mc29c = 0.0
  mc29w = 0.0
  scale_jac = .false.
!$OMP END PARALLEL

! Set the value for the maximum allowed error in the call to MA48CD
  maxerr = 1.0d-04
  
  Read(600) ridx,cidx,pb
  
! Build a compressed row format version of the identity matrix
  Do i=1,lval
    If (ridx(i)==cidx(i)) sident(i)=1.0
  EndDo
  
  Read(600) l1s,l2s,l3s
  
! Build  arrays for direct sparse representation Jacobian build
  Allocate(ns11(l1s))
  Allocate(ns21(l2s),ns22(l2s))
  Allocate(ns31(l3s),ns32(l3s),ns33(l3s))
  
  ns11 = 0
  ns21 = 0
  ns22 = 0
  ns31 = 0
  ns32 = 0
  ns33 = 0
  
  Read(600) ns11,ns21,ns22
  Read(600) ns31
  Read(600) ns32
  Read(600) ns33
  Close(600)  
  
End Subroutine read_jacobian_data

Subroutine jacobian_build(diag,mult)
!===============================================================================
! This routine calculates the Jacobian matrix dYdot/dY, and and augments 
! by multiplying all elements by mult and adding diag to the diagonal elements.
!===============================================================================
  Use controls
  Use conditions
  Use abundances
  Use reac_rate_data
  Use jacobian_data
  Use timers
  
  Real(8), Intent(in) :: diag, mult
  Integer :: i,j,kstep,i0,la1,le1,la2,le2,la3,le3,j1,l1,l2,l3
  Integer :: ls1,ls2,ls3
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_jacob = timer_jacob - start_timer  

! The quick solution to taking advantage of sparseness is to create a values array that has the maximum
! number of non-zero elements as well as a 2-by-#non-zero elements matrix and the other vectors
! required by the sparse solver.  The second matrix will contain the ordered pairs that map to a 
! particular place in the Jacobian.  Build the Jacobian as it is built now, Then pull the values 
! from it using the ordered pairs to fill the values array. 

! Build the Jacobian, species by species
  tvals = 0.0

  Do i0=1,ny
    la1=la(1,i0)
    le1=le(1,i0)
    Do j1=la1,le1
      ls1 = ns11(j1) ! ns11(j1) gives the index effected reaction j1 by in the compressed row storage scheme
      tvals(ls1)=tvals(ls1)+b1(j1)
    EndDo
    la2=la(2,i0) 
    le2=le(2,i0)  
    Do j1=la2,le2
      ls1=ns21(j1) ! ns21(j1) gives the first index effected reaction j1 by in the compressed row storage scheme
      ls2=ns22(j1) ! ns22(j1) gives the second index effected reaction j1 by in the compressed row storage scheme
      l1=n21(j1)   ! n21(k) gives the index of first reactant in reaction mu2(k)
      l2=n22(j1)   ! n22(k) gives the index of second reactant in reaction mu2(k)
      tvals(ls1)=tvals(ls1)+b2(j1)*yt(l2)
      tvals(ls2)=tvals(ls2)+b2(j1)*yt(l1)
    EndDo
    la3=la(3,i0)
    le3=le(3,i0)
    Do j1=la3,le3
      ls1=ns31(j1) ! ns31(j1) gives the first index effected reaction j1 by in the compressed row storage scheme
      ls2=ns32(j1) ! ns32(j1) gives the second index effected reaction j1 by in the compressed row storage scheme
      ls3=ns33(j1) ! ns33(j1) gives the third index effected reaction j1 by in the compressed row storage scheme
      l1=n31(j1)   ! n21(k) gives the index of first reactant in reaction mu2(k)
      l2=n32(j1)   ! n22(k) gives the index of second reactant in reaction mu2(k)
      l3=n33(j1)   ! n22(k) gives the index of third reactant in reaction mu2(k)
      tvals(ls1)=tvals(ls1)+b3(j1)*yt(l2)*yt(l3)
      tvals(ls2)=tvals(ls2)+b3(j1)*yt(l1)*yt(l3)
      tvals(ls3)=tvals(ls3)+b3(j1)*yt(l1)*yt(l2)
    EndDo
  EndDo

! Augment matrix with externally provided factors  
  tvals = mult * tvals
  tvals = tvals + sident * diag 
  
  If(idiag>=5) Then
    Write(lun_diag,"(a9,2es14.7)") 'JAC_build',diag,mult
    Write(lun_diag,"(14es9.1)") tvals
  EndIf

! Stop timer
  stop_timer = xnet_wtime()
  timer_jacob = timer_jacob + stop_timer

  Return   
End Subroutine jacobian_build
  
Subroutine jacobian_solve(kstep,rhs,dy) 
!===============================================================================
! This routine solves the system of abundance equations composed of the jacobian
! matrix and rhs vector.
!===============================================================================
  Use controls
  Use jacobian_data
  Use timers
  Integer, Intent(in)  :: kstep
  Real(8), Intent(IN)  :: rhs(ny)
  Real(8), Intent(out) ::  dy(ny)
  
  Call jacobian_decomp(kstep)
  Call jacobian_bksub(rhs,dy)
  
! Diagnostic output
  If(idiag>=4) Then
    Write(lun_diag,"(a)") 'JAC_SOLV'
    Write(lun_diag,"(14es10.3)") dy
  EndIf

  Return
End Subroutine jacobian_solve                                                                       
  
Subroutine jacobian_decomp(kstep) 
!===============================================================================
! This routine performs a matrix decomposition for the jacobian
!===============================================================================
  Use controls
  Use jacobian_data
  Use timers
  Integer, Intent(in)  :: kstep
  Integer :: i,j,kdecomp 
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_solve = timer_solve - start_timer  
  
! Use Harwell MA48 sparse solver (this one is portable)
  If (kstep == 1) Then
    jobA = 1
    jobB = 1
    vals = 0.0
    jcn = 0
    irn = 0
    keep = 0
    Do i=1,lval
      jcn(i) = cidx(i)
      irn(i) = ridx(i)
    EndDo
    vals(1:lval) = tvals

    If (scale_jac) Then
      ifail = 0
      Call MC29AD(ny,ny,lval,vals,irn,jcn,mc29r,mc29c,mc29w,lun_diag,ifail)
      If(ifail<0) Write(lun_diag,"(1x,a,i2)") 'Error during MC29AD.  Code: ',ifail
      mc29r(1:ny) = exp(mc29r(1:ny))
      mc29c(1:ny) = exp(mc29c(1:ny))
      vals(1:lval) = vals(1:lval)*mc29r(irn(1:lval))*mc29c(jcn(1:lval))
    EndIf

! Perform symbolic decomposition
    iwA = 0
    info = 0
    rinfo = 0.0
    Call MA48AD(ny,ny,lval,jobA,lia,vals,irn,jcn,keep,cntl,icntl,iwA,info,rinfo)
    If(info(1)/=0) Write(lun_diag,"(1x,a,i5,a,i2)") 'kstep=',kstep,', Warning during MA48AD.  Code: ',info(1)

    Do kdecomp = 1,10
      If (info(1) == 0) Exit

      If (info(3) > lia) Then
        Write(lun_diag,"(1x,3(a,i7))") 'info(3) = ',info(3),' > lia = ',lia,', info(4) = ',info(4)
        Write(lun_diag,*) 'Trying to reAllocate MA48 arrays'
        Deallocate(irn,jcn,vals)
        lia = info(4)
        Allocate(irn(lia),jcn(lia),vals(lia))
      ElseIf (info(2) > 10) Then
        Write(lun_diag,"(1x,a,i3,a,i7)") 'MA48 arrays too small, info(2) = ',info(2),' > 10, info(4) = ',info(4)
        Write(lun_diag,*) 'Trying to reAllocate MA48 arrays'
        Deallocate(irn,jcn,vals)
        lia = int(1.2*info(4))
        Allocate(irn(lia),jcn(lia),vals(lia))
      ElseIf (info(2) <= 1) Then
        Write(lun_diag,"(1x,a,i3,a,i7)") 'MA48 arrays too large, info(2) = ',info(2),' <= 1, info(4) = ',info(4)
        Write(lun_diag,*) 'Trying to reAllocate MA48 arrays'
        Deallocate(irn,jcn,vals)
        lia = max(int(0.5*lia),int(1.2*info(4)))
        Allocate(irn(lia),jcn(lia),vals(lia))
      EndIf  

! Rebuild vals array
      vals = 0.0
      irn = 0
      jcn = 0
      Do i=1,lval
        jcn(i) = cidx(i)
        irn(i) = ridx(i)
      EndDo
      vals(1:lval) = tvals

      If (scale_jac) Then
        vals(1:lval) = vals(1:lval)*mc29r(irn(1:lval))*mc29c(jcn(1:lval))
      ElseIf (kdecomp>2) Then
        scale_jac = .true.
        Write(lun_diag,*) 'Trying to scale input matrix'
        ifail = 0
        Call MC29AD(ny,ny,lval,vals,irn,jcn,mc29r,mc29c,mc29w,lun_diag,ifail)
        If(ifail<0) Write(lun_diag,"(1x,a,i2)") 'Error during MC29AD.  Code: ',ifail
        mc29r(1:ny) = exp(mc29r(1:ny))
        mc29c(1:ny) = exp(mc29c(1:ny))
        vals(1:lval) = vals(1:lval)*mc29r(irn(1:lval))*mc29c(jcn(1:lval))
      EndIf

      If(idiag>=1) Write(lun_diag,"(a,i6)") 'Recomputing symbolic decomposition at step ',kstep
      iwA = 0
      info = 0
      rinfo = 0.0
      Call MA48AD(ny,ny,lval,jobA,lia,vals,irn,jcn,keep,cntl,icntl,iwA,info,rinfo)
      If(info(1)/=0) Write(lun_diag,"(1x,a,i5,a,i2)") 'kstep=',kstep,', Warning during MA48AD.  Code: ',info(1)
    EndDo
    lia_min = info(4)
  Else
! Build sparse value array    
    vals(1:lval) = tvals
    If(scale_jac) vals(1:lval) = vals(1:lval)*mc29r(irn(1:lval))*mc29c(jcn(1:lval))
  EndIf

! Perform numerical decomposition using previously determined symbolic decomposition
  wB = 0.0
  iwB = 0
  info = 0
  rinfo = 0.0
  Call MA48BD(ny,ny,lval,jobB,lia,vals,irn,jcn,keep,cntl,icntl,wB,iwB,info,rinfo)
  If(info(1)==0 .or. ( info(1)==-3 .and. icntl(10)/=1)) lia_min = info(4)

  Select Case (info(1))
    Case(1)
      Write(lun_diag,"(1x,a,i5,a,i2)") 'kstep=',kstep,', Warning during MA48BD. Code: ',info(1)
      Write(lun_diag,"(1x,a)")         '  Switched from JOB=2 to JOB=1'
    Case(2)
      Write(lun_diag,"(1x,a,i5,a,i2)") 'kstep=',kstep,', Warning during MA48BD. Code: ',info(1)
      Write(lun_diag,"(1x,a,i5)")      '  Matrix is rank deficient.  info(5) = ',info(5)
    Case(3)
      Write(lun_diag,"(1x,a,i5,a,i2)") 'kstep=',kstep,', Warning during MA48BD. Code: ',info(1)
      Write(lun_diag,"(1x,a)")         '  Switched from JOB=2 to JOB=1'
      Write(lun_diag,"(1x,a,i5)")      '  Matrix is rank deficient.  info(5) = ',info(5)
    Case(:-1)
      Write(lun_diag,"(1x,a,i5,a,i2)") 'kstep=',kstep,', Error during MA48BD. Code: ',info(1)
  End Select

  If( info(6)>0 .or. info(1)/=0 ) Then
    jobB = 1 ! If entries have been dropped, the next call to MA48BD needs to be a normal call
  ElseIf( icntl(8)==0 ) Then
    jobB = 2 ! Otherwise, unless using special case of icntl(8)/=0, use the "fast" MA48BD call
  Else
    jobB = 3 ! For the case of icntl(8)/=0, use the "intermediate" MA48BD all
  EndIf
      
! Perform check to see if an error occured in MA48BD.  Most likely a singularity error
! caused by incorrect symbolic matrix, so run MA48AD again.  Also, make sure that 
! workspaces are large enough.
  Do kdecomp = 1,10
    If (info(1) == 0) Exit

    If (info(1) == -3 .or. lia_min > lia) Then
      Write(lun_diag,"(1x,a,i7,a,i2,a,i7)")'MA48 arrays too small, lia = ',lia,', info(4) = ',info(4),' lia_min =',lia_min
      Write(lun_diag,*)'Trying to reAllocate MA48 arrays'
      Deallocate(irn,jcn,vals)
      lia = max(int(1.2*lia_min),info(4))
      Allocate(irn(lia),jcn(lia),vals(lia))
    EndIf
      
! Rebuild vals array
    vals = 0.0
    irn = 0
    jcn = 0
    Do i=1,lval
      jcn(i) = cidx(i)
      irn(i) = ridx(i)
    EndDo
    vals(1:lval)=tvals

    If (scale_jac) Then
      vals(1:lval) = vals(1:lval)*mc29r(irn(1:lval))*mc29c(jcn(1:lval))
    ElseIf (kdecomp>2) Then
      scale_jac = .true.
      Write(lun_diag,*)'Trying to scale input matrix'
      ifail = 0
      Call MC29AD(ny,ny,lval,vals,irn,jcn,mc29r,mc29c,mc29w,lun_diag,ifail)
      If(ifail<0) Write(lun_diag,"(1x,a,i2)")'Error during MC29AD.  Code: ',ifail
      mc29r(1:ny) = exp(mc29r(1:ny))
      mc29c(1:ny) = exp(mc29c(1:ny))
      vals(1:lval) = vals(1:lval)*mc29r(irn(1:lval))*mc29c(jcn(1:lval))
    EndIf

    jobA = 1
    jobB = 1
    keep = 0
    iwA = 0
    info = 0
    rinfo = 0.0
    Call MA48AD(ny,ny,lval,jobA,lia,vals,irn,jcn,keep,cntl,icntl,iwA,info,rinfo)
    If(idiag>=1) Write(lun_diag,"(a,i6)") 'Recomputing symbolic decomposition at step ',kstep
    If(info(1)/=0) Write(lun_diag,"(1x,a,i5,a,i3,a,i2)") 'kstep=',kstep,',',kdecomp+1,': Warning during MA48AD.  Code: ',info(1)
    lia_min = info(4)

    If( info(1)==0 ) Then
      wB = 0.0
      iwB = 0
      info = 0
      rinfo = 0.0
      Call MA48BD(ny,ny,lval,jobB,lia,vals,irn,jcn,keep,cntl,icntl,wB,iwB,info,rinfo)
      If(info(1)==0 .or. ( info(1)==-3 .and. icntl(10)/=1)) lia_min = info(4)
      If(info(1)/=0) Write(lun_diag,"(1x,a,i5,a,i3,a,i2)") 'kstep=',kstep,',',kdecomp+1,': Warning during MA48BD.  Code: ',info(1)

      If( info(6)>0 .or. info(1)/=0 ) Then
        jobB = 1 ! If entries have been dropped, the next call to MA48BD needs to be a normal call
      ElseIf( icntl(8)==0 ) Then
        jobB = 2 ! Otherwise, unless using special case of icntl(8)/=0, use the "fast" MA48BD call
      Else
        jobB = 3 ! For the case of icntl(8)/=0, use the "intermediate" MA48BD all
      EndIf
    EndIf

  EndDo
    
! Stop timer
  stop_timer = xnet_wtime()
  timer_solve = timer_solve + stop_timer

  Return
End Subroutine jacobian_decomp                                                                       
  
Subroutine jacobian_bksub(rhs,dy) 
!===============================================================================
! This routine performs back-substitution for a previously factored matrix and 
! the vector rhs.   
!===============================================================================
  Use controls
  Use jacobian_data
  Use timers
  Real(8), Intent(inout) :: rhs(ny)
  Real(8), Intent(out) ::  dy(ny)
  Real(8) :: relerr(3)
  
! Initiate timer
  start_timer = xnet_wtime()
  timer_solve = timer_solve - start_timer  

  If(scale_jac) rhs(1:ny) = rhs(1:ny)*mc29r(1:ny)

! Perform back substitution 
  jobC = 2
  wC = 0.0
  iwC = 0
  relerr = 0.0
  info = 0
  Call MA48CD(ny,ny,trans,jobC,lia,vals,irn,keep,cntl,icntl,rhs,dy,relerr,wC,iwC,info)
  If(info(1)/=0) Write(lun_diag,"(1x,a,i2)") 'Warning during MA48CD.  Code: ',info(1)

! If the relative error becomes sufficiently large, regenerate data structures in the next call to MA48BD
  If( maxval(relerr)>yacc*maxerr ) jobB = 1

  If(scale_jac) dy(1:ny) = dy(1:ny)*mc29c(1:ny)
    
! Diagnostic output
  If(idiag>=4) Then
    Write(lun_diag,"(a)") 'BKSUB'
    Write(lun_diag,"(14es10.3)") dy
  EndIf

! Stop timer
  stop_timer = xnet_wtime()
  timer_solve = timer_solve + stop_timer

  Return
End Subroutine jacobian_bksub
