!*******************************************************************************
! Jacobian Broadcast for PARDISO matrix solver, part of XNet 6, 1/31/12
!
! Needed for MPI execution with PARDISO matrix solver.
! This routine broadcasts the jacobian data between MPI tasks.
!
!*******************************************************************************

Subroutine jacobian_bcast(data_dir)
!===============================================================================
! This routine distributes sparse Jacobian data for PARDISO solver
!===============================================================================
  Use controls
  Use file_data
  Use jacobian_data
  Use mpi
  Character (LEN=80) , INTENT(IN):: data_dir  
  Integer :: sparse_int(4),i,ierr

! For PARDISO solver, parameters for the compressed row storage must be
! read and broadcast, along with allocations being performed 
! On PE0 ...  
  If(myid==0) Then
    Call read_jacobian_data(data_dir)

! Pack sparse array size passing arrays
    sparse_int(1)  = lval
    sparse_int(2)  = l1s
    sparse_int(3)  = l2s
    sparse_int(4)  = l3s

  EndIf

! Broadcast Jacobian array sizes
  Call mpi_bcast(sparse_int,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! Unpack Jacobian array sizes and allocate on PE /=0
  If(myid/=0) Then

! Unpack sparse array size passing arrays
    lval = sparse_int(1)
    l1s = sparse_int(2)
    l2s = sparse_int(3)
    l3s = sparse_int(4)
    lia = 8*lval

! Allocate sparse data arrays   
    Allocate(ridx(lval),cidx(lval),sident(lval))
    Allocate(ns11(l1s),ns21(l2s),ns22(l2s))
    Allocate(ns31(l3s),ns32(l3s),ns33(l3s))

  EndIf ! PE/=0      

! Broadcast Jacobian data arrays
  call mpi_bcast(ridx,lval,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(cidx,lval,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns11, l1s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns21, l2s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns22, l2s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns31, l3s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns32, l3s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ns33, l3s,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  
! Complete MA48 Initilization on PE /=0
  If(myid/=0) Then

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
!   icntl(4) = 0

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
!   icntl(6) = ny

! Set option to move each column for which iwA(j)=0 to the end of the pivot sequence within its block (default=0)
!   icntl(8) = 1

! Set option to automatically reset to job=1 if MA48BD fails for job=2
!   icntl(11) = 1

! Set the pivoting threshold (near 0.0 emphasizes sparsity; near 1.0 emphasizes stability)
!   cntl(2) = 1.0

! Allocate variables used by the MA48 solver
    Allocate(tvals(lval),vals(lia),jcn(lia),irn(lia))
    Allocate(wB(ny),wC(4*ny))
    Allocate(iwA(9*ny),iwB(4*ny),iwC(ny))
    If (icntl(8) == 0) Then
      Allocate(keep(max(ny/icntl(6),1)))
    Else
      Allocate(keep(6*ny+4*ny/icntl(6)+7))
    EndIf

! These are variables to be used by MC29 if scaling is deemed necessary
    Allocate(mc29r(ny),mc29c(ny),mc29w(5*ny))
    mc29r = 0.0d0
    mc29c = 0.0d0
    mc29w = 0.0d0
    scale_jac = .false.
!$OMP End PARALLEL

! Set the value for the maximum allowed error in the call to MA48CD
  maxerr = 1.0d-04

! Build a compressed row format version of the identity matrix
    Do i=1,lval
      If (ridx(i)==cidx(i)) sident(i)=1.0d0
    Enddo

  Endif ! PE/=0      

End Subroutine jacobian_bcast
