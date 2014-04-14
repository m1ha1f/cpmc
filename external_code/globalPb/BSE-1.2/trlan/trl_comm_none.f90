! $Id: trl_comm_none.f90,v 1.1 2004/05/04 06:29:35 fowlkes Exp $
!!!
! All communication routines used by TRLAN are in this file.
! This file contains routines to be used on sequential or shared memory
! environments.   No actual data exchange is performed.
!!!
! Initialization routine -- initializes a TRL_INFO variable.
! *********
! MUST BE CALLED before calling any other user level routine in TRLAN
! package
! *********
!
! Arguments:
! info   -- the information package to be used by TRLAN
! nrow   -- local dimension of the problem
! mxlan  -- maximum number of basis vectors to be used
! lh     -- which end of the spectrum to compute
!           <0 : lower end, the smallest eigenvalues
!           >0 : high end, the largest eigenvalues
!            0 : either lower and or high end, whoever converges
!                first 
! ned    -- number of wanted eigenvalues and eigenvectors
! tol    -- (optional) tolerance on residual norm,
!           default: sqrt(epsilon)
! trestart-- (optional) thick-restarting scheme, 1-4,
!           default: 1
! mxmv   -- (optional) maximum number of matrix-vector multiplications
!           allowed
!           default: info%ntot*info%ned
! mpicom -- (optional) the MPI communicator,
!           default: a duplicate of MPI_COMM_WORLD
!           in sequential case, this is set to 0.
!!!
Subroutine trl_init_info(info, nrow, mxlan, lh, ned, tol, trestart,&
     & mxmv, mpicom)
  Use trl_info
  Implicit None 
  Integer, Intent(in) :: nrow, mxlan, lh, ned
  Integer, Intent(in), Optional :: mpicom, mxmv, trestart
  Real(8), Intent(in), Optional :: tol
  Type(TRL_INFO_T), Intent(out) :: info
  !
  info%maxlan = mxlan
  If (mxlan <= ned) Then
     info%maxlan = ned + Max(ned, 6)
  End If
  info%lohi = lh
  info%ned = ned
  info%mpicom = -Huge(ned)
  info%nloc = nrow
  info%ntot = nrow
  info%guess = 0
  If (Present(mxmv)) Then
     info%maxmv = mxmv
  Else
     info%maxmv = Min(Max(info%ntot, 1000), 1000*info%ned)
  End If
  If (Present(trestart)) Then
     info%restart = trestart
  Else
     info%restart = 0
  End If
  If (Present(tol)) Then
     If (tol .Le. Tiny(tol)) Then
        info%tol = Epsilon(tol)
     Else If (tol .Ge. 1D0) Then
        info%tol = Min(0.1D0, 1D0/tol)
     Else
        info%tol = tol
     End If
  Else
     info%tol = Sqrt(Epsilon(info%tol))
  End If
  info%nec = 0
  info%locked = info%nec
  info%matvec = 0
  info%nloop = 0
  info%north = 0
  info%nrand = 0
  info%flop = 0
  info%rflp = 0
  info%flop_h = 0
  info%rflp_h = 0
  info%flop_r = 0
  info%rflp_r = 0
  info%clk_rate = 0
  info%clk_max = -1
  info%clk_tot = 0
  info%clk_op = 0
  info%clk_orth = 0
  info%clk_res = 0
  info%tick_t = 0
  info%tick_o = 0
  info%tick_h = 0
  info%tick_r = 0
  info%clk_in = 0
  info%clk_out = 0
  info%wrds_in = 0
  info%wrds_out = 0
  info%verbose = 0
  info%log_io = 99
  info%log_file = 'TRL_LOG_'
  info%stat = 0
  info%anrm = 0
  info%tmv = -1
  info%trgt = - Huge(info%anrm)
  info%tres = 0.0D0
  info%crat = 1.0D0
  info%my_pe = 0
  info%npes = 1
  info%cpflag = 0
  info%cpio = 98
  info%cpfile = 'TRL_CHECKPOINT_'
  info%oldcpf = ''
  Return
End Subroutine trl_init_info
!****** The rest are lower level working routines ******
!!!
! trl_g_sum performs global sum in the parallel environment, nothing
! is done here
!!!
Subroutine trl_g_sum(mpicom, nelm, x, y)
  Implicit None
  Integer, Intent(in) :: mpicom, nelm
  Real(8) :: x(nelm), y(*)
  Return
End Subroutine trl_g_sum
!!!
! trl_sync_flag
! given an integer value, this function returns the minimum value of
! all the PEs
!!!
Function trl_sync_flag(mpicom, inflag) Result(outflag)
  Implicit None
  Integer :: outflag
  Integer, Intent(in) :: mpicom, inflag
  outflag = inflag
End Function trl_sync_flag
!!!
! trl_g_dot implements a distributed version of BLAS routine dgemv
! which is used to compute dot-products by TRLAN
! this function performs (in matlab notation) wrk = [V1, V2]'*rr
!!!
! Arguments:
! mpicom  -- MPI communicator
! nrow    -- number of rows on the local processor
! v1      -- the first part of the matrix
! ld1     -- leading dimension of v1
! m1      -- number of columns in v1
! v2      -- the second part of the matrix
! ld2     -- the leading dimension of v2
! m2      -- number of columns in v2
! rr      -- the vector to be multiplied
! wrk     -- results of this operation.  !! size not checked !!
Subroutine trl_g_dot(mpicom, nrow, v1, ld1, m1, v2, ld2, m2, rr, wrk)
  Implicit None
  Integer, Intent(in) :: mpicom, nrow, ld1, ld2, m1, m2
  Real(8), Intent(in) :: v1(ld1,m1), v2(ld2,m2), rr(nrow)
  Real(8), Intent(out) :: wrk(m1+m2)
  !
  ! local variables
  Character, Parameter :: trans='T'
  Integer :: i, m1p1, nd
  Real(8), Parameter :: zero=0.0D0, one=1.0D0
  !
  ! BLAS routines used here
  External dgemv
  !
  nd = m1 + m2
  ! nothing to do if both m1 and m2 are zero
  If (nd.Le.0) Return
  ! make sure the array sizes are correct
  If (ld1.Lt.nrow .Or. ld2.Lt.nrow) Then
     Stop 'trl_g_dot: incorrect array sizes'
  End If
  m1p1 = m1 + 1
  If (m1 .Gt. 2) Then
     Call dgemv(trans, nrow, m1, one, v1, ld1, rr, 1, zero, wrk, 1)
  Else If (m1 .Eq. 2) Then
     wrk(1) = zero
     wrk(2) = zero
     Do i = 1, nrow
        wrk(1) = wrk(1) + v1(i,1)*rr(i)
        wrk(2) = wrk(2) + v1(i,2)*rr(i)
     End Do
  Else If (m1 .Eq. 1) Then
     wrk(1) = Dot_product(v1(1:nrow,1), rr(1:nrow))
  End If
  If (m2 .Gt. 2) Then
     Call dgemv(trans, nrow, m2, one, v2, ld2, rr, 1, zero,&
          & wrk(m1p1), 1)
  Else If (m2 .Eq. 2) Then
     wrk(m1p1) = zero
     wrk(nd) = zero
     Do i = 1, nrow
        wrk(m1p1) = wrk(m1p1) + v2(i,1) * rr(i)
        wrk(nd) = wrk(nd) + v2(i,2) * rr(i)
     End Do
  Else If (m2 .Eq. 1) Then
     wrk(m1p1) = Dot_product(v2(1:nrow,1), rr(1:nrow))
  End If
End Subroutine trl_g_dot
