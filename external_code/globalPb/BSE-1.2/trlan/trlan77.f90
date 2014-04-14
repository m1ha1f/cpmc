! $Id: trlan77.f90,v 1.1 2004/05/04 06:29:35 fowlkes Exp $
!!!
!     a Fortran 77 wrapper for TRLAN
!
!     mapping of ipar to TRL_INFO_T
!  O  ipar(1)  = stat
! I   ipar(2)  = lohi
! I   ipar(3)  = ned
! IO  ipar(4)  = nec
! I   ipar(5)  = maxlan
! I   ipar(6)  = restart
! I   ipar(7)  = maxmv
! I   ipar(8)  = mpicom
! I   ipar(9)  = verbose
! I   iapr(10) = log_io
! I   ipar(11) = iguess
! I   ipar(12) = cpflag
! I   ipar(13) = cpio
! I   ipar(14) = mvflop
!  O  ipar(24) = locked
!  O  ipar(25) = matvec
!  O  ipar(26) = nloop
!  O  ipar(27) = north
!  O  ipar(28) = nrand
!  O  ipar(29) = total time in milliseconds
!  O  ipar(30) = MATVEC time in milliseconds
!  O  ipar(31) = re-orthogonalization time in milliseconds
!  O  ipar(32) = time spent in restarting
!
!     on Entry:
!     wrk(1) = tol
!     on Return:
!     wrk(1:nec) = residual norms
!     wrk(nec+1) = crat
!
! ** NOTE **
! Since TRLINFO is not available to the outside world, 'verbose' is
! interpreted slightly differently as in TRLAN.  In addition to its
! normal used in TRLAN, if 'verbose' is not a negative number, than
! the content of TRLINFO will be printed to the log files using
! TRL_PRINT_INFO.   The user should provide ipar(14) to indicate the
! number of floating-point operations performed on each PE while
! completing one matrix-vector multiplication (applying the operator on
! one vector).  This is used to compute MFLOPS rate of the program.
!!!
      Subroutine trlan77(op, ipar, nrow, mev, eval, evec, lde, wrk, &
     &     lwrk)
      Use trl_info
      Use trl_interface
      Implicit None
                                ! arguments
      Integer :: ipar(32), nrow, mev, lde, lwrk
      Real(8) :: eval(mev), evec(lde, mev), wrk(lwrk), tmp
      External op
                                ! local variables
      Type(TRL_INFO_T) :: trlinfo
      Integer :: j
                                ! executable statements
      If (lwrk.Le.mev) Then
         Print *, 'TRLANf77: should have at MEV+1 elements in wrk.'
         Return
      End If
                                ! initialize trlinfo
      Call trl_init_info(trlinfo, nrow, ipar(5), ipar(2), ipar(3), &
     &    wrk(1), ipar(6), ipar(7), ipar(8))
      Call trl_set_iguess(trlinfo, ipar(4), ipar(11))
      If (ipar(9) .Gt. 0) Then
         Call trl_set_debug(trlinfo, ipar(9), ipar(10))
      End If
      trlinfo%cpflag = ipar(12)
      trlinfo%cpio = ipar(13)
                                ! calling Fortran 90 version of TRLAN
      Call trlan(op, trlinfo, nrow, mev, eval, evec, lde, wrk, lwrk)
                                ! copy the parameters back to ipar
      ipar(1) = trlinfo%stat
      ipar(4) = trlinfo%nec
      ipar(24) = trlinfo%locked
      ipar(25) = trlinfo%matvec
      ipar(26) = trlinfo%nloop
      ipar(27) = trlinfo%north
      ipar(28) = trlinfo%nrand
      tmp = 1.0D3 / trlinfo%clk_rate
      ipar(29) = (trlinfo%tick_t+trlinfo%clk_tot)*tmp
      ipar(30) = (trlinfo%tick_o+trlinfo%clk_op)*tmp
      ipar(31) = (trlinfo%tick_h+trlinfo%clk_orth)*tmp
      ipar(32) = (trlinfo%tick_r+trlinfo%clk_res)*tmp
! print the summary information to the debug file if one is created
! else print summary information to the standard output
      Call trl_print_info(trlinfo, mvflop=ipar(14))
      If (ipar(4).Gt.0) Then
         j = ipar(4)
      Else
         j = Min(mev-1, ipar(3))
      End If
      Call trl_check_ritz(op, trlinfo, nrow, evec(:,1:j),&
           & eval, beta=wrk(1:j), wrk=wrk(j+1:))
!!
!! try one of the currection schemes
!!!$      Call trl_rayleigh_quotients(op, trlinfo, evec(:,1:j), wrk,&
!!!$           & evec(:,mev))
!      Call trl_ritz_projection(op, trlinfo, evec, wrk(1:j+j+2),&
!           & wrk(j+j+3:))
!      If (trlinfo%stat .Eq. 0) Then
!         Call trl_check_ritz(op, trlinfo, nrow, evec(:,1:j),&
!              & wrk(1:j), beta=wrk(j+1:j+j), wrk=wrk(j+j+1:))
!         wrk(1:j) = wrk(j+1:j+j)
!      End If
      wrk(j+1) = trlinfo%crat
      ipar(1) = trlinfo%stat
      Return
      End Subroutine trlan77
