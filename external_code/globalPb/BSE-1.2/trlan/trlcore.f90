! $Id: trlcore.f90,v 1.1 2004/05/04 06:29:35 fowlkes Exp $
!!!
!\Description:
! the actual work routine of restarted Lanczos program for real
! symmetric eigenvalue problems
!
! user may directly invoke this sunroutine but she/he is responsible
! for input correct arguments as needed
!
! 1) info needs to be initialized
! 2) if info%nec>0, the first nec elements of eval and first nec
!    columns of evec must contain valid eigenpairs
! 3) workspace of currect size
!    eval(mev)
!    evec(lde, mev) (lde >= nrow, mev >= ned)
!    base(ldb, info%maxlan-mev+1) (ldb>=nrow, not used if mev>maxlan)
!    wrk(lwrk) minimum amount of memory required by TRLANCZOS is
!    maxlan*(maxlan+10)
! 4) if log files are to be written, the user need to open files on IO
!    unit log_io so that the log gile may be written correctly.
! 5) the user must set the timing variable info%clk_tot and
!    info%clk_max using system_clock function call in order for this
!    subroutine to track timing results correctly
!
! the operator that defines the eigenvalue problem is expected to have
! the following interface
!  Subroutine OP(nrow, ncol, xin, ldx, yout, ldy)
!    Integer, Intent(in) :: nrow, ncol, ldx, ldy
!    Real(8), Dimension(ldx*ncol), Intent(in) :: xin
!    Real(8), Dimension(ldy*ncol), Intent(out) :: yout
!  End Subroutine OP
!
!\Algorithm
!  0. initialize input vector
!  1. do while (more eigenvalues to compute .and. more MATVEC allowed)
!  2.    first step
!     o   alpha(k+1) = dot_product(q_{k+1}, Aq_{k+1})
!     o   rr = A*q_{k+1}-alpha(k+1)*q_{k+1}-\sum_{i=1}^{k} beta(i)q_i
!     o   (re-orthogonalization)
!  3.    do j = k+2, m
!     o     rr = Aq_j
!     o     alpha(j) = dot_product(q_j, rr)
!     o     rr = rr - alpha(j)*q_j - beta(j-1)*q_{j-1}
!     o     (re-orthogonalization)
!        end do j = k+2, m
!  4.    restarting
!     o   call dstqrb to decompose the tridiagonal matrix
!     o   perform convergence test
!     o   determine what and how many Ritz pairs to save
!     o   compute the Ritz pairs to be saved
!     end do while
!
! The re-orthogonalization procedure is implemented in trl_orth.  it
! produces a normalized vector rr that is guaranteed to be orthogonal
! to the current basis.  An error will be reported if it can not
! achieve its goal.
!!!
Subroutine trlanczos(op, info, nrow, mev, eval, evec, lde, base, ldb,&
     & nbas, wrk, lwrk)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T) :: info
  Integer, Intent(in) :: nrow, ldb, lde, mev, nbas, lwrk
  Real(8) :: eval(mev)
  Real(8), Target :: evec(lde,mev), base(1:ldb,nbas), wrk(lwrk)
  Interface
     Subroutine op(nrow, ncol, xx, ldx, yy, ldy)
       Implicit None
       Integer, Intent(in) :: nrow, ncol, ldx, ldy
       Real(8), Intent(in) :: xx(ldx*ncol)
       Real(8), Intent(out) :: yy(ldy*ncol)
     End Subroutine op
     ! this routine normalizes rr and make sure it is orthogonal to all
     ! colums of v1 and v2
     Subroutine trl_orth(nrow, v1, ld1, m1, v2, ld2, m2, rr, kept,&
          & alpha, beta, wrk, lwrk, info)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T) :: info
       Integer, Intent(in) :: nrow, ld1, ld2, m1, m2, kept, lwrk
       Real(8), Intent(in), Target :: v1(ld1,m1), v2(ld2,m2)
       Real(8) :: rr(nrow), alpha(m1+m2), beta(m1+m2), wrk(lwrk)
     End Subroutine trl_orth
     Subroutine trl_check_orth(info, nrow, v1, ld1, j1, v2, ld2, j2, wrk, lwrk)
       Use trl_info
       Implicit None
       Integer, Intent(in) :: nrow, ld1, ld2, j1, j2, lwrk
       Type(TRL_INFO_T), Intent(in) :: info
       Real(8), Intent(in) :: v1(ld1,j1), v2(ld2,j2)
       Real(8) :: wrk(lwrk)
     End Subroutine trl_check_orth
     ! compute the errors in the Lanczos recurrence
     Subroutine trl_check_recurrence(op, info, nrow, v1, ld1, m1, v2,&
          & ld2, m2, kept, alpha, beta, wrk, lwrk)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: nrow, ld1, ld2, m1, m2, kept, lwrk
       Real(8), Intent(in), Target :: v1(ld1,m1), v2(ld2,m2)
       Real(8), Intent(in) :: alpha(m1+m2-1), beta(m1+m2-1)
       Real(8), Target :: wrk(lwrk)
       Interface
          Subroutine op(nrow, ncol, xx, ldx, yy, ldy)
            Implicit None
            Integer, Intent(in) :: nrow, ncol, ldx, ldy
            Real(8), Intent(in) :: xx(ldx*ncol)
            Real(8), Intent(out) :: yy(ldy*ncol)
          End Subroutine op
       End Interface
     End Subroutine trl_check_recurrence
     Subroutine trl_initial_guess(nrow, evec, lde, mev, base, ldb, nbas,&
          & alpha, beta, j1, j2, info, wrk, lwrk)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T) :: info
       Integer, Intent(in) :: nrow, lde, ldb, mev, nbas, lwrk
       Integer, Intent(out) :: j1, j2
       Real(8) :: wrk(lwrk), base(ldb,nbas), evec(lde,mev)
       Real(8) :: alpha(mev+nbas-1), beta(mev+nbas-1)
     End Subroutine trl_initial_guess
     Subroutine trl_g_sum(mpicom, n, xx, wrk)
       Implicit None
       Integer, Intent(in) :: mpicom, n
       Real(8), Intent(out) :: wrk(n)
       Real(8), Intent(inout) :: xx(n)
     End Subroutine trl_g_sum
     Subroutine trl_tridiag(nd, alpha, beta, rot, alfrot, betrot, wrk,&
          & lwrk, ierr)
       Implicit None
       Integer, Intent(in) :: nd, lwrk
       Integer, Intent(out) :: ierr
       Real(8), Intent(in), Dimension(nd) :: alpha, beta
       Real(8) :: rot(nd*nd), alfrot(nd), betrot(nd), wrk(lwrk)
     End Subroutine trl_tridiag
     Subroutine trl_get_eval(nd, locked, alpha, beta, lambda, res, wrk,&
          & lwrk, ierr)
       Implicit None
       Integer, Intent(in) :: nd, locked, lwrk
       Integer :: ierr
       Real(8), Intent(in) :: alpha(nd), beta(nd)
       Real(8), Intent(out) :: lambda(nd), res(nd), wrk(lwrk)
     End Subroutine trl_get_eval
     Subroutine trl_convergence_test(nd, lambda, res, info, wrk)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T) :: info
       Integer, Intent(in) :: nd
       Real(8), Intent(in) :: lambda(nd), res(nd)
       Real(8), Intent(out) :: wrk(nd+nd)
     End Subroutine trl_convergence_test
     Subroutine trl_print_progress(info)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
     End Subroutine trl_print_progress
     Subroutine trl_shuffle_eig(nd, mnd, lambda, res, info, kept)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: nd, mnd
       Integer :: kept
       Real(8) :: lambda(nd), res(nd)
     End Subroutine trl_shuffle_eig
     Subroutine trl_get_tvec(nd, alpha, beta, irot, nrot, rot, nlam,&
          & lambda, yy, iwrk, wrk, lwrk, ierr)
       Implicit None
       Integer, Intent(in) :: nd, nrot, nlam, irot, lwrk
       Integer :: ierr, iwrk(nd,4)
       Real(8), Intent(in) :: alpha(nd), beta(nd), lambda(nlam), rot(nrot*nrot)
       Real(8) :: wrk(lwrk), yy(nd*nlam)
     End Subroutine trl_get_tvec
     Subroutine trl_get_tvec_a(nd, kept, alpha, beta, nlam, lambda,&
          & yy, wrk, lwrk, iwrk, ierr)
       Implicit None
       Integer, Intent(in) :: kept, nd, nlam, lwrk
       Integer :: ierr, iwrk(nd)
       Real(8), Intent(in) :: lambda(nd)
       Real(8) :: alpha(nd), beta(nd), yy(nd*nd), wrk(lwrk)
     End Subroutine trl_get_tvec_a
     Subroutine trl_ritz_vectors(nrow, lck, ny, yy, ldy, vec1, ld1, m1,&
          & vec2, ld2, m2, wrk, lwrk)
       Implicit None
       Integer, Intent(in) :: nrow, ny, lck, ldy, ld1, m1, m2, ld2, lwrk
       Real(8), Intent(in) :: yy(ldy, ny)
       Real(8) :: vec1(ld1,m1), vec2(ld2,m2), wrk(lwrk)
     End Subroutine trl_ritz_vectors
  End Interface
  !
  ! local variables
  !
  Character, Parameter :: notrans = 'N'
  Real(8), Parameter :: zero=0.0D0, one=1.0D0
  Character(132) :: title
  Integer :: i1, i2, j1, j2, jnd, jml, j1n, j2n, kept, prek
  Integer :: next_test, clk1, lwrk2, chkpnt, locked
  Integer, Dimension(:), Pointer :: iwrk
  Real(8), Dimension(:), Pointer :: alpha, beta, rr, rot,&
       & alfrot, betrot, lambda, res, yy, qa, qb, wrk2
  External dgemv
  !
  ! alpha, beta, alfrot and betrot have fixed locations in wrk
  !
  Nullify(alpha, beta, rr, rot, alfrot, betrot, lambda, res, yy, qa,&
       & qb, wrk2)
  title = ''
  clk1 = 0
  i1 = info%maxlan + 1
  i2 = info%maxlan + info%maxlan
  alpha => wrk(1:info%maxlan)
  beta  => wrk(i1:i2)
  i1 = i2 + 1
  i2 = i2 + info%maxlan
  alfrot => wrk(i1:i2)
  i1 = i2 + 1
  i2 = i2 + info%maxlan
  betrot => wrk(i1:i2)
  Allocate(iwrk(info%maxlan*4))
  iwrk = 0
  If (info%cpflag.Le.0) Then
     chkpnt = info%maxmv + info%maxlan
  Else
     chkpnt = info%maxmv / info%cpflag
  End If
  locked = info%nec
  !
  ! assign values to alpha, beta
  ! uses the assumption that the content of eval(1:nec) are eigenvalues
  ! and their residual norms are zero
  !
  locked = info%nec
  If (locked.Gt.0) Then
     alpha(1:locked) = eval(1:locked)
     beta(1:locked) = zero
  End If
  ! get valid initial guess for the Lanczos iterations
  wrk2 => wrk(i2+1:lwrk)
  lwrk2 = lwrk - i2
  Call trl_initial_guess(nrow, evec, lde, mev, base, ldb, nbas, alpha,&
       & beta, j1, j2, info, wrk2, lwrk2)
  If (info%stat .Ne. 0) Goto 888
  jnd = j1 + j2
  kept = jnd
  ! we will perform the first convergence test after next_test
  ! matrix-vector multiplications
  i1 = info%ned - jnd
  next_test = i1 + Min(i1, 6, info%ned/2)
  !!*************************************************!!
  !! the restart loop --                             !!
  !! restart if there is more eigenvalues to compute !!
  !! and more MATVEC allowed to do it                !!
  !!*************************************************!!
  Do While (info%matvec.Lt.info%maxmv .And. info%nec.Lt.info%ned)
     jnd = jnd + 1
     i2 = lwrk - (jnd-locked)**2
     rot => wrk(i2+1:lwrk)
     i1 = 4*info%maxlan + 1
     wrk2 => wrk(i1:i2)              ! workspace for orthogonalization
     lwrk2 = i2 - i1 + 1
     i1 = Max(5*info%maxlan-3*locked, 4*info%maxlan)
     If (lwrk2 .Lt. i1) Then
        info%stat = -11
        Return
     End If
     !
     ! first step
     !
     If (j1 .Lt. mev) Then
        j1 = j1 + 1
        qa => evec(1:nrow, j1)
     Else
        j2 = j2+1
        qa => base(1:nrow, j2)
     End If
     If (j1 .Lt. mev) Then
        j1n = j1 + 1
        j2n = 0
     Else
        j1n = mev
        j2n = j2 + 1
     End If
     Nullify(qb)
     If (j1.Lt.mev) Then
        rr => evec(1:nrow, j1n)
     Else
        rr => base(1:nrow, j2n)
     End If
     ! perform matrix-vector multiplication
     ! record the total time and the time in MATVEC
     Call System_clock(clk1)
     If (clk1 .Le. info%clk_tot) Then
        info%tick_t = info%tick_t + (info%clk_max-info%clk_tot)
        info%tick_t = info%tick_t + clk1
        info%clk_tot = clk1
     End If
!!$print *, 'input to OP #', info%matvec + 1
!!$print *, qa
!!$print *, 'should be the same as'
!!$if (j2 .eq. 0) then
!!$print *, evec(1:nrow,j1)
!!$else
!!$print *, base(1:nrow,j2)
!!$end if
     Call OP(nrow, 1, qa, nrow, rr, nrow)
!!$print *, 'return of OP #', info%matvec + 1
!!$print *, rr
     Call add_clock_ticks(info%clk_op, info%tick_o)
!!$     Call hpcheck(info%stat)
!!$     Print *, 'hpcheck (trlcore.f90:298) returns ', info%stat,&
!!$          & ', nloop=', info%nloop, ', nec=', info%nec, ', jnd=',jnd
     info%matvec = info%matvec+1
     ! perform the Lanczos orthogonalization
     alpha(jnd) = Dot_product(qa, rr)
     Call trl_g_sum(info%mpicom, 1, alpha(jnd:jnd), wrk2)
     beta(jnd) = alpha(jnd)
     info%flop = info%flop + nrow + nrow
     If (j1 .Gt. 2) Then
        Call dgemv(notrans, nrow, j1, -one, evec, lde,&
             & beta, 1, one, rr, 1)
        info%flop = info%flop + 2*j1*nrow
     Else If (j1 .Eq. 1) Then
        rr = rr - alpha(1)*qa
        info%flop = info%flop + nrow + nrow
     Else If (j1 .Eq. 2) Then
        rr = rr - beta(1)*evec(1:nrow,1) - beta(2)*evec(1:nrow,2)
        info%flop = info%flop + 4*nrow
     End If
     If (j2 .Gt. 2) Then
        Call dgemv(notrans, nrow, j2, -one, base, ldb,&
             & beta(j1+1), 1, one, rr, 1)
        info%flop = info%flop + 2*j2*nrow
     Else If (j2 .Eq. 1) Then
        rr = rr - beta(jnd)*qa
        info%flop = info%flop + nrow + nrow
     Else If (j2 .Eq. 2) Then
        rr = rr - beta(j1+1)*base(1:nrow,1) - beta(jnd)*base(1:nrow,2)
        info%flop = info%flop + 4*nrow
     End If
     ! perform re-orthogonalization
     info%flop_h = info%flop_h - info%flop
     Call System_clock(clk1)
     Call trl_orth(nrow, evec, lde, j1, base, ldb, j2, rr,&
          & kept, alpha, beta, wrk2, lwrk2, info)
     If (info%verbose.Gt.8 .Or. info%stat.Ne.0) Then
        ! check orthogonality after the initilization step
        Call trl_check_orth(info, nrow, evec, lde, j1n, base, ldb, j2n,&
             & wrk2, lwrk2)
     End If
     Call add_clock_ticks(info%clk_orth, info%tick_h)
     info%flop_h = info%flop_h + info%flop
     If (info%stat .Ne. 0) Goto 888
     If (info%verbose .Gt. 5) Call print_alpha_beta
     !
     ! transform the matrix formed by alpha and beta into a
     ! tridiagonal matrix, rot stores the transformation matrix
     !
     alfrot(1:locked) = alpha(1:locked)
     betrot(1:locked) = zero
     i1 = jnd-locked
     Call trl_tridiag(i1, alpha(locked+1:locked+i1), beta(locked+1:locked+i1),&
          & rot, alfrot(locked+1:locked+i1), betrot(locked+1:locked+i1), wrk2,&
          & lwrk2, info%stat)
     info%flop = info%flop + 8*i1*i1*i1/3  ! Golub:1996:MC, P415
     If (info%stat .Ne. 0) Goto 888
     betrot(jnd) = beta(jnd)
     !!******************************************************!!
     !! regular iterations of Lanczos algorithm (inner loop) !!
     !!******************************************************!!
     Do While (jnd .Lt. info%klan .And. info%nec .Lt. info%ned)
        qb => qa
        qa => rr
        j1 = j1n
        j2 = j2n
        jnd = jnd + 1
        If (j1n .Lt. mev) Then
           j1n = j1n + 1
        Else
           j2n = j2n + 1
        End If
        If (j1.Lt.mev) Then
           rr => evec(1:nrow, j1n)
        Else
           rr => base(1:nrow, j2n)
        End If
        ! MATVEC
        Call System_clock(clk1)
        Call OP(nrow, 1, qa, nrow, rr, nrow)
        Call add_clock_ticks(info%clk_op, info%tick_o)
        info%matvec = info%matvec+1
        ! compute alpha(jnd)
        alpha(jnd) = Dot_product(qa, rr)
        ! the Lanczos orthogonalization
        Call trl_g_sum(info%mpicom, 1, alpha(jnd:jnd), wrk2)
        rr = rr - alpha(jnd)*qa - beta(jnd-1)*qb
        info%flop = info%flop + 6*nrow
        ! re-orthogonalization, compute beta(jnd)
        info%flop_h = info%flop_h - info%flop
        Call System_clock(clk1)
        Call trl_orth(nrow, evec, lde, j1, base, ldb, j2, rr, kept,&
             & alpha, beta, wrk2, lwrk2, info)
        If (info%verbose.Gt.8 .Or. info%stat.Ne.0) Then
           ! check orthogonality after the initilization step
           Call trl_check_orth(info, nrow, evec, lde, j1n, base, ldb, j2n,&
                & wrk2, lwrk2)
        End If
        Call add_clock_ticks(info%clk_orth, info%tick_h)
        info%flop_h = info%flop_h + info%flop
        If (info%stat.Ne.0) Goto 888
        alfrot(jnd) = alpha(jnd)
        betrot(jnd) = beta(jnd)
        If (info%verbose .Gt. 4) Call print_alpha_beta
        !
        ! perform convergence test once in a while
        If (info%matvec.Ge.next_test .And. (beta(jnd).Ne.zero .Or. &
             & info%matvec.Gt.info%maxlan)) Then
           If (info%verbose .Gt. 5) Call print_all_alpha_beta
           lambda => wrk2(1:jnd)
           res => wrk2(jnd+1:jnd+jnd)
           Call trl_get_eval(jnd, locked, alfrot, betrot, lambda, res,&
                & wrk2(jnd+jnd+1:lwrk2), lwrk2-jnd-jnd, info%stat)
           If (info%stat .Ne. 0) Goto 888
           If (info%verbose .Gt. 2) Call print_lambda_res
           i1 = Min(mev, jnd)
           eval(1:i1) = wrk2(1:i1)
           Call trl_convergence_test(jnd, lambda, res, info,&
                & wrk2(jnd+jnd+1:4*jnd))
           ! decide when to perform the next test
           If (info%nec .Lt. info%ned .And. info%nec.Gt.0) Then
              next_test = Dble(info%ned * info%matvec) / Dble(info%nec)
           Else If (info%nec.Eq.0) Then
              next_test = next_test + next_test
              If (info%maxlan .Eq. info%ntot) Then
                 next_test = Ceiling(0.5*(info%maxlan + info%matvec))
              End If
           End If
           If (info%verbose .Gt. 0) Call trl_print_progress(info)
        End If
     !!**********************************************************!!
     End Do ! inner (regular Lanczos three-term recurrence) loop !!
     !!**********************************************************!!
     !
     ! error checking for debugging use
     !
     lambda => wrk2(1:jnd)
     res => wrk2(jnd+1:jnd+jnd)
     If (info%verbose.Gt.6) Then
        wrk2 => wrk2(jnd+jnd+1:lwrk2)
        i2 = lwrk2 - jnd - jnd
        Call trl_check_orth(info, nrow, evec, lde, j1n, base, ldb, j2n,&
             & wrk2, i2)
        If (info%verbose.Gt.7) Then
           Call trl_check_recurrence(op, info, nrow, evec, lde, j1n,&
                & base, ldb, j2n, kept, alpha, beta, wrk2, i2)
        End If
     End If
     !
     ! convert the integer counters to floating-point counters
     !
     i2 = info%clk_max / 4
     If (info%flop .Gt. i2) Then
        info%rflp = info%rflp + info%flop
        info%flop = 0
     End If
     If (info%flop_h .Gt. i2) Then
        info%rflp_h = info%rflp_h + info%flop_h
        info%flop_h = 0
     End If
     If (info%flop_r .Gt. i2) Then
        info%rflp_r = info%rflp_r + info%flop_r
        info%flop_r = 0
     End If
     If (info%clk_op .Gt. i2) Then
        info%tick_o = info%tick_o + info%clk_op
        info%clk_op = 0
     End If
     If (info%clk_orth .Gt. i2) Then
        info%tick_h = info%tick_h + info%clk_orth
        info%clk_orth = 0
     End If
     If (info%clk_res .Gt. i2) Then
        info%tick_r = info%tick_r + info%clk_res
        info%clk_res = 0
     End If
     info%flop_r = info%flop_r - info%flop
     !
     ! *** Determine whether to restart ***
     ! compute the Ritz values and Ritz vectors if they are not up to
     ! date
     !
     Call System_clock(clk1)
     prek = kept
     jml = jnd - locked
     i2 = kept - locked + 1
     If (info%nec .Lt. info%ned) Then
        ! need to compute the updated Ritz values and residual norms
        wrk2 => wrk(4*info%maxlan+2*jnd+1:lwrk-i2*i2)
        lwrk2 = lwrk-i2*i2-4*info%maxlan-2*jnd
        If (lwrk2 .Lt. 3*jnd) Then
           info%stat = -12
           Goto 888
        End If
        If (info%verbose .Gt. 5) Call print_all_alpha_beta
        Call trl_get_eval(jnd, locked, alfrot, betrot,&
             & lambda, res, wrk2, lwrk2, info%stat)
        If (info%stat .Ne. 0) Goto 888
        If (info%verbose .Gt. 2) Call print_lambda_res
        Call trl_convergence_test(jnd, lambda, res, info, wrk2)
        If (info%verbose .Gt. 0) Call trl_print_progress(info)
     End If
     !
     ! compute the Ritz vectors
     ! decide how many vectors to save if restart
     If (info%nec.Lt.info%ned .And. info%matvec.Lt.info%maxmv) Then
        ! prepare to restart
        Call trl_shuffle_eig(jml, info%klan-locked, lambda(locked&
             &+1:locked+jml), res(locked+1:locked+jml), info, kept)
        ! compute eigenvectors using dstein (inverse interations)
        If (kept*3 .Lt. jml) Then
           i1 = 4*info%maxlan+kept*jml+jnd
           yy => wrk(4*info%maxlan+jnd+1:i1)
           wrk2 => wrk(i1+1:lwrk-i2*i2)
           lwrk2 = lwrk - i1 - i2*i2
           Call trl_get_tvec(jml, alfrot(locked+1:locked+jml),&
                & betrot(locked+1:locked+jml), 0, i2, rot, kept,&
                & lambda(locked+1:locked+jml), yy, iwrk, wrk2, lwrk2,&
                & info%stat)
           If (info%stat.Eq.0)&
                & alpha(1:locked+kept) = lambda(1:locked+kept)
        End If
        ! compute eigenvectors using dsyev (QR)
        If (kept*3.Ge.jml .Or. info%stat.Ne.0) Then
!0918        If (info%stat.Ne.0) Then
           alfrot(1:locked+kept) = lambda(1:locked+kept)
           i1 = 4*info%maxlan+jml*jml
           yy => wrk(4*info%maxlan+1:i1)
           wrk2 => wrk(i1+1:lwrk)
           lwrk2 = lwrk - i1
           Call trl_get_tvec_a(jml, prek-locked, alpha(locked+1:locked&
                &+jml), beta(locked+1:locked+jml), kept, alfrot(locked&
                &+1:locked+jml), yy, wrk2, lwrk2, iwrk, info%stat)
        End If
        If (info%stat .Ne. 0) Goto 888
        beta(locked+1:locked+kept) = yy(jml:jml*kept:jml)*betrot(jnd)
        If (jml.Gt.info%ned+(info%ned/5+6)) Then
           Call trl_set_locking(jml, kept, alpha(locked+1:locked+kept),&
                & beta(locked+1:locked+kept), yy, info%anrm, i2)
        Else
           i2 = 0
        End If
        !
        ! generate Ritz vectos, reclaim the space pointed by ROT
        i1 = 4*info%maxlan + kept*jml + jnd
        wrk2 => wrk(i1+1:lwrk)
        lwrk2 = lwrk - i1
        Call trl_ritz_vectors(nrow, locked, kept, yy, jml,&
             & evec, lde, j1, base, ldb, j2, wrk2, lwrk2)
        info%flop = info%flop + 2*nrow*jml*kept
        If (info%verbose .Gt. 0) Call print_restart_state
        ! reset the counters and indices to the correct values for
        ! restarting
        kept = kept + locked
        locked = locked + i2
        info%locked = locked
        jnd = kept
        If (jnd .Le. mev) Then
           j1 = jnd
           j2 = 0
        Else
           j1 = mev
           j2 = jnd - mev
        End If
        If (info%nec .Gt. 0) Then
           next_test = Int(Dble(info%matvec*info%ned)/Dble(info%nec))
        Else
           next_test = next_test + info%maxlan
        End If
        i1=Min(mev, jnd)
        eval(1:i1) = lambda(1:i1)
        If (jnd .Lt. mev) Then
           j1n = j1 + 1
           j2n = 0
           evec(1:nrow, j1n) = rr
        Else
           j1n = mev
           j2n = j2 + 1
           base(1:nrow, j2n) = rr
        End If
        ! write checkpoint files
        If (info%matvec .Ge. chkpnt) Then
           Call write_checkpoint
           chkpnt = chkpnt + info%maxmv/info%cpflag
        End If
     Else
        !
        ! all wanted eigenpairs converged or maximum MATVEC used
        ! sort the eigenvalues in final output order
        kept = Min(info%nec, Max(info%ned,mev-1))
        info%nec = kept
        If (kept.Eq.0) kept = Min(mev-1, info%ned)
        Call trl_sort_eig(jnd, info%lohi, kept, lambda, res)
        eval(1:kept) = lambda(1:kept)
        If (kept*3 .Lt. jnd) Then
           ! eigenvectors of the projection matrix (try inverse
           ! interations)
           i1 = kept*jnd + 4*info%maxlan
           yy => wrk(4*info%maxlan+1:i1)
           wrk2 => wrk(i1+1:lwrk-i2*i2)
           lwrk2 = lwrk - i1 - i2*i2
           Call trl_get_tvec(jnd, alfrot, betrot, locked, i2, rot,&
                & kept, eval, yy, iwrk, wrk2, lwrk2, info%stat)
        End If
        If (kept*3.Ge.jnd .Or. info%stat .Ne. 0) Then
           ! too many eigenvectors or inverse iterations have failed,
           ! try QR
           i1 = 4*info%maxlan+jnd*jnd
           yy => wrk(4*info%maxlan+1:i1)
           wrk2 => wrk(i1+1:lwrk)
           lwrk2 = lwrk - i1
           Call trl_get_tvec_a(jnd, prek, alpha, beta, kept, eval, yy,&
                & wrk2, lwrk2, iwrk, info%stat)
           If (info%stat .Ne. 0) Goto 888
        End If
        alpha(1:kept) = eval(1:kept)
        beta(1:kept) = betrot(jnd)*yy(jnd:kept*jnd:jnd)
        !
        ! generate eigenvectos, reclaim the space pointed by ROT
        i1 = kept*jnd + 4*info%maxlan
        wrk2 => wrk(i1+1:lwrk)
        lwrk2 = lwrk - i1
        Call trl_ritz_vectors(nrow, 0, kept, yy, jnd,&
             & evec, lde, j1, base, ldb, j2, wrk2, lwrk2)
        info%flop = info%flop + 2*nrow*jml*kept
        If (info%verbose .Gt. 1) Call print_final_state
        ! reset the counters and indices to be used by check_orth and
        ! check_recurrence
        jnd = kept
        j1 = kept
        j2 = 0
        If (kept .Lt. mev) Then
        End If
        If (j1 .Lt. mev) Then
           j1n = j1 + 1
           j2n = 0
           evec(1:nrow, j1n) = rr
        Else
           j1n = mev
           j2n = 1
           base(1:nrow, 1) = rr
        End If
        ! write checkpoint files
        If (info%cpflag .Gt. 0) Call write_checkpoint
     End If
     !
     ! check the orthogonality of the basis vectors before restarting
     !
     If (info%verbose.Gt.6) Then
        Call trl_check_orth(info, nrow, evec, lde, j1n, base, ldb, j2n,&
             & wrk2, lwrk2)
        If (info%verbose.Gt.7) Then
           Call trl_check_recurrence(op, info, nrow, evec, lde, j1n,&
                & base, ldb, j2n, kept, alpha, beta, wrk2, lwrk2)
        End If
     End If
     Call add_clock_ticks(info%clk_res, info%tick_r)
     info%flop_r = info%flop_r + info%flop
     info%nloop = info%nloop + 1
  End Do
  !!*********************!!
  !! end of restart_loop !!
  !!*********************!!
  ! write the estimated residual norms to the beginning of WRK
  wrk(1:j1) = Abs(beta(1:j1))
888 If (info%stat.Lt.0 .And. (info%verbose.Gt.0.Or.info%my_pe.Eq.0))&
       & Call log_error_state
  Deallocate(iwrk)
  Return
  ! internal subroutines for printing etc
Contains
  ! add clock ticks to a integer variable, if there is potential for
  ! overflow convert it to floating-point numbers
  Subroutine add_clock_ticks(time, rtime)
    Integer :: time, clk2
    Real(8) :: rtime
    Call System_clock(clk2)
    If (clk2 .Ge. clk1) Then
       clk2 = clk2 - clk1
    Else
       clk2 = clk2 + (info%clk_max - clk1)
    End If
    If (clk2+time .Ge. time) Then
       time = clk2 + time
    Else
       rtime = (rtime + time) + clk2
       time = 0
    End If
  End Subroutine add_clock_ticks
  ! print the current alpha(i) and beta(i)
  Subroutine print_alpha_beta
    Interface
       ! print real array for debugging
       Subroutine trl_print_real(info, title, array)
         Use trl_info
         Implicit None
         Type(TRL_INFO_T), Intent(in) :: info
         Character(*), Intent(in) :: title
         Real(8), Dimension(:), Intent(in) :: array
       End Subroutine trl_print_real
    End Interface
    Write(title, *) 'alpha(',jnd,') ='
    Call trl_print_real(info, title, alpha(jnd:jnd))
    Write(title, *) ' beta(',jnd,') ='
    Call trl_print_real(info, title, beta(jnd:jnd))
  End Subroutine print_alpha_beta
  ! print all computed alpha and beta
  Subroutine print_all_alpha_beta
    Interface
       ! print real array for debugging
       Subroutine trl_print_real(info, title, array)
         Use trl_info
         Implicit None
         Type(TRL_INFO_T), Intent(in) :: info
         Character(*), Intent(in) :: title
         Real(8), Dimension(:), Intent(in) :: array
       End Subroutine trl_print_real
    End Interface
    Write(title, *) 'alfrot(1:', jnd, ')..'
    Call trl_print_real(info, title, alfrot(1:jnd))
    title(1:3) = 'bet'
    Call trl_print_real(info, title, betrot(1:jnd))
  End Subroutine print_all_alpha_beta
  ! print lambda and residual norm
  Subroutine print_lambda_res
    Interface
       ! print real array for debugging
       Subroutine trl_print_real(info, title, array)
         Use trl_info
         Implicit None
         Type(TRL_INFO_T), Intent(in) :: info
         Character(*), Intent(in) :: title
         Real(8), Dimension(:), Intent(in) :: array
       End Subroutine trl_print_real
    End Interface
    title = 'Current eigenvalues..'
    Call trl_print_real(info, title, lambda)
    title = 'Current residual norms..'
    Call trl_print_real(info, title, res)
  End Subroutine print_lambda_res
  ! print current solution status
  Subroutine print_restart_state
    Interface
       ! print integer array for debugging
       Subroutine trl_print_int(info, title, array)
         Use trl_info
         Implicit None
         Type(TRL_INFO_T), Intent(in) :: info
         Character(*), Intent(in) :: title
         Integer, Dimension(:), Intent(in) :: array
       End Subroutine trl_print_int
       ! print real array for debugging
       Subroutine trl_print_real(info, title, array)
         Use trl_info
         Implicit None
         Type(TRL_INFO_T), Intent(in) :: info
         Character(*), Intent(in) :: title
         Real(8), Dimension(:), Intent(in) :: array
       End Subroutine trl_print_real
    End Interface
    iwrk(1) = kept+locked
    iwrk(2) = locked+i2
    title = 'Number of saved and locked Ritz pairs..'
    Call trl_print_int(info, title, iwrk(1:2))
    If (info%verbose.Gt.2) Then
       If (iwrk(2) .Eq. 0) Then
          title = 'Ritz values saved (ascending ordered)..'
       Else
          title = 'Ritz values saved (may not be ordered)..'
       End If
       Call trl_print_real(info, title, alpha(1:kept+locked))
       title = 'Residual norms of the saved Ritz pairs..'
       betrot(1:kept+locked) = Abs(beta(1:kept+locked))
       Call trl_print_real(info, title, betrot(1:kept+locked))
    End If
    If (info%verbose .Gt. 7) Then
       Do j1 = 1, Min(kept, info%verbose)
          Do j2 = 1, j1
             wrk2(j2) = Dot_product(yy(j2*jml-jml+1:j2*jml),&
                  & yy(j1*jml-jml+1:j1*jml))
          End Do
          wrk2(j1) = wrk2(j1) - 1
          Write(title, *) 'Orthogonality level of y(', j1, ') ..'
          Call trl_print_real(info, title, wrk2(1:j1))
       End Do
    End If
    If (info%verbose .Gt. 10) Then
       Do j1 = 1, Min(kept, info%verbose)
          Write(title,*) 'eigenvector ', j1, ' of Q''AQ..'
          Call trl_print_real(info, title, yy((j1-1)*jml+1:j1*jml))
       End Do
    End If
    If (info%verbose .Gt. 10) Then
       j1n = Min(nrow, info%verbose)
       Do j1 = 1, Min(kept+locked, mev)
          Write(title,*) 'Ritz vector ', j1, ' (1:', j1n, ') ..'
          Call trl_print_real(info, title, evec(1:j1n,j1))
       End Do
    End If
  End Subroutine print_restart_state
  ! print the final state
  Subroutine print_final_state
    Interface
       ! print real array for debugging
       Subroutine trl_print_real(info, title, array)
         Use trl_info
         Implicit None
         Type(TRL_INFO_T), Intent(in) :: info
         Character(*), Intent(in) :: title
         Real(8), Dimension(:), Intent(in) :: array
       End Subroutine trl_print_real
    End Interface
    title = 'Final eigenvalues  (in ascending order)..'
    Call trl_print_real(info, title, eval(1:kept))
    If (info%verbose .Gt. 4) Then
       title = 'Final residual norms..'
       Call trl_print_real(info, title, beta(1:kept))
    End If
    If (info%verbose .Gt. 8) Then
       Do j1 = 1, Min(kept, info%verbose)
          Write(title,*) 'Eigenvector ', j1, ' of Q''AQ ..'
          Call trl_print_real(info, title, yy((j1-1)*jml+1:j1*jml))
       End Do
    End If
    If (info%verbose .Gt. 10) Then
       j1n = Min(nrow, info%verbose)
       Do j1 = 1, Min(kept, mev)
          Write(title,*) 'Ritz vector ', j1, ' (1:', j1n, ') ..'
          Call trl_print_real(info, title, evec(1:j1n,j1))
       End Do
    End If
  End Subroutine print_final_state
  ! do check pointing
  Subroutine write_checkpoint
    Interface
       ! write the checkpoint files
       Subroutine trl_write_checkpoint(cp_io, filename, nrow, alpha,&
            & beta, evec, lde, me, base, ldb, nb, ierr)
         Implicit None
         Character(*), Intent(in) :: filename
         Integer, Intent(in) :: cp_io, nrow, ldb, lde, me, nb
         Real(8), Intent(in) :: alpha(me+nb-1), beta(me+nb-1)
         Real(8), Intent(in) :: evec(lde,me), base(ldb,nb)
         Integer, Intent(out) :: ierr
       End Subroutine trl_write_checkpoint
       ! find the minimum of the flag on all PEs
       Function trl_sync_flag(mpicom, inflag) Result(outflag)
         Implicit None
         Integer :: outflag
         Integer, Intent(in) :: mpicom, inflag
       End Function trl_sync_flag
       Function trl_pe_filename(base, my_rank, npe) Result(filename)
         Implicit None
         Integer, Intent(in) :: my_rank, npe
         Character(*), Intent(in) :: base
         Character(132) :: filename
       End Function trl_pe_filename
    End Interface
    Integer :: ii, c1, c2
    title = trl_pe_filename(info%cpfile, info%my_pe, info%npes)
    Call System_clock(c1)
    Call trl_write_checkpoint(info%cpio, title, nrow, alpha,&
         & beta, evec, lde, j1n, base, ldb, j2n, ii)
    Call System_clock(c2)
    If (c2.Gt.c1) Then
       info%clk_out = info%clk_out + (c2 - c1)
    Else
       info%clk_out = info%clk_out + ((info%clk_max-c1)+c2)
    End If
    info%wrds_out = info%wrds_out + jnd*(nrow+nrow+2) + nrow + 2
    info%stat = trl_sync_flag(info%mpicom, ii)
  End Subroutine write_checkpoint
  ! dump important variables when return due to error
  Subroutine log_error_state
    Interface
       ! print real array for debugging
       Subroutine trl_print_real(info, title, array)
         Use trl_info
         Type(TRL_INFO_T), Intent(in) :: info
         Character(*), Intent(in) :: title
         Real(8), Dimension(:), Intent(in) :: array
       End Subroutine trl_print_real
       ! print integer array for debugging
       Subroutine trl_print_int(info, title, array)
         Use trl_info
         Type(TRL_INFO_T), Intent(in) :: info
         Character(*), Intent(in) :: title
         Integer, Dimension(:), Intent(in) :: array
       End Subroutine trl_print_int
       ! print a time stamp
       Subroutine trl_time_stamp(iou)
         Implicit None
         Integer, Intent(in) :: iou
       End Subroutine trl_time_stamp
       ! a short of the print info function that can be called by
       ! indivadual PE
       Subroutine trl_terse_info(info, iou)
         Use trl_info
         Implicit None
         Type(TRL_INFO_T), Intent(in) :: info
         Integer, Intent(in) :: iou
       End Subroutine trl_terse_info
    End Interface
    Call trl_time_stamp(info%log_io)
    title = 'Dumping the content of the variables on error..'
    iwrk(1) = info%stat
    Call trl_print_int(info, title, iwrk(1:1))
    i2 = info%log_io
    Call trl_terse_info(info, i2)
    Write (i2,*) 'This Lanczos iteration started With ', kept,&
         & ' vectors.'
    Write (i2,*) 'There are ', jnd, '(', j1, ', ', j2,&
         & ') Lanczos vectors currently.'
    If (jnd .Ne. j1+j2) jnd = j1 + j2
    If (jnd.Lt.0 .Or. jnd.Gt.info%klan) jnd = 0
    title = 'Content of eval ..'
    Call trl_print_real(info, title, eval)
    If (jnd .Gt. 0) Then
       Write (title, *) 'Alpha(1:', jnd, ') .. '
       Call trl_print_real(info, title, alpha(1:jnd))
       title(1:5) = ' Beta'
       Call trl_print_real(info, title, beta(1:jnd))
       Write (title, *) 'Alfrot(1:', jnd, ') .. '
       Call trl_print_real(info, title, alfrot(1:jnd))
       title(1:3) = 'bet'
       Call trl_print_real(info, title, betrot(1:jnd))
    End If
    If (j1 .Gt. 0) Then
       title = 'the First row of evec ..'
       Call trl_print_real(info, title, evec(1, 1:j1))
       Write (title, *) 'row ', nrow, ' of evec ..'
       Call trl_print_real(info, title, evec(nrow, 1:j1))
    End If
    If (j2 .Gt. 0) Then
       title = 'the First row of base ..'
       Call trl_print_real(info, title, base(1,1:j2))
       Write (title, *) 'row ', nrow, ' of base ..'
       Call trl_print_real(info, title, base(nrow,1:j2))
    End If
    If (Associated(qb)) Then
       Write (title, *) 'Content of qb (q_', jnd-1, ') ..'
       Call trl_print_real(info, title, qb)
    End If
    If (Associated(qa)) Then
       Write (title, *) 'Content of qa (q_', jnd, ') ..'
       Call trl_print_real(info, title, qa)
    End If
    If (Associated(rr)) Then
       Write (title, *) 'Content of rr (residual == q_', jnd+1, ') ..'
       Call trl_print_real(info, title, rr)
    End If
    If (info%my_pe.Eq.0 .And. info%log_io.Ne.6) Then
       Write (*, *) 'TRLanczos returned with error ', info%stat
       If (info%verbose .Gt. 0) Then
          Write (*, *) 'Contents of most variables are dumped to ',&
               & 'log file ', Trim(info%log_file)
       Else
          Write (*, *) 'Contents of most variables are dumped to ',&
               &'IO unit # ', info%log_io
       End If
    End If
  End Subroutine log_error_state
End Subroutine trlanczos
!!!
! check to make sure the initial guess vector contains valid nonzero
! numbers if not fill with random numbers
! this routine will also read the checkpoint files to restore the
! previous state of the Lancozs iterations
!!!
Subroutine trl_initial_guess(nrow, evec, lde, mev, base, ldb, nbas,&
     & alpha, beta, j1, j2, info, wrk, lwrk)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T) :: info
  Integer, Intent(in) :: nrow, lde, ldb, mev, nbas, lwrk
  Integer, Intent(out) :: j1, j2
  Real(8) :: wrk(lwrk), base(ldb,nbas), evec(lde,mev)
  Real(8) :: alpha(mev+nbas-1), beta(mev+nbas-1)
  Interface
     ! read the checkpointfile
     Subroutine trl_read_checkpoint(cp_io, filename, nrow, evec, lde,&
          & mev, j1, base, ldb, nbas, j2, alpha, beta, ierr)
       Implicit None
       Character(*), Intent(in) :: filename
       Integer, Intent(in) :: cp_io, nrow, lde, mev, ldb, nbas
       Integer, Intent(out) :: j1, j2, ierr
       Real(8), Intent(out) :: alpha(mev+nbas-1), beta(mev+nbas-1)
       Real(8), Intent(out) :: base(ldb,nbas), evec(lde,mev)
     End Subroutine trl_read_checkpoint
     ! the orthogonalization routine -- Classical Gram-Schmidt
     Subroutine trl_cgs(mpicom, myid, nrow, v1, ld1, m1, v2, ld2,&
          & m2, rr, rnrm, alpha, north, nrand, wrk, ierr)
       Integer, Intent(inout) :: north, nrand, ierr
       Integer, Intent(in) :: mpicom, myid, nrow, ld1, ld2, m1, m2
       Real(8), Intent(inout) :: rnrm, alpha
       Real(8), Intent(in) :: v1(ld1,m1), v2(ld2*m2)
       Real(8), Intent(inout) :: rr(nrow), wrk(m1+m2+m1+m2)
     End Subroutine trl_cgs
     ! generate filenames for each PE
     Function trl_pe_filename(base, my_rank, npe) Result(filename)
       Implicit None
       Integer, Intent(in) :: my_rank, npe
       Character(*), Intent(in) :: base
       Character(132) :: filename
     End Function trl_pe_filename
     ! find the minimum of the flag on all PEs
     Function trl_sync_flag(mpicom, inflag) Result(outflag)
       Implicit None
       Integer :: outflag
       Integer, Intent(in) :: mpicom, inflag
     End Function trl_sync_flag
     Subroutine trl_g_sum(mpicom, n, xx, wrk)
       Implicit None
       Integer, Intent(in) :: mpicom, n
       Real(8), Intent(out) :: wrk(n)
       Real(8), Intent(inout) :: xx(n)
     End Subroutine trl_g_sum
     Subroutine trl_check_orth(info, nrow, v1, ld1, j1, v2, ld2, j2, wrk, lwrk)
       Use trl_info
       Implicit None
       Integer, Intent(in) :: nrow, ld1, ld2, j1, j2, lwrk
       Type(TRL_INFO_T), Intent(in) :: info
       Real(8), Intent(in) :: v1(ld1,j1), v2(ld2,j2)
       Real(8) :: wrk(lwrk)
     End Subroutine trl_check_orth
  End Interface
  !
  ! local variable
  !
  Integer :: i, ii, j, nran, north
  Real(8) :: tmp, rnrm
  Character(132) :: file
  Integer, Dimension(:), Allocatable :: rseed
  !
  ! generate random seeds based on current clock ticks
  Call Random_seed(SIZE=ii)    ! find out the size of seed
  Allocate(rseed(ii))
  Call Random_seed(GET=rseed)  ! get the current seeds
  i = Mod(info%my_pe, ii) + 1
  j = Max(Maxval(rseed), info%nloc)
  Call System_clock(rseed(i))  ! get the current clock ticks
  ! just in case the clocks are absolutely synchronized
  If (info%my_pe .Gt. 0) Then
     rseed(i) = rseed(i) - Sign(Int(info%my_pe*Sqrt(Dble(j))), rseed(i))
  End If
  Call Random_seed(PUT=rseed)  ! initialize the random number generator
  Deallocate(rseed)
  !
  j = info%nec+1
  If (info%guess .Gt. 1) Then
     ! retrieve a check-point file
     i = info%cpio
     If (info%oldcpf .Ne. '') Then
        file = trl_pe_filename(info%oldcpf, info%my_pe, info%npes)
     Else
        file = trl_pe_filename(info%cpfile, info%my_pe, info%npes)
     End If
     Call System_clock(ii)
     Call trl_read_checkpoint(info%cpio, file, nrow, evec(1,j), lde,&
          & mev-info%nec, j1, base, ldb, nbas, j2, alpha(j), beta(j),&
          & i)
     info%stat = trl_sync_flag(info%mpicom, i)
     Call System_clock(i)
     If (i.Gt.ii) Then
        info%clk_in = i - ii
     Else
        info%clk_in = (info%clk_max - ii) + i
     End If
     info%wrds_in = (j1+j2)*(nrow+nrow+2) + nrow + 2
     j1 = j1 + info%nec
     If (info%stat .Ne. 0) Return
  Else
     If (info%guess .Le. 0) Then
        ! generate an arbitrary initial starting vector
        ! if (info%guess == 0), use the vector [1, 1, ...]^T
        ! else perturb some random locations of the above vector
        evec(1:nrow,j) = 1.0D0
        nran = Min(1-info%guess, lwrk)
        nran = 2*(nran/2)
        If (nran.Gt.0 .And. nran.Lt.nrow) Then
           Call Random_number(wrk(1:nran))
           Do i = 1, nran-1, 2
              ii = Int(nrow*wrk(i))+1
              evec(ii,j) = evec(ii,j) + wrk(i+1) - 0.5D0
           End Do
           info%flop = info%flop + nran + nran
        Else If (nran .Ge. nrow) Then
           Call random_vector
        End If
     End If
     j1 = info%nec
     j2 = 0
  End If
  tmp = 0.0
  ! make sure the norm of the next vector can be computed
  wrk(1) = Dot_product(evec(1:nrow,j), evec(1:nrow,j))
  Call trl_g_sum(info%mpicom, 1, wrk(1), wrk(2))
  info%flop = info%flop + nrow+nrow
  If (wrk(1).Ge.Tiny(tmp) .And. wrk(1).Le.Huge(tmp)) Then
     ! set rnrm to let trl_CGS normalize evec(1:nrow, j)
     rnrm = Sqrt(wrk(1))
  Else
     Call random_vector
  End If
  !
  ! orthogonalize initial guess against all existing vectors
  !
  i = 0
  tmp = 1.0D0
  nran = info%nrand
  north = info%north
  If (j1 .Lt. mev) Then
     Call trl_cgs(info%mpicom, info%my_pe, nrow, evec, lde, j1,&
          & base, ldb, 0, evec(1,j1+1), rnrm, tmp, i,&
          & info%nrand, wrk, info%stat)
  Else If (j2 .Le. 0) Then
     Call trl_cgs(info%mpicom, info%my_pe, nrow, evec, lde, j1,&
          & evec, lde, 0, base, rnrm, tmp, i, info%nrand, wrk,&
          & info%stat)
  Else
     Call trl_cgs(info%mpicom, info%my_pe, nrow, evec, lde, j1,&
          & base, ldb, j2, base(1,j2+1), rnrm, tmp, i, info%nrand,&
          & wrk, info%stat)
  End If
  info%flop = info%flop + 4*nrow*((info%north-north)*(j1+j2) + &
       & info%nrand-nran) + nrow
  If (info%verbose.Gt.6) Then
     If (j1 .Lt. mev) Then
        i = j1 + 1
        ii = j2
     Else
        i = j1
        ii = j2 + 1
     End If
     Call trl_check_orth(info, nrow, evec, lde, j1, base, ldb, ii, wrk,&
          & lwrk)
  End If
  Return
Contains
  ! an internal subroutine to generate random vector
  Subroutine random_vector
    Call Random_number(evec(1:nrow, j))
    evec(1:nrow, info%nec+1) = evec(1:nrow, info%nec+1) +&
         & evec(1:nrow, info%nec+1) +&
         & Cshift(evec(1:nrow, info%nec+1), 1) +&
         & Cshift(evec(1:nrow, info%nec+1), -1)
    info%nrand = info%nrand + 1
    info%flop = info%flop + 4*nrow
  End Subroutine random_vector
End Subroutine trl_initial_guess
!!!
! full re-orthogonalization for TRLANCZOS
!
!\Algorithm
!
! 1. if (global re-orthogonalization is needed)
!      call trl_cgs
!    else
!      perform extended local re-reorthogonalization
!    endif
! 2. perform normalization
!
! workspace requirement: lwrk >= 2*(m1+m2)
!!!
Subroutine trl_orth(nrow, v1, ld1, m1, v2, ld2, m2, rr, kept, alpha,&
     & beta, wrk, lwrk, info)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T) :: info
  Integer, Intent(in) :: nrow, ld1, ld2, m1, m2, kept, lwrk
  Real(8), Intent(in), Target :: v1(ld1,m1), v2(ld2,m2)
  Real(8) :: rr(nrow), alpha(m1+m2), beta(m1+m2), wrk(lwrk)
  !
  ! local variables
  !
  Real(8), Parameter :: zero=0.0D0, one=1.0D0
  Integer :: i, jnd, jm1, no, nr
  Real(8) :: tmp
  Real(8), Dimension(:), Pointer :: qa, qb
  Interface
     Subroutine trl_cgs(mpicom, myid, nrow, v1, ld1, m1, v2, ld2,&
          & m2, rr, rnrm, alpha, north, nrand, wrk, ierr)
       Integer, Intent(inout) :: north, nrand, ierr
       Integer, Intent(in) :: mpicom, myid, nrow, ld1, ld2, m1, m2
       Real(8), Intent(inout) :: rnrm, alpha
       Real(8), Intent(in) :: v1(ld1,m1), v2(ld2*m2)
       Real(8), Intent(inout) :: rr(nrow), wrk(m1+m2+m1+m2)
     End Subroutine trl_cgs
     Subroutine trl_g_sum(mpicom, n, xx, wrk)
       Implicit None
       Integer, Intent(in) :: mpicom, n
       Real(8), Intent(out) :: wrk(n)
       Real(8), Intent(inout) :: xx(n)
     End Subroutine trl_g_sum
  End Interface
  !
  ! check for workspace size
  !
  jnd = m1 + m2
  jm1 = jnd - 1
  tmp = zero
  If (ld1.Ge.nrow .And. ld2.Ge.nrow .And. lwrk.Ge.Max(4,jnd+jnd)) Then
     info%stat = 0
  Else
     info%stat = -101
     Return
  End If
!!$Print *, 'Entering TRL_ORTH', m1, m2 
!!$Do i = 1, m1
!!$Print *, 'v(:,',i, ') '
!!$Print *, v1(1:nrow,i)
!!$End Do
!!$Do i = 1, m2
!!$Print *, 'v(:,',m1+i,') '
!!$Print *, v2(1:nrow, i)
!!$End Do
!!$Print *, 'rr '
!!$Print *, rr
  !
  ! compute the norm of the vector RR
  !
  wrk(1) = Dot_product(rr, rr)
  Call trl_g_sum(info%mpicom, 1, wrk(1), wrk(2))
  If (.Not.(wrk(1).Ge.zero) .Or. .Not.(wrk(1).Le.Huge(tmp))) Then
     info%stat = -102
     Return
  End If
  beta(jnd) = Sqrt(wrk(1))
  tmp = alpha(jnd)*alpha(jnd)
  If (jm1 .Gt. kept) Then
     tmp = tmp + beta(jm1)*beta(jm1)
     info%flop = info%flop + 2*nrow + 4
  Else If (kept.Gt.0) Then
     tmp = tmp + Dot_product(beta(1:jm1), beta(1:jm1))
     info%flop = info%flop + 2*(nrow + kept + 2)
  End If
  !
  ! decide whether to perform global or local reothogonalization
  If (jnd .Ge. info%ntot) Then
     ! very specfial case: do nothing more
     info%stat = 0
     beta(jnd) = zero
  Else If (tmp.Ge.wrk(1) .Or. wrk(1).Le.Tiny(tmp) .Or. &
       & beta(jm1).Le.Tiny(tmp) .Or. jm1.Eq.kept) Then
     ! perform global re-orthogonalization
     nr = info%nrand
     no = info%north
     Call trl_cgs(info%mpicom, info%my_pe, nrow, v1, ld1, m1, v2, &
          & ld2, m2, rr, beta(jnd), alpha(jnd), info%north, &
          & info%nrand, wrk, info%stat)
     info%flop = info%flop + 4*nrow*((info%north-no)*jnd +&
          & info%nrand-nr) + nrow
  Else If (jnd .Gt. 1) Then
     ! perform local re-orthogonalization against two previous vectors
     If (m2.Gt.1) Then
        qa => v2(1:nrow, m2)
        qb => v2(1:nrow, m2-1)
     Else If (m2.Eq.1) Then
        qa => v2(1:nrow, 1)
        qb => v1(1:nrow, m1)
     Else
        qa => v1(1:nrow, m1)
        qb => v1(1:nrow, jm1)
     End If
     wrk(1) = zero
     wrk(2) = zero
     Do i = 1, nrow
        wrk(1) = wrk(1) + qa(i)*rr(i)
        wrk(2) = wrk(2) + qb(i)*rr(i)
     End Do
     Call trl_g_sum(info%mpicom, 2, wrk(1:2), wrk(3:4))
     alpha(jnd) = alpha(jnd) + wrk(1)
     rr = rr - wrk(1)*qa - wrk(2)*qb
     tmp = one / beta(jnd)
     rr = tmp*rr
     info%flop = info%flop + 9*nrow
  Else
     ! perform local re-orthogonalization against the only vector
     If (m1.Eq.1) Then
        qa => v1(1:nrow, 1)
     Else
        qa => v2(1:nrow, 1)
     End If
     wrk(1) = Dot_product(qa, rr)
     Call trl_g_sum(info%mpicom, 1, wrk(1:1), wrk(2:2))
     alpha(jnd) = alpha(jnd) + wrk(1)
     rr = rr - wrk(1)*qa
     tmp = one / beta(jnd)
     rr = tmp*rr
     info%flop = info%flop + 5*nrow
  End If
  ! when beta(jnd) is exceedingly small, it should be treated as zero
  If (info%stat .Eq. 0) Then
     If (beta(jnd) .Le. Epsilon(tmp)*Abs(alpha(jnd))) Then
        beta(jnd) = zero
     Else If (jnd .Ge. info%ntot) Then
        beta(jnd) = zero
     End If
  End If
End Subroutine trl_orth
!!!
! transforms an real symmetric arrow matrix into a
! symmetric tridiagonal matrix
!!!
Subroutine trl_tridiag(nd, alpha, beta, rot, alfrot, betrot, wrk, lwrk,&
     & ierr)
  Implicit None
  Integer, Intent(in) :: nd, lwrk
  Integer, Intent(out) :: ierr
  Real(8), Intent(in), Dimension(nd) :: alpha, beta
  Real(8) :: rot(nd*nd), alfrot(nd), betrot(nd), wrk(lwrk)
  External dsytrd, dorgtr
  !  Interface
  !     Subroutine dsytrd(uplo, n, a, lda, d, e, tau, wrk, lwrk, ierr)
  !       Character, Intent(in) :: uplo
  !       Integer, Intent(in) :: n, lda, lwrk
  !       Integer :: ierr
  !       Real(8) :: a(lda,n), d(n), e(n), tau(n), wrk(lwrk)
  !     End Subroutine dsytrd
  !     Subroutine dorgtr(uplo, n, a, lda, tau, wrk, lwrk, ierr)
  !       Character, Intent(in) :: uplo
  !       Integer, Intent(in) :: n, lda, lwrk
  !       Integer :: ierr
  !       Real(8) :: a(lda,n), tau(n), wrk(lwrk)
  !     End Subroutine dorgtr
  !  End Interface
  !
  ! local variables
  !
  Character, Parameter :: upper='U'
  Integer :: lwrk2
  !
  ! special case
  !
  If (nd .Le. 1) Then
     rot=1.0D0
     alfrot = alpha
     betrot = beta
     ierr = 0
     Return
  End If
  If (lwrk.Ge.nd+nd) Then
     ierr = 0
  Else
     ierr = -111
     Return
  End If
  !
  ! first form the array matrix as a full matrix in rot
  !
  rot = 0.0D0
  rot(1:nd*nd:nd+1) = alpha(1:nd)
  rot((nd-1)*nd+1:nd*nd-1) = beta(1:nd-1)
  rot(nd:nd*(nd-1):nd) = beta(1:nd-1)
  lwrk2 = lwrk - nd
  !
  ! call LAPACK routines to reduce the matrix into tridiagonal form
  ! and generate the rotation matrix
  !
  Call dsytrd(upper, nd, rot, nd, alfrot, betrot, wrk, wrk(nd+1), &
       &  lwrk2, ierr)
  If (ierr .Ne. 0) Then
     ierr = -112
     Return
  End If
  betrot(nd) = beta(nd)
  Call dorgtr(upper, nd, rot, nd, wrk, wrk(nd+1), lwrk2, ierr)
  If (ierr .Ne. 0) Then
     ierr = -113
     Return
  End If
End Subroutine trl_tridiag
!!!
! evaluates the eigenvalues and their corresponding residual norms of a
! real symmetric tridiagonal matrix
!
! it returns eigenvalues in two sections
! the first section is the locked eigenvalues, their residual norms are
! zero
! the second section contains the new Ritz values, they are in
! ascending order.
! res will contain corresponding residual norms
!!!
Subroutine trl_get_eval(nd, locked, alpha, beta, lambda, res, wrk,&
     & lwrk, ierr)
  Implicit None
  Integer, Intent(in) :: nd, locked, lwrk
  Integer :: ierr
  Real(8), Intent(in) :: alpha(nd), beta(nd)
  Real(8), Intent(out) :: lambda(nd), res(nd), wrk(lwrk)
  External dstqrb
  !  Interface
  !     Subroutine dstqrb(n, d, e, z, wrk, ierr)
  !       Integer, Intent(in) :: n
  !       Integer :: ierr
  !       Real(8) :: d(n), e(n), z(n), wrk(n+n-2)
  !     End Subroutine dstqrb
  !  End Interface
  If (lwrk .Ge. 3*nd) Then
     ierr = 0
  Else
     ierr = -121
     Return
  End If
  lambda = alpha
  wrk(1:nd-locked) = beta(locked+1:nd)
  Call dstqrb(nd-locked, lambda(locked+1:nd), wrk(1:nd-locked),&
       & res(locked+1:nd), wrk(1+nd:), ierr)
  If (ierr .Eq. 0) Then
     res(1:locked) = 0.0D0
     res(locked+1:nd) = beta(nd) * Abs(res(locked+1:nd))
  Else
     ierr = -122
  End If
End Subroutine trl_get_eval
!!!
! count the numer of wanted eigenvalues that have small residual norms
! eigenvalues are assumed to be order from small to large
!!!
Subroutine trl_convergence_test(nd, lambda, res, info, wrk)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T) :: info
  Integer, Intent(in) :: nd
  Real(8), Intent(in) :: lambda(nd), res(nd)
  Real(8), Intent(out) :: wrk(nd+nd)
  Real(8) :: bnd
  Integer :: i, j, ncl, ncr
  External dsort2
  !
  ! copy lambda and res to wrk, sort them in ascending order of lambda
  wrk(nd+1:nd+nd) = lambda(1:nd)
  wrk(1:nd) = Abs(res(1:nd))
  Call dsort2(nd, wrk(nd+1:nd+nd), wrk(1:nd))
  !
  ! determine the convergence rate of the previous target
  If (info%tmv.Gt.0 .And. info%matvec.Gt.info%tmv) Then
     j = 1
     bnd = Abs(lambda(j)-info%trgt)
     Do i = 1, nd
        If (Abs(lambda(i)-info%trgt) .Lt. bnd) Then
           bnd = Abs(lambda(i)-info%trgt)
           j = i
        End If
     End Do
     If (info%tres .Gt. res(j)) Then
        bnd = res(j) / info%tres
        If (bnd .Gt. 0.0D0) Then
           info%crat = Log(bnd)/Dble(info%matvec-info%tmv)
        Else
           info%crat = 0.0D0
        End If
     Else
        info%crat = 1.0D0
     End If
  End If
  !
  ! find out who has converged at the lower end of the spectrum
  info%anrm = Max(info%anrm, Abs(wrk(nd+1)), Abs(wrk(nd+nd)))
  bnd = Tiny(info%anrm) + info%tol*info%anrm
  ncl = 0
  ncr = nd+1
  If (info%lohi .Le. 0) Then
     ncl = nd
     i = 1
     Do While (i.Le.nd)
        If (wrk(i) .Lt. bnd) Then
           i = i + 1
        Else
           ncl = i - 1
           i = nd + 1
        End If
     End Do
  End If
  ! find out who has converged at the high end of the spectrum
  If (info%lohi .Ge. 0) Then
     ncr = 1
     i = nd
     Do While (i.Ge.1)
        If (wrk(i) .Lt. bnd) Then
           i = i - 1
        Else
           ncr = i + 1
           i = 0
        End If
     End Do
  End If
  ! determine the number of wanted eigenvalues that have converged
  ! compute the next target
  info%tmv = info%matvec
  If (info%lohi .Lt. 0) Then
     info%nec = ncl
     info%trgt = wrk(nd+Min(nd,ncl+1))
     info%tres = wrk(Min(nd,ncl+1))
  Else If (info%lohi .Gt. 0) Then
     info%nec = nd-ncr+1
     info%trgt = wrk(nd+Max(1,ncr-1))
     info%tres = wrk(Max(1,ncr-1))
  Else
     If (ncr .Le. ncl) Then
        ncl = nd/2
        ncr = ncl+1
        info%trgt = wrk(nd+(nd+1)/2)
        info%tres = wrk((nd+1)/2)
     Else If (wrk(ncl+1).Le.wrk(ncr-1)) Then
        info%trgt = wrk(nd+ncl+1)
        info%tres = wrk(ncl+1)
     Else
        info%trgt = wrk(nd+ncr-1)
        info%tres = wrk(ncr-1)
     End If
     info%nec = ncl + nd - ncr + 1
     Do i = ncl+1, ncr-1
        If (wrk(i) .Lt. bnd) info%nec = info%nec + 1
     End Do
  End If
End Subroutine trl_convergence_test
!!!
! sort the eigenvalues so that the wanted eigenvalues are ouputed to the
! user in the front of the arrays
! the final Ritz values are in ascending order so that
! DSTEIN can be used to compute the eigenvectors
!!!
Subroutine trl_sort_eig(nd, lohi, nec, lambda, res)
  Implicit None
  Integer, Intent(in) :: nd, lohi, nec
  Real(8) :: lambda(nd), res(nd)
  External dsort2, dsort2a
  Integer :: i, j
  If (lohi .Eq. 0) Then
     ! sort the eigenvalues according to their absolute values
     Call dsort2a(nd, lambda, res)
     ! sort the first nec eigenvalue in the order of lambda
     Call dsort2(nec, lambda, res)
  Else
     ! sort the eigenvalues and residual norms in ascending order of the
     ! eigenvalues
     Call dsort2(nd, lambda, res)
     If (lohi .Gt. 0) Then
        ! move the largest ones to the front (still ascending order)
        j = nd - nec + 1
        Do i = 1, nec
           res(i) = res(j)
           lambda(i) = lambda(j)
           j = j + 1
        End Do
     End If
  End If
End Subroutine trl_sort_eig
!!!
! generating eigenvectors of the projected eigenvalue problem
! corresponding to the given Ritz values
! using LAPACK routine DSTEIN (inverse iterations)
!
! workspace size:
! iwrk: dimension(4*nd)
! wrk:  lwrk >= 5*nd
!!!
Subroutine trl_get_tvec(nd, alpha, beta, irot, nrot, rot, nlam,&
     & lambda, yy, iwrk, wrk, lwrk, ierr)
  Implicit None
  Integer, Intent(in) :: nd, nrot, nlam, irot, lwrk
  Integer :: ierr, iwrk(nd,4)
  Real(8), Intent(in) :: alpha(nd), beta(nd), lambda(nlam), rot(nrot*nrot)
  Real(8) :: wrk(lwrk), yy(nd*nlam)
  !
  ! local variables
  Character, Parameter :: notrans='N'
  Real(8), Parameter :: zero=0.0D0, one=1.0D0
  Integer :: i, j, k, ncol, ii, ioff
  !
  ! conventional external subprograms
  External dstein, dgemm, dgemv
  If (nlam .Le. 0) Return
  If (lwrk.Ge.5*nd) Then
     ierr = 0
  Else
     ierr = -131
     Return
  End If
  !
  ! set up IBLOCK and ISPLIT for calling dstein
  iwrk(1:nd,1) = 1
  iwrk(1:nd,2) = nd
  Call dstein(nd, alpha, beta, nlam, lambda, iwrk(1,1), iwrk(1,2), yy,&
       & nd, wrk, iwrk(1,3), iwrk(1,4), ierr)
  If (ierr .Ne. 0) Then
     Print *, 'TRL_GET_TVEC: dstein failed with error code ', ierr
     ierr = -132
     Return
  End If
  !
  ! apply the rotations to the IROT+1:IROT+NROT rows of YY
  ! generates results 'NCOL' columns at a time
  If (nrot.Gt.1) Then
     ncol = lwrk / nrot
     Do i = 1, nlam, ncol
        j = Min(nlam, i+ncol-1)
        k = j - i + 1
        If (k .Gt. 1) Then
           Call dgemm(notrans, notrans, nrot, k, nrot, one, rot, nrot,&
                & yy((i-1)*nd+irot+1), nd, zero, wrk, nrot)
           Do ii = i-1, j-1
              ioff = (ii-i+1)*nrot
              yy(ii*nd+irot+1:ii*nd+irot+nrot) = wrk(ioff+1:ioff+nrot)
           End Do
        Else
           Call dgemv(notrans, nrot, nrot, one, rot, nrot,&
                & yy((i-1)*nd+irot+1), 1, zero, wrk, 1)
           yy((i-1)*nd+irot+1:(i-1)*nd+irot+nrot) = wrk(1:nrot)
        End If
     End Do
  End If
End Subroutine trl_get_tvec
!!!
! compute all eigenvalues and eigenvectors of the projected matrix
! use LAPACK routine DSYEV
! The eigenvectors corresponding to lambda(1:nlam) are placed at the
! first nlam*nd locations of yy on exit.
!
! On entry, Alpha and Beta defines the arrayhead plus tridiagonal
! matrix.
! On return, Alpha will store the Ritz values in the same order as
! lambda.
!
! lwrk >= 3*nd
!!!
Subroutine trl_get_tvec_a(nd, kept, alpha, beta, nlam, lambda,&
     & yy, wrk, lwrk, iwrk, ierr)
  Implicit None
  Integer, Intent(in) :: kept, nd, nlam, lwrk
  Integer :: ierr, iwrk(nd)
  Real(8), Intent(in) :: lambda(nd)
  Real(8) :: alpha(nd), beta(nd), yy(nd*nd), wrk(lwrk)
  !
  ! local variables
  Integer :: i, j, i2, j2, ii
  Real(8) :: tmp
  External dsyev
  !
  ! fill yy with the projection matrix, then call DSYEV
  If (nlam .Le. 0) Return
  If (lwrk.Ge.nd+nd+nd) Then
     ierr = 0
  Else
     ierr = -141
     Return
  End If
  yy = 0.0D0
  yy(1:nd*nd:nd+1) = alpha(1:nd)
  If (kept.Gt.0) yy(kept*nd+1:kept*nd+kept) = beta(1:kept)
  yy((kept+1)*(nd+1):nd*nd:nd+1) = beta(kept+1:nd-1)
  Call dsyev('V', 'U', nd, yy, nd, alpha, wrk, lwrk, ierr)
  If (ierr .Ne. 0) Then
     ierr = -142
     Return
  End If
  If (nlam .Ge. nd) Return
  !
  ! reorder the eigenvectors
  ! both lambda(1:kept) and alpha are in ascending order
  !
  tmp = Max(alpha(nd)-alpha(1), Abs(alpha(nd)), Abs(alpha(1)))
  tmp = Epsilon(tmp)*tmp*nd
  j = 1
  i = 1
  Do While (i .Le. nlam)
     ! move j so that alpha(j) is within tmp distance away
     ii = j
     j = nd
     Do While (ii .Le. nd)
        If (alpha(ii) .Lt. lambda(i)-tmp) Then
           ii = ii + 1
        Else
           j = ii
           ii = nd+1
        End If
     End Do
     If (alpha(j) .Gt. lambda(i)+tmp) Then
        ierr = -143
        Return
     End If
     ! identify the group size in lambda
     ii = i + 1
     i2 = nlam
     Do While (ii .Le. nlam)
        If (lambda(ii) .Le. lambda(i)+tmp) Then
           ii = ii + 1
        Else
           i2 = ii - 1
           ii = nd+1
        End If
     End Do
     ! identify the group size in alpha
     ii = j + 1
     j2 = nd
     Do While (ii .Le. nd)
        If (alpha(ii) .Le. lambda(i)+tmp) Then
           ii = ii + 1
        Else
           j2 = ii - 1
           ii = nd+1
        End If
     End Do
     ! assign the index values
     If (j2.Eq.j .And. i2.Eq.i) Then
        iwrk(i) = j
     Else If (j2-j .Eq. i2-i) Then
        iwrk(i:i2) = (/(ii, ii=j, j2)/)
     Else If (j2-j .Gt. i2-i) Then
        j2 = j + i2 - i
        iwrk(i:i2) = (/(ii, ii=j, j2)/)
     Else If (j2 .Lt. nd) Then
        i2 = i + j2 - j
        iwrk(i:i2) = (/(ii, ii=j, j2)/)
     Else
        ierr = -144
        Return
     End If
     i = i2 + 1
     j = j2 + 1
  End Do
  ! perform the actual copying
  Do i = 1, nlam
     j = iwrk(i)
     If (j.Gt.i) Then
        alpha(i) = alpha(j)
        yy((i-1)*nd+1:i*nd) = yy((j-1)*nd+1:j*nd)
     End If
  End Do
End Subroutine trl_get_tvec_a
!!!
! move the Ritz pairs with extremely small residual norms to the front
! of the arrays so that locking can be performed cleanly
!!!
Subroutine trl_set_locking(jnd, nlam, lambda, res, yy, anrm, locked)
  Implicit None
  Integer, Intent(in) :: jnd, nlam
  Integer :: locked
  Real(8), Intent(in) :: anrm
  Real(8) :: lambda(nlam), res(nlam), yy(jnd*nlam)
  !
  Real(8), Parameter :: zero = 0.0D0
  Integer :: i, j, ii, ioff
  Real(8) :: tmp, eps, small
  Logical :: ti, tj
  small(tmp, eps) = eps*Max(Abs(tmp), eps*anrm)
  eps = Epsilon(zero)
  i = 1
  j = nlam
  ti = (Abs(res(i)) .Le. small(lambda(i),eps))
  tj = (Abs(res(j)) .Le. small(lambda(j),eps))
  Do While (i .Lt. j)
     If (ti) Then
        res(i) = zero
        i = i + 1
        If (i.Le.j) Then
           ti = (Abs(res(i)) .Le. small(lambda(i),eps))
        Else
           ti = .False.
        End If
     Else
        If (tj) Then
           tmp = lambda(i)
           lambda(i) = lambda(j)
           lambda(j) = tmp
           res(j) = res(i)
           res(i) = zero
           ioff = (j-i)*jnd
           Do ii = i*jnd-jnd+1, i*jnd
              tmp = yy(ii)
              yy(ii) = yy(ii+ioff)
              yy(ii+ioff) = tmp
           End Do
           i = i + 1
           If (i.Le.j) Then
              ti = (Abs(res(i)) .Le. small(lambda(i),eps))
           Else
              ti = .False.
           End If
        End If
        j = j - 1
        If (j.Gt.i) Then
           tj = (Abs(res(j)) .Le. small(lambda(j),eps))
        Else
           tj = .False.
        End If
     End If
  End Do
  If (ti) Then
     locked = i
  Else
     locked = i - 1
  End If
End Subroutine trl_set_locking
!!!
! compute the Ritz vectors from the basis vectors and the eigenvectors
! of the projected system
! the basis vectors may be stored in two separete arrays
! the result need to be stored back in them
!
! lwrk should be no less than ny (lwrk>=ny) ***NOT checked inside***
!!!
Subroutine trl_ritz_vectors(nrow, lck, ny, yy, ldy, vec1, ld1, m1,&
     & vec2, ld2, m2, wrk, lwrk)
  Implicit None
  Integer, Intent(in) :: nrow, ny, lck, ldy, ld1, m1, m2, ld2, lwrk
  Real(8), Intent(in) :: yy(ldy, ny)
  Real(8) :: vec1(ld1,m1), vec2(ld2,m2), wrk(lwrk)
  !
  ! local variables
  !
  Character, Parameter :: notrans='N'
  Real(8), Parameter :: zero=0.0D0, one=1.0D0
  Integer :: i,j,k,stride,ii,jl1,jl2,il1,il2,kv1
  !
  ! conventional external procedures
  External dgemm, dgemv
  !
  If (lck .Le. m1) Then
     il1 = lck + 1
     jl1 = m1 - lck
     il2 = 1
     jl2 = m2
  Else
     il1 = m1+1
     jl1 = 0
     il2 = lck - m1 + 1
     jl2 = m1 + m2 - lck
  End If
  If (jl1.Eq.0 .And. jl2.Eq.0) Return
  kv1 = Min(m1-il1+1, ny)
  wrk = zero
  If (ny .Gt. 1) Then
     stride = lwrk / ny
     Do i = 1, nrow, stride
        j = Min(nrow, i+stride-1)
        k = j - i + 1
        If (jl1 .Gt. 0) Then
           Call dgemm(notrans, notrans, k, ny, jl1, one, vec1(i,il1),&
                & ld1, yy, ldy, zero, wrk, k)
        Else
           wrk = zero
        End If
        If (jl2 .Gt. 0) Call dgemm(notrans, notrans, k, ny, jl2, one,&
             & vec2(i,il2), ld2, yy(jl1+1,1), ldy, one, wrk, k)
        Do ii = 0, kv1-1
           vec1(i:j,ii+il1) = wrk(ii*k+1:ii*k+k)
        End Do
        Do ii = 0, ny-kv1-1
           vec2(i:j,ii+il2) = wrk((kv1+ii)*k+1:(kv1+ii)*k+k)
        End Do
     End Do
  Else If (ny.Eq.1) Then
     stride = lwrk
     Do i = 1, nrow, stride
        j = Min(nrow, i+stride-1)
        k = j - i + 1
        If (jl1 .Gt. 0) Then
           Call dgemv(notrans, k, jl1, one, vec1(i,il1), ld1, yy, 1,&
                & zero, wrk, 1)
           If (jl2 .Gt. 0) Call dgemv(notrans, k, jl2, one,&
                & vec2(i,il2), ld2, yy(jl1+1,1), 1, one, wrk, 1)
        Else
           Call dgemv(notrans, k, jl2, one, vec2(i,il2),&
                & ld2, yy(jl1+1,1), 1, zero, wrk, 1)
        End If
        If (kv1 .Gt. 0) Then
           vec1(i:j,il1) = wrk(1:k)
        Else
           vec2(i:j,il2) = wrk(1:k)
        End If
     End Do
  End If
End Subroutine trl_ritz_vectors
!!!
! full Gram-Schmidt routine -- orthogonalize a new vector against
! all existing vectors
!
! On entry:
! rnrm and alpha are expected to have the correct values.
!!!
Subroutine trl_cgs(mpicom, myid, nrow, v1, ld1, m1, v2, ld2, m2, rr,&
     & rnrm, alpha, north, nrand, wrk, ierr)
  Implicit None
  Integer, Intent(in) :: mpicom, myid, nrow, ld1, m1, ld2, m2
  Integer, Intent(inout) :: north, nrand, ierr
  Real(8), Intent(inout) :: rnrm, alpha
  Real(8), Intent(in) :: v1(ld1,m1), v2(ld2,m2)
  Real(8), Intent(inout) :: rr(nrow), wrk(m1+m2+m1+m2)
  Interface
     ! global dot-product routine, calls dgemv
     Subroutine trl_g_dot(mpicom, nrow, v1, ld1, m1, v2, ld2, m2, rr, wrk)
       Implicit None
       Integer, Intent(in) :: mpicom, nrow, ld1, ld2, m1, m2
       Real(8), Intent(in) :: v1(ld1,m1), v2(ld2,m2), rr(nrow)
       Real(8), Intent(out) :: wrk(m1+m2+m1+m2)
     End Subroutine trl_g_dot
     Subroutine trl_g_sum(mpicom, n, xx, wrk)
       Implicit None
       Integer, Intent(in) :: mpicom, n
       Real(8), Intent(out) :: wrk(n)
       Real(8), Intent(inout) :: xx(n)
     End Subroutine trl_g_sum
  End Interface
  !
  ! local variables
  !
  Character, Parameter :: notrans='N'
  Integer, Parameter :: maxorth = 5
  Real(8), Parameter :: one=1.0D0, zero=0.0D0
  Integer :: i, j, nold, irnd, cnt
  Real(8) :: tmp
  Logical :: sane
  External dgemv
  !
  nold = m1+m2
  If (ld1.Ge.nrow .And. (ld2.Ge.nrow.Or.m2.Le.0)) Then
     ierr = 0
  Else
     ierr = -201
     Return
  End If
  irnd = 0
  If (nold .Gt. 0) Then
     ! repeated performing Gram-Schmidt procedure to ensure
     ! orthogonality
     cnt = 0
     Do While (cnt .Le. maxorth)
        Call trl_g_dot(mpicom, nrow, v1, ld1, m1, v2, ld2, m2, rr, wrk)
!!$print *, 'TRL_CGS', cnt, wrk(1:nold)
        If (m1 .Gt. 1) Then
           Call dgemv(notrans, nrow, m1, -one, v1, ld1, wrk, 1,&
                & one, rr, 1)
        Else If (m1 .Eq. 1) Then
           rr = rr - wrk(1)*v1(1:nrow,1)
        End If
        If (m2 .Gt. 1) Then
           Call dgemv(notrans, nrow, m2, -one, v2, ld2, wrk(m1+1),&
                & 1, one, rr, 1)
        Else If (m2 .Eq. 1) Then
           rr = rr - wrk(nold)*v2(1:nrow,1)
        End If
        sane = .True.
        If (irnd .Eq. 0) Then
           tmp = (ld1*nold)*Epsilon(one)*Max(Abs(alpha), rnrm)
           If (Abs(wrk(nold)).Gt.tmp .And. tmp.Gt.zero) Then
              Print *, 'TRL_CGS(', nold, ') caution: small eigenvalues ', &
                   & 'may not have any correct digits'
              sane = .False.
           End If
           alpha = alpha + wrk(nold)
        End If
        north = north + 1
        cnt = cnt + 1
        tmp = Dot_product(wrk(1:nold), wrk(1:nold))
        wrk(1) = Dot_product(rr,rr)
        Call trl_g_sum(mpicom, 1, wrk, wrk(2))
        rnrm = Sqrt(wrk(1))
        !
        ! decisions about whether to re-orthogonalize is based on
        ! relative size of tmp and wrk(1) (R norm square)
        If (wrk(1) .Gt. 1e4*tmp) Then
           ! no need for more orthogonalization
           cnt = maxorth + 1
        Else If (((.Not.(wrk(1) .Gt. Epsilon(tmp)*tmp) .And. cnt.Gt.1) .Or.&
             & .Not.(sane)) .And. irnd.Lt.maxorth) Then
           ! the input vector is so nearly linear dependent on the
           ! existing vectors, we need to perturb it in order to
           ! generate a new vector that is orthogonal to the existing
           ! ones
           ! the perturbation is done in two stages:
           ! -- perturbing only one number
           ! -- call random_number to generate a whole random vector
           cnt = 0
           irnd = irnd + 1
           nrand = nrand + 1
           If (irnd.Eq.1) Then
              ! modify only one element
              Call Random_number(tmp)
              i = Int(nrow*tmp) + 1
              If (rnrm .Lt. Epsilon(tmp) .And. rnrm.Gt.Tiny(tmp)) Then
                 tmp = one / rnrm
                 rr = tmp*rr
                 rnrm = one
              Else If (rnrm .Le. Tiny(tmp)) Then
                 rr = zero
                 rnrm = one
              End If
              tmp = 0.5D0
              Do While (Abs(tmp-0.5D0) .Le. Epsilon(tmp))
                 Call Random_number(tmp)
                 rr(i) = rr(i) + rnrm * (tmp - 0.5D0)
              End Do
           Else
              ! fill with random numbers produced by intrinsic function
              Do i = 0, myid
                 Call Random_number(tmp)
              End Do
              i = Int(nrow*tmp) + 1
              Call Random_number(tmp)
              j = Int(nrow*tmp) + 1
              If (i .Lt. j) Then
                 Call Random_number(rr(i:j))
              Else If (i .Gt. j) Then
                 Call Random_number(rr(j:i))
              Else
                 Call Random_number(rr)
              End If
           End If
           rr = rr + rr + Cshift(rr, 1) + Cshift(rr, -1)
        End If
     End Do
     ! failed to reduce the dot-products between the new vector
     ! and the old vectors to small number -- usually an indication of
     ! problem with orthogonality of the old vectors
     If (.Not. (wrk(1).Ge.tmp)) Then
        ierr = -203
     End If
  End If
  !
  ! normalization
  !
  If (ierr.Eq.0) Then
     If (rnrm.Gt.Tiny(rnrm)) Then
        tmp = one / rnrm
        rr = tmp * rr
     Else
        ierr = -204
     End If
  End If
  If (irnd.Gt.0) rnrm = zero
End Subroutine trl_cgs
