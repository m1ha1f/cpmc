! $Id: trlan.f90,v 1.1 2004/05/04 06:29:35 fowlkes Exp $
!!! Top (user) level routines
!\Description:
! A thick-restart Lanczos routine for computing eigenvalues and
! eigenvectors of a real symmetric operator/matrix (A).
! -- only accept one input vector, the input vector is expected
!    to be stored in the (nec+1)st column of EVEC.
! -- it extends the Lanczos basis one vector at a time.
! -- orthogonality among the Lanczos vectors are maintained using
!    full re-orthogonalization when necessary.
!
!\Requirements:
! 1) User supplies OP with the specified interface
! 2) If (info%nec>0), evec(1:nrow, 1:info%nec) must contain valid
! eigenvectors and eval(1:nec) must be the corresponding eigenvalues.
! These eigenpairs are assumed to have zero residual norm inside TRLAN.
! 3) lde >= nrow
! 4) The arrays evec and eval must be large enough to store the
! solutions, i.e., mev >= info%ned and mev >= info%nec
! 5) The array wrk may be of arbitrary size.  Internally, the workspace
! size is
!
!        nrow*max(0,info%ned-size(evec,2))+maxlan*(maxlan+10)
!
! If wrk is smaller than this, trlan routine will allocate additional
! workspace to accommodate.
!
!\Arguments:
! op   -- the operator. Given a set of vectors X, it returns
!         op(X) == A*X
! info -- data structure to store the information about the eigenvalue
!         problem and the progress of TRLAN
! nrow -- number of rows that is on this processor (local problem size)
! mev  -- the number of Ritz pairs can be stored in eval and evec.
! eval -- array to store the eigenvalues
! evec -- array to store the eigenvectors
! lde  -- the leading dimension of the array evec (lde >= nrow)
! wrk  -- (optional) workspace, if it is provided and there is enough
!         space, the residual norm of the converged eigenpairs will be
!         stored at wrk(1:info%nec) on exit.
! lwrk -- (optional) the size of WRK.  When both WRK and LWRK are
!         present, than LWRK should correctly indicate the size of WRK.
!         If WRK is present by not LWRK, the size of WRK is assumed to
!         be MEV which is only enough to store the residual norms on
!         exit.  If WRK is not present, LWRK is not used even if it is
!         present.
!
! NOTE ON USING EXPLICIT SIZE ARRARY ARGUMENTS !
! This interface can be shorten to use assumed shape array arguments. 
! However, because it may cause extra copying of arrays when calling
! BLAS and LAPACK routines, we have decided to use only explicit size
! arrary arguments for most of the computing procedures.  Only the
! printing routines and debugging routines use assumed shape array
! arguments.
!
! the operator that defines the eigenvalue problem is expected to have
! the following interface
!  Subroutine OP(nrow, ncol, xin, ldx, yout, ldy)
!    Integer, Intent(in) :: nrow, ncol, ldx, ldy
!    Real(8), Dimension(ldx*ncol), Intent(in) :: xin
!    Real(8), Dimension(ldy*ncol), Intent(out) :: yout
!  End Subroutine OP
! Where the arguments are
! nrow -- number of rows in Xin and Yout on this PE (local problem
!         size)
! ncol -- number of columns in Xin and Yout, i.e., the number of
!         vectors to be multiplied
! xin  -- input vector to be multiplied.  It consists of Ncol column
!         vectors with each column stored in consecutive order
! ldx  -- the leading dimension of the array Xin.  The i-th column
!         vector starts with element (i-1)*ldx+1 and ends with element
!         (i-1)*ldx+nrow in Xin.
! yout -- the result array.  It stores the result of matrix-vector
!         multiplications.
! ldy  -- the leading dimension of the array yout.  The i-th column
!         vector starts with element (i-1)*ldy+1 in Yout and ends with
!         element (i-1)*ldy+nrow.
!
! Since it does not have an interface definition in Fortran 90, it is
! assumed to have an old Fortran 77 interface.  The user may declear
! the arrays Xin and Yout as two-dimension Fortran arrays in their
! implementation of the operator.
!!!
Subroutine trlan(op, info, nrow, mev, eval, evec, lde, wrk, lwrk)
  Use trl_info
  Implicit None
  External op
  Type(TRL_INFO_T) :: info
  Integer, Intent(in) :: nrow, mev, lde
  Real(8) :: eval(mev), evec(lde,mev)
  Integer, Intent(in), Optional :: lwrk
  Real(8), Target, Optional :: wrk(*)
  Interface
     !  ! the operator that defined the eigenvalue problem
     !  Subroutine OP(nrow, ncol, xin, ldx, yout, ldy)
     !    Integer, Intent(in) :: nrow, ncol, ldx, ldy
     !    Real(8), Dimension(ldx*ncol), Intent(in) :: xin
     !    Real(8), Dimension(ldy*ncol), Intent(out) :: yout
     !  End Subroutine OP
     ! the actual workhorse of the restarted Lanczos routine
     Subroutine trlanczos(op, info, nrow, mev, eval, evec, lde, base,&
          & ldb, nbas, wrk, lwrk)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T) :: info
       Integer, Intent(in) :: nrow, mev, lde, ldb, nbas, lwrk
       Real(8) :: eval(1:mev)
       Real(8), Target :: evec(1:lde,1:mev), base(ldb*nbas), wrk(1:lwrk)
       External op
     End Subroutine trlanczos
     Subroutine trl_clear_counter(info, nrow, mev, lde)
       Use trl_info
       Implicit None 
       Type(TRL_INFO_T) :: info
       Integer, Intent(in) :: nrow, mev, lde
     End Subroutine trl_clear_counter
     Subroutine trl_print_setup(info, lbas, lmis, lwrk)
       Use trl_info
       Implicit None 
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: lbas, lmis
       Integer, Intent(in), Optional :: lwrk
     End Subroutine trl_print_setup
     ! find the minimum flag value of all PEs
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
     Subroutine trl_time_stamp(iou)
       Implicit None
       Integer, Intent(in) :: iou
     End Subroutine trl_time_stamp
  End Interface
  Character(132) :: filename
  !
  ! The main task of this subroutine is prepare workspace required and
  ! check arguments
  !
  Integer :: ii, nbas, nmis, ibas, imis, ldb, lwrk0
  !
  ! local arrays -- all declared as 1-D array though they may
  ! accutually be used as 2-D arrays.
  !
  Real(8), Dimension(:), Pointer :: base, misc
  imis = -1    ! if this routine allocated misc, imis will be 0
  ibas = -1    ! if this routine allocated base, ibas will be 0
  Nullify(base, misc)
  Call System_clock(COUNT=ii, COUNT_RATE=info%clk_rate,&
       & COUNT_MAX=info%clk_max)
  info%clk_tot = ii
  !
  If (info%ned > mev) Then
     Print *, 'info%ned (', info%ned, ') is larger than mev (', mev,&
          & ')', 'reducing info%ned'
     info%ned = mev
  End If
  !
  ! there is nothing to do if there is no more eigenvalue to compute
  If (info%ned .Le. info%nec .Or. info%ned .Le. 0) Goto 888
  !
  ! determine the workspace size and leading dimensions
  If (Present(wrk) .And. Present(lwrk)) Then
     lwrk0 = lwrk
  Else If (Present(wrk)) Then
     lwrk0 = mev
  Else
     lwrk0 = 0
  End If
  info%stat = 0
  ldb = ((nrow+3)/4)*4
  If (Mod(ldb,4096) .Eq. 0) ldb = ldb + 4
  Call trl_clear_counter(info, nrow, mev, lde)
  If (info%stat .Ne. 0) Goto 888
  !
  ! Internally, the workspace is broken into two parts
  ! one to store (maxlan+1) Lanczos vectors, and the other to
  ! store all others (size maxlan*(maxlan+ned+14))
  ! The next If-block decides how the two arrays are mapped.
  !
  nbas = Max(1, info%maxlan-mev+1)
  ii = nbas * ldb
  nmis = info%maxlan*(info%maxlan+10)
  If (lwrk0 .Gt. Min(ii, nmis)) Then
     ! use wrk either as base or misc or both depending its size
     If (lwrk0 .Ge. ii+nmis) Then
        ! WRK is large enough for both arrays
        base => wrk(1:ii)
        misc => wrk(ii+1:lwrk0)
        nmis = lwrk0 - ii
     Else If (lwrk0 .Ge. Max(ii, nmis)) Then
        ! associate the larger one of base and misc to WRK
        If (ii .Ge. nmis) Then
           base => wrk(1:ii)
           Allocate(misc(nmis), stat=imis)
           If (imis .Ne. 0) info%stat = -4
        Else
           misc => wrk(1:lwrk0)
           nmis = lwrk0
           Allocate(base(ii), stat=ibas)
           If (ibas .Ne. 0) info%stat = -5
        End If
     Else If (ii .Le. nmis) Then
        ! base is smaller, associate base with WRK
        base => wrk(1:ii)
        Allocate(misc(nmis), stat=imis)
        If (imis .Ne. 0) info%stat = -4
     Else
        ! misc is smaller, associate misc with WRK
        misc => wrk(1:lwrk0)
        nmis = lwrk0
        Allocate(base(ii), stat=ibas)
        If (ibas .Ne. 0) info%stat = -5
     End If
  Else
     ! have to allocate both base and misc
     Allocate(base(ii), stat=ibas)
     If (ibas .Ne. 0) info%stat = -5
     Allocate(misc(nmis), stat=imis)
     If (imis .Ne. 0) info%stat = -4
  End If
  !
  ! make sure every process is successful so far
  ii = trl_sync_flag(info%mpicom, info%stat)
  info%stat = ii
  If (ii .Ne. 0) Goto 888
  !
  ! open debug files if desired
  If (info%verbose.Gt.0 .And. info%log_io.Ne.6) Then
     filename = trl_pe_filename(info%log_file, info%my_pe, info%npes)
     Open(info%log_io, file=filename, status='replace',&
          & action='write', position='rewind', iostat=ii)
     If (ii .Ne. 0) info%log_io = 6
  End If
  If (info%verbose .Gt. 0) Then
     Call trl_time_stamp(info%log_io)
     Call trl_print_setup(info, nbas*ldb, nmis, lwrk0)
  End If
  !
  ! call trlanczos to do the real work
  base = 0.0D0
  misc = 0.0D0
  Call trlanczos(op, info, nrow, mev, eval, evec, lde, base, ldb, nbas,&
       & misc, nmis)
  If (info%verbose.Gt.0 .And. info%log_io.Ne.6) Close(info%log_io)
  ii = Max(info%nec, info%ned)
  If (lwrk0 .Ge. ii) Then
     wrk(1:ii) = misc(1:ii)
  End If
  !
  ! DONE, reclaim the space allocated
888 If (imis .Eq. 0) Deallocate(misc)
  If (ibas .Eq. 0) Deallocate(base)
  Call System_clock(count=ii)
  If (ii .Ge. info%clk_tot) Then
     info%clk_tot = ii - info%clk_tot
  Else
     info%clk_tot = ii + (info%clk_max - info%clk_tot)
  End If
  Return
End Subroutine trlan
!!!
! clear the counters inside info and performs a minimal check on the
! input parameters
!!!
Subroutine trl_clear_counter(info, nrow, mev, lde)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T) :: info
  Integer, Intent(in) :: nrow, mev, lde
  Integer, External :: trl_sync_flag
  Integer :: ntmp
  !!
  info%stat = 0
  If (nrow.Ne.info%nloc .Or. nrow.Gt.info%ntot) Then
     Print *, 'TRLAN: ''info'' not setup for this problem.'
     Print *, '       Please reset ''info'' before calling TRLAN.'
     info%stat = -1
  End If
  If (info%ned+info%ned >= info%ntot) Then
     Print *, 'TRLAN: info%ned (', info%ned, ') is large relative ', &
          & 'to the matrix dimension (', info%ntot, ')'
     Print *, '       it might be more appropriate to use LAPACK ', &
          & 'dsyev/ssyev'
  End If
  If (info%nec .Lt. 0) info%nec = 0
  If (lde .Lt. nrow) Then
     Print *, 'TRLAN: leading dimension of EVEC to small.'
     info%stat = -2
  End If
  If (.Not. (info%tol .Lt. 1D0)) Then
     info%tol = Sqrt(Epsilon(info%tol))
  Else If (.Not. (info%tol .Gt. Tiny(info%tol))) Then
     info%tol = Epsilon(info%tol)
  End If
  info%ned = Min(info%ntot, info%ned)
  If (mev .Lt. info%ned) Then
     Print *, 'TRLAN: array EVAL and EVEC can not hold wanted', &
          & ' number of eigenpairs.'
     info%stat = -3
  End If
  info%nrand = info%stat
  info%stat = trl_sync_flag(info%mpicom, info%nrand)
  !
  ! decide what is a good maximum basis size to use
  !
  If (info%maxlan .Le. Min(info%ned+3, info%ntot)) Then
     info%maxlan = info%ned+Min(info%ned,20)+Int(Log(Dble(info%ntot)))
     info%maxlan = Min(info%maxlan, info%ntot)
  End If
  If (info%maxlan < mev) Then
     ntmp = Min(info%ntot, Max(100+info%ned, 10*info%ned))
     If (mev < ntmp) Then
        info%maxlan = mev
     Else
        info%maxlan = ntmp;
     End If
  End If
  !!! clear regular counters
  info%tmv = -1
  info%klan = Min(info%maxlan, info%ntot)
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
  info%tick_t = 0.0D0
  info%clk_op = 0
  info%tick_o = 0.0D0
  info%clk_orth = 0
  info%tick_h = 0.0D0
  info%clk_res = 0
  info%tick_r = 0.0D0
  info%clk_in = 0
  info%clk_out = 0
  info%wrds_in = 0
  info%wrds_out = 0
  Return
End Subroutine trl_clear_counter
!!!
! print the definition of the eigenvalue proble
!!!
Subroutine trl_print_setup(info, lbas, lmis, lwrk)
  Use trl_info
  Implicit None 
  Type(TRL_INFO_T), Intent(in) :: info
  Integer, Intent(in) :: lbas, lmis
  Integer, Intent(in), Optional :: lwrk
  !
  ! print the problem parameters
  If (info%lohi .Gt. 0) Then
     Write (info%log_io, FMT=100) info%ned, 'largest'
  Else If (info%lohi .Lt. 0) Then
     Write (info%log_io, FMT=100) info%ned, 'smallest'
  Else
     Write (info%log_io, FMT=100) info%ned, 'first converged'
  Endif
  Write(info%log_io, FMT=150) info%nloc, info%my_pe, info%ntot
  Write(info%log_io, FMT=200) 'Maximum basis size:', info%maxlan
  Write(info%log_io, FMT=200) 'Dynamic restarting scheme:', info%restart
  Write(info%log_io, FMT=200) 'Maximum applications of the operator:',&
       & info%maxmv
  Write(info%log_io, FMT=300) 'Relative convergence tolerance:', info%tol
100 Format('TRLAN is to compute ', I6, 1X, A, ' eigenpair(s).')
150 Format('Problem dimension: ', I9, '(PE:', I4, '),', I12, '(Global)')
200 Format(A, T40, I10)
300 Format(A, T40, 1PE10.3)
  ! initial guess
  If (info%guess .Eq. 1) Then
     Write(info%log_io, *) 'User provided the starting vector.'
  Else If (info%guess .Eq. 0) Then
     Write(info%log_io, *) 'TRLAN uses [1,1,...] as starting vctor.'
  Else If (info%guess .Lt. 0) Then
     Write(info%log_io, *) 'TRLAN generates a random starting vector.'
  Else If (info%oldcpf .Ne. '') Then
     Write(info%log_io, *)&
          & 'Restarting with existing checkpoint files ',&
          & Trim(info%oldcpf), '####'
  Else
     Write(info%log_io, *)&
          & 'Restarting with existing checkpoint files ',&
          & Trim(info%cpfile), '####'
  End If
  If (info%cpflag .Gt. 0) Then
     Write(info%log_io, *) 'TLRAN will write about ', info%cpflag,&
          &' sets of checkpointing files ', Trim(info%cpfile), '####.'
  End If
  !
  ! print the workspace size parameters
  Write(info%log_io, *) '(required) array BASE size is ', lbas
  Write(info%log_io, *) '(required) array MISC size is ', lmis
  If (Present(lwrk)) Then
     If (lwrk .Gt. 0) Then
        Write(info%log_io, *) 'Caller has supplied a work array with ',&
          & lwrk, ' elements.'
     Else
        Write(info%log_io, *) 'Caller did not supply work array.'
     End If
  Else
     Write(info%log_io, *) 'Caller did not supply work array.'
  End If
End Subroutine trl_print_setup
!!!
! set information related to debugging, the initialization routine
! trl_init_info sets the parameters to not allow no debug information
! to be printed.
!
! arguments:
! info   -- the info to be modified
! msglvl -- the new message level
!           < 0 : nothing printed
!           1:10 -- the high the level, the more debug information is
!                printed
! iou    -- Fortran I/O unit to be used to write out the debug
!           information
! file   -- leading part of the debug file name.
!           Each processor will generate its own debug file with the
!           file name formed by appending the processor number to
!           the string FILE
! The actual name is composed by the routine TRL_PE_FILENAME which is
! in trlaux.f90
!!!
Subroutine trl_set_debug(info, msglvl, iou, file)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T) :: info
  Integer, Intent(in) :: msglvl, iou
  Character(*), Optional :: file
  info%verbose = msglvl
  info%log_io = iou
  If (Present(file)) Then
     info%log_file = file
  End If
End Subroutine trl_set_debug
!!!
! set up the information related to check-pointing
!!!
Subroutine trl_set_checkpoint(info, cpflag, cpio, file)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T) :: info
  Integer, Intent(in) :: cpflag, cpio
  Character(*), Optional :: file
  !
  ! copy the flags
  info%cpflag = cpflag
  info%cpio = cpio
  If (Present(file)) Then
     info%cpfile = file
  End If
End Subroutine trl_set_checkpoint
!!!
! set up parameters related to initial guesses of the Lanczos iterations
!
! It sets the number of eigenvector already converged (initially
! assumed to be zero) and whether the user has provided initial guess
! vector (initially assumed no).  It can also tell TRLan to read
! checkpoint file that is different from the default name.
!
! It uses info%cpio to open the checkpoint files.  If the default value
! of info%cpio is in use, make sure that the routine is called after
! trl_set_checkpoint function is called.
!
! If oldcpf is not specified, TRLAN will use the string info%cpfile as
! the name of the checkpoint file to use.  (INFO%cpfile is the name of
! the output checkpoint files.)
!!!
Subroutine trl_set_iguess(info, nec, iguess, oldcpf)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T) :: info
  Integer, Intent(in) :: iguess, nec
  Character(*), Intent(in), Optional :: oldcpf
  Interface
     ! function to synchronize the flags at output
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
  Character(132) :: cpf
  ! assign nec and iguess flags to info
  info%nec = nec
  info%guess = iguess
  If (Present(oldcpf)) Then
     info%oldcpf = oldcpf
  Else
     info%oldcpf = ''
  End If
  If (info%oldcpf .Ne. '' .And. info%guess.Gt.1) Then
     ! check to make sure the files exist
     cpf = trl_pe_filename(info%oldcpf, info%my_pe, info%npes)
     Open(info%cpio, file=cpf, status='OLD', form='UNFORMATTED',&
          & iostat=info%stat)
     If (info%stat .Eq. 0) Then
        Close(info%cpio, iostat=info%stat)
        If (info%stat.Ne.0) info%stat = -9
     Else
        info%stat = -8
     End If
     info%stat = trl_sync_flag(info%mpicom, info%stat)
  Else
     info%stat = 0
  End If
End Subroutine trl_set_iguess
!!!
! next function provides an uniform way of printing information
! stored in TRL_INFO_T.  It needs to be called by all PEs.
!
! Arguments
! info  -- the TRL_INFO_T variable to be printed
! mvflop-- the number of floating-operations per MATVEC per PE.  This
!          information has to be supplied by user, otherwise related
!          entries are left blank in the print out.
!
! In parallel environment, when writing to standard outputd device, only
! PE0 will actually write its local summary information.
!
! *** NOTE *** MUST be called on all PEs at the same time ***
!!!
Subroutine trl_print_info(info, mvflop)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Integer, Intent(in), Optional :: mvflop
  Interface
     Function trl_pe_filename(base, my_rank, npe) Result(filename)
       Implicit None
       Integer, Intent(in) :: my_rank, npe
       Character(*), Intent(in) :: base
       Character(132) :: filename
     End Function trl_pe_filename
     Subroutine trl_time_stamp(iou)
       Implicit None
       Integer, Intent(in) :: iou
     End Subroutine trl_time_stamp
     Subroutine trl_g_sum(mpicom, n, xx, wrk)
       Implicit None
       Integer, Intent(in) :: mpicom, n
       Real(8), Intent(out) :: wrk(n)
       Real(8), Intent(inout) :: xx(n)
     End Subroutine trl_g_sum
  End Interface
  Real(8), Parameter :: zero=0.0D0, one=1.0D0
  Integer :: rate, iot, ierr
  Real(8) :: t_tot, t_op, t_orth, t_res, t_in, t_out, rinv
  Real(8) :: r_tot, r_op, r_orth, r_res, r_in, r_out
  Real(8) :: tmp1(12), tmp2(12)
  Character(132) :: filename
  If (info%clk_rate.Gt.0) Then
     rinv = one / Dble(info%clk_rate)
  Else
     ! get clock rate
     Call System_clock(count_rate=rate)
     rinv = one / Dble(rate)
  Endif
  iot = info%log_io
  If (iot.Le.0 .Or. info%verbose.Le.0) iot = 6
  t_op  = (info%tick_o+info%clk_op) * rinv
  t_tot = (info%tick_t+info%clk_tot) * rinv
  t_res  = (info%tick_r+info%clk_res) * rinv
  t_orth = (info%tick_h+info%clk_orth) * rinv
  t_in = info%clk_in * rinv
  t_out = info%clk_out * rinv
  If (t_op.Ne.zero .And. Present(mvflop)) Then
     If (mvflop .Gt. 0) Then
        r_op = mvflop * info%matvec
     Else
        r_op = zero
     End If
  Else
     r_op = zero
  End If
  If (t_orth .Ne. zero) Then
     r_orth = info%rflp_h + info%flop_h
  Else
     r_orth = zero
  End If
  If (t_res .Ne. zero) Then
     r_res = info%rflp_r+info%flop_r
  Else
     r_res = zero
  End If
  If (r_op .Gt. zero) Then
     r_tot = Dble(mvflop)*Dble(info%matvec)+info%rflp+&
          &Dble(info%flop)
  Else
     r_tot = zero
  End If
  If (info%clk_in .Gt. 0) Then
     r_in = 8D0 * info%wrds_in
  Else
     r_in = zero
  End If
  If (info%clk_out .Gt. 0) Then
     r_out = 8D0 * info%wrds_out
  Else
     r_out = zero
  End If
  tmp1(1) = t_tot
  tmp1(2) = t_op
  tmp1(3) = t_orth
  tmp1(4) = t_res
  tmp1(5) = t_in
  tmp1(6) = t_out
  tmp1(7) = r_tot
  tmp1(8) = r_op
  tmp1(9) = r_orth
  tmp1(10) = r_res
  tmp1(11) = r_in
  tmp1(12) = r_out
  Call trl_g_sum(info%mpicom, 12, tmp1, tmp2)
  If (iot.Eq.info%log_io .And. iot.Ne.6) Then
     filename = trl_pe_filename(info%log_file, info%my_pe, info%npes)
     Open(iot, file=filename, status='UNKNOWN', position='APPEND',&
          & action='WRITE', iostat=ierr)
     If (ierr .Ne. 0) iot = 6
  End If
  If (iot.Eq.6 .And. info%my_pe.Gt.0) Return
  Call trl_time_stamp(iot)
  If (info%npes .Gt. 1) Then
     Write (iot, *) 'TRLAN execution summary (exit status =', info%stat,&
          & ') on PE ', info%my_pe
  Else
     Write (iot, *) 'TRLAN execution summary (exit status =', info%stat,&
          & ')'
  End If
  If (info%lohi.Gt.0) Then
     Write (iot, FMT=100) 'LARGEST', info%nec, info%ned
  Else If (info%lohi.Lt.0) Then
     Write (iot, FMT=100) 'SMALLEST', info%nec, info%ned
  Else
     Write (iot, FMT=100) 'EXTREME', info%nec, info%ned
  End If
  Write (iot, FMT=200) 'Times the operator is applied:', info%matvec,&
       & info%maxmv
  Write (iot, FMT=300) info%nloc, info%my_pe, info%ntot
  Write (iot, FMT=400) 'Convergence tolerance:', info%tol,&
       & info%tol*info%anrm
  Write (iot, FMT=500) 'Maximum basis size:', info%maxlan
  Write (iot, FMT=500) 'Restarting scheme:', info%restart
  Write (iot, FMT=500) 'Number of re-orthogonalizations:', info%north
  Write (iot, FMT=500) 'Number of (re)start loops:', info%nloop
  If (info%nrand .Gt. 0) Then
     Write (iot, FMT=500) 'Number of random vectors used:', info%nrand
  Endif
  If (info%npes .Gt. 1) Then
     Write (iot, FMT=500) 'Number of MPI processes:', info%npes
  Endif
  Write (iot, FMT=500) 'Number of eigenpairs locked:', info%locked
  If (t_op .Gt. zero) Then
     Write (iot, FMT=610) 'OP(MATVEC):', t_op, r_op/t_op, &
          & r_op
  Else
     Write (iot, FMT=600) 'time in OP:', t_op
  End If
  If (t_orth .Gt. zero) Then
     Write (iot, FMT=610) 'Re-Orthogonalization:', t_orth, &
          & r_orth/t_orth, r_orth
  Else
     Write (iot, FMT=600) 'time in orth:', t_orth
  End If
  If (t_res .Gt. zero) Then
     Write (iot, FMT=610) 'Restarting:', t_res, r_res/t_res,&
          & r_res
  Else
     Write (iot, FMT=600) 'time in restarting:', t_res
  End If
  If (t_tot .Gt. zero) Then
     Write (iot, FMT=610) 'TRLAN on this PE:', t_tot, r_tot/t_tot, &
          & r_tot
  Else
     Write (iot, FMT=600) 'total time in TRLAN:', t_tot
  End If
  If (info%verbose.Gt.0 .And. info%log_io.Ne.iot) Then
     Write(iot, *) 'Debug infomation written to files ',&
          & info%log_file(1:Index(info%log_file, ' ')-1), '####'
  Endif
  If (info%guess.Gt.1 .And. info%wrds_in.Gt.0) Then
     If (info%oldcpf .Ne. '') Then
        Write(iot,*) 'TRLAN restarted with checkpoint files ',&
             & info%oldcpf(1:Index(info%oldcpf,' ')-1), '####'
     Else
        Write(iot,*) 'TRLAN restarted with checkpoint files ',&
             & info%cpfile(1:Index(info%cpfile,' ')-1), '####'
     End If
     Write(iot,700) 'read', r_in, t_in, r_in/t_in
  End If
  If (info%clk_out.Gt.0 .And. info%wrds_out.Gt.0) Then
     Write(iot,*) 'Checkpoint files are ', &
          & info%cpfile(1:Index(info%cpfile, ' ')-1), '####'
     Write(iot,700) 'written', r_out, t_out, r_out/t_out
  End If
100 Format('Number of ', A, ' eigenpairs:', T32, I10, &
         &' (computed),', I11, ' (wanted)')
200 Format(A, T32, I10, ' (MAX:', I17, ' ).')
300 Format('Problem size: ', T32, I10, ' (PE:', I4, '),', I12,&
         & ' (Global)')
400 Format(A, T32, 1P, E10.3, ' (rel),', E16.3, ' (abs)')
500 Format(A, T32, I10)
600 Format(A, T25, 1P, E17.5, ' sec')
610 Format(A, T23, 1P, E12.4, ' sec,', E12.4, ' FLOP/S (', E11.4, ' FLOP)')
700 Format('Bytes ', A7, ': ', 1P, E12.5, ', Time(sec): ', E12.5,&
         & ', Rate(B/s): ', E12.5)
  If (info%npes .Gt. 1) Then
     ! write global performance information
     rinv = one / info%npes
     tmp1(1:12) = tmp1(1:12) * rinv
     Where (tmp1(1:6) .Gt. 0)
        tmp1(7:12) = tmp1(7:12) / tmp1(1:6)
     Elsewhere
        tmp1(7:12) = 0.0D0
     End Where
     If (tmp1(5).Eq.tmp1(6) .And. tmp1(5).Eq.zero) Then
        Write(iot, FMT=810)
        Write(iot, FMT=820) tmp1(1:4)
        Write(iot, FMT=830) tmp1(7:10)
     Else
        Write(iot, FMT=800)
        Write(iot, FMT=820) tmp1(1:6)
        Write(iot, FMT=830) tmp1(7:12)
     End If
  End If
800 Format(' -- Global summary --', /&
         &13X, 'Overall', 5X, 'MATVEC', 4X, 'Re-orth', 4X, 'Restart',&
         & 7X, 'Read', 6X, 'Write')
810 Format(' -- Global summary --', /&
         &13X, 'Overall', 5X, 'MATVEC', 4X, 'Re-orth', 4X, 'Restart')
820 Format('Time(ave)', 1P, 6E11.4)
830 Format('Rate(tot)', 1P, 6E11.4)
  If (iot.Eq.info%log_io .And. iot.Ne.6) Close(iot)
  Return
End Subroutine trl_print_info
!!!
! a more compact version of trl_print_info
!
! this is a local routine, indivadual PE can call it without regard of
! whether other PEs do the same and the output may be written to a
! different I/O unit number than the log_io
!!!
Subroutine trl_terse_info(info, iou)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Integer, Intent(in) :: iou
  Integer :: rate
  Real(8) :: t_tot, t_op, t_orth, t_res
  If (info%clk_rate.Gt.0) Then
     t_op  = (info%tick_o+info%clk_op) / Dble(info%clk_rate)
     t_tot = (info%tick_t+info%clk_tot) / Dble(info%clk_rate)
     t_res  = (info%tick_r+info%clk_res) / Dble(info%clk_rate)
     t_orth = (info%tick_h+info%clk_orth) / Dble(info%clk_rate)
  Else
     ! get clock rate
     Call System_clock(count_rate=rate)
     t_op  = (info%tick_o+info%clk_op) / Dble(rate)
     t_tot = (info%tick_t+info%clk_tot) / Dble(rate)
     t_res  = (info%tick_r+info%clk_res) / Dble(rate)
     t_orth = (info%tick_h+info%clk_orth) / Dble(rate)
  Endif
  If (info%lohi.Gt.0) Then
     Write (iou, FMT=100) info%maxlan, info%restart, '+', info%ned,&
          & info%nec
  Else If (info%lohi.Lt.0) Then
     Write (iou, FMT=100) info%maxlan, info%restart, '-', info%ned,&
          & info%nec
  Else
     Write (iou, FMT=100) info%maxlan, info%restart, '0', info%ned,&
          & info%nec
  End If
  Write (iou, FMT=200) info%matvec, info%north, info%nloop, info%locked
  If (t_tot.Gt.1.0D-3 .And. Max(t_tot,t_op,t_res,t_orth).Lt.1.0D3) Then
     Write (iou, FMT=400) t_tot, t_op, t_orth, t_res
  Else
     Write (iou, FMT=300) t_tot, t_op, t_orth, t_res
  End If
100 Format('MAXLAN:', I10, ', Restart:', I10, ',   NED:', 2X, A1,&
         & I7, ',     NEC:', I10)
200 Format('MATVEC:', I10, ',  Reorth:', I10, ', Nloop:', I10,&
         & ', Nlocked:', I10)
300 Format('Ttotal:', 1PE10.3, ',    T_op:', 1PE10.3, ', Torth:',&
         & 1PE10.3, ',  Tstart:', 1PE10.3)
400 Format('Ttotal:', F10.6, ',    T_op:', F10.6, ', Torth:',&
         & F10.6, ',  Tstart:', F10.6)
  Return
End Subroutine trl_terse_info
!!!
! trl_check_ritz
! performs a standard check on the computed Ritz pairs
!
! Arguments:
! op    -- the operator, aka the matrix-vector multiplication routine
! nrow  -- the problem size
! rvec  -- the array storing the Ritz vectors
! alpha -- The ritz values
! beta  -- The residual norms returned from a Lanczos routine (optional)
! eval  -- The actual eigenvalues (optional)
! wrk   -- workspace (optional)
!!!
Subroutine trl_check_ritz(op, info, nrow, rvec, alpha, beta, eval, wrk)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Integer, Intent(in) :: nrow
  Real(8), Dimension(:,:), Intent(in) :: rvec
  Real(8), Dimension(:), Intent(in) :: alpha
  Real(8), Dimension(:), Intent(in), Optional :: beta, eval
  Real(8), Dimension(:), Optional, Target :: wrk
  Interface
     Function trl_pe_filename(base, my_rank, npe) Result(filename)
       Implicit None
       Integer, Intent(in) :: my_rank, npe
       Character(*), Intent(in) :: base
       Character(132) :: filename
     End Function trl_pe_filename
     Subroutine op(nrow, ncol, xx, ldx, yy, ldy)
       Implicit None
       Integer, Intent(in) :: nrow, ncol, ldx, ldy
       Real(8), Intent(in) :: xx(ldx*ncol)
       Real(8), Intent(out) :: yy(ldy*ncol)
     End Subroutine op
     Subroutine trl_g_sum(mpicom, n, xx, wrk)
       Implicit None
       Integer, Intent(in) :: mpicom, n
       Real(8), Intent(out) :: wrk(n)
       Real(8), Intent(inout) :: xx(n)
     End Subroutine trl_g_sum
  End Interface
  ! local variables
  ! aq -- store the result of matrix-vector multiplication, size nrow
  ! rq -- store the Rayleigh-Quotient and the residual norms
  ! gsumwrk -- workspace left over for trl_g_sum to use
  Character(132) :: filename
  Integer :: i,lde,nev,lwrk,iaq,irq,iou
  Real(8) :: gapl, gapr
  Real(8), Dimension(:), Pointer :: aq, rq, gsumwrk, res, err
  ! dimension of the input arrays
  lde = Size(rvec, 1)
  nev = Size(rvec, 2)
  If (nev .Le. 0) Return  ! there is nothing to do
  If (nev .Gt. Size(alpha,1)) Then
     Print *, 'TRL_CHECK_RITZ: Rvec and Alpha array size do not match.'
     Return
  End If
  ! figure out whether it is necessary to allocate new workspace
  iaq = nrow
  irq = nev+nev+nev
  If (Present(wrk)) Then
     ! WRK provide by caller, try to use it
     lwrk = Size(wrk)
     If (lwrk .Ge. iaq+irq) Then
        aq => wrk(1:iaq)
        rq => wrk(iaq+1:iaq+irq)
        gsumwrk => wrk(iaq+irq+1:lwrk)
     Else If (lwrk .Ge. iaq) Then
        aq => wrk(1:iaq)
        gsumwrk => wrk(iaq+1:lwrk)
        Allocate(rq(nev+nev+nev), stat=irq)
        If (irq .Ne. 0) Then
           Print *, 'TRL_CHECK_RITZ: Failed to allocated workspace RQ.'
           Return
        End If
     Else If (lwrk .Ge. irq) Then
        rq => wrk(1:irq)
        gsumwrk => wrk(irq+1:lwrk)
        Allocate(aq(nrow), stat=iaq)
        If (iaq .Ne. 0) Then
           Print *, 'TRL_CHECK_RITZ: Failed to allocated workspace AQ.'
           Return
        End If
     Else
        gsumwrk => wrk
        Allocate(aq(nrow),  stat=iaq)
        If (iaq .Ne. 0) Then
           Print *, 'TRL_CHECK_RITZ: Failed to allocated workspace AQ.'
           Return
        End If
        Allocate(rq(nev+nev+nev), stat=irq)
        If (irq .Ne. 0) Then
           Print *, 'TRL_CHECK_RITZ: Failed to allocated workspace RQ.'
           Deallocate(aq)
           Return
        End If
     End If
  Else
     ! WRK not provided -- allocate space for AQ and RQ,
     ! gsumwrk points to the last third of RQ
     Allocate(aq(nrow),  stat=iaq)
     If (iaq .Ne. 0) Then
        Print *, 'TRL_CHECK_RITZ: Failed to allocated workspace AQ.'
        Return
     End If
     Allocate(rq(nev+nev+nev), stat=irq)
     If (irq .Ne. 0) Then
        Print *, 'TRL_CHECK_RITZ: Failed to allocated workspace RQ.'
        Deallocate(aq)
        Return
     End If
     gsumwrk => rq(nev+nev+1:nev+nev+nev)
  End If
  aq = 0
  rq = 0
  gsumwrk = 0
  !
  ! go through each Ritz pair one at a time, compute Rayleigh
  ! quotient and the corresponding residual norm
  res => rq(nev+1:nev+nev)
  Do i = 1, nev
     Call op(nrow, 1, rvec(:,i), lde, aq, nrow)
     ! Rayleigh quotient -- assuming rvec(:,i) has unit norm
     rq(i) = Dot_product(rvec(1:nrow,i), aq)
     Call trl_g_sum(info%mpicom, 1, rq(i:i), gsumwrk)
     aq = aq - rq(i)*rvec(1:nrow,i)
     res(i) = Dot_product(aq, aq)
  End Do
  Call trl_g_sum(info%mpicom, nev, rq(nev+1:nev+nev), gsumwrk)
  !
  ! compute the error estimate based on computed residual norms and the
  ! Ritz values
  err => rq(nev+nev+1:nev+nev+nev)
  res = Sqrt(res)
  gapl = alpha(nev) - alpha(1)
  Do i = 1, nev-1
     gapr = alpha(i+1) - alpha(i)
     gapl = Min(gapl, gapr)
     If (res(i).Ge.gapl) Then
        err(i) = res(i)
     Else
        err(i) = res(i) * res(i) / gapl
     End If
     gapl = gapr
  End Do
  If (res(nev).Ge.gapl) Then
     err(nev) = res(nev)
  Else
     err(nev) = res(nev) * res(nev) / gapl
  End If
  !
  ! if a log file was opened before, open it again to append the new
  ! stuff to the end
  iou = info%log_io
  If (iou.Le.0 .Or. info%verbose.Le.0) iou = 6
  If (iou .Ne. 6) Then
     filename = trl_pe_filename(info%log_file, info%my_pe, info%npes)
     Open(iou, file=filename, status='OLD', position='APPEND',&
          & action='WRITE', iostat=i)
     If (i .Ne. 0) iou = 6
  End If
  ! if writing to IO unit 6, only PE 0 does it
  If (iou.Eq.6 .And. info%my_pe.Gt.0) go to 888
  ! print out the information
  Write (iou, FMT=100)
  If (Present(beta) .And. Present(eval)) Then
     Do i = 1, nev
        Write(iou, FMT=200) alpha(i), res(i), beta(i)-res(i),&
             & err(i), rq(i)-alpha(i), eval(i)-alpha(i)
     End Do
  Else If (Present(beta)) Then
     Do i = 1, nev
        Write(iou, FMT=200) alpha(i), res(i),&
             & beta(i)-res(i), err(i), rq(i)-alpha(i)
     End Do
  Else If (Present(eval)) Then
     Do i = 1, nev
        Write(iou, FMT=300) alpha(i), res(i), err(i), rq(i)-alpha(i),&
             & eval(i)-alpha(i)
     End Do
  Else
     Do i = 1, nev
        Write(iou, FMT=300) alpha(i), res(i), err(i), rq(i)-alpha(i)
     End Do
  End If
100 Format('TRL_CHECK_RITZ:', /,&
         & 11X, 'Ritz value', 7X, 'res norm', '   res diff',&
         & '  est error', '  diff w rq', ' act. error')
200 Format(1P, G25.17, 5E11.3)
300 Format(1P, G25.17, E11.3, 11X, 3E11.3)
  If (iou.Ne.6) Close(iou)
  ! deallocate workspace allocated in this routine
888 If (irq .Eq. 0) Deallocate(rq)
  If (iaq .Eq. 0) Deallocate(aq)
End Subroutine trl_check_ritz
!!!
! A separate Rayleigh-Ritz projection routine
! Given a set of approximately orthonormal vectors (V), this routine
! performs the following operations
! (1) V'*V ==> G
! (2) R'*R :=  G
! (3) V'*A*V => H1, inv(R')*H1*inv(R) => H
! (4) Y*D*Y' := H
! (5) V*inv(R)*Y => V, diag(D) => lambda,
!     r(i) = ||A*V(:,i)-lambda(i)*V(:,i)||
!
! Arguments:
! op   -- the operator (matrix-vector multiplication routine)
! info -- the structure for storing information about TRLAN
! evec -- the Ritz vectors to be manipulated
! eres -- the array to store new Ritz values and residual norms
! base -- the (optional) workspace to store the result of MATVEC
! wrk  -- the (optional) workspace to store projection matrix, etc..
!!!
Subroutine trl_ritz_projection(op, info, evec, eres, wrk, base)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T) :: info
  Real(8), Dimension(:), Optional, Target :: wrk, base
  Real(8), Dimension(:, :), Target :: evec
  Real(8), Dimension(:) :: eres
  !! local variables
  Real(8), Parameter :: one=1.0D0, zero=0.0D0
  Integer :: i, ierr, lde, lwrk, mev, nev, nsqr, nrow, iuau, irvv
  Real(8), Dimension(:), Pointer :: rvv, uau, wrk2, avec
  External dgemm, dpotrf, dstev, dtrtrs
  Interface
     Subroutine op(nrow, ncol, xx, ldx, yy, ldy)
       Implicit None
       Integer, Intent(in) :: nrow, ncol, ldx, ldy
       Real(8), Intent(in) :: xx(ldx*ncol)
       Real(8), Intent(out) :: yy(ldy*ncol)
     End Subroutine op
     Subroutine trl_g_sum(mpicom, n, xx, wrk)
       Implicit None
       Integer, Intent(in) :: mpicom, n
       Real(8), Intent(out) :: wrk(n)
       Real(8), Intent(inout) :: xx(n)
     End Subroutine trl_g_sum
     Subroutine trl_ritz_vectors(nrow, lck, ny, yy, ldy, vec1, ld1, m1,&
          & vec2, ld2, m2, wrk, lwrk)
       Implicit None
       Integer, Intent(in) :: nrow, ny, lck, ldy, ld1, m1, m2, ld2, lwrk
       Real(8), Intent(in) :: yy(ldy, ny)
       Real(8) :: vec1(ld1,m1), vec2(ld2,m2), wrk(lwrk)
     End Subroutine trl_ritz_vectors
  End Interface
  !
  lde = Size(evec,1)
  mev = Size(evec,2)
  nrow = info%nloc
  If (info%nec .Gt. 0) Then
     nev = info%nec + 1
  Else
     nev = Min(info%ned, mev-1)
     If (info%lohi.Ne.0) nev = nev + 1
  End If
  nsqr = nev*nev
  If (Present(wrk)) Then
     lwrk = Size(wrk, 1)
  Else
     lwrk = 0
  End If
  If (Present(base)) Then
     avec => base(1:nrow)
  Else If (mev .Gt. nev) Then
     avec => evec(1:nrow, mev)
  Else
     Allocate(avec(nrow))
  End If
  ! memory allocation -- need 3*nev*nev elements, will allocate them
  ! in two consecutive blocks, uau(nev*nev), rvv(2*nev*nev)
  ! in actual use, rvv is further split in two until the last operation
  iuau = nsqr
  irvv = nsqr+nsqr
  If (lwrk .Ge. iuau+irvv) Then
     uau => wrk(1:nsqr)
     rvv => wrk(nsqr+1:nsqr+nsqr)
     wrk2 => wrk(nsqr+nsqr+1:lwrk)
  Else If (lwrk .Ge. irvv) Then
     rvv => wrk(1:nsqr)
     wrk2 => wrk(nsqr+1:lwrk)
     Allocate(uau(nsqr), stat=iuau)
     If (iuau .Ne. 0) Then
        info%stat = -231
        Goto 888
     End If
  Else If (lwrk .Ge. iuau) Then
     uau => wrk(1:nsqr)
     Allocate(rvv(nsqr+nsqr), stat=irvv)
     If (irvv .Ne. 0) Then
        info%stat = -232
        Goto 888
     End If
     wrk2 => rvv(nsqr+1:nsqr+nsqr)
  Else
     Allocate(uau(nsqr), stat=iuau)
     If (iuau .Ne. 0) Then
        info%stat = -231
        Goto 888
     End If
     Allocate(rvv(nsqr+nsqr), stat=irvv)
     If (irvv .Ne. 0) Then
        info%stat = -232
        Goto 888
     End If
     wrk2 => rvv(nsqr+1:nsqr+nsqr)
  End If
  ! step (1) : V'*V ==> G
  Call dgemm('T', 'N', nev, nev, nrow, one, evec, lde, evec, lde, zero,&
       & rvv, nev)
  Call trl_g_sum(info%mpicom, nsqr, rvv(1:nsqr), wrk2)
  ! step (2) : Choleskey factorization of G
  Call dpotrf('U', nev, rvv, nev, ierr)
  If (ierr .Ne. 0) Then
     info%stat = -234
     Goto 888
  End If
  ! step (3) : compute H_1 = V'*A*V
  ! use the first nrow elements of avec to store the results of
  ! matrix-vector multiplication
  wrk2 = zero
  Do i = 1, nev
     Call op(nrow, 1, evec(1:nrow, i), lde, avec, nrow)
     Call dgemv('T', nrow, i, one, evec, lde, avec, 1, zero,&
          & wrk2((i-1)*nev+1), 1)
  End Do
  Call trl_g_sum(info%mpicom, nsqr, wrk2(1:nsqr), uau)
  Do i = 2, nev
     wrk2(i:(i-1)*nev:nev) = wrk2((i-1)*nev+1:(i-1)*nev+i-1)
  End Do
  ! compute solution of R^T H_2 = H_1
  Call dtrtrs('U', 'T', 'N', nev, nev, rvv, nev, wrk2, nev, ierr)
  If (ierr .Ne. 0) Then
     info%stat = -235
     Goto 888
  End If
  ! compute solution of R^T H = H_2^T
  Do i = 1, nev
     uau(i:nsqr:nev) = wrk2((i-1)*nev+1:i*nev)
  End Do
  Call dtrtrs('U', 'T', 'N', nev, nev, rvv, nev, uau, nev, ierr)
  If (ierr .Ne. 0) Then
     info%stat = -236
     Goto 888
  End If
  ! solve the small symmetric eigenvalue problem
  Call dsyev('V', 'U', nev, uau, nev, eres, wrk2, nsqr, ierr)
  If (ierr .Ne. 0) Then
     info%stat = -237
     Goto 888
  End If
  ! solve R Y = Y to prepare for multiplying with V
  Call dtrtrs('U', 'N', 'N', nev, nev, rvv, nev, uau, nev, ierr)
  If (ierr .Ne. 0) Then
     info%stat = -238
     Goto 888
  End If
  ! call trl_ritz_vector to do the final multiplication
  If (lwrk .Ge. 3*nsqr) Then
     wrk2 => wrk(nsqr+1:lwrk)
  Else If (lwrk .Ge. nsqr+nsqr) Then
     wrk2 => wrk
  Else
     wrk2 => rvv
  End If
  i = Size(wrk2)
  Call trl_ritz_vectors(nrow, 0, nev, uau, nev,&
       & evec, lde, nev, avec, nrow, 0, wrk2, i)
  ! compute the residual norms
  Do i = 1, nev
     Call op(nrow, 1, evec(1:nrow,i), lde, avec, nrow)
     avec = avec - eres(i)*evec(1:nrow,i)
     eres(nev+i) = Dot_product(avec, avec)
  End Do
  Call trl_g_sum(info%mpicom, nev, eres(nev+1:nev+nev), avec)
  Do i = nev+1, nev+nev
     If (eres(i) .Gt. 0.0D0) Then
        eres(i) = Sqrt(eres(i))
     Else
        eres(i) = -Huge(eres(1))
     End If
  End Do
  If (info%lohi .Lt. 0) Then
     Do i = nev, nev+nev-2
        eres(i) = eres(i+1)
     End Do
  Else If (info%lohi .Gt. 0) Then
     Do i = 1, nev-1
        eres(i) = eres(i+1)
        evec(:,i) = evec(:,i+1)
     End Do
     Do i = nev, nev+nev-2
        eres(i) = eres(i+2)
     End Do
  End If
888 If (iuau .Eq. 0) Deallocate(uau)
  If (irvv .Eq. 0) Deallocate(rvv)
  If ((.Not. Present(base)) .And. mev.Le.nev) Deallocate(avec)
End Subroutine trl_ritz_projection
!!!
! a subroutine to compute Rayleigh quotients --
! Given a set of Ritz vectors and Ritz values, normalize the Ritz
! vectors and compute their Rayleigh quotients to replace the Ritz
! values.
!
! Arguments:
! op   -- the matrix-vector multiplication routine
! info -- the data structure that stores infomation for TRLAN
! evec -- the array to store the portion of eigenvectors on this PE
! base -- the workspace used to store results of MATVEC
! eres -- store new Ritz values and new Rresidual norms
!         if there are NEV Ritz pairs, eres(1:NEV) stores the new
!         Rayleigh quotient and eres(nev+1:nev+nev) stores the
!         new residual norms.
!!!
Subroutine trl_rayleigh_quotients(op, info, evec, eres, base)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T) :: info
  Real(8), Dimension(:) :: eres
  Real(8), Dimension(:,:) :: evec
  Real(8), Dimension(:), Target, Optional :: base
  Interface
     Subroutine op(nrow, ncol, xx, ldx, yy, ldy)
       Implicit None
       Integer, Intent(in) :: nrow, ncol, ldx, ldy
       Real(8), Intent(in) :: xx(ldx*ncol)
       Real(8), Intent(out) :: yy(ldy*ncol)
     End Subroutine op
     Subroutine trl_g_sum(mpicom, n, xx, wrk)
       Implicit None
       Integer, Intent(in) :: mpicom, n
       Real(8), Intent(out) :: wrk(n)
       Real(8), Intent(inout) :: xx(n)
     End Subroutine trl_g_sum
  End Interface
  ! local variables
  Integer :: i, nev, nrow
  Real(8), Dimension(4) :: wrk
  Real(8), Dimension(:), Pointer :: avec
  nev = Size(evec, 2)
  nrow = info%nloc
  If (nev .Le. 0) Return
  If (Present(base)) Then
     avec => base(1:nrow)
  Else
     Allocate(avec(nrow))
  End If
  avec = 0
  ! loop through each vector to normalize the vector, compute Rayleigh
  ! quotient and compute residual norm of the new Ritz pairs
  Do i = 1, nev
     wrk(1) = Dot_product(evec(1:nrow,i), evec(1:nrow,i))
     Call op(nrow, 1, evec(1:nrow,i), nrow, avec, nrow)
     wrk(2) = Dot_product(evec(1:nrow,i), avec)
     Call trl_g_sum(info%mpicom, 2, wrk(1:2), wrk(3:4))
     info%matvec = info%matvec + 1
     info%flop = info%flop + 4*nrow
     If (wrk(1) .Gt. 0.0D0) Then
        eres(i) = wrk(2)/wrk(1)
        avec = avec - eres(i)*evec(1:nrow,i)
        wrk(2) = Dot_product(avec, avec)
        Call trl_g_sum(info%mpicom, 1, wrk(2), wrk(3))
        wrk(1) = 1.0D0 / Sqrt(wrk(1))
        eres(nev+i) = wrk(1) * Sqrt(wrk(2))
        evec(1:nrow,i) = wrk(1) * evec(1:nrow,i)
        info%flop = info%flop + 6*nrow + 3
     Else
        eres(i) = -Huge(wrk(1))
        eres(nev+i) = -Huge(wrk(1))
     End If
  End Do
  If (.Not. Present(base)) Deallocate(avec)
End Subroutine trl_rayleigh_quotients
