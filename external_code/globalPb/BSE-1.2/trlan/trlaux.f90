! $Id: trlaux.f90,v 1.1 2004/05/04 06:29:35 fowlkes Exp $
!!!
! This file contains most auxillary routines which are not extensively
! used in during the normal operations of the TRLan, e.g., printing,
! error checking, etc.
!!!
! print an integer array for debugging
!!!
Subroutine trl_print_int(info, title, array)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Character(*), Intent(in) :: title
  Integer, Dimension(1:), Intent(in) :: array
  If (Size(array).Gt.3) Then
     Write (info%log_io, *) 'PE', info%my_pe, ': ',&
          & Trim(title)
     Write (info%log_io, FMT=100) array
  Else
     Write (info%log_io, *) 'PE', info%my_pe, ': ',&
          & Trim(title), array
  End If
100 Format(8I10)
End Subroutine trl_print_int
!!!
! print a real(8) array for debugging
!!!
Subroutine trl_print_real(info, title, array)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Character(*), Intent(in) :: title
  Real(8), Dimension(1:), Intent(in) :: array
  If (Size(array).Gt.3) Then
     Write (info%log_io, *) 'PE', info%my_pe, ': ',&
          & Trim(title)
     Write (info%log_io, FMT=100) array
  Else If (Size(array).Gt.1) Then
     Write(info%log_io, FMT=*) 'PE', info%my_pe, ': ',&
          & Trim(title)
     Write (info%log_io, FMT=100) array
  Else
     Write(info%log_io, FMT=*) 'PE', info%my_pe, ': ',&
          & Trim(title), array
  End If
100 Format(1P,8E10.3)
End Subroutine trl_print_real
!!!
! current progress of eigenvalue solution
!!!
Subroutine trl_print_progress(info)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Write(info%log_io, FMT=100) info%matvec, info%nloop, info%nec,&
       & info%north, info%nrand, info%stat
  Write(info%log_io, FMT=200) info%trgt, info%tres, info%crat
100 Format('MATVEC:', I10, ',    Nloop:', I10, ',     Nec:', I10, /,&
       & 'Reorth:', I10, ',    Nrand:', I10, ',    Ierr:', I10)
200 Format('Target:', 1P,E10.3, ',   ResNrm:', 1P,E10.3, ',   CFact:',&
         & 1P,E10.3)
End Subroutine trl_print_progress
!!!
! check orthogonality of the basis given
! -- used to debug of TRLANCZOS
!!!
Subroutine trl_check_orth(info, nrow, v1, ld1, j1, v2, ld2, j2, wrk, lwrk)
  Use trl_info
  Implicit None
  Integer, Intent(in) :: nrow, ld1, ld2, j1, j2, lwrk
  Type(TRL_INFO_T), Intent(in) :: info
  Real(8), Intent(in) :: v1(ld1,j1), v2(ld2,j2)
  Real(8) :: wrk(lwrk)
  Interface
     Subroutine trl_g_dot(mpicom, nrow, v1, ld1, m1, v2, ld2, m2, rr, wrk)
       Implicit None
       Integer, Intent(in) :: mpicom, nrow, ld1, ld2, m1, m2
       Real(8), Intent(in) :: v1(ld1,m1), v2(ld2,m2), rr(nrow)
       Real(8), Intent(out) :: wrk(m1+m2+m1+m2)
     End Subroutine trl_g_dot
  End Interface
  !
  ! local variables
  !
  Real(8), Parameter :: one=1.0D0, zero=0.0D0
  Integer :: i, j, jnd
  Real(8) :: nrmfro, nrminf
  !
  jnd = j1 + j2
  nrmfro = zero
  nrminf = zero
  If (jnd .Le. 0) Return
  If (lwrk .Lt. (jnd+jnd)) Then
     Write(info%log_io, *) 'TRL_CHECK_ORTH: workspace too small.'
     Return
  End If
  Write(info%log_io, *) 'TRL_CHECK_ORTH: check orthogonality of ',&
       & jnd, ' basis vectors.'
  !
  ! check orthognality of the basis vectors
  !
  Do i = 1, Min(j1, info%ntot) ! check no more than ntot vectors
     Call trl_g_dot(info%mpicom, nrow, v1, ld1, i, v2, ld2, 0, &
          & v1(1,i), wrk)
     wrk(i) = wrk(i) - one
     If (info%verbose.Gt.7) Then
        Write(info%log_io, *) 'Orthogonality level of v(', i,&
             & ') ..'
        Write(info%log_io, '(1P,8E10.3)') wrk(1:i)
     End If
     nrmfro = nrmfro + 2*Dot_product(wrk(1:i-1), wrk(1:i-1)) + &
          & wrk(i)*wrk(i)
     wrk(i+1) = Maxval(Abs(wrk(1:i)))
     nrminf = Max(nrminf, wrk(i+1))
  End Do
  Do i = 1, Min(j2, info%ntot-j1) ! check no more than ntot vectors
     j = j1 + i
     Call trl_g_dot(info%mpicom, nrow, v1, ld1, j1, v2, ld2, i,&
          & v2(1,i), wrk)
     wrk(j) = wrk(j) - one
     If (info%verbose.Gt.7) Then
        Write(info%log_io, *) 'Orthogonality level of v(', j, ') ..'
        Write(info%log_io, '(1P,8E10.3)') wrk(1:j)
     End If
     nrmfro = nrmfro + 2*Dot_product(wrk(1:j-1), wrk(1:j-1)) + &
          & wrk(j)*wrk(j)
     nrminf = Max(nrminf, Maxval(Abs(wrk(1:j))))
  End Do
  Write(info%log_io, FMT=100) info%matvec, jnd, Sqrt(nrmfro)
  Write(info%log_io, FMT=200) nrminf
100 Format('Frobenius norm of orthogonality level ', I10, I4, 1P, E14.5)
200 Format('Maximum abs. value of orthogonality level is ', 1P, E14.5)
End Subroutine trl_check_orth
!!!
! check Lanczos recurrence relation for debug purpose
!!!
Subroutine trl_check_recurrence(op, info, nrow, v1, ld1, m1, v2, ld2,&
     & m2, kept, alpha, beta, wrk, lwrk)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Integer, Intent(in) :: nrow, ld1, ld2, m1, m2, kept, lwrk
  Real(8), Intent(in), Target :: v1(ld1,m1), v2(ld2,m2)
  Real(8), Intent(in) :: alpha(m1+m2-1), beta(m1+m2-1)
  Real(8), Target :: wrk(lwrk)
  External op
  Interface
     Subroutine trl_print_real(info, title, array)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Character*(*), Intent(in) :: title
       Real(8), Dimension(1:), Intent(in) :: array
     End Subroutine trl_print_real
  End Interface
  !
  ! local variables
  !
  Integer :: i, ii, j, j1, j2, jnd, mv1
  Character*80 :: title
  Real(8), Dimension(:), Pointer :: aq, qkp1, cs, alf, bet
  Real(8), Parameter :: zero=0.0D0, one=1.0d0
  !
  mv1 = m1
  If (m2 .Gt. 0) Then
     j2 = m2 -1
     j1 = m1
  Else
     j2 = 0
     j1 = m1 - 1
  End If
  jnd = j1 + j2
  If (lwrk .Lt. jnd*4+Max(jnd*4,nrow)) Then
     Write(info%log_io, *)&
          & 'TRL_CHECK_RECURRENCE: not enough workspace.'
     Return
  End If
  If (lwrk .Ge. jnd*4+nrow) Then
     aq => wrk(lwrk-nrow+1:lwrk)
  Else If (lwrk .Ge. jnd*4) Then
     Allocate(aq(nrow), stat=i)
     If (i.Ne.0) Then
        Write(info%log_io, *)&
             & 'TRL_CHECK_RECURRENCE: failed to allcoate workspace.'
        Return
     End If
  End If
  wrk(1:jnd*4) = zero
  cs => wrk(jnd+1:jnd+jnd)
  alf => wrk(jnd+jnd+1:jnd+jnd+jnd)
  bet => wrk(jnd+jnd+jnd+1:jnd*4)
  !
  ! first type of relation
  ! A q_i = Alpha_i q_i + Beta_i q_{k+1}
  !
  If (kept.Lt.Size(v1,2)) Then
     qkp1 => v1(1:nrow, kept+1)
  Else
     qkp1 => v2(1:nrow, kept-j1+1)
  End If
  Do i = 1, Min(j1, kept)
     Call op(nrow, 1, v1(1,i), nrow, aq, nrow)
     Do ii = 1, nrow
        alf(i) = alf(i) + aq(ii) * v1(ii,i)
        aq(ii) = aq(ii) - alpha(i) * v1(ii,i)
        bet(i) = bet(i) + aq(ii) * aq(ii)
        cs(i)  = cs(i)  + aq(ii) * qkp1(ii)
        aq(ii) = aq(ii) - beta(i) * qkp1(ii)
        wrk(i) = wrk(i) + aq(ii) * aq(ii)
     End Do
  End Do
  Do i = 1, kept-j1
     j = i + j1
     Call op(nrow, 1, v2(1,i), nrow, aq, nrow)
     Do ii = 1, nrow
        alf(j) = alf(j) + aq(ii) * v2(ii, i)
        aq(ii) = aq(ii) - alpha(j) * v2(ii, i)
        bet(j) = bet(j) + aq(ii) * aq(ii)
        cs(j)  = cs(j)  + aq(ii) * qkp1(ii)
        aq(ii) = aq(ii) - beta(j) * qkp1(ii)
        wrk(j) = wrk(j) + aq(ii) * aq(ii)
     End Do
  End Do
  !
  ! the (k+1)st base vector need to orthogonalize against all previous
  ! vectors
  !
  If (jnd .Gt. kept) Then
     Call op(nrow, 1, qkp1, nrow, aq, nrow)
     alf(kept+1) = Dot_product(aq, qkp1)
     aq = aq - alpha(kept+1) * qkp1
     Do i = 1, Min(j1,kept)
        aq = aq - beta(i) * v1(1:nrow,i)
     End Do
     Do i = 1, kept-j1
        j = j1 + i
        aq = aq - beta(j) * v2(1:nrow,i)
     End Do
     bet(kept+1) = Dot_product(aq, aq)
     If (kept+2 .Le. j1) Then
        cs(kept+1) = Dot_product(aq, v1(1:nrow, kept+2))
        aq = aq - beta(kept+1)*v1(1:nrow, kept+2)
     Else
        cs(kept+1) = Dot_product(aq, v2(1:nrow, kept+2-j1))
        aq = aq - beta(kept+1)*v2(1:nrow, kept+2-j1)
     End If
     wrk(kept+1) = Dot_product(aq, aq)
  End If
  !
  ! the third kind of relation -- normal three term recurrence
  ! depending the fact that if the lower-bound of loop is less than
  ! upper bound, the look should not be executed
  !
  Do i = kept+2, j1
     Call op(nrow, 1, v1(1:nrow,i), nrow, aq, nrow)
     If (i.Lt.mv1) Then
        Do ii = 1, nrow
           alf(i) = alf(i) + aq(ii) * v1(ii,i)
           aq(ii) = aq(ii) - alpha(i)*v1(ii,i) - beta(i-1)*v1(ii,i-1)
           bet(i) = bet(i) + aq(ii) * aq(ii)
           cs(i)  = cs(i)  + aq(ii) * v1(ii, i+1)
           aq(ii) = aq(ii) - beta(i) * v1(ii, i+1)
           wrk(i) = wrk(i) + aq(ii) * aq(ii)
        End Do
     Else
        Do ii = 1, nrow
           alf(i) = alf(i) + aq(ii) * v1(ii,i)
           aq(ii) = aq(ii) - alpha(i)*v1(ii,i) - beta(i-1)*v1(ii,i-1)
           bet(i) = bet(i) + aq(ii) * aq(ii)
           cs(i)  = cs(i)  + aq(ii) * v2(ii,1)
           aq(ii) = aq(ii) - beta(i) * v2(ii,1)
           wrk(i) = wrk(i) + aq(ii) * aq(ii)
        End Do
     End If
  End Do
  Do i = Max(1, kept-j1+2), j2
     j = i + j1
     Call op(nrow, 1, v2(1,i), nrow, aq, nrow)
     If (i.Gt.1) Then
        Do ii = 1, nrow
           alf(j) = alf(j) + aq(ii) * v2(ii,i)
           aq(ii) = aq(ii) - beta(j-1)*v2(ii,i-1) - alpha(j)&
                &*v2(ii,i)
           bet(j) = bet(j) + aq(ii) * aq(ii)
           cs(j)  = cs(j)  + aq(ii) * v2(ii,i+1)
           aq(ii) = aq(ii) - beta(j) * v2(ii,i+1)
           wrk(j) = wrk(j) + aq(ii) * aq(ii)
        End Do
     Else
        Do ii = 1, nrow
           alf(j) = alf(j) + aq(ii) * v2(ii,1)
           aq(ii) = aq(ii) - beta(j-1)*v1(ii,j1) - alpha(j)*v2(ii,1)
           bet(j) = bet(j) + aq(ii) * aq(ii)
           cs(j)  = cs(j)  + aq(ii) * v2(ii,2)
           aq(ii) = aq(ii) - beta(j) * v2(ii,2)
           wrk(j) = wrk(j) + aq(ii) * aq(ii)
        End Do
     End If
  End Do
  !
  Call trl_g_sum(info%mpicom, jnd*4, wrk, wrk(jnd*4+1))
  aq(1) = Sqrt(Sum(wrk(1:jnd)))
  wrk(1:jnd) = Sqrt(wrk(1:jnd))
  Do ii =1, jnd
     If (bet(ii) .Gt. zero) Then
        bet(ii) = Sign(Sqrt(bet(ii)), beta(ii))
        cs(ii) = cs(ii) / bet(ii)
     Else
        bet(ii) = zero
     End If
  End Do
  title = 'Alpha computed by TRLAN ..'
  Call trl_print_real(info, title, alpha(1:jnd))
  title = 'Alpha computed explicitly in TRL_CHECK_RECURRENCE ..'
  Call trl_print_real(info, title, alf)
  title = 'Differences in alpha ..'
  alf = alf - alpha(1:jnd)
  Call trl_print_real(info, title, alf)
  title = 'Beta computed by TRLAN ..'
  Call trl_print_real(info, title, beta(1:jnd))
  title = 'Beta computed explicitly in TRL_CHECK_RECURRENCE ..'
  Call trl_print_real(info, title, bet)
  title = 'Differences in beta ..'
  bet = bet - beta(1:jnd)
  Call trl_print_real(info, title, bet)
  title = 'Error in Lanczos recurrence (overall) ='
  Call trl_print_real(info, title, aq(1:1))
  If (info%verbose.Gt.7) Then
     title = '|| A q_i - alpha_i q_i - beta_{i-1} q_{i-1} - beta_i q_{i&
          &+1} ||..'
     Call trl_print_real(info, title, wrk(1:jnd))
     title = '(A q_i - alpha_i q_i - beta_{i-1} q_{i-1})*q_{i+1}/beta_i ..'
     Call trl_print_real(info, title, cs(1:jnd))
     title = 'Sine of the angles ..'
     cs = cs*cs
     Do ii = 1, jnd
        If (cs(ii) .Lt. one) Then
           cs(ii) = Sqrt(one - cs(ii))
        Else
           cs(ii) = -one
        End If
     End Do
     Call trl_print_real(info, title, cs(1:jnd))
  End If
  If (lwrk.Lt.jnd*4+nrow) Deallocate(aq)
End Subroutine trl_check_recurrence
!!!
! write a check-point file
!!!
Subroutine trl_write_checkpoint(cp_io, filename, nrow, alpha, beta,&
     & evec, lde, me, base, ldb, nb, ierr)
  Implicit None
  Character(*), Intent(in) :: filename
  Integer, Intent(in) :: cp_io, nrow, ldb, lde, me, nb
  Real(8), Intent(in) :: alpha(me+nb-1), beta(me+nb-1)
  Real(8), Intent(in) :: evec(lde,me), base(ldb,nb)
  Integer, Intent(out) :: ierr
  ! local variables
  Integer :: jnd, i
  !
  ierr = 0
  jnd = me + nb - 1
  Open(cp_io, file=filename, form='UNFORMATTED', position='REWIND',&
       & status='REPLACE', action='WRITE', iostat=ierr)
  If (ierr.Ne.0) Then
     Print *, 'TRL_WRITE_CHECKPOINT: failed to open file: ',&
          & Trim(filename), '(', ierr, ')'
     ierr = -221
     Return
  End If
  Write(cp_io, err=100) nrow, jnd
  Write(cp_io, err=100) alpha(1:jnd)
  Write(cp_io, err=100) beta(1:jnd)
  Do i = 1, me
     Write(cp_io, err=100) evec(1:nrow,i)
  End Do
  Do i = 1, nb
     Write(cp_io, err=100) base(1:nrow,i)
  End Do
  Close(cp_io, err=200)
  Return
  ! handling write error
100 ierr = -222
  Close(cp_io, err=200)
  Return
  ! handle close error -- do nothing
200 If (ierr .Eq. 0) ierr = -223
  Return
End Subroutine trl_write_checkpoint
!!!
! read check-point file
!!!
Subroutine trl_read_checkpoint(cp_io, filename, nrow, evec, lde, mev,&
     & j1, base, ldb, nbas, j2, alpha, beta, ierr)
  Implicit None
  Character(*), Intent(in) :: filename
  Integer, Intent(in) :: cp_io, nrow, lde, mev, ldb, nbas
  Integer, Intent(out) :: j1, j2, ierr
  Real(8), Intent(out) :: alpha(mev+nbas-1), beta(mev+nbas-1)
  Real(8), Intent(out) :: base(ldb,nbas), evec(lde,mev)
  ! local variables
  Integer :: i
  ! array size parameters
  If (lde.Ge.nrow .And. ldb.Ge.nrow) Then
     ierr = 0
  Else
     Print *, 'TRL_READ_CHECKPOINT: leading dimensions too small.'
     ierr = -211
     Return
  End If
  ! open file
  Open(cp_io, file=filename, status='OLD', action='READ',&
       & form='UNFORMATTED', position='REWIND', iostat=ierr)
  If (ierr.Ne.0) Then
     Print *, 'TRL_READ_CHECKPOINT: failed to open check-point file ',&
          & Trim(filename), ' (', ierr, ')'
     ierr = -212
     Return
  End If
  ! read size information
  Read(cp_io, err=100) j1, j2
  If (j1 .Ne. nrow) Then
     Print *, 'TRL_READ_CHECKPOINT: Nrow mismatch.'
     ierr = -213
     Return
  End If
  If (j2 .Gt. Min(Size(alpha), Size(beta), mev+nbas-1)) Then
     Print *, 'TRL_READ_CHECKPOINT: MAXLAN too small.'
     ierr = -214
     Return
  End If
  ! can continue read all data
  Read(cp_io, err=100) alpha(1:j2)
  Read(cp_io, err=100) beta(1:j2)
  j1 = Min(mev, j2)
  j2 = j2 - j1
  If (j1 .Lt. mev) Then
     Do i = 1, j1+1
        Read(cp_io, err=100) evec(1:nrow,i)
     End Do
  Else
     Do i = 1, j1
        Read(cp_io, err=100) evec(1:nrow, i)
     End Do
     Do i = 1, j2+1
        Read(cp_io, err=100) base(1:nrow, i)
     End Do
  End If
  Close(cp_io, err=200)
  Return
  ! handle read error
100 ierr = -215
  Close(cp_io, err=200)
  Return
  ! handle close error
200 If (ierr.Eq.0) ierr = -216
  Return
End Subroutine trl_read_checkpoint
!!!
! a function to generate file name from a base name and the PE number
!
! the resulting filename is limited to 132 character long.  The
! processor number takes four spaces
!!!
Function trl_pe_filename(base, my_rank, npe) Result(filename)
  Implicit None
  Integer, Intent(in) :: my_rank, npe
  Character(*), Intent(in) :: base
  Character(132) :: filename
  Character(8) :: fmtint
  !! local variable
  Integer :: lead, ndig
  ndig = 1
  lead = npe
  Do While (lead .Gt. 9)
     lead = lead / 10
     ndig = ndig + 1
  End Do
  lead = Min(133-ndig, Index(base, ' '))
  filename = base(1:lead-1)
  If (ndig.Gt.0 .And. ndig.Lt.10) Then
     Write(fmtint, FMT=100) ndig, ndig
  Else If (ndig.Ge.10 .And. ndig.Lt.100) Then
     Write (fmtint, FMT=200) ndig, ndig
  Else
     Stop 'TRL_PE_FILENAME: to many PEs'
  End If
  Write(filename(lead:), FMT=fmtint) my_rank
100 Format('(I', I1, '.', I1,')')
200 Format('(I', I2, '.', I2,')')
End Function trl_pe_filename
!!!
! print the current date and time
!!!
Subroutine trl_time_stamp(iou)
  Implicit None
  Integer, Intent(in) :: iou
  ! local variables
  Character*(10) :: date, time, zone
  !
  Call Date_and_time(date, time, zone)  ! F90 intrinsic function
100 Format(T40, A5, '/', A2, '/', A2, 1X, A2, ':', A2, ':', A6,&
         & ' (', A3, ':', A2, ')')
  Write (iou, 100) date(1:4), date(5:6), date(7:8), time(1:2),&
       & time(3:4), time(5:10), zone(1:3), zone(4:5)
  Return
End Subroutine trl_time_stamp
