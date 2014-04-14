! $Id: restart.f90,v 1.1 2004/05/04 06:29:35 fowlkes Exp $
!!!  TRLAN Low level utility routines
! This file contains a number of routines related to decisions of how
! many Ritz pairs to save when restarting the Thick-Restart Lanczos
! method.  The subroutine trl_shuffle_eig is the main access point.
!!!
! The subroutine trl_shuffle_eig accomplishes two tasks:
! (1) sort the Ritz values so that those to be saved after
! restart are ordered in the front of the array,
! (2) determine the number of Ritz values to be saved.
! On return from this subroutine, the Ritz values are in ascending
! order so that DSTEIN can be used to compute the eigenvectors
!!!
Subroutine trl_shuffle_eig(nd, mnd, lambda, res, info, kept)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Integer, Intent(in) :: nd, mnd
  Integer :: kept
  Real(8) :: lambda(nd), res(nd)
  !
  Integer :: i, ncl, ncr, kl, kr, tind, minsep
  Real(8) :: bnd
  External dsort2
  Interface
     Subroutine trl_restart_fixed(nd, mnd, tind, lambda, res, info, kl, kr)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: nd, mnd, tind
       Integer, Intent(inout) :: kl, kr
       Real(8), Intent(in) :: lambda(nd), res(nd)
     End Subroutine trl_restart_fixed
     Subroutine trl_restart_scan(nd, res, info, kept, kl, kr)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: nd, kept
       Integer, Intent(inout) :: kl, kr
       Real(8), Intent(in) :: res(nd)
     End Subroutine trl_restart_scan
     Subroutine trl_restart_small_res(nd, mnd, tind, lambda, res, info, kl, kr)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: nd, mnd, tind
       Integer, Intent(inout) :: kl, kr
       Real(8), Intent(in) :: lambda(nd), res(nd)
     End Subroutine trl_restart_small_res
     Subroutine trl_restart_max_gap_ratio(nd, tind, kept, lambda, res,&
          & info, kl, kr)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: kept, nd, tind
       Integer, Intent(inout) :: kl, kr
       Real(8), Intent(in) :: lambda(nd), res(nd)
     End Subroutine trl_restart_max_gap_ratio
     Subroutine trl_restart_max_progress(nd, tind, kept, lambda, res,&
          & info, kl, kr)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: kept, nd, tind
       Integer, Intent(inout) :: kl, kr
       Real(8), Intent(in) :: lambda(nd), res(nd)
     End Subroutine trl_restart_max_progress
     Subroutine trl_restart_max_reduction(nd, tind, kept, lambda, res,&
          & info, kl, kr)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: kept, nd, tind
       Integer, Intent(inout) :: kl, kr
       Real(8), Intent(in) :: lambda(nd), res(nd)
     End Subroutine trl_restart_max_reduction
  End Interface
  !
  ! very small basis -- save the half with the smallest residual norms
  If (nd .Le. 5) Then
     Call dsort2(nd, res, lambda)
     If (nd.Gt.3) Then
        kept = 2
     Else If (nd.Gt.0) Then
        kept = 1
     Else
        kept = 0
     End If
     If (kept.Gt.1) Call dsort2(kept, lambda, res)
     Return
  End If
  !
  ! preparation for normal case, first determine what are converged
  ! ncl are converged from the left and ncr are converged from the
  ! right
  bnd = Min(info%tol, Epsilon(info%tol))*info%anrm
  ncr = 1
  ncl = nd
  i = nd
  Do While (i .Gt. 0) ! determine how many has converged from the right
     If (res(i) .Le. bnd) Then
        i = i - 1
     Else
        ncr = i + 1
        i = 0
     End If
  End Do
  i = 1
  Do While (i .Le. nd) ! determine how many has converged from the left
     If (res(i) .Le. bnd) Then
        i = i + 1
     Else
        ncl = i - 1
        i = nd + 1
     End If
  End Do
  kl = Max(1, ncl)
  kr = Min(nd, ncr)
  If (ncr .Gt. ncl) Then
     ! find the one that is closest to info%trgt
     tind = (kl+kr)/2
     Do While (lambda(tind).Ne.info%trgt .And. kr.Gt.kl)
        If (lambda(tind) .Lt. info%trgt) Then
           kl = tind + 1
           tind = (kl + kr) / 2
        Else If (lambda(tind) .Gt. info%trgt) Then
           kr = tind - 1
           tind = (kl + kr) / 2
        Else
           kl = tind
           kr = tind
        End If
     End Do
     ! assign kl to the largest index of lambda that is smaller than
     ! info%trgt
     If (lambda(tind) .Eq. info%trgt) Then
        kl = tind - 1
10      If (kl .Gt. 0) Then
           If (lambda(kl) .Eq. info%trgt) Then
              kl = kl - 1
              Goto 10
           End If
        End If
     ! assign kr to the smallest index of lambda that is greater than
     ! info%trgt
        kr = tind + 1
20      If (kr .Le. nd) Then
           If (lambda(kr) .Eq. info%trgt) Then
              kr = kr + 1
              Goto 20
           End If
        End If
     Else
        kl = tind - 1
        kr = tind + 1
     End If
     ! initial assignment of kl and kr
     If (info%lohi .Gt. 0) Then
        kr = kl
        kl = Min(ncl, Max(0, nd-info%ned))
     Else If (info%lohi .Lt. 0) Then
        kl = kr
        kr = Max(ncr, Min(nd-info%nec,info%ned+1))
     Else If (ncr-tind .Gt. tind-ncl) Then
        kl = kr
        kr = Max(ncr, Min(nd-info%nec,info%ned+1))
     Else
        kr = kl
        kl = Min(ncl, Max(0, nd-info%ned))
     End If
  Else
     ! all have converged, keep all -- should not happen
     kept = nd
     Return
  End If
  !
  ! normal cases, call subroutines to complete the tasks
  ! [1 .. kl] and [kr .. nd] are saved for later
  ! the initial values of kl and kr are simply ncl and ncr
  ! they are further analyzed according to the restarting strategy
  ! requested
  Select Case (info%restart)
  Case (1)
     ! fixed number beyond the currently converged ones
     Call trl_restart_fixed(nd, mnd, tind, lambda, res, info, kl, kr)
  Case (2)
     ! add the ones with smallest residual nroms
     Call trl_restart_small_res(nd, mnd, tind, lambda, res, info, kl, kr)
  Case (3)
     If (info%nloop .Gt. 0) Then
        ! maximize the gap ratio
        Call trl_restart_max_gap_ratio(nd, tind, kept, lambda, res, info,&
             & kl, kr)
     Else
        ! this is the first restart -- use trl_restart_fixed instead
        Call trl_restart_fixed(nd, mnd, tind, lambda, res, info, kl, kr)
     End If
  Case (4)
     If (info%nloop .Gt. 0) Then
        ! maximize [gap-ratio * (m-k)]
        Call trl_restart_max_progress(nd, tind, kept, lambda, res,&
             & info, kl, kr)
     Else
        ! this is the first restart -- use trl_restart_fixed instead
        Call trl_restart_fixed(nd, mnd, tind, lambda, res, info, kl, kr)
     End If
  Case (5)
     If (info%nloop .Gt. 0) Then
        ! maximize [sqrt(gap tatio) * (m-k)]
        Call trl_restart_max_reduction(nd, tind, kept, lambda, res,&
             & info, kl, kr)
     Else
        ! this is the first restart -- use trl_restart_fixed instead
        Call trl_restart_fixed(nd, mnd, tind, lambda, res, info, kl, kr)
     End If
  Case (6)
     ! progressively vary the thickness
     Call trl_restart_scan(nd, res, info, kept, kl, kr)
  Case default
     If (info%restart .Le. -info%ned) Then
        If (info%lohi .Ge. 0) Then
           kl = 0
           kr = Max(2, nd+info%restart)+1
        Else If (info%lohi .Lt. 0) Then
           kl = Min(-info%restart, nd-3)
           kr = nd+1
        Else
           kl = Min(nd-3, -info%restart)/2
           kr = nd - kl
        End If
     Else
        Call trl_restart_fixed(nd, mnd, tind, lambda, res, info, kl, kr)
     End If
  End Select
  !
  ! make sure kr > kl+minsep
  minsep = Max(3, nd/6, nd-6*info%ned)
  If (kr.Le.kl+minsep .Or. kl+nd-kr+minsep.Gt.mnd) Then
     If (ncl.Lt.kl .And. kl.Lt.kr .And. kr.Lt.ncr) Then
        kl = kl - 1
        kr = kr + 1
     Else If (info%lohi .Gt. 0) Then
        kr = Max(minsep, Min(nd/3, ncr-1))
        kl = 0
     Else If (info%lohi .Lt. 0) Then
        kl = Min(nd-minsep, Max((nd+nd)/3, ncl+1))
        kr = nd+1
     Else
        kl = (nd-minsep)/2-1
        kr = (nd-minsep+1)/2+1
     End If
  End If
  ! copy the (kr:nd) elements to (kl+1:kl+nd-kr+2)
  ! kr is temporarily decreased by 1 to make indices easier to compute
  kr = kr - 1
  Do i = 1, nd-kr
     lambda(kl+i) = lambda(kr+i)
     res(kl+i) = res(kr+i)
  End Do
  kept = kl + Max(0, nd-kr)
  Return
End Subroutine trl_shuffle_eig
!!!
! save fixed number of extra Ritz pairs
! January 99 -- added modification to ensure a minimal gap ratio is
! maintained
!!!
Subroutine trl_restart_fixed(nd, mnd, tind, lambda, res, info, kl, kr)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Integer, Intent(in) :: nd, mnd, tind
  Integer, Intent(inout) :: kl, kr
  Real(8), Intent(in) :: lambda(nd), res(nd)
  Interface
     Function trl_min_gap_ratio(info, nd, tind, res) Result(gamma)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: nd, tind
       Real(8), Intent(in) :: res(nd)
       Real(8) :: gamma
     End Function trl_min_gap_ratio
  End Interface
  !
  ! local variables
  Real(8), Parameter :: ten=1.0D1
  Integer :: extra, i, kl0, kr0, minsep
  Real(8) :: gmin
  ! the number of extra Ritz pairs to be saved
  kl0 = kl
  kr0 = kr
  extra = Nint((mnd-info%nec)*(0.4D0+0.1D0*info%ned/Dble(mnd)))
  If (extra.Gt.info%ned+info%ned .And. extra.Gt.5) Then
     gmin = Dble(mnd)/Dble(info%ned)
     extra = Nint((extra + (Log(gmin)*info%ned) * gmin) / (1.0D0 +&
          & gmin))
  End If
  minsep = Max(3, nd/5, nd-4*info%ned)
  gmin = trl_min_gap_ratio(info, nd, tind, res)
  If (info%lohi .Gt. 0) Then
     kr = Min(tind-1, kr) - extra
     kl = 0
  Else If (info%lohi .Lt. 0) Then
     kl = Max(tind+1, kl) + extra
     kr = nd+1
  Else
     extra = extra + 1
     kl = kl + extra/2
     kr = kr - extra/2
     i = 1
     Do While (kl.Gt.kl0 .And. kr.Lt.kr0 .And. i.Eq.1)
        If (ten*res(kl).Lt.res(kr)) Then
           If (res(kl+1).Lt.res(kr+1)) Then
              kl = kl + 1
              kr = kr + 1
           Else
              i = 0
           End If
        Else If (ten*res(kr).Lt.res(kl)) Then
           If (res(kr-1).Lt.res(kl-1)) Then
              kr = kr - 1
              kl = kl - 1
           Else
              i = 0
           End If
        Else
           i = 0
        End If
     End Do
  End If
  ! adjust kl and kr until the minimal gap ratio is satisfied
  Do While (kl+minsep.Lt.kr .And. gap_ratio(Max(1,kl),Min(kr,nd)).Lt.gmin)
     If (info%lohi .Lt. 0) Then
        kl = kl + 1
     Else If (info%lohi .Gt. 0) Then
        kr = kr - 1
     Else If (res(kl) .Lt. res(kr)) Then
        kl = kl + 1
     Else
        kr = kr + 1
     End If
  End Do
  ! make sure nearly degenerate Ritz pairs are included
  ! lambda(kl)+r(kl) > lambda(j) and 
  !                lambda(kl)-r(kl) > lambda(j)-r(j)
  If (info%lohi .Gt. 0) Then
     i = kr-1
     Do While (i.Gt.kl+minsep .And. lambda(kr)-res(kr).Lt.lambda(i) .And.&
          & lambda(kr)+res(kr).Lt.lambda(i)+res(i))
        i = i - 1
     End Do
     kr = i+1
  Else
     kl0 = kl
     i = kl + 1
     Do While (i.Lt.kr-minsep .And. lambda(kl)+res(kl).Gt.lambda(i) .And.&
          & lambda(kl)-res(kl).Gt.lambda(i)-res(i))
        i = i + 1
     End Do
     kl = i-1
  End If
Contains
  Function gap_ratio(i,j) Result(gamma)
    Integer, Intent(in) :: i, j
    Real(8) :: gamma
    gamma = (lambda(i) - lambda(tind)) / (lambda(j) - lambda(tind))
    Return
  End Function gap_ratio
End Subroutine trl_restart_fixed
!!!
! progressively vary the thickness to scan all possible values
!
! This subroutine determines the number Ritz vectors to be saved by
! adding some quantity to the current thickness.  If the thickness is
! larger than nd-2, it is reset to something smaller.
!!!
Subroutine trl_restart_scan(nd, res, info, kept, kl, kr)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Integer, Intent(in) :: nd, kept
  Integer, Intent(inout) :: kl, kr
  Real(8), Intent(in) :: res(nd)
  !
  ! local variables
  Real(8), Parameter :: ten=1.0D1
  Integer :: extra, i, kl0, kr0
  ! three cases have to be dealt with separately
  If (info%lohi .Lt. 0) Then
     kr = nd + 1
     kl = kept + Min(Max(info%nec,1), (nd-kept)/2)
     If (kl .Le. 1) Then
        If (nd.Gt.6) Then
           kl = nd/2
        Else If (nd.Gt.2) Then
           kl = 2
        End If
     Else If (kl+3 .Ge. nd) Then
        kl = info%nec + Min(info%ned, 10, (nd-info%ned)/2)
     End If
  Else If (info%lohi .Gt. 0) Then
     kl = 0
     kr = kept + Min(Max(info%nec,1), (nd-kept)/2)
     If (kr .Le. 1) Then
        If (nd .Gt. 6) Then
           kr = nd / 2
        Else If (nd .Gt. 2) Then
           kr = 2
        End If
     Else If (kr+3 .Gt. nd) Then
        kr = info%nec + Min(info%ned, 10, (nd-info%ned)/2)
     End If
     kr = nd - kr + 1
  Else
     kl0 = kl
     kr0 = kr
     extra = kept + Min(info%nec, (nd-kept)/2) + 1
     If (extra .Le. 1) Then
        If (nd .Gt. 6) Then
           extra = nd / 2
        Else If (nd .Gt. 2) Then
           extra = 2
        End If
     Else If (extra+3 .Gt. nd) Then
        extra = info%nec + Min(info%ned, 10, (nd-info%ned)/2)
     End If
     kl = Max(kl, extra/2)
     kr = Min(kr, nd+1-extra/2)
     i = 1
     Do While (kl.Gt.kl0 .And. kr.Lt.kr0 .And. i.Eq.1)
        If (ten*res(kl).Lt.res(kr)) Then
           If (res(kl+1).Lt.res(kr+1)) Then
              kl = kl + 1
              kr = kr + 1
           Else
              i = 0
           End If
        Else If (ten*res(kr).Lt.res(kl)) Then
           If (res(kr-1).Lt.res(kl-1)) Then
              kr = kr - 1
              kl = kl - 1
           Else
              i = 0
           End If
        Else
           i = 0
        End If
     End Do
  End If
End Subroutine trl_restart_scan
!!!
! save those that are close to converge (as close as the target)
!!!
Subroutine trl_restart_small_res(nd, mnd, tind, lambda, res, info, kl, kr)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Integer, Intent(in) :: nd, mnd, tind
  Integer, Intent(inout) :: kl, kr
  Real(8), Intent(in) :: lambda(nd), res(nd)
  Interface
     Function trl_min_gap_ratio(info, nd, tind, res) Result(gamma)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: nd, tind
       Real(8), Intent(in) :: res(nd)
       Real(8) :: gamma
     End Function trl_min_gap_ratio
  End Interface
  !
  ! local variables
  Integer :: extra, i, j, ii, kl0, kr0, minsep
  Real(8) :: acpt, resmax, gmin
  ! the number of extra Ritz pairs to be saved
  minsep = Max(3, nd/5, nd-4*info%ned)
  extra = Nint((mnd-info%nec)*(0.4D0+0.1D0*info%ned/Dble(mnd)))
  If (extra.Gt.info%ned+info%ned .And. extra.Gt.5) Then
     gmin = Dble(mnd)/Dble(info%ned)
     extra = Nint((extra + (Log(gmin)*info%ned) * gmin) / (1.0D0 +&
          & gmin))
  End If
  kl0 = kl
  kr0 = kr
  resmax = Maxval(res)
  acpt = resmax / res(tind)
  !
  ! determine the number of Ritz pairs that has to be saved
  If (info%lohi .Gt. 0) Then
     If (acpt.Lt.0.999D0 .And. acpt.Ge.0.0D0) Then
        ii = tind - 1
        acpt = Max(Sqrt(acpt)*res(tind), res(ii)+res(ii))
        acpt = Min(acpt, resmax)
        kr = ii - 1
        Do While (res(kr).Lt.acpt .And. kr.Gt.kl+3)
           kr = kr - 1
        End Do
     Else
        kr = kr0 - extra
     End If
     kr = Max(3, kr)
     kl = Min(kl, kr-2)
  Else If (info%lohi .Lt. 0) Then
     If (acpt.Lt.0.999D0 .And. acpt.Ge.0.0D0) Then
        ii = tind + 1
        acpt = Max(Sqrt(acpt)*res(tind), res(ii)+res(ii))
        acpt = Min(acpt, resmax)
        kl = ii + 1
        Do While (res(kl).Lt.acpt .And. kl.Lt.kr-3)
           kl = kl + 1
        End Do
     Else
        kl = kl0 + extra
     End If
     kl = Min(nd-3, kl)
     kr = Max(kr, kl+2)
  Else
     ! save whoever has smaller residual norms
     i = kl + 1
     j = kr - 1
     Do ii = 1, extra
        If (res(i) .Lt. res(j)) Then
           kl = i
           i = i + 1
        Else If (res(i) .Gt. res(j)) Then
           kr = j
           j = j - 1
        Else If (i .Le. nd-j) Then
           kl = i
           i = i + 1
        Else
           kr = j
           j = j - 1
        End If
     End Do
  End If
  ! adjust kl and kr until the minimal gap ratio is satisfied
  kl0 = kl
  kr0 = kr
  gmin = trl_min_gap_ratio(info, nd, tind, res)
  Do While (kl+minsep.Lt.kr .And. gap_ratio(Max(1,kl),Min(nd,kr)).Lt.gmin)
     If (info%lohi .Lt. 0) Then
        kl = kl + 1
     Else If (info%lohi .Gt. 0) Then
        kr = kr - 1
     Else If (res(kl) .Lt. res(kr)) Then
        kl = kl + 1
     Else
        kr = kr + 1
     End If
  End Do
  ! make sure nearly degenerate Ritz pairs are included
  ! lambda(kl)+r(kl) > lambda(j) and 
  !                lambda(kl)-r(kl) > lambda(j)-r(j)
  If (info%lohi .Gt. 0) Then
     i = kr0-1
     Do While (i.Gt.kl+minsep .And. lambda(kr)-res(kr).Lt.lambda(i) .And.&
          & lambda(kr)+res(kr).Lt.lambda(i)+res(i))
        i = i - 1
     End Do
     kr = Min(kr, i+1)
  Else
     i = kl0 + 1
     Do While (i.Lt.kr-minsep .And. lambda(kl)+res(kl).Gt.lambda(i) .And.&
          & lambda(kl)-res(kl).Gt.lambda(i)-res(i))
        i = i + 1
     End Do
     kl = Max(kl, i-1)
  End If
Contains
  Function gap_ratio(i,j) Result(gamma)
    Integer, Intent(in) :: i, j
    Real(8) :: gamma
    gamma = (lambda(i) - lambda(tind)) / (lambda(j) - lambda(tind))
    Return
  End Function gap_ratio
End Subroutine trl_restart_small_res
!!!
! search throught all pairs of (kl, kr) for the one with the maximum
! gap ratio for the next Ritz pair(target)
!
! This is an optimized version of the original version.  It only search
! through nd number once. (Only single loop!)
!!!
Subroutine trl_restart_max_gap_ratio(nd, tind, kept, lambda, res, info, kl, kr)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Integer, Intent(in) :: kept, nd, tind
  Integer, Intent(inout) :: kl, kr
  Real(8), Intent(in) :: lambda(nd), res(nd)
  Interface
     Subroutine trl_restart_search_range(nd, lambda, res, info, ncl,&
          & ncr, lohi, tind, klm, krm)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: nd, ncl, ncr, tind
       Integer, Intent(out) :: klm, krm, lohi
       Real(8), Intent(in) :: lambda(nd), res(nd)
     End Subroutine trl_restart_search_range
     Function trl_min_gap_ratio(info, nd, tind, res) Result(gamma)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: nd, tind
       Real(8), Intent(in) :: res(nd)
       Real(8) :: gamma
     End Function trl_min_gap_ratio
  End Interface
  !
  ! local variables
  Integer :: i, j, lohi, klm, krm, igap
  Real(8) :: bnd, tmp, gratio
  ! statement function for computing gap ratio
  gratio(i,j) = (lambda(i)-info%trgt)/(lambda(j)-info%trgt)
  ! determine the search range
  Call trl_restart_search_range(nd, lambda, res, info, kl, kr, lohi,&
       & tind, klm, krm)
  kl = klm
  kr = krm
  igap = Max(Min(nd-info%ned, Nint((krm-klm)*0.4D0)), 2)
  If (igap.Gt.2 .And. info%matvec.Gt.info%maxlan) Then
     If (info%clk_op+info%tick_o.Gt.1.0D1*(info%clk_orth&
          &+info%tick_h+info%clk_res+info%tick_r)) Then
        igap = Max(2, nd-kept-1)
     Else
        bnd = trl_min_gap_ratio(info, nd, tind, res)
        If (info%crat .Lt. bnd) igap = Max(2, nd-kept-1)
     End If
  End If
  If (kl+igap .Gt. kr) Return
  ! two cases depending on lohi
  If (lohi .Gt. 0) Then
     ! target is at the high end of spectrum
     bnd = gratio(kr, kl)
     Do i = klm, krm-igap
        j = i + igap
        tmp = gratio(j,i)
        If (tmp .Gt. bnd) Then
           kl = i
           kr = j
           bnd = tmp
        End If
     End Do
  Else
     bnd = gratio(kl, kr)
     Do i = klm, krm-igap
        j = i + igap
        tmp = gratio(i,j)
        If (tmp .Gt. bnd) Then
           kl = i
           kr = j
           bnd = tmp
        End If
     End Do
  End If
End Subroutine trl_restart_max_gap_ratio
!!!
! search for a pair (kl, kr) such that the reduction in residual norm
! of the target (info%trgt) will be the largest before next restart
! The merit function is [gap-ratio * (m-k)]
!!!
Subroutine trl_restart_max_progress(nd, tind, kept, lambda, res, info, kl, kr)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Integer, Intent(in) :: kept, nd, tind
  Integer, Intent(inout) :: kl, kr
  Real(8), Intent(in) :: lambda(nd), res(nd)
  Interface
     Subroutine trl_restart_search_range(nd, lambda, res, info, ncl,&
          & ncr, lohi, tind, klm, krm)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: nd, ncl, ncr, tind
       Integer, Intent(out) :: klm, krm, lohi
       Real(8), Intent(in) :: lambda(nd), res(nd)
     End Subroutine trl_restart_search_range
     Function trl_min_gap_ratio(info, nd, tind, res) Result(gamma)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: nd, tind
       Real(8), Intent(in) :: res(nd)
       Real(8) :: gamma
     End Function trl_min_gap_ratio
  End Interface
  !
  ! local variables
  Integer :: i, j, lohi, klm, krm, igap
  Real(8) :: tmp, ss, merit
  ! merit measure the factor of residual norm reduction
  merit(i,j) = (lambda(i)-info%trgt) * Abs(j-i) / (lambda(j)-info%trgt)
  ! determine the search range
  Call trl_restart_search_range(nd, lambda, res, info, kl, kr, lohi,&
       & tind, klm, krm)
  ! perform the brute-force search
  kl = klm
  kr = krm
  igap = Max(Min(nd-info%ned, Nint((krm-klm)*0.4D0)), 2)
  If (igap.Gt.2 .And. igap+kept.Gt.nd .And. info%crat.Gt.0.0D0) Then
     ss = trl_min_gap_ratio(info, nd, tind, res)
     If (ss .Gt. info%crat) igap = Max(2, nd-kept-1)
  End If
  If (lohi .Gt. 0) Then
     ss = merit(kr, kl)
     Do i = klm, krm-igap
        Do j = i+igap, krm
           tmp = merit(j, i)
           If (tmp.Gt.ss) Then
              ss = tmp
              kl = i
              kr = j
           End If
        End Do
     End Do
  Else
     ss = merit(kl, kr)
     Do i = klm, krm-igap
        Do j = i+igap, krm
           tmp = merit(i,j)
           If (tmp.Gt.ss) Then
              ss = tmp
              kl = i
              kr = j
           End If
        End Do
     End Do
  End If
End Subroutine trl_restart_max_progress
!!!
! search for a pair (kl, kr) such that the reduction in residual norm
! of the target (info%trgt) will be the largest before next restart
! the merit function is [sqrt(gap ratio) * (m-k)]
!!!
Subroutine trl_restart_max_reduction(nd, tind, kept, lambda, res, info, kl, kr)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Integer, Intent(in) :: kept, nd, tind
  Integer, Intent(inout) :: kl, kr
  Real(8), Intent(in) :: lambda(nd), res(nd)
  Interface
     Subroutine trl_restart_search_range(nd, lambda, res, info, ncl,&
          & ncr, lohi, tind, klm, krm)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: nd, ncl, ncr, tind
       Integer, Intent(out) :: klm, krm, lohi
       Real(8), Intent(in) :: lambda(nd), res(nd)
     End Subroutine trl_restart_search_range
     Function trl_min_gap_ratio(info, nd, tind, res) Result(gamma)
       Use trl_info
       Implicit None
       Type(TRL_INFO_T), Intent(in) :: info
       Integer, Intent(in) :: nd, tind
       Real(8), Intent(in) :: res(nd)
       Real(8) :: gamma
     End Function trl_min_gap_ratio
  End Interface
  !
  ! local variables
  Integer :: i, j, lohi, klm, krm, igap
  Real(8) :: tmp, ss, merit
  ! merit measure the factor of residual norm reduction
  merit(i,j) = Sqrt((lambda(i)-info%trgt)/(lambda(j)-info%trgt)) * Abs(j-i)
  ! determine the search range
  Call trl_restart_search_range(nd, lambda, res, info, kl, kr, lohi,&
       & tind, klm, krm)
  ! perform the brute-force search
  kl = klm
  kr = krm
  igap = Max(Min(nd-info%ned, Nint((krm-klm)*0.4D0)), 2)
  If (igap.Gt.2 .And. igap+kept.Gt.nd .And. info%crat.Gt.0.0D0) Then
     ss = trl_min_gap_ratio(info, nd, tind, res)
     If (ss .Gt. info%crat) igap = Max(2, nd-kept-1)
  End If
  If (lohi .Gt. 0) Then
     ss = merit(kr, kl)
     Do i = klm, krm-igap
        Do j = i+igap, krm
           tmp = merit(j, i)
           If (tmp.Gt.ss) Then
              ss = tmp
              kl = i
              kr = j
           End If
        End Do
     End Do
  Else
     ss = merit(kl, kr)
     Do i = klm, krm-igap
        Do j = i+igap, krm
           tmp = merit(i,j)
           If (tmp.Gt.ss) Then
              ss = tmp
              kl = i
              kr = j
           End If
        End Do
     End Do
  End If
End Subroutine trl_restart_max_reduction
!!!
! determine the search range -- used by the schemes that performs brute
! -force search
!!!
Subroutine trl_restart_search_range(nd, lambda, res, info, ncl, ncr,&
     & lohi, tind, klm, krm)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Integer, Intent(in) :: nd, ncl, ncr, tind
  Integer, Intent(out) :: klm, krm, lohi
  Real(8), Intent(in) :: lambda(nd), res(nd)
  !
  ! local variables
  Integer :: j
  Real(8) :: bnd
  klm = Max(ncl,1)
  krm = Min(ncr,nd)
  bnd = info%tol * info%anrm
  lohi = info%lohi
  ! make sure there is some distance between the boundary and the
  ! target Ritz value
  If (info%lohi .Gt. 0) Then
     ! save high end
     krm = Min(Max(info%maxlan-info%ned, (info%maxlan+info%nec)/2),&
          & krm, tind-1)
     Do While (krm+krm.Ge.ncl+ncr .And. res(krm).Le.bnd)
        krm = krm - 1
     End Do
  Else If (info%lohi .Lt. 0) Then
     ! save lower end
     klm = Max(Min(info%ned, (info%maxlan+info%nec)/2), tind+1, klm)
     Do While (klm+klm.Le.ncl+ncr .And. res(klm).Le.bnd)
        klm = klm + 1
     End Do
  Else
     ! save both ends
     If (tind-klm .Lt. krm-tind) Then
        lohi = -1
        klm = tind + 1
     Else
        lohi = 1
        krm = tind - 1
     End If
     j = info%locked + klm + nd - krm + 1
     If (j .Gt. 0) Then
        j = j / 2
        klm = klm - j
        krm = krm + j
     End If
  End If
End Subroutine trl_restart_search_range
!!!
! try to determine the minimal gap ratio need to compute all wanted
! eigenvalues 
!!!
Function trl_min_gap_ratio(info, nd, tind, res) Result(gamma)
  Use trl_info
  Implicit None
  Type(TRL_INFO_T), Intent(in) :: info
  Integer, Intent(in) :: nd, tind
  Real(8), Intent(in) :: res(nd)
  Real(8) :: gamma
  !
  gamma = info%maxmv * (info%nec + 1.0D0) / info%ned - info%matvec
  If (gamma .Lt. info%klan) Then
     gamma = Max((info%maxmv-info%matvec)/Dble(info%ned-info%nec), 2.0D0)
  End If
  gamma = Min(Log(res(tind)/(info%tol*info%anrm))/gamma, 0.5D0)
  Return
End Function trl_min_gap_ratio
