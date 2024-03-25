 
        subroutine ums2so (n, job, transc, lvalue, lindex, value,
     $          index, keep, b, x, w, cntl, icntl, info, rinfo)
c
cc UMS2SO solves a linear system associated with a factored sparse matrix.
c
        integer n, job, lvalue, lindex, index (lindex), keep (20),
     $          icntl (20), info (40)
        real
     $          value (lvalue), b (n), x (n), w (*), cntl (10),
     $          rinfo (20)
        logical transc
 
c=== ums2so ============================================================
c
c  unsymmetric-pattern multifrontal package (umfpack). version 2.0s.
c  copyright (c) 1995, timothy a. davis, university of florida, usa.
c  joint work with iain s. duff, rutherford appleton laboratory, uk.
c  september 1995. work supported by the national science foundation
c  (dms-9223088 and dms-9504974) and the state of florida; and by cray
c  research inc. through the allocation of supercomputing resources.
 
c***********************************************************************
c* notice:  "the umfpack package may be used solely for educational,   *
c* research, and benchmarking purposes by non-profit organizations and *
c* the u.s. government.  commericial and other organizations may make  *
c* use of umfpack solely for benchmarking purposes only.  umfpack may  *
c* be modified by or on behalf of the user for such use but at no time *
c* shall umfpack or any such modified version of umfpack become the    *
c* property of the user.  umfpack is provided without warranty of any  *
c* kind, either expressed or implied.  neither the authors nor their   *
c* employers shall be liable for any direct or consequential loss or   *
c* damage whatsoever arising out of the use or misuse of umfpack by    *
c* the user.  umfpack must not be sold.  you may make copies of        *
c* umfpack, but this notice and the copyright notice must appear in    *
c* all copies.  any other use of umfpack requires written permission.  *
c* your use of umfpack is an implicit agreement to these conditions."  *
c*                                                                     *
c* the ma38 package in the harwell subroutine library (hsl) has        *
c* equivalent functionality (and identical calling interface) as       *
c* umfpack.  it is available for commercial use.   technical reports,  *
c* information on hsl, and matrices are available via the world wide   *
c* web at http://www.cis.rl.ac.uk/struct/arcd/num.html, or by          *
c* anonymous ftp at seamus.cc.rl.ac.uk/pub.  also contact john         *
c* harding, harwell subroutine library, b 552, aea technology,         *
c* harwell, didcot, oxon ox11 0ra, england.                            *
c* telephone (44) 1235 434573, fax (44) 1235 434340,                   *
c* email john.harding@aeat.co.uk, who will provide details of price    *
c* and conditions of use.                                              *
c***********************************************************************
 
c=======================================================================
c  user-callable.
 
c=======================================================================
c  description:
c=======================================================================
c
c  given lu factors computed by ums2fa or ums2rf, and the
c  right-hand-side, b, solve a linear system for the solution x.
c
c  this routine handles all permutations, so that b and x are in terms
c  of the original order of the matrix, a, and not in terms of the
c  permuted matrix.
c
c  if iterative refinement is done, then the residual is returned in w,
c  and the sparse backward error estimates are returned in
c  rinfo (7) and rinfo (8).  the computed solution x is the
c  exact solution of the equation (a + da)x = (b + db), where
c    da (i,j)  <= max (rinfo (7), rinfo (8)) * abs (a(i,j))
c  and
c    db (i) <= max (rinfo (7) * abs (b (i)),
c                   rinfo (8) * maxnorm (a) * maxnorm (x computed))
c  note that da has the same sparsity pattern as a.
c  the method used to compute the sparse backward error estimate is
c  described in m. arioli, j. w. demmel, and i. s. duff, "solving
c  sparse linear systems with sparse backward error," siam j. matrix
c  analysis and applications, vol 10, 1989, pp. 165-190.
 
c=======================================================================
c  arguments:
c=======================================================================
 
c           ------------------------------------------------------------
c  n:       an integer variable.
c           must be set by caller on input (not modified).
c           must be the same as passed to ums2fa or ums2rf.
 
c           ------------------------------------------------------------
c  job:     an integer variable.
c           must be set by caller on input (not modified).
c           what system to solve (see the transc argument below).
c           iterative refinement is only performed if job = 0,
c           icntl (8) > 0, and only if the original matrix was
c           preserved (job = 1 in ums2fa or ums2rf).
 
c           ------------------------------------------------------------
c  transc:  an integer variable.
c           must be set by caller on input (not modified).
c           solve with l and u factors or with l' and u', where
c           transa was passed to ums2fa or ums2rf.
c
c           if transa = false, then paq = lu was performed,
c           and the following systems are solved:
c
c                               transc = false          transc = true
c                               ----------------        ----------------
c                  job = 0      solve ax = b            solve a'x = b
c                  job = 1      solve p'lx = b          solve l'px = b
c                  job = 2      solve uq'x = b          solve qu'x = b
c
c           if transa = true, then a was transformed prior to lu
c           factorization, and p(a')q = lu
c
c                               transc = false          transc = true
c                               ----------------        ----------------
c                  job = 0      solve a'x = b           solve ax = b
c                  job = 1      solve p'lx = b          solve l'px = b
c                  job = 2      solve uq'x = b          solve qu'x = b
c
c           other values of job are treated as zero.  iterative
c           refinement can be done only when solving ax=b or a'x=b.
c
c           the comments below use matlab notation, where
c           x = l \ b means x = (l^(-1)) * b, premultiplication by
c           the inverse of l.
 
c           ------------------------------------------------------------
c  lvalue:  an integer variable.
c           must be set by caller on input (not modified).
c           the size of value.
 
c           ------------------------------------------------------------
c  lindex:  an integer variable.
c           must be set by caller on input (not modified).
c           the size of index.
 
c           ------------------------------------------------------------
c  value:   a real array of size lvalue.
c           must be set by caller on input (normally from last call to
c           ums2fa or ums2rf) (not modified).
c           the lu factors, in value (keep (1) ... keep (2)).
c           the entries in value (1 ... keep (1) - 1) and in
c           value (keep (2) + 1 ... lvalue) are not accessed.
 
c           ------------------------------------------------------------
c  index:   an integer array of size lindex.
c           must be set by caller on input (normally from last call to
c           ums2fa or ums2rf) (not modified).
c           the lu factors, in index (keep (3) ... keep (5)).
c           the entries in index (1 ... keep (3) - 1) and in
c           index (keep (5) + 1 ... lindex) are not accessed.
 
c           ------------------------------------------------------------
c  keep:    an integer array of size 20.
c
c           keep (1..5): must be set by caller on input (normally from
c               last call to ums2fa or ums2rf) (not modified).
c               layout of the lu factors in value and index
 
c           ------------------------------------------------------------
c  b:       a real array of size n.
c           must be set by caller on input (not modified).
c           the right hand side, b, of the system to solve.
 
c           ------------------------------------------------------------
c  w:       a real array of size 2*n or 4*n.
c           need not be set by caller on input.  modified on output.
c           workspace of size w (1..2*n) if icntl (8) = 0, which
c           is the default value.  if iterative refinement is
c           performed, and w must be of size w (1..4*n) and the
c           residual b-ax (or b-a'x) is returned in w (1..n).
 
c           ------------------------------------------------------------
c  x:       a real array of size n.
c           need not be set by caller on input.  modified on output.
c           the solution, x, of the system that was solved.  valid only
c           if info (1) is greater than or equal to 0.
 
c           ------------------------------------------------------------
c  cntl:    a real array of size 10.
c           must be set by caller on input (not modified).
c           real control parameters, see ums2in for a description,
c           which sets the defaults.
 
c           ------------------------------------------------------------
c  icntl:   an integer array of size 20.
c           must be set by caller on input (not modified).
c           integer control parameters, see ums2in for a description,
c           which sets the defaults.  in particular, icntl (8) is
c           the maximum number of steps of iterative refinement to be
c           performed.
 
c           ------------------------------------------------------------
c  info:    an integer array of size 40.
c           need not be set by caller on input.  modified on output.
c           it contains information about the execution of ums2so.
c
c           info (1) is the error flag.  if info (1) is -7, then
c           the lu factors are uncomputed, or have been corrupted since
c           the last call to ums2fa or ums2rf.  no system is solved,
c           and x (1..n) is not valid on output.  if info (1) is 8,
c           then iterative refinement was requested but cannot be done.
c           to perform iterative refinement, the original matrix must be
c           preserved (job = 1 in ums2fa or ums2rf) and ax=b or a'x=b
c           must be solved (job = 0 in ums2so).  info (24) is the
c           steps of iterative refinement actually taken.
 
c           ------------------------------------------------------------
c  rinfo:   a real array of size 20.
c           need not be set by caller on input.  modified on output.
c           it contains information about the execution of ums2so.
c
c           if iterative refinement was performed then
c           rinfo (7) is the sparse error estimate, omega1, and
c           rinfo (8) is the sparse error estimate, omega2.
 
c=======================================================================
c  to be preserved between calls to ums2fa, ums2rf, ums2so:
c=======================================================================
c
c  the following must be unchanged since the call to ums2fa or ums2rf
c  that computed the lu factors:
c
c       n
c       value (keep (1) ... keep (2))
c       index (keep (3) ... keep (5))
c       keep (1 ... 20)
 
c## end of user documentation ##########################################
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   user routine
c       subroutines called:     ums2er, ums2p1, ums2s2
c       functions called:       max
        intrinsic max
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer nblks, offip, offxp, n1, nz, ne, offpp, blkpp, lublpp,
     $          app, an, anz, on, lui1, lui2, lux1, lux2, aip, axp,
     $          cpermp, rpermp, nzoff, irstep, yp, ly, lw, sp, ip1, ip2,
     $          xp1, luir1, io, prl
        logical presrv, badlu
        real
     $          zero
        parameter (zero = 0.0)
 
c  printing control:
c  -----------------
c  io:      i/o unit for diagnostic messages
c  prl:     printing level
c
c  location and status of lu factors:
c  ----------------------------------
c  lui1:    integer part of lu factors start in index (lui1...)
c  luir1:   index (luir1 ... lui2) is needed for a call to ums2rf
c  lui2:    integer part of lu factors end in index (..lui2)
c  lux1:    real part of lu factors start in value (lux1...)
c  lux2:    real part of lu factors end in value (...lux1)
c  ip1:     pointer into leading part of lu factors in index
c  ip2:     pointer into trailing part of lu factors in index
c  xp1:     pointer into leading part of lu factors in value
c  badlu:   if true, then lu factors are corrupted or not computed
c
c  arrays and scalars allocated in lu factors (in order):
c  ------------------------------------------------------
c  app:     ap (1..n+1) array located in index (app...app+n)
c  axp:     ax (1..nz) array located in value (axp...axp+nz-1)
c  aip:     ai (1..nz) array located in index (aip...aip+nz-1)
c  an:      n if a is preserved, 1 otherwise
c  anz:     nz if a is preserved, 1 otherwise
c  offip:   offi (1..nzoff) array loc. in index (offip...offip+nzoff-1)
c  offxp:   offx (1..nzoff) array loc. in value (offxp...offxp+nzoff-1)
c  ...      lu factors of each diagonal block located here
c  lublpp:  lublkp (1..nblks) array in index (lublpp..lublpp+nblks-1)
c  blkpp:   blkp (1..nblks+1) array loc. in index (blkpp...blkpp+nblks)
c  offpp:   offp (1..n+1) array located in index (offpp...offpp+n)
c  on:      size of offp (1..n+1):  n if nblks > 1, 1 otherwise
c  cpermp:  cperm (1..n) array located in index (cpermp...cpermp+n-1)
c  rpermp:  rperm (1..n) array located in index (rpermp...rpermp+n-1)
c  ...      seven scalars in index (lui2-6...lui2):
c  nzoff:   number of entries in off-diagonal part
c  nblks:   number of diagonal blocks
c  presrv:  true if original matrix was preserved when factorized
c  nz:      entries in a
c  n1:      n argument in ums2fa or ums2rf when matrix factorized
c  ne:      ne argument in ums2fa or ums2rf when matrix factorized
c
c  arrays allocated from w work array:
c  -----------------------------------
c  lw:      size of w
c  yp:      y (1..n) located in w (yp...yp+n-1)
c  sp:      s (1..n) located in w (sp...sp+n-1)
c  ly:      size of y and s
c
c  other:
c  ------
c  irstep:  maximum number of iterative refinement steps to take
 
c=======================================================================
c  executable statements:
c=======================================================================
 
        io = icntl (2)
        prl = icntl (3)
 
c-----------------------------------------------------------------------
c  clear informational output
c-----------------------------------------------------------------------
 
        info (1) = 0
        info (24) = 0
        rinfo (7) = zero
        rinfo (8) = zero
 
c-----------------------------------------------------------------------
c  print input arguments if requested
c-----------------------------------------------------------------------
 
        irstep = max (0, icntl (8))
        if (irstep .eq. 0) then
           lw = 2*n
        else
           lw = 4*n
        endif
        call ums2p1 (3, 1,
     $          n, ne, job, transc, lvalue, lindex, value,
     $          index, keep, cntl, icntl, info, rinfo,
     $          b, x, n, w, lw)
 
c-----------------------------------------------------------------------
c  get pointers to lu factors
c-----------------------------------------------------------------------
 
        lux1 = keep (1)
        lux2 = keep (2)
        lui1 = keep (3)
        luir1 = keep (4)
        lui2 = keep (5)
        badlu = luir1 .le. 0 .or. lui2-6 .lt. luir1
     $     .or. lui2 .gt. lindex
     $     .or. lux1 .le. 0 .or. lux1 .gt. lux2 .or. lux2 .gt. lvalue
     $     .or. lui1 .le. 0 .or. luir1 .lt. lui1 .or. luir1 .gt. lui2
        if (badlu) then
           call ums2er (3, icntl, info, -7, 0)
c          error return, lu factors are corrupted:
           go to 9000
        endif
 
c-----------------------------------------------------------------------
c  get seven scalars (transa, nzoff, nblks, presrv, nz, n, ne) from lu
c-----------------------------------------------------------------------
 
        ne = index (lui2)
        n1 = index (lui2-1)
        nz = index (lui2-2)
        presrv = index (lui2-3) .ne. 0
        nblks = index (lui2-4)
        nzoff = index (lui2-5)
c       transa = index (lui2-6) .ne. 0, we don't actually need this here
 
c-----------------------------------------------------------------------
c  get pointers to permutation vectors
c-----------------------------------------------------------------------
 
        rpermp = (lui2-6) - n
        cpermp = rpermp - n
        ip2 = cpermp - 1
        xp1 = lux1
        ip1 = lui1
 
c-----------------------------------------------------------------------
c  get pointers to preserved column-oriented copy of input matrix
c-----------------------------------------------------------------------
 
        if (presrv) then
 
c          -------------------------------------------------------------
c          original matrix preserved in index (lui1..lui1+nz+n) and
c          value (lux1..lux1+nz-1)
c          -------------------------------------------------------------
 
           app = ip1
           aip = app + n+1
           ip1 = aip + nz
           axp = xp1
           xp1 = axp + nz
           an = n
           anz = nz
 
        else
 
c          -------------------------------------------------------------
c          original matrix not preserved, pass dummy argument to ums2s2
c          -------------------------------------------------------------
 
           app = 1
           aip = 1
           axp = 1
           an = 1
           anz = 1
 
        endif
 
c-----------------------------------------------------------------------
c  get pointers to block-triangular information, if btf was used
c-----------------------------------------------------------------------
 
        if (nblks .gt. 1) then
 
c          -------------------------------------------------------------
c          get pointers to off-diagonal nonzeros, and btf arrays
c          -------------------------------------------------------------
 
           offip = ip1
           ip1 = ip1 + nzoff
           offxp = xp1
           xp1 = xp1 + nzoff
           offpp = cpermp - (n+1)
           blkpp = offpp - (nblks+1)
           lublpp = blkpp - (nblks)
           ip2 = lublpp - 1
           on = n
 
        else
 
c          -------------------------------------------------------------
c          matrix was factorized as a single block, pass dummy arg.
c          -------------------------------------------------------------
 
           offip = 1
           offxp = 1
           offpp = 1
           blkpp = 1
           lublpp = 1
           on = 1
 
        endif
 
        badlu = n .ne. n1 .or. nz .le. 0 .or. luir1 .gt. ip2 .or.
     $     nblks .le. 0 .or. nblks .gt. n .or.
     $     xp1 .gt. lux2 .or. nzoff .lt. 0 .or. ip1 .ne. luir1
        if (badlu) then
           call ums2er (3, icntl, info, -7, 0)
c          error return, lu factors are corrupted:
           go to 9000
        endif
 
c-----------------------------------------------------------------------
c  get the number of steps of iterative refinement
c-----------------------------------------------------------------------
 
        if (irstep .gt. 0 .and. .not. presrv) then
c          original matrix not preserved (ums2fa/ums2rf job .ne. 1)
           call ums2er (3, icntl, info, 8, 0)
           irstep = 0
        endif
        if (irstep .gt. 0 .and. (job .eq. 1 .or. job .eq. 2)) then
c          iterative refinement for ax=b and a'x=b only (job = 0)
           call ums2er (3, icntl, info, 8, 1)
           irstep = 0
        endif
        if (irstep .eq. 0) then
c          pass a dummy argument as y, which is not accessed in ums2s2
           yp = 1
           ly = 1
           sp = 1
           lw = 2*n
        else
c          pass w (yp ... yp+n-1) as y (1..n) to ums2s2
           yp = 2*n+1
           ly = n
           sp = 3*n+1
           lw = 4*n
        endif
 
c-----------------------------------------------------------------------
c  solve; optional iterative refinement and sparse backward error
c-----------------------------------------------------------------------
 
        call ums2s2 (n, job, transc, lux2-xp1+1, value (xp1),
     $     ip2-luir1+1, index (luir1), b, x,
     $     w, w (n+1), ly, w (yp), w (sp),
     $     cntl, icntl, info, rinfo, index (cpermp), index (rpermp),
     $     presrv, an, anz, index (app), index (aip), value (axp),
     $     on, max (1, nzoff), index (offpp), index (offip),
     $     value (offxp), nblks, index (lublpp), index (blkpp), irstep)
 
c-----------------------------------------------------------------------
c  print output arguments if requested
c-----------------------------------------------------------------------
 
c       error return label:
9000    continue
        call ums2p1 (3, 2,
     $          n, ne, job, transc, lvalue, lindex, value,
     $          index, keep, cntl, icntl, info, rinfo,
     $          b, x, n, w, lw)
        return
        end
