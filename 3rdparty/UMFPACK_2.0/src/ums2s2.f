 
        subroutine ums2s2 (n, job, transc, luxsiz, lux,
     $          luisiz, lui, b, x, r, z, ly, y, s, cntl, icntl, info,
     $          rinfo, cperm, rperm, presrv, an, anz, ap, ai, ax, on,
     $          nzoff, offp, offi, offx, nblks, lublkp, blkp, irstep)
c
cc UMS2S2 is a utility routine which solves a factored system.
c
        integer n, job, luxsiz, luisiz, lui (luisiz), ly, irstep,
     $          icntl (20), info (40), cperm (n), rperm (n), an,
     $          anz, ap (an+1), ai (anz), on, nzoff, offp (on+1),
     $          offi (nzoff), nblks, lublkp (nblks), blkp (nblks+1)
        logical transc, presrv
        real
     $          lux (luxsiz), b (n), x (n), r (n), z (n), y (ly),
     $          s (ly), cntl (10), rinfo (20), ax (anz),
     $          offx (nzoff)
 
c=== ums2s2 ============================================================
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
c  not user-callable.
 
c=======================================================================
c  description:
c=======================================================================
c
c  solve a system, given lu factors, permutation arrays, original
c  matrix (if preserved), and off-diagonal blocks (if btf was used).
 
c=======================================================================
c  input:
c=======================================================================
c
c       n:              order of matrix
c       job:            0: solve ax=b, 1: solve lx=b, 2: solve ux=b
c       transc:         if true, solve with transposed factors instead
c       luxsiz:         size of lux
c       lux (1..luxsiz) real values in lu factors for each block
c       luisiz:         size of lui
c       lui (1..luisiz) integers in lu factors for each block
c       b (1..n):       right-hand-side
c       ly:             size of y and s, ly=n if y and s are used
c       cntl:           real control parameters, see ums2in
c       icntl:          integer control parameters, see ums2in
c       cperm (1..n):   q, column permutation array
c       rperm (1..n):   p, row permutation array
c       presrv:         if true, then original matrix was preserved
c       nblks:          number of diagonoal blocks (1 if no btf)
c       irstep:         maximum number of steps of iterative refinement
c
c       if presrv then
c           an:                 order of preserved matrix, n
c           anz:                number of entries in preserved matrix
c           ap (1..an+1):       column pointers of preserved matrix
c           ai (1..anz):        row indices of preserved matrix
c           ax (1..anz):        values of preserved matrix
c           an, anz, ap, ai, ax:        not accessed
c
c       if nblks > 1 then
c           on:                 n
c           nzoff:              number of off-diagonoal entries
c           offp (1..n+1)       row pointers for off-diagonal part
c           offi (1..nzoff):    column indices for off-diagonal part
c           offx (1..nzoff):    values of off-diagonal part
c           lublkp (1..nblks):  pointers to lu factors of each block
c           blkp (1..nblks+1):  index range of each block
c       else
c           on, nzoff, offp, offi, offx, lublkp, blkp:  not accessed
 
c=======================================================================
c  workspace:
c=======================================================================
c
c       r (1..n), z (1..n)
c       y (1..ly), s (1..ly):   unaccessed if no iterative refinement
 
c=======================================================================
c  output:
c=======================================================================
c
c       x (1..n):       solution
c       info:           integer informational output, see ums2in
c       rinfo:          real informational output, see ums2in
c
c       if irsteps > 0 and presrv is true then
c           w (1..n):           residual
c           rinfo (7):  sparse error estimate, omega1
c           rinfo (8):  sparse error estimate, omega2
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   ums2so
c       subroutines called:     ums2er, ums2sl, ums2lt, ums2su, ums2ut
c       functions called:       isamax, abs, max
        intrinsic abs, max
        integer isamax
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer nlu, i, blk, k1, k2, kn, p, step, npiv, j
        real
     $          zero, one, xnorm, tau, nctau, omega1, omega2, d1,
     $          d2, omega, omlast, om1lst, om2lst, two, eps, maxeps,
     $          thosnd, a, axx, r2, x2, y2, z2
        parameter (zero = 0.0, one = 1.0, two = 2.0,
     $          maxeps = two ** (-15), thosnd = 1000.0)
 
c  lu factors:
c  -----------
c  blk:     current diagonal block
c  k1,k2:   current diagonal block is a (k1..k2, k1..k2)
c  kn:      size of diagonal block (= k2-k1+1)
c  nlu:     number of elements in the lu factors of a single diag block
c  npiv:    number of pivots in the lu factors of a single diag block
c
c  iterative refinement and sparse backward error:
c  -----------------------------------------------
c  step:    number of steps of iterative refinement taken
c  xnorm:   ||x|| maxnorm of solution vector, x
c  tau:     threshold for selecting which estimate to use (1 or 2)
c  nctau:   1000*n*eps
c  eps:     largest positive value such that fl (1.0 + eps) = 1.0
c  omega1:  current sparse backward error estimate 1
c  omega2:  current sparse backward error estimate 2
c  d1:      divisor for omega1
c  d2:      divisor for omega2
c  omega:   omega1 + omega2
c  omlast:  value of omega from previous step
c  om1lst:  value of omega1 from previous step
c  om2lst:  value of omega2 from previous step
c  maxeps:  2**(-16), maximum value that eps is allowed to be
c  a:       value of an entry in a, a_ij
c  axx:     a_ij * x_j
c
c  other:
c  ------
c  i,j:     loop indices
c  p:       pointer
c  r2:      r (i)
c  x2:      x (i)
c  y2:      y (i)
c  z2:      z (i)
 
c=======================================================================
c  executable statements:
c=======================================================================
 
c-----------------------------------------------------------------------
c  initializations for sparse backward error
c-----------------------------------------------------------------------
 
        omega = zero
        omega1 = zero
        omega2 = zero
        eps = cntl (3)
        if (eps .le. zero .or. eps .gt. maxeps) then
c          eps is too small or too big: set to a large default value
           eps = maxeps
        endif
        nctau = thosnd * n * eps
 
c-----------------------------------------------------------------------
c  get information on lu factorization if btf was not used
c-----------------------------------------------------------------------
 
        if (nblks .eq. 1) then
c          p is 1, and lui (p) is 1
           nlu = lui (2)
           npiv = lui (3)
        endif
 
c-----------------------------------------------------------------------
        if (job .eq. 1) then
c-----------------------------------------------------------------------
 
c          -------------------------------------------------------------
           if (.not. transc) then
c          -------------------------------------------------------------
 
c             ----------------------------------------------------------
c             solve p'lx=b:  x = l \ pb
c             ----------------------------------------------------------
 
              do i = 1, n
                 x (i) = b (rperm (i))
              end do
              if (nblks .eq. 1) then
                 call ums2sl (nlu, npiv, n, lui(6), lui(6+nlu), lux,x,z)
              else
                 do 20 blk = 1, nblks
                    k1 = blkp (blk)
                    k2 = blkp (blk+1) - 1
                    kn = k2-k1+1
                    if (kn .gt. 1) then
                       p = lublkp (blk)
                       nlu = lui (p+1)
                       npiv = lui (p+2)
                       call ums2sl (nlu, npiv, kn, lui (p+5),
     $                    lui (p+5+nlu), lux (lui (p)), x (k1), z)
                    endif
20               continue
              endif
 
c          -------------------------------------------------------------
           else
c          -------------------------------------------------------------
 
c             ----------------------------------------------------------
c             solve l'px=b:  x = p' (l' \ b)
c             ----------------------------------------------------------
 
              do 30 i = 1, n
                 r (i) = b (i)
30            continue
              if (nblks .eq. 1) then
                 call ums2lt (nlu, npiv, n, lui(6), lui(6+nlu), lux,r,z)
              else
                 do 40 blk = 1, nblks
                    k1 = blkp (blk)
                    k2 = blkp (blk+1) - 1
                    kn = k2-k1+1
                    if (kn .gt. 1) then
                       p = lublkp (blk)
                       nlu = lui (p+1)
                       npiv = lui (p+2)
                       call ums2lt (nlu, npiv, kn, lui (p+5),
     $                    lui (p+5+nlu), lux (lui (p)), r (k1), z)
                    endif
40               continue
              endif
              do 50 i = 1, n
                 x (rperm (i)) = r (i)
50            continue
 
c          -------------------------------------------------------------
           endif
c          -------------------------------------------------------------
 
c-----------------------------------------------------------------------
        else if (job .eq. 2) then
c-----------------------------------------------------------------------
 
c          -------------------------------------------------------------
           if (transc) then
c          -------------------------------------------------------------
 
c             ----------------------------------------------------------
c             solve qu'x=b:  x = u' \ q'b
c             ----------------------------------------------------------
 
              do 60 i = 1, n
                 x (i) = b (cperm (i))
60            continue
              if (nblks .eq. 1) then
                 call ums2ut (nlu, npiv, n, lui(6), lui(6+nlu), lux,x,z)
              else
                 do 100 blk = 1, nblks
                    k1 = blkp (blk)
                    k2 = blkp (blk+1) - 1
                    kn = k2-k1+1
                    if (kn .eq. 1) then
                       x (k1) = x (k1) / lux (lublkp (blk))
                       r (k1) = x (k1)
                    else
                       p = lublkp (blk)
                       nlu = lui (p+1)
                       npiv = lui (p+2)
                       call ums2ut (nlu, npiv, kn, lui (p+5),
     $                    lui (p+5+nlu), lux (lui (p)), x (k1), z)
                       do 70 i = k1, k2
                          r (i) = x (i)
70                     continue
                       call ums2lt (nlu, npiv, kn, lui (p+5),
     $                    lui (p+5+nlu), lux (lui (p)), r (k1), z)
                    endif
                    do 90 i = k1, k2
                       r2 = r (i)
                       do 80 p = offp (i), offp (i+1)-1
                          x (offi (p)) = x (offi (p)) - offx (p) * r2
80                     continue
90                  continue
100              continue
              endif
 
c          -------------------------------------------------------------
           else
c          -------------------------------------------------------------
 
c             ----------------------------------------------------------
c             solve uq'x=b:  x = q (u \ b)
c             ----------------------------------------------------------
 
              if (nblks .eq. 1) then
                 do 110 i = 1, n
                    r (i) = b (i)
110              continue
                 call ums2su (nlu, npiv, n, lui(6), lui(6+nlu), lux,r,z)
              else
                 do 150 blk = nblks, 1, -1
                    k1 = blkp (blk)
                    k2 = blkp (blk+1) - 1
                    kn = k2-k1+1
                    do 130 i = k1, k2
                       x2 = zero
                       do 120 p = offp (i), offp (i+1)-1
                          x2 = x2 + offx (p) * r (offi (p))
120                    continue
                       x (i) = x2
130                 continue
                    if (kn .eq. 1) then
                       r (k1) = (b (k1) - x (k1)) / lux (lublkp (blk))
                    else
                       p = lublkp (blk)
                       nlu = lui (p+1)
                       npiv = lui (p+2)
                       call ums2sl (nlu, npiv, kn, lui (p+5),
     $                    lui (p+5+nlu), lux (lui (p)), x (k1), z)
                       do 140 i = k1, k2
                          r (i) = b (i) - x (i)
140                    continue
                       call ums2su (nlu, npiv, kn, lui (p+5),
     $                    lui (p+5+nlu), lux (lui (p)), r (k1), z)
                    endif
150              continue
              endif
              do 160 i = 1, n
                 x (cperm (i)) = r (i)
160           continue
 
c          -------------------------------------------------------------
           endif
c          -------------------------------------------------------------
 
c-----------------------------------------------------------------------
        else
c-----------------------------------------------------------------------
 
           do 450 step = 0, irstep
 
c             ----------------------------------------------------------
c             if transa was true in ums2fa or ums2rf, then c = a'.
c             otherwise c = a.  in both cases, the factorization is
c             pcq = lu, and c is stored in column-form in ai,ax,ap if
c             it is preserved.
c             ----------------------------------------------------------
 
c             ----------------------------------------------------------
              if (.not. transc) then
c             ----------------------------------------------------------
 
c                -------------------------------------------------------
c                solve cx=b (step 0):
c                   x = q (u \ l \ pb)
c                and then perform iterative refinement (step > 0):
c                   x = x + q (u \ l \ p (b-cx))
c                -------------------------------------------------------
 
                 if (step .eq. 0) then
                    do 170 i = 1, n
                       r (i) = b (rperm (i))
170                 continue
                 else
                    do 180 i = 1, n
                       z (i) = b (i)
180                 continue
                    do 200 i = 1, n
                       x2 = x (i)
                       do 190 p = ap (i), ap (i+1) - 1
                          z (ai (p)) = z (ai (p)) - ax (p) * x2
190                    continue
200                 continue
                    do 210 i = 1, n
                       r (i) = z (rperm (i))
210                 continue
                 endif
                 if (nblks .eq. 1) then
                    call ums2sl (nlu, npiv, n,lui(6),lui(6+nlu),lux,r,z)
                    call ums2su (nlu, npiv, n,lui(6),lui(6+nlu),lux,r,z)
                 else
                    do 240 blk = nblks, 1, -1
                       k1 = blkp (blk)
                       k2 = blkp (blk+1) - 1
                       kn = k2-k1+1
                       do 230 i = k1, k2
                          r2 = r (i)
                          do 220 p = offp (i), offp (i+1)-1
                             r2 = r2 - offx (p) * r (offi (p))
220                       continue
                          r (i) = r2
230                    continue
                       if (kn .eq. 1) then
                          r (k1) = r (k1) / lux (lublkp (blk))
                       else
                          p = lublkp (blk)
                          nlu = lui (p+1)
                          npiv = lui (p+2)
                          call ums2sl (nlu, npiv, kn, lui (p+5),
     $                       lui (p+5+nlu), lux (lui (p)), r (k1), z)
                          call ums2su (nlu, npiv, kn, lui (p+5),
     $                       lui (p+5+nlu), lux (lui (p)), r (k1), z)
                       endif
240                 continue
                 endif
                 if (step .eq. 0) then
                    do 250 i = 1, n
                       x (cperm (i)) = r (i)
250                 continue
                 else
                    do 260 i = 1, n
                       x (cperm (i)) = x (cperm (i)) + r (i)
260                 continue
                 endif
 
c             ----------------------------------------------------------
              else
c             ----------------------------------------------------------
 
c                -------------------------------------------------------
c                solve c'x=b (step 0):
c                   x = p' (l' \ u' \ q'b)
c                and then perform iterative refinement (step > 0):
c                   x = x + p' (l' \ u' \ q' (b-c'x))
c                -------------------------------------------------------
 
                 if (step .eq. 0) then
                    do 270 i = 1, n
                       r (i) = b (cperm (i))
270                 continue
                 else
                    do 280 i = 1, n
                       z (i) = b (i)
280                 continue
                    do 300 i = 1, n
                       z2 = z (i)
                       do 290 p = ap (i), ap (i+1) - 1
                          z2 = z2 - ax (p) * x (ai (p))
290                    continue
                       z (i) = z2
300                 continue
                    do 310 i = 1, n
                       r (i) = z (cperm (i))
310                 continue
                 endif
                 if (nblks .eq. 1) then
                    call ums2ut (nlu, npiv, n,lui(6),lui(6+nlu),lux,r,z)
                    call ums2lt (nlu, npiv, n,lui(6),lui(6+nlu),lux,r,z)
                 else
                    do 340 blk = 1, nblks
                       k1 = blkp (blk)
                       k2 = blkp (blk+1) - 1
                       kn = k2-k1+1
                       if (kn .eq. 1) then
                          r (k1) = r (k1) / lux (lublkp (blk))
                       else
                          p = lublkp (blk)
                          nlu = lui (p+1)
                          npiv = lui (p+2)
                          call ums2ut (nlu, npiv, kn, lui (p+5),
     $                       lui (p+5+nlu), lux (lui (p)), r (k1), z)
                          call ums2lt (nlu, npiv, kn, lui (p+5),
     $                       lui (p+5+nlu), lux (lui (p)), r (k1), z)
                       endif
                       do 330 i = k1, k2
                          r2 = r (i)
                          do 320 p = offp (i), offp (i+1)-1
                             r (offi (p)) = r (offi (p)) - offx (p) * r2
320                       continue
330                    continue
340                 continue
                 endif
                 if (step .eq. 0) then
                    do 350 i = 1, n
                       x (rperm (i)) = r (i)
350                 continue
                 else
                    do 360 i = 1, n
                       x (rperm (i)) = x (rperm (i)) + r (i)
360                 continue
                 endif
 
c             ----------------------------------------------------------
              endif
c             ----------------------------------------------------------
 
c             ----------------------------------------------------------
c             sparse backward error estimate
c             ----------------------------------------------------------
 
              if (irstep .gt. 0) then
 
c                xnorm = ||x|| maxnorm
                 xnorm = abs (x (isamax (n, x, 1)))
 
c                r (i) = (b-ax)_i, residual (or a')
c                z (i) = (|a||x|)_i
c                y (i) = ||a_i||, maxnorm of row i of a (or a')
                 do 370 i = 1, n
                    r (i) = b (i)
                    z (i) = zero
                    y (i) = zero
370              continue
 
                 if (.not. transc) then
 
c                   ----------------------------------------------------
c                   sparse backward error for cx=b, c stored by column
c                   ----------------------------------------------------
 
                    do 390 j = 1, n
                       x2 = x (j)
cfpp$ nodepchk l
                       do 380 p = ap (j), ap (j+1) - 1
                          i = ai (p)
                          a = ax (p)
                          axx = a * x2
                          r (i) = r (i) -     (axx)
                          z (i) = z (i) + abs (axx)
                          y (i) = y (i) + abs (a)
380                    continue
390                 continue
 
                 else
 
c                   ----------------------------------------------------
c                   sparse backward error for c'x=b, c' stored by row
c                   ----------------------------------------------------
 
                    do 410 i = 1, n
                       r2 = r (i)
                       z2 = z (i)
                       y2 = y (i)
cfpp$ nodepchk l
                       do 400 p = ap (i), ap (i+1) - 1
                          j = ai (p)
                          a = ax (p)
                          axx = a * x (j)
                          r2 = r2 -     (axx)
                          z2 = z2 + abs (axx)
                          y2 = y2 + abs (a)
400                    continue
                       r (i) = r2
                       z (i) = z2
                       y (i) = y2
410                 continue
 
                 endif
 
c                -------------------------------------------------------
c                save the last iteration in case we need to reinstate it
c                -------------------------------------------------------
 
                 omlast = omega
                 om1lst = omega1
                 om2lst = omega2
 
c                -------------------------------------------------------
c                compute sparse backward errors: omega1 and omega2
c                -------------------------------------------------------
 
                 omega1 = zero
                 omega2 = zero
                 do 420 i = 1, n
                    tau = (y (i) * xnorm + abs (b (i))) * nctau
                    d1 = z (i) + abs (b (i))
                    if (d1 .gt. tau) then
                       omega1 = max (omega1, abs (r (i)) / d1)
                    else if (tau .gt. zero) then
                       d2 = z (i) + y (i) * xnorm
                       omega2 = max (omega2, abs (r (i)) / d2)
                    endif
420              continue
                 omega = omega1 + omega2
                 rinfo (7) = omega1
                 rinfo (8) = omega2
 
c                -------------------------------------------------------
c                stop the iterations if the backward error is small
c                -------------------------------------------------------
 
                 info (24) = step
                 if (one + omega .le. one) then
c                   further iterative refinement will no longer improve
c                   the solution
                    return
                 endif
 
c                -------------------------------------------------------
c                stop if insufficient decrease in omega
c                -------------------------------------------------------
 
                 if (step .gt. 0 .and. omega .gt. omlast / two) then
                    if (omega .gt. omlast) then
c                      last iteration better than this one, reinstate it
                       do 430 i = 1, n
                          x (i) = s (i)
                          rinfo (7) = om1lst
                          rinfo (8) = om2lst
430                    continue
                    endif
                    info (24) = step - 1
                    return
                 endif
 
c                -------------------------------------------------------
c                save current solution in case we need to reinstate
c                -------------------------------------------------------
 
                 do 440 i = 1, n
                    s (i) = x (i)
440              continue
 
              endif
 
450        continue
 
c-----------------------------------------------------------------------
        endif
c-----------------------------------------------------------------------
 
        return
        end
