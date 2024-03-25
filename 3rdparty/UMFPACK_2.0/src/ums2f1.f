 
        subroutine ums2f1 (cp, n, cperm, rperm, nzoff,
     $          itail, xtail, xx, xsize, xuse, ii, isize, iuse,
     $          icntl, cntl, info, rinfo, nblks,
     $          ap, ai, ax, presrv, k1, an, anz, pr, keep,
     $          rmax, cmax, totnlu, xrmax, xruse, iout, xout)
c
cc UMS2F1 is a utility routine which factors part of the matrix.
c
        integer xsize, isize, n, icntl (20), info (40), xuse, iuse,
     $          itail, xtail, ii (isize), cp (n+1), cperm (n), nzoff,
     $          an, anz, rperm (n), ai (anz), ap (an+1), k1, pr (an),
     $          nblks, keep (20), rmax, cmax, totnlu, xrmax, xruse
        logical presrv, iout, xout
        real
     $          xx (xsize), cntl (10), rinfo (20), ax (anz)
 
c=== ums2f1 ============================================================
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
c  ums2f1 factorizes the n-by-n column-form matrix at the head of ii/xx
c  or in ap/ai/ax, and places its lu factors at the tail of ii/xx.  the
c  input matrix overwritten if it is located in ii/xx on input.  if
c  block-triangular-form (btf) is in use, this routine factorizes a
c  single diagonal block.
 
c=======================================================================
c  input:
c=======================================================================
c
c       n:              order of matrix (or order of diagonal block
c                       if btf is in use).
c       cp (1..n+1):    column pointers for input matrix
c       nblks:          number of diagonal blocks in btf form
c       isize:          size of ii
c       xsize:          size of xx
c       k1:             first index of this matrix (1 if btf not used)
c       icntl:          integer control parameters, see ums2in
c       cntl:           real control parameters, see ums2in
c       keep (6..8):    integer control parameters, see ums2in
c       iuse:           memory usage in index
c       xuse:           memory usage in value
c       rmax:           maximum ludegr seen so far (see ums2f2 for info)
c       cmax:           maximum ludegc seen so far (see ums2f2 for info)
c       totnlu:         total number of lu arrowheads constructed so far
c
c       if nblks>1 then:
c          cperm (1..n):        col permutation to btf
c          rperm (1..n):        row permutation to btf
c       else
c          cperm (1..n):        undefined on input
c          rperm (1..n):        undefined on input
c
c
c       presrv:         if true then input matrix is preserved
c
c       if presrv is true then:
c           an:                 order of preserved matrix (all blocks)
c           anz:                entries in preserved matrix
c           ap (1..an+1):       column pointers for preserved matrix
c           ai (1..anz):        row indices of preserved matrix
c           ax (1..anz):        values of preserved matrix
c                               the preserved matrix is not in btf form;
c                               it is in the orginal order.
c           if nblks > 1:
c               pr (1..n):      inverse row permutations to btf form
c               nzoff           entries in off-diagonal blocks
c                               seen so far
c           else
c               pr (1..n):      undefined on input
c
c           ii (1..isize):      undefined on input
c           xx (1..xsize):      undefined on input
c           cp (1..n+1):        undefined on input
c
c       else, if presrv is false:
c           an:                         1
c           anz:                        1
c           ii (1..cp (1) - 1):         unused
c           ii (cp (1) ... cp (n+1)-1): row indices of matrix to factor,
c                                       will be overwritten on output
c           ii (cp (n+1) ... isize):    unused on input
c
c           xx (1..cp (1) - 1):         unused
c           xx (cp (1) ... cp (n+1)-1): values of matrix to factorize,
c                                       will be overwritten on output
c           xx (cp (n+1) ... xsize):    unused on input
c                       if btf is in use, then ii and xx contain a
c                       single diagonal block.
 
c=======================================================================
c  output:
c=======================================================================
c
c       xx (xtail ... xsize), xtail,  ii (itail ... isize), itail:
c
c                       the lu factors of a single diagonal block.
c                       see ums2f2 for a description.
c
c       ii (cp1 ... itail-1):   undefined on output
c       xx (cp1 ... xtail-1):   undefined on output,
c                       where cp1 is equal to the value of cp (1)
c                       if presrv is false, or cp1 = 1 if presrv is
c                       true.
c
c       info:           integer informational output, see ums2fa
c       rinfo:          real informational output, see ums2fa
c       cperm (1..n):   the final col permutations, including btf
c       rperm (1..n):   the final row permutations, including btf
c
c       iuse:           memory usage in index
c       xuse:           memory usage in value
c       rmax:           maximum ludegr seen so far (see ums2f2 for info)
c       cmax:           maximum ludegc seen so far (see ums2f2 for info)
c       totnlu:         total number of lu arrowheads constructed so far
c
c       if nblks>1 and presrv:
c           nzoff       entries in off-diagonal blocks seen so far
c
c       iout:           true if ran out of integer memory in ums2f1
c       xout:           true if ran out of real memory in ums2f1
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutine:   ums2f0
c       subroutines called:     ums2f2
c       functions called:       max, sqrt
        intrinsic max, sqrt
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        integer cp1, pc, pend, pcol, cdeg, col, csiz, nz, xp, ip, is, p,
     $          dn, dsiz, wrksiz, i, clen, d1, d2, n2, row, cscal
        parameter (cscal = 9)
        real
     $          xn
 
c  original and expanded column-form:
c  ----------------------------------
c  cp1:     = value of cp (1) on input
c  pc:      pointer to integer part of expanded column-form matrix
c  pend:    column col ends here in the input column-form matrix
c  pcol:    column col starts here in the input column-form matrix
c  cdeg:    degree (number of entries) in a column
c  clen:    number of original entries in a column (= degree, here,
c           but decreases in ums2f2)
c  csiz:    size of the integer part of an expanded column (cdeg+cscal)
c  cscal:   = 9, the number of scalars in column data structure
c  nz:      number of entries in the diagonal block being factorized
c  xp:      pointer to real part of the expanded column
c  ip:      pointer to integer part of the expanded column
c
c  memory usage:
c  -------------
c  wrksiz:  size of integer workspace needed by ums2f2
c
c  "dense" columns: (converted to a prior, or artificial, frontal mat.):
c  ----------------
c  d1:      = keep (7), dense column control
c  d2:      = keep (8), dense column control
c  dn:      number of "dense" columns
c  dsiz:    a column is "dense" if it has more than dsiz entries
c  xn:      = sqrt (real (n))
c  n2:      = int (sqrt (real (n)))
c
c  other:
c  ------
c  row:     row index
c  col:     a column index
c  i:       loop index
c  p:       pointer
 
c=======================================================================
c  executable statements:
c=======================================================================
 
        iout = .false.
        xout = .false.
 
c-----------------------------------------------------------------------
c  count "dense" columns (they are treated as a priori frontal matrices)
c-----------------------------------------------------------------------
 
c       a column is "dense" if it has more than dsiz entries
        d1 = keep (7)
        d2 = keep (8)
        xn = n
        xn = sqrt (xn)
        n2 = xn
        dsiz = max (0, d1, d2 * n2)
        dn = 0
        if (presrv) then
           if (nblks .eq. 1) then
              do 10 col = 1, n
                 if (ap (col+1) - ap (col) .gt. dsiz) then
c                   this is a "dense" column
                    dn = dn + 1
                    endif
10            continue
           else
              do 40 col = 1, n
c                if col might be dense, check more carefully:
                 cdeg = ap (cperm (col) + 1)- ap (cperm (col))
                 if (cdeg .gt. dsiz) then
                    cdeg = 0
                    do 20 p = ap (cperm (col)), ap (cperm (col) + 1) -1
                       row = pr (ai (p))
                       if (row .ge. k1) then
                          cdeg = cdeg + 1
                          if (cdeg .gt. dsiz) then
c                            this is a "dense" column, exit out of loop
                             dn = dn + 1
                             go to 30
                          endif
                       endif
20                  continue
c                   loop exit label:
30                  continue
                 endif
40            continue
           endif
        else
           do col = 1, n
              if (cp (col+1) - cp (col) .gt. dsiz) then
c                this is a "dense" column
                 dn = dn + 1
              endif
           end do
        endif
 
c-----------------------------------------------------------------------
c  get size of workspaces to allocate from ii
c-----------------------------------------------------------------------
 
c       workspaces: wir (n), wic (n), wpr (n), wpc (n),
c       wm (n), head (n), rp (n+dn), wc (n+dn), wr (n+dn), wj (n)
        if (nblks .eq. 1) then
c          rperm (1..n) is used as wir (1..n), and
c          cperm (1..n) is used as wic (1..n) in ums2f2
           wrksiz = 8*n + 3*dn
        else
           wrksiz = 10*n + 3*dn
        endif
 
c-----------------------------------------------------------------------
c  construct the expanded column-form of the matrix or the diag. block
c-----------------------------------------------------------------------
 
        if (presrv) then
 
c          -------------------------------------------------------------
c          allocate space for wrksiz workspace and nz+cscal*n
c          integers and nz reals for the expanded column-form matrix.
c          -------------------------------------------------------------
 
           cp1 = 1
           xp = 1
           ip = 1 + wrksiz
           if (nblks .eq. 1) then
 
c             ----------------------------------------------------------
c             construct copy of entire matrix
c             ----------------------------------------------------------
 
              nz = anz
              is = nz + wrksiz + cscal*n
              iuse = iuse + is
              xuse = xuse + nz
              info (18) = max (info (18), iuse)
              info (19) = max (info (19), iuse)
              info (20) = max (info (20), xuse)
              info (21) = max (info (21), xuse)
              iout = is .gt. isize
              xout = nz .gt. xsize
              if (iout .or. xout) then
c                error return, if not enough integer and/or real memory:
                 go to 9000
              endif
 
              pc = ip
              do 70 col = 1, n
                 cp (col) = pc - wrksiz
                 cdeg = ap (col+1) - ap (col)
                 clen = cdeg
                 csiz = cdeg + cscal
                 ii (pc) = csiz
                 ii (pc+1) = cdeg
                 ii (pc+5) = 0
                 ii (pc+6) = clen
                 ii (pc+7) = 0
                 ii (pc+8) = 0
                 ii (pc+2) = xp
                 xp = xp + cdeg
                 pc = pc + cscal
                 p = ap (col)
                 do 60 i = 0, cdeg - 1
                    ii (pc + i) = ai (p + i)
60               continue
                 pc = pc + cdeg
70            continue
              do p = 1, nz
                 xx (p) = ax (p)
              end do
 
           else
 
c             ----------------------------------------------------------
c             construct copy of a single block in btf form
c             ----------------------------------------------------------
 
c             check for memory usage during construction of block
              do 100 col = 1, n
                 pc = ip
                 cp (col) = pc - wrksiz
                 ip = ip + cscal
                 iout = ip .gt. isize
                 if (iout) then
c                   error return, if not enough integer memory:
                    go to 9000
                 endif
                 ii (pc+2) = xp
                 cdeg = ip
                 do 90 p = ap (cperm (col)), ap (cperm (col)+1)-1
                    row = pr (ai (p))
                    if (row .ge. k1) then
                       iout = ip .gt. isize
                       xout = xp .gt. xsize
                       if (iout .or. xout) then
c                         error return, if not enough memory
                          go to 9000
                       endif
                       ii (ip) = row - k1 + 1
                       xx (xp) = ax (p)
                       ip = ip + 1
                       xp = xp + 1
                    else
c                      entry in off-diagonal part
                       nzoff = nzoff + 1
                    endif
90               continue
                 cdeg = ip - cdeg
                 clen = cdeg
                 csiz = cdeg + cscal
                 ii (pc) = csiz
                 ii (pc+1) = cdeg
                 ii (pc+5) = 0
                 ii (pc+6) = clen
                 ii (pc+7) = 0
                 ii (pc+8) = 0
100           continue
 
              nz = xp - 1
              is = nz + wrksiz + cscal*n
              iuse = iuse + is
              xuse = xuse + nz
              info (18) = max (info (18), iuse)
              info (19) = max (info (19), iuse)
              info (20) = max (info (20), xuse)
              info (21) = max (info (21), xuse)
 
           endif
 
c          -------------------------------------------------------------
c          get memory usage for next call to ums2rf
c          -------------------------------------------------------------
 
           xruse = xruse + nz
           xrmax = max (xrmax, xruse)
 
        else
 
c          -------------------------------------------------------------
c          allocate space for wrksiz workspace and additional cscal*n
c          space for the expanded column-form of the matrix.
c          -------------------------------------------------------------
 
           cp1 = cp (1)
           nz = cp (n+1) - cp1
           pc = cp1 + wrksiz + (nz+cscal*n)
           iuse = iuse + wrksiz + cscal*n
           info (18) = max (info (18), iuse)
           info (19) = max (info (19), iuse)
           iout = pc .gt. isize+1
           if (iout) then
c             error return, if not enough integer memory:
              go to 9000
           endif
 
c          -------------------------------------------------------------
c          expand the column form in place and make space for workspace
c          -------------------------------------------------------------
 
           xp = nz + 1
           ip = nz + cscal*n + 1
           pend = cp (n+1)
           do 120 col = n, 1, -1
              pcol = cp (col)
              do p = pend-1, pcol, -1
                 pc = pc - 1
                 ii (pc) = ii (p)
              end do
              pc = pc - cscal
              cdeg = pend - pcol
              clen = cdeg
              pend = pcol
              csiz = cdeg + cscal
              ip = ip - csiz
              cp (col) = ip
              ii (pc) = csiz
              ii (pc+1) = cdeg
              ii (pc+5) = 0
              ii (pc+6) = clen
              ii (pc+7) = 0
              ii (pc+8) = 0
              xp = xp - cdeg
              ii (pc+2) = xp
120        continue
        endif
 
c-----------------------------------------------------------------------
c  factorize the expanded column-form, with allocated workspaces
c-----------------------------------------------------------------------
 
        xp = cp1
        ip = cp1 + wrksiz
 
        if (nblks .eq. 1) then
 
c          pass rperm and cperm as the wir and wic arrays in ums2f2:
           call ums2f2 (cp, nz, n, 1, 1, 1, itail, xtail,
     $          xx (xp), xsize-xp+1, ii (ip), isize-ip+1, icntl, cntl,
     $          info, rinfo, .false., iuse, xuse,
     $          rperm, cperm, ii (cp1), ii (cp1+n),
     $          ii (cp1+2*n), ii (cp1+3*n), ii (cp1+4*n), ii (cp1+5*n),
     $          ii (cp1+6*n+dn), ii (cp1+7*n+2*dn),
     $          dn, dsiz, keep, rmax, cmax, totnlu, xrmax, xruse)
 
        else
 
c          pass cperm, rperm, wic and wir as separate arrays, and
c          change cperm and rperm from the btf permutations to the
c          final permutations (including btf and numerical pivoting).
           call ums2f2 (cp, nz, n, n, cperm, rperm, itail, xtail,
     $          xx (xp), xsize-xp+1, ii (ip), isize-ip+1, icntl, cntl,
     $          info, rinfo, .true., iuse, xuse,
     $          ii (cp1), ii (cp1+n), ii (cp1+2*n), ii (cp1+3*n),
     $          ii (cp1+4*n), ii (cp1+5*n), ii (cp1+6*n), ii (cp1+7*n),
     $          ii (cp1+8*n+dn), ii (cp1+9*n+2*dn),
     $          dn, dsiz, keep, rmax, cmax, totnlu, xrmax, xruse)
 
        endif
 
        if (info (1) .lt. 0) then
c          error return, if error occured in ums2f2:
           return
        endif
 
c-----------------------------------------------------------------------
c  adjust tail pointers, and save pointer to numerical part of lu
c-----------------------------------------------------------------------
 
c       head = cp1
        iuse = iuse - wrksiz
        itail = itail + ip - 1
        xtail = xtail + xp - 1
        ii (itail) = xtail
        return
 
c=======================================================================
c  error return
c=======================================================================
 
c       error return label:
9000    continue
        return
        end
