 
        subroutine ums2p1 (who, where,
     $          n, ne, job, trans, lvalue, lindex, value,
     $          index, keep, cntl, icntl, info, rinfo,
     $          b, x, lx, w, lw)
c
cc UMS2P1 is a utility function which prints arguments for several other routines.
c
        integer who, where, n, ne, job, lvalue, lindex, index (lindex),
     $          keep (20), icntl (20), info (40), lx, lw
        real
     $          value (lvalue), cntl (10), rinfo (20), b (lx),
     $          x (lx), w (lw)
        logical trans
 
c=== ums2p1 ============================================================
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
c  print input/output arguments for ums2fa, ums2rf, and ums2so
 
c=======================================================================
c  installation note:
c=======================================================================
c
c  this routine can be deleted on installation (replaced with a dummy
c  routine that just returns without printing) in order to completely
c  disable the printing of all input/output parameters.  to completely
c  disable all i/o, you can also replace the ums2p2 routine with a
c  dummy subroutine.  if you make this modification, please do
c  not delete any original code - just comment it out instead.  add a
c  comment and date to your modifications.
 
c=======================================================================
c  input:
c=======================================================================
c
c       who:            what routine called ums2p1:
c                       1: ums2fa, 2: ums2rf, 3: ums2so
c       where:          called from where:
c                       1: entry of routine, else exit of routine
c       icntl (3):      if < 3 then print nothing, if 3 then print
c                       terse output, if >= 4 print everything
c       icntl (2):      i/o unit on which to print.  no printing
c                       occurs if < 0.
c
c       parameters to print, see ums2fa, ums2rf, or ums2so for
c       descriptions:
c
c           n, ne, job, trans, lvalue, lindex, value, index, keep,
c           icntl, info, rinfo, b, x, lx, w, lw
 
c=======================================================================
c  output:
c=======================================================================
c
c       on icntl (2) i/o unit only
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutines:  ums2fa, ums2rf, ums2so
c       functions called:       min
        intrinsic min
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        logical transa, transc, prlu, badlu, sglton, presrv, symbol
        integer io, prl, prn, k, lui1, lui2, lux1, lux2, row, col,
     $          facne, facn, nz, facjob, nblks, nzoff, factra, cpermp,
     $          rpermp, app, axp, aip, offip, offxp, lublpp, offpp,
     $          blkpp, p1, p2, p, blk, k1, k2, kn, luiip, luxxp, npiv,
     $          nlu, e, luk, lupp, luip, luxp, ludegr, ludegc, lunson,
     $          lusonp, lucp, lurp, i, j, nzcol, nzrow, uxp, son,
     $          prmax, ludimr, ludimc, maxdr, maxdc, luir1, ip1, ip2,
     $          xp1
        parameter (prmax = 10)
 
c  printing control:
c  -----------------
c  io:      i/o unit for diagnostic messages
c  prl:     printing level
c  prn:     number of entries printed so far
c  prmax:   maximum number of entries to print if prl = 3
c  prlu:    true if printing lu factors
c
c  location and status of lu factors:
c  ----------------------------------
c  transc:  transc argument in ums2so
c  transa:  transa argument in ums2fa or ums2rf when matrix factorized
c  badlu:   true if lu factors uncomputed or corrupted
c  presrv:  true if original matrix was preserved when factorized
c  symbol:  true if only symbolic part of lu factors needed on input
c  lui1:    integer part of lu factors start in index (lui1...)
c  luir1:   index (luir1 ... lui2) is needed for a call to ums2rf
c  lui2:    integer part of lu factors end in index (..lui2)
c  lux1:    real part of lu factors start in value (lux1...)
c  lux2:    real part of lu factors end in value (...lux1)
c  ip1:     pointer into leading part of lu factors in index
c  ip2:     pointer into trailing part of lu factors in index
c  xp1:     pointer into leading part of lu factors in value
c
c  arrays and scalars allocated in lu factors (in order):
c  ------------------------------------------------------
c  app:     ap (1..n+1) array located in index (app...app+n)
c  axp:     ax (1..nz) array located in value (axp...axp+nz-1)
c  aip:     ai (1..nz) array located in index (aip...aip+nz-1)
c  offip:   offi (1..nzoff) array loc. in index (offip...offip+nzoff-1)
c  offxp:   offx (1..nzoff) array loc. in value (offxp...offxp+nzoff-1)
c  ...      lu factors of each diagonal block located here
c  lublpp:  lublkp (1..nblks) array in index (lublpp..lublpp+nblks-1)
c  blkpp:   blkp (1..nblks+1) array loc. in index (blkpp...blkpp+nblks)
c  offpp:   offp (1..n+1) array located in index (offpp...offpp+n)
c  cpermp:  cperm (1..n) array located in index (cpermp...cpermp+n-1)
c  rpermp:  rperm (1..n) array located in index (rpermp...rpermp+n-1)
c  ...      seven scalars in index (lui2-6...lui2):
c  factra:  0/1 if transa argument was false/true in ums2fa or ums2rf
c  nzoff:   number of entries in off-diagonal part
c  nblks:   number of diagonal blocks
c  facjob:  job argument in ums2fa or ums2rf when matrix factorized
c  nz:      entries in a
c  facn:    n argument in ums2fa or ums2rf when matrix factorized
c  facne:   ne argument in ums2fa or ums2rf when matrix factorized
c
c  a single diagonal block and its lu factors:
c  -------------------------------------------
c  blk:     current diagonal block
c  k1,k2:   current diagonal is a (k1..k2, k1..k2)
c  kn:      order of current diagonal block (= k2-k1+1)
c  sglton:  true if current diagonal block is 1-by-1 (a singleton)
c  luiip:   lu factors of a diagonal block start in index (luiip...)
c  luxxp:   lu factors of a diagonal block start in value (luxxp...)
c  npiv:    number of pivots in a diagonal block (0 <= npiv <= kn)
c  nlu:     number of elements in a diagonal block
c  lupp:    lup (1..nlu) array located in index (lupp...lupp+nlu-1)
c
c  an element in the lu factors of a single diagonal block:
c  --------------------------------------------------------
c  e:       element
c  luk:     number of pivots in element e
c  luip:    integer part of element is in index (luip...)
c  luxp:    real part of element e is in value (luxp...)
c  ludegr:  row degree (number of columns) of u2 block in element e
c  ludegc:  column degree (number of rows) of l2 block in element e
c  lunson:  number of sons of element e in the assembly dag
c  lusonp:  list of sons of element e in index(lusonp...lusonp+lunson-1)
c  lucp:    column pattern (row indices) of l2 block in index (lucp..)
c  lurp:    row pattern (column indices) of u2 block in index (lurp..)
c  nzcol:   entries in a column of l, including unit diagonal
c  nzrow:   entries in a row of u, including non-unit diagonal
c  uxp:     a row of the u2 block located in value (uxp...)
c  son:     a son of the element e
c  ludimr:  row dimension (number of columns) in frontal matrix
c  ludimc:  column dimension (number of rows) in frontal matrix
c  maxdr:   largest ludimr for this block
c  maxdc:   largest ludimc for this block
c
c  other:
c  ------
c  row:     row index
c  col:     column index
c  k:       kth pivot, and general loop index
c  i, j:    loop indices
c  p:       pointer
c  p1:      column of a starts ai/ax (p1...), or row offi/x (p1...)
c  p2:      column of a ends in ai/ax (...p2), or row offi/x (...p2)
 
c=======================================================================
c  executable statements:
c       if (printing disabled on installation) return
c=======================================================================
 
c-----------------------------------------------------------------------
c  get printing control parameters
c-----------------------------------------------------------------------
 
        io = icntl (2)
        prl = icntl (3)
        if (prl .lt. 3 .or. io .lt. 0) then
c          printing has not been requested
           return
        endif
 
c-----------------------------------------------------------------------
c  who is this, and where.  determine if lu factors are to be printed
c-----------------------------------------------------------------------
 
        if (who .eq. 1) then
           if (where .eq. 1) then
              write (io, 200) 'ums2fa input:       '
              prlu = .false.
           else
              write (io, 200) 'ums2fa output:      '
              prlu = .true.
           endif
        else if (who .eq. 2) then
           if (where .eq. 1) then
              write (io, 200) 'ums2rf input:       '
              prlu = .true.
           else
              write (io, 200) 'ums2rf output:      '
              prlu = .true.
           endif
        else if (who .eq. 3) then
           if (where .eq. 1) then
              write (io, 200) 'ums2so input:       '
              prlu = .true.
           else
              write (io, 200) 'ums2so output:      '
              prlu = .false.
           endif
        endif
 
c-----------------------------------------------------------------------
c  print scalar input arguments: n, ne, job, trans, lvalue, lindex
c-----------------------------------------------------------------------
 
        if (where .eq. 1) then
           write (io, *)  '   scalar arguments:'
           write (io, *)  '      n:      ', n, ' : order of matrix a'
           if (who .eq. 3) then
c             ums2so:
c             was a or a^t factorized?
              lui2 = keep (5)
              transa = .false.
              if (lui2-6 .ge. 1 .and. lui2-6 .le. lindex) then
                 transa = index (lui2-6) .ne. 0
              endif
              transc = trans
              if (.not. transc) then
                 if (job .eq. 1) then
                    write (io, *)'      job:    ',job,' : solve p''lx=b'
                 else if (job .eq. 2) then
                    write (io, *)'      job:    ',job,' : solve uq''x=b'
                 else if (.not. transa) then
                    write (io, *)'      job:    ',job,' : solve ax=b',
     $                           ' (paq=lu was factorized)'
                 else
                    write (io, *)'      job:    ',job,' : solve a''x=b',
     $                           ' (pa''q=lu was factorized)'
                 endif
              else
                 if (job .eq. 1) then
                    write (io, *)'      job:    ',job,' : solve l''px=b'
                 else if (job .eq. 2) then
                    write (io, *)'      job:    ',job,' : solve qu''x=b'
                 else if (.not. transa) then
                    write (io, *)'      job:    ',job,' : solve a''x=b',
     $                           ' (paq=lu was factorized)'
                 else
                    write (io, *)'      job:    ',job,' : solve ax=b',
     $                           ' (pa''q=lu was factorized)'
                 endif
              endif
              if (transc) then
                 write (io, *)
     $                  '      transc:   .true. : see job above '
              else
                 write (io, *)
     $                  '      transc:   .false. : see job above '
              endif
           else
c             ums2fa or ums2rf:
              write (io, *)'      ne:     ', ne,' : entries in matrix a'
              if (job .eq. 1) then
                 write (io, *)
     $                  '      job:    ',job,' : matrix a preserved'
              else
                 write (io, *)
     $                  '      job:    ',job,' : matrix a not preserved'
              endif
              transa = trans
              if (transa) then
                 write (io, *)
     $                  '      transa:   .true. : factorize a transpose'
              else
                 write (io, *)
     $                  '      transa:   .false. : factorize a'
              endif
           endif
           write (io, *)
     $  '      lvalue: ',lvalue,' : size of value array'
           write (io, *)
     $  '      lindex: ',lindex,' : size of index array'
        endif
 
c-----------------------------------------------------------------------
c  print control parameters:  icntl, cntl, and keep (6..8)
c-----------------------------------------------------------------------
 
        if (where .eq. 1) then
           write (io, *)
     $  '   control parameters, normally initialized by ums2in:'
           write (io, *)
     $  '      icntl (1...8): integer control parameters'
           write (io, *)
     $  '      icntl (1): ',icntl (1),' : i/o unit for error and',
     $                    ' warning messages'
           write (io, *)
     $  '      icntl (2): ',io,' : i/o unit for diagnostics'
           write (io, *)
     $  '      icntl (3): ',prl,' : printing control'
           if (who .eq. 1) then
              if (icntl (4) .eq. 1) then
                 write (io, *)
     $  '      icntl (4): ',icntl (4),' : use block triangular',
     $                    ' form (btf)'
              else
                 write (io, *)
     $  '      icntl (4): ',icntl (4),' : do not permute to block',
     $                    ' triangular form (btf)'
              endif
              write (io, *)
     $  '      icntl (5): ',icntl (5),' : columns examined during',
     $                    ' pivot search'
              if (icntl (6) .ne. 0) then
                 write (io, *)
     $  '      icntl (6): ',icntl (6),' : preserve symmetry'
              else
                 write (io, *)
     $  '      icntl (6): ',icntl (6),' : do not preserve symmetry'
              endif
           endif
           if (who .ne. 3) then
              write (io, *)
     $  '      icntl (7): ',icntl (7),' : block size for dense matrix',
     $                    ' multiply'
           else
              write (io, *)
     $  '      icntl (8): ',icntl (8),' : maximum number',
     $                    ' of iterative refinement steps'
           endif
           if (who .eq. 1) then
              write (io, *)
     $  '      cntl (1...3): real control parameters'
              write (io, *)
     $  '      cntl (1):  ',cntl (1),' : relative pivot tolerance'
              write (io, *)
     $  '      cntl (2):  ',cntl (2),' : frontal matrix',
     $                    ' growth factor'
              write (io, *)
     $  '      keep (6...8): integer control parameters',
     $                          ' not normally modified by user'
              write (io, *)
     $  '      keep (6):  ',keep(6),' : largest positive integer'
              write (io, *)
     $  '      keep (7):  ',keep(7),' : dense row/col control, d1'
              write (io, *)
     $  '      keep (8):  ',keep(8),' : dense row/col control, d2'
           else if (who .eq. 3) then
              write (io, *)
     $  '      cntl (1...3): real control parameters'
              write (io, *)
     $  '      cntl (3):  ',cntl(3),' : machine epsilon'
           endif
        endif
 
c-----------------------------------------------------------------------
c  print the informational output
c-----------------------------------------------------------------------
 
        if (where .ne. 1) then
           write (io, *)
     $  '   output information:'
           write (io, *)
     $  '      info (1...24): integer output information'
           if (info (1) .lt. 0) then
              write (io, *)
     $  '      info (1):  ',info (1),' : error occurred!'
           else if (info (1) .gt. 0) then
              write (io, *)
     $  '      info (1):  ',info (1),' : warning occurred'
           else
              write (io, *)
     $  '      info (1):  ',info (1),' : no error or warning',
     $                    ' occurred'
           endif
           if (who .ne. 3) then
              write (io, *)
     $  '      info (2):  ',info (2),' : duplicate entries in a'
              write (io, *)
     $  '      info (3):  ',info (3),' : invalid entries in a',
     $                          ' (indices not in 1..n)'
              write (io, *)
     $  '      info (4):  ',info (4),' : invalid entries in a',
     $                          ' (not in prior pattern)'
              write (io, *)
     $  '      info (5):  ',info (5),' : entries in a after adding'
              write (io, *)
     $  '                       duplicates and removing invalid entries'
              write (io, *)
     $  '      info (6):  ',info (6),' : entries in diagonal',
     $                    ' blocks of a'
              write (io, *)
     $  '      info (7):  ',info (7),' : entries in off-diagonal',
     $                    ' blocks of a'
              write (io, *)
     $  '      info (8):  ',info (8),' : 1-by-1 diagonal blocks',
     $                    ' in a'
              write (io, *)
     $  '      info (9):  ',info (9),' : diagonal blocks in a',
     $                    ' (>1 only if btf used)'
              write (io, *)
     $  '      info (10): ',info (10),' : entries below diagonal in l'
              write (io, *)
     $  '      info (11): ',info (11),' : entries above diagonal in u'
              write (io, *)
     $  '      info (12): ',info (12),' : entries in l + u +',
     $                    ' offdiagonal blocks of a'
              write (io, *)
     $  '      info (13): ',info (13),' : frontal matrices'
              write (io, *)
     $  '      info (14): ',info (14),' : integer garbage',
     $                    ' collections'
              write (io, *)
     $  '      info (15): ',info (15),' : real garbage collections'
              write (io, *)
     $  '      info (16): ',info (16),' : diagonal pivots chosen'
              write (io, *)
     $  '      info (17): ',info (17),' : numerically valid pivots',
     $                    ' found in a'
              write (io, *)
     $  '      info (18): ',info (18),' : memory used in index'
              write (io, *)
     $  '      info (19): ',info (19),' : minimum memory needed in',
     $                    ' index'
              write (io, *)
     $  '      info (20): ',info (20),' : memory used in value'
              write (io, *)
     $  '      info (21): ',info (21),' : minimum memory needed in',
     $                    ' value'
              write (io, *)
     $  '      info (22): ',info (22),' : memory needed in',
     $                    ' index for next call to ums2rf'
              write (io, *)
     $  '      info (23): ',info (23),' : memory needed in',
     $                    ' value for next call to ums2rf'
           else
              write (io, *)
     $  '      info (24): ',info (24),' : steps of iterative',
     $                    ' refinement taken'
           endif
           if (who .ne. 3) then
              write (io, *)
     $  '      rinfo (1...8): real output information'
              write (io, *)
     $  '      rinfo (1): ',rinfo (1),' : total blas flop count'
              write (io, *)
     $  '      rinfo (2): ',rinfo (2),' : assembly flop count'
              write (io, *)
     $  '      rinfo (3): ',rinfo (3),' : pivot search flop count'
              write (io, *)
     $  '      rinfo (4): ',rinfo (4),' : level-1 blas flop count'
              write (io, *)
     $  '      rinfo (5): ',rinfo (5),' : level-2 blas flop count'
              write (io, *)
     $  '      rinfo (6): ',rinfo (6),' : level-3 blas flop count'
           else if (lw .eq. 4*n) then
              write (io, *)
     $  '      rinfo (1...8): real output information'
              write (io, *)
     $  '      rinfo (7): ',rinfo (7),' : sparse error estimate',
     $                    ' omega1'
              write (io, *)
     $  '      rinfo (8): ',rinfo (8),' : sparse error estimate',
     $                    ' omega2'
           endif
        endif
 
c-----------------------------------------------------------------------
c  print input matrix a, in triplet form, for ums2fa and ums2rf
c-----------------------------------------------------------------------
 
        if (where .eq. 1 .and. who .ne. 3) then
 
           if (transa) then
              write (io, *) '   the input matrix a transpose:'
              write (io, *)
     $  '      value (1 ... ',ne,' ): numerical values'
              write (io, *)
     $  '      index (1 ... ',ne,' ): column indices'
              write (io, *)
     $  '      index (',ne+1,' ... ',2*ne,' ): row indices'
              write (io, *) '   entries in the matrix a transpose',
     $                      ' (entry number: row, column, value):'
           else
              write (io, *) '   the input matrix a:'
              write (io, *)
     $  '      value (1 ... ',ne,' ): numerical values'
              write (io, *)
     $  '      index (1 ... ',ne,' ): row indices'
              write (io, *)
     $  '      index (',ne+1,' ... ',2*ne,' ): column indices'
              write (io, *) '   entries in the matrix a',
     $                      ' (entry number: row, column, value):'
           endif
 
           prn = min (prmax, ne)
           if (prl .ge. 4) then
              prn = ne
           endif
           do k = 1, prn
              if (transa) then
                 row = index (k+ne)
                 col = index (k)
              else
                 row = index (k)
                 col = index (k+ne)
              endif
              write (io, *) '      ', k, ': ',row,' ',col,' ', value (k)
           end do
           if (prn .lt. ne) then
              write (io, 220)
           endif
        endif
 
c-----------------------------------------------------------------------
c  print the lu factors:  ums2fa output, ums2rf input/output,
c                         and ums2so input
c-----------------------------------------------------------------------
 
        if (prlu .and. info (1) .lt. 0) then
           write (io, *) '   lu factors not printed because of error',
     $                   ' flag, info (1) = ', info (1)
           prlu = .false.
        endif
 
        if (prlu) then
 
c          -------------------------------------------------------------
c          description of what must be preserved between calls
c          -------------------------------------------------------------
 
           lux1 = keep (1)
           lux2 = keep (2)
           lui1 = keep (3)
           luir1 = keep (4)
           lui2 = keep (5)
 
           xp1 = lux1
           ip1 = lui1
           ip2 = lui2
 
c          -------------------------------------------------------------
c          on input to ums2rf, only the symbol information is used
c          -------------------------------------------------------------
 
           symbol = who .eq. 2 .and. where .eq. 1
 
           if (symbol) then
              write (io, *)
     $  '   keep (4...5) gives the location of lu factors'
              write (io, *)
     $  '      which must be preserved for calls to ums2rf: '
           else
              write (io, *)
     $  '   keep (1...5) gives the location of lu factors'
              write (io, *)
     $  '      which must be preserved for calls to ums2so: '
              write (io, *)
     $  '         value ( keep (1): ', lux1,' ... keep (2): ', lux2,' )'
              write (io, *)
     $  '         index ( keep (3): ', lui1,' ... keep (5): ', lui2,' )'
              write (io, *)
     $  '      and for calls to ums2rf: '
           endif
           write (io, *)
     $  '         index ( keep (4): ',luir1,' ... keep (5): ', lui2,' )'
 
           badlu = luir1 .le. 0 .or. lui2-6 .lt. luir1 .or.
     $        lui2 .gt. lindex
           if (.not. symbol) then
              badlu = badlu .or. lux1 .le. 0 .or.
     $        lux1 .gt. lux2 .or. lux2 .gt. lvalue .or. lui1 .le. 0 .or.
     $        luir1 .lt. lui1 .or. luir1 .gt. lui2
           endif
 
c          -------------------------------------------------------------
c          get the 7 scalars, and location of permutation vectors
c          -------------------------------------------------------------
 
           if (badlu) then
c             pointers are bad, so these values cannot be obtained
              facne  = 0
              facn   = 0
              nz     = 0
              facjob = 0
              nblks  = 0
              nzoff  = 0
              factra = 0
           else
              facne  = index (lui2)
              facn   = index (lui2-1)
              nz     = index (lui2-2)
              facjob = index (lui2-3)
              nblks  = index (lui2-4)
              nzoff  = index (lui2-5)
              factra = index (lui2-6)
           endif
 
           presrv = facjob .ne. 0
           transa = factra .ne. 0
           rpermp = (lui2-6) - (facn)
           cpermp = rpermp - (facn)
           ip2 = cpermp - 1
 
           if (symbol) then
              write (io, *)'   layout of lu factors in index:'
           else
              write (io, *)'   layout of lu factors in value and index:'
           endif
 
c          -------------------------------------------------------------
c          print location of preserved input matrix
c          -------------------------------------------------------------
 
           if (presrv) then
c             preserved column-form of original matrix
              app = ip1
              aip = app + (facn+1)
              ip1 = aip + (nz)
              axp = xp1
              xp1 = xp1 + (nz)
              if (.not. symbol) then
                 write (io, *)'      preserved copy of original matrix:'
                 write (io, *)
     $  '         index ( ',app,' ... ', aip-1,' ): column pointers'
                 write (io, *)
     $  '         index ( ',aip,' ... ',ip1-1,' ): row indices'
                 write (io, *)
     $  '         value ( ',axp,' ... ',xp1-1,' ): numerical values'
              endif
           else
              if (.not. symbol) then
                 write (io, *) '      original matrix not preserved.'
              endif
           endif
 
           badlu = badlu .or.
     $          n .ne. facn .or. nz .le. 0 .or. luir1 .gt. ip2 .or.
     $          nblks .le. 0 .or. nblks .gt. n
           if (.not. symbol) then
              badlu = badlu .or. xp1 .gt. lux2 .or. nzoff .lt. 0
           endif
           if (badlu) then
              nblks = 0
           endif
 
           if (nblks .le. 1) then
 
c             ----------------------------------------------------------
c             single block (or block triangular form not used),
c             or lu factors are corrupted
c             ----------------------------------------------------------
 
              write (io, *)
     $  '      collection of elements in lu factors:'
              write (io, *)
     $  '         (an "element" contains one or columns of l and'
              write (io, *)
     $  '         rows of u with similar nonzero pattern)'
              write (io, *)
     $  '         index ( ',luir1,' ... ', ip2,
     $                                      ' ): integer data'
              if (.not. symbol) then
                 write (io, *)
     $  '         value ( ',xp1,' ... ', lux2,' ): numerical values'
              endif
 
           else
 
c             ----------------------------------------------------------
c             block triangular form with more than one block
c             ----------------------------------------------------------
 
              offip = ip1
              ip1 = ip1 + (nzoff)
              offxp = xp1
              xp1 = xp1 + (nzoff)
              offpp = cpermp - (n+1)
              blkpp = offpp - (nblks+1)
              lublpp = blkpp - (nblks)
              ip2 = lublpp - 1
              badlu = badlu .or. luir1 .gt. ip2
              if (.not. symbol) then
                 badlu = badlu .or. ip1 .gt. ip2 .or.
     $           xp1 .gt. lux2 .or. luir1 .ne. ip1
              endif
              write (io, *)
     $  '      matrix permuted to upper block triangular form.'
              if (nzoff .ne. 0 .and. .not. symbol) then
                 write (io, *) '      entries not in diagonal blocks:'
                 write (io, *)
     $  '         index ( ',offip,' ... ',luir1-1,' ): row indices'
                 write (io, *)
     $  '         value ( ',offxp,' ... ',xp1-1,' ): numerical values'
              endif
              write (io, *)
     $  '      collection of elements in lu factors of diagonal blocks:'
              write (io, *)
     $  '         (an "element" contains one or columns of l and'
              write (io, *)
     $  '         rows of u with similar nonzero pattern)'
              if (luir1 .le. lublpp-1) then
                 write (io, *)
     $  '         index ( ',luir1,' ... ', ip2,
     $                                         ' ): integer data'
              endif
              if (xp1 .le. lux2 .and. .not. symbol) then
                 write (io, *)
     $  '         value ( ',xp1,' ... ', lux2,' ): numerical values'
              endif
              write (io, *) '      other block triangular data:'
              write (io, *)
     $  '         index ( ',lublpp,' ... ',blkpp-1,' ):',
     $                          ' pointers to block factors'
              write (io, *)
     $  '         index ( ', blkpp,' ... ',offpp-1,' ):',
     $                          ' index range of blocks'
              if (.not. symbol) then
                 write (io, *)
     $  '         index ( ', offpp,' ... ',   lui2-7,' ):',
     $          ' row pointers for off-diagonal part'
              endif
           endif
 
c          -------------------------------------------------------------
c          print location of permutation vectors and 7 scalars at tail
c          -------------------------------------------------------------
 
           write (io, *)
     $  '      permutation vectors (start at keep(4)-2*n-6):'
           write (io, *)
     $  '         index ( ', cpermp,' ... ',rpermp-1,' ):',
     $                          ' column permutations'
           write (io, *)
     $  '         index ( ', rpermp,' ... ',lui2-7,' ):',
     $                          ' row permutations'
           write (io, *) '      other data in index: '
           write (io, *)
     $  '         index ( ',lui2-6,' ): ', factra,' :',
     $                                  ' transa ums2fa/ums2rf argument'
           write (io, *)
     $  '         index ( ',lui2-5,' ): ', nzoff,' :',
     $                                  ' entries in off-diagonal part'
           write (io, *)
     $  '         index ( ',lui2-4,' ): ', nblks,' :',
     $                                  ' number of diagonal blocks'
           write (io, *)
     $  '         index ( ',lui2-3,' ): ', facjob,' :',
     $                                  ' job ums2fa/ums2rf argument'
           write (io, *)
     $  '         index ( ',lui2-2,' ): ', nz,' :',
     $                                  ' entries in original matrix'
           write (io, *)
     $  '         index ( ',lui2-1,' ): ', facn,' :',
     $                                  ' n ums2fa/ums2rf argument'
           write (io, *)
     $  '         index ( ',lui2  ,' ): ', facne,' :',
     $                                  ' ne ums2fa/ums2rf argument'
 
           if (.not. symbol) then
              badlu = badlu .or. ip1 .ne. luir1
           endif
           ip1 = luir1
           if (badlu) then
              write (io, *) '   lu factors uncomputed or corrupted!'
              presrv = .false.
              nblks = 0
           endif
 
c          -------------------------------------------------------------
c          copy of original matrix in column-oriented form
c          -------------------------------------------------------------
 
           if (presrv .and. .not. symbol) then
              write (io, 230)
              write (io, *) '   preserved copy of original matrix:'
              do 20 col = 1, n
                 p1 = index (app-1 + col)
                 p2 = index (app-1 + col+1) - 1
                 write (io, *) '      col: ', col, ' nz: ', p2-p1+1
                 if (prl .eq. 3) then
                    p2 = min (prmax, p2)
                 endif
                 write (io, *) (index (aip-1 + p), p = p1, p2)
                 write (io, 210) (value (axp-1 + p), p = p1, p2)
                 if (prl .eq. 3 .and. p2 .ge. prmax) then
c                   exit out of loop if done printing:
                    go to 30
                 endif
20            continue
c             loop exit label:
30            continue
              if (prl .eq. 3 .and. nz .gt. prmax) then
                 write (io, 220)
              endif
           endif
 
c          -------------------------------------------------------------
c          entries in off-diagonal blocks, in row-oriented form
c          -------------------------------------------------------------
 
           if (nblks .gt. 1 .and. .not. symbol) then
              write (io, 230)
              write (io, *)'   entries not in diagonal blocks:'
              if (nzoff .eq. 0) then
                 write (io, *)'      (none)'
              endif
              do 40 row = 1, n
                 p1 = index (offpp-1 + row)
                 p2 = index (offpp-1 + row+1) - 1
                 if (p2 .ge. p1) then
                    write (io, *) '      row: ', row, ' nz: ',p2-p1+1
                    if (prl .eq. 3) then
                       p2 = min (prmax, p2)
                    endif
                    write (io, *) (index (offip-1 + p), p = p1, p2)
                    write (io, 210) (value (offxp-1 + p), p = p1, p2)
                 endif
                 if (prl .eq. 3 .and. p2 .ge. prmax) then
c                   exit out of loop if done printing:
                    go to 50
                 endif
40            continue
c             loop exit label:
50            continue
              if (prl .eq. 3 .and. nz .gt. prmax) then
                 write (io, 220)
              endif
           endif
 
c          -------------------------------------------------------------
c          lu factors of each diagonal block
c          -------------------------------------------------------------
 
           write (io, 230)
           if (nblks .gt. 0) then
              write (io, *) '   lu factors of each diagonal block:'
           endif
           prn = 0
           do 140 blk = 1, nblks
 
c             ----------------------------------------------------------
c             print the factors of a single diagonal block
c             ----------------------------------------------------------
 
              if (nblks .gt. 1) then
                 k1 = index (blkpp-1 + blk)
                 k2 = index (blkpp-1 + blk+1) - 1
                 kn = k2-k1+1
                 sglton = kn .eq. 1
                 if (sglton) then
c                   this is a singleton
                    luxxp = xp1-1 + index (lublpp-1 + blk)
                 else
                    luiip = ip1-1 + index (lublpp-1 + blk)
                 endif
              else
                 sglton = .false.
                 k1 = 1
                 k2 = n
                 kn = n
                 luiip = ip1
              endif
 
              write (io, 240)
              if (sglton) then
 
c                -------------------------------------------------------
c                this is a singleton
c                -------------------------------------------------------
 
                 write (io, *)'   singleton block: ', blk,
     $                        ' at index : ', k1
                 if (.not. symbol) then
                    write (io, *)
     $           '   located in value ( ', luxxp,' ): ', value (luxxp)
                 endif
                 if (prl .eq. 3 .and. prn .gt. prmax) then
c                   exit out of loop if done printing:
                    go to 150
                 endif
                 prn = prn + 1
 
              else
 
c                -------------------------------------------------------
c                this block is larger than 1-by-1
c                -------------------------------------------------------
 
                 luxxp = xp1-1 + index (luiip)
                 nlu = index (luiip+1)
                 npiv = index (luiip+2)
                 maxdc = index (luiip+3)
                 maxdr = index (luiip+4)
                 lupp = luiip+5
                 write (io, *) '   block: ',blk,' first index: ',k1,
     $                         ' last index: ',k2, '   order: ', kn
                 write (io, *) '   elements: ', nlu, '   pivots: ', npiv
                 write (io, *) '   largest contribution block: ',
     $                         maxdc, ' by ', maxdr
                 write (io, *) '   located in index ( ',luiip,' ... )'
                 if (.not. symbol) then
                    write (io, *) '   and in value ( ',luxxp,' ... )'
                 endif
                 luiip = lupp + nlu
 
c                note: the indices of the lu factors of the block range
c                from 1 to kn, even though the kn-by-kn block resides in
c                a (k1 ... k2, k1 ... k2).
                 k = 0
 
                 do 130 e = 1, nlu
 
c                   ----------------------------------------------------
c                   print a single element
c                   ----------------------------------------------------
 
                    luip = luiip-1 + index (lupp-1 + e)
                    luxp = luxxp-1 + index (luip)
                    luk  = index (luip+1)
                    ludegr = index (luip+2)
                    ludegc = index (luip+3)
                    lunson = index (luip+4)
                    ludimr = index (luip+5)
                    ludimc = index (luip+6)
                    lucp = luip + 7
                    lurp = lucp + ludegc
                    lusonp = lurp + ludegr
                    write (io, *) '      e: ',e, ' pivots: ', luk,
     $                  ' children in dag: ', lunson,
     $                  ' frontal matrix: ', ludimr, ' by ', ludimc
 
c                   ----------------------------------------------------
c                   print the columns of l
c                   ----------------------------------------------------
 
                    p = luxp
                    do 80 j = 1, luk
                       col = k+j
                       nzcol = luk-j+1+ludegc
                       write (io, *) '         col: ',col,' nz: ',nzcol
c                      l is unit diagonal
                       prn = prn + 1
                       row = col
                       if (symbol) then
                          write (io, *) '            ', row
                       else
                          write (io, *) '            ', row, '  1.0'
                       endif
                       p = p + 1
c                      pivot block
                       do 60 i = j+1, luk
                          if (prl.eq.3 .and. prn.gt.prmax) then
c                            exit out of loop if done printing:
                             go to 150
                          endif
                          prn = prn + 1
                          row = k+i
                          if (symbol) then
                             write (io, *)'            ', row
                          else
                             write (io, *)'            ', row, value (p)
                          endif
                          p = p + 1
60                     continue
c                      l block
                       do 70 i = 1, ludegc
                          if (prl.eq.3 .and. prn.gt.prmax) then
c                            exit out of loop if done printing:
                             go to 150
                          endif
                          prn = prn + 1
                          row = index (lucp-1+i)
                          if (symbol) then
                             write (io, *)'            ', row
                          else
                             write (io, *)'            ', row, value (p)
                          endif
                          p = p + 1
70                     continue
                       p = p + j
80                  continue
 
c                   ----------------------------------------------------
c                   print the rows of u
c                   ----------------------------------------------------
 
                    uxp = luxp + luk*(ludegc+luk)
                    do 110 i = 1, luk
                       row = k+i
                       nzrow = luk-i+1+ludegr
                       write (io, *) '         row: ',row,' nz: ',nzrow
                       p = luxp + (i-1) + (i-1) * (ludegc+luk)
c                      pivot block
                       do 90 j = i, luk
                          if (prl.eq.3 .and. prn.gt.prmax) then
c                            exit out of loop if done printing:
                             go to 150
                          endif
                          prn = prn + 1
                          col = k+j
                          if (symbol) then
                             write (io, *)'            ', col
                          else
                             write (io, *)'            ', col, value (p)
                          endif
                          p = p + (ludegc+luk)
90                     continue
                       p = uxp
c                      u block
                       do 100 j = 1, ludegr
                          if (prl.eq.3 .and. prn.gt.prmax) then
c                            exit out of loop if done printing:
                             go to 150
                          endif
                          prn = prn + 1
                          col = index (lurp-1+j)
                          if (symbol) then
                             write (io, *)'            ', col
                          else
                             write (io, *)'            ', col, value (p)
                          endif
                          p = p + luk
100                    continue
                       uxp = uxp + 1
110                 continue
 
c                   ----------------------------------------------------
c                   print the sons of the element in the assembly dag
c                   ----------------------------------------------------
 
                    do 120 i = 1, lunson
                       prn = prn + 1
                       son = index (lusonp-1+i)
                       if (son .le. kn) then
c                         an luson
                          write (io, *) '         luson: ', son
                       else if (son .le. 2*kn) then
c                         a uson
                          write (io, *) '         uson:  ', son-kn
                       else
c                         an lson
                          write (io, *) '         lson:  ', son-2*kn
                       endif
120                 continue
 
c                   ----------------------------------------------------
c                   increment count of pivots within this block
c                   ----------------------------------------------------
 
                    k = k + luk
130              continue
              endif
140        continue
c          if loop was not exited prematurely, do not print "..." :
           go to 160
c          loop exit label:
150        continue
           write (io, 220)
160        continue
 
c          -------------------------------------------------------------
c          row and column permutations
c          -------------------------------------------------------------
 
           if (.not. badlu) then
              write (io, 230)
              write (io, *) '      column permutations'
              if (prl .ge. 4 .or. n .le. prmax) then
                 write (io, *) (index (cpermp+i-1), i = 1, n)
              else
                 write (io, *) (index (cpermp+i-1), i = 1, prmax)
                 write (io, 220)
              endif
 
              write (io, 230)
              write (io, *) '      row permutations'
              if (prl .ge. 4 .or. n .le. prmax) then
                 write (io, *) (index (rpermp+i-1), i = 1, n)
              else
                 write (io, *) (index (rpermp+i-1), i = 1, prmax)
                 write (io, 220)
              endif
           endif
 
        endif
 
c-----------------------------------------------------------------------
c  print b (on input) or w and x (on output) for ums2so
c-----------------------------------------------------------------------
 
        if (who .eq. 3) then
           write (io, 230)
           prn = min (prmax, n)
           if (prl .ge. 4) then
c             print all of b, or w and x
              prn = n
           endif
           if (where .eq. 1) then
              write (io, *) '   w (1 ... ',lw,' ), work vector:',
     $                      ' not printed'
              write (io, *) '   b (1 ... ',n,' ), right-hand side: '
              do 170 i = 1, prn
                 write (io, *) '      ', i, ': ', b (i)
170           continue
              if (prn .lt. n) then
                 write (io, 220)
              endif
           else
              if (info (1) .lt. 0) then
                 write (io, *) '   w (1 ... ',lw,' ), work vector, and'
                 write (io, *) '   x (1 ... ',n,' ), solution,'
                 write (io, *) '      not printed because of error',
     $                         ' flag, info (1) = ', info (1)
              else
                 if (lw .eq. 4*n) then
c                   ums2so did iterative refinement
                    write (io, *) '   w (1 ... ',n,' ), residual: '
                    do 180 i = 1, prn
                       write (io, *) '      ', i, ': ', w (i)
180                 continue
                    if (prn .lt. n) then
                       write (io, 220)
                    endif
                    write (io, *) '   w (',n+1,' ... ',lw,' )',
     $                            ', work vector: not printed'
                 else
c                   no iterative refinement
                    write (io, *) '   w (1 ... ',lw,' ),',
     $                            ' work vector: not printed'
                 endif
                 write (io, *) '   x (1 ... ',n,' ), solution: '
                 do 190 i = 1, prn
                    write (io, *) '      ', i, ': ', x (i)
190              continue
                 if (prn .lt. n) then
                    write (io, 220)
                 endif
              endif
           endif
        endif
 
c-----------------------------------------------------------------------
c  who is this, and where:
c-----------------------------------------------------------------------
 
        if (who .eq. 1) then
           if (where .eq. 1) then
              write (io, 200) 'end of ums2fa input '
           else
              write (io, 200) 'end of ums2fa output'
           endif
        else if (who .eq. 2) then
           if (where .eq. 1) then
              write (io, 200) 'end of ums2rf input '
           else
              write (io, 200) 'end of ums2rf output'
           endif
        else if (who .eq. 3) then
           if (where .eq. 1) then
              write (io, 200) 'end of ums2so input '
           else
              write (io, 200) 'end of ums2so output'
           endif
        endif
 
        return
 
c-----------------------------------------------------------------------
c  format statments
c-----------------------------------------------------------------------
 
200     format (60('='), a20)
210     format (5e16.8)
220     format ('        ...')
230     format ('   ', 77 ('-'))
240     format ('   ', 77 ('.'))
        end
