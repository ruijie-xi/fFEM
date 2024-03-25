 
        subroutine ums2p2 (who, error, i, j, x, io)
c
cc UMS2P2 is a utility which prints error and warning messages for several functions.
c
        integer who, error, i, j, io
        real
     $          x
 
c=== ums2p2 ============================================================
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
c  not user-callable
 
c=======================================================================
c  description:
c=======================================================================
c
c  print error and warning messages for for ums2fa, ums2rf, and ums2so.
 
c=======================================================================
c  installation note:
c=======================================================================
c
c  this routine can be deleted on installation (replaced with a dummy
c  routine that just returns without printing) in order to completely
c  disable the printing of all error and warning messages.  the error
c  and warning return flag (info (1)) will not be affected.  to
c  completely disable all i/o, you can also replace the ums2p1 routine
c  with a dummy subroutine.  if you make this modification, please do
c  not delete any original code - just comment it out instead.  add a
c  comment and date to your modifications.
 
c=======================================================================
c  input:
c=======================================================================
c
c       who:            what user-callable routine called ums2p2:
c                       1: ums2fa, 2: ums2rf, 3: ums2so
c       i, j, x:        the relevant offending value(s)
c       io:             i/o unit on which to print.  no printing
c                       occurs if < 0.
c       error:          the applicable error (<0) or warning (>0)
c                       errors (<0) cause the factorization/solve to
c                       be terminated.  if an error occurs, a prior
c                       warning status is overwritten with the error
c                       status.
c
c  the following error codes are returned in info (1) by ums2er.
c  these errors cause the factorization or solve to terminate:
c
c  where**      error   description
c
c  fa rf  -     -1      n < 1 or n > maximum value
c  fa rf  -     -2      ne < 1 or ne > maximum value
c  fa rf  -     -3      lindex too small
c  fa rf  -     -4      lvalue too small
c  fa rf  -     -5      both lindex and lvalue are too small
c   - rf  -     -6      prior pivot ordering no longer acceptable
c   - rf so     -7      lu factors are uncomputed, or are corrupted
c
c  the following warning codes are returned in info (1) by ums2er.
c  the factorization or solve was able to complete:
c
c  fa rf  -     1       invalid entries
c  fa rf  -     2       duplicate entries
c  fa rf  -     3       invalid and duplicate entries
c  fa rf  -     4       singular matrix
c  fa rf  -     5       invalid entries, singular matrix
c  fa rf  -     6       duplicate entries, singular matrix
c  fa rf  -     7       invalid and duplicate entries, singular matrix
c   -  - so     8       iterative refinement cannot be done
c
c  the following are internal error codes (not returned in info (1))
c  for printing specific invalid or duplicate entries.  these codes are
c  for ums2co, ums2of, and ums2r2.  invalid and duplicate entries are
c  ignored during factorization.  warning levels (1..7) will be set
c  later by ums2er, above.
c
c  fa rf  -     99      invalid entry, out of range 1..n
c  fa rf  -     98      duplicate entry
c   - rf  -     97      invalid entry:  within a diagonal block, but not
c                       in the pattern of the lu factors of that block.
c   - rf  -     96      invalid entry:  below the diagonal blocks.  can
c                       only occur if the matrix has been ordered into
c                       block-upper-triangular form.
c   - rf  -     95      invalid entry:  matrix is singular.  the
c                       remaining rank 0 submatrix yet to be factorized
c                       is replaced with the identity matrix in the lu
c                       factors.  any entry that remains is ignored.
 
c ** fa: ums2fa, rf: ums2rf, so: ums2so
 
c=======================================================================
c  output:
c=======================================================================
c
c  error or warning message printed on i/o unit
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutines:  ums2er, ums2co, ums2r0, ums2r2
 
c=======================================================================
c  executable statements:
c       if (printing disabled on installation) return
c=======================================================================
 
        if (io .lt. 0) then
c          printing of error / warning messages has not been requested
           return
        endif
 
        if (who .eq. 1) then
 
c          -------------------------------------------------------------
c          ums2fa error messages
c          -------------------------------------------------------------
 
           if (error .eq. -1) then
              if (i .lt. 0) then
                 write (io, *)'ums2fa: n less than one!'
              else
                 write (io, *)'ums2fa: n too large, must be <= ', i,'!'
              endif
           else if (error .eq. -2) then
              write (io, *)'ums2fa: ne less than one!'
           else if (error .eq. -3) then
              write (io, *)'ums2fa: insufficient integer workspace! ',
     $                     ' lindex must be >= ', i, '.'
           else if (error .eq. -4) then
              write (io, *)'ums2fa: insufficient real workspace!    ',
     $                     ' lvalue must be >= ', i, '.'
 
c          -------------------------------------------------------------
c          ums2fa cumulative warning messages
c          -------------------------------------------------------------
 
           else if (error .eq. 1) then
              write (io, *)'ums2fa: ', i,' invalid entries ignored',
     $                     ' (out of range 1..n).'
           else if (error .eq. 2) then
              write (io, *)'ums2fa: ', i,' duplicate entries summed.'
           else if (error .eq. 4) then
              write (io, *)'ums2fa: matrix is singular.  only ', i,
     $                     ' pivots found.'
 
c          -------------------------------------------------------------
c          ums2fa non-cumulative warning messages (internal error codes)
c          -------------------------------------------------------------
 
           else if (error .eq. 99) then
              write (io, *)'ums2fa: invalid entry (out of range 1..n):'
              write (io, *)'        row: ',i,' col: ',j,' ',x
           else if (error .eq. 98) then
              write (io, *)'ums2fa: duplicate entry summed:'
              write (io, *)'        row: ',i,' col: ',j,' ',x
           endif
 
        else if (who .eq. 2) then
 
c          -------------------------------------------------------------
c          ums2rf error messages
c          -------------------------------------------------------------
 
           if (error .eq. -1) then
              write (io, *)'ums2rf: n less than one!'
           else if (error .eq. -2) then
              if (i .lt. 0) then
                 write (io, *)'ums2rf: ne less than one!'
              else
                 write (io, *)'ums2rf: ne too large, must be <= ',i, '!'
              endif
           else if (error .eq. -3) then
              write (io, *)'ums2rf: insufficient integer workspace! ',
     $                       ' lindex must be >= ', i, '.'
           else if (error .eq. -4) then
              write (io, *)'ums2rf: insufficient real workspace!    ',
     $                       ' lvalue must be >= ', i, '.'
           else if (error .eq. -6) then
              write (io, *)'ums2rf: pivot order from ums2fa failed!'
           else if (error .eq. -7) then
              write (io, *)'ums2rf: lu factors uncomputed,',
     $                     ' or corrupted!'
 
c          -------------------------------------------------------------
c          ums2rf cumulative warning messages
c          -------------------------------------------------------------
 
           else if (error .eq. 1) then
              if (i .gt. 0) then
                 write (io, *)'ums2rf: ', i,' invalid entries ignored',
     $                        ' (out of range 1..n).'
              else
                 write (io, *)'ums2rf: ',-i,' invalid entries ignored',
     $                        ' (not in prior pattern).'
              endif
           else if (error .eq. 2) then
              write (io, *)'ums2rf: ', i,' duplicate entries summed.'
           else if (error .eq. 4) then
              write (io, *)'ums2rf: matrix is singular.  only ', i,
     $                     ' pivots found.'
 
c          -------------------------------------------------------------
c          ums2rf non-cumulative warning messages (internal error codes)
c          -------------------------------------------------------------
 
           else if (error .eq. 99) then
              write (io, *)'ums2rf: invalid entry (out of range 1..n):'
              write (io, *)'        row: ',i,' col: ',j,' ',x
           else if (error .eq. 98) then
              write (io, *)'ums2rf: duplicate entry summed:'
              write (io, *)'        row: ',i,' col: ',j,' ',x
           else if (error .eq. 97) then
              write (io, *)'ums2rf: invalid entry (not in pattern of',
     $                     ' lu factor(s) of diagonal block(s)):'
              write (io, *)'        row: ',i,' col: ',j,' ',x
           else if (error .eq. 96) then
              write (io, *)'ums2rf: invalid entry (below diagonal',
     $                     ' blocks):'
              write (io, *)'        row: ',i,' col: ',j,' ',x
           else if (error .eq. 95) then
              write (io, *)'ums2rf: invalid entry (because',
     $                     ' matrix factorized by ums2fa was singular):'
              write (io, *)'        row: ',i,' col: ',j,' ',x
           endif
 
        else if (who .eq. 3) then
 
c          -------------------------------------------------------------
c          ums2so error messages
c          -------------------------------------------------------------
 
           if (error .eq. -7) then
              write (io, *)'ums2so: lu factors uncomputed,',
     $                     ' or corrupted!'
 
c          -------------------------------------------------------------
c          ums2so non-cumulative warning messages
c          -------------------------------------------------------------
 
           else if (error .eq. 8) then
              if (i .eq. 0) then
                 write (io, *)'ums2so: iterative refinement requested',
     $           ' but original matrix not preserved.'
              else
                 write (io, *)'ums2so: iterative refinement requested',
     $           ' but only available for ax=b or a''x=b.'
              endif
           endif
 
        endif
        return
        end
