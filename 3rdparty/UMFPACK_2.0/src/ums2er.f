 
        subroutine ums2er (who, icntl, info, error, s)
c
cc UMS2ER prints error and warning messages and sets error flags.
c
        integer who, icntl (20), info (40), error, s
 
c=== ums2er ============================================================
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
c  print error and warning messages, and set error flags.
 
c=======================================================================
c  input:
c=======================================================================
c
c       who             which user-callable routine called:
c                       1: ums2fa, 2: ums2rf, 3: ums2so
c       icntl (1):      i/o unit for error and warning messages
c       icntl (3):      printing level
c       info (1):       the error/warning status
c       error:          the applicable error (<0) or warning (>0).
c                       see ums2p2 for a description.
c       s:              the relevant offending value
 
c=======================================================================
c  output:
c=======================================================================
c
c       info (1):       the error/warning status
 
c=======================================================================
c  subroutines and functions called / called by:
c=======================================================================
c
c       called by subroutines:  ums2co, ums2fa, ums2f0, ums2f1, ums2f2,
c                               ums2rf, ums2r0, ums2r2, ums2so, ums2s2
c       subroutines called:     ums2p2
c       functions called:       mod
        intrinsic mod
 
c=======================================================================
c  local scalars:
c=======================================================================
 
        logical both
        real
     $          zero
        parameter (zero = 0.0)
        integer ioerr, prl
 
c  ioerr:   i/o unit for error and warning messages
c  prl:     printing level
c  both:    if true, then combine errors -3 and -4 into error -5
 
c=======================================================================
c  executable statements:
c=======================================================================
 
        ioerr = icntl (1)
        prl = icntl (3)
        if (error .lt. 0) then
c          this is an error message
           both = (info (1) .eq. -3 .and. error .eq. -4) .or.
     $            (info (1) .eq. -4 .and. error .eq. -3)
           if (both) then
c             combine error -3 (out of integer memory) and error -4
c             (out of real memory)
              info (1) = -5
           else
              info (1) = error
           endif
           if (prl .ge. 1) then
              call ums2p2 (who, error, s, 0, zero, ioerr)
           endif
        else if (error .gt. 0) then
c          this is a warning message
           if (info (1) .ge. 0) then
c             do not override a prior error setting, sum up warnings
              if (mod (info (1) / error, 2) .eq. 0) then
                 info (1) = info (1) + error
              endif
           endif
           if (prl .ge. 2) then
              call ums2p2 (who, error, s, 0, zero, ioerr)
           endif
        endif
        return
        end
