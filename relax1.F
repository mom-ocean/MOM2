      subroutine relax1 (npt, variable, bc_symm
     &,                  guess, dpsi, forc, res
     &,                  cf
     &,                  sor, mxscan, mscan, crit
     &,                  imask, iperm, jperm, iofs, nisle, nippts
     &,                  map
     &,                  converged
     &,                  estimated_error
     &                   )
c
c=======================================================================
c     MOM 2 Relax using symmetric coefficients as input, but
c     normalizes them as in MOM 1.
c     Normalized coefficients are cfn2, cfs2, etc.
c     Uses parallelization trick to get Gauss/Seidel update
c=======================================================================
c
c                          O L D   R E L A X
c
c      solve:
c
c             A * dpsi = forc
c
c      for "dpsi" with dirichlet boundary conditions (dpsi=const on
c      each component of the boundary) by a "hypergrid" version of
c      Gauss-Seidel iteration.  In this version, the grid is
c      decomposed into 4 sets, each with the same values of
c      (i mod 2, j mod 2).  All calculations within a set may be
c      done in parallel.
c
c      inputs:
c              npt   = 5 or 9 (active coefficients)
c              variable = character string identifying solution variable
c              bc_symm = equatorial symmetry type (used only when the
c                        symmetry option is on. otherwise ignore it)
c              guess = initial approximation to solution
c              A     = linear operator (assumed symmetric)
c                      typically A is  grad{(1/h)*grad(dpsi)} -
c                      2dt*acor*{grad(f/h) x grad(dpsi)}
c                      using 5 or 9 pt discretizations
c              cf    = imt x jmt x 3 x 3 array of coefficients of A
c              sor   = over-relaxation multiplier
c              forc  = the sum of all terms evaluated at times tau
c                      or tau-1
c              epsilon = convergence criterion
c              max_iterations = maximum number of iterations
c              imask = shows which land masses have perimeter equations
c              iperm = i coordinate of island perimeter points
c              jperm = j coordinate of island perimeter points
c              iofs  = offset in iperm, jperm for start of perimeter
c                      of land_mass(isle)
c              nisle = actual number of land_masses
c              nippts = number of perimeter ocean points for a land_mass
c      output:
c              dpsi   = answer
c              iterations = actual number of iterations performed
c              converged = logical value
c              estimated_error = estimated maximum error in solution
c                          based on step sizes and convergence rate
c
c=======================================================================
c
c      more specifically, the equations to be solved are
c
c             sum (A(ij,i'j') * dpsi(i'j')) = forc(ij)
c
c      where the subscripts ij and i'j' range over all "free ocean"
c      T cells ij=(i,j) that are not adjacent to land T cells,
c      and one ij=isle for each boundary component of the ocean.
c
c      with this choice of variables, in the absence of coriolis terms
c      (acor=0), the operator A is symmetric, i.e.,
c
c             A(ij,i'j') = A(i'j',ij)
c
c=======================================================================
c
c     previous versions by:     m. cox, b. semtner, r.c. pacanowski
c     author:              c.h. goldberg    e-mail => chg@gfdl.gov
c
c=======================================================================
c
c     dimensions of local arrays
#include "size.h"
c
      logical imask(-mnisle:mnisle)
      dimension dpsi(imt,jmt), forc(imt,jmt), res(imt,jmt)
      dimension cf(imt,jmt,-1:1,-1:1)
      dimension relmsk(imt,jmt), guess(imt,jmt)
      dimension nippts(mnisle)
      dimension iofs(mnisle), iperm(maxipp), jperm(maxipp)
      dimension map(imt,jmt)
      dimension rncfdiag(imt,jmt)
      dimension cfn2(imt,jmt)
      dimension cfs2(imt,jmt)
      dimension cfe2(imt,jmt)
      dimension cfw2(imt,jmt)
      dimension forc2(imt,jmt)
      dimension diagsum(mnisle)
      logical converged
      character * 16 variable
      character*(*) bc_symm
c
#ifdef debug_relax1
c
c-----------------------------------------------------------------------
c     verify that domain(forc) is a subset of domain(A)
c-----------------------------------------------------------------------
c
      call subset (forc, cf(1,1,0,0), nerror)
      if (nerror .ne. 0) then
        print '(a,a,a,i5,a,a)'
     &,   'fatal error: ',variable(1:len_trim(variable))
     &,   ' forcing is nonzero at ',nerror,' points '
     &,   ' with no equation, i.e., cf(i,j,0,0)=0.'
c       stop '=>relax1'
      end if
#endif
c
c-----------------------------------------------------------------------
c     the parallelization tricks used in relax1 work only for 5 pt
c     operators.  do not use relax1 with 9 point operators.
c-----------------------------------------------------------------------
c
      if (npt .ne. 5) then
        print '(a)', 'WARNING:  relax1 works only with 5 pt operators'
        mscan = 0
        converged = .false.
        stop '=>relax1'
      end if
c
c-----------------------------------------------------------------------
c     set locally needed constants
c-----------------------------------------------------------------------
c
      c0    = 0.0
      c1    = 1.0
c
c-----------------------------------------------------------------------
c     "normalize" coefficients for "oldrelax" method as in MOM1
c     relmsk is now a locally computed array
c     it is 1 on mid-ocean points, and 0 elsewhere
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do i=1,imt
          if (map(i,j) .eq. 0) then
            relmsk(i,j) = c1
          else
            relmsk(i,j) = c0
          end if
        end do
      end do
c
c-----------------------------------------------------------------------
c     initialize arrays
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do i=1,imt
          cfn2(i,j)=c0
          cfs2(i,j)=c0
          cfe2(i,j)=c0
          cfw2(i,j)=c0
          rncfdiag(i,j) = c1
        end do
      end do
c
      do isle=1,nisle
        diagsum(isle) = c0
      end do
c
      do j=2,jmt-1
        do i=2,imt-1
          if (map(i,j) .eq. 0) then
            rncfdiag(i,j) =
     &         c1/(cf(i,j,0,1)+cf(i,j,0,-1)+cf(i,j,1,0)+cf(i,j,-1,0))
c
c           normalize coefficients (mid ocean)
c
            cfn2(i,j) = cf(i,j, 0, 1)*rncfdiag(i,j)
            cfs2(i,j) = cf(i,j, 0,-1)*rncfdiag(i,j)
            cfe2(i,j) = cf(i,j, 1, 0)*rncfdiag(i,j)
            cfw2(i,j) = cf(i,j,-1, 0)*rncfdiag(i,j)
          end if
c
c         sum diagonal coefficients on island boundary
c
          if (map(i,j) .le. -1) then
            isle = -map(i,j)
            if (imask(isle)) then
              diagsum(isle) = diagsum(isle)+cf(i,j,0,0)
            end if
          end if
        end do
      end do
c
c-----------------------------------------------------------------------
c     normalize coefficients on island boundaries
c-----------------------------------------------------------------------
c
      do isle=1,nisle
        if (imask(isle)) then
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            rncfdiag(i,j) = -c1/diagsum(isle)
c
c           normalize coefficients (island boundary)
c
              cfn2(i,j) = cf(i,j, 0, 1)*rncfdiag(i,j)
              cfs2(i,j) = cf(i,j, 0,-1)*rncfdiag(i,j)
              cfe2(i,j) = cf(i,j, 1, 0)*rncfdiag(i,j)
              cfw2(i,j) = cf(i,j,-1, 0)*rncfdiag(i,j)
          end do
        end if
      end do
c
c-----------------------------------------------------------------------
c     pre-multiply all coefficients by sor
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do i=1,imt
          cfn2(i,j) = cfn2(i,j)*sor
          cfs2(i,j) = cfs2(i,j)*sor
          cfe2(i,j) = cfe2(i,j)*sor
          cfw2(i,j) = cfw2(i,j)*sor
        end do
      end do
c
#ifdef symmetry
c
c-----------------------------------------------------------------------
c     check for allowable symmetry conditions at equator
c-----------------------------------------------------------------------
c
      if (bc_symm .ne. 't even' .and. bc_symm .ne. 't odd') then
        write (*,*) 
     &'=>Error: type of bc_symm is ',bc_symm,' in hyper3.F'
     &, '... Only allowed type is "t even" or "t odd"'
        stop
      endif
#endif
c
c-----------------------------------------------------------------------
c     impose boundary conditions on guess
c     dpsi(0) = guess
c-----------------------------------------------------------------------
c
      call border(guess, bc_symm)
c
c-----------------------------------------------------------------------
c     set residuals to zero and normalize forcing
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do i=1,imt
          res(i,j)  = c0
          forc2(i,j) = forc(i,j)*rncfdiag(i,j)
          dpsi(i,j) = guess(i,j)
        end do
      end do
c
c-----------------------------------------------------------------------
c     begin iteration loop
c-----------------------------------------------------------------------
c
      do mscan=1,mxscan
c
c-----------------------------------------------------------------------
c       compute residuals without using updated "dpsi" values to get
c       vector of maximum length
c-----------------------------------------------------------------------
c
        do j=2,jmt-1
          do i=2,imt-1
          res(i,j) = (cfn2(i,j)*dpsi(i,j+1) +
     &                cfs2(i,j)*dpsi(i,j-1) +
     &                cfe2(i,j)*dpsi(i+1,j) +
     &                cfw2(i,j)*dpsi(i-1,j) -
     &                sor*(dpsi(i,j)+forc2(i,j)))*relmsk(i,j)
          end do
        end do
c
        call border(res, bc_symm)
c
c-----------------------------------------------------------------------
c       correct southern point using updated "dpsi" to get vectors on "i"
c-----------------------------------------------------------------------
c
        do j=2,jmt-1
          do i=2,imt-1
            res(i,j) = res(i,j) + cfs2(i,j)*res(i,j-1)*relmsk(i,j)
          end do
c
c---------------------------------------------------------------------
c       correct western point using updated "dpsi" to get vectors on "j"
c---------------------------------------------------------------------
c
          do i=2,imt-1
            res(i,j) = res(i,j) + cfw2(i,j)*res(i-1,j)*relmsk(i,j)
          end do
        end do
c
        call border(res, bc_symm)
c
c---------------------------------------------------------------------
c       make a correction to dpsi based on the residuals
c---------------------------------------------------------------------
c
        do j=2,jmt-1
          do i=1,imt
            res(i,j)  = res(i,j)*relmsk(i,j)
            dpsi(i,j) = dpsi(i,j) + res(i,j)
          end do
        end do
c
c---------------------------------------------------------------------
c       find the maximum absolute residual to determine convergence
c---------------------------------------------------------------------
c
        resmax = absmax(res)
c
c-----------------------------------------------------------------------
c       do a line integral around each island
c---------------------------------------------------------------------
c
        do isle=1,nisle
          if (imask(isle)) then
            resis = c0
            do n=1,nippts(isle)
              i = iperm(iofs(isle)+n)
              j = jperm(iofs(isle)+n)
              resis = resis +  cfn2(i,j)*dpsi(i  ,j+1)
     &                        +cfs2(i,j)*dpsi(i  ,j-1)
     &                        +cfe2(i,j)*dpsi(i+1,j  )
     &                        +cfw2(i,j)*dpsi(i-1,j  )
     &                        -sor*(          forc2(i,j))
            end do
            resis = resis - sor*dpsi(i,j)
c
            resmax = max(abs(resis),resmax)
c
            do n=1,nippts(isle)
              i = iperm(iofs(isle)+n)
              j = jperm(iofs(isle)+n)
              dpsi(i,j) = dpsi(i,j) + resis
            end do
          end if
        end do
c
        call border(dpsi, bc_symm)
c
c-----------------------------------------------------------------------
c       test for convergence of the relaxation.
c-----------------------------------------------------------------------
c
        step = resmax
c
c-----------------------------------------------------------------------
c       the solver is deemed to have converged when the estimated
c       maximum sum of all future corrections does not exceed
c       crit at any point.
c-----------------------------------------------------------------------
c
        if (mscan .eq. 1) then
          step1 = step
          estimated_error = step
          if (step .lt. crit) goto 1001
        else if (step .lt. crit) then
          cfactor = log(step/step1)
          convergence_rate = exp(cfactor/(mscan-1))
          estimated_error = step*convergence_rate/(1.0-convergence_rate)
#ifdef debug_relax1
          print '(a,i6,4(a,e13.4))'
     &,      'iteration=', mscan, ', correction=', step
     &,      ', convergence rate=', convergence_rate
     &,      ', estimated error=', estimated_error
#endif
          if (estimated_error  .lt. crit)  goto 1001
        end if
      end do
c
c---------------------------------------------------------------------
c     end of iteration loop
c---------------------------------------------------------------------
c
1001  continue
      if (mscan .lt. mxscan) then
        converged = .true.
      else
         converged = .false.
      end if
c
c---------------------------------------------------------------------
c     return the last increment to dpsi in the argument res
c-----------------------------------------------------------------------
c
      do i=1,imt
        do j=1,jmt
          res(i,j) = res(i,j)
        end do
      end do
c
      return
      end
