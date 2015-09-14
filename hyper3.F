      subroutine hyper3 (npt, variable, bc_symm 
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
c     MOM 2 "hypergrid" relax using symmetric coefficients
c     It does not normalize symmetric coefficients ala MOM 1.
c     Hyper3 does checkerboard updating as 4 loops of constant stride
c
c=======================================================================
c
c                          H Y P E R G R I D
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
#include "size.h"
c
      logical imask(-mnisle:mnisle)
      dimension dpsi(imt,jmt), forc(imt,jmt), res(imt,jmt)
      dimension cf(imt,jmt,-1:1,-1:1)
      dimension relmsk(imt,jmt), guess(imt,jmt)
      dimension nippts(mnisle)
      dimension iofs(mnisle), iperm(maxipp), jperm(maxipp)
      dimension map(imt,jmt)
      dimension rcfdiag(imt,jmt)
      dimension diagsum(mnisle)
      dimension resmi(jmt)
      logical converged
      character * 16 variable
      character*(*) bc_symm
c
#ifdef debug_hyper3
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
c       stop '=>hyper3'
      end if
#endif
c
c-----------------------------------------------------------------------
c     set locally needed constants
c-----------------------------------------------------------------------
c
      c0    = 0.0
      c1    = 1.0
c
c-----------------------------------------------------------------------
c     calculate "normalized" coefficients used in MOM1
c
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
      do isle=1,nisle
        diagsum(isle) = c0
      end do
c
      do j=2,jmt-1
        do i=2,imt-1
	  if (cf(i,j,0,0) .eq. 0.0) then
	    rcfdiag(i,j) = c0
          else if (map(i,j) .eq. 0) then
            rcfdiag(i,j) = c1/cf(i,j,0,0)
          else
            rcfdiag(i,j) = c0
          end if
c
c         sum diagonal coefficients on island boundary
c
          isle = -map(i,j)
          if (isle .gt. 0 .and. imask(isle)) then
            diagsum(isle) = diagsum(isle)+cf(i,j,0,0)
          end if
        end do
      end do
c
      do isle=1,nisle
        if (imask(isle)) then
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            rcfdiag(i,j) = c1/diagsum(isle)
          end do
        end if
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
c     set residuals to zero and initialize dpsi
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do i=1,imt
          res(i,j)  = c0
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
        do j=2,jmt-1
          resmi(j) = c0
        end do
c
c-----------------------------------------------------------------------
c       consider the arrays as being defined on the squares of an
c       "imt by jmt" checkerboard. take four passes: first solve the
c       equation on the black squares in even columns, then on the red 
c       squares in even columns, then red squares in odd columns, and
c       finally on black squares in odd columns..
c-----------------------------------------------------------------------
c
        if (npt .eq. 5) then
c
c         5 point calculation
c
          do i1=0,1
            do j1=0,1
              do j=2+j1,jmt-1,2
                do i=2+i1,imt-1,2
                  res(i,j) =  relmsk(i,j) * 
     &                      ((forc(i,j)
     &                       -cf(i,j, 0, 1)*dpsi(i,j+1)
     &                       -cf(i,j, 0,-1)*dpsi(i,j-1)
     &                       -cf(i,j, 1, 0)*dpsi(i+1,j)
     &                       -cf(i,j,-1, 0)*dpsi(i-1,j)
     &                       )*rcfdiag(i,j) - dpsi(i,j) )
                end do
c
                call border(res, bc_symm)
c
c               make a correction to dpsi based on the residuals
c
                do i=2+i1,imt,2
                  dpsi(i,j) = dpsi(i,j) + sor * res(i,j)
                end do
c
c               find the maximum absolute residual to determine convergence
c
                do i=2+i1,imt,2
                  resmi(j) = max(abs(res(i,j)),resmi(j))
                end do
              end do
            end do
            call border(res, bc_symm)
          end do
        else
c
c         9 point calculation
c
          do i1=0,1
            do j1=0,1
              do j=2+j1,jmt-1,2
                do i=2+i1,imt-1,2
                  res(i,j) =  relmsk(i,j) * 
     &                      ((forc(i,j)
     &                       -cf(i,j, 0, 1)*dpsi(i,j+1)
     &                       -cf(i,j, 0,-1)*dpsi(i,j-1)
     &                       -cf(i,j, 1, 0)*dpsi(i+1,j)
     &                       -cf(i,j,-1, 0)*dpsi(i-1,j)
     &                       -cf(i,j, 1, 1)*dpsi(i+1,j+1)
     &                       -cf(i,j,-1, 1)*dpsi(i-1,j+1)
     &                       -cf(i,j, 1,-1)*dpsi(i+1,j-1)
     &                       -cf(i,j,-1,-1)*dpsi(i-1,j-1)
     &                       )*rcfdiag(i,j) - dpsi(i,j) )
                end do
c
                call border(res, bc_symm)
c
c               make a correction to dpsi based on the residuals
c
                do i=2+i1,imt,2
                  dpsi(i,j) = dpsi(i,j) + sor * res(i,j)
                end do
c
c               find the maximum absolute residual to determine convergence
c
                do i=2+i1,imt,2
                  resmi(j) = max(abs(res(i,j)),resmi(j))
                end do
              end do
            end do
            call border(res, bc_symm)
          end do
        end if
c
c-----------------------------------------------------------------------
c       find maximum residual
c-----------------------------------------------------------------------
c
        resmax = c0
        do j=2,jmt-1
          resmax = max(resmi(j),resmax)
        end do
c
c-----------------------------------------------------------------------
c       do integration around each island
c-----------------------------------------------------------------------
c
        do isle=1,nisle
          if (imask(isle)) then
            resis = c0
            do n=1,nippts(isle)
              i = iperm(iofs(isle)+n)
              j = jperm(iofs(isle)+n)
              resis = resis +  forc(i,j)
     &                      -  cf(i,j, 0, 1)*dpsi(i,j+1)
     &                      -  cf(i,j, 0,-1)*dpsi(i,j-1)
     &                      -  cf(i,j, 1, 0)*dpsi(i+1,j)
     &                      -  cf(i,j,-1, 0)*dpsi(i-1,j)
     &                      -  cf(i,j, 1, 1)*dpsi(i+1,j+1)
     &                      -  cf(i,j,-1, 1)*dpsi(i-1,j+1)
     &                      -  cf(i,j, 1,-1)*dpsi(i+1,j-1)
     &                      -  cf(i,j,-1,-1)*dpsi(i-1,j-1)
            end do
            resis =  resis / diagsum(isle) - dpsi(i,j)
c 
            resmax = max(abs(resis),resmax)
c
            do n=1,nippts(isle)
              i = iperm(iofs(isle)+n)
              j = jperm(iofs(isle)+n)
              dpsi(i,j) = dpsi(i,j) + sor * resis
            end do
          end if
        end do
c
        call border(dpsi, bc_symm)
c
c-----------------------------------------------------------------------
c       test for convergence of the relaxation.
c       the solver is deemed to have converged when the estimated
c       maximum sum of all future corrections does not exceed
c       crit at any point.
c-----------------------------------------------------------------------
c
        step = sor * resmax
        if (mscan .eq. 1) then
          step1 = step
	  estimated_error = step
          if (step .lt. crit) goto 1001
        else if (step .lt. crit) then
          cfactor = log(step/step1)
          convergence_rate = exp(cfactor/(mscan-1))
	  estimated_error = step*convergence_rate/(1.0-convergence_rate)
#ifdef debug_hyper3
          print '(a,i6,4(a,e13.4))'
     &,      'iteration=', mscan, ', correction=', step
     &,      ', convergence rate=', convergence_rate
     &,      ', estimated error=', estimated_error
#endif
          if (estimated_error .lt. crit)  goto 1001
        end if
      end do
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
          res(i,j) = sor * res(i,j)
        end do
      end do
c
      return
      end
