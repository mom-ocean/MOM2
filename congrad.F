      subroutine congr  (npt, variable, bc_symm
     &,                  guess, dpsi, forc, res
     &,                  cf
     &,                  max_iterations, iterations, epsilon
     &,                  imask, iperm, jperm, iofs, nisle, nippts
     &,                  converged
     &,                  estimated_error
     &                  )
c
c=======================================================================
c
c                            C O N G R A D
c
c      solve:
c
c             A * dpsi = forc
c
c      for "dpsi" with dirichlet boundary conditions (dpsi=const on
c      each component of the boundary) by a preconditioned conjugate
c      gradient algorithm.
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
c
c      author:      Charles Goldberg        e-mail: chg@gfdl.gov
c
c      based on the preconditioned conjugate gradient algorithm given
c      in:
c
c     A Reformulation and Implementation of the Bryan-Cox-Semtner
c     Ocean Model on the Connection Machine
c     J.K. Dukowicz, R.D. Smith, and R.C. Malone
c     Journal of Atmospheric and Oceanic Technology
c     Vol 10. No. 2 April 1993
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
c
c      the algorithm (essentially executable in Fortran 90) is...
c
c      subroutine congrad (A, guess, forc, dpsi, iterations)
c
c      use matrix_module
c
c      intent (in)     :: A, guess, forc
c      intent (out)    :: dpsi, iterations
c
c      type(dpsi_type) :: guess, dpsi, Zres, s
c      type(res_type)  :: res, As, forc
c      type(operator)  :: A
c      type(inv_op)    :: Z
c      dimension (0:max_iterations) :: dpsi, res, s, As, beta, alpha
c
c      dpsi(0) = guess
c      res(0)  = forc - A * dpsi(0)
c      beta(0) = 1
c      s(0)    = zerovector()
c      do k = 1 to max_iterations
c        Zres(k-1) = Z * res(k-1)
c        beta(k)   = res(k-1) * Zres(k-1)
c        s(k)      = Zres(k-1) + (beta(k)/beta(k-1)) * s(k-1)
c        As(k)     = A * s(k)
c        alpha(k)  = beta(k) / (s(k) * As(k))
c        dpsi(k)   = dpsi(k-1) + alpha(k) * s(k)
c        res(k)    = res(k-1) - alpha(k) * As(k)
c        estimated_error = err_est(k, alpha(k), s(k))
c        if (estimated_error) < epsilon) exit
c      end do
c      if (k > max_iterations) then
c        print *, 'did not converge in ',k,' iterations'
c        stop '=>congrad'
c      end if
c
c      iterations = k
c      dpsi = dpsi(k)
c
c      end
c
c      where...
c
c      the "vector" and "operator" types used in conjugate gradient
c      are mapped to ordinary 2-dimensional fortran arrays as follows:
c
c      type(dpsi_type) :: guess, dpsi, Zres, s
c          if ij=(i,j) is a mid-ocean point, map dpsi(ij)-->dpsi(i,j)
c          if ij=isle is an ocean boundary subscript, replicate the
c          value dpsi(isle) in dpsi(i,j) for each (i,j) in the ocean
c          perimeter of land_mass(isle).  the arrays iperm(isle) and
c          jperm(isle), along with iofs(isle) locate these ocean
c          perimeter T cells.
c      type(res_type)  :: res, As, forc
c          if ij=(i,j) is a mid-ocean point, res(ij)-->res(i,j)
c          if ij=isle is an ocean boundary subscript, the value of
c          res(isle) = sum (res(i,j))
c          where the sum is taken over all (i,j) in the ocean perimeter
c          of land_mass(isle).  sometimes, the computed values
c          res(i,j) represent contributions of T cell (i,j) to the
c          component res(isle), and sometimes the values are balanced
c          so that res(i,j)=res(isle)/nippts(isle).  note that, even
c          when balanced, the relation between type(res_type) variables
c          res(isle) and res(i,j) differs from that of type(dpsi_type)
c          variables dpsi(isle) and dpsi(i,j) on T cells in the ocean
c          perimeter.
c      type(operator)  :: A
c          the nearly diagonal quality of the operators used
c          permits a representation as a small collection of
c          2-dimensional arrays.
c          the diagonal, A(ij,ij), is stored in an array cfdiag(i,j)
c          as follows:
c          if ij=(i,j) is a mid-ocean point, A(ij,ij) = cfdiag(i,j)
c          if ij=isle is an ocean boundary subscript,
c          A(isle, isle) = sum (cfdiag(i,j))
c          where the sum is taken over all (i,j) in the ocean perimeter
c          of land_mass(isle).  each cfdiag(i,j) represents the contribution
c          of T cell (i,j) to the island variable diagonal coefficient.
c          the off-diagonal terms A(ij,i`j`) are stored in 4 arrays
c          cfn, cfs, cfe, and cfw if A is a 5-point operator, and in
c          these and 4 additional arrays, cfne, cfnw, cfse, cfsw, if
c          A is a 9-point operator.  For example, if i`=i and j`=j+1,
c          then A(ij,i`j`) is stored in cfn(i,j).
c          if ij=(i,j) is a mid-ocean point and i`j`=isle` is and ocean
c          perimeter subscript, with i`=i and j`=j+1, then
c          cfn(i,j)=A(ij,isle`) is the coefficient of the island
c          variable dpsi(isle`) in the equation for mid-ocean point
c          dpsi(ij)=dpsi(i,j).
c          if ij=isle is an ocean perimeter point and i`j`=(i`,j`) is
c          a mid-ocean point, with i`=i and j`=j-1, then
c          cfs(i,j)=A(isle,i`j`) is the coefficient of the mid-ocean
c          variable dpsi(i`j`)=dpsi(i,j) in the equation for the island
c          variable dpsi(isle).  note that equations for island
c          variables dpsi(isle) are "non-local" in the sense that
c          they usually contain more than 5 or 9 terms, some of which
c          involve values dpsi(i`j`) outside of a compact 5-point
c          or 9-point neighborhood.
c      type(inv_op)    :: Z
c          the approximate inverse operator Z used at present is a
c          diagonal operator Z(ij,ij) = 1/A(ij,ij).
c          if ij=(i,j) is a mid-ocean point,
c            then Z(i,j)=Z(ij)=1/A(ij)=1/cfdiag(i,j)
c          if ij=isle is an ocean perimeter point, then
c          Z(isle) is replicated at each ocean perimeter T cell
c          bordering land_mass(isle).
c            Z(i,j)=Z(isle)=1/A(isle)=1/sum(A(i,j))
c
c=======================================================================
c
#include "size.h"
c
      logical imask, converged
      character * 16 variable
      character *(*) bc_symm
      dimension guess(imt,jmt), dpsi(imt,jmt), Zres(imt,jmt)
      dimension s(imt,jmt)
      dimension res(imt,jmt), As(imt,jmt), forc(imt,jmt)
      dimension cf(imt,jmt,-1:1,-1:1)
      dimension Z(imt,jmt)
      dimension imask(-mnisle:mnisle)
      dimension iperm(maxipp), jperm(maxipp), iofs(mnisle)
      dimension nippts(mnisle)
c
#ifdef debug_congrad
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
c       stop '=>congr'
      end if
#endif
c
#ifdef symmetry
c
c-----------------------------------------------------------------------
c     check for allowable symmetry conditions at equator
c-----------------------------------------------------------------------
c
      if (bc_symm .ne. 't even' .and. bc_symm .ne. 't odd') then
        write (*,*) 
     &'=>Error: type of bc_symm is ',bc_symm,' in congrad.F'
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
      do i=1,imt
        do j=1,jmt
          dpsi(i,j) = guess(i,j)
        end do
      end do
c
c-----------------------------------------------------------------------
c     make approximate inverse operator Z (always even symmetry)
c-----------------------------------------------------------------------
c
      call make_inv (cf, Z,
     &               imask, iperm, jperm, iofs, nisle, nippts)
      call border(Z, 't even')
c
c-----------------------------------------------------------------------
c     res(0)  = forc - A * dpsi(0)
c     impose cyclic and/or symmetry conditions on res(i,j)
c-----------------------------------------------------------------------
c
      if (npt .eq. 5) then
        call op5_vec(cf, dpsi, res)
      else
        call op9_vec(cf, dpsi, res)
      end if
      do i=2,imt-1
        do j=2,jmt-1
          res(i,j) = forc(i,j) - res(i,j)
        end do
      end do
c
      call border(res, bc_symm)
c
c-----------------------------------------------------------------------
c     Zres(k-1) = Z * res(k-1)
c     see if guess is a solution, bail out to avoid division by zero
c-----------------------------------------------------------------------
c
        k = 0
        call inv_op(Z, res, Zres,
     &              imask, iperm, jperm, iofs, nisle, nippts)
c
c       set borders of Zres using cyclic/symmetry, if defined.
c
        call border(Zres, bc_symm)
c
#ifdef debug_congrad
        Zresmax = absvecmax(Zres,imax,jmax)
        print '(a,a,i5,a,e14.7,a,i3,a,i3,a)'
     &,      variable(1:len_trim(variable))
     &,     ' iteration',k,'   max Zres    = ',Zresmax
     &,          ' at (',imax,',',jmax,')'
#else
        Zresmax = absmax(Zres)
#endif
c
c       Assume convergence rate of 0.99 to extrapolate error
c
        if (100.0 * Zresmax .lt. epsilon) then
	  estimated_error = 100.0 * Zresmax 
	  goto 101
	endif
c
c-----------------------------------------------------------------------
c     beta(0) = 1
c     s(0)    = zerovector()
c-----------------------------------------------------------------------
c
      betakm1 = 1.0
      call zero_vec(s)
c
c-----------------------------------------------------------------------
c     begin iteration loop
c-----------------------------------------------------------------------
c
      do k = 1,max_iterations
c
c-----------------------------------------------------------------------
c       Zres(k-1) = Z * res(k-1)
c-----------------------------------------------------------------------
c
        call inv_op(Z, res, Zres,
     &              imask, iperm, jperm, iofs, nisle, nippts)
c
c       set borders of Zres using cyclic/symmetry, if defined.
c
        call border(Zres, bc_symm)
c
c-----------------------------------------------------------------------
c       beta(k)   = res(k-1) * Zres(k-1)
c-----------------------------------------------------------------------
c
        betak = dot2(Zres, res)
        if (k .eq. 1) then
          betak_min = abs(betak)
	elseif (k .gt. 2) then
          betak_min = min(betak_min, abs(betak))
          if (abs(betak) .gt. 100.0*betak_min) then
	    write (*,'(/a/a/)')
     &      'WARNING: conjugate gradient terminated because correction'
     &,     '         steps are diverging. Probable cause...roundoff'
            smax = absmax(s)
            step = abs(alpha) * smax
	    estimated_error=step*convergence_rate/(1.0-convergence_rate)
	    go to 101
	  endif
	endif
c
c-----------------------------------------------------------------------
c       s(k)      = Zres(k-1) + (beta(k)/beta(k-1)) * s(k-1)
c-----------------------------------------------------------------------
c
        betaquot = betak/betakm1
        do i=1,imt
          do j=1,jmt
            s(i,j) = Zres(i,j) + betaquot * s(i,j)
          end do
        end do
c
c-----------------------------------------------------------------------
c       As(k)     = A * s(k)
c-----------------------------------------------------------------------
c
        if (npt .eq. 5) then
          call op5_vec(cf, s, As)
        else
          call op9_vec(cf, s, As)
        end if
c
        call border(As, bc_symm)
c
c-----------------------------------------------------------------------
c       If s=0 then the division for alpha(k) gives a float exception.
c       Assume convergence rate of 0.99 to extrapolate error.
c       Also assume alpha(k) ~ 1.
c-----------------------------------------------------------------------
c
        s_dot_As = dot2(s, As)
        if (abs(s_dot_As) .lt. abs(betak)*1.e-10) then
          smax = absmax(s)
#ifdef debug_congrad
          print *, 'Congrad exited on smax=',smax
#endif
	  estimated_error = 100.0 * smax 
	  goto 101
	endif
c
c-----------------------------------------------------------------------
c       alpha(k)  = beta(k) / (s(k) * As(k))
c-----------------------------------------------------------------------
c
        alpha = betak / s_dot_As
#ifdef debug_congrad
        print *, 'alpha=',alpha,' beta(k)=',betak
#endif
c
c-----------------------------------------------------------------------
c       update values:
c       dpsi(k)   = dpsi(k-1) + alpha(k) * s(k)
c       res(k)    = res(k-1) - alpha(k) * As(k)
c-----------------------------------------------------------------------
c
        do i=1,imt
          do j=1,jmt
            dpsi (i,j) = dpsi(i,j) + alpha * s(i,j)
            res  (i,j) = res (i,j) - alpha * As(i,j)
          end do
        end do
        call avg_dist (res,
     &          imask, iperm, jperm, iofs, nisle, nippts)
c
        call border(res, bc_symm)
c
#ifdef debug_congrad
        smax = absvecmax(s, imax, jmax)
        print '(a,i5,a,e14.7,a,i3,a,i3,a)',
     &           'iteration',k,'   max correct = ',alpha*smax,
     &           ' at (',imax,',',jmax,')'
#else
        smax = absmax(s)
#endif
c
c-----------------------------------------------------------------------
c       test for convergence
c       if (estimated_error) < epsilon) exit
c-----------------------------------------------------------------------
c
        step = abs(alpha) * smax
        if (k .eq. 1) then
          step1 = step
          estimated_error = step
          if (step .lt. epsilon) goto 101
        else if (step .lt. epsilon) then
          cfactor = log(step/step1)
          convergence_rate = exp(cfactor/(k-1))
	  estimated_error = step*convergence_rate/(1.0-convergence_rate)
#ifdef debug_congrad
          print '(a,i6,4(a,e13.4))'
     &,      'iteration=', k, ', correction=', step
     &,      ', convergence rate=', convergence_rate
     &,      ', estimated error=', estimated_error
#endif
          if (estimated_error .lt. epsilon) goto 101
        end if
c
        betakm1 = betak
c
      end do
c
c-----------------------------------------------------------------------
c     end of iteration loop
c-----------------------------------------------------------------------
c
  101 continue
      if (k .gt. max_iterations) then
        cfactor = log(step/step1)
        convergence_rate = exp(cfactor/(k-1))
	estimated_error = step*convergence_rate/(1.0-convergence_rate)
#ifdef debug_congrad
        print '(a,i4,a,a,e10.3,a,e10.3)',
     & '****> congrad did not converge in ',k,' iterations.'
     &,' maximum dpsi is ',absmax(dpsi)
     &,' estimated error is ',estimated_error
#endif
        converged = .false.
      else
        converged = .true.
      end if

      iterations = k
c
c-----------------------------------------------------------------------
c     return the last increment of dpsi in the argument res
c-----------------------------------------------------------------------
c
      if (iterations .eq. 0) then
        do i=1,imt
          do j=1,jmt
            res(i,j) = Zres(i,j)
          end do
        end do
      else
        do i=1,imt
          do j=1,jmt
            res(i,j) = alpha * s(i,j)
          end do
        end do
      endif
#ifdef debug_congrad
c
c     print out convergence information each call
c
        print '(a,l1,a,i4,a,a,a,a,e10.3,a,e10.3)',
     & '      congrad converged = ',converged,' in ',k,' iterations.'
     &,' maximum ',variable(1:len_trim(variable)),' is ',absmax(dpsi)
     &,' estimated error is ',estimated_error
#endif
c
      return
      end



c
c=======================================================================
c
c     M A T R I X   M O D U L E   F O R   C O N G R A D
c
c=======================================================================
c
      subroutine zero_vec (v)
#include "size.h"
      dimension v(imt,jmt)
      do i=1,imt
        do j=1,jmt
          v(i,j) = 0.0
        end do
      end do
      return
      end


      subroutine add_vec (v,w,vpw)
#include "size.h"
      dimension v(imt,jmt), w(imt,jmt), vpw(imt,jmt)
      do i=1,imt
        do j=1,jmt
          vpw(i,j) = v(i,j) + w(i,j)
        end do
      end do
      return
      end


      subroutine sub_vec (v,w,vmw)
#include "size.h"
      dimension v(imt,jmt), w(imt,jmt), vmw(imt,jmt)
      do i=1,imt
        do j=1,jmt
          vmw(i,j) = v(i,j) - w(i,j)
        end do
      end do
      return
      end


      subroutine mult_vec(v,w,vtw)
#include "size.h"
      dimension v(imt,jmt), w(imt,jmt), vtw(imt,jmt)
      do i=1,imt
        do j=1,jmt
          vtw(i,j) = v(i,j) * w(i,j)
        end do
      end do
      return
      end


      subroutine div_vec(v,w,vdw)
#include "size.h"
      dimension v(imt,jmt), w(imt,jmt), vdw(imt,jmt)
      do i=1,imt
        do j=1,jmt
	  if (w(i,j) .ne. 0) then
            vdw(i,j) = v(i,j) / w(i,j)
	  else
	    vdw(i,j) = 0.0
	  end if
        end do
      end do
      return
      end


      subroutine scalar_vec (scalar,w,sw)
#include "size.h"
      dimension w(imt,jmt), sw(imt,jmt)
      do i=1,imt
        do j=1,jmt
          sw(i,j) = scalar * w(i,j)
        end do
      end do
      return
      end


      subroutine neg_vec (v)
#include "size.h"
      dimension v(imt,jmt)
      do i=1,imt
        do j=1,jmt
          v(i,j) = -v(i,j)
        end do
      end do
      return
      end



      function dot2 (dp_vec, res_vec)
c
c     this dot product produces the correct answers because for
c     ocean perimeter subscripts, ij=isle, the value on a
c     type(dpsi_type) vector, dp_vec(isle)=dp_vec(i,j), i.e., the true
c     value is replicated, and for a type(res_type) vector,
c     res_vec(isle) = sum (res_vec(i,j)), i.e., the true value is the
c     accumulation of the distributed values.
c
#include "size.h"
      dimension dp_vec(imt,jmt), res_vec(imt,jmt)
      dimension rowsum (jmt)
c
      do j=2,jmt-1
        rowsum(j) = 0.0
        do i=2,imt-1
          rowsum(j) = rowsum(j) + dp_vec(i,j) * res_vec(i,j)
        end do
      end do
c
      dot2 = 0.0
      do j=2,jmt-1
        dot2 = dot2 + rowsum(j)
      end do
c
      return
      end


      subroutine op5_vec(cf, dpsi, res)
c
c                       res = A * dpsi
c
c     this subroutine does not collect the terms of the true value
c     of res(isle) = sum (res(i,j)).  the contributions to the sum
c     remain distributed among the T cells (i,j) that form the
c     ocean perimeter of land_mass(isle).
c
c     at present, borders are not computed [i=1 or imt] [j=1 or jmt]
c
#include "size.h"
      dimension cf(imt,jmt,-1:1,-1:1)
      dimension dpsi(imt,jmt)
      dimension res(imt,jmt)
c
      do j=2,jmt-1
        do i=2,imt-1
          res(i,j) = cf(i,j, 0, 0) * dpsi(i,j)   +
     &               cf(i,j, 0, 1) * dpsi(i,j+1) +
     &               cf(i,j, 0,-1) * dpsi(i,j-1) +
     &               cf(i,j, 1, 0) * dpsi(i+1,j) +
     &               cf(i,j,-1, 0) * dpsi(i-1,j)
        end do
      end do
c
      return
      end


      subroutine op9_vec(cf, dpsi, res)
c
c                       res = A * dpsi
c
c     this subroutine does not collect the terms of the true value
c     of res(isle) = sum (res(i,j)).  the contributions to the sum
c     remain distributed among the T cells (i,j) that form the
c     ocean perimeter of land_mass(isle).
c
c     at present, borders are not computed [i=1 or imt] [j=1 or jmt]
c
#include "size.h"
      dimension cf(imt,jmt,-1:1,-1:1)
      dimension dpsi(imt,jmt)
      dimension res(imt,jmt)
c
      do j=2,jmt-1
        do i=2,imt-1
          res(i,j) = cf(i,j, 0, 0) * dpsi(i  ,j  ) +
     &               cf(i,j, 0, 1) * dpsi(i  ,j+1) +
     &               cf(i,j, 0,-1) * dpsi(i  ,j-1) +
     &               cf(i,j, 1, 0) * dpsi(i+1,j  ) +
     &               cf(i,j,-1, 0) * dpsi(i-1,j  ) +
     &               cf(i,j, 1, 1) * dpsi(i+1,j+1) +
     &               cf(i,j,-1, 1) * dpsi(i-1,j+1) +
     &               cf(i,j, 1,-1) * dpsi(i+1,j-1) +
     &               cf(i,j,-1,-1) * dpsi(i-1,j-1)
        end do
      end do
c
      return
      end


      subroutine subset (a, b, nerror)
c
c     verifies that the set of subscripts for which a(i,j) .ne. 0.0
c     is a subset of the set of subscripts for which b(i,j) .ne. 0.0
c
#include "size.h"
      dimension a(imt,jmt), b(imt,jmt)
c
      nerror = 0
      do i=2,imt-1
        do j=2,jmt-1
          if (a(i,j) .ne. 0.0 .and. b(i,j) .eq. 0.0) then
            nerror = nerror + 1
            print '(a,i3,a,i3,a,a)', '(',i,',',j,')'
     &                           ,' forcing is reset to zero'
c           set forcing (i.e., a(i,j)) to zero
            a(i,j) = 0.0
          end if
        end do
      end do
c
      return
      end


      subroutine inv_op(Z, res, Zres,
     &              imask, iperm, jperm, iofs, nisle, nippts)
c
c
c     apply and approximate inverse Z or the operator A
c
c     res is type(res_type), i.e., perimeter values res(isle)
c         are the sum of the distributed contributions res(i,j)
c     Zres is type(dpsi_type), i.e., perimeter values Zres(isle)
c         must be replicated at each perimeter point Zres(i,j)
c
c     borders  of Zres [i=1 or imt] [j=1 or jmt] must be defined
c     and must satisfy cyclic and/or symmetry, if defined.
c
c     currently, Z is diagonal:  Z(ij) = 1/A(ij)
c     and is stored in type(dpsi_type) format, i.e., Z(isle) is
c     replicated and stored in each Z(i,j) in the perimeter of
c     land_mass(isle).
c
#include "size.h"
c
      logical imask
      dimension imask (-mnisle:mnisle)
      dimension iperm(maxipp), jperm(maxipp)
      dimension iofs (mnisle)
      dimension nippts(mnisle)
c
      dimension Z(imt,jmt)
      dimension res(imt,jmt)
      dimension Zres(imt,jmt)
c
      do i=1,imt
        do j=1,jmt
          Zres(i,j) = Z(i,j) * res(i,j)
        end do
      end do
c
c     sum contributions to Zres(isle)
c     distribute Zres(isle) to all perimeter points
c
      call sum_dist (Zres,
     &        imask, iperm, jperm, iofs, nisle, nippts)
c
      return
      end


      function absvecmax(res, imax, jmax)
c
#include "size.h"
      dimension res(imt,jmt)
c
      absvecmax = 0.0
      do i=2,imt-1
        do j=2,jmt-1
          if (abs(res(i,j)) .gt. absvecmax) then
            absvecmax = abs(res(i,j))
            imax = i
            jmax = j
          end if
        end do
      end do
      return
      end


      function absmax (f)
c
#include "size.h"
      dimension f(imt,jmt)
c
      amax = 0.0
      do i=1,imt
        do j=1,jmt
	  amax = max(amax, abs(f(i,j)))
        end do
      end do
      absmax = amax
      return
      end


      function absmin (f)
c
#include "size.h"
      dimension f(imt,jmt)
c
      amin = 1.0e37
      do i=1,imt
        do j=1,jmt
          if (f(i,j) .ne. 0 .and. abs(f(i,j)) .lt. amin) then
            amin = abs(f(i,j))
          end if
        end do
      end do
      absmin = amin
      return
      end


      subroutine make_inv (cf, Z,
     &              imask, iperm, jperm, iofs, nisle, nippts)
c
c     construct an approximate inverse Z to A
c
c     Z will be diagonal:  Z(ij) = 1/A(ij)
c     and values for ocean perimeter entries Z(isle) will be replicated
c     at all T cells Z(i,j) in the ocean perimeter of land_mass(isle).
c
c     T cells (i,j) for which there is no diagonal coefficient
c     i.e., A(ij)=A(i,j)=0, are masked off by assigning Z(i,j)=0.
c     there are effectively no equations and no variables dpsi(i,j)
c     at these points.
c
#include "size.h"
      logical imask
      dimension imask (-mnisle:mnisle)
      dimension iperm(maxipp), jperm(maxipp)
      dimension iofs (mnisle)
      dimension nippts(mnisle)
c
      dimension cf(imt,jmt,-1:1,-1:1)
      dimension Z(imt,jmt)
c
c     copy diagonal coefficients of A to Z
c
      do i=2,imt-1
        do j=2,jmt-1
          Z(i,j) = cf(i,j,0,0)
        end do
      end do
c
c     for each land_mass(isle),
c     sum the contributions to cfdiag(isle)=A(isle,isle)
c     now stored in Z(i,j) at ocean perimeter T cells and replicate
c     the sum in all Z(i,j) for which (i,j) is in ocean perimeter
c     of land_mass(isle).
c
      call sum_dist (Z,
     &        imask, iperm, jperm, iofs, nisle, nippts)
c
c     now invert Z
c
      do i=2,imt-1
        do j=2,jmt-1
          if (Z(i,j) .ne. 0.0) then
            Z(i,j) = 1/Z(i,j)
          else
            Z(i,j) = 0.0
          end if
        end do
      end do
c
c     make inverse zero on island perimeters that are not integrated
c
      do isle=1,nisle
        if (.not. imask(isle)) then
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            Z(i,j) = 0.0
          end do
        end if
      end do
c
      return
      end



      subroutine sum_dist (Zres,
     &              imask, iperm, jperm, iofs, nisle, nippts)
c
c
c     sum contributions to Zres(isle)
c     distribute Zres(isle) to all perimeter points
c
c     this subroutine converts a type(res_type) vector with
c     distributed contributions to perimeter values
c        Zres(isle) = sum (Zres(i,j))
c     into a type (dpsi_type) vector with replicated values
c     for land_mass perimeters
c        Zres(isle) = Zres(i,j)
c     for all (i,j) in the ocean perimeter of land_mass(isle).
c
#include "size.h"
      logical imask
      dimension imask (-mnisle:mnisle)
      dimension iperm(maxipp), jperm(maxipp)
      dimension iofs (mnisle)
      dimension nippts(mnisle)
c
      dimension Zres(imt,jmt)
      dimension Zresisle(mnisle)
c
c     sum contributions to Zres(isle)
c
#ifdef timing
      call tic ('isle','loop1')
#endif
      do isle=1,nisle
        if (imask(isle)) then
          Zresisle(isle) = 0.0
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            Zresisle(isle) = Zresisle(isle) + Zres(i,j)
          end do
        end if
      end do
#ifdef timing
      call toc ('isle','loop1')
#endif
c
c     distribute Zres(isle) to all perimeter points
c
#ifdef timing
      call tic ('isle','loop2')
#endif
      do isle=1,nisle
        if (imask(isle)) then
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            Zres(i,j) = Zresisle(isle)
          end do
        end if
      end do
#ifdef timing
      call toc ('isle','loop2')
#endif
c
      return
      end


#ifdef debug_congrad
      subroutine examine_isles (cf,forc,dpsi,res,
     &              imask, iperm, jperm, iofs, nisle, nippts)
c
c     sum contributions to res(isle)
c     distribute res(isle) to all perimeter points
c
c     this subroutine is for debugging printouts of
c     island values.
c
#include "size.h"
      logical imask
      dimension imask (-mnisle:mnisle)
      dimension iperm(maxipp), jperm(maxipp)
      dimension iofs (mnisle)
      dimension nippts(mnisle)
c
      dimension res(imt,jmt)
      dimension resisle(mnisle)
c
      dimension cf(imt,jmt,-1:1,-1:1)
      dimension forc  (imt,jmt)
      dimension dpsi  (imt,jmt)
c
c     sum contributions to res(isle)
c
#ifdef timing
      call tic ('isle','loop5')
#endif
      do isle=1,nisle
        if (imask(isle)) then
          print *,'examine: isle=',isle,' nisle=',nisle
          resisle(isle) = 0.0
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            resisle(isle) = resisle(isle) + res(i,j)
            if (1 .ne. 1) then
            print '(2i3,a,e12.5)',i,j,' res  =',res(i,j)
            print '(2i3,a,e12.5)',i,j,' forc =',forc(i,j)
            print '(2i3,a,e12.5)',i,j,' dpsi =',dpsi(i,j)
            print '(2i3,a,e12.5)',i,j,' cfdiag=',cf(i,j,0,0)
            print '(2i3,a,e12.5)',i,j,' cfn  =',cf(i,j,0,1)
            print '(2i3,a,e12.5)',i,j,' cfs  =',cf(i,j,0,-1)
            print '(2i3,a,e12.5)',i,j,' cfe  =',cf(i,j,1,0)
            print '(2i3,a,e12.5)',i,j,' cfw  =',cf(i,j,-1,0)
            print *,'  '
            end if
          end do
          print *, 'resisle = ', resisle(isle)
          print *,' '
        end if
      end do
#ifdef timing
      call toc ('isle','loop5')
#endif
c
      return
      end
#endif



      subroutine avg_dist (Zres,
     &              imask, iperm, jperm, iofs, nisle, nippts)
c
c     avg contributions to Zres(isle)
c     distribute Zres(isle) to all perimeter points
c
c     this subroutine converts a type(res_type) vector with
c     distributed contributions to perimeter values
c        Zres(isle) = avg (Zres(i,j))
c     into a type (dpsi_type) vector with replicated values
c     for land_mass perimeters
c        Zres(isle) = Zres(i,j)
c     for all (i,j) in the ocean perimeter of land_mass(isle).
c
#include "size.h"
      logical imask
      dimension imask (-mnisle:mnisle)
      dimension iperm(maxipp), jperm(maxipp)
      dimension iofs (mnisle)
      dimension nippts(mnisle)
c
      dimension Zres(imt,jmt)
      dimension Zresisle(mnisle)
c
c     avg contributions to Zres(isle)
c
#ifdef timing
      call tic ('isle','loop3')
#endif
      do isle=1,nisle
        if (imask(isle)) then
c         print *,'isle=',isle,' nisle=',nisle
          Zresisle(isle) = 0.0
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            Zresisle(isle) = Zresisle(isle) + Zres(i,j)
          end do
        end if
      end do
#ifdef timing
      call toc ('isle','loop3')
#endif
c
c     distribute Zres(isle) to all perimeter points
c
#ifdef timing
      call tic ('isle','loop4')
#endif
      do isle=1,nisle
        if (imask(isle)) then
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            Zres(i,j) = Zresisle(isle)/nippts(isle)
          end do
        end if
      end do
#ifdef timing
      call toc ('isle','loop4')
#endif
c
      return
      end
