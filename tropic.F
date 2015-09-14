#ifndef test_poisson
      subroutine tropic (c2dtsf, acor, f, itt)
c
c=======================================================================
c
c        S O L V E   T H E   B A R O T R O P I C   E Q U A T I O N  
c
c   There are several significant changes made in MOM_2 in the
c   calculation of the vertically averaged velocities.
c
c   1. Based on a 1994 rederivation of the finite difference equations
c     for the stream function by Charles Goldberg (GFDL/Princeton
c     University/Trenton State College), the coefficients in
c     the Poisson equations for dpsi differ slighty from those
c     used in MOM_1.  Designed to reduce the residuals in the
c     finite difference momentum equations, these new equations
c     for dpsi seem to be less stiff.  Tests on a variety of 
c     geometries and topographies indicate that the new equations
c     converge more rapidly to the same tolerances in all solvers 
c     at all resolutions.
c
c   2. In all three Poisson solvers for stream function, the MOM_2
c     default mode is that the values of psi and dpsi are not held
c     fixed at zero on the boundary of the "main continent".  They
c     are allowed to float in the same way that other "island
c     boundary" values float.  Although this requires an island
c     integral to be computed on the boundary of every land mass, it
c     also makes the iterative system less stiff, and again results
c     in fewer iterations to converge to the same tolerances
c     in all solvers at all resolutions.
c
c     Tests indicate that except on architectures where island 
c     integrals are prohibitively expensive and the stream function
c     formulation therefore computationally infeasible, the cost
c     of the extra island integrals is small and the savings due
c     to reduced numbers of iterations significant. 
c
c     The user may specify either a land mass on whose boundary the
c     stream function is later normalized to zero or that no post
c     solver normalization is to take place.  Surface pressures are
c     always normalized to have mean zero.  Options are also provided 
c     for turning off the island integrals on the boundaries of
c     selected land masses; however, this practice is not recommended.
c
c   3. The convergence criterion in all three solvers has been
c     changed to the following:  "stop when the predicted maximum error
c     in the solved variable (dpsi or surface pressure) is less than
c     the user specified tolerance."  Convergence tolerances are now in
c     the same units as the variable solved for, and tell the user how
c     many digits of the answer are correct.Thus, if one expects a
c     maximum dpsi of 1.0e12 and desires convergence to 5 significants
c     digits, one chooses a tolerance of 1.0e7. Note that the tolerance
c     used in MOM 2 is NOT the same as "crit" in MOM 1.
c
c
c     The maximum error in the solution is predicted as follows:
c     First, convergence of the solver is assumed to be "geometric" in
c     the sense that the maximum absolute correction added to the
c     solved variable in iteration k is modeled as
c
c	   step(k) = step(1)*(convergence_rate)**(k-1)
c
c     The estimated maximum error in the solution after k iterations
c     is then bounded by the sum of the missing terms in the geometric
c     series truncated after k terms: 
c     
c     sum {step(i)} = step(k)*convergence_rate/(1.0 - convergence_rate)
c       i=k+1,infinity
c       
c		 
c     Experimental evidence indicates that the convergence rate 
c     of an iterative solver remains essentially stable over many
c     iterations, and that the maximum errors when the solvers are
c     stopped by this criterion as compared to solutions obtained
c     by allowing the solvers to run to machine precision are indeed
c     less than the stated tolerances.
c  
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
c
# include "param.h"
# include "emode.h"
# include "grdvar.h"
# include "iounit.h"
# include "mw.h"
# include "switch.h"
      character*8 bc_symm
      dimension f(jmt)
c#define chg_surface
# ifdef chg_surface
      dimension ubarp1(imt,jmt,2), psp1(imt,jmt)
# endif
c
# ifdef timing
      call tic ('ocean', 'tropic')
# endif
c
      call setbcx (zu(1,1,1), imt, jmt)
      call setbcx (zu(1,1,2), imt, jmt)
c
# if defined rigid_lid_surface_pressure || defined implicit_free_surface
c
c-----------------------------------------------------------------------
c     solve for the "tau+1" surface pressure and barotropic velocities
c-----------------------------------------------------------------------
c
#  ifdef chg_surface
c
c     for testing only ....
c
      call surface (ubarm1, ubar, ubarp1
     &,             ps(1,1,2), ps(1,1,1), psp1, zu 
     &,             c2dtsf/2.0, cori
     &,             map
     &,             imask, iperm, jperm, iofs, nislsp, nippts
     &,             res, cf, mxscan, mscan, h, tolrsp, tolrfs
     &,             divf)
c
c     shift quantities back 1 time step
c
      call copy_vec2 (ubar, ubarm1)
      call copy_vec2 (ubarp1, ubar)
      call copy_vec (ps(1,1,1), ps(1,1,2))
      call copy_vec (psp1, ps(1,1,1))
#  else
      call bardiv
#  endif
c
c-----------------------------------------------------------------------
c     test accuracy of solving for change in surface pressure
c-----------------------------------------------------------------------
c
      if (.not.converged) then
        write (stdout,'(a,i5,3(a,1pe10.3))')
     &  ' WARNING: SOLVER DID NOT CONVERGE in ',mscan
     &, ' scans. max(ps)='
     &, absmax(ps(1,1,2)), ' max(ptd)=',absmax(ptd)
     &, ' estimated max(err)=', esterr
      endif
# endif
c
# ifdef stream_function
c
c-----------------------------------------------------------------------
c     construct the forcing for the stream function equation
c-----------------------------------------------------------------------
c
      call sfforc (zu, dxu, dyu, csu, ztd) 

c
#  if defined fourfil || defined firfil
c
c     filter forcing at high latitudes
c
      call filz (ztd, cf)
#  endif
c
c-----------------------------------------------------------------------
c     read in solution from 2 previous timesteps for the
c     purpose of computing an initial guess for the present solution.
c-----------------------------------------------------------------------
c
      luptdb = nkflds - mod(itt,2)
      luptd  = nkflds - 1 + mod(itt,2)
      call oget (kflds, nwds, luptdb, res)
      call oget (kflds, nwds, luptd, ptd)
c
      fxa=c1
      if (.not. leapfrog) fxa=p5
      do jrow=1,jmt
        do i=1,imt
          ptd(i,jrow) = fxa*(c2*ptd(i,jrow)-res(i,jrow))
        enddo
      enddo
#  ifdef cyclic
      do jrow=2,jmtm1
        ptd(1,jrow)   = ptd(imtm1,jrow)
        ptd(imt,jrow) = ptd(2,jrow)
      enddo
#  endif
c
c-----------------------------------------------------------------------
c     choose 5 or 9 point numerics
c-----------------------------------------------------------------------
c
#  ifdef sf_5_point
c
c     initialize coefficients using 5 point numerics
c
      call sfc5pt (acor, f, c2dtsf, dxu, dyu, csu, hr, cf)
      npt = 5
#  endif
#  ifdef sf_9_point
c
c     initialize coefficients using 9 point numerics
c
      call sfc9pt (acor, f, c2dtsf, dxu, dyu, csu, hr, cf)
      npt = 9
#  endif
c
c-----------------------------------------------------------------------
c     choose a method for solving for the "tau+1" stream function change
c-----------------------------------------------------------------------
c
      type   = 'dpsi'
      bc_symm = 't odd'
c
#  if defined conjugate_gradient
#   ifdef timing
      call tic ('congrad', 'stream function (tropic)')
#   endif
      call congr (npt, type, bc_symm, ptd, ptd, ztd, res 
     &,           cf
     &,           mxscan, mscan, tolrsf
     &,           imask, iperm, jperm, iofs, nisle, nippts
     &,           converged, esterr)
#   ifdef timing
      call toc ('congrad', 'stream function (tropic)')
#   endif
#  endif
c
#  ifdef oldrelax
c
c     use sequential over-relaxation to solve the 5 pt laplacian
c     as in the codes of Cox (1984) and Semtner (1974).
c
#   ifdef timing
      call tic ('relax1', 'stream function (tropic)')
#   endif
      call relax1 (npt, type, bc_symm, ptd, ptd, ztd, res 
     &,            cf
     &,            sor, mxscan, mscan, tolrsf
     &,            imask, iperm, jperm, iofs, nisle, nippts
     &,            map
     &,            converged
     &,            esterr
     &            )
#   ifdef timing
      call toc ('relax1', 'stream function (tropic)')
#   endif
#  endif
c
#  ifdef hypergrid
c
c     use sequential over-relaxation to solve the 5 or 9 pt laplacian
c     along diagonals.
c
#   ifdef timing
      call tic ('hyper3', 'stream function (tropic)')
#   endif
      call hyper3 (npt, type, bc_symm, ptd, ptd, ztd, res
     &,            cf
     &,            sor, mxscan, mscan, tolrsf
     &,            imask, iperm, jperm, iofs, nisle, nippts
     &,            map
     &,            converged
     &,            esterr
     &            )
#   ifdef timing
      call toc ('hyper3', 'stream function (tropic)')
#   endif
#  endif
c
c     correct for drifting dpsi on land mass "imain"
c
      if (imain .gt. 0) then
        dpsi1 = ptd(iperm(iofs(imain)+1), jperm(iofs(imain)+1))
        call con_adjust (ptd, dpsi1, map)
      end if
c
c-----------------------------------------------------------------------
c     test accuracy of solving for change in stream function
c-----------------------------------------------------------------------
c
      if (.not.converged) then
        write (stdout,'(a,i5,3(a,1pe10.3))')
     &  ' WARNING: SOLVER DID NOT CONVERGE in ',mscan
     &, ' scans. max(psi)='
     &, absmax(psi(1,1,2)), ' max(dpsi)=',absmax(ptd)
     &, ' estimated max(err)=', esterr
      endif
c
c-----------------------------------------------------------------------
c     update the stream function based upon the solution
c-----------------------------------------------------------------------
c
      if (euler2) then
        do jrow=1,jmt
          do i=1,imt
            psi(i,jrow,1) = psi(i,jrow,2) + ptd(i,jrow)
          enddo
	enddo
      else
        do jrow=1,jmt
          do i=1,imt
            res(i,jrow)    = psi(i,jrow,2) + ptd(i,jrow)
            psi(i,jrow,2)  = psi(i,jrow,1)
            psi(i,jrow,1)  = res(i,jrow)
          enddo
	enddo
      endif
c
c-----------------------------------------------------------------------
c     save ptd to compute 1st guess for relaxation next timestep
c     (..note.. on 1st pass of euler backward timestep, bypass this
c            save, since it will be done on the 2nd pass)
c     (..note.. on a mixing timestep, alter ptd to be consistent with
c            normal, leap-frog stepping)
c-----------------------------------------------------------------------
c
      if (.not. euler1) then
c
        if (.not. leapfrog) then
          do jrow=1,jmt
            do i=1,imt
              ptd(i,jrow)=c2*ptd(i,jrow)
            enddo
	  enddo
        endif
c
        call oput (kflds, nwds, luptdb, ptd)
c
      endif
# endif
c
# ifdef timing
      call toc ('ocean', 'tropic')
# endif
c
      return
      end
#endif




      subroutine sfforc (zu, dxu, dyu, csu, forc)
c
c=======================================================================
c
c           S T R E A M   F U N C T I O N   F O R C I N G
c
c  
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
c
#include "size.h"
      dimension zu(imt,jmt,2), dxu(imt), dyu(jmt), csu(jmt) 
      dimension forc(imt,jmt)
      dimension ustuff(imt,jmt), vstuff(imt,jmt)
      dimension cddxu(0:1,0:1), cddyu(0:1,0:1)
      dimension cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)
      parameter (p5=0.5, c0=0.0)
c
#ifdef timing
      call tic ('sfforc', 'sfforc')
#endif
c
c-----------------------------------------------------------------------
c     initialize the forcing
c-----------------------------------------------------------------------
c
      do i=1,imt
        do jrow=1,jmt
          forc(i,jrow) = c0
        end do
      end do
c
c-----------------------------------------------------------------------
c     construct coefficients for partial differences. a partial
c     difference in "x" is defined as an "x" difference of a quantity 
c     which is averaged in "y". (and symmetrically for "y" differences).
c     Note that this is an x difference and NOT an x derivitive.
c     partial differences of quantities on the "t" grid are defined on
c     the "u" grid and visa versa.
c     therefore partial differences at: 
c     u/v points (i,j), involve nearby t/s points with subscripts: 
c        (i  ,j+1)    (i+1,j+1)
c        (i  ,j  )    (i+1,j  )
c     t/s points (i,j), involve nearby u/v points with subscripts:
c        (i-1,j  )    (i  ,j  )
c        (i-1,j-1)    (i  ,j-1)
c     thus if qu(i,j) is defined on u/v points, its partial
c     difference ddxqt = ddxt(qu) is defined on t/s points and has the
c     value 
c     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
c                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
c-----------------------------------------------------------------------
c
      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5
c
      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5
c
      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5
c
      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5
c
c-----------------------------------------------------------------------
c     multiply the u eqn by dx*cos, the v eqn by dy, then subtract their 
c     partial differences to eliminate the unknown surface pressure from 
c     the resulting equation
c-----------------------------------------------------------------------
c
      do i=1,imt-1
        do jrow=1,jmt-1
          ustuff(i,jrow) = zu(i,jrow,1)*dxu(i)*csu(jrow)
          vstuff(i,jrow) = zu(i,jrow,2)*dyu(jrow)
        end do
      end do
c
      do i1=-1,0
        do j1=-1,0
          do jrow=2,jmt-1
            do i=2,imt-1
              forc(i,jrow) = forc(i,jrow) 
     &                     - cddyt(i1,j1)*ustuff(i+i1,jrow+j1)
     &                     + cddxt(i1,j1)*vstuff(i+i1,jrow+j1)
            end do
          end do
        end do
      end do
c
#ifdef timing
      call toc ('sfforc', 'sfforc')
#endif
c
      return
      end




      subroutine sfc5pt (acor, f, c2dtsf, dxu, dyu, csu, hr, coef)
c
c=======================================================================
c
c     5  P T    C O E F F I C I E N T   I N I T I A L I A Z A T I O N
c
c     coefficient initialization for 5 point elliptic solvers
c
c     inputs:
c
c     acor   = implicit coriolis factor (0.0 => 1.0)
c     f      = 2*omega*sin(phi(j))
c     c2dtsf = twice the time step (seconds)
c     dxu    = width of "u" grid cell (cm)
c     dyu    = height of "u" grid cell (cm)
c     csu    = cosine of "u" grid cell
c     hr     = 1/depth at "u" cells (cm)
c
c     outputs:
c
c     coeff   = 3 x 3 array of coefficients at each (i,j) point
c
c     author: c.h. goldberg         e-mail => chg@gfdl.gov
c=======================================================================
c
#include "size.h"
      parameter (p5=0.5, c0=0.0)
      dimension f(jmt), csu(jmt), dxu(imt), dyu(jmt), hr(imt,jmt)
      dimension coef(imt,jmt,-1:1,-1:1)
      dimension ustuff(imt,jmt), vstuff(imt,jmt)
      dimension cddxu(0:1,0:1), cddyu(0:1,0:1)
      dimension cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)
c
#ifdef timing
      call tic ('sfc5pt', 'sfc5pt')
#endif
c
c-----------------------------------------------------------------------
c     initialize the coefficients
c-----------------------------------------------------------------------
c
      do jj=-1,1
        do ii=-1,1
          do j=1,jmt
            do i=1,imt
              coef(i,j,ii,jj) = c0
            end do
          end do
        end do
      end do
c
c-----------------------------------------------------------------------
c     construct coefficients for partial differences. a partial
c     difference in "x" is defined as an "x" difference of a quantity 
c     which is averaged in "y". (and symmetrically for "y" differences).
c     Note that this is an x difference and NOT an x derivitive.
c     partial differences of quantities on the "t" grid are defined on
c     the "u" grid and visa versa.
c     therefore partial differences at: 
c     u/v points (i,j), involve nearby t/s points with subscripts: 
c        (i  ,j+1)    (i+1,j+1)
c        (i  ,j  )    (i+1,j  )
c     t/s points (i,j), involve nearby u/v points with subscripts:
c        (i-1,j  )    (i  ,j  )
c        (i-1,j-1)    (i  ,j-1)
c     thus if qu(i,j) is defined on u/v points, its partial
c     difference ddxqt = ddxt(qu) is defined on t/s points and has the
c     value 
c     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
c                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
c-----------------------------------------------------------------------
c
      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5
c
      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5
c
      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5
c
      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5
c
c-----------------------------------------------------------------------
c     compute coefficients for all points
c-----------------------------------------------------------------------
c
      do i=1,imt-1
        do j=1,jmt-1
          ustuff(i,j) = dxu(i)*csu(j)*hr(i,j) / (c2dtsf*dyu(j))
          vstuff(i,j) = dyu(j)*hr(i,j) / (c2dtsf*dxu(i)*csu(j))
        end do
      end do
c
c-----------------------------------------------------------------------
c     calculate 5 point coefficients
c
c     note that ne [and nw] coefficient adds to n coefficient in 
c     ustuff term, but ne [and se] coefficient adds to e coefficient in
c     vstuff term for the 5 point operator.
c-----------------------------------------------------------------------
c
      do i1=0,1
        do j1=0,1
          do i2=-1,0
            do j2=-1,0
              do j=2,jmt-1
                do  i=2,imt-1
                  coef(i,j,0,j1+j2) = coef(i,j,0,j1+j2) +
     &               cddyu(i1,j1)*cddyt(i2,j2)*ustuff(i+i2,j+j2)
                  coef(i,j,i1+i2,0) = coef(i,j,i1+i2,0) +
     &               cddxu(i1,j1)*cddxt(i2,j2)*vstuff(i+i2,j+j2)
                end do
	      end do
            end do
          end do
        end do
      end do
c
c-----------------------------------------------------------------------
c     augment coefficients for implicit treatment of coriolis term
c     all coefficients are calculated, but corner ones are zero.
c-----------------------------------------------------------------------
c
      if (acor .ne. 0.0) then
        do i=1,imt-1
          do j=1,jmt-1
            ustuff(i,j) = acor*hr(i,j)*(-f(j))
            vstuff(i,j) = acor*hr(i,j)*( f(j))
          end do
        end do
        do i1=0,1
          do j1=0,1
            do i2=-1,0
              do j2=-1,0
                do j=2,jmt-1
                  do  i=2,imt-1
                    coef(i,j,i1+i2,j1+j2) = coef(i,j,i1+i2,j1+j2)
     &                 - cddxu(i1,j1)*cddyt(i2,j2)*ustuff(i+i2,j+j2)
                    coef(i,j,i1+i2,j1+j2) = coef(i,j,i1+i2,j1+j2)
     &                 - cddyu(i1,j1)*cddxt(i2,j2)*vstuff(i+i2,j+j2)
                  end do
                end do
              end do
            end do
          end do
        end do
      end if
c
#ifdef timing
      call toc ('sfc5pt', 'sfc5pt')
#endif
c
      return
      end



      subroutine sfc9pt  (acor, f, c2dtsf, dxu, dyu, csu, hr, coef)
c
c=======================================================================
c
c     9  P T    C O E F F I C I E N T   I N I T I A L I A Z A T I O N
c
c     coefficient initialization for 9 point elliptic solvers
c
c     inputs:
c
c     acor   = implicit coriolis factor (0.0 => 1.0)
c     f      = 2*omega*sin(phi(j))
c     c2dtsf = twice the time step (seconds)
c     dxu    = width of "u" grid cell (cm)
c     dyu    = height of "u" grid cell (cm)
c     csu    = cosine of "u" grid cell
c     hr     = 1/depth at "u" cells (cm)
c
c     outputs:
c
c     coeff   = 3 x 3 array of coefficients at each (i,j) point
c
c     author: c.h. goldberg         e-mail => chg@gfdl.gov
c=======================================================================
c
#include "size.h"
      parameter (c0=0.0, p5=0.5)
      dimension f(jmt), csu(jmt), dxu(imt), dyu(jmt), hr(imt,jmt)
      dimension coef(imt,jmt,-1:1,-1:1)
      dimension ustuff(imt,jmt), vstuff(imt,jmt)
      dimension cddxu(0:1,0:1), cddyu(0:1,0:1)
      dimension cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)
c
#ifdef timing
      call tic ('sfc9pt', 'sfc9pt')
#endif
c
c-----------------------------------------------------------------------
c     initialize the work area
c-----------------------------------------------------------------------
c
      do jj=-1,1
        do ii=-1,1
          do j=1,jmt
            do i=1,imt
              coef(i,j,ii,jj) = c0
            end do
          end do
        end do
      end do
c
c-----------------------------------------------------------------------
c     generate arrays of coefficients
c
c     construct coefficients for partial differences. a partial
c     difference in "x" is defined as an "x" difference of a quantity 
c     which is averaged in "y". (and symmetrically for "y" differences).
c     Note that this is an x difference and NOT an x derivitive.
c     partial differences of quantities on the "t" grid are defined on
c     the "u" grid and visa versa.
c     therefore partial differences at: 
c     u/v points (i,j), involve nearby t/s points with subscripts: 
c        (i  ,j+1)    (i+1,j+1)
c        (i  ,j  )    (i+1,j  )
c     t/s points (i,j), involve nearby u/v points with subscripts:
c        (i-1,j  )    (i  ,j  )
c        (i-1,j-1)    (i  ,j-1)
c     thus if qu(i,j) is defined on u/v points, its partial
c     difference ddxqt = ddxt(qu) is defined on t/s points and has the
c     value 
c     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
c                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
c-----------------------------------------------------------------------
c
      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5
c
      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5
c
      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5
c
      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5
c
c-----------------------------------------------------------------------
c     compute coefficients for all points
c-----------------------------------------------------------------------
c
      do i=1,imt-1
        do j=1,jmt-1
          ustuff(i,j) = dxu(i)*csu(j)*hr(i,j) / (c2dtsf*dyu(j))
          vstuff(i,j) = dyu(j)*hr(i,j) / (c2dtsf*dxu(i)*csu(j))
        end do
      end do
c
c---------------------------------------------------------------------
c     calculate 9 point coefficients
c---------------------------------------------------------------------
c
      do i1=0,1
        do j1=0,1
          do i2=-1,0
            do j2=-1,0
              do j=2,jmt-1
                do  i=2,imt-1
                  coef(i,j,i1+i2,j1+j2) = coef(i,j,i1+i2,j1+j2) +
     &               cddyu(i1,j1)*cddyt(i2,j2)*ustuff(i+i2,j+j2)  +
     &               cddxu(i1,j1)*cddxt(i2,j2)*vstuff(i+i2,j+j2)
                end do
              end do
            end do
          end do
        end do
      end do
c
c-----------------------------------------------------------------------
c     augment coefficients for implicit treatment of coriolis term
c     all coefficients are calculated, but corner ones are zero.
c-----------------------------------------------------------------------
c
      if (acor .ne. 0.0) then
        do i=1,imt-1
          do j=1,jmt-1
            ustuff(i,j) = acor*hr(i,j)*(-f(j))
            vstuff(i,j) = acor*hr(i,j)*( f(j))
          end do
        end do
        do i1=0,1
          do j1=0,1
            do i2=-1,0
              do j2=-1,0
                do j=2,jmt-1
                  do  i=2,imt-1
                    coef(i,j,i1+i2,j1+j2) = coef(i,j,i1+i2,j1+j2)
     &               - cddxu(i1,j1)*cddyt(i2,j2)*ustuff(i+i2,j+j2)
     &               - cddyu(i1,j1)*cddxt(i2,j2)*vstuff(i+i2,j+j2)
                  end do
                end do
              end do
            end do
          end do
        end do
      end if
c
#ifdef timing
      call toc ('sfc9pt', 'sfc9pt')
#endif
c
      return
      end


      subroutine spforc (zu, dxu, dyu, csu, h, forc)
c
c=======================================================================
c
c           S U R F A C E   P R E S S U R E   F O R C I N G
c
c  
c     author:      c. h. goldberg      e-mail=> chg@gfdl.gov
c=======================================================================
c
#include "size.h"
      dimension zu(imt,jmt,2)
      dimension dxu(imt), dyu(jmt), csu(jmt)
      dimension h(imt,jmt)
      dimension forc(imt,jmt)
      dimension ustuff(imt,jmt), vstuff(imt,jmt)
      dimension cddxu(0:1,0:1), cddyu(0:1,0:1)
      dimension cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)
      parameter (p5=0.5)
c
#ifdef timing
      call tic ('spforc', 'spforc')
#endif
c
c-----------------------------------------------------------------------
c     generate arrays of coefficients
c
c     construct coefficients for partial differences. a partial
c     difference in "x" is defined as an "x" difference of a quantity 
c     which is averaged in "y". (and symmetrically for "y" differences).
c     Note that this is an x difference and NOT an x derivitive.
c     partial differences of quantities on the "t" grid are defined on
c     the "u" grid and visa versa.
c     therefore partial differences at: 
c     u/v points (i,j), involve nearby t/s points with subscripts: 
c        (i  ,j+1)    (i+1,j+1)
c        (i  ,j  )    (i+1,j  )
c     t/s points (i,j), involve nearby u/v points with subscripts:
c        (i-1,j  )    (i  ,j  )
c        (i-1,j-1)    (i  ,j-1)
c     thus if qu(i,j) is defined on u/v points, its partial
c     difference ddxqt = ddxt(qu) is defined on t/s points and has the
c     value 
c     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
c                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
c-----------------------------------------------------------------------
c
      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5
c
      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5
c
      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5
c
      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5
c
c     weight "zu" and "zv" by the cell area and take the divergence
c
      do i=1,imt-1
	do j=1,jmt-1
	  ustuff(i,j) = h(i,j)*zu(i,j,1)*dyu(j)
	  vstuff(i,j) = h(i,j)*zu(i,j,2)*dxu(i)*csu(j)
	end do
      end do
c
      do i=1,imt
	do j=1,jmt
	  forc(i,j) = 0.0
	end do
      end do
c
      do i1=-1,0
	do j1=-1,0
	  do i=2,imt-1
	    do j=2,jmt-1
	      forc(i,j) = forc(i,j) + cddxt(i1,j1)*ustuff(i+i1,j+j1)
     &                              + cddyt(i1,j1)*vstuff(i+i1,j+j1)
	    end do
	  end do
	end do
      end do
c
#ifdef timing
      call toc ('spforc', 'spforc')
#endif
c
      return
      end



      subroutine spc9pt (dxu, dyu, csu, h, coef)
c
c=======================================================================
c
c     S U R F A C E   P R E S S U R E    C O E F F I C I E N T
c
c                I N I T I A L I A Z A T I O N
c
c     inputs:
c
c     dxu    = width of "u" grid cell (cm)
c     dyu    = height of "u" grid cell (cm)
c     csu    = cosine of "u" grid cell
c     h      = depth at "u,v" cells (cm)
c
c     outputs:
c
c
c     coeff   = 3 x 3 array of coefficients at each (i,j) point
c
c     author: c.h. goldberg         e-mail => chg@gfdl.gov
c=======================================================================
c
#include "size.h"
      parameter (c0=0.0, p5=0.5)
      dimension csu(jmt), dxu(imt), dyu(jmt), h(imt,jmt)
      dimension coef(imt,jmt,-1:1,-1:1)
      dimension ustuff(imt,jmt), vstuff(imt,jmt)
      dimension cddxu(0:1,0:1), cddyu(0:1,0:1)
      dimension cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)
c
#ifdef timing
      call tic ('spc9pt', 'spc9pt')
#endif
c
c-----------------------------------------------------------------------
c     generate arrays of coefficients
c
c     construct coefficients for partial differences. a partial
c     difference in "x" is defined as an "x" difference of a quantity 
c     which is averaged in "y". (and symmetrically for "y" differences).
c     Note that this is an x difference and NOT an x derivitive.
c     partial differences of quantities on the "t" grid are defined on
c     the "u" grid and visa versa.
c     therefore partial differences at: 
c     u/v points (i,j), involve nearby t/s points with subscripts: 
c        (i  ,j+1)    (i+1,j+1)
c        (i  ,j  )    (i+1,j  )
c     t/s points (i,j), involve nearby u/v points with subscripts:
c        (i-1,j  )    (i  ,j  )
c        (i-1,j-1)    (i  ,j-1)
c     thus if qu(i,j) is defined on u/v points, its partial
c     difference ddxqt = ddxt(qu) is defined on t/s points and has the
c     value 
c     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
c                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
c-----------------------------------------------------------------------
c
      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5
c
      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5
c
      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5
c
      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5
c
c-----------------------------------------------------------------------
c     compute coefficients for all points
c-----------------------------------------------------------------------
c
c     initialize all 9 coefficients to zero
c
      do i1=-1,1
        do j1=-1,1
          do i=1,imt
            do j=1,jmt
              coef(i,j,i1,j1) = c0
            end do
          end do
	end do
      end do
c
      do j=1,jmt
        do i=1,imt
	  ustuff(i,j) = 0.0
	  vstuff(i,j) = 0.0
	enddo
      enddo
      do i=1,imt-1
        do j=1,jmt-1
          ustuff(i,j) = h(i,j)*dyu(j)/(dxu(i)*csu(j))
          vstuff(i,j) = h(i,j)*dxu(i)*csu(j)/dyu(j)
        end do
      end do
c
c     calculate divergence = ddx (ddx (ustuff)) + ddy( ddy (vstuff))
c
      do i1=0,1
        do j1=0,1
          do i2=-1,0
            do j2=-1,0
              do i=2,imt-1
                do j=2,jmt-1
                  coef(i,j,i1+i2,j1+j2) = coef(i,j,i1+i2,j1+j2)
     &                 + cddxu(i1,j1) * cddxt(i2,j2) * ustuff(i+i2,j+j2)
     &                 + cddyu(i1,j1) * cddyt(i2,j2) * vstuff(i+i2,j+j2)
                end do
              end do
            end do
          end do
        end do
      end do
c
#ifdef timing
      call toc ('spc9pt', 'spc9pt')
#endif
      return
      end



      subroutine filz (fext, cf)
#if defined fourfil || defined firfil
c
c=======================================================================
c     subroutine filz sets up input needed for fourier filtering
c     (when the "fourfil" ifdef is defined) -or- symmetric finite
c     impulse response filtering (when the "firfil" ifdef is defined)
c     of "fext" at the specified high latitudes. "fext" is forcing for
c     the external mode.
c=======================================================================
c
# include "param.h"
# include "emode.h"
# include "grdvar.h"
# include "index.h"
# include "levind.h"
      dimension fext(imt,jmt)
      dimension cf(imt,jmt,3)
# ifdef fourfil
      dimension temp(imt)
# endif
# ifdef firfil
      dimension jtof(jmt)
# endif
c=======================================================================
c
# ifdef fourfil
      do jrow=jfrst,jmtm1
        if (jrow.le.jft1 .or. jrow.ge.jft2) then
c
c---------------------------------------------------------------------
c         fourier filter fext at high latitudes
c---------------------------------------------------------------------
c
            jj = jrow - jfrst + 1
            if (jrow .ge. jft2) jj = jj - jskpt + 1
            do l=1,lsegf
              is = iszf(jj,l)
              if (is .ne. 0) then
                ie = iezf(jj,l)
                do ii=is,ie
                  i = mod(ii-2,imtm2) + 2
                  temp(ii+1-is) = fext(i,jrow)
                enddo
                im = ie-is+1
#  ifdef cyclic
                if (im .ne. imtm2) then
                   m = 1
                   n = nint(im*cst(jrow)*cstr(jft0))
                else
                   m = 3
                   n = nint(im*cst(jrow)*cstr(jft0)*p5)
                endif
#  else
                m = 1
                n = nint(im*cst(jrow)*cstr(jft0))
#  endif
c
                call filtr (temp(1), im, m ,n, 0)
c
                do ii=is,ie
                  i = mod(ii-2,imtm2)+2
                  fext(i,jrow) = temp(ii+1-is)
                enddo
              endif
	    enddo
        endif
      enddo
# endif
# ifdef firfil
c
c-----------------------------------------------------------------------
c     filter external mode forcing at high latitudes with
c     symmetric finite impulse response filter
c-----------------------------------------------------------------------
c
c     construct mask and number of jrows to filter "jtof"
c
      do jrow=1,jmt
	jtof(jrow) = 0
        do i=1,imt
          cf(i,jrow,1) = c0
        enddo
      enddo
      do jrow=jfrst,jmt-1
        if (jrow.le.jft1 .or. jrow.ge.jft2) then
            jj = jrow - jfrst + 1
            if (jrow .ge. jft2) jj = jj - jskpt + 1
            jtof(jrow) = numflt(jj)
c
c           construct appropriate land/sea mask 
c           (1,0) for (filtered,non filtered) points
c
            do l=1,lseg
              isv = iszf(jj,l)
              if (isv .ne. 0) then
                iev = iezf(jj,l)
                do ii=isv,iev
                  i = mod(ii-2,imt-2) + 2
                  cf(i,jrow,1) = c1
                enddo
              endif
            enddo
c
	endif
      enddo
      call setbcx (cf(1,1,1), imt, jmt)
c
c     select points to filter (non filtered points = zero)
c
      do jrow=1,jmt
	do i=2,imtm1
	  cf(i,jrow,2) = fext(i,jrow)*cf(i,jrow,1)
	enddo
      enddo
      call setbcx (cf(1,1,2), imt, jmt)
c
# ifdef timing
      call tic ('filtering', 'filfir (finite impulse)')
# endif
c
c     each filtering consists of a double pass
c
      do j=2,jmtm1
        num = jtof(j)
	do n=1,num 
          do i=2,imtm1
            cf(i,j,3) = cf(i,j,1)*(p25*(cf(i-1,j,2) + cf(i+1,j,2)) +
     &               cf(i,j,2)*(c1 - p25*(cf(i-1,j,1) + cf(i+1,j,1))))
          enddo
# ifdef cyclic
          cf(1,j,3)   = cf(imtm1,j,3)
          cf(imt,j,3) = cf(2,j,3)
# else
          cf(1,j,3)   = c0
          cf(imt,j,3) = c0
# endif
          do i=2,imtm1
            cf(i,j,2) = cf(i,j,1)*(p25*(cf(i-1,j,3) + cf(i+1,j,3)) +
     &               cf(i,j,3)*(c1 - p25*(cf(i-1,j,1) + cf(i+1,j,1))))
          enddo
# ifdef cyclic
          cf(1,j,2)   = cf(imtm1,j,2)
          cf(imt,j,2) = cf(2,j,2)
# else
          cf(1,j,2)   = c0
          cf(imt,j,2) = c0
# endif
        enddo
      enddo
c
# ifdef timing
      call toc ('filtering', 'filfir (finite impulse)')
# endif
c
c     restore filtered "fext" on ocean points
c
      do jrow=2,jmtm1
	do i=1,imt
	  if (cf(i,jrow,1) .ne. c0) fext(i,jrow) = cf(i,jrow,2)
	enddo
      enddo
      call setbcx (fext, imt, jmt)      
# endif
c
#endif
      return
      end




