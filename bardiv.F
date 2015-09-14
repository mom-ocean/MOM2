      subroutine bardiv
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
c
c=======================================================================
c
c     compute uncorrected barotropic velocities and their divergence
c     for r.h.s. of surface pressure equation
c
c     Parallel Ocean General Circulation Modeling
c     R.D. Smith, J.K. Dukowicz and R.C. Malone
c     Physica D 60 (1992) 38-61
c     North-Holland
c     Elsevier Science Publishers B.V.
c
c     Implicit Free-Surface Method for the Bryan-Cox-Semtner Ocean
c     Model.
c     J.K. Dukowicz, R.D. Smith
c     Submitted to J. Geophysical Research   June 1993
c
c
c     author:       R. Smith             e-mail   rds@acl.lanl.gov
c                   A. Rosati            e-mail   ar@gfdl.gov
c=======================================================================
c
# include "param.h"
# include "coord.h"
# include "cprnts.h"
# include "diag.h"
# include "emode.h"
# include "grdvar.h"
# include "index.h"
# include "levind.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# include "tmngr.h"
      character *8 bc_symm
c
# ifdef timing
      call tic ('tropic', 'bardiv')
# endif
c
c-----------------------------------------------------------------------
c      construct uncorrected barotropic velocities
c      based on gradients of surface pressure
c-----------------------------------------------------------------------
c
      do i=1,imt
        ubarm1(i,1,1)   = c0
        ubarm1(i,jmt,1) = c0
        ubarm1(i,1,2)   = c0
        ubarm1(i,jmt,2) = c0
      enddo
c
      do jrow=2,jmtm1
	fx = acor*c2dtsf*cori(jrow,1)
	fy = c1/(c1 + fx**2)
        factu = p5*c2dtsf*csur(jrow)
        factv = p5*c2dtsf*dyur(jrow)
c
c       leapfrog
c
        if (leapfrog) then 
	  do i=2,imtm1
	    d1 = gam*(ps(i+1,jrow+1,1) - ps(i,jrow,1))
     &       + (c1-gam)*(ps(i+1,jrow+1,2) - ps(i,jrow,2))
	    d2 = gam*(ps(i+1,jrow,1) - ps(i,jrow+1,1))
     &       + (c1-gam)*(ps(i+1,jrow,2) - ps(i,jrow+1,2))
            utwid = zu(i,jrow,1)*c2dtsf - factu*(d1 + d2)*dxur(i)
            vtwid = zu(i,jrow,2)*c2dtsf - factv*(d1 - d2)
c
c           calculate uncorrected velocities at time tau+1
c
            uhat(i,jrow,1) = fy*(utwid + fx*vtwid) + ubarm1(i,jrow,1)
# ifdef implicit_free_surface
     &                    + ubar(i,jrow,1)
# endif
            uhat(i,jrow,2) = fy*(vtwid - fx*utwid) + ubarm1(i,jrow,2)  
# ifdef implicit_free_surface
     &                   + ubar(i,jrow,2)
# endif
          enddo
        elseif (euler1 .or. forward) then
c
c         forward or 1st pass of euler
c
	  do i=2,imtm1
	    d1 = ps(i+1,jrow+1,1) - ps(i,jrow,1)
	    d2 = ps(i+1,jrow,1) - ps(i,jrow+1,1)
            utwid = zu(i,jrow,1)*c2dtsf - factu*(d1 + d2)*dxur(i)
            vtwid = zu(i,jrow,2)*c2dtsf - factv*(d1 - d2)
c
c           calculate uncorrected velocities at time tau+1
c
            uhat(i,jrow,1) = fy*(utwid + fx*vtwid)
     &	     + ubarm1(i,jrow,1)
            uhat(i,jrow,2) = fy*(vtwid - fx*utwid)
     &       + ubarm1(i,jrow,2)
          enddo
        elseif (euler2) then
c
c         euler 2nd pass
c
	  do i=2,imtm1
	    d1 = theta*(pguess(i+1,jrow+1) - pguess(i,jrow))
     &         + (c1-theta)*(ps(i+1,jrow+1,1) - ps(i,jrow,1))
	    d2 = theta*(pguess(i+1,jrow) - pguess(i,jrow+1))
     &       + (c1-theta)*(ps(i+1,jrow,1) - ps(i,jrow+1,1))
            utwid = zu(i,jrow,1)*c2dtsf - factu*(d1 + d2)*dxur(i)
            vtwid = zu(i,jrow,2)*c2dtsf - factv*(d1 - d2)
c
c           calculate uncorrected velocities at time tau+1
c
            uhat(i,jrow,1) = fy*(utwid + fx*vtwid)
     &                    + ubarm1(i,jrow,1)
            uhat(i,jrow,2) = fy*(vtwid - fx*utwid)
     &                    + ubarm1(i,jrow,2)
          enddo
        else
	  write (stdout,*) '=>Error: leapfrog, euler1, forward, euler2='
     &,   leapfrog, euler1, forward, euler2
          stop '=>bardiv'
        endif
c
        do i=2,imtm1
	  if (kmu(i,jrow) .eq. 0 )then
	    uhat(i,jrow,1) = c0
	    uhat(i,jrow,2) = c0
          endif
        enddo
      enddo
c
      call border (uhat(1,1,1), 'u even')
      call border (uhat(1,1,2), 'u odd')
c
# if defined fourfil || defined firfil
c
c-----------------------------------------------------------------------
c       filter uhat and vhat at high latitudes
c-----------------------------------------------------------------------
c
        call filhat (uhat(1,1,1))
        call filhat (uhat(1,1,2))
# endif
c
c-----------------------------------------------------------------------
c     form divergence of uncorrected barotropic velocities for
c     r.h.s. of surface pressure eqn
c-----------------------------------------------------------------------
c
c     compute the divergence of the sea surface pressure gradients
c
      call spforc (uhat, dxu, dyu, csu, h, divf)
c
      fx = c1/(apgr*c2dtsf)
      do jrow=2,jmtm1
        do i=2,imtm1
	  divf(i,jrow) = fx*divf(i,jrow)
	enddo
      enddo
c
# ifdef implicit_free_surface
      if (euler2) then
c
c       correct r.h.s on 2nd pass euler backward
c
        fx = c1/(grav*dtsf*c2dtsf*apgr)
        do jrow=2,jmtm1
	  fy = dyt(jrow)*cst(jrow)*fx
	  do i=1,imtm1
	    divf(i,jrow) = divf(i,jrow) 
     &                   + fy*dxt(i)*(pguess(i,jrow) - ps(i,jrow,1))
          enddo
        enddo
      endif
# endif
c
        call border (divf, 't even')
c
c-----------------------------------------------------------------------
c     calculate guess for change in pressure from guess for pressure 
c-----------------------------------------------------------------------
c
      if (leapfrog) then
        do jrow=1,jmt
          do i=1,imt
            ptd(i,jrow) = pguess(i,jrow)-ps(i,jrow,2)
          enddo
        enddo
      elseif (euler1 .or. forward) then
        do jrow=1,jmt
          do i=1,imt
            ptd(i,jrow) = pguess(i,jrow)-ps(i,jrow,1)
          enddo
        enddo
      else
        do jrow=1,jmt
          do i=1,imt
            ptd(i,jrow) = c0 
          enddo
        enddo
      endif
c
c-----------------------------------------------------------------------
c     use conjugate gradient 9 point laplacian to solve directly for
c     the change in surface pressure.
c-----------------------------------------------------------------------
c
c     initialize coefficients for the conjugate gradient solver
c
      call spc9pt (dxu, dyu, csu, h, cf)
c
c-----------------------------------------------------------------------
c     add diag term to central coeff in the implicit free surface
c-----------------------------------------------------------------------
c
# ifdef implicit_free_surface
      do jrow=2,jmtm1
	fx = cst(jrow)*dyt(jrow)/(apgr*c2dtsf*dtsf*grav)
	do i=2,imtm1
          if (map(i,jrow) .le. 0) then
	    cf(i,jrow,0,0) = cf(i,jrow,0,0) - fx*dxt(i)
	  endif
        enddo
      enddo
# endif
      npt    = 9
      variable   = 'd(surf press)'
      nislsp = 0
      bc_symm = 't even'
# ifdef implicit_free_surface
      tolr   = tolrfs
# else
      tolr   = tolrsp
# endif
c
# ifdef timing
      call tic ('congrad', 'surface pressure (bardiv)')
# endif
      call congr (npt, variable, bc_symm, ptd, ptd, divf, res 
     &,           cf
     &,           mxscan, mscan, tolr
     &,           imask, iperm, jperm, iofs, nislsp, nippts
     &,           converged, esterr)
# ifdef timing
      call toc ('congrad', 'surface pressure (bardiv)')
# endif
c
# ifndef implicit_free_surface
c
c     remove null space from ptd (rigid lid only)
c
      call checkerboard (ptd, map)
      call border (ptd, bc_symm)
c
c     remove mean
c
      call zero_level (ptd, 'surf press', map, dxt, dyt, cst)
      call border (ptd, bc_symm)
# endif
c
c-----------------------------------------------------------------------
c     correct barotropic velocities with change in surface pressure
c     gradient. this completes the calculation of barotropic
c     velocity at tau+1
c     also update barotropic velocities for next time step
c-----------------------------------------------------------------------
c
      do jrow=2,jmtm1
        factu = p5*apgr*c2dtsf*csur(jrow)
        factv = p5*apgr*c2dtsf*dyur(jrow)
        do i=2,imtm1
          d1 = ptd(i+1,jrow+1) - ptd(i,jrow)
          d2 = ptd(i+1,jrow) - ptd(i,jrow+1) 
c
c         tau + 1
c
	  tempu = uhat(i,jrow,1) - factu*(d1 + d2)*dxur(i)
	  tempv = uhat(i,jrow,2) - factv*(d1 - d2)
	  if (leapfrog) then
# ifdef implicit_free_surface
	    tempu = tempu - ubar(i,jrow,1)
	    tempv = tempv - ubar(i,jrow,2)
# endif
c
c           tau - 1  <= tau
c
	    ubarm1(i,jrow,1) = ubar(i,jrow,1)
	    ubarm1(i,jrow,2) = ubar(i,jrow,2)
	  endif
c
c         tau  <= tau + 1
c
          ubar(i,jrow,1) = tempu
          ubar(i,jrow,2) = tempv
        enddo
        do i=2,imtm1
          if (kmu(i,jrow) .eq. 0 ) then
            ubar(i,jrow,1)   = c0
            ubar(i,jrow,2)   = c0
            ubarm1(i,jrow,1) = c0
            ubarm1(i,jrow,2) = c0
          endif
        enddo
      enddo
c
      call border (ubar(1,1,1), 'u even')
      call border (ubar(1,1,2), 'u odd')
      call border (ubarm1(1,1,1), 'u even')
      call border (ubarm1(1,1,2), 'u odd')
c
c-----------------------------------------------------------------------
c     update the surface pressure based upon the relaxation solution
c-----------------------------------------------------------------------
c
      if (leapfrog) then
c
c       leapfrog
c
	do jrow=1,jmt
	  do i=1,imt
            pnew           = ptd(i,jrow) + ps(i,jrow,2)
            pguess(i,jrow) = c3*(pnew - ps(i,jrow,1)) + ps(i,jrow,2)
            ps(i,jrow,2)   = ps(i,jrow,1)
            ps(i,jrow,1)   = pnew
          enddo
	enddo
      elseif (euler1) then
c
c       eb 1st pass
c
        do jrow=1,jmt
          do i=1,imt
# ifdef implicit_free_surface
	    pnew = ptd(i,jrow) + ps(i,jrow,1)
	    pguess(i,jrow) = pnew
            ps(i,jrow,2)   = ps(i,jrow,1)
            ps(i,jrow,1)   = pnew
# else
	    pnew = ptd(i,jrow) + ps(i,jrow,1)
	    pguess(i,jrow) = pnew
# endif
          enddo
	enddo
      elseif (forward) then
c
c       forward
c
        do jrow=1,jmt
          do i=1,imt
	    pnew           = ptd(i,jrow) + ps(i,jrow,1)
	    pguess(i,jrow) = c3*(pnew - ps(i,jrow,1)) + ps(i,jrow,2)
            ps(i,jrow,2)   = ps(i,jrow,1)
            ps(i,jrow,1)   = pnew
          enddo
	enddo
      elseif (euler2) then
c
c       eb 2nd pass
c
        do jrow=1,jmt
          do i=1,imt
	    pnew           = ptd(i,jrow) + pguess(i,jrow)
	    pguess(i,jrow) = c3*(pnew - ps(i,jrow,1)) + ps(i,jrow,2)
            ps(i,jrow,2)   = ps(i,jrow,1)
            ps(i,jrow,1)   = pnew
          enddo
	enddo
      endif
c
# ifdef remove_ps_checkerboard
#  if !defined implicit_free_surface
c
c     test accumulation of residual checkerboard 
c
      call checkerboard(ps(1,1,1), map)
      call border (ps(1,1,1), bc_symm)
      call zero_level (ps(1,1,1), 'surf press', map, dxt, dyt, cst)
      call border (ps(1,1,1), bc_symm)
#  endif
# endif
# ifdef timing
      call toc ('tropic', 'bardiv')
# endif
#endif
      return
      end

      subroutine filhat (uvhat)
#if defined fourfil || defined firfil
c
c=======================================================================
c     subroutine filhat sets up input needed for fourier filtering
c     (when the "fourfil" ifdef is defined) or anti-symmetric finite
c     impulse response filtering (when the "firfil" ifdef is defined)
c     of "uvhat" at the specified high latitudes. "uvhat" is forcing for
c     the external mode.
c
c     author: r.c.pacanowski     e-mail => rcp@gfdl.gov
c=======================================================================
c
# include "param.h"
# include "emode.h"
# include "grdvar.h"
# include "index.h"
# include "levind.h"
      dimension uvhat(imt,jmt)
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
        if (jrow.le.jfu1 .or. jrow.ge.jfu2) then
c
c---------------------------------------------------------------------
c         fourier filter uvhat at high latitudes
c---------------------------------------------------------------------
c
            jj = jrow - jfrst + 1
            if (jrow .ge. jfu2) jj = jj - jskpu + 1
            do l=1,lsegf
              is = iszf(jj,l)
              if (is .ne. 0) then
                ie = iezf(jj,l)
                do ii=is,ie
                  i = mod(ii-2,imtm2) + 2
                  temp(ii+1-is) = uvhat(i,jrow)
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
                  uvhat(i,jrow) = temp(ii+1-is)
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
#  ifdef timing
      call tic ('filtering', 'filfir (finite impulse)')
#  endif
      do jrow=jfrst,jmt-1
        if (jrow.le.jfu1 .or. jrow.ge.jfu2) then
          jj = jrow - jfrst + 1
          if (jrow .ge. jfu2) jj = jj - jskpu + 1
          jtof(jrow) = numflt(jj)
c
          do l=1,lseg
            isv = iszf(jj,l)
            iev = iezf(jj,l)
            if (isv .ne. 0) then
              call filter_it (uvhat(1,jrow), isv, iev, numflt(jj))
            endif
          enddo
	endif
      enddo
#  ifdef timing
      call toc ('filtering', 'filfir (finite impulse)')
#  endif
      call setbcx (uvhat, imt, jmt)      
# endif
c
      return
      end


      subroutine filter_it (a, iss, ies, num)
c
c=======================================================================
c     each filtering consists of a double pass using a anti-symmetric
c     finite impulse response filter.
c
c     author: r.c.pacanowski     e-mail => rcp@gfdl.gov
c=======================================================================
c
# include "size.h"
      dimension a(imt), b(0:imt), c(0:imt)
c
      imtm2 = imt-2
      len = ies-iss+1
c
c     extract the strip and apply anti-symmetry or cyclic condition
c
      do ii=iss,ies
        i = mod(ii-2,imtm2) + 2
        b(ii-iss+1) = a(i)
      enddo
c            
      b(0)     = 0.0
      b(len+1) = 0.0
      c(0)     = 0.0
      c(len+1) = 0.0
# ifdef cyclic
      if (len .eq. imtm2) then
	b(0)     = b(len)
	b(len+1) = b(1)
      endif
# endif
      do n=1,num 
        do i=1,len
          c(i) = 0.25*b(i-1) + 0.5*b(i) + 0.25*b(i+1)
        enddo
# ifdef cyclic
        if (len .eq. imtm2) then
	  c(0)     = c(len)
	  c(len+1) = c(1)
        endif
# endif
        do i=1,len
          b(i) = 0.25*c(i-1) + 0.5*c(i) + 0.25*c(i+1)
        enddo
# ifdef cyclic
        if (len .eq. imtm2) then
	  b(0)     = b(len)
	  b(len+1) = b(1)
        endif
# endif
      enddo      
c
c     restore the filtered strip
c
      do ii=iss,ies
        i = mod(ii-2,imtm2) + 2
        a(i) = b(ii-iss+1)
      enddo
#endif
      return
      end
