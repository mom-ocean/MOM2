      subroutine clinic (joff, js, je, is, ie)
c
c=======================================================================
c     compute internal mode velocity components for rows js through je
c     in the MW.
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
      parameter (istrt=2, iend=imt-1)
#include "coord.h"
#include "csbc.h"
#include "grdvar.h"
#include "hmixc.h"
#include "emode.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "vmixc.h"
      dimension tempik(imt,km,jsmw:jmw)
      dimension baru(imt,jsmw:jemw,2)
#include "fdifm.h"
# ifdef timing
      call tic ('ocean', 'clinic')
      call tic ('clinic', 'misc')
# endif
c
c-----------------------------------------------------------------------
c     limit the longitude indices based on those from the argument list
c     Note: these are currently bypassed. istrt and iend are set as
c           parameters to optimize performance
c-----------------------------------------------------------------------
c
c      istrt = max(2,is)
c      iend  = min(imt-1,ie)
c
c
c-----------------------------------------------------------------------
c     build coefficients to minimize advection and diffusion computation 
c-----------------------------------------------------------------------
c
#ifdef biharmonic
      do j=js,je+1
#else
      do j=js,je
#endif
        jrow = j + joff
        do i=istrt-1,iend
	  rho0csudxur(i,j) = rho0r*csur(jrow)*dxur(i)
	  csudxu2r(i,j)    = csur(jrow)*dxur(i)*p5
#ifdef consthmix
	  am_csudxtr(i,j)  = am*csur(jrow)*dxtr(i+1)
#endif
	enddo
      enddo
c
c-----------------------------------------------------------------------
c     construct the hydrostatic pressure gradients: 1 = dp/dx; 2 = dp/dy
c-----------------------------------------------------------------------
c
c     compute horizontal pressure gradient at the first level
c
      grav_rho0r = grav*rho0r
      do j=js,je
        jrow = j + joff
        fxa  = grav_rho0r*dzw(0)*csur(jrow)
        fxb  = grav_rho0r*dzw(0)*dyu2r(jrow)
        do i=istrt-1,iend
          t1              = rho(i+1,1,j+1) - rho(i  ,1,j)
          t2              = rho(i  ,1,j+1) - rho(i+1,1,j)
          grad_p(i,1,j,1) = (t1-t2)*fxa*dxu2r(i)
          grad_p(i,1,j,2) = (t1+t2)*fxb
        enddo
      enddo
c
c     compute the change in pressure gradient between levels
c
      do j=js,je+1
        do k=2,km
          do i=istrt-1,iend+1
	    tempik(i,k,j) = rho(i,k-1,j) + rho(i,k,j)
          enddo
        enddo
      enddo
c
      do j=js,je
        jrow = j + joff
        fxa = grav_rho0r*csur(jrow)*p5
        fxb = grav_rho0r*dyu4r(jrow)
        do k=2,km
          do i=istrt-1,iend
            t1              = tempik(i+1,k,j+1) - tempik(i  ,k,j)
            t2              = tempik(i  ,k,j+1) - tempik(i+1,k,j)
            grad_p(i,k,j,1) = fxa*(t1-t2)*dzw(k-1)*dxu2r(i)
            grad_p(i,k,j,2) = fxb*(t1+t2)*dzw(k-1)
          enddo
        enddo
      enddo
c
c     integrate downward from the first level
c
      do j=js,je
        do k=1,kmm1
          do i=istrt-1,iend
            grad_p(i,k+1,j,1) = grad_p(i,k,j,1) + grad_p(i,k+1,j,1)
            grad_p(i,k+1,j,2) = grad_p(i,k,j,2) + grad_p(i,k+1,j,2)
          enddo
        enddo
      enddo
c
      do j=js,je
        call setbcx (grad_p(1,1,j,1), imt, km)
        call setbcx (grad_p(1,1,j,2), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     solve for one component of velocity at a time
c     n = 1 => zonal component
c     n = 2 => meridional component
c-----------------------------------------------------------------------
c
      do n=1,2
c
#if !defined linearized_advection
c
c-----------------------------------------------------------------------
c       calculate 2*advective flux (for speed) across east face of
c       "u" cells.
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=1,km
            do i=istrt-1,iend
	      adv_fe(i,k,j) = adv_veu(i,k,j)*(u(i,  k,j,n,tau) + 
     &                                        u(i+1,k,j,n,tau))
            enddo
          enddo
        enddo
c
c-----------------------------------------------------------------------
c       2*advective flux across northern face of "u" cells is built
c       into ADV_Uy. (It's done this way for performance issues)
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c       diffusive flux across east face of "u" cell
c       diffusive flux across north face of "u" cell
c-----------------------------------------------------------------------
c
# if defined consthmix && !defined biharmonic
c
c       build diffusive flux on eastern face of "u" cells
c
        do j=js,je
          do k=1,km
            do i=istrt-1,iend              
              diff_fe(i,k,j) = am_csudxtr(i,j)*
     &                        (u(i+1,k,j,n,taum1) - u(i,k,j,n,taum1))
            enddo
          enddo
        enddo
c
c       diffusive flux on northern face of "u" cells is built
c       into DIFF_Uy
c
# endif
#endif
#if defined consthmix && defined biharmonic
c
c       calculate diffusive flux on eastern and northern faces of
c      "u" cells
c
        call delsqu (joff, js, je+1, istrt, iend, n)
#endif
#ifdef smagnlmix
c
c       calculate diffusive flux on eastern and northern faces of
c      "u" cells
c
        call smagnlm (joff, js, je, istrt, iend, n)
#endif
c
#ifdef free_slip
c
c       no viscous transfer of tangential stress to boundary
c
        do j=js,je
	  do k=1,km
	    do i=istrt-1,iend
	      diff_fe(i,k,j) = diff_fe(i,k,j)*
     &            (umask(i,k,j)*(c1-umask(i,k,j-1)) + umask(i,k,j-1))
	    enddo
	  enddo
        enddo
# if !defined consthmix || defined biharmonic
        do j=js-1,je
          do k=1,km
	    do i=istrt-1,iend
	      diff_fn(i,k,j) = diff_fn(i,k,j)*
     &            (umask(i,k,j)*(c1-umask(i-1,k,j)) + umask(i-1,k,j))
	    enddo
	  enddo
        enddo
# endif
#endif
c
#if !defined linearized_advection
c
c-----------------------------------------------------------------------
c       calculate 2*advective flux (for speed) on bottom face of
c       "u" cell. also diffusive flux on bottom face of "u" cell
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=1,kmm1
            do i=istrt,iend
	      adv_fb(i,k,j) = adv_vbu(i,k,j)*(u(i,k,  j,n,tau) +
     &                                        u(i,k+1,j,n,tau))
# if defined constvmix 
              diff_fb(i,k,j) = kappa_m_dzwr(k)*
     &                         (u(i,k,j,n,taum1) - u(i,k+1,j,n,taum1))
# else
              diff_fb(i,k,j) = visc_cbu(i,k,j)*dzwr(k)*
     &                         (u(i,k,j,n,taum1) - u(i,k+1,j,n,taum1))
# endif
            enddo
          enddo
        enddo
#endif
c
c-----------------------------------------------------------------------
c       set surface and bottom vert b.c. on "u" cells for mixing
c       set surface and bottom vert b.c. on "u" cells for advection
c-----------------------------------------------------------------------
c
        do j=js,je
          jrow = j + joff
          do i=istrt,iend
	    kb              = kmu(i,jrow)
            diff_fb(i,0,j)  = smf(i,j,n)
            diff_fb(i,kb,j) = bmf(i,j,n)
	    adv_fb(i,0,j)   = adv_vbu(i,0,j)*(u(i,1,j,n,tau) +
     &                                        u(i,1,j,n,tau)) 
	    adv_fb(i,km,j)  = c0 
          enddo
        enddo
c
#ifdef source_term
c
c-----------------------------------------------------------------------
c       set source term for "u" cell
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=1,km
            do i=istrt,iend
	      source(i,k,j) = c0
            enddo
          enddo
        enddo
#endif
c
c-----------------------------------------------------------------------
c       solve for the internal mode part of du/dt at center of 
c       "u" cells by neglecting the surface pressure gradients. use
c       statement functions to represent each component of the 
c       calculation. 
c-----------------------------------------------------------------------
c
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
	      u(i,k,j,n,taup1) = 
     &            (DIFF_Ux(i,k,j) + DIFF_Uy(i,k,j) + DIFF_Uz(i,k,j) 
     &             + DIFF_metric(i,k,j,n)
#if !defined linearized_advection
     &             - ADV_Ux(i,k,j) - ADV_Uy(i,k,j) - ADV_Uz(i,k,j)
     &             + ADV_metric(i,k,j,n) 
#endif
     &             - grad_p(i,k,j,n) + CORIOLIS(i,k,j,n)
#ifdef source_term
     &             + source(i,k,j)
#endif
     &              )*umask(i,k,j)
            enddo
          enddo
        enddo
#ifdef implicitvmix
c
c-----------------------------------------------------------------------
c       add in du/dt component due to implicit vertical diffusion
c-----------------------------------------------------------------------
c
        call ivdifu (joff, js, je, istrt, iend, n)
#endif
c
c-----------------------------------------------------------------------
c       construct diagnostics associated with velocity component "n" 
c-----------------------------------------------------------------------
c
        call diagc1 (joff, js, je, istrt, iend, n)
c
c-----------------------------------------------------------------------
c       construct the vertical average of du/dt for forcing
c       the barotropic equation
c-----------------------------------------------------------------------
c
        do j=js,je
          jrow = j + joff
          do i=istrt,iend
            zu(i,jrow,n) = c0
          enddo
        enddo
        do j=js,je
          jrow = j + joff
          do k=1,km
            fx = dzt(k)
            do i=istrt,iend
              zu(i,jrow,n) = zu(i,jrow,n) + u(i,k,j,n,taup1)*fx
            enddo
          enddo
        enddo
c
        do j=js,je
          jrow = j + joff
          do i=istrt,iend
            zu(i,jrow,n) = zu(i,jrow,n)*hr(i,jrow)
          enddo
        enddo
#ifdef symmetry
        do j=js,je
          jrow = j + joff
          if (jrow .eq. jmtm1 .and. n .eq. 2) then
            do i=istrt,iend
              zu(i,jrow,2)   =  c0
	      zu(i,jrow+1,2) = -zu(i,jrow-1,2)
	      zu(i,jrow+1,1) =  zu(i,jrow-1,1)
            enddo
          endif
        enddo
#endif
c
c-----------------------------------------------------------------------
c       end of velocity component "n" loop
c-----------------------------------------------------------------------
c
      enddo
c
c-----------------------------------------------------------------------
c     compute "tau+1" velocities accounting for implicit part of the
c     coriolis term if treated implicitly. velocities are in error by an
c     arbitrary constant related to neglecting the unknown surface
c     pressure gradients
c-----------------------------------------------------------------------
c
      if (acor .ne. c0) then
	do j=js,je
          jrow  = j + joff
          fx    = c2dtuv*acor*c2*omega*sine(jrow)
          detmr = c1/(c1 + fx*fx)
          do k=1,km
            do i=istrt,iend
              t1 = (u(i,k,j,1,taup1) + fx*u(i,k,j,2,taup1))*detmr
              t2 = (u(i,k,j,2,taup1) - fx*u(i,k,j,1,taup1))*detmr
              u(i,k,j,1,taup1) = u(i,k,j,1,taum1) + c2dtuv*t1
              u(i,k,j,2,taup1) = u(i,k,j,2,taum1) + c2dtuv*t2
            enddo
	  enddo
	enddo
      else
        do n=1,2
	  do j=js,je
            do k=1,km
              do i=istrt,iend
                u(i,k,j,n,taup1) = u(i,k,j,n,taum1)
     &                            + c2dtuv*u(i,k,j,n,taup1)
              enddo
	    enddo
	  enddo
	enddo
      endif
c
c-----------------------------------------------------------------------
c     subtract incorrect vertical means (related to ignoring horizontal
c     gradients of the surface pressure) to get pure internal modes.
c-----------------------------------------------------------------------
c
      do n=1,2
        do j=js,je
          do i=istrt,iend
            baru(i,j,n) = c0
          enddo
        enddo
        do j=js,je
          do k=1,km
            do i=istrt,iend
              baru(i,j,n) = baru(i,j,n) + u(i,k,j,n,taup1)*dzt(k)
            enddo
	  enddo
        enddo
        do j=js,je
          jrow  = j + joff
          do i=istrt,iend
            baru(i,j,n) = baru(i,j,n)*hr(i,jrow)
	  enddo
        enddo
        do j=js,je
          do k=1,km
            do i=istrt,iend
              u(i,k,j,n,taup1) = u(i,k,j,n,taup1)
     &                          - umask(i,k,j)*baru(i,j,n)
            enddo
	  enddo
          call setbcx (u(1,1,j,n,taup1), imt, km)
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     construct diagnostics involving internal mode velocity at "tau+1"
c-----------------------------------------------------------------------
c
      call diagc2 (joff, js, je, is, ie)
c
#if defined fourfil || defined firfil
c
c-----------------------------------------------------------------------
c     filter velocity components at high latitudes
c-----------------------------------------------------------------------
c
      if (istrt .eq. 2 .and. iend .eq. imt-1) then
        call filuv (joff, js, je)
      else
        write (stdout,'(a)') 
     &  'Error: filtering requires is=2 and ie=imt-1 in clinic'
        stop '=>clinic'
      endif
#endif
      do j=js,je
        call setbcx (u(1,1,j,1,taup1), imt, km)
        call setbcx (u(1,1,j,2,taup1), imt, km)
      enddo
c
#if !defined simple_sbc
c
c-----------------------------------------------------------------------
c     if needed, construct the Atmos S.B.C.(surface boundary conditions)
c     averaged over this segment      
c     reference: Pacanowski, R.C., Effect of Equatorial Currents
c                on Surface Stress (JPO, Vol 17, No. 6, June 1987)
c-----------------------------------------------------------------------
c
      iu = mapsbc(8)
      iv = mapsbc(9)
      if (iu .ne. 0 .and. iv .ne. 0) then
        call asbcu (joff, js, je, istrt, iend, iu, iv)
      endif
#endif
c
#ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     & "=> In clinic: js=",js," je=",je," joff=",joff
     &," jrows=",js+joff," to ",je+joff
#endif
#ifdef timing
      call toc ('clinic', 'misc')
      call toc ('ocean', 'clinic')
#endif
      return
      end




      subroutine diagc1 (joff, js, je, is, ie, n)
c
c-----------------------------------------------------------------------
c     construct diagnostics which don`t require internal mode velocity
c     at "tau+1" for each velocity component "n"
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       n    = (1,2) = (u,v) velocity component
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
#include "param.h"
#include "coord.h"
#include "diag.h"
#include "diaga.h"
#include "grdvar.h"
#include "hmixc.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
      dimension temp(imt,km)
c
#ifdef timing
        if (n .eq. 1) then
	  call tic ('clinic', 'diagnostics')
	else
	  call ticr ('clinic', 'diagnostics')
	endif
#endif
c
#if defined save_mixing_coeff
c
c-----------------------------------------------------------------------
c     diagnostic: estimate mixing coefficients on east, north, and
c                 bottom face of U cells from the flux
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (cmixts .and. n .eq. 1 .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'mixing coeff')
# endif
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=2,imt-1
              dudx = (u(i+1,k,j,1,taum1)-u(i,k,j,1,taum1))
     &                *csur(jrow)*dxtr(i+1) + epsln
              ce(i,k,j,1) = diff_fe(i,k,j)/dudx
              dudy = (u(i,k,j+1,1,taum1)-u(i,k,j,1,taum1))
     &                *dytr(jrow+1) + epsln
# if !defined consthmix || defined biharmonic || defined isopycmix
              cn(i,k,j,1) = diff_fn(i,k,j)/dudy
# else
              cn(i,k,j,1) = am
# endif
            enddo
          enddo
        enddo
        do j=js,je
          jrow = j + joff
          do k=1,km-1
            do i=2,imt-1
              dudz = (u(i,k,j,1,taum1)-u(i,k+1,j,1,taum1))
     &                *dzwr(k) + epsln
              cb(i,k,j,1) = diff_fb(i,k,j)/dudz
            enddo
          enddo
          do i=2,imt-1
            cb(i,km,j,1) = 0.0
          enddo
        enddo
c
        do j=js,je
          call setbcx (ce(1,1,j,1), imt, km)
          call setbcx (cn(1,1,j,1), imt, km)
          call setbcx (cb(1,1,j,1), imt, km)
        enddo
# ifdef timing
        call toc ('diagnostic', 'mixing coeff')
# endif
      endif
#endif
c
#ifdef time_step_monitor
c
c-----------------------------------------------------------------------
c     diagnostic: accumulate global kinetic energy on "tau" velocity
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c               (based on diagnostic by M. Cox)
c-----------------------------------------------------------------------
c
      if (tsits .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'time step monitor')
# endif
        do j=js,je
          jrow = j + joff
          fx = rho0*p5*csu(jrow)*dyu(jrow)
# ifdef symmetry
          if (jrow .eq. jmtm1) fx = fx*p5
# endif
          do k=1,km
            do i=is,ie
	      weight    = fx*dzt(k)*dxu(i)
              temp(i,k) = u(i,k,j,n,tau)**2*weight
            enddo
            do i=is,ie
              ektot(k,jrow) = ektot(k,jrow) + temp(i,k)
            enddo
          enddo
        enddo
# ifdef timing
        call toc ('diagnostic', 'time step monitor')
# endif
      endif
#endif
c
#ifdef energy_analysis
c
c-----------------------------------------------------------------------
c     diagnostic: integrate work done by the r.h.s. terms in the
c                  momentum equations.
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c               (this is not done the same way as in MOM 1)
c-----------------------------------------------------------------------
c
      if (glents .and. eots) call ge1 (joff, js, je, is, ie, n)
#endif
c
#ifdef term_balances
c
c-----------------------------------------------------------------------
c     diagnostic: integrate r.h.s. terms in the momentum equations
c                 over specified regional volumes
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (trmbts .and. eots) call utb1 (joff, js, je, is, ie, n)
#endif
#ifdef xbts
c
c-----------------------------------------------------------------------
c     diagnostic: accumulate r.h.s terms in the momentum equation
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (eots) call uxbt1 (joff, js, je, n)
#endif
#ifdef timing
        call toc ('clinic', 'diagnostics')
#endif
      return
      end






      subroutine diagc2 (joff, js, je, is, ie)
c
c-----------------------------------------------------------------------
c     construct diagnostics requiring internal mode velocity at "tau+1"
c     and those not dependent on velocity component fluxes.
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
#include "param.h"
#include "levind.h"
#include "scalar.h"
#include "switch.h"
c
#ifdef timing
      call ticr ('clinic', 'diagnostics')
#endif
#ifdef energy_analysis
c
c-----------------------------------------------------------------------
c     diagnostic: integrate work done by du/dt in the momentum equations
c                 the external mode part of "u" at "tau+1" will be 
c                 accounted for after the external mode is solved.
c                 also, integrate the work done by buoyancy.
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c               (this is not done the same way as in MOM 1)
c-----------------------------------------------------------------------
c
      if (glents .and. eots) then
        call ge2 (joff, js, je, is, ie, kmt, kmu, c2dtuv, grav, rho0r)
      endif
#endif
c
#ifdef term_balances
c
c-----------------------------------------------------------------------
c     diagnostic: add du/dt and implicit coriolis terms to the integrals
c                 over specified volumes. the external mode parts will
c                 be accounted for after the external mode is solved.                 
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (trmbts .and. eots) then
        call utb2 (joff, js, je, is, ie, c2dtuv, acor)
      endif
#endif
c
#ifdef xbts
c
c-----------------------------------------------------------------------
c     diagnostic: accumulate du/dt and implicit coriolis terms from the 
c                 momentum equations
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (eots) call uxbt2 (joff, js, je, c2dtuv, acor)
#endif
#ifdef timing
      call toc ('clinic', 'diagnostics')
#endif
      return
      end



#if !defined simple_sbc
      subroutine asbcu (joff, js, je, is, ie, iu, iv)
c
c-----------------------------------------------------------------------
c     construct the Atmos S.B.C.(surface boundary conditions)
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       iu   = index for u component 
c       iv   = index for v component 
c
c     reference: Pacanowski, R.C., Effect of Equatorial Currents
c                on Surface Stress (JPO, Vol 17, No. 6, June 1987)
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
# include "param.h"
# include "csbc.h"
# include "levind.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# ifdef timing
      call tic ('clinic', 'constructing SBC')
# endif
c
c     initialize S.B.C. at the beginning of each ocean segment
c     (do not alter values in land)
c
      if (eots .and. osegs .and. iu .ne. 0 .and. iv .ne. 0) then
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) then
              sbcocn(i,jrow,iu) = c0
              sbcocn(i,jrow,iv) = c0
            endif
          enddo
        enddo
      endif
c
c     accumulate surface currents for the Atmos S.B.C. every time step
c
      if (eots .and. iu .ne. 0 .and. iv .ne. 0) then
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            sbcocn(i,jrow,iu) = sbcocn(i,jrow,iu) + p25*(
     &                          u(i,1,j,1,tau) + u(i-1,1,j,1,tau) 
     &                        + u(i,1,j-1,1,tau) + u(i-1,1,j-1,1,tau))
            sbcocn(i,jrow,iv) = sbcocn(i,jrow,iv) + p25*(
     &                          u(i,1,j,2,tau) + u(i-1,1,j,2,tau) 
     &                        + u(i,1,j-1,2,tau) + u(i-1,1,j-1,2,tau))
          enddo
        enddo
      endif
c
c     average the surface currents for the Atmos S.B.C. at the end
c     of each ocean segment. (do not alter values in land)
c
      if (eots .and. osege .and. iu .ne. 0 .and. iv .ne. 0) then
        rts = c1/ntspos
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) then
              sbcocn(i,jrow,iu) = rts*sbcocn(i,jrow,iu)
	      sbcocn(i,jrow,iv) = rts*sbcocn(i,jrow,iv)
            endif
          enddo
        enddo
      endif
# ifdef timing
      call toc ('clinic', 'constructing SBC')
# endif
      return
      end
#endif





#ifdef implicitvmix
      subroutine ivdifu (joff, js, je, is, ie, n)
c
c-----------------------------------------------------------------------
c     solve vertical diffusion of velocity implicitly
c
c     input:
c       joff  = offset relating "j" in the MW to latitude "jrow"
c       js    = starting row in the MW
c       je    = ending row in the MW
c       is    = starting longitude index in the MW
c       ie    = ending longitude index in the MW
c       n     = velocity component
c       twodt = (2*dtuv, dtuv) on (leapfrog, mixing) time steps
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
# include "size.h"
# include "levind.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
# include "vmixc.h"
      dimension twodt_rho0r(km)
c
# ifdef timing
        call tic ('clinic', 'implicit vert diff')
# endif
c
c     set some constants
c
      c1 = 1.0      
c
c     construct the "tau+1" velocity without implicit vertical diffusion
c
      do j=js,je
        do k=1,km
          do i=is,ie
            u(i,k,j,n,taup1) = u(i,k,j,n,taum1)+c2dtuv*u(i,k,j,n,taup1)
          enddo
        enddo
      enddo
c
c     store terms to compute implicit vertical diffusion on
c     diagnostic time steps
c
# ifdef xbts
      if (eots) then
# else
      if ((glents .or. trmbts) .and. eots) then
# endif
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = u(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif
c
c     add in the implicit vertical diffusion
c
      do k=1,km
        twodt_rho0r(k) = c2dtuv*rho0r
      enddo
# if !defined constvmix
#  ifdef tcvmix
      call invtri (up1(1,1,1,n), smf(1,jsmw,n), bmf(1,jsmw,n), vvca(???)
     &,            twodt_rho0r, kmu, umask(1,1,1), is, ie, joff, js, je)
      if (imt-1 .gt. 1) then
        print *,'Error: ivdif.F is not converted for ifvef tvcmi x'
	stop '=>ivdif'
      endif
#  else
      call invtri (u(1,1,1,n,taup1), smf(1,jsmw,n), bmf(1,jsmw,n)
     &, visc_cbu(1,1,jsmw), twodt_rho0r, kmu, umask(1,1,1), is, ie
     &, joff, js, je)
#  endif
# endif
c
c     compute residual implicit vertical diffusion
c
      r2dtuv = c1/c2dtuv
# ifdef xbts
      if (eots) then
# else
      if ((glents .or. trmbts) .and. eots) then
# endif
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = r2dtuv*(u(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
          enddo
        enddo
      endif
c
c     convert back to time change of velocity
c
      do j=js,je
        do k=1,km
          do i=is,ie
            u(i,k,j,n,taup1) =r2dtuv*(u(i,k,j,n,taup1)-u(i,k,j,n,taum1))
          enddo
        enddo
      enddo
# ifdef timing
        call toc ('clinic', 'implicit vert diff')
# endif
      return
      end
#endif




#ifdef biharmonic
      subroutine delsqu (joff, js, je, is, ie, n)
c
c=======================================================================
c     compute del**2 of prognostic variables on MW rows "js" ... "je"
c     and gradients of del**2 quantities for variable "n" in the
c     momentum equation
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       n    = (1,2) = (u,v) velocity component
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#include "grdvar.h"
#include "hmixc.h"
#include "mw.h"
#include "fdifm.h"
c
# ifdef timing
        call tic ('clinic', 'biharmonic mixing')
# endif
c
c       set -del**2 = 0 on southern wall
c
        j    = js
        jrow = j + joff
        if (jrow .eq. 2) then
          do k=1,km
            do i=is-1,ie+1
	      del2(i,k,j-1,n) = c0
            enddo
          enddo
	endif
c
c       set -del**2 = 0 on northern wall
c
        j    = je
        jrow = j + joff
        if (jrow .eq. jmt-1) then
          do k=1,km
            do i=is-1,ie+1
	      del2(i,k,j+1,n) = c0
            enddo
          enddo
	endif
c
c-----------------------------------------------------------------------
c       diffusive flux across east face of "u" cell
c       diffusive flux across north face of "u" cell
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=1,km
          jrow = j + joff
	  ambi_csur = visc_ceu*csur(jrow)
            do i=is-1,ie              
              diff_fe(i,k,j) = ambi_csur*dxtr(i+1)*
     &                         (u(i+1,k,j,n,taum1) - u(i,k,j,n,taum1))
            enddo
          enddo
        enddo
        do j=js-1,je
          jrow = j + joff
	  ambi_cst_dytr = visc_cnu*cst(jrow+1)*dytr(jrow+1)
          do k=1,km
            do i=is,ie
              diff_fn(i,k,j) = ambi_cst_dytr* 
     &                        (u(i,k,j+1,n,taum1) - u(i,k,j,n,taum1))
            enddo
          enddo
        enddo
c
c       compute -am*del**2
c
        do j=js,je
          jrow   = j + joff
          do k=1,km
            do i=is,ie
	      del2(i,k,j,n) = (-DIFF_Ux(i,k,j) - DIFF_Uy(i,k,j) 
     &                      - am3(jrow)*u(i,k,j,n,taum1) - am4(jrow,n)*
     &         (u(i+1,k,j,3-n,taum1) - u(i-1,k,j,3-n,taum1))*dxmetr(i)
     &                        )*umask(i,k,j)

            enddo
          enddo
	  call setbcx (del2(1,1,j,n), imt, km)
        enddo
c
c-----------------------------------------------------------------------
c       diffusive flux across eastern face of "u" cell
c       diffusive flux across northern face of "u" cell
c-----------------------------------------------------------------------
c
        do j=js,je-1
          jrow = j + joff
	  ambi_csur = visc_ceu*csur(jrow)
          do k=1,km
            do i=is-1,ie
              diff_fe(i,k,j) = ambi_csur*dxtr(i+1)*
     &                        (del2(i+1,k,j,n)-del2(i,k,j,n))
            enddo
          enddo
        enddo
        do j=js-1,je-1
          jrow = j + joff
	  ambi_cst_dytr = visc_cnu*cst(jrow+1)*dytr(jrow+1)
          do k=1,km
            do i=is,ie
              diff_fn(i,k,j) = ambi_cst_dytr*
     &                        (del2(i,k,j+1,n) - del2(i,k,j,n))
            enddo
          enddo
        enddo
# ifdef timing
        call toc ('clinic', 'biharmonic mixing')
# endif
        return
	end
#endif

