      subroutine tracer (joff, js, je, is, ie)
c
c=======================================================================
c     compute tracers at "tau+1" for rows js through je in the MW.
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
#include "accel.h"
#include "coord.h"
#include "cregin.h"
#include "csbc.h"
#include "emode.h"
#include "grdvar.h"
#include "hmixc.h"
#ifdef isopycmix
# include "isopyc.h"
#endif
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "vmixc.h"
      dimension twodt(km)
#include "fdift.h"
# ifdef timing
      call tic ('ocean', 'tracer')
      call tic ('tracer', 'misc')
# endif
c
c-----------------------------------------------------------------------
c     limit the longitude indices based on those from the argument list
c     Note: this is currently bypassed. istrt and iend are set as
c           parameters to optimize performance
c-----------------------------------------------------------------------
c
c      istrt = max(2,is)
c      iend  = min(imt-1,ie)
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
	  cstdxtr(i,j)    = cstr(jrow)*dxtr(i)
	  cstdxt2r(i,j)   = cstr(jrow)*dxtr(i)*p5
#if defined smagnlmix || defined bryan_lewis_horizontal
	  cstdxur(i,j) = cstr(jrow)*dxur(i)
#endif
#ifdef consthmix
	  ah_cstdxur(i,j) = ah*cstr(jrow)*dxur(i)
#endif
	enddo
      enddo
c
c-----------------------------------------------------------------------
c     solve for one tracer at a time
c     n = 1 => temperature
c     n = 2 => salinity
c     n > 2 => other tracers (if applicable)
c-----------------------------------------------------------------------
c
      do n=1,nt
c
#if !defined linearized_advection
c
c-----------------------------------------------------------------------
c       calculate 2*advective flux across eastern face of "T" cells.
c       (It's done this way for performance issues)
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=1,km
            do i=istrt-1,iend
	      adv_fe(i,k,j) = adv_vet(i,k,j)*(t(i,  k,j,n,tau) + 
     &                                        t(i+1,k,j,n,tau))
            enddo
          enddo
        enddo
c
c-----------------------------------------------------------------------
c       2*advective flux across northern face of "T" cells is built
c       into ADV_Ty. (It's done this way for performance issues)
c-----------------------------------------------------------------------
c
#endif
c
c-----------------------------------------------------------------------
c       calculate 2*advective flux across bottom face of "T" cells.
c       (It's done this way for performance issues)
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=1,km-1
            do i=istrt,iend
#ifdef linearized_advection
	      adv_fb(i,k,j)  = adv_vbt(i,k,j)*(tbarz(k) + tbarz(k+1))
#else
	      adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,  j,n,tau) +
     &                                         t(i,k+1,j,n,tau))
#endif
            enddo
          enddo
        enddo
c
c-----------------------------------------------------------------------
c       calculate diffusive flux across eastern and northern faces 
c       of "T" cells due to various parameterizations for diffusion.
c-----------------------------------------------------------------------
c
#if defined consthmix && !defined biharmonic
c
c       diffusive flux on eastern face of "T" cells
c
        do j=js,je
          do k=1,km
            do i=istrt-1,iend              
              diff_fe(i,k,j) = 
# ifdef bryan_lewis_horizontal
     &                         diff_cet(k)*cstdxur(i,j)*
# else
     &                         ah_cstdxur(i,j)*
# endif
     &             (t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
# if defined isopycmix
c
c       diffusive flux on northern face of "T" cells is calculated for
c       use as a background for isopycnal mixing
c
        do j=js-1,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) =
#  ifdef bryan_lewis_horizontal
     &                         diff_cnt(k)*
#  else
     &                         diff_cnt*
#  endif
     &           csu_dyur(jrow)*(t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
# else
c
c       diffusive flux on northern face of "T" cells is not calculated.
c       instead, it is incorporated into DIFF_Ty for performance issues.
c
# endif
#endif
#if defined consthmix && defined biharmonic
c
c       diffusive flux on eastern and northern faces of "T" cells
c      
        call delsqt (joff, js, je+1, istrt, iend, n)
#endif
#ifdef smagnlmix
c
c       diffusive flux on eastern and northern faces of "T" cells
c      
        do j=js,je
          do k=1,km
            do i=istrt-1,iend              
              diff_fe(i,k,j) = diff_cet(i,k,j)*cstdxur(i,j)*
     &                         (t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
        do j=js-1,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) = diff_cnt(i,k,j)*csu_dyur(jrow)*
     &                     (t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
#endif
c
c-----------------------------------------------------------------------
c       calculate diffusive flux across bottom face of "T" cells
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=1,km-1
            do i=istrt,iend
              diff_fb(i,k,j) = 
# if defined constvmix && !defined isopycmix
     &                         kappa_h_dzwr(k)*
# else
     &                         diff_cbt(i,k,j)*dzwr(k)*
# endif
     &                         (t(i,k,j,n,taum1) - t(i,k+1,j,n,taum1))
            enddo
          enddo
        enddo
c
c-----------------------------------------------------------------------
c       set boundary conditions on "t" cells for vertical diffusion
c       and vertical advection
c-----------------------------------------------------------------------
c
        do j=js,je
          jrow   = j + joff
          do i=istrt,iend
	    kb              = kmt(i,jrow)
            diff_fb(i,0,j)  = stf(i,j,n)
            diff_fb(i,kb,j) = btf(i,j,n)
	    adv_fb(i,0,j)   = adv_vbt(i,0,j)*(t(i,1,j,n,tau) +
     &                                       t(i,1,j,n,tau)) 
	    adv_fb(i,km,j)  = c0 
          enddo
        enddo
c
#ifdef isopycmix
c
c-----------------------------------------------------------------------
c       add isopycnal mixing components to previously computed
c       diffusive flux on eastern and northern faces of "T" cells.
c       xz and yz isopycnal diffusive flux at bottom of "t" cells is
c       in "diff_fbiso". All of these components are solved explicitly.
c       The zz component will be solved implicitly.
c-----------------------------------------------------------------------
c
        call isoflux (joff, js, je, is, ie, n)
#endif
c
# ifdef source_term
c
c-----------------------------------------------------------------------
c       set source term for "T" cells
c-----------------------------------------------------------------------
c
        do j=js,je
          do k=1,km
            do i=istrt,iend
	      source(i,k,j) = c0
            enddo
          enddo
        enddo
c
# ifdef sponges
c
c-----------------------------------------------------------------------
c       add newtonian damping term for sponge layers to source
c-----------------------------------------------------------------------
c
        call sponge1 (joff, js, je, istrt, iend, n, t(1,1,1,1,taum1)
     &,               source)
# endif
c
# if defined shortwave && !defined simple_sbc
c
c-----------------------------------------------------------------------
c       incorporate short wave penetration into source
c-----------------------------------------------------------------------
c
        if (n .eq. 1) then
	  call swflux0 (joff, js, je, istrt, iend, source)
	endif
# endif
#endif
c
c-----------------------------------------------------------------------
c       solve for "tau+1" tracer at center of "T" cells using statement
c       functions to represent each component of the calculation 
c-----------------------------------------------------------------------
c
        do j=js,je
          jrow   = j + joff
          do k=1,km
	    twodt(k) = c2dtts*dtxcel(k)
            do i=istrt,iend
	      t(i,k,j,n,taup1) = t(i,k,j,n,taum1) + twodt(k)*(
     &                  DIFF_Tx(i,k,j) + DIFF_Ty(i,k,j) + DIFF_Tz(i,k,j)
     &                 - ADV_Tx(i,k,j) - ADV_Ty(i,k,j) - ADV_Tz(i,k,j)
#if defined isopycmix && gent_mcwilliams
     &                 - ADV_Txiso(i,k,j) - ADV_Tyiso(i,k,j)
     &                 - ADV_Tziso(i,k,j)
#endif
#ifdef source_term
     &                 + source(i,k,j)
#endif
     &                           )*tmask(i,k,j)
            enddo
          enddo
        enddo
c
#if defined implicitvmix || defined isopycmix
c
c-----------------------------------------------------------------------
c       add dT/dt component due to implicit vertical diffusion
c-----------------------------------------------------------------------
c
        call ivdift (joff, js, je, istrt, iend, n, twodt)
#endif
c
        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
	enddo
c
c-----------------------------------------------------------------------
c       construct diagnostics associated with tracer "n" 
c-----------------------------------------------------------------------
c
        call diagt1 (joff, js, je, istrt, iend, n, twodt)
c
c-----------------------------------------------------------------------
c       end of tracer component "n" loop
c-----------------------------------------------------------------------
c
      enddo
c
#if !defined implicitvmix
c
c-----------------------------------------------------------------------
c     explicit convection: adjust column if gravitationally unstable
c-----------------------------------------------------------------------
c
# ifdef fullconvect
      call convct2 (t(1,1,1,1,taup1), joff, js, je, istrt, iend, kmt)
# else
      call convct (t(1,1,1,1,taup1), ncon, joff, js, je, istrt, iend
     &,            kmt)
# endif
#endif
c
c-----------------------------------------------------------------------
c     construct diagnostics after convection 
c-----------------------------------------------------------------------
c
      idiag = 10
      call diagt2 (joff, js, je, istrt, iend, idiag)
c
#if defined fourfil || defined firfil
c
c-----------------------------------------------------------------------
c     filter tracers at high latitudes
c-----------------------------------------------------------------------
c
      if (istrt .eq. 2 .and. iend .eq. imt-1) then
        call filt (joff, js, je)
      else
        write (stdout,'(a)') 
     &  'Error: filtering requires is=2 and ie=imt-1 in tracer'
        stop '=>tracer'
      endif
#endif
      do n=1,nt
        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     construct diagnostics after filtering (for total dT/dt) 
c-----------------------------------------------------------------------
c
      idiag = 1
      call diagt2 (joff, js, je, istrt, iend, idiag)
c
#if !defined simple_sbc
c
c-----------------------------------------------------------------------
c     if needed, construct the Atmos S.B.C.(surface boundary conditions)
c     averaged over this segment
c     eg: SST and possibly SSS
c-----------------------------------------------------------------------
c
      isst = mapsbc(6)
      isss = mapsbc(7)
      if (isst .ne. 0 .or. isss .ne. 0) then
        call asbct (joff, js, je, istrt, iend, isst, isss)
      endif
#endif
c
#ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     & "=> In tracer: js=",js," je=",je," joff=",joff
     &," jrows=",js+joff," to ",je+joff
#endif
#ifdef timing
      call toc ('tracer', 'misc')
      call toc ('ocean', 'tracer')
#endif
      return
      end





      subroutine diagt1 (joff, js, je, is, ie, n, twodt)
c
c-----------------------------------------------------------------------
c     construct diagnostics associated with tracer component "n" 
c
c     input:
c       joff  = offset relating "j" in the MW to latitude "jrow"
c       js    = starting row in the MW
c       je    = ending row in the MW
c       is    = starting longitude index in the MW
c       ie    = ending longitude index in the MW
c       n     = (1,2) = (u,v) velocity component
c       twodt = (2*dtts,dtts) on (leapfrog,mixing) time steps
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
#include "param.h"
#include "accel.h"
#include "coord.h"
#include "cregin.h"
#include "csbc.h"
#ifdef tracer_averages
# include "ctavg.h"
#endif
#include "diag.h"
#include "diaga.h"
#include "emode.h"
#include "grdvar.h"
#include "hmixc.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "vmixc.h"
#ifdef meridional_tracer_budget
# include "ctmb.h"
      dimension stor(imt,km), div(imt,km), sorc(imt,km), dif(imt,km)
      dimension t1(imt), t2(imt), t3(imt), t4(imt)
#endif
#ifdef time_step_monitor
      dimension temp1(imt,km), temp2(imt,km), temp3(imt,km)
#endif
      dimension twodt(km)
#include "fdift.h"
c
#ifdef timing
      if (n .eq. 1) then
	call tic ('tracer', 'diagnostics')
      else
	call ticr ('tracer', 'diagnostics')
      endif
#endif
c
#if defined save_mixing_coeff
c
c-----------------------------------------------------------------------
c     diagnostic: estimate mixing coefficients on east, north, and
c                 bottom face of T cells from the flux
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
              dtdx = (t(i+1,k,j,1,taum1)-t(i,k,j,1,taum1))
     &                *tmask(i+1,k,j)*tmask(i,k,j)
     &                *cstr(jrow)*dxur(i) + epsln
              ce(i,k,j,2) = diff_fe(i,k,j)/dtdx
              dtdy = (t(i,k,j+1,1,taum1)-t(i,k,j,1,taum1))
     &                *tmask(i,k,j+1)*tmask(i,k,j)
     &                *dyur(jrow) + epsln
# if !defined consthmix || defined biharmonic || defined isopycmix
              cn(i,k,j,2) = diff_fn(i,k,j)/dtdy
# else
              cn(i,k,j,2) = ah
# endif
            enddo
          enddo
        enddo
        do j=js,je
          jrow = j + joff
          do k=1,km-1
            do i=2,imt-1
              dtdz = (t(i,k,j,1,taum1)-t(i,k+1,j,1,taum1))
     &                *tmask(i,k,j)*tmask(i,k+1,j)
     &                *dzwr(k) + epsln
              cb(i,k,j,2) = diff_fb(i,k,j)/dtdz
            enddo
          enddo
          do i=2,imt-1
            cb(i,km,j,2) = 0.0
          enddo
        enddo
c
        do j=js,je
          call setbcx (ce(1,1,j,2), imt, km)
          call setbcx (cn(1,1,j,2), imt, km)
          call setbcx (cb(1,1,j,2), imt, km)
        enddo
# ifdef timing
        call toc ('diagnostic', 'mixing coeff')
# endif
      endif
#endif
c
# if defined save_convection
c
c-----------------------------------------------------------------------
c     diagnostic: initialize temerature before convection
c
c     author:   S. Griffies       e-mail smg@gfdl.gov
c-----------------------------------------------------------------------
c
      if (exconvts .and. n .eq. 1 .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'explicit convection')
# endif
        do j=js,je
          do k=1,km
            do i=1,imt
              excnv0(i,k,j) = t(i,k,j,1,taup1)
            enddo
          enddo
        enddo
# ifdef timing
        call toc ('diagnostic', 'explicit convection')
# endif
      endif
c
#endif
c
#ifdef time_step_monitor
c
c-----------------------------------------------------------------------
c     diagnostic: integrate |d(tracer)/dt|  and tracer variance on "tau"
c                 globally
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
          r2dt    = c1/c2dtts
          cosdyt  = cst(jrow)*dyt(jrow)
          do k=1,km
            fx = r2dt/dtxcel(k)
            do i=is,ie
              darea      = dzt(k)*dxt(i)*cosdyt*tmask(i,k,j)
              temp3(i,k) = t(i,k,j,n,tau)*darea
              temp1(i,k) = t(i,k,j,n,tau)**2*darea
              temp2(i,k) = abs(t(i,k,j,n,taup1)-t(i,k,j,n,taum1))*
     &                     darea*fx
            enddo
            do i=is,ie
              tbar(k,n,jrow)   = tbar(k,n,jrow) + temp3(i,k)
              travar(k,n,jrow) = travar(k,n,jrow) + temp1(i,k)
              dtabs(k,n,jrow)  = dtabs(k,n,jrow) + temp2(i,k)
            enddo
          enddo
        enddo
# ifdef timing
        call toc ('diagnostic', 'time step monitor')
# endif
      endif
#endif
c
#ifdef tracer_averages
c
c-----------------------------------------------------------------------
c     diagnostic: accumulate tracers for averages under horizontal
c                 regions (use units of meters, rather than cm) 
c
c     author:   K. Dixon       e-mail kd@gfdl.gov
c-----------------------------------------------------------------------
c
      if (tavgts .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'tracer averages')
# endif
        do j=js,je
          jrow = j + joff
          do i=is,ie
            mask = mskhr(i,jrow)
            if (mask .ne. 0) then
              boxar = cst(jrow)*dxt(i)*dyt(jrow)*tmask(i,1,j)*0.0001
              sumbf(mask,n) = sumbf(mask,n) + stf(i,j,n)*boxar
              do k=1,km
                sumbk(mask,k,n) = sumbk(mask,k,n) + t(i,k,j,n,tau)
     &                             *boxar*dzt(k)*tmask(i,k,j)*0.01
             enddo
            endif
          enddo
        enddo
# ifdef timing
        call toc ('diagnostic', 'tracer averages')
# endif
      endif
#endif
c
#ifdef meridional_tracer_budget
c
c----------------------------------------------------------------------
c     diagnostic: integrate equations in depth, lon, and time for 
c                 each basin and jrow. basin # 0 is used to catch
c                 values in land areas
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c----------------------------------------------------------------------
c
      if (tmbint .ge. c0 .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'meridional tracer budget')
# endif
        do j=js,je
          jrow = j + joff
          if (jrow .eq. jmtm1 .and. n .eq. 1) numtmb = numtmb + 1
          cosdyt = cst(jrow)*dyt(jrow)
	  do k=1,km
	    do i=is,ie
              stor(i,k) = 0.0
              div(i,k)  = 0.0
              sorc(i,k) = 0.0
              dif(i,k)  = 0.0
            enddo
	  enddo
          do k=1,km
            rtwodt = c1/twodt(k)
            do i=is,ie
              dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
              dxdydz          = dxdy*dzt(k) 
              stor(i,k) = stor(i,k) + rtwodt*dxdydz*
     &                            (t(i,k,j,n,taup1) - t(i,k,j,n,taum1))
              div(i,k)  = div(i,k)  - dxdydz*ADV_Ty(i,k,j)
# ifdef source_term
              sorc(i,k) = sorc(i,k) + dxdydz*source(i,k,j)
# endif
              dif(i,k)  = dif(i,k)  + dxdydz*DIFF_Ty(i,k,j)
            enddo
          enddo
c
c         the accumulation is done this way for issues of speed
c
          do i=is,ie
            t1(i) = stor(i,1)
            t2(i) = div(i,1)
            t3(i) = sorc(i,1)
            t4(i) = dif(i,1)
          enddo
c
	  do k=2,km
	    do i=is,ie
              t1(i) = t1(i) + stor(i,k)
              t2(i) = t2(i)  + div(i,k)
              t3(i) = t3(i) + sorc(i,k)
              t4(i) = t4(i)  + dif(i,k)
            enddo
	  enddo
c
	  do i=is,ie
            m               = msktmb(i,jrow)
            tstor(jrow,n,m) = tstor(jrow,n,m) + t1(i)
            tdiv(jrow,n,m)  = tdiv(jrow,n,m)  + t2(i)
            tsorc(jrow,n,m) = tsorc(jrow,n,m) + t3(i)
            tdif(jrow,n,m)  = tdif(jrow,n,m)  + t4(i)
          enddo
c
	  k = 1
	  do i=is,ie
            m               = msktmb(i,jrow)
            dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
            tflux(jrow,n,m) = tflux(jrow,n,m) + dxdy*stf(i,j,n)
          enddo
          if (n .eq. 1) then
	    do k=1,km
	      do i=is,ie
                m               = msktmb(i,jrow)
                dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
                dxdydz          = dxdy*dzt(k) 
                smdvol(jrow,m)  = smdvol(jrow,m)  + dxdydz
              enddo
	    enddo
          endif
        enddo
# ifdef timing
        call toc ('diagnostic', 'meridional tracer budget')
# endif
      endif
#endif
c
#ifdef gyre_components
c
c-----------------------------------------------------------------------
c     diagnostic: compute the northward transport components of
c                 each tracer
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c               (based on diagnostic by K. Bryan)
c-----------------------------------------------------------------------
c
      if (gyrets .and. eots)  call gyre (joff, js, je, is, ie, n)
#endif
c
#ifdef term_balances
c
c-----------------------------------------------------------------------
c     diagnostic: integrate r.h.s. terms in the tracer equations  
c                 over specified regional volumes.
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (trmbts .and. eots)  call ttb1 (joff, js, je, is, ie, n)
#endif
c
#ifdef xbts
c
c-----------------------------------------------------------------------
c     diagnostic: accumulate r.h.s. terms in the tracer equations
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (eots) call txbt1 (joff, js, je, n)
#endif
#ifdef timing
      call toc ('tracer', 'diagnostics')
#endif
      return
      end




      subroutine diagt2 (joff, js, je, is, ie, idiag)
c
c-----------------------------------------------------------------------
c     construct d(tracer)/dt diagnostics 
c
c     input:
c       joff  = offset relating "j" in the MW to latitude "jrow"
c       js    = starting row in the MW
c       je    = ending row in the MW
c       is    = starting longitude index in the MW
c       ie    = ending longitude index in the MW
c       idiag = 1  => total tracer change
c       idiag = 10 => change of tracer due to filtering(also convection)
c-----------------------------------------------------------------------
c
#include "param.h"
#include "coord.h"
#include "diaga.h"
#include "iounit.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"
c
#ifdef timing
      call ticr ('tracer', 'diagnostics')
#endif
c
#if defined save_convection
c
c-----------------------------------------------------------------------
c     diagnostic: save tracer time change due to explicit convection
c     idiag = 10 signifies diagnostics immediately after convection
c
c     author:   S. Griffies       e-mail smg@gfdl.gov
c-----------------------------------------------------------------------
c
      if (exconvts .and. idiag .eq. 10 .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'explicit convection')
# endif
        rdt = c1/c2dtts
c
c       excnv1 = epsln over land cells and d(convect)/dt over ocean
c
        do j=js,je
          do k=1,km
            do i=1,imt
              excnv1(i,k,j) = tmask(i,k,j)*
     &         (t(i,k,j,1,taup1)-excnv0(i,k,j))*rdt
     &                       + (c1-tmask(i,k,j))*epsln
            enddo
          enddo
        enddo
c
        if (joff .eq. 0) then
          write (stdout,*) ' =>Writing explicit convection at ts=',itt
     &  , ' ',stamp
	  call getunit (iocv, 'cvct.dta','u s a ieee')
c
          period = 0.0
	  reltim = relyr
          iotext = 'read(iocv) reltim, period, imt, jmt, km, flag'
          write (iocv) stamp, iotext, expnam
          write (iocv) reltim, period, imt, jmt, km, epsln
c
          iotext = 'read(iocv) (xt(i),i=1,imt)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, xt, imt)
c
          iotext = 'read(iocv) (yt(j),j=1,jmt)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, yt, jmt)
c
          iotext = 'read(iocv) (zt(k),k=1,km)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, zt, km)
c
          call relunit (iocv)
        endif
c
        call getunit (iocv, 'cvct.dta','u s a ieee')   
c
        do j=js,je
          jrow = j+joff
c
          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocv) (convect(i,k),i=1,imt),k=1,km)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, excnv1(1,1,j), imt*km)
c
        enddo
        call relunit(iocv)
# ifdef timing
        call toc ('diagnostic', 'explicit convection')
# endif
      endif
#endif
c
#ifdef term_balances
c
c-----------------------------------------------------------------------
c     diagnostic: integrate d/dt(tracer) over specified regional volumes 
c                  after convection and filtering
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (trmbts .and. eots) call ttb2 (joff, js, je, is, ie, idiag)
#endif
c
#ifdef xbts
c
c-----------------------------------------------------------------------
c     diagnostic: accumulate d/dt(tracer) after convection
c                 and filtering
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (eots) call txbt2 (joff, js, je, idiag)
#endif
#ifdef timing
      call toc ('tracer', 'diagnostics')
#endif
      return
      end



#if !defined simple_sbc
      subroutine asbct (joff, js, je, is, ie, isst, isss)
c
c-----------------------------------------------------------------------
c     construct the Atmos S.B.C. (surface boundary conditions) SST
c     and possibly SSS
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       isst   = index for sst 
c       isss   = index for sss 
c
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
      call tic ('tracer', 'constructing SBC')
# endif
c
c     initialize the Atmos S.B.C. at the start of each ocean segment
c     (do not alter values in land)
c
      if (eots .and. osegs) then
        do j=js,je
          jrow   = j + joff
          if (isst .ne. 0) then
            do i=is,ie
              if (kmt(i,jrow) .ne. 0) sbcocn(i,jrow,isst) = c0
            enddo
          endif
          if (isss .ne. 0) then
            do i=is,ie
              if (kmt(i,jrow) .ne. 0) sbcocn(i,jrow,isss) = c0
            enddo
          endif
        enddo
      endif
c
c     accumulate surface tracers for the Atmos S.B.C. every time step
c
      if (eots) then
        do j=js,je
          jrow   = j + joff
	  if (isst .ne. 0) then
            do i=is,ie
              sbcocn(i,jrow,isst) = sbcocn(i,jrow,isst)+t(i,1,j,1,taup1)
            enddo
          endif
          if (isss .ne. 0) then
            do i=is,ie
              sbcocn(i,jrow,isss) = sbcocn(i,jrow,isss)+t(i,1,j,2,taup1)
            enddo
	  endif
        enddo
      endif
c
c     average the surface tracers for the Atmos S.B.C. at the end of
c     each ocean segment. (do not alter values in land)
c
      if (eots .and. osege) then
        rts = c1/ntspos
        do j=js,je
          jrow   = j + joff
          if (isst .ne. 0) then
            do i=is,ie
              if (kmt(i,jrow) .ne. 0) sbcocn(i,jrow,isst) = 
     &                                        rts*sbcocn(i,jrow,isst)
            enddo
          endif
          if (isss .ne. 0) then
            do i=is,ie
              if (kmt(i,jrow) .ne. 0) sbcocn(i,jrow,isss) = 
     &	                                      rts*sbcocn(i,jrow,isss)
            enddo
          endif
        enddo
      endif
# ifdef timing
      call toc ('tracer', 'constructing SBC')
# endif
      return
      end
#endif


#if defined implicitvmix || defined isopycmix
      subroutine ivdift (joff, js, je, is, ie, n, twodt)
c
c-----------------------------------------------------------------------
c     solve vertical diffusion of tracers implicitly
c
c     input:
c       joff  = offset relating "j" in the MW to latitude "jrow"
c       js    = starting row in the MW
c       je    = ending row in the MW
c       is    = starting longitude index in the MW
c       ie    = ending longitude index in the MW
c       n     = tracer component
c       twodt = (2*dtts, dtts) on (leapfrog, mixing) time steps
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
# include "size.h"
# include "levind.h"
# include "mw.h"
# include "switch.h"
# include "vmixc.h"
      dimension twodt(km)
c
c     store terms to compute implicit vertical mixing on
c     diagnostic time steps
c
# ifdef timing
      call tic ('tracer', 'implicit vert diff')
# endif
c
c     set some constants
c
      c1 = 1.0
c
# ifdef xbts
      if (eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = t(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif
# else
#  ifdef term_balances
      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = t(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif
#  endif
# endif
c
# if !defined constvmix || defined isopycmix
#  ifdef tcvmix
      call invtri (tp1(1,1,1,n), stf(1,jsmw,n), btf(1,jsmw,n), vdca(???)
     &,              twodt, kmt, tmask(1,1,1), is, ie, joff, js, je)
      if (imt-1 .gt. 1) then
        print *,'Error: ivdif.F is not converted for ifvef tvcmi x'
	stop '=>ivdif'
      endif
#  else
      call invtri (t(1,1,1,n,taup1), stf(1,jsmw,n), btf(1,jsmw,n)
     &, diff_cbt(1,1,jsmw), twodt, kmt, tmask(1,1,1), is, ie
     &, joff, js, je)
#  endif
# endif
c
c     compute residual implicit vertical mixing
c
# ifdef xbts
      if (eots) then
        do j=js,je
          do k=1,km
            rc2dt = c1/twodt(k)
            do i=is,ie
              zzi(i,k,j) = rc2dt*(t(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
	  enddo
	enddo
      endif
# else
#  ifdef term_balances
      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            rc2dt = c1/twodt(k)
            do i=is,ie
              zzi(i,k,j) = rc2dt*(t(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
	  enddo
	enddo
      endif
#  endif
# endif
# ifdef timing
      call toc ('tracer', 'implicit vert diff')
# endif
      return
      end
#endif



#ifdef biharmonic
      subroutine delsqt (joff, js, je, is, ie, n)
c
c=======================================================================
c     compute del**2 of prognostic variables on MW rows "js" ... "je"
c     and gradients of del**2 quantities for variable "n" in the
c     tracer equation
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       n    = tracer component
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#include "grdvar.h"
#include "hmixc.h"
#include "mw.h"
#include "fdift.h"
c
# ifdef timing
        call tic ('tracer', 'biharmonic mixing')
# endif
c
c       set index to skip over the two components of velocity in "del2"
c
        m = n + 2
c
c       set -del**2 = 0 on southern wall
c
        j    = js
        jrow = j + joff
        if (jrow .eq. 2) then
          do k=1,km
            do i=is-1,ie+1
	      del2(i,k,j-1,m) = c0
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
	      del2(i,k,j+1,m) = c0
            enddo
          enddo
	endif
c
c-----------------------------------------------------------------------
c       diffusive flux across eastern face of "T" cell
c       diffusive flux across northern face of "T" cell
c-----------------------------------------------------------------------
c
        do j=js,je
          jrow = j + joff
	  ahbi_cstr = diff_cet*cstr(jrow)
          do k=1,km
            do i=is-1,ie              
              diff_fe(i,k,j) = ahbi_cstr*dxur(i)*
     &                         (t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
        do j=js-1,je
          jrow = j + joff
	  ahbi_csu_dyur = diff_cnt*csu(jrow)*dyur(jrow)
          do k=1,km
            do i=is,ie
              diff_fn(i,k,j) = ahbi_csu_dyur* 
     &                     (t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
c
c       compute -ah*del**2 of tracer
c
        do j=js,je
          jrow   = j + joff
          do k=1,km
            do i=is,ie
	      del2(i,k,j,m) = -DIFF_Tx(i,k,j) - DIFF_Ty(i,k,j)
            enddo
          enddo
	  call setbcx (del2(1,1,j,m), imt, km)
        enddo
c
c-----------------------------------------------------------------------
c       diffusive flux across eastern face of "T" cell
c       diffusive flux across northern face of "T" cell
c-----------------------------------------------------------------------
c
        do j=js,je-1
          jrow = j + joff
	  ahbi_cstr = diff_cet*cstr(jrow)
          do k=1,km
            do i=is-1,ie
              diff_fe(i,k,j) = ahbi_cstr*dxur(i)*
     &                         (del2(i+1,k,j,m)-del2(i,k,j,m))        
            enddo
          enddo
        enddo
        do j=js-1,je-1
          jrow = j + joff
	  ahbi_csu_dyur = diff_cnt*csu(jrow)*dyur(jrow)
          do k=1,km
            do i=is,ie
              diff_fn(i,k,j) = ahbi_csu_dyur*
     &                        (del2(i,k,j+1,m) - del2(i,k,j,m))
            enddo
          enddo
        enddo
c
# ifdef timing
        call toc ('tracer', 'biharmonic mixing')
# endif
        return
	end
#endif


#ifdef sponges
      subroutine sponge1 (joff, js, je, is, ie, n, tm1, source)
c
c=======================================================================
c     newtonian damping variables for sponge regions adjacent to
c     artificial southern and northern boundaries in limited domain
c     basins. data must be prepared using the "sponge" routines 
c     included in the programs for working with the MOM dataset.
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       n    = (1,2) = tracer component
c       tm1  = tracer at "tau-1"
c
c     output:
c       source = newtonian damping term
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#include "iounit.h"
#include "sponge.h"
#include "switch.h"
#include "tmngr.h"
c
      dimension source(imt,km,jsmw:jemw)
      dimension tm1(imt,km,jmw,nt)
c
# ifdef timing
        call tic ('tracer', 'sponge boundary')
# endif
      if (joff .eq. 0 .and. n .eq. 1 .and. is .eq. 2) then
c
c-----------------------------------------------------------------------
c       decide whether to read sponge data or not
c-----------------------------------------------------------------------
c
        begtim = (realdays(initial) - 1.0) + realdays(imodeltime)
        method = 3
        call timeinterp (begtim, indxsp, tspng, spgdpm, 12, .true.
     &,           method, inextd, iprevd, wprev, readsp, inext, iprev)
c
c       read in the next data record from disk when needed
c
        if (readsp) then
          call getunit (ionew2, 'sponges', opt_sponge) 
          read (ionew2, rec=inextd) stnext, spdpmn, im, kk, jm, spngs
     &,     spngn, spbuf1 
          do in=1,4
            do k=1,km
              do i=1,imt
                spbuf(i,k,in,inext) = spbuf1(i,k,in)
              enddo
           enddo
          enddo
          write (stdout,'(/a,i3,a,i2,a,i2,a,a/)')
     &    '=> read sponge record =',inextd,' into buffer =', inext
     &,   ' method #',method,' at ', stamp
          call relunit (ionew2)
        endif
      endif
      if (n .le. 2) then
c
c-----------------------------------------------------------------------
c       construct newtonian damping term using sponge data
c-----------------------------------------------------------------------
c
	do j=js,je
	  jrow = joff + j
	  if (spngs(jrow) .ne. c0) then
	    tnext = c1-wprev
	    do k=1,km
	      do i=is,ie
	        data = tnext*spbuf(i,k,n,inext) 
     &               + wprev*spbuf(i,k,n,iprev)
                source(i,k,j) = source(i,k,j) - 
     &                            spngs(jrow)*(tm1(i,k,j,n) - data)
              enddo
	    enddo
	  endif
c
	  if (spngn(jrow) .ne. c0) then
	    tnext = c1-wprev
	    do k=1,km
	      do i=is,ie
	        data = tnext*spbuf(i,k,n+2,inext) 
     &               + wprev*spbuf(i,k,n+2,iprev)
                source(i,k,j) = source(i,k,j) - 
     &                            spngn(jrow)*(tm1(i,k,j,n) - data)
              enddo
	    enddo
	  endif
        enddo
      endif
# ifdef timing
        call toc ('tracer', 'sponge boundary')
# endif
      return
      end
#endif

#if defined shortwave && !defined simple_sbc
      subroutine swflux0 (joff, js, je, is, ie, source)
c
c=======================================================================
c       incorporate short wave penetration via the "source" 
c       term. note that the divergence of shortwave for the first
c       level "divpen(1)" is compensating for the effect of having
c       the shortwave component already included in the total
c       surface tracer flux "stf(i,j,1)"
c
c       incorporating the penetrative shortwave radiative flux into
c       the vertical diffusive flux term directly is not correct when
c       vertical diffusion is solved implicitly.
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c
c     output:
c       source = source term penetrative short wave
c
c     author:   r.c.pacanowski       e-mail rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#include "csbc.h"      
#include "cshort.h"      
c
      dimension source(imt,km,jsmw:jemw)
# ifdef timing
        call tic ('tracer', 'shortwave penetration')
# endif
c
      isw = mapsbc(5)
      if (isw .ne. 0) then
        do j=js,je
          jrow   = j + joff
	  do k=1,km
	    do i=is,ie
              source(i,k,j) = source(i,k,j)
     &                       + sbcocn(i,jrow,isw)*divpen(k) 
            enddo
          enddo
        enddo
      endif
# ifdef timing
        call toc ('tracer', 'shortwave penetration')
# endif
      return
      end
#endif

