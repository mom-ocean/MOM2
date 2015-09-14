      subroutine loadmw (joff, js, je, is, ie, num1, nu)
c
c=======================================================================
c     load variables into the MW for rows "js" through "je"
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       num1 = "tau-1" latitude disk unit
c       nu   = "tau" latitude disk unit
c     
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#include "emode.h"
#include "grdvar.h"
#include "iounit.h"
#include "levind.h"
#include "mw.h"
#include "switch.h"
#include "tmngr.h"
      save taum1old
c
#ifdef timing
      call tic ('loadmw', 'loading the MW')
#endif
c
c-----------------------------------------------------------------------
c     limit the longitude indices
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
c-----------------------------------------------------------------------
c     for all MW`s after the first, move the MW northward by moving
c     data from the last two rows into the first two rows.
c     (last 3 rows into first 3 rows if using the biharmonic option)
c-----------------------------------------------------------------------
c
      if (joff .ne. 0) then
        call movemw (istrt-1, iend+1)
      endif
c
c-----------------------------------------------------------------------
c     construct "t" cell and "u" cell land/sea masks
c-----------------------------------------------------------------------
c
#ifdef timing
      call tic ('loadmw', 'getvar+mask+ext mode')
#endif
      if (.not. wide_open_mw .or. (wide_open_mw .and. first)) then
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt-1,iend+1
              if (kmt(i,jrow) .ge. k) then
                tmask(i,k,j) = c1
              else
                tmask(i,k,j) = c0
              endif
              if (kmu(i,jrow) .ge. k) then
                umask(i,k,j) = c1
              else
                umask(i,k,j) = c0
              endif
            enddo
          enddo
        enddo
      endif
c
c-----------------------------------------------------------------------
c     read data from "tau-1" and "tau" disk units into the MW
c     or if MW is wide open, move data to proper time level and add
c     external mode to "tau" only since "tau-1" was done last timestep
c-----------------------------------------------------------------------
c
      if (wide_open_mw) then
        if (leapfrog) then
          call add_ext_mode (joff, js, je, istrt, iend, 'tau')
	  if (first) then
	    call add_ext_mode (joff, js, je, istrt, iend, 'tau-1')
	  endif
	elseif (euler2) then
          tau   = taup1
	  taup1 = taum1old
          call add_ext_mode (joff, js, je, istrt, iend, 'tau')
	elseif (euler1 .or. forward) then
          taum1old = taum1
	  taum1    = tau
          call add_ext_mode (joff, js, je, istrt, iend, 'tau')
	endif
      else
        if (leapfrog .or. euler2) then
c
c         read "tau","tau-1" disk data into "tau","tau-1" MW positions
c
          call getvar (joff, js, je, istrt-1, iend+1, num1, nu)
	elseif (forward .or. euler1) then
c
c         read "tau" disk data into "tau" and "tau-1" MW positions
c
          call getvar (joff, js, je, istrt-1, iend+1, nu, nu)
	endif
c
c       add external mode to both since only internal modes are on disk
c
        call add_ext_mode (joff, js, je, istrt, iend, 'tau-1')
        call add_ext_mode (joff, js, je, istrt, iend, 'tau')
      endif
#ifdef timing
      call toc ('loadmw', 'getvar+mask+ext mode')
#endif
c
#ifdef symmetry
c
c-----------------------------------------------------------------------
c     set symmetry condition at latitude of "u" cell row jmtm1 
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        if (jrow .eq. jmt) then
          do k=1,km
	    do i=istrt-1,iend+1
	      u(i,k,j,1,tau)     = u(i,k,j-2,1,tau)
	      u(i,k,j,1,taum1)   = u(i,k,j-2,1,taum1)
	      u(i,k,j-1,2,tau)   = c0
	      u(i,k,j-1,2,taum1) = c0
	      u(i,k,j,2,tau)     = -u(i,k,j-2,2,tau)
	      u(i,k,j,2,taum1)   = -u(i,k,j-2,2,taum1)
	    enddo
	  enddo
	  do n=1,nt
            do k=1,km
	      do i=istrt-1,iend+1
	        t(i,k,j,n,tau)   = t(i,k,j-1,n,tau)
	        t(i,k,j,n,taum1) = t(i,k,j-1,n,taum1)
	      enddo
	    enddo
	  enddo
        endif
      enddo
#endif
c
c-----------------------------------------------------------------------
c     compute density  at "t" cell centers
c-----------------------------------------------------------------------
c
#ifdef timing
      call tic ('loadmw', 'state')
#endif
      call state (t(1,1,1,1,tau), t(1,1,1,2,tau), rho(1,1,jsmw)
     &,           max(jsmw,js), je, istrt-1, iend+1)
#ifdef timing
      call toc ('loadmw', 'state')
#endif
c
#ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     & "=> In loadmw: loading rows ",js+joff," through ",je+joff
     &," from disk into rows ",js," through ",je
     &," in the Memory window. joff= ",joff
#endif
#ifdef timing
      call toc ('loadmw', 'loading the MW')
#endif
c
      return
      end      



      subroutine movemw (iss, iee)
c
c=======================================================================
c     move the MW up (northward) by moving data from the last two rows
c     into the first two rows. (last 3 rows if using biharmonic option)
c
c     input:
c      is = starting longitude index in the MW
c      ie = ending longitude index in the MW
c      note: iss,iee bypassed to optimize performance
c     
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#include "mw.h"
#if defined biharmonic && !defined constvmix
# include "vmixc.h"
#endif
      parameter (is=1, ie=imt)
c
#ifdef timing
      call tic ('loadmw', 'moving MW northward')
#endif
c
      nrows = jmw - ncrows
      do move=1,nrows
        jfrom = jmw - (nrows - move)
	jto   = move
c
c-----------------------------------------------------------------------
c       move quantities with rows dimensioned (1:jmw)
c-----------------------------------------------------------------------
c
        do k=1,km
	  do i=is,ie
	    do n=1,nt
	      t(i,k,jto,n,taum1) = t(i,k,jfrom,n,taum1)
	      t(i,k,jto,n,tau)   = t(i,k,jfrom,n,tau)
	    enddo
	    do n=1,2
	      u(i,k,jto,n,taum1) = u(i,k,jfrom,n,taum1)
	      u(i,k,jto,n,tau)   = u(i,k,jfrom,n,tau)
	    enddo
	    tmask(i,k,jto)   = tmask(i,k,jfrom)
	    umask(i,k,jto)   = umask(i,k,jfrom)
	    adv_vnt(i,k,jto) = adv_vnt(i,k,jfrom)
          enddo
        enddo
c
#ifdef biharmonic
        do n=1,nvarbh
          do k=1,km
	    do i=is,ie
	      del2(i,k,jto,n) = del2(i,k,jfrom,n)
            enddo
          enddo
        enddo
#endif
c
c-----------------------------------------------------------------------
c       move quantities with rows dimensioned (1:jemw)
c-----------------------------------------------------------------------
c
        if (jfrom .le. jemw) then
          do k=1,km
	    do i=is,ie
	      adv_vnu(i,k,jto) = adv_vnu(i,k,jfrom)
            enddo
          enddo
# ifdef trajectories
          do k=0,km
	    do i=is,ie
	      adv_vbu(i,k,jto) = adv_vbu(i,k,jfrom)
            enddo
          enddo
# endif
	endif
c
c-----------------------------------------------------------------------
c       move quantities with rows dimensioned (jsmw:jmw)
c-----------------------------------------------------------------------
c
        if (jto .ge. jsmw) then
          do k=1,km
	    do i=is,ie
	      adv_vet(i,k,jto) = adv_vet(i,k,jfrom)
	      rho(i,k,jto)  = rho(i,k,jfrom)
            enddo
          enddo
          do k=0,km
	    do i=is,ie
	      adv_vbt(i,k,jto) = adv_vbt(i,k,jfrom)
            enddo
          enddo
	endif
#ifdef biharmonic
c
c-----------------------------------------------------------------------
c       move quantities with rows dimensioned (jsme:jemw)
c-----------------------------------------------------------------------
c
        if (jto .eq. jsmw) then
          do k=1,km
	    do i=is,ie
#if !defined constvmix
	      visc_cbu(i,k,jto)  = visc_cbu(i,k,jfrom)
	      diff_cbt(i,k,jto)  = diff_cbt(i,k,jfrom)
#endif
	      adv_veu(i,k,jto) = adv_veu(i,k,jfrom)
            enddo
          enddo
          do k=0,km
	    do i=is,ie
	      adv_vbu(i,k,jto) = adv_vbu(i,k,jfrom)
            enddo
          enddo
	endif
#endif
c
# ifdef trace_indices
        write (stdout,'(4x,2(a,i4))')
     & "=> In movemw: moving variables on row ",jfrom," to row ",jto
# endif
      enddo
#ifdef timing
      call toc ('loadmw', 'moving MW northward')
#endif
      return
      end



      subroutine getvar (joff, js, je, is, ie, num1, nu)
c
c=======================================================================
c     read prognostic quantities from disk units "num1" (tau-1) and
c     "nu" (tau) into the MW for rows "js" through "je"
c
c     input:
c       joff = offset between "j" in the MW and latitude "jrow"
c       js = starting row in the MW
c       je = ending row in the MW
c       is = starting longitude index in the MW
c       ie = ending longitude index in the MW
c       num1 = "tau-1" latitude disk unit
c       nu   = "tau" latitude disk unit
c     
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#include "mw.h"
c
#ifdef timing
      call tic ('getvar', 'from disk')
#endif
c
      do j=js,je
        jrow   = j + joff
c
c       read "tau-1" data into the "tau-1" portion of the MW
c
        call getrow (num1, nslab, jrow, u(1,1,j,1,taum1)
     &,                                 t(1,1,j,1,taum1))
        if (nu .ne. num1) then
c
c         read "tau" data into the "tau" portion of the MW 
c
          call getrow (nu, nslab, jrow, u(1,1,j,1,tau)
     &,                                 t(1,1,j,1,tau))
        else
c
c         copy "tau" data into "tau-1" portion of the MW
c
	  do k=1,km
	    do i=is,ie
	      u(i,k,j,1,tau) = u(i,k,j,1,taum1)
	      u(i,k,j,2,tau) = u(i,k,j,2,taum1)
	    enddo
	  enddo
	  do n=1,nt
	    do k=1,km
	      do i=is,ie
	        t(i,k,j,n,tau) = t(i,k,j,n,taum1)
	      enddo
	    enddo
	  enddo
        endif
      enddo
c
#ifdef timing
      call toc ('getvar', 'from disk')
#endif
c
      return
      end


      subroutine add_ext_mode (joff, js, je, is, ie, timelev)
c
c=======================================================================
c     add external mode to velocity for time level "timelev"
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c       timelev = "tau" or "tau-1"
c     
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#include "emode.h"
#include "grdvar.h"
#include "iounit.h"
#include "mw.h"
#include "switch.h"
#ifdef stream_function
      dimension ext(imt,1:jmw,2)
#endif
      character*(*) timelev
c
      if (timelev .eq. 'tau') then
c
c-----------------------------------------------------------------------
c       add external mode to "tau" velocity
c-----------------------------------------------------------------------
c
        do j=js,je
          jrow = j + joff
          if (jrow .lt. jmt) then
#ifdef stream_function
            do i=is,ie
              diag1       = psi(i+1,jrow+1,1) - psi(i  ,jrow,1)
              diag0       = psi(i  ,jrow+1,1) - psi(i+1,jrow,1)
	      ext(i,j,1)  = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
	      ext(i,j,2)  =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)
     &                       *csur(jrow)
            enddo
            do n=1,2
              do k=1,km
                do i=is,ie
                  u(i,k,j,n,tau) = (u(i,k,j,n,tau) + ext(i,j,n))
     &                             *umask(i,k,j)
                enddo
              enddo
            enddo
#endif
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
            do n=1,2
	      do k=1,km
                do i=is,ie
                  u(i,k,j,n,tau) = (u(i,k,j,n,tau) + ubar(i,jrow,n))
     &                             *umask(i,k,j)
                enddo
              enddo
            enddo
#endif
            do n=1,2
	      call setbcx (u(1,1,j,n,tau), imt, km)
	    enddo
	  endif
	enddo
c
      elseif (timelev .eq. 'tau-1') then
c
c-----------------------------------------------------------------------
c       add external mode to "tau-1" velocity
c-----------------------------------------------------------------------
c
        do j=js,je
          jrow = j + joff
          if (jrow .lt. jmt) then
#ifdef stream_function
            do i=is,ie
              diag1       = psi(i+1,jrow+1,2) - psi(i  ,jrow,2)
              diag0       = psi(i  ,jrow+1,2) - psi(i+1,jrow,2)
	      ext(i,j,1)  = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
	      ext(i,j,2)  =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)
     &                       *csur(jrow)
            enddo
            do n=1,2
              do k=1,km
                do i=is,ie
                  u(i,k,j,n,taum1) = (u(i,k,j,n,taum1) + ext(i,j,n))
     &                               *umask(i,k,j)
                enddo
              enddo
            enddo
#endif
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
            do n=1,2
	      do k=1,km
                do i=is,ie
                  u(i,k,j,n,taum1) = (u(i,k,j,n,taum1)+ubarm1(i,jrow,n))
     &                           *umask(i,k,j)
                enddo
              enddo
            enddo
#endif
            do n=1,2
	      call setbcx (u(1,1,j,n,taum1), imt, km)
	    enddo
	  endif
	enddo
      else
        write (stdout,'(/a,a/)') 'Error: timelev = ',timelev
	stop '=>loadmw'
      endif
c
      return
      end
      


      subroutine putmw (joff, js, je, nup1)
c
c=======================================================================
c     write prognostic quantities from MW to disk unit "nup1" (tau+1)
c     for rows "js" to "je".
c     
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#include "mw.h"
c
#ifdef timing
      call tic ('putmw', 'to disk')
#endif
c
c     write all newly computed quantities to disk "tau+1"
c
      do j=js,je
        jrow   = j + joff
        call putrow (nup1, nslab, jrow, u(1,1,j,1,taup1)
     &,                                 t(1,1,j,1,taup1))
      enddo
c
#ifdef timing
      call toc ('putmw', 'to disk')
#endif
c
      return
      end


      
