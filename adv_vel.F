      subroutine adv_vel (joff, js, je, is, ie)
c
c=======================================================================
c     calculate advection velocities for momentum and tracer equations     
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c
c     output:
c       adv_vet = advection velocity on east face of "t" cell
c       adv_vnt = advection velocity on north face of "t" cell
c       adv_vbt = advection velocity on bottom face of "t" cell
c       adv_veu = advection velocity on east face of "u" cell
c       adv_vnu = advection velocity on north face of "u" cell
c       adv_vbu = advection velocity on bottom face of "u" cell
c
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#include "coord.h"
#include "grdvar.h"
#include "levind.h"
#include "mw.h"
#ifdef implicit_free_surface
# include "emode.h"
# include "scalar.h"
# include "switch.h"
#endif
# ifdef timing
      call tic ('ocean', 'advection velocity')
# endif
c
c-----------------------------------------------------------------------
c     limit the longitude indices
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
c-----------------------------------------------------------------------
c     advection velocity on northern face of "T" cells. Note the
c     imbedded cosine.
c     adv_vnt = WT_AVG_X(u(1,1,1,2,tau))
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        do k=1,km
          do i=istrt,iend
            adv_vnt(i,k,j) = (u(i,k,j,2,tau)*dxu(i) +   
     &                     u(i-1,k,j,2,tau)*dxu(i-1))*csu(jrow)*dxt2r(i)
          enddo
        enddo
        call setbcx (adv_vnt(1,1,j), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     advection velocity on the eastern face of "T" cells
c     adv_vnt = WT_AVG_Y(u(1,1,1,1,tau))
c-----------------------------------------------------------------------
c
      jstbe = max(js,jsmw)
      do j=jstbe,je
        jrow = j + joff
        do k=1,km
          do i=istrt-1,iend+1
            adv_vet(i,k,j) = (u(i,k,j,1,tau)*dyu(jrow) + 
     &                     u(i,k,j-1,1,tau)*dyu(jrow-1))*dyt2r(jrow)
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     construct vertical velocity on the bottom face of "T" cells
c-----------------------------------------------------------------------
c
      do j=jstbe,je
        jrow = j + joff
c
c       set "adv_vbt" at surface to 0.0 (rigid-lid) or dh/dt (free surf)
c
        do i=istrt,iend
#ifdef implicit_free_surface
	  if (euler2) then
            adv_vbt(i,0,j) = (pguess(i,jrow) - ps(i,jrow,2))/(grav*dtsf)
          else
            adv_vbt(i,0,j) = (ps(i,jrow,1) - ps(i,jrow,2))/(grav*dtsf)
	  endif 
#else
          adv_vbt(i,0,j)   = c0
#endif
        enddo     
c
c       construct divergence of advection velocity * level thickness
c
        do k=1,km
          do i=istrt,iend
            adv_vbt(i,k,j) =
     &                   ((adv_vet(i,k,j) - adv_vet(i-1,k,j))*dxtr(i)
     &                   +(adv_vnt(i,k,j) - adv_vnt(i,k,j-1))*dytr(jrow)
     &                   )*cstr(jrow)*dzt(k)
          enddo
        enddo
c
c       integrate downward to define "adv_vbt" at the bottom of levels
c
        do k=1,km
          do i=istrt,iend
            adv_vbt(i,k,j) = adv_vbt(i,k,j) + adv_vbt(i,k-1,j)
          enddo
        enddo
c
	call setbcx (adv_vbt(1,0,j), imt, km+1)
c
      enddo
c
#ifdef linearized_advection
c
c-----------------------------------------------------------------------
c     Advective velocities for U cells are to remain zero. Only the
c     vertical advective velocity on T cells will be retained
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=1,km
          do i=istrt-1,iend+1
            adv_vnt(i,k,j) = c0
          enddo
        enddo
      enddo
      do j=jstbe,je
        do k=1,km
          do i=istrt-1,iend+1
            adv_vet(i,k,j) = c0
          enddo
        enddo
      enddo
      
#endif
c
c-----------------------------------------------------------------------
c     construct advection velocity on the northern face of "u" cells by
c     averaging advection velocity on northern face of "t" cells 
c     note: je-1 is used instead of jemw to account for possible non
c           integral number of MW`s in jmt 
c     adv_vnu = LINEAR_INTRP_Y(WT_AVG_X(adv_vnt))
c-----------------------------------------------------------------------
c
      jsun = max(js,jsmw)-1
      do j=jsun,je-1
        jrow = j + joff
	dyr  = dytr(jrow+1)
        do k=1,km
          do i=istrt,iend
#ifdef linearized_advection
	    adv_vnu(i,k,j) = c0 
#else
	    adv_vnu(i,k,j) = ((adv_vnt(i,k,j)*duw(i) 
     &                       + adv_vnt(i+1,k,j)*due(i)
     &                        )*dus(jrow+1) +
     &                        (adv_vnt(i,k,j+1)*duw(i) 
     &                       + adv_vnt(i+1,k,j+1)*due(i)
     &                        )*dun(jrow))*dyr*dxur(i)
#endif
          enddo
        enddo
	call setbcx (adv_vnu(1,1,j), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     construct advection velocity on the eastern face of "u" cells by
c     averaging advection velocity on eastern face of "t" cells 
c     note: take special care of zonal b.c. on this term.
c     adv_veu = LINEAR_INTRP_X(WT_AVG_Y(adv_vet))
c-----------------------------------------------------------------------
c
      jsube = max(js-1,jsmw)
      do j=jsube,je-1
        jrow = j + joff
	dyr  = dyur(jrow)
        do k=1,km
          do i=istrt-1,iend
#ifdef linearized_advection
	    adv_veu(i,k,j) = c0
#else
	    adv_veu(i,k,j) = ((adv_vet(i,k,j)*dus(jrow)
     &                       + adv_vet(i,k,j+1)*dun(jrow)
     &                        )*duw(i+1) +    
     &                        (adv_vet(i+1,k,j)*dus(jrow)
     &                       + adv_vet(i+1,k,j+1)*dun(jrow)
     &                        )*due(i))*dyr*dxtr(i+1)
#endif
          enddo
        enddo
#ifdef cyclic
	call setbcx (adv_veu(1,1,j), imt, km)
#else
        do k=1,km
	  adv_veu(imt,k,j) = c0
	enddo
#endif
      enddo
c
c-----------------------------------------------------------------------
c     construct advection velocity on the bottom face of "u" cells by
c     averaging advection velocity on bottom face of "t" cells 
c-----------------------------------------------------------------------
c
      do j=jsube,je-1
        jrow = j + joff
	dyn  = dun(jrow)*cst(jrow+1)
	dys  = dus(jrow)*cst(jrow)
	dyr  = dyur(jrow)*csur(jrow)
        do k=0,km
          do i=istrt,iend
	    asw = duw(i)*dys
	    anw = duw(i)*dyn
	    ase = due(i)*dys
	    ane = due(i)*dyn
#ifdef linearized_advection
	    adv_vbu(i,k,j) = c0
#else
	    adv_vbu(i,k,j) = dyr*dxur(i)*(
     &                    adv_vbt(i,k,j)*asw + adv_vbt(i+1,k,j)*ase
     &                  + adv_vbt(i,k,j+1)*anw + adv_vbt(i+1,k,j+1)*ane)
#endif
          enddo
        enddo
c
	call setbcx (adv_vbu(1,0,j), imt, km+1)
c
      enddo
# ifdef trajectories
c
c     insure advection velocity on bottom face of 1st latitude cell = 0
c
      if (joff .eq. 0) then
        j = 1
        do k=0,km
          do i=istrt-1,iend+1
	    adv_vbu(i,k,j) = c0
          enddo
        enddo
      endif
# endif
# ifdef trace_indices
      write (stdout,'(2x,3(a,i4),6(/4x,4(a,i4),a))')
     & "=> In adv_vel: js=",js," je=",je," joff=",joff
     &, " adv_vnt calculated on jrows ",js+joff," through ",je+joff
     &, " (rows ",js," through ",je," in the MW)"
     &, " adv_vet calculated on jrows ",jstbe+joff," through ",je+joff
     &, " (rows ",jstbe," through ",je," in the MW)"
     &, " adv_vbt calculated on jrows ",jstbe+joff," through ",je+joff
     &, " (rows ",jstbe," through ",je," in the MW)"
     &, " adv_vnu calculated on jrows ",jsun+joff," through ",je-1+joff
     &, " (rows ",jsun," through ",je-1," in the MW)"
     &, " adv_veu calculated on jrows ",jsube+joff," through ",je-1+joff
     &, " (rows ",jsube," through ",je-1," in the MW)"
     &, " adv_vbu calculated on jrows ",jsube+joff," through ",je-1+joff
     &, " (rows ",jsube," through ",je-1," in the MW)"
# endif
# ifdef timing
      call toc ('ocean', 'advection velocity')
# endif
c
c#define debug_adv_vel
#ifdef debug_adv_vel
c
c-----------------------------------------------------------------------
c     to inspect the divergence of "t" and "u" cells at point (ipt,jpt)
c     for all levels
c-----------------------------------------------------------------------
c
      sml  = 0
      ipt = 47
      jpt = 43
      do j=jstrt,je-1
        jrow = j + joff
        if (jrow .eq. jpt) then
	  do i=ipt,ipt
          do k=1,km
	    divgt =
     &       ((adv_vet(i,k,j) - adv_vet(i-1,k,j))*cstr(jrow)*dxtr(i)
     & +      (adv_vnt(i,k,j) - adv_vnt(i,k,j-1))*cstr(jrow)*dytr(jrow)
     & +      (adv_vbt(i,k-1,j) - adv_vbt(i,k,j))*dztr(k))*tmask(i,k,j)
	    divgu =
     &        ((adv_veu(i,k,j) - adv_veu(i-1,k,j))*csur(jrow)*dxur(i)
     & +      (adv_vnu(i,k,j) - adv_vnu(i,k,j-1))*csur(jrow)*dyur(jrow)
     & +      (adv_vbu(i,k-1,j) - adv_vbu(i,k,j))*dztr(k))*umask(i,k,j)
            if (abs(divgt) .ge. sml .or. abs(divgu) .ge. sml) then
              write (stdout,98) i,k,divgt
     &,                         adv_vet(i,k,j), adv_vet(i-1,k,j)
     &,                         adv_vnt(i,k,j), adv_vnt(i,k,j-1)
     &,                         adv_vbt(i,k-1,j), adv_vbt(i,k,j)
              write (stdout,99) i,k,divgu
     &,                         adv_veu(i,k,j), adv_veu(i-1,k,j)
     &,                         adv_vnu(i,k,j), adv_vnu(i,k,j-1)
     &,                         adv_vbu(i,k-1,j), adv_vbu(i,k,j)
              write (stdout,*) ' tmask=',tmask(i,k,j), ', umask='
     &,                        umask(i,k,j),' kmt=',kmt(i,jrow),' kmu='
     &,                        kmu(i,jrow), ', jrow=',jrow
              write (stdout,*) ' '
	    endif
          enddo
	  write (stdout,*) 'adv_vbt(i,0,j) = ',adv_vbt(i,0,j)
	  write (stdout,*) 'adv_vbu(i,0,j) = ',adv_vbu(i,0,j)
	  write (stdout,'(///)')
	  enddo
        endif
      enddo
98      format (1x,'tp: i=',i3,', k=',i2,'divt=',7e14.7)
99      format (1x,'tp: i=',i3,', k=',i2,'divu=',7e14.7)
#endif
      return
      end
      
