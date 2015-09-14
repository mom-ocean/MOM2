      subroutine vmixc (joff, js, je, is, ie)
c
c=======================================================================
c     set viscosity coefficient on bottom face of "u" cells
c     set diffusion coefficient on bottom face of "t" cells
c
c     input:
c       joff = offset relating "j" in the MW to latitude "jrow"
c       js   = starting row in the MW
c       je   = ending row in the MW
c       is   = starting longitude index in the MW
c       ie   = ending longitude index in the MW
c
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#include "coord.h"
#ifdef held_larichev
# include "hmixc.h"
#endif
#include "mw.h"
#include "switch.h"
#include "vmixc.h"
#ifdef isopycmix
# include "isopyc.h"
#endif
c
c-----------------------------------------------------------------------
c     limit the longitude and latitude indices
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
      jstrt = max(js-1,jsmw)
      jend  = je-1
c
c-----------------------------------------------------------------------
c     set all vertical mixing coefficients
c-----------------------------------------------------------------------
c
#if defined constvmix 
c
c     for momentum
c
      if (first) then
        visc_cbu = kappa_m
        do k=1,km
          kappa_m_dzwr(k) = kappa_m/dzw(k)
        enddo
      endif
#endif
c
#if defined constvmix && !defined isopycmix
c
c     for tracers
c
      if (first) then
# ifdef bryan_lewis_vertical
        do k=1,km
	  diff_cbt(k)     = Ahv(k)
          kappa_h_dzwr(k) = Ahv(k)/dzw(k)
        enddo
# else
        diff_cbt = kappa_h
        do k=1,km
          kappa_h_dzwr(k) = kappa_h/dzw(k)
        enddo
# endif
      endif
#endif
c
#if defined constvmix && defined isopycmix
c
c     for tracers
c
      do j=jstrt,jend
        do i=istrt,iend
          do k=1,km
# ifdef bryan_lewis_vertical
	    diff_cbt(i,k,j) = Ahv(k)
# else
            diff_cbt(i,k,j) = kappa_h
# endif
          enddo
        enddo
      enddo
#endif
c
#ifdef ppvmix
c
c-----------------------------------------------------------------------
c     for momentum and tracers based on the pacanowski & philander
c     richardson mixing scheme (JPO vol 11, #11, 1981).
c-----------------------------------------------------------------------
c
      call ppmix (joff, js, je, istrt, iend)
#endif
c
#if defined isopycmix
c
c-----------------------------------------------------------------------
c     add K33 component of isopycnal slope to vertical diffusion coeff
c-----------------------------------------------------------------------
c
      do j=jstrt,jend
        do i=istrt,iend
          do k=1,km-1
# ifdef held_larichev
            ahisop = hl_b(i,j)
# endif
            diff_cbt(i,k,j) = (ahisop*K3(i,k,j,3)
# if defined constvmix && !defined bryan_lewis_vertical
     &                         + kappa_h
#endif
# if defined constvmix && defined bryan_lewis_vertical
     &                         + Ahv(k)
#endif
# if defined ppvmix || defined tcvmix
     &                         + diff_cbt(i,k,j)
#endif
     &                        )*tmask(i,k+1,j)
          enddo
          diff_cbt(i,km,j) = c0
        enddo
      enddo
#endif
#ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     & "=> In vmixc: js=",js," je=",je," joff=",joff
     &," jstrt=",jstrt," jend=",jend
#endif
      return
      end
