      subroutine hmixc (joff, js, je, is, ie)
c
c=======================================================================
c     set horizontal mixing coeffs on north and east face of "t" and
c     "u" cells.
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
#include "grdvar.h"
#include "hmixc.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
c
c-----------------------------------------------------------------------
c     limit the latitude indices
c-----------------------------------------------------------------------
c
      jstrt = max(1,js-1)
      jend  = je-1
c
c-----------------------------------------------------------------------
c     set all horizontal mixing coefficients
c-----------------------------------------------------------------------
c
#if defined consthmix
c
c     for momentum
c
      if (first) then
# ifdef biharmonic
        visc_cnu = sqrt(abs(ambi))
        visc_ceu = sqrt(abs(ambi))
# else
        visc_cnu  = am
        visc_ceu  = am
# endif
c
        visc_cnu_rho0r = visc_cnu/rho0
        do j=jstrt,jend
	  jrow = j + joff
	  jm1 = max(1,jrow-1)
	  jp1 = min(jmt,jrow+1)
	  amc_north(jrow) = visc_cnu_rho0r*cst(jp1)*dytr(jp1)
     &                              *csur(jrow)*dytr(jrow)
	  amc_south(jrow) = visc_cnu_rho0r*cst(jrow)*dytr(jrow)
     &                              *csur(jrow)*dytr(jrow)
        enddo
      endif
#endif
c
#if defined consthmix
c
c     for tracers
c
      if (first) then
# ifdef bryan_lewis_horizontal
        do k=1,km
          diff_cet(k) = Ahh(k)
          diff_cnt(k) = Ahh(k)
        enddo
        do j=jstrt,jend
	  jrow = j + joff
	  jm1 = max(1,jrow-1)
	  jp1 = min(jmt,jrow+1)
          do k=1,km
	    ahc_north(jrow,k) = diff_cnt(k)*csu(jrow)*dyur(jrow)
     &                                     *cstr(jrow)*dytr(jrow)
	    ahc_south(jrow,k) = diff_cnt(k)*csu(jm1)*dyur(jm1)
     &                                     *cstr(jrow)*dytr(jrow)
          enddo
        enddo
# else
#  ifdef biharmonic
        diff_cet = sqrt(abs(ahbi))
        diff_cnt = sqrt(abs(ahbi))
#  else
        diff_cnt  = ah
	diff_cet  = ah
#  endif
c
        do j=jstrt,jend
	  jrow = j + joff
	  jm1 = max(1,jrow-1)
	  jp1 = min(jmt,jrow+1)
          ahc_north(jrow) = diff_cnt*csu(jrow)*dyur(jrow)*cstr(jrow)
     &                              *dytr(jrow)
          ahc_south(jrow) = diff_cnt*csu(jm1)*dyur(jm1)*cstr(jrow)
     &                              *dytr(jrow)
        enddo
# endif
      endif
#endif
c
#ifdef smagnlmix
c
c-----------------------------------------------------------------------
c     limit the longitude indices and calculate deformation rates
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
      call smagnlc (joff, js, je, istrt, iend)
#endif
c
#ifdef held_larichev
c
c-----------------------------------------------------------------------
c     calculate Tracer mixing coefficient based on held_larichev
c     (Note: this is a hybrid scheme used only for tracers)
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
      call hlmix (joff, js, je, istrt, iend)
#endif
#ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     & "=> In hmixc: js=",js," je=",je," joff=",joff," jstrt=",jstrt
     &," jend=",jend
#endif
c
      return
      end
