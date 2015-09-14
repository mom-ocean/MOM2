      subroutine calc_psgrad (psgrad, uext, vext, js, je, is, ie)
c
c-----------------------------------------------------------------------
c     compute the surface pressure gradients
c
c     inputs:
c
c     js   = index of starting row
c     je   = index of ending row
c     is   = index of starting longitude
c     ie   = index of ending longitude
c
c     outputs:
c
c     psgrad  = grad(surf press)
c     uext = external mode u (tau+1) for point (ie,je) only
c     vext = external mode v (tau+1) for point (ie,je) only
c-----------------------------------------------------------------------
c
#include "param.h"
#include "emode.h"
#include "grdvar.h"
#include "levind.h"
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"
      dimension psgrad(is:ie,js:je,2)
c
      js1 = max(js,2)
      je1 = min(je,jmt-1)
      is1 = max(is,2)
      ie1 = min(ie,imtm1)
c
c     on mixing time steps "ptd" has been multiplied by a factor of
c     two and the time step has to be adjusted also.
c
      if (mod(itt,nmix) .eq. 1) then
        fxa = p5
      else
        fxa = c1
      endif
      r2dtuv = c1/c2dtuv
      do jrow=js1,je1
        atosp = acor*c2*omega*sine(jrow)
        f2    = atosp*c2dtuv
c       f3    = c2dtuv/c2dtsf
        f3    = c2dtuv
	do i=is1,ie1
          kz = kmu(i,jrow)
          if (kz .ne. 0) then
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
            uext         = ubar(i,jrow,1)
	    vext         = ubar(i,jrow,2)
            d1           = ps(i+1,jrow+1,1) - ps(i,jrow,1)
            d2           = ps(i+1,jrow,1) - ps(i,jrow+1,1) 
	    psgrad(i,jrow,1)  = (d1 + d2)*dxu2r(i)*csur(jrow)
	    psgrad(i,jrow,2)  = (d1 - d2)*dyu2r(jrow)
#endif
#ifdef stream_function
            diag1        = psi(i+1,jrow+1,1)-psi(i  ,jrow,1)
            diag0        = psi(i  ,jrow+1,1)-psi(i+1,jrow,1)
            uext         = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
            vext         =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)*csur(jrow)
            diag3        = fxa*(ptd(i+1,jrow+1)-ptd(i  ,jrow))
            diag4        = fxa*(ptd(i  ,jrow+1)-ptd(i+1,jrow))
            dubdt        = (diag3+diag4)*dyu2r(jrow)*hr(i,jrow)
            dvbdt        = (diag3-diag4)*dxu2r(i)*hr(i,jrow)*csur(jrow)
            psgrad(i,jrow,1)=r2dtuv*(dubdt + f3*zu(i,jrow,1) + f2*dvbdt)
            psgrad(i,jrow,2)=r2dtuv*(-dvbdt+ f3*zu(i,jrow,2) + f2*dubdt)
#endif
          endif
        enddo
      enddo
      return
      end
