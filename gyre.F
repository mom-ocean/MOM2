      subroutine gyre (joff, js, je, is, ie, n)
c
c-----------------------------------------------------------------------
c     compute the northward transport components of each tracer
c
c     input:
c       joff  = offset relating "j" in the MW to latitude "jrow"
c       js    = starting row in the MW
c       je    = ending row in the MW
c       is    = starting longitude index in the MW
c       ie    = ending longitude index in the MW
c       n     = tracer component
c
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c     based on work by M.Cox and B. Bryan
c     modifications for basins: gokhan danabasoglu
c-----------------------------------------------------------------------
c
#ifdef gyre_components
# include "param.h"
# include "coord.h"
# include "cregin.h"
# include "diag.h"
# include "grdvar.h"
# include "hmixc.h"
# ifdef isopycmix
#  include "isopyc.h"
# endif
# include "mw.h"
# include "scalar.h"
c
# ifdef timing
      call tic ('diagnostic', 'gyre components')
# endif
c
      do j=js,je
        jrow = j + joff
        if (jrow .lt. jmtm1) then
	  small = 1.e-10
          do k=1,km
            totdxn = small
            totdxs = small
            vbr    = c0
            tbrs   = c0
            tbrn   = c0
            do i=is,ie
              totdxn = totdxn + dxt(i)*tmask(i,k,j+1)
              totdxs = totdxs + dxt(i)*tmask(i,k,j)
              vbr    = vbr  + u(i,k,j,2,tau)*dxu(i)*csu(jrow)
              tbrn   = tbrn + t(i,k,j+1,n,tau)*tmask(i,k,j+1)*dxt(i)
              tbrs   = tbrs + t(i,k,j,n,tau)*tmask(i,k,j)*dxt(i)
	    enddo
            tbrn          = tbrn/totdxn
            tbrs          = tbrs/totdxs
            ttn(1,jrow,n) = ttn(1,jrow,n) + vbr*p5*(tbrn+tbrs)*dzt(k)
            do i=is,ie
	      tempdiff_fn = 
#if defined consthmix && !defined biharmonic && !defined isopycmix
# ifdef bryan_lewis_horizontal
     &                      ahc_north(jrow,k)*
     &                      (t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))*
# else
     &                      ahc_north(jrow)*
     &                      (t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))*
# endif
#else
     &                      diff_fn(i,k,j)*
#endif
     &                           tmask(i,k,j+1)*tmask(i,k,j)*
     &                           dxt(i)*dzt(k)*csu(jrow)
	      tempadv_fn       = p5*adv_vnt(i,k,j)*(t(i,k,j,n,tau) +
     &                           t(i,k,j+1,n,tau))*dxt(i)*dzt(k)
              ttn(6,jrow,n)    = ttn(6,jrow,n) + tempadv_fn
              ttn(7,jrow,n)    = ttn(7,jrow,n) - tempdiff_fn
              ttn2(6,jrow,n,0) = ttn2(6,jrow,n,0) + tempadv_fn
              ttn2(7,jrow,n,0) = ttn2(7,jrow,n,0) - tempdiff_fn
              if (mskhr(i,jrow) .ne. 0) then
                ttn2(6,jrow,n,mskhr(i,jrow)) = 
     &                   ttn2(6,jrow,n,mskhr(i,jrow)) + tempadv_fn
                ttn2(7,jrow,n,mskhr(i,jrow)) = 
     &                   ttn2(7,jrow,n,mskhr(i,jrow)) - tempdiff_fn
              endif
# if defined isopycmix && defined gent_mcwilliams
              fxb = adv_vntiso(i,k,j)*
     &              (t(i,k,j,n,taum1)+t(i,k,j+1,n,taum1))
     &              *p5*dxt(i)*dzt(k)*csu(jrow)
              ttn2(9,jrow,n,0) = ttn2(9,jrow,n,0) + fxb
              if (mskhr(i,jrow) .ne. 0) then
                ttn2(9,jrow,n,mskhr(i,jrow)) = 
     &                   ttn2(9,jrow,n,mskhr(i,jrow)) + fxb
              endif
# endif
            enddo
          enddo
c
          if (sine(jrow) .eq. c0 .and. jrow .gt. 1) then
            factor = 8.0*omega*sine(jrow-1)
          else
            factor = 8.0*omega*sine(jrow)
          endif
          do i=is,ie
            totz = c0
            vbrz = c0
            tbrz = c0
            do k=1,km
	      mask = tmask(i,k,j)*tmask(i,k,j+1)
              vbrz = vbrz + adv_vnt(i,k,j)*dxt(i)*dzt(k)
              tbrz = tbrz +mask*(t(i,k,j,n,tau)+t(i,k,j+1,n,tau))*dzt(k)
              totz = totz + mask*dzt(k)
            enddo
            if (totz .ne. c0) then
              tbrz = tbrz/totz
              ttn(3,jrow,n) = ttn(3,jrow,n) + vbrz*tbrz*p5
              ttn(5,jrow,n) = ttn(5,jrow,n) - (smf(i,j,1)*dxu(i) + 
     &                        smf(i-1,j,1)*dxu(i-1))*(t(i,1,j,n,tau)
     &                        +t(i,1,j+1,n,tau)-tbrz)
     &         *csu(jrow)/factor
            endif
          enddo
          ttn(2,jrow,n) = ttn(6,jrow,n)-ttn(1,jrow,n)
          ttn(4,jrow,n) = ttn(6,jrow,n)-ttn(3,jrow,n)-ttn(5,jrow,n)
          ttn(8,jrow,n) = ttn(6,jrow,n)+ttn(7,jrow,n)
# if defined isopycmix && defined gent_mcwilliams
     &                  + ttn2(9,jrow,n,0)
# endif
          do ll=0,nhreg
            ttn2(8,jrow,n,ll) = ttn2(6,jrow,n,ll)+ttn2(7,jrow,n,ll)
# if defined isopycmix && defined gent_mcwilliams
     &                         +ttn2(9,jrow,n,ll) 
# endif
          enddo
        endif
      enddo
# ifdef timing
      call toc ('diagnostic', 'gyre components')
# endif
#endif
      return
      end
