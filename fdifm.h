c====================== include file "fdifm.h" =========================
c
c     finite difference numerics for momentum
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
c-----------------------------------------------------------------------
c     advective terms
c-----------------------------------------------------------------------
c
#ifdef linearized_advection
      ADV_Ux(i,k,j) = 0.0
      ADV_Uy(i,k,j) = 0.0
      ADV_Uz(i,k,j) = 0.0
      ADV_metric(i,k,j,n)=0.0
#else
      ADV_Ux(i,k,j) = (adv_fe(i,k,j) - adv_fe(i-1,k,j))*csudxu2r(i,j)
      ADV_Uy(i,k,j) = (adv_vnu(i,k,j)*(u(i,k,j,n,tau)+u(i,k,j+1,n,tau))
     &              - adv_vnu(i,k,j-1)*(u(i,k,j-1,n,tau)+u(i,k,j,n,tau))
     &                )*csudyu2r(jrow)
      ADV_Uz(i,k,j) = (adv_fb(i,k-1,j) - adv_fb(i,k,j))*dzt2r(k)
      ADV_metric(i,k,j,n)=advmet(jrow,n)*u(i,k,j,1,tau)*u(i,k,j,3-n,tau)
#endif
c
c-----------------------------------------------------------------------
c     diffusive terms
c-----------------------------------------------------------------------
c
      DIFF_Ux(i,k,j) = (diff_fe(i,k,j) - diff_fe(i-1,k,j))
     &                 *rho0csudxur(i,j)
      DIFF_Uz(i,k,j) = (diff_fb(i,k-1,j) - diff_fb(i,k,j))*rho0dztr(k)
#ifdef implicitvmix
     &                *(c1-aidif)
#endif
#if defined consthmix && !defined biharmonic
# ifdef free_slip
      DIFF_Uy(i,k,j) = amc_north(jrow)*
     &               (u(i,k,j+1,n,taum1) - u(i,k,j,n,taum1))
     &              *(umask(i,k,j)*(c1-umask(i-1,k,j)) + umask(i-1,k,j))
     &             - amc_south(jrow)*
     &               (u(i,k,j,n,taum1) - u(i,k,j-1,n,taum1))
     &      *(umask(i,k,j-1)*(c1-umask(i-1,k,j-1)) + umask(i-1,k,j-1))
# else
      DIFF_Uy(i,k,j) = amc_north(jrow)*
     &                (u(i,k,j+1,n,taum1) - u(i,k,j,n,taum1))
     &              - amc_south(jrow)*
     &                (u(i,k,j,n,taum1) - u(i,k,j-1,n,taum1))
# endif
#else
      DIFF_Uy(i,k,j) = (diff_fn(i,k,j) - diff_fn(i,k,j-1))
     &                 *rho0csudyur(jrow)
#endif
c
c-----------------------------------------------------------------------
c     metric term
c-----------------------------------------------------------------------
c
#if defined consthmix
# if defined  biharmonic
      DIFF_metric(i,k,j,n) = am3(jrow)*del2(i,k,j,n) + am4(jrow,n)*
     &              (del2(i+1,k,j,3-n) - del2(i-1,k,j,3-n))*dxmetr(i)
# else
#  ifdef free_slip
      DIFF_metric(i,k,j,n) = am3(jrow)*u(i,k,j,n,taum1) 
     &                     + am4(jrow,n)*dxmetr(i)
     &             *(
     &               (n-1)*(u(i+1,k,j,3-n,taum1) - u(i-1,k,j,3-n,taum1))
     &             + (2-n)*(u(i+1,k,j,3-n,taum1)*umask(i+1,k,j) +
     &                      u(i  ,k,j,3-n,taum1)*(c1-umask(i+1,k,j)) -
     &                     (u(i-1,k,j,3-n,taum1)*umask(i-1,k,j)
     &                    + u(i  ,k,j,3-n,taum1)*(c1-umask(i-1,k,j))))
     &              )
#  else
      DIFF_metric(i,k,j,n) = am3(jrow)*u(i,k,j,n,taum1)
     &                     + am4(jrow,n)*dxmetr(i)
     &                    *(u(i+1,k,j,3-n,taum1) - u(i-1,k,j,3-n,taum1))
#  endif
# endif
#endif
#if defined smagnlmix
      DIFF_metric(i,k,j,n) = smag_metric(i,k,j)
#endif
c
c-----------------------------------------------------------------------
c     coriolis term
c-----------------------------------------------------------------------
c
      CORIOLIS(i,k,j,n) = cori(jrow,n)*( 
     &             gcor*u(i,k,j,3-n,tau) + (c1-gcor)*u(i,k,j,3-n,taum1))
c
