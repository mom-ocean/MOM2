c====================== include file "fdift.h" =========================
c
c     finite difference numerics for tracers
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
c-----------------------------------------------------------------------
c     advective terms
c-----------------------------------------------------------------------
c
#ifdef linearized_advection
      ADV_Tx(i,k,j) = 0.0
      ADV_Ty(i,k,j) = 0.0
      ADV_Tz(i,k,j) = (adv_fb(i,k-1,j) - adv_fb(i,k,j))*dzt2r(k)
#else
      ADV_Tx(i,k,j) = (adv_fe(i,k,j) - adv_fe(i-1,k,j))*cstdxt2r(i,j)
      ADV_Ty(i,k,j) = (adv_vnt(i,k,j)*(t(i,k,j,n,tau)+t(i,k,j+1,n,tau))
     &              - adv_vnt(i,k,j-1)*(t(i,k,j-1,n,tau)+t(i,k,j,n,tau))
     &                )*cstdyt2r(jrow)
      ADV_Tz(i,k,j) = (adv_fb(i,k-1,j) - adv_fb(i,k,j))*dzt2r(k)
#endif
#if defined gent_mcwilliams
c
c     gent_mcwilliams isopycnal advective terms simulating the effect
c     of eddies on the isopycnals
c
      ADV_Txiso(i,k,j) = cstdxt2r(i,j)*
     &                   (adv_vetiso(i,k,j)*
     &                    (t(i+1,k,j,n,taum1) + t(i,k,j,n,taum1))
     &                   -adv_vetiso(i-1,k,j)*
     &                    (t(i,k,j,n,taum1) + t(i-1,k,j,n,taum1)))
      ADV_Tyiso(i,k,j) = cstdyt2r(jrow)*
     &                   (adv_vntiso(i,k,j)*
     &                   (t(i,k,j+1,n,taum1) + t(i,k,j,n,taum1))
     &                  -adv_vntiso(i,k,j-1)*
     &                   (t(i,k,j,n,taum1) + t(i,k,j-1,n,taum1)))   
      ADV_Tziso(i,k,j) = dzt2r(k)*(adv_fbiso(i,k-1,j)-adv_fbiso(i,k,j))
#endif
c
c-----------------------------------------------------------------------
c     diffusive terms
c-----------------------------------------------------------------------
c
c     zonal component
c
      DIFF_Tx(i,k,j) = (diff_fe(i,  k,j)*tmask(i+1,k,j)
     &                - diff_fe(i-1,k,j)*tmask(i-1,k,j))*cstdxtr(i,j)
c
c     meridional component
c
#if defined consthmix && !defined biharmonic && !defined isopycmix
# ifdef bryan_lewis_horizontal
      DIFF_Ty(i,k,j) = ahc_north(jrow,k)*tmask(i,k,j+1)*
     &                 (t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
     &               - ahc_south(jrow,k)*tmask(i,k,j-1)*
     &                 (t(i,k,j,n,taum1) - t(i,k,j-1,n,taum1))
# else
      DIFF_Ty(i,k,j) = ahc_north(jrow)*tmask(i,k,j+1)*
     &                 (t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
     &               - ahc_south(jrow)*tmask(i,k,j-1)*
     &                 (t(i,k,j,n,taum1) - t(i,k,j-1,n,taum1))
# endif
#else
      DIFF_Ty(i,k,j) = (diff_fn(i,k,j  )*tmask(i,k,j+1)
     &                - diff_fn(i,k,j-1)*tmask(i,k,j-1))*cstdytr(jrow)
#endif
c
c     vertical component
c
      DIFF_Tz(i,k,j) = (diff_fb(i,k-1,j) - diff_fb(i,k,j))*dztr(k)
#if defined implicitvmix || defined isopycmix
     &                 *(c1-aidif)
#endif
#if defined isopycmix
     &               + (diff_fbiso(i,k-1,j) - diff_fbiso(i,k,j))*dztr(k)
c
c   The K31 and K32 terms are in "diff_fbiso and treated explicitly.
c   The K33 is treated semi-implicitly in diff_fb
c
#endif
