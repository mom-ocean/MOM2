c====================== include file "dens.h" ==========================
c
c-----------------------------------------------------------------------
c     statement function
c-----------------------------------------------------------------------
c
#ifdef linearized_density
c
c     approximating rho = 1.035*(1-alpha*tq). The one is removed by
c     gradients and the 1.035 is absorbed into alpha.
c
      dens(tq,sq,k) = -2.e-4*tq
#else
      dens (tq, sq, k) = (c(k,1) + (c(k,4) + c(k,7)*sq)*sq +
     &                   (c(k,3) + c(k,8)*sq + c(k,6)*tq)*tq)*tq +
     &                   (c(k,2) + (c(k,5) + c(k,9)*sq)*sq)*sq
#endif
