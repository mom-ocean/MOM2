      subroutine invtri (z, topbc, botbc, dcb, tdt, kmz, mask, is, ie
     &,                  joff, js, je)
#if defined implicitvmix || defined isopycmix
c
c=======================================================================
c     solve the vertical diffusion equation implicitly using the
c     method of inverting a tridiagonal matrix as described in
c     Richtmyer and Morton, 1967, Difference Methods for Initial Value 
c     Problems.
c     this routine assums that the variables are defined at grid points
c     and the top and bottom b.c. are flux conditions.
c
c     inputs:
c     z         = right hand side terms
c     topbc     = top boundary condition
c     botbc     = bottom boundary condition
c     dcb       = vertical mixing coeff
c     tdt       = 2 * timestep
c     kmz       = level indicator
c     mask      = land - sea mask
c     is        = index of starting longitude
c     ie        = index of ending longitude
c     js        = starting latitude row in MW
c     je        = ending latitude row in MW
c     joff      = offset between jrow on disk and j in the MW
c
c     outputs:
c     z         = returned solution
c
c=======================================================================
c
#include "param.h"
#include "grdvar.h"
#include "vmixc.h"
      dimension z(imt,km,jmw)
      dimension topbc(imt,jsmw:jemw), botbc(imt,jsmw:jemw)
      dimension dcb(imt,km,jsmw:jemw)
      dimension kmz(imt,jmt), tdt(km)
      real mask(imt,km,1:jmw)
c
      dimension a(imt,km,jsmw:jemw), b(imt,km,jsmw:jemw)
      dimension c(imt,km,jsmw:jemw), d(imt,km,jsmw:jemw)
      dimension e(imt,0:km,jsmw:jemw), f(imt,0:km,jsmw:jemw)
      dimension g(imt,jsmw:jemw)
c
      do j=js,je
        do k=2,km
          do i=is,ie
            a(i,k,j)   = dcb(i,k-1,j)*dztur(k)*tdt(k)*aidif
            c(i,k,j)   = dcb(i,k,j)*dztlr(k)*tdt(k)*aidif
            b(i,k,j)   = c1 + a(i,k,j) + c(i,k,j)
            d(i,k,j)   = z(i,k,j)
            e(i,k-1,j) = c0
            f(i,k-1,j) = c0
          enddo
        enddo
      enddo
c
c     b. c. at top
c
      k = 1
      do j=js,je
        do i=is,ie
          a(i,k,j)   = dztr(k)*tdt(k)*aidif
          c(i,k,j)   = dcb(i,k,j)*dztlr(k)*tdt(k)*aidif
          b(i,k,j)   = c1 + c(i,k,j)
          d(i,k,j)   = z(i,k,j)
          e(i,k-1,j) = c0
          f(i,k-1,j) = c0
        enddo
      enddo
c
c     b. c. at bottom
c
      do j=js,je
        jrow = j + joff
        do i=is,ie
          kz = kmz(i,jrow)
          if (kz .ne. 0) then
            b(i,kz,j) = c1 + a(i,kz,j)
            c(i,kz,j) = dztr(kz)*tdt(kz)*aidif
            e(i,kz,j) = c0
            f(i,kz,j) = -botbc(i,j)
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     now invert
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        do k=km,1,-1
          do i=is,ie
            if (k .le. kmz(i,jrow)) then
              g(i,j)     = c1/(b(i,k,j)-c(i,k,j)*e(i,k,j))
              e(i,k-1,j) = a(i,k,j)*g(i,j)
              f(i,k-1,j) = (d(i,k,j)+c(i,k,j)*f(i,k,j))*g(i,j)
            endif
          enddo
        enddo
      enddo
c
c     b.c. at surface
c
      do j=js,je
        do i=is,ie
          z(i,1,j) = (e(i,0,j)*topbc(i,j) + f(i,0,j))*mask(i,1,j)
        enddo
      enddo
c
      do j=js,je
        do k=2,km
          do i=is,ie
            z(i,k,j) = (e(i,k-1,j)*z(i,k-1,j) + f(i,k-1,j))*mask(i,k,j)
          enddo
        enddo
      enddo
#endif
      return
      end
