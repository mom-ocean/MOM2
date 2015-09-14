      subroutine filfir (t, f, s, kl, kind, jtof, jsf, jef)
#ifdef firfil
c
c=======================================================================
c     simple finite impulse response filter with [.25, .5, .25] weights
c     modified for assymetric and symmetric boundary conditions
c
c     input:
c             t    = array of quantity to be filtered along
c                    the first dimension.
c                    note: t(i,k) must be zero where f(i,k) = zero
c                    for this filter to work.
c             f    = mask of zeroes & ones to indicate land
c                    and ocean. zero indicates a land point
c             s    = scratch array
c             kl   = number of vertical levels to filter
c             kind = (0,1) = (symmetric, asymmetric) boundary condition
c                    symmetric is appropriate for tracers & vorticity
c                    asymmetric is appropriate for velocities
c             jtof = number of filter passes per row
c             jsf  = starting row
c             jef  = ending row
c
c     output:
c             t    = (imt,km) array of filtered quantities
c
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
c
      dimension t(imt,kl,jsmw:jemw), f(imt,kl,jsmw:jemw)
      dimension s(imt,kl,jsmw:jemw)
      dimension jtof(jmw)
c
      call tic ('filtering', 'filfir (finite impulse)')
      do j=jsf,jef
        call setbcx (t(1,1,j), imt, kl)
      enddo
      if (kind .eq. 0) then
c
c-----------------------------------------------------------------------
c       apply the filter "num" times using a symmetric (no flux)
c       boundary condition
c-----------------------------------------------------------------------
c
        do j=jsf,jef
	  num = jtof(j)
	  do n=1,num
            do k=1,kl
              do i=2,imtm1
                s(i,k,j) = f(i,k,j)*(p25*(t(i-1,k,j) + t(i+1,k,j)) +
     &                   t(i,k,j)*(c1 - p25*(f(i-1,k,j) + f(i+1,k,j))))
              enddo
            enddo
            call setbcx (s(1,1,j), imt, kl)
            do k=1,kl
              do i=2,imtm1
                t(i,k,j) = f(i,k,j)*(p25*(s(i-1,k,j) + s(i+1,k,j)) +
     &                   s(i,k,j)*(c1 - p25*(f(i-1,k,j) + f(i+1,k,j))))
              enddo
            enddo
            call setbcx (t(1,1,j), imt, kl)
          enddo
        enddo
      elseif (kind .eq. 1) then
c
c----------------------------------------------------------------------
c       apply the filter "num" times using an asymmetric (flux)
c       boundary condition
c----------------------------------------------------------------------
c
        do j=jsf,jef
	  num = jtof(j)
	  do n=1,num
            do k=1,kl
              do i=2,imtm1
                s(i,k,j) = f(i,k,j)*(p25*t(i-1,k,j) + p5*t(i,k,j) + 
     &                               p25*t(i+1,k,j))
              enddo
            enddo
            call setbcx (s(1,1,j), imt, kl)
            do k=1,kl
              do i=2,imtm1
                t(i,k,j) = f(i,k,j)*(p25*s(i-1,k,j) + p5*s(i,k,j)
     &                              +p25*s(i+1,k,j))
              enddo
            enddo
            call setbcx (t(1,1,j), imt, kl)
          enddo
        enddo
      else
        write (stdout,'(/a,i10,a)') ' error=> kind =', kind,' in filfir'
        stop '=>filfir'
      endif
c
      call toc ('filtering', 'filfir (finite impulse)')
#endif
      return
      end
