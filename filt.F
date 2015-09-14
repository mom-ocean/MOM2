      subroutine filt (joff, js, je)
#if defined fourfil || defined firfil
c
c=======================================================================
c     subroutine filt sets up input needed for fourier filtering
c     (when the "fourfil" option is defined) -or- symmetric finite
c     impulse response filtering (when the "firfil" option is defined)
c     of tracers at the specifiied high latitude row "jrow".
c=======================================================================
c
# include "param.h"
# include "grdvar.h"
# include "index.h"
# include "levind.h"
# include "mw.h"
# ifdef fourfil
      dimension tempik(imt,km)
# endif
# ifdef firfil
      dimension temp(imt,km,jsmw:jemw)
      dimension jtof(jmw)
# endif
c
# ifdef timing
      call tic ('tracer', 'filtering')
# endif
      do n=1,nt
        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo
c
# ifdef fourfil
c
c---------------------------------------------------------------------
c     fourier filter tracers at high latitudes
c---------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
      if((jrow.gt.jft1.and.jrow.lt.jft2) .or. jrow.lt.jfrst) goto 101
      jj = jrow-jfrst+1
      if (jrow .ge. jft2) jj = jj-jskpt+1
c
c    if previous strips were of same length, do not recompute
c    fourier coeffs
c
      isave = 0
      ieave = 0
      do l=1,lsegf
        do k=1,km
          if (istf(jj,l,k) .ne. 0) then
            is    = istf(jj,l,k)
            ie    = ietf(jj,l,k)
            iredo = 0
            if (is.ne.isave .or. ie.ne.ieave) then
              iredo = -1
              isave = is
              ieave = ie
              im = ie-is+1
#  ifdef cyclic
              if (im.ne.imtm2 .or. kmt(1,jrow).lt.k) then
                m = 1
                n = nint(im*cst(jrow)*cstr(jft0))
              else
                m = 3
                n = nint(im*cst(jrow)*cstr(jft0)*0.5)
              endif
#  else
              m = 1
              n = nint(im*cst(jrow)*cstr(jft0))
#  endif
            endif
            do mm=1,nt
              idx  = iredo+mm
              ism1 = is-1
              iea  = ie
              if (ie .ge. imt) iea = imtm1
              do i=is,iea
                tempik(i-ism1,k) = t(i,k,j,mm,taup1)
              enddo
              if (ie .ge. imt) then
                ieb = ie-imtm2
                ii  = imtm1-is
                do i=2,ieb
                  tempik(i+ii,k) = t(i,k,j,mm,taup1)
                enddo
              endif
c
              call filtr (tempik(1,k), im, m, n, idx)
c
              do i=is,iea
                t(i,k,j,mm,taup1) = tempik(i-ism1,k)
              enddo
              if (ie .ge. imt) then
                do i=2,ieb
                  t(i,k,j,mm,taup1) = tempik(i+ii,k)
                enddo
              endif
            enddo
          endif
        enddo
      enddo
101   continue
      enddo
# endif
c
# ifdef firfil
c
c     build starting and ending rows to filter
c
      jsf = 0
      jef = 0
      do j=js,je
        jrow = j + joff
        if ((jrow.le.jft1 .or. jrow .ge. jft2) .and. jrow .ge. jfrst)
     &    then
          if (jsf .eq. 0) jsf = j
	  jef = j 
	endif
      enddo
      if (jsf .eq. 0) goto 101
c
      do j=jsf,jef
        jrow = j + joff
        jj = jrow - jfrst + 1
        if (jrow .ge. jft2) jj = jj - jskpt + 1
        if ((jrow.le.jft1 .or. jrow .ge. jft2) .and. jrow .ge. jfrst)
     &    then
          jtof(j) = numflt(jj)
	else
          jtof(j) = 0
	endif
      enddo
c
c
c-----------------------------------------------------------------------
c     filter tracers at high latitudes with symmetric finite impulse
c     response filter
c-----------------------------------------------------------------------
c
        do mm=1,nt
          call filtrb (t(1,1,jsmw,mm,taup1), tmask(1,1,jsmw)
     &,                temp(1,1,jsmw), km, jtof, jsf, jef)
        enddo
101   continue
# endif
# ifdef timing
      call toc ('tracer', 'filtering')
# endif
#endif
      return
      end

#ifdef firfil
      subroutine filtrb (t, f, s, kl, jtof, jsf, jef)
c
c=======================================================================
c     simple finite impulse response filter with [.25, .5, .25] weights
c     using symmetric boundary conditions on each latitude row. this
c     filter does an entire row at a time.
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
# include "param.h"
c
      dimension t(imt,kl,jsmw:jemw), f(imt,kl,jsmw:jemw)
      dimension s(imt,kl,jsmw:jemw)
      dimension jtof(jmw)
c
# ifdef timing
      call tic ('filtering', 'filfir (finite impulse)')
# endif
c
      do j=jsf,jef
        do k=1,kl
# ifdef cyclic
          t(1,k,j)   = t(imtm1,k,j)
          t(imt,k,j) = t(2,k,j)
# else
          s(1,k,j)   = c0
          s(imt,k,j) = c0
# endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     apply the filter "num" times using a symmetric (no flux)
c     boundary condition
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
# ifdef cyclic
          do k=1,kl
            s(1,k,j)   = s(imtm1,k,j)
            s(imt,k,j) = s(2,k,j)
          enddo
# endif
          do k=1,kl
            do i=2,imtm1
              t(i,k,j) = f(i,k,j)*(p25*(s(i-1,k,j) + s(i+1,k,j)) +
     &                   s(i,k,j)*(c1 - p25*(f(i-1,k,j) + f(i+1,k,j))))
            enddo
          enddo
# ifdef cyclic
          do k=1,kl
            t(1,k,j)   = t(imtm1,k,j)
            t(imt,k,j) = t(2,k,j)
          enddo
# endif
        enddo
      enddo
c
# ifdef timing
      call toc ('filtering', 'filfir (finite impulse)')
# endif
      return
      end
#endif
