      subroutine diag (joff, js, je, is, ie)
c
c=======================================================================
c     calculate diagnostics
c
c     input:
c
c      joff   = offset between row j in the MW and latitude jrow on disk
c               "joff" > 0 moves variables
c      js     = starting row for calculations
c      je     = ending row for calculations
c      is     = starting longitude index for calculations
c      ie     = ending longitude index for calculations
c=======================================================================
c
#include "param.h"
#include "coord.h"
#ifdef matrix_sections
# include "cprnts.h"
#endif
#include "diag.h"
#include "diaga.h"
#include "docnam.h"
#include "grdvar.h"
#include "iounit.h"
#include "isopyc.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"
#include "vmixc.h"
c
#ifdef meridional_overturning
      dimension vbarx(km)
#endif
c
c-----------------------------------------------------------------------
c     limit longitudes
c-----------------------------------------------------------------------
c
      istrt  = max(2,is)
      iend   = min(imt-1,ie)
c
      do j=js,je
        jrow = joff + j
c
#ifdef time_averages
c
c-----------------------------------------------------------------------
c       diagnostic: accumulate "tau" data for time means 
c       author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c-----------------------------------------------------------------------
c
        if (.not. euler2) then
	  if (istrt .ne. 2 .and. iend .ne. imt-1) then
	    write (stdout,*) '=>Error: istrt = ',istrt,' and iend ='
     &,     iend,' are not allowed when calling "avgvar"'
	    stop '=>diag'
          else
            call avgvar (j, jrow, adv_vbt(1,1,j), u(1,1,1,1,tau) 
     &,                  t(1,1,1,1,tau), stf, smf)
	  endif
        endif
#endif
c
#ifdef snapshots
c
c-----------------------------------------------------------------------
c       diagnostic: write instantaneous sample of MOM data
c       author:    r.c.pacanowski   e-mail  rcp@gfdl.gov
c                  a. rosati        e-mail  ar@gfdl.gov
c-----------------------------------------------------------------------
c
        if (snapts .and. .not. euler2) then
	  if (istrt .ne. 2 .and. iend .ne. imt-1) then
	    write (stdout,*) '=>Error: istrt = ',istrt,' and iend ='
     &,     iend,' are not allowed when calling "snapit"'
            stop '=>diag'
	  else
	    call snapit (j, jrow)
	  endif
        endif
#endif
c
#ifdef stability_tests
c
c-----------------------------------------------------------------------
c       diagnostic: compute stability diagnostics
c       author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c-----------------------------------------------------------------------
c
        if (stabts .and. eots) then
	  if (istrt .ne. 2 .and. iend .ne. imt-1) then
	    write (stdout,*) '=>Error: istrt = ',istrt,' and iend ='
     &,     iend,' are not allowed when calling "stab"'
            stop '=>diag'
	  else
	    call stab (j, jrow)
	  endif
        endif
#endif
c
#ifdef trajectories
c
c-----------------------------------------------------------------------
c       diagnostic: integrate particle trajectories
c       author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c-----------------------------------------------------------------------
c
        if (eots) then
	  if (istrt .ne. 2 .and. iend .ne. imt-1) then
	    write (stdout,*) '=>Error: istrt = ',istrt,' and iend ='
     &,     iend,' are not allowed when calling "ptraj"'
            stop '=>diag'
	  else
	    call ptraj (j, jrow)
          endif
        endif
#endif
c
#ifdef meridional_overturning
c
c-----------------------------------------------------------------------
c       construct meridional overturning of mass
c       author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c-----------------------------------------------------------------------
c
        if (jrow .lt. jmtm1 .and. vmsfts .and. eots) then
# ifdef timing
          call tic ('diagnostic', 'meridional overturning')
# endif
          do k=1,km
            vbarx(k) = c0
          enddo
c
          do k=1,km
            do i=istrt,iend
              vbarx(k) = vbarx(k) + u(i,k,j,2,tau)*csu(jrow)*dxu(i)
            enddo
            if (k .eq. 1) then
	      vsf(jrow,k) = vbarx(k)*dzt(k)
	    else
              vsf(jrow,k) = vsf(jrow,k-1) + vbarx(k)*dzt(k)
	    endif
          enddo
# ifdef timing
          call toc ('diagnostic', 'meridional overturning')
# endif
        endif
#endif
c
#ifdef show_zonal_mean_of_sbc
c
c-----------------------------------------------------------------------
c       construct zonal mean of surface b.c. and related items
c       author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c-----------------------------------------------------------------------
c
        if (zmbcts .and. eots) then
# ifdef timing
          call tic ('diagnostic', 'show zonal mean of s.b.c.')
# endif
          zmau(jrow) = c0
          zmat(jrow) = c0
          do i=istrt,iend
            zma1            = umask(i,1,j)*csu(jrow)*dxu(i)*dyu(jrow)
            zmau(jrow)      = zmau(jrow)    + zma1
            zmsmf(jrow,1)   = zmsmf(jrow,1) + zma1*smf(i,j,1)
            zmsmf(jrow,2)   = zmsmf(jrow,2) + zma1*smf(i,j,2)
            zmsm(jrow,1)    = zmsm(jrow,1)  + zma1*u(i,1,j,1,tau)
            zmsm(jrow,2)    = zmsm(jrow,2)  + zma1*u(i,1,j,2,tau)
            zma2            = tmask(i,1,j)*cst(jrow)*dxt(i)*dyt(jrow)
            zmat(jrow)      = zmat(jrow) + zma2
            do n=1,nt
              zmstf(jrow,n) = zmstf(jrow,n) + zma2*stf(i,j,n)
              zmst(jrow,n)  = zmst(jrow,n)  + zma2*t(i,1,j,n,tau)
            enddo
          enddo
# ifdef timing
          call toc ('diagnostic', 'show zonal mean of s.b.c.')
# endif
        endif
#endif
c
#ifdef matrix_sections
c
c-----------------------------------------------------------------------
c       print "tau" (not "tau+1") variables  at specified latitudes
c
c       author:    r.c.pacanowski   e-mail  rcp@gfdl.gov
c                  a. rosati        e-mail  ar@gfdl.gov
c       based on previous code by M. Cox
c-----------------------------------------------------------------------
c
        if (prxzts .and. eots) then
# ifdef timing
          call tic ('diagnostic', 'matrix sections')
# endif
          reltim = relyr
          do jlat=1,nlatpr
            jj = indp (prlat(jlat), yt, jmt)
            if (jj .eq. jrow .and. prlat(jlat) .le. yt(jmt))
     &        then
              is = indp (prslon(jlat), xt, imt)
              ie = indp (prelon(jlat), xt, imt)
              ks = indp (prsdpt(jlat), zt, km)
              ke = indp (predpt(jlat), zt, km)
              fx = 1.0e-2
              if (jlat .eq. 1) write(stdout,9000)
              do m=1,nt
                scl = c1
                if (m .eq. 2) scl=1.e-3
                if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
                  write (stdout,9100) trname(m), itt, jrow 
     &,           yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
                  call matrix (t(1,1,j,m,tau), imt, is, ie, ks, ke, scl)
                endif
                if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
c
	          call getunit (io, 'sections.dta','u s a ieee')
c
                  write (stdout,*) ' =>Zonal ',trname(m), ' slice: lat='
     &,           yt(jrow), ' written unformatted to file sections.dta'
     &,           ' on ts=', itt, stamp
                  write (stdout,'(///)')
                  iotext = ' read (ioprxz) imt, km, m, nt, reltim'
                  write (io) stamp, iotext, expnam
                  write (io) imt, km, m, nt, reltim
                  write(iotext,'(a10,i4,a4,i2)') ' for jrow=',jrow
     &,           '  m=',m
                  iotext(18:)=':read(ioprxz)((t(i,k,m),i=1,imt),k=1,km)'
                  write (io) stamp, iotext, expnam
                  call wrufio (io, t(1,1,j,m,tau), imt*km)
c
                  call relunit (io)
c
                endif
              enddo
              scl = 1.e-3
              if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
	        i1 = max(is,2)
		i2 = min(ie,imtm1)
                write (stdout,9100)  'adv_vbu ', itt, jrow
     &,         yt(jrow), xt(i1), xt(i2), fx*zw(ks), fx*zw(ke), scl
                call matrix (adv_vbu(1,1,j), imt, i1, i2, ks, ke, scl)
                write (stdout,9100)  'adv_vbt ', itt, jrow
     &,         yt(jrow), xt(i1), xt(i2), fx*zw(ks), fx*zw(ke), scl
                call matrix (adv_vbt(1,1,j), imt, i1, i2, ks, ke, scl)
              endif
              if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
c
	        call getunit (io, 'sections.dta','u s a ieee')
c
                write (stdout,*) ' => Zonal adv_vbu slice: lat='
     &,         yt(jrow), ' written unformatted to file sections.dta'
     &,         ' on ts=', itt, stamp
                write (stdout,'(///)')
                iotext = ' read (ioprxz) imt, km, reltim'
                write (io) stamp, iotext, expnam
                write (io) imt, km, reltim
                write(iotext,'(a10,i4)') ' for jrow=',jrow
                iotext(12:)=
     &          ': read(ioprxz)((adv_vbu(i,k),i=1,imt),k=0,km)'
                write (io) stamp, iotext, expnam
                call wrufio (io, adv_vbu(1,0,j), imt*(km+1))
c
                write (stdout,*) ' => Zonal adv_vbt slice: lat='
     &,         yt(jrow), ' written unformatted to file sections.dta'
     &,         ' on ts=', itt, stamp
                write (stdout,'(///)')
                iotext = ' read (ioprxz) imt, km, reltim'
                write (io) stamp, iotext, expnam
                write (io) imt, km, reltim
                write(iotext,'(a10,i4)') ' for jrow=',jrow
                iotext(12:)=
     &          ': read(ioprxz)((adv_vbt(i,k),i=1,imt),k=1,km)'
                write (io) stamp, iotext, expnam
                call wrufio (io, adv_vbt(1,0,j), imt*(km+1))
c
                call relunit (io)
c
              endif
c
              scl = c1
              if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
                write (stdout,9100) 'u velocity', itt
     &,        jrow, yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
                call matrix (u(1,1,j,1,tau), imt, is, ie, ks, ke, scl)
              endif
              if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
c
	        call getunit (io, 'sections.dta','u s a ieee')
c
                write (stdout,*) ' => Zonal u velocity slice: lat='
     &,         yt(jrow), ' written unformatted to file sections.dta'
     &,         ' on ts=', itt, stamp
                write (stdout,'(///)')
                iotext = ' read (ioprxz) imt, km, reltim'
                write (io) stamp, iotext, expnam
                write (io) imt, km, reltim
                write(iotext,'(a10,i4)') ' for jrow=',jrow
                iotext(12:) = ': read (ioprxz)((u(i,k),i=1,imt),k=1,km)'
                write (io) stamp, iotext, expnam
                call wrufio (io, u(1,1,j,1,tau), imt*km)
c
                call relunit (io)
c
              endif
c
              scl = c1
              if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
                write (stdout,9100) 'v velocity', itt, jrow
     &,         yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
                call matrix (u(1,1,j,2,tau), imt, is, ie, ks, ke, scl)
              endif
              if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
c
	        call getunit (io, 'sections.dta','u s a ieee')
c
                write (stdout,*) ' => Meridional v velocity slice: lat='
     &,         yt(jrow),' written unformatted to file sections.dta'
     &,          ' on ts=', itt, stamp
                write (stdout,'(///)')
                iotext = ' read (ioprxz) imt, km, reltim'
                write (io) stamp, iotext, expnam
                write (io) imt, km, reltim
                write(iotext,'(a10,i4)') ' for jrow=',jrow
                iotext(12:) = ': read (ioprxz)((v(i,k),i=1,imt),k=1,km)'
                write (io) stamp, iotext, expnam
                call wrufio (io, u(1,1,j,2,tau), imt*km)
c
                call relunit (io)
c
              endif
            endif
          enddo
# ifdef timing
          call toc ('diagnostic', 'matrix sections')
# endif
        endif
9000    format(/' Zonal section printouts at specified latitudes:'/)
9100    format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     &,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     &,', scaling=',1pg10.3)
#endif
      enddo
c
#if defined save_mixing_coeff
c
c-----------------------------------------------------------------------
c     diagnostic: save estimated mixing coefficients on east, north, and
c                 bottom face of T and U cells
c
c     author:   R. C. Pacanowski       e-mail rcp@gfdl.gov
c-----------------------------------------------------------------------
c
      if (cmixts .and. eots) then
# ifdef timing
        call tic ('diagnostic', 'mixing coeff')
# endif
c
        if (joff .eq. 0) then
          write (stdout,*) ' =>Writing mixing coefficients at ts=',itt
     &  , ' ',stamp
	  call getunit (iocm, 'cmix.dta','u s a ieee')
c
          period = 0.0
	  reltim = relyr
          iotext = 'read(iocm) reltim, period, imt, jmt, km'
          write (iocm) stamp, iotext, expnam
          write (iocm) reltim, period, imt, jmt, km
c
          iotext = 'read(iocm) (xt(i),i=1,imt)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, xt, imt)
c
          iotext = 'read(iocm) (yt(j),j=1,jmt)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, yt, jmt)
c
          iotext = 'read(iocm) (zt(k),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, zt, km)
c
          iotext = 'read(iocm) (xu(i),i=1,imt)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, xu, imt)
c
          iotext = 'read(iocm) (yu(j),j=1,jmt)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, yu, jmt)
c
          iotext = 'read(iocm) (zw(k),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, zw, km)
c
          call relunit (iocm)
        endif
c
        call getunit (iocm, 'cmix.dta','u s a ieee')   
c
        do j=js,je
          jrow = j+joff
c
          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (diff_ceu(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, ce(1,1,j,1), imt*km)
c
          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (diff_cnu(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, cn(1,1,j,1), imt*km)
c
          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (diff_cbu(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, cb(1,1,j,1), imt*km)
c
          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (diff_cet(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, ce(1,1,j,2), imt*km)
c
          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (diff_cnt(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, cn(1,1,j,2), imt*km)
c
          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (diff_cbt(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, cb(1,1,j,2), imt*km)
c
#ifdef isopycmix
c
          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (K13(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, K1(1,1,j,3), imt*km)
c
          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (K23(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, K2(1,1,j,3), imt*km)
c
          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (K31(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, K3(1,1,j,1), imt*km)
c
          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (K32(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, K3(1,1,j,2), imt*km)
c
          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (K33(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, K3(1,1,j,3), imt*km)
#endif
        enddo
        call relunit(iocm)
# ifdef timing
        call toc ('diagnostic', 'mixing coeff')
# endif
      endif
#endif
      return
      end
