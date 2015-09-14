#ifdef makeic
# define driver_only
# include "iomngr.F"
# include "grids.F"
# include "size_check.F"
# include "util.F"

      program ic
c
c-----------------------------------------------------------------------
c     prepare initial temperature & salt for MOM based on
c     monthly Levitus climatology from MOM database
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
# include "size.h"
# include "stdunits.h"
c
      character*12 tempin, tempout, saltin, saltout
      dimension tempin(13), tempout(13), saltin(13), saltout(13)
      dimension iotin(13), iotout(13), iosin(13), iosout(13)
c
      data tempout/
     &               'jan.temp.mom', 'feb.temp.mom', 'mar.temp.mom'
     &,              'apr.temp.mom', 'may.temp.mom', 'jun.temp.mom'
     &,              'jul.temp.mom', 'aug.temp.mom', 'sep.temp.mom'
     &,              'oct.temp.mom', 'nov.temp.mom', 'dec.temp.mom'
     &,              'ann.temp.mom' /
c
      data saltout/
     &               'jan.salt.mom', 'feb.salt.mom', 'mar.salt.mom'
     &,              'apr.salt.mom', 'may.salt.mom', 'jun.salt.mom'
     &,              'jul.salt.mom', 'aug.salt.mom', 'sep.salt.mom'
     &,              'oct.salt.mom', 'nov.salt.mom', 'dec.salt.mom'
     &,              'ann.salt.mom' /
c
c
      data tempin /
     &               'jan.temp', 'feb.temp', 'mar.temp'
     &,              'apr.temp', 'may.temp', 'jun.temp'
     &,              'jul.temp', 'aug.temp', 'sep.temp'
     &,              'oct.temp', 'nov.temp', 'dec.temp'
     &,              'ann.temp' /
c
      data saltin /
     &               'jan.salt', 'feb.salt', 'mar.salt'
     &,              'apr.salt', 'may.salt', 'jun.salt'
     &,              'jul.salt', 'aug.salt', 'sep.salt'
     &,              'oct.salt', 'nov.salt', 'dec.salt'
     &,              'ann.salt' /
c
c     intrp  = 1 when MOM grid resolution is less than Levitus
c              2 when MOM grid resolution is greater than Levitus
c     intrpv = 1 when MOM grid resolution is less than Levitus in the
c              vertical
c              2 when MOM grid resolution is greater than Levitus in the
c              vertical
c
      data intrp /1/, intrpv /1/
c
      namelist /input/   intrp, intrpv
c
# ifdef cyclic
      write (stdout,*) '  ==> cyclic option set'
# else
      write (stdout,*) '  ==> cyclic option NOT set'
# endif
      read  (stdin, input)
      write (stdout,input)
c
      if (intrp .eq. 1) then
        print *,' Assuming MOM horizontal resolution is < than Levitus'
      elseif (intrp .eq. 2) then
        print *,' Assuming MOM horizontal resolution is > than Levitus'
      else
        print *,' intrp must be 1 or 2'
	stop
      endif
c
      if (intrpv .eq. 1) then
        print *,' Assuming MOM vertical resolution is < than Levitus'
      elseif (intrpv .eq. 2) then
        print *,' Assuming MOM vertical resolution is > than Levitus'
      else
        print *,' intrpv must be 1 or 2'
	stop
      endif
c
c-----------------------------------------------------------------------
c     initialize the i/o units
c-----------------------------------------------------------------------
c
      call ioinit
      do n=1,13
        call getunit (iotin(n), tempin(n), 'u s r ieee')
        call getunit (iotout(n), tempout(n), 'u s r ieee')
        call getunit (iosin(n), saltin(n), 'u s r ieee')
        call getunit (iosout(n), saltout(n), 'u s r ieee')
      enddo
      call getunit (iosst, 'sst.mom', 'u s r ieee')
      call getunit (iosss, 'sss.mom', 'u s r ieee')
c
c
c-----------------------------------------------------------------------
c     interpolate levitus data (T & S) to MOM "t" grid
c     temperature varies monthly down to level 19
c     salinity varies monthly down to level 24
c-----------------------------------------------------------------------
c
      months = 13
      lev = 19
      it = 0
      write (stdout,'(/a/)') 'Converting Levitus Temperature'
      do n=1,months
        call doit (lev, iotin, iotout, iosst, n, intrp, intrpv, it)
      enddo
c
      lev = 24
      write (stdout,'(/a/)') 'Converting Levitus Salinity'
      do n=1,months
        call doit (lev, iosin, iosout, iosss, n, intrp, intrpv, it)
      enddo
c
      stop
      end


      subroutine doit (lev, ioin, ioout, iosurf, month, intrp, intrpv
     &,                it)
# include "size.h"
# include "stdunits.h"
      parameter (imd=362, jmd=180, kd=33, lenw=imd*20)
      character*32 stamp, stamp2
      character*80 iotext
# include "coord.h"
c
      common /bigone/ lons, lats, depth, cx(imd), cy(jmd)
      common /bigone/ dd(kd), data(imd,jmd), tdp(imd), cxt(imd)
c
      dimension block(imt,jmt,kd), tm(imt,km), tij(imt,jmt), tdv(imt,kd)
      common /icwrk/ work1(lenw), work2(lenw)

      dimension ioin(13), ioout(13)
c
      it = it + 1
      if (it .eq. 1) then
c
c-----------------------------------------------------------------------
c       construct the grid.
c-----------------------------------------------------------------------
c
        write (stdout,'(//1x,a,//)') 'Construct the MOM grid'
c
        maxlen = max(imt,jmt,km)
        call gcoord (maxlen, imt2, jmt2, km2, dxtdeg, dytdeg, dxudeg
     &,              dyudeg, dzt, dzw, xt, xu, yt, yu, zt, zw)
c
        call size_check (imt2, jmt2, km2, 'ic', 'stop')
c
c-----------------------------------------------------------------------
c       bring in the grid definition for levitus climatology
c-----------------------------------------------------------------------
c
        read (ioin(13))
        read (ioin(13)) stamp, dpm, lons, lats, kk, kl, n, cx
     &,                  cy, dd, data
# ifdef cyclic
        call tranlon (data, imd, imd-2, jmd, tdp, cx, xt(2), imt-2, cxt)
# else
        call tranlon (data, imd, imd-2, jmd, tdp, cx, xt, imt, cxt)
# endif
c
        write(stdout,'(///,1x, a,/)') 'Levitus grid definition:'
        write (stdout,'(//,1x,a,/)') '            latitudes:'
        write (stdout,'(1x,10f10.5)') cy
        write (stdout,'(//,1x,a,/)') '  original longitudes:'
        write (stdout,'(1x,10f10.5)') cx
        write (stdout,'(//,1x,a,/)') 'translated longitudes:'
        write (stdout,'(1x,10f10.5)') cxt
        write (stdout,'(//,1x,a,/)') '               depths:'
        write (stdout,'(1x,10g10.2)') dd
        write (stdout,'(///)')
      endif
c
      rewind ioin(13)
      do k=1,kd
c
c       read in Levitus data at standard level k
c
        read (ioin(13))
        read (ioin(13)) stamp2, dpm2, lons, lats, kk, kl, n, cx
     &,                   cy, dd, data
        if (month .ne. 13) then
	  if (k .le. lev) then
            read (ioin(month))
            read (ioin(month)) stamp, dpm, lons, lats, kk, kl, n, cx
     &,                        cy, dd, data
          endif
	else
	  stamp = stamp2
	  dpm   = dpm2
	endif
c
c       translate levitus data on standard level k and interpolate 
c       to MOM "t" grid
c
        if (month .eq. 1 .and. k .eq. 1) then
	  nnn = 1 
        else
	  nnn = 2
	endif
# ifdef cyclic
        call tranlon (data, imd, imd-2, jmd, tdp, cx, xt(2), imt-2
     &,               cxt)
        if (intrp .eq. 1) then
          call ftc (data, imd, jmd, cxt, cy, block(1,1,k), imt, jmt, 2
     &,         imt-1, 1, jmt, xt, yt, nnn, work1, lenw)
        else if (intrp .eq. 2) then
          call ctf (data, imd, jmd, cxt, cy, block(1,1,k), imt, jmt, 2
     &,         imt-1, 1, jmt, xt, yt, nnn, work1, lenw)
        else
          write (stdout,*) ' Error: intrp must = 1 or 2'
	 stop
        endif
c
        do j=1,jmt
          block(1,j,k)   = block(imt-1,j,k)
          block(imt,j,k) = block(2,j,k)
        enddo
# else
        call tranlon (data, imd, imd-2, jmd, tdp, cx, xt, imt, cxt)
        if (intrp .eq. 1) then
          call ftc (data, imd, jmd, cxt, cy, block(1,1,k), imt, jmt, 1
     &,       imt, 1, jmt, xt, yt, nnn, work1, lenw)
        else if (intrp .eq. 2) then
          call ctf (data, imd, jmd, cxt, cy, block(1,1,k), imt, jmt, 1
     &,       imt, 1, jmt, xt, yt, nnn, work1, lenw)
        else
          write (stdout,*) ' Error: intrp must = 1 or 2'
	  stop
        endif
# endif
c
c       normalize salilnity to model units
c
        if (lev .eq. 24) then
          do j=1,jmt
            do i=1,imt
              block(i,j,k) = 1.e-3*block(i,j,k)-.035
            enddo
	  enddo
        endif
c
      enddo
c
c-----------------------------------------------------------------------
c    interpolate vertically to MOM levels for each latitude row
c-----------------------------------------------------------------------
c
      do j=1,jmt
        do k=1,kd
          do i=1,imt
            tdv(i,k) = block(i,j,k)
          enddo
	enddo
c
        if (intrpv .eq. 1) then
          call ftc (tdv, imt, kd, xt, dd, tm, imt, km, 1
     &,         imt, 1, km, xt, zt, j, work2, lenw)
        else if (intrpv .eq. 2) then
          call ctf (tdv, imt, kd, xt, dd, tm, imt, km, 1
     &,         imt, 1, km, xt, zt, j, work2, lenw)
        else
          write (stdout,*) ' Error: intrp must = 1 or 2'
	  stop
        endif
c
c       save a slice of surface data
c
        do i=1,imt
          tij(i,j) = tm(i,1)
        enddo
c
        iotext =
     &   'dim x(i),z(k),d(i,k);read(nu)stamp,avg,i,k,j,y,mon,x,z,d'
        write (ioout(month)) iotext
        write (ioout(month)) stamp, dpm, imt, km, j, yt(j), month
     &, xt, zt, tm
        print *,' ==> writing data(imt,km) for jrow =',j, ' at lat='
     &, yt(j),' on ',stamp, ' averaged over preceeding '
     &, dpm,' days'
c
	if (j .eq. jmt/2) then
	  write (6,'(//"j=",i2/10(1pe10.3,1x))')
     &          j,(tm(i,1),i=1,imt)
	  write (6,'(//)')
	endif
      enddo
c
c-----------------------------------------------------------------------
c     store surface data on MOM "t" grid for each month (annual is 13)
c-----------------------------------------------------------------------
c
      k = 1
      iotext =
     & 'dim x(i),y(j),d(i,j);read(nu)stamp,avg,i,j,k,z,mon,x,y,d'
      write (iosurf) iotext
      write (iosurf) stamp, dpm, imt, jmt, k, zt(1), month, xt, yt, tij 
      write (stdout, '(/,a,a,a,g14.7,a,g14.7,a//,a,i3,a,//)')
     & ' ==> writing surface data for ', stamp, ' at depth = ',zt(1)
     &, 'cm,   averaged over previous ', dpm, ' days'
     &,' ====    E N D   OF   M O N T H  ',month,'    ===='
      return
      end
#else
      subroutine ic
      return
      end
#endif
