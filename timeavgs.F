#ifdef drive_timeavgs
      program driver
c
c=======================================================================
c
c              T I M E    M E A N    G R I D    M O D U L E
c
c To generate a grid for constructing time averages: 
c     
c  1) set the model grid in the "grids.F" module.
c
c  2) set the grid specifications in the USER INPUT section below.
c
c  3) compile and run using the run script "run_timeavgs"
c
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
#include "size.h"
#include "coord.h"
      integer cvxz, cvx, cvy, cvz
      parameter (imkmax=imt*km)
      dimension cvxz(imkmax), cvx(imt), cvy(jmt), cvz(km)
      dimension javgr(jmt)
c
      write (stdout,'(//1x,a,//)')  'Constructing MOM grid'
c
c-----------------------------------------------------------------------
c     construct the grid.
c-----------------------------------------------------------------------
c
      maxlen = max(imt, jmt, km)
      call gcoord (maxlen, imt2, jmt2, km2
     &,            dxtdeg, dytdeg, dxudeg, dyudeg
     &,            dzt, dzw, xt, xu, yt, yu, zt, zw)
c
c     check that returned grid sizes match those in file "size.h"
c
      call size_check (imt2, jmt2, km2, 'driver for timemeans', 'stop')
c
c-----------------------------------------------------------------------
c     generate the averaging grid
c-----------------------------------------------------------------------
c
      call avgset (xt, xu, yt, yu, zt, zw, imt*km
     &,            cvxz, cvx, cvy, cvz, javgr, imav, jmav, levav)
c
c-----------------------------------------------------------------------
c     remind user to change imav, jmav, and lenav in "timeavgs.h"
c-----------------------------------------------------------------------
c
      write (stdout,'(//1x,a/)')
     & '==>TIMEMEAN GRID INSTALLATION DIRECTIONS'
c
      write (stdout,'(//,1x,a,/,6x,a,i4,a,i4,a,i3,a,/)')
     & '==> To use in MOM, change parameter in "timeavgs.h" to:'
     &,'parameter (imtav=',imav,', jmtav=',jmav,', kmav=',levav,')'
c
      write (stdout,'(/a/)')
     &'(To change resolution, see the USER INPUT section in timeavgs.F)'
      stop
      end

# define driver_only
# include "grids.F"
# include "iomngr.F"
# include "size_check.F"
# include "util.F"

#endif


      subroutine avgset (xt, xu, yt, yu, zt, zw, imkmax
     &,         cvxz, cvx, cvy, cvz, javgr, imav, jmav, levav)
c
c-----------------------------------------------------------------------
c     setup the "averaging" grid for use with option "time_averages" in
c     MOM for accumulating and saving time averaged data.
c
c     Warning: whenever the definition of this time  averaging
c              grid is changed below, it must be re-installed by 
c              running the "timeavgs.F" module. 
c
c     input:
c       imt   = number of longitudes in the MOM grid
c       jmt   = number of latitudes in the MOM grid
c       km    = number of levels in the MOM grid
c       xt    = longitudinal coordinates of "t" points (deg)
c       xu    = longitudinal coordinates of "u" points (deg)
c       yt    = latitudinal coordinates of "t" points (deg)
c       yu    = latitudinal coordinates of "u" points (deg)
c       zt    = depth of "t" points (cm)
c       zw    = depth of "w" points (bottom of "t" cells) (cm)
c       imkmax= imt*km for dimensioning purposes
c
c     output:
c       cvxz  = indicates points in (i,k) plane for time averaging data
c       cvx   = indicates points in longitude for time averaging data
c       cvy   = indicates points in latitude for time averaging data
c       cvz   = indicates points in depth for time averaging data
c       javgr = number of the latitude row containing averaging points
c       imav  = number of longitudinal points for averaging data
c       jmav  = number of latitudinal points for averaging data
c       levav = number of depth points for averaging data
c     
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
# include "size.h"
# include "stdunits.h"
      logical savex, savey, savez
      integer cvxz, cvx, cvy, cvz
      dimension xt(imt), xu(imt), yt(jmt), yu(jmt), zt(km), zw(km)
      dimension cvxz(imkmax), cvx(imt), cvy(jmt), cvz(km), javgr(jmt)
      dimension savex(imt), savey(jmt), savez(km)
      dimension xtav(imt), xuav(imt), ytav(jmt), yuav(jmt)
      dimension ztav(km), zwav(km)
c
c     initialize all points as not being on the "averaging" grid
c
c
c     savex(i) = (t,f) implies that the model grid point longitude
c                corresponding to grid point index "i" (is, is not)
c                on the "averaging" grid
c     savey(j) = (t,f) implies that the model grid point latitude
c                corresponding to grid point index "j" (is, is not)
c                on the "averaging" grid
c     savez(k) = (t,f) implies that the model grid point depth
c                correspondingto grid point index "k" (is, is not)
c                on the "averaging" grid
c
c     a particular grid point with index (i,j,k) is on the "averaging"
c     grid when "savex(i)", "savey(j)", & "savez(k)" are all true.
c
      do jrow=1,jmt
        savey(jrow) = .false.
      enddo
c
      do i=1,imt
        savex(i) = .false.
      enddo
c
      do k=1,km
        savez(k) = .false.
      enddo
c
c----------------------------------------------------------------------
c     USER INPUT ===> set up the "averaging" grid
c----------------------------------------------------------------------
c
c     for the test case, choose MOM grid points that are nearest to
c     every "dely" degress of latitude. {any (or all) latitudes from
c     jrow=2,jmt-1 on MOM`s grid may be selected)
c
      dely = 10.0
      widlat = 180.0
      mlat = widlat/dely + 1
      do jrow=1,mlat
        m = indp(-90.0 + (jrow-1)*dely, yt, jmt)
        savey(m) = .true.
      enddo
c
c     for the test case, choose MOM grid points that are nearest to
c     every "delx" degress of longitude. {any (or all) longitudes from
c     i=2,imt-1 on MOM`s grid may be selected)
c
      delx = 10.0
      widlon = 720.0
      mlon = widlon/delx + 1
      do i=1,mlon
        m = indp((i-1)*delx, xt, imt)
        savex(m) = .true.
      enddo
c
c     for the test case, choose MOM grid points down to "sdepth" cm.
c     {any (or all) depths from k=1,km on MOM`s grid may be selected)     
c
      sdepth = 500.e2
      nmax   = indp (sdepth, zt, km)
      do n=1,nmax
        savez(n) = .true.
      enddo
c
c----------------------------------------------------------------------
c     ===> end USER INPUT for choosing grid
c----------------------------------------------------------------------
c
      savey(1)   = .false.
      savey(jmt) = .false.
      savex(1)   = .false.
      savex(imt) = .false.
c
c     save indices for jrows which define the "averaging" grid
c
      j = 0
      do jrow=1,jmt
        javgr(jrow) = 0
        if (savey(jrow)) then
          j     = j + 1
          javgr(jrow) = j
        endif
      enddo
c
c----------------------------------------------------------------------
c     calculate the size of the "averaging" grid
c----------------------------------------------------------------------
c
      imav = 0
      do i=1,imt
        if (savex(i)) imav = imav+1
      enddo
c
      jmav = 0
      do jrow=1,jmt
        if (savey(jrow)) jmav = jmav + 1
      enddo
c
      levav = 0
      do k=1,km
        if (savez(k)) levav = levav + 1
      enddo
c
c----------------------------------------------------------------------
c     establish control vectors for gathering data from the MOM grid
c     onto the "averaging" grid
c----------------------------------------------------------------------
c
      lcvxz = 0
      do k=1,km
        if (savez(k)) then
          do i=1,imt
            if (savex(i)) then
              lcvxz = lcvxz + 1
              cvxz(lcvxz) = (k-1)*imt + i
            endif
          enddo
        endif
      enddo
c
      lcvx = 0
      do i=1,imt
        if (savex(i)) then
          lcvx = lcvx + 1
          cvx(lcvx) = i
        endif
      enddo
c
      lcvy = 0
      do jrow=1,jmt
        if (savey(jrow)) then
          lcvy = lcvy + 1
          cvy(lcvy) = jrow
        endif
      enddo
c
      lcvz = 0
      do k=1,km
        if (savez(k)) then
          lcvz = lcvz + 1
          cvz(lcvz) = k
        endif
      enddo
c
      do i=1,imav
        xtav(i) = xt(cvx(i))
        xuav(i) = xu(cvx(i))
      enddo
      do j=1,jmav
        ytav(j) = yt(cvy(j))
        yuav(j) = yu(cvy(j))
      enddo
      do k=1,levav
        ztav(k) = zt(cvz(k))
        zwav(k) = zw(cvz(k))
      enddo
c
c----------------------------------------------------------------------
c     show the "averaging" grid coordinates
c----------------------------------------------------------------------
c
      write (stdout,'(//,20x,a)')
     & 'T I M E   A V E R A G E S   G R I D   I N I T I A L I Z E D'
      write (stdout,'(//1x,a,/)')
     & ' Time averages will be taken at the following grid points:'
      write (stdout,9900) imav,' grid "t" longitudes:'
     &,                   (xtav(i),i=1,imav)
      write (stdout,9900) imav,' grid "v" longitudes:'
     &,                   (xuav(i),i=1,imav)
      write (stdout,9900) jmav,' grid "t" latitudes:'
     &,                   (ytav(j),j=1,jmav)
      write (stdout,9900) jmav,' grid "v" latitudes:'
     &,                   (yuav(j),j=1,jmav)
      write (stdout,9900) levav,' "t" grid depths:'
     &,                   (ztav(k),k=1,levav)
      write (stdout,9900) levav,' "w" grid depths:'
     &,                   (zwav(k),k=1,levav)
c
c---------------------------------------------------------------------
c     compute a grid checksum
c---------------------------------------------------------------------
c
      cksum = 0.0
      cksum = cksum + checksum (xtav, imav, 1)
      cksum = cksum + checksum (ytav, jmav, 1)
      cksum = cksum + checksum (ztav, levav, 1)
      cksum = cksum + checksum (xuav, imav, 1)
      cksum = cksum + checksum (yuav, jmav, 1)
      cksum = cksum + checksum (zwav, levav, 1)
      write (stdout,'(/)')
      write (stdout,*) 'Time average grid checksum = ',cksum
      write (stdout,'(/)')
      return
9900  format (1x,i4,a,/ (1x,10f10.2))
      end

#ifndef drive_timeavgs
      subroutine avgvar (j, jrow, w, varu, vart, stf, smf)
c
c-----------------------------------------------------------------------
c
c     "avgvar" produces a four-dimensional dataset of time averaged
c      quantities on a grid defined by the "timeavgs.F" module.
c      this grid may be the entire MOM grid, a region of it, or any
c      coarser subset of it.  For example:
c      if MOM uses a 2 deg x 2 deg grid, the "averaging" grid
c      could be a 6 deg x 8 deg subset of part or the entire domain.
c      spatial resolution of the "averaging" grid is controlled in
c      the USER INPUT section of "avgset".
c
c      how to use:
c
c      1) set up the "averaging" grid in "avgset" (in this module)
c      2) run this "timeavgs.F" module and follow its directions
c      3) enable the "time_averages" option
c          (possibly the "time_averages_disk" also and set the 
c           averaging period "timavgint" (See "switch.h")
c     
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
#ifdef time_averages
# include "param.h"
# include "timeavgs.h"
c
c     the following parameters are needed to define the workspace area
c
c     n3dvar  = # of 3d fields (# of prognostic variables + one for w)
c     n2dvar  = # of 2d fields ("stf(nt)","smf(2)", and "psi")
c
c     (note: the data structure is similar to "snapshots")
c
# ifdef tcvmix
      parameter (n3dvar=nt+2+2+1, n3dtra=nt+2)
# else
      parameter (n3dvar=nt+2+1, n3dtra=nt)
# endif
      parameter (n3dvam=n3dvar-1, n2dvar=nt+2+1)
      parameter (imkmav=imtav*kmav, imjmav=imtav*jmtav
     &,          lenrow=imkmav*n3dvar+imtav*n2dvar)
c
      common /avgblk/ xtav(imtav), xuav(imtav), ytav(jmtav)
      common /avgblk/ yuav(jmtav), ztav(kmav), zwav(kmav)
      common /avgblk/ avg3d(imkmav,n3dvar), avg2d(imtav,n2dvar)
# ifndef time_averages_disk
      common /avgblk/ spbuf(imkmav,n3dvar,jmtav)
      common /avgblk/ spbuf2(imtav,n2dvar,jmtav)
# endif
      character*12 name3d, name2d
      common /avgblc/ name3d(n3dvar), name2d(n2dvar)
c
c     imav    =  number of points in longitude on the "averaging" grid
c     jmav    =  number of points in latitude on the "averaging" grid
c     levav   =  number of points in depth on the "averaging" grid
c
c     javgr   =  "averaging" grid "jrow" number. it maps "jrow" on the
c                MOM grid into a "jrow" on the "averaging" grid.
c
c     cvxz    =  control vector used in gathering data from array
c                dimensioned (imt,km) into one dimensioned (imav,levav)
c     cvx     =  control vector used in gathering data from array
c                dimensioned (imt) into one dimensioned (imav)
c     cvy     =  control vector used in gathering data from array
c                dimensioned (jmt) into one dimensioned (jmav)
c     cvz     =  control vector used in gathering data from array
c                dimensioned (km) into one dimensioned (levav)
c
c     navgts  =  a counter for tracking the number of time steps
c                 within an averaging period.
c     navgp   =  a counter for tracking the number of averaging
c                 periods.
c
      integer cvxz, cvx, cvy, cvz
      common /avgbli/ cvxz(imtkm), cvx(imt), cvy(jmt), cvz(km)
      common /avgbli/ javgr(jmt), imav, jmav, levav, navgts, navgp
c
# include "coord.h"
# include "docnam.h"
# include "csbc.h"
# include "emode.h"
# include "grdvar.h"
# include "iounit.h"
# include "scalar.h"
# include "switch.h"
# include "tmngr.h"
c
      dimension varu(imtkm,jmw,2), w(imtkm)
      dimension vart(imtkm,jmw,n3dtra)
      dimension stf(imt,jsmw:jemw,nt), smf(imt,jsmw:jemw,2)
c
# ifdef timing
      call tic ('diagnostic', 'time means')
# endif
c
c-----------------------------------------------------------------------
c     only process those "jrows" that are on the "averaging" grid
c-----------------------------------------------------------------------
c
      if (javgr(jrow) .ne. 0) then
c
        jav = javgr(jrow)
c
# ifdef time_averages_disk
c
c-----------------------------------------------------------------------
c       bring in data from SSD for the "averaging" grid row jav
c-----------------------------------------------------------------------
c
        call oget (iotim1, lenrow, jav, avg3d)
c
c-----------------------------------------------------------------------
c       integrate three dimensional data only on those longitudes
c       that are included in the "averaging" grid
c-----------------------------------------------------------------------
c
        do n=1,n3dtra
          do i=1,imkmav
            avg3d(i,n) = avg3d(i,n) + vart(cvxz(i),j,n)
          enddo
	enddo
        do n=1,2
          do i=1,imkmav
	    nvel = n + n3dtra
            avg3d(i,nvel) = avg3d(i,nvel) + varu(cvxz(i),j,n)
          enddo
	enddo
        n = n3dvar
        do i=1,imkmav
          avg3d(i,n) = avg3d(i,n) + w(cvxz(i))
	enddo
c
c-----------------------------------------------------------------------
c       integrate two dimensional fields only on those longitudes
c       that are included in the "averaging" grid
c-----------------------------------------------------------------------
c
        do n=1,n2dvar
          if (n .le. nt) then
            do i=1,imtav
              avg2d(i,n) = avg2d(i,n) + stf(cvx(i),j,n)
            enddo
          elseif (n .lt. n2dvar) then
            do i=1,imtav
              avg2d(i,n) = avg2d(i,n) + smf(cvx(i),j,n-nt)
            enddo
          else
            do i=1,imtav
# if defined rigid_lid_surface_pressure || defined implicit_free_surface
              avg2d(i,n) = avg2d(i,n) + ps(cvx(i),jrow,1)
# endif
# ifdef stream_function
              avg2d(i,n) = avg2d(i,n) + psi(cvx(i),jrow,1)
# endif
            enddo
          endif
        enddo
c
c-----------------------------------------------------------------------
c       write out data to SSD for the "averaging" grid row j
c-----------------------------------------------------------------------
c
        call oput (iotim1, lenrow, jav, avg3d)
# else
c
c-----------------------------------------------------------------------
c       integrate three dimensional data only on those longitudes
c       that are included in the "averaging" grid
c-----------------------------------------------------------------------
c
        do n=1,n3dtra
          do i=1,imkmav
            spbuf(i,n,jav) = spbuf(i,n,jav) + vart(cvxz(i),j,n)
          enddo
	enddo
        do n=1,2
          do i=1,imkmav
	    nvel = n + n3dtra
            spbuf(i,nvel,jav) = spbuf(i,nvel,jav) + varu(cvxz(i),j,n)
          enddo
	enddo
        n = n3dvar
        do i=1,imkmav
          spbuf(i,n,jav) = spbuf(i,n,jav) + w(cvxz(i))
	enddo
c
c-----------------------------------------------------------------------
c       integrate two dimensional fields only on those longitudes
c       that are included in the "averaging" grid
c-----------------------------------------------------------------------
c
        do n=1,n2dvar
          if (n .le. nt) then
            do i=1,imtav
              spbuf2(i,n,jav) = spbuf2(i,n,jav) + stf(cvx(i),j,n)
            enddo
          elseif (n .lt. n2dvar) then
            do i=1,imtav
              spbuf2(i,n,jav) = spbuf2(i,n,jav) + smf(cvx(i),j,n-nt)
            enddo
          else
            do i=1,imtav
# if defined rigid_lid_surface_pressure || defined implicit_free_surface
              spbuf2(i,n,jav) = spbuf2(i,n,jav) + ps(cvx(i),jrow,1)
# endif
#  ifdef stream_function
              spbuf2(i,n,jav) = spbuf2(i,n,jav) + psi(cvx(i),jrow,1)
#  endif
            enddo
          endif
        enddo
# endif
c
c-----------------------------------------------------------------------
c       update integration counter once per time step on the last row
c-----------------------------------------------------------------------
c
        if (jav .eq. jmtav) navgts = navgts + 1
c
      endif
c
# ifdef timing
      call toc ('diagnostic', 'time means')
# endif
c
      return
c
      entry avgout
c
# ifdef timing
      call tic ('diagnostic', 'time means')
# endif
c
c-----------------------------------------------------------------------
c     save the "averaging" grid data at the end of the averaging period
c-----------------------------------------------------------------------
c
      call getunit (iotim, 'time_mean.dta','u s a ieee')
c
      do l=1,imtav
        xtav(l) = xt(cvx(l))
        xuav(l) = xu(cvx(l))
      enddo
      do l=1,jmtav
        ytav(l) = yt(cvy(l))
        yuav(l) = yu(cvy(l))
      enddo
      do l=1,kmav
        ztav(l)  = zt(cvz(l))
        zwav(l)  = zw(cvz(l))
      enddo
      lensp   = jmtav*lenrow + 2*(imtav+jmtav+kmav) + 8
      navgp   = navgp + 1
c
      write (stdout,9999)  navgp, itt, stamp, lensp
9999  format (/1x
     & ,'===> Time means from averaging period #',i4
     &,        ' written to file time_mean.dta',/15x,'on ts = ',i10
     &,        ', ', a, 2x, '(',i10,' words written) ',/)
      reltim = prelyr
      period = dtts*navgts
c
      iotext ='read (iotim) reltim, period, imt, jmt, km, nt, n3d, n2d'
      write (iotim) pstamp, iotext, expnam
      write (iotim) reltim, period, imtav, jmtav, kmav, nt, n3dvar
     &,             n2dvar
c
      iotext = 'read (iotim) (xtav(i),i=1,imtav)'
      write (iotim) pstamp, iotext, expnam
      call wrufio (iotim, xtav, imtav)
c
      iotext = 'read (iotim) (xuav(i),i=1,imtav)'
      write (iotim) pstamp, iotext, expnam
      call wrufio (iotim, xuav, imtav)
c
      iotext = 'read (iotim) (ytav(j),j=1,jmtav)'
      write (iotim) pstamp, iotext, expnam
      call wrufio (iotim, ytav, jmtav)
c
      iotext = 'read (iotim) (yuav(j),j=1,jmtav)'
      write (iotim) pstamp, iotext, expnam
      call wrufio (iotim, yuav, jmtav)
c
      iotext = 'read (iotim) (ztav(k),k=1,kmav)'
      write (iotim) pstamp, iotext, expnam
      call wrufio (iotim, ztav, kmav)
c
      iotext = 'read (iotim) (zwav(k),k=1,kmav)'
      write (iotim) pstamp, iotext, expnam
      call wrufio (iotim, zwav, kmav)
c
      do n=1,nt
        name3d(n) = trname(n)
        if (n .eq. 1) name2d(n) = ' heat flux'
        if (n .eq. 2) name2d(n) = ' salt flux'
        if (n .gt. 2) then
          write (name2d(n),'(a5,i2,a5)') 'trac#',n,' flux'
        endif
      enddo
c
      n = nt+1
      name3d(n) = ' u velocity '
      name2d(n) = ' taux       '
      n = n + 1
      name3d(n) = ' v velocity '
      name2d(n) = ' tauy       '
      n = n + 1
      name3d(n) = ' adv_vbt    '
# if defined stream_function
      name2d(n) = ' psi        '
#endif
# if defined rigid_lid_surface_pressure || defined implicit_free_surface
      name2d(n) = ' surf press '
#endif
c
# ifdef tcvmix
      n = n + 1
      name3d(n) = ' q2         '
      n = n + 1
      name3d(n) = ' vdc        '
# endif
c
      iotext = 'read (iotim) (name3d(1:12)(n),n=1,n3d),'
      iotext(41:) = ' (name2d*(1:12)(n),n=1,n2d)'
      write (iotim) pstamp, iotext, expnam
      write (iotim) name3d, name2d
c
c-----------------------------------------------------------------------
c     construct time mean quantities on the "averaging" grid
c-----------------------------------------------------------------------
c
      rnavgt = c1/navgts
      do jj=1,jmtav
# ifdef time_averages_disk
c
c       bring in data from disk on the "averaging" grid for row jj
c
        call oget (iotim1, lenrow, jj, avg3d)
c
c       construct the time average
c
        do n=1,n3dvar
          do l=1,imkmav
            avg3d(l,n) = rnavgt*avg3d(l,n)
          enddo
	enddo
c
        do n=1,n2dvar
          do l=1,imtav
            avg2d(l,n) = rnavgt*avg2d(l,n)
          enddo
	enddo
# else
c
c       construct the time average
c
        do n=1,n3dvar
          do l=1,imkmav
            avg3d(l,n) = rnavgt*spbuf(l,n,jj)
          enddo
	enddo
c
        do n=1,n2dvar
          do l=1,imtav
            avg2d(l,n) = rnavgt*spbuf2(l,n,jj)
          enddo
	enddo
# endif
c
c-----------------------------------------------------------------------
c       write out data using unformatted i/o
c-----------------------------------------------------------------------
c
        do m=1,nt
          write(iotext,'(a13,i4,a4,i2)') ' for jrowav =',jj,', m=',m
          iotext(24:) =': read (iotim) ((t(i,k,m),i=1,imtav),k=1,kmav)'
          write (iotim) pstamp, iotext, expnam
          call wrufio (iotim, avg3d(1,m), imkmav)
        enddo
c
        write(iotext,'(a13,i4)') ' for jrowav =',jj
        iotext(18:) =': read (iotim) ((u(i,k),i=1,imtav),k=1,kmav)'
        write (iotim) pstamp, iotext, expnam
        call wrufio (iotim, avg3d(1,nt+1), imkmav)
c
        write(iotext,'(a13,i4)') ' for jrowav =',jj
        iotext(18:) =': read (iotim) ((v(i,k),i=1,imtav),k=1,kmav)'
        write (iotim) pstamp, iotext, expnam
        call wrufio (iotim, avg3d(1,nt+2), imkmav)
c
        write(iotext,'(a13,i4)') ' for jrowav =',jj
        iotext(18:) =': read (iotim)((adv_vbt(i,k),i=1,imtav),k=1,kmav)'
        write (iotim) pstamp, iotext, expnam
        call wrufio (iotim, avg3d(1,nt+3), imkmav)
# ifdef tcvmix
c
        write(iotext,'(a13,i4)') ' for jrowav =',jj
        iotext(18:) =': read (iotim) ((q2(i,k),i=1,imtav),k=1,kmav)'
        write (iotim) pstamp, iotext, expnam
        call wrufio (iotim, avg3d(1,nt+4), imkmav)
c
        write(iotext,'(a13,i4)') ' for jrowav =',jj
        iotext(18:) =': read (iotim) ((vdc(i,k),i=1,imtav),k=1,kmav)'
        write (iotim) pstamp, iotext, expnam
        call wrufio (iotim, avg3d(1,nt+5), imkmav)
# endif
c
        write(iotext,'(a13,i4)') ' for jrowav =',jj
        iotext(18:) = ': read (iotim) ((stf(i,n),i=1,imtav),n=1,nt)'
        write (iotim) pstamp, iotext, expnam
        call wrufio (iotim, avg2d(1,1), imtav*nt)
c
        write(iotext,'(a13,i4)') ' for jrowav =',jj
        iotext(18:) = ': read (iotim) ((smf(i,n),i=1,imtav),n=1,2)'
        write (iotim) pstamp, iotext, expnam
        call wrufio (iotim, avg2d(1,nt+1), imtav*2)
c
        write(iotext,'(a13,i4)') ' for jrowav =',jj
# ifdef stream_function 
        iotext(18:) = ': read (iotim) (psi(i,jrow),i=1,imtav)'
# endif
# if defined rigid_lid_surface_pressure || defined implicit_free_surface
        iotext(18:) = ': read (iotim) (ps(i,jrow),i=1,imtav)'
# endif
        write (iotim) pstamp, iotext, expnam
        call wrufio (iotim, avg2d(1,nt+3), imtav)
c
c-----------------------------------------------------------------------
c       zero out the "averaging" data the for the next averaging period
c-----------------------------------------------------------------------
c
# ifdef time_averages_disk
c
        do n=1,n3dvar
          do l=1,imkmav
            avg3d(l,n)  = c0
          enddo
	enddo
c
        do n=1,n2dvar
          do l=1,imtav
            avg2d(l,n) = c0
          enddo
	enddo
c
c       save out the "zeroed" data to SSD on "averaging" grid row jj
c
        call oput (iotim1, lenrow, jj, avg3d)
# else
        do n=1,n3dvar
          do l=1,imkmav
            spbuf(l,n,jj) = c0
          enddo
	enddo
c
        do n=1,n2dvar
          do l=1,imtav
            spbuf2(l,n,jj) = c0
          enddo
	enddo
# endif
      enddo
c
c-----------------------------------------------------------------------
c     zero out the "averaging" counter for the next averaging period
c-----------------------------------------------------------------------
c
      navgts = 0
c
      call relunit (iotim)
# ifdef time_averages_disk
       if (eorun) call oclose (iotim1)
# endif
# ifdef timing
      call toc ('diagnostic', 'time means')
# endif
c
      return
c
      entry avgi
c
#ifdef timing
      call tic ('diagnostic', 'initialization')
#endif
c
c-----------------------------------------------------------------------
c     initialize counters for tracking number of time steps within
c     an averaging period and the number of averaging periods.
c-----------------------------------------------------------------------
c
      navgts = 0
      navgp  = 0
c
c-----------------------------------------------------------------------
c     setup the "averaging" grid (data will be time averaged only on
c     these grid cells)
c-----------------------------------------------------------------------
c
      call avgset (xt, xu, yt, yu, zt, zw, imtkm
     &,         cvxz, cvx, cvy, cvz, javgr, imav, jmav, levav)
c
c-----------------------------------------------------------------------
c     verify that the "averaging" grid size matches parameters set
c     in "timeavgs.h"
c-----------------------------------------------------------------------
c
      if (imav .ne. imtav .or. jmav .ne. jmtav .or. levav .ne. kmav)then
        write (stdout,*) '=> Error: number of grid points in averaging'
     &,' grid does not match the parameter setting in timeavgs.h'
        write (stdout,*) ' imtav=',imtav, ': avgset returns', imav
        write (stdout,*) ' jmtav=',jmtav, ': avgset returns', jmav
        write (stdout,*) '  kmav=',kmav,  ': avgset returns', levav
        stop '=>avgi'
      endif
c
c-----------------------------------------------------------------------
c     initialize time mean averages to "zero"
c-----------------------------------------------------------------------
c
# ifdef time_averages_disk
c
c-----------------------------------------------------------------------
c     open unit for "averaging" grid data
c     nbuf = number of i/o buffers (currently not used) 
c-----------------------------------------------------------------------
c
# ifdef ramdrive
      print '(/a/a)'
     &,'=>Error: "time_averages_disk" is not allowed with "ramdrive"'
     &,'          response: disable the "time_averages_disk" option'
      stop '=>avgi'
# endif
      write (stdout,'(/a/)')
     &     'Space for "time_averages_disk" file "iotim1" follows:'
      nbuf = 2
      call ostart (iotim1, 'iotim1', jmtav*lenrow, lenrow, nbuf)
c
      do n=1,n3dvar
        do l=1,imkmav
          avg3d(l,n) = c0
        enddo
      enddo
c
      do n=1,n2dvar
        do l=1,imtav
          avg2d(l,n) = c0
        enddo
      enddo
c
      do jj=1,jmtav
        call oput (iotim1, lenrow, jj, avg3d)
      enddo
# else
      do jj=1,jmtav
        do n=1,n3dvar
          do l=1,imkmav
            spbuf(l,n,jj) = c0
          enddo
	enddo
c
        do n=1,n2dvar
          do l=1,imtav
            spbuf2(l,n,jj) = c0
          enddo
        enddo
      enddo
# endif
#endif
c
#ifdef timing
      call toc ('diagnostic', 'initialization')
#endif
c
      return
      end
#endif
