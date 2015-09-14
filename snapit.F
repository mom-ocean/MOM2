      subroutine snapit (j, jrow)
#ifdef snapshots 
# if !defined netcdf 
c
c=======================================================================
c     save data for an instantaneous snapshot. preceed each data record
c     by a character record describing the data record. note that the
c     latitude and depth range for saving data is in "csnap.h" and
c     is initialized via "blkdta.F"
c     
c     stucture of data saved is similar to "time_averages" 
c
c     author:        r.c.pacanowski    e-mail=> rcp@gfdl.gov
c                    a.rosati          e-mail=> ar@gfdl.gov
c=======================================================================
c
c
# include "param.h"
# include "coord.h"
# include "csnap.h"
# include "docnam.h"
# include "emode.h"
# include "iounit.h"
# ifdef isopycmix
#  include "isopyc.h"
# endif
# include "mw.h"
# include "tmngr.h"
c
c     n3d = number of 3d fields (# of prognostic variables + one for
c           adv_vbt
c     n2d = number of 2d fields (stf(nt) + smf(2) + psi)
c
# ifdef tcvmix
      parameter (n3d=nt+2+2+1)
# else
      parameter (n3d=nt+2+1)
# endif
      parameter (n2d=nt+2+1)
      character*12 snap3d, snap2d
      common /snapc/ snap3d(n3d), snap2d(n2d)
c
c     jssnap = nearest jrow corresponding to "snapls"
c     jesnap = nearest jrow corresponding to "snaple"
c     kesnap = nearest k level corresponding to "snapde"
c     jl     = number of jrows to be saved
c     kl     = number of depths to be saved
c
      common /csnapi/ jssnap, jesnap, kesnap, jl, kl
c
# ifdef timing
      call tic ('diagnostic', 'snapshots')
# endif
c
      call getunit (ios, 'snapshots.dta','u s a ieee')
c
      if (jrow .eq. 2) then
c
c       calculate range of "jrows" and "k" levels for snapshots
c       row "jmt" isn`t allowed
c
        jssnap = indp (snapls, yt, jmt)
        jesnap = indp (snaple, yt, jmt)
        if (jesnap .eq. jmt) jesnap = jesnap-1
        kesnap = indp (snapde, zt, km)
        snapls = yt(jssnap)
        snaple = yt(jesnap)
        snapde = zt(kesnap)
        write (stdout,*)
     &       ' ==>  Writing a snapshot of MOM data in file '
     &,            'snapshots.dta, ts = ', itt, stamp
        write (stdout,*)
     &       '      between latitudes:',snapls,' and ', snaple
     &,      ' and from the surface down to depth  ',snapde   
c
        reltim = prelyr
        period = c0
        jl     = jesnap - jssnap + 1
        kl     = kesnap
c
        iotext =
     & 'read(ios) reltim, period, imt, jmt, km, nt, n3d, n2d'
        write (ios) pstamp, iotext, expnam
        write (ios) reltim, period, imt, jl, kl, nt, n3d, n2d
c
        iotext = 'read(ios) (xt(i),i=1,imt)'
        write (ios) pstamp, iotext, expnam
        call wrufio (ios, xt, imt)
c
        iotext = 'read(ios) (xu(i),i=1,imt)'
        write (ios) pstamp, iotext, expnam
        call wrufio (ios, xu, imt)
c
        iotext = 'read(ios) (yt(j),j=1,jmt)'
        write (ios) pstamp, iotext, expnam
        call wrufio (ios, yt(jssnap), jl)
c
        iotext = 'read(ios) (yu(j),j=1,jmt)'
        write (ios) pstamp, iotext, expnam
        call wrufio (ios, yu(jssnap), jl)
c
        iotext = 'read(ios) (zt(k),k=1,km)'
        write (ios) pstamp, iotext, expnam
        call wrufio (ios, zt, kl)
c
        iotext = 'read(ios) (zw(k),k=1,km)'
        write (ios) pstamp, iotext, expnam
        call wrufio (ios, zw, kl)
c
        do n=1,nt
          snap3d(n) = trname(n)
          if (n .eq. 1) snap2d(n) = ' heat flux'
          if (n .eq. 2) snap2d(n) = ' salt flux'
          if (n .gt. 2) then
            write (snap2d(n),'(a5,i2,a5)') 'trac#',n,' flux'
          endif
        enddo
c
        n = nt+1
        snap3d(n) = ' u velocity '
        snap2d(n) = ' taux       '
        n = n + 1
        snap3d(n) = ' v velocity '
        snap2d(n) = ' tauy       '
        n = n + 1
        snap3d(n) = ' adv_vbt    '
# if defined stream_function
        snap2d(n) = ' psi        '
#endif
# if defined rigid_lid_surface_pressure || defined implicit_free_surface
        snap2d(n) = ' surf press '
#endif
c
# ifdef tcvmix
        n = n + 1
        snap3d(n) = ' q2         '
        n = n + 1
        snap3d(n) = ' vdc        '
# endif
c
        iotext = 'read (ios) (name3d(1:12)(n),n=1,n3d),'
        iotext(41:) = ' (name2d*(1:12)(n),n=1,n2d)'
        write (ios) pstamp, iotext, expnam
        write (ios) snap3d, snap2d
      endif
c
      if (jrow .lt. jssnap .or. jrow .gt. jesnap) goto 1001
c
c     save data from this row to depth given by "kl"
c
      do m=1,nt
        write(iotext,'(a10,i4,a4,i2)') ' for jrow=',jrow,', m=',m
        iotext(21:) = ': read (ios) (t(i,k,m),i=1,imt),k=1,km)'
        write (ios) pstamp, iotext, expnam
        call wrufio (ios, t(1,1,j,m,tau), imt*kl)
      enddo
c
      write(iotext,'(a10,i4)') ' for jrow=',jrow
      iotext(15:) = ': read (ios) (u(i,k),i=1,imt),k=1,km)'
      write (ios) pstamp, iotext, expnam
      call wrufio (ios, u(1,1,j,1,tau), imt*kl)
c
      write(iotext,'(a10,i4)') ' for jrow=',jrow
      iotext(15:) = ': read (ios) (v(i,k),i=1,imt),k=1,km)'
      write (ios) pstamp, iotext, expnam
      call wrufio (ios, u(1,1,j,2,tau), imt*kl)
c
      write(iotext,'(a10,i4)') ' for jrow=',jrow
      iotext(15:) = ': read (ios) (adv_vbt(i,k),i=1,imt),k=1,km)'
      write (ios) pstamp, iotext, expnam
      call wrufio (ios, adv_vbt(1,1,j), imt*kl)
# if defined isopycmix && defined gent_mcwilliams
c
      write(iotext,'(a10,i4)') ' for jrow=',jrow
      iotext(15:) = ': read (ios) (adv_vetiso(i,k),i=1,imt),k=1,km)'
      write (ios) pstamp, iotext, expnam
      call wrufio (ios, adv_vetiso(1,1,j), imt*kl)
c
      write(iotext,'(a10,i4)') ' for jrow=',jrow
      iotext(15:) = ': read (ios) (adv_vntiso(i,k),i=1,imt),k=1,km)'
      write (ios) pstamp, iotext, expnam
      call wrufio (ios, adv_vntiso(1,1,j), imt*kl)
c
      write(iotext,'(a10,i4)') ' for jrow=',jrow
      iotext(15:) = ': read (ios) (adv_vbtiso(i,k),i=1,imt),k=1,km)'
      write (ios) pstamp, iotext, expnam
      call wrufio (ios, adv_vbtiso(1,1,j), imt*kl)
# endif
# ifdef tcvmix
c
      write(iotext,'(a10,i4)') ' for jrow=',jrow
      iotext(15:) = ': read (ios) (q2(i,k),i=1,imt),k=1,km)'
      write (ios) pstamp, iotext, expnam
      call wrufio (ios, q2(1,1,j), imt*kl)
c
      write(iotext,'(a10,i4)') ' for jrow=',jrow
      iotext(15:) = ': read (ios) (vdc(i,k),i=1,imt),k=1,km)'
      write (ios) pstamp, iotext, expnam
      call wrufio (ios, vdc(1,1,j), imt*kl)
# endif
c
      write(iotext,'(a10,i4)') ' for jrow=',jrow
      iotext(15:) = ': read (ios) ((stf(i,n),i=1,imt),n=1,nt)'
      write (ios) pstamp, iotext, expnam
      write (ios) ((stf(l,j,m),l=1,imt),m=1,nt)
c      call wrufio (ios, stf(1,j,n), imt*nt)
c
      write(iotext,'(a10,i4)') ' for jrow=',jrow
      iotext(15:) = ': read (ios) ((smf(i,n),i=1,imt),n=1,2)'
      write (ios) pstamp, iotext, expnam
      write (ios) ((smf(l,j,m),l=1,imt),m=1,2)
c      call wrufio (ios, smf, imt*2)
c
# if defined rigid_lid_surface_pressure || defined implicit_free_surface
      write(iotext,'(a10,i4)') ' for jrow=',jrow
      iotext(15:) = ': read (ios) (ps(i,jrow),i=1,imt)'
      write (ios) pstamp, iotext, expnam
c
c     ps is at "tau" since the external mode has not been updated yet.
c
      call wrufio (ios, ps(1,jrow,1), imt)
# endif
# ifdef stream_function 
      write(iotext,'(a10,i4)') ' for jrow=',jrow
      iotext(15:) = ': read (ios) (psi(i,jrow),i=1,imt)'
      write (ios) pstamp, iotext, expnam
c
c     psi is at "tau" since the external mode has not been updated yet.
c
      call wrufio (ios, psi(1,jrow,1), imt)
# endif
c
1001  continue
      call relunit (ios)
# ifdef timing
      call toc ('diagnostic', 'snapshots')
# endif
# endif
# ifdef netcdf
c
c=======================================================================
c     save data for an instantaneous snapshot.  Using netcdf format.
c     note that the latitude and depth range for saving data is in
c     "csnap.h" and is initialized via "blkdta.F"
c
c    author:    Matt Harrison    e-mail =>   mjh@gfdl.gov
c=======================================================================
c
c
# include "param.h"
# include "coord.h"
# include "csnap.h"
# include "docnam.h"
# include "emode.h"
# include "iounit.h"
# include "mw.h"
# include "tmngr.h"
c
c     n3d = number of 3d fields (# of prognostic variables + one for w)
c     n2d = number of 2d fields (stf(nt) + smf(2) + psi)
c
# ifdef tcvmix
      parameter (n3d=nt+2+2+1)
# else
      parameter (n3d=nt+2+1)
# endif
      parameter (n2d=nt+2+1)
      character*12 snap3d, snap2d
      character*4 cyy
      character*2 cmm,cdd
      common /snapc/ snap3d(n3d), snap2d(n2d)
c
c     jssnap = nearest jrow corresponding to "snapls"
c     jesnap = nearest jrow corresponding to "snaple"
c     kesnap = nearest k level corresponding to "snapde"
c     jl     = number of jrows to be saved
c     kl     = number of depths to be saved
c
      common /csnapi/ jssnap, jesnap, kesnap, jl, kl
c
c--------------------------------------------------------------------------
c     ndims = number of axes (t grid and u,v grid)
c     mxdimin,mxdimout = maximum length of (model,snapshot file) grids
c     nvars = number of variables (not including axis coordinates) being
c     written
c--------------------------------------------------------------------------
c
      parameter (ndims=6)
      parameter (mxdimin=400,mxdimout=400)
      parameter (nvars = n3d + n2d)
c
# include "nc_stuff.h"
c
c--------------------------------------------------------------------------
c     maxndims = maximum number of axes in netcdf file (>=ndims)
c     maxnvars = maximum number of variables in netcdf file (>=nvars)
c     cdfid    = unit number of netcdf file (not using i/o manager for
c     this routine ndims1,id_dim,nvars1,id_var,irecdim = internal
c     variables for netcdf routines
c--------------------------------------------------------------------------
c
      parameter (maxndims=maxncdim, maxnvars=maxncvar)
      common /ncidstuff/ cdfid, ndims1, id_dim(maxndims), 
     &                            nvars1, id_var(maxnvars),
     &			          irecdim
c
c
# ifdef timing
      call tic ('diagnostic', 'snapshots')
# endif
c
      if (jrow .eq. 2) then
c
c---------------------------------------------------------------------------
c       Determine all of the variables needed to define NetCDF file
c       for jrow = 2, define snapshot filename and time header 
c       information from "stamp" variable.  Currently writing netcdf
c       files using one time level per file.
c---------------------------------------------------------------------------
c
c---------------------------------------------------------------------------
c       Name netcdf snapshot file
c---------------------------------------------------------------------------
c
        cyy=stamp(13:16)
        cmm=stamp(7:8)
        cdd=stamp(10:11)
        fname = 'snapshot.'//cyy(NON_BLANK(cyy):4)//'.'
     &           //cmm(NON_BLANK(cmm):2)//'.'
     &          //cdd(NON_BLANK(cdd):2)//'.nc'
c
c---------------------------------------------------------------------------
c       Define time header information using numbers derived from
c       "stamp"
c---------------------------------------------------------------------------
c
        irecdimtyp = 2
        call rdstmp(stamp,iyy,imm,idd,ihh,imin,isec)
c
c---------------------------------------------------------------------------
c       lclobber =  Rewrite (Do not write to) existing netcdf files
c       gtitle = generic name of netcdf file
c       gspval = global undefined value
c---------------------------------------------------------------------------
c
        lclobber = .true.
        gtitle = '3d and 2d mom_2 snapshot data'
        lgspval = .false.
        gspval = -1.e+10
        cdatim = '  '
c
c--------------------------------------------------------------------------
c       Define time-independent dimension information
c--------------------------------------------------------------------------
c
        cdimnam(1) = 'Longitude_t'
        cdimnam(2) = 'Longitude_u'
        cdimnam(3) = 'Depth_t'
        cdimnam(4) = 'Depth_w'
        cdimnam(5) = 'Latitude_t'
        cdimnam(6) = 'Latitude_u'
c
        cdimlnam(1) = 'Longitude on T grid '
        cdimlnam(2) = 'Longitude on U grid '
        cdimlnam(3) = 'Depth of T grid points '
        cdimlnam(4) = 'Depth at bottom of T boxes'
        cdimlnam(5) = 'Latitude on T grid '
        cdimlnam(6) = 'Latitude on U grid '
c
        cdimunits(1) = 'degrees_W'
        cdimunits(2) = 'degrees_W'
        cdimunits(3) = 'centimeters'
        cdimunits(4) = 'centimeters'
        cdimunits(5) = 'degrees_N'
        cdimunits(6) = 'degrees_N'
c 
c--------------------------------------------------------------------------
c       Define Model variable information.  See ncsetup.F for variable
c       definitions.
c--------------------------------------------------------------------------
c
c--------------------------------------------------------------------------
c      Tracers
c--------------------------------------------------------------------------
c
        do n=1,nt 
          cvarnam(n) = trname(n)
          cvarlnam(n) = trname(n)
          if (n.eq.1) then
            cvarunits(n) = 'deg C'
            nvdims(n) = 3
            idimindx(1,n) = 1
            idimindx(2,n) = 3
            idimindx(3,n) = 5
            validr(1,n) = -5.0
            validr(2,n) = 50.0
            lspval(n) = .false.
            spval(n) = 0.0
          endif
          if (n.eq.2) then
            cvarunits(n) = 'PSU'
            nvdims(n) = 3
            idimindx(1,n) = 1
            idimindx(2,n) = 3
            idimindx(3,n) = 5
            validr(1,n) = -100.0
            validr(2,n) = 100.0
            lspval(n) = .false.
            spval(n) = 0.0
          endif
        enddo
c
c--------------------------------------------------------------------------
c       Velocities
c--------------------------------------------------------------------------
c 
        cvarnam(nt+1) = 'u_velocity ' 
        cvarlnam(nt+1) = 'u_velocity '   
        cvarunits(nt+1) = 'cm/sec'
        nvdims(nt+1) = 3
        idimindx(1,nt+1) = 2
        idimindx(2,nt+1) = 4
        idimindx(3,nt+1) = 6
        validr(1,nt+1) = -500.0
        validr(2,nt+1) = 500.0
        lspval(nt+1) = .false.
        spval(nt+1) = 0.0

        cvarnam(nt+2) = 'v_velocity ' 
        cvarlnam(nt+2) = 'v_velocity ' 
        cvarunits(nt+2) = 'cm/sec'
        nvdims(nt+2) = 3
        idimindx(1,nt+2) = 2
        idimindx(2,nt+2) = 4
        idimindx(3,nt+2) = 6
        validr(1,nt+2) = -500.0
        validr(2,nt+2) = 500.0
        lspval(nt+2) = .false.
        spval(nt+2) = 0.0

        cvarnam(nt+3) = 'adv_vbt    '
        cvarlnam(nt+3) = 'adv_vbt    '
        cvarunits(nt+3) = 'cm/sec'
        nvdims(nt+3) = 3
        idimindx(1,nt+3) = 2
        idimindx(2,nt+3) = 4
        idimindx(3,nt+3) = 6
        validr(1,nt+3) = -500.0
        validr(2,nt+3) = 500.0
        lspval(nt+3) = .false.
        spval(nt+3) = 0.0
c
c--------------------------------------------------------------------------
c     Additional three dimensional valiables
c--------------------------------------------------------------------------
c
        n = nt+3
# ifdef tcvmix
        n = n + 1
        cvarnam(n) = ' q2         '
        cvarlnam(n) = ' q2         '
        cvarunits(n) = '  '
        nvdims(n) = 3
        idimindx(1,n) = 2
        idimindx(2,n) = 4
        idimindx(3,n) = 6 
        validr(1,n) = 0.0
        validr(2,n) = 0.0
        lspval(n) = .false.
        spval(n) = 0.0
c
        n = n + 1
        cvarnam(n) = ' vdc        '
        cvarlnam(n) = ' vdc         '
        cvarunits(n) = '  '
        nvdims(n) = 3
        idimindx(1,n) = 2
        idimindx(2,n) = 4
        idimindx(3,n) = 6
        validr(1,n) = 0.0
        validr(2,n) = 0.0
        lspval(n) = .false.
        spval(n) = 0.0
# endif
c
c--------------------------------------------------------------------------
c     Two-dimensional variables
c--------------------------------------------------------------------------
c
        n = n + 1 
        cvarnam(n) = 'heat_flux'
        cvarlnam(n) = 'heat_flux'
        cvarunits(n) = 'cal/cm**2/sec'
        nvdims(n) = 2
        idimindx(1,n) = 1
        idimindx(2,n) = 5
        validr(1,n) = 0.0
        validr(2,n) = 0.0
        lspval(n) = .false.
        spval(n) = 0.0
c
        n = n + 1
        cvarnam(n) = 'salt_flux'
        cvarlnam(n) = 'salt_flux'
        cvarunits(n) = 'gm/cm**2/sec'
        nvdims(n) = 2
        idimindx(1,n) = 1
        idimindx(2,n) = 5
        validr(1,n) = 0.0
        validr(2,n) = 0.0
        lspval(n) = .false.
        spval(n) = 0.0
c
        n = n + 1
        cvarnam(n) = 'taux       '
        cvarlnam(n) = 'taux       '
        cvarunits(n) = 'dynes/cm**2'
        nvdims(n) = 2
        idimindx(1,n) = 2
        idimindx(2,n) = 6
        validr(1,n) = 0.0
        validr(2,n) = 0.0
        lspval(n) = .false.
        spval(n) = 0.0

        n = n + 1
        cvarnam(n) = 'tauy       '
        cvarlnam(n) = 'tauy       '
        cvarunits(n) = 'dynes/cm**2' 
        nvdims(n) = 2
        idimindx(1,n) = 2
        idimindx(2,n) = 6
        validr(1,n) = 0.0
        validr(2,n) = 0.0
        lspval(n) = .false.
        spval(n) = 0.0
c
# if defined stream_function
        n = n + 1
        cvarnam(n) = 'psi        '
        cvarlnam(n) = 'psi        '
        cvarunits(n) = 'cm**3/sec'
        nvdims(n) = 2
        idimindx(1,n) = 1
        idimindx(2,n) = 5
        validr(1,n) = 0.0
        validr(2,n) = 0.0
        lspval(n) = .false.
        spval(n) = 0.0
#endif
# if defined rigid_lid_surface_pressure || defined implicit_free_surface
        n = n + 1
        cvarnam(n) = ' surf press '
        cvarlnam(n) = ' surf press '
        cvarunits(n) = 'centimeters'
        nvdims(n) = 2
        idimindx(1,n) = 1
        idimindx(2,n) = 5
        validr(1,n) = 0.0
        validr(2,n) = 0.0
        lspval(n) = .false.
        spval(n) = 0.0
#endif
c
c---------------------------------------------------------------------------
c    Finished defining variable information for netcdf header
c---------------------------------------------------------------------------
c
c       calculate range of "jrows" and "k" levels for snapshots
c       row "jmt" isn`t allowed
c
        jssnap = indp (snapls, yt, jmt)
        jesnap = indp (snaple, yt, jmt)
        if (jesnap .eq. jmt) jesnap = jesnap-1
        if (jssnap .eq. 1) jssnap = 2
        kesnap = indp (snapde, zt, km)
        snapls = yt(jssnap)
        snaple = yt(jesnap)
        snapde = zt(kesnap) 
c
        jl   = jesnap - jssnap + 1
        kl   = kesnap  
c
c       assign values for lendims array
c
        lendims(1) = imt
        lendims(2) = imt 
        lendims(3) = kl
        lendims(4) = kl
        lendims(5) = jl
        lendims(6) = jl
c
c       define array of axis values
c
        do i = 1,imt
          dimvals(i,1) = xt(i)
          dimvals(i,2) = xu(i)
        enddo
c
        do k = 1,kl
          dimvals(k,3) = zt(k)
          dimvals(k,4) = zw(k)
        enddo
c
        do jj = jssnap,jesnap
          dimvals(jj-jssnap+1,5) = yt(jj)
          dimvals(jj-jssnap+1,6) = yu(jj)
        enddo
c
c------------------------------------------------------------------------
c      After defining variables, call routine to setup NetCDF file 
c------------------------------------------------------------------------
c
        call  ncsetup ( fname, lclobber, gtitle, lgspval,gspval,
     &			irecdimtyp, cdatim,
     &                  iyy,imm,idd,ihh,imin,isec,
     &                  iyr_refmin, imo_refmin,  ida_refmin, 
     &                  ihr_refmin, imin_refmin, isec_refmin,
     &                  ndims, lendims, cdimnam, cdimlnam,
     &                  cdimunits, dimvals, mxdimout,
     &                  nvars, cvarnam, cvarlnam, cvarunits, 
     &			nvdims, idimindx, ndims,
     &                  validr, lspval, vspval )
c
        write (stdout,*)
     &       ' ==>  Writing a snapshot of MOM data on unit '
     &,            cdfid, ', ts = ', itt, stamp
        write (stdout,*)
     &       '      between latitudes:',snapls,' and ', snaple
     &,      ' and from the surface down to depth  ',snapde
      endif
c
c
      if (jrow .lt. jssnap .or. jrow .gt. jesnap) goto 1001
c
      lendims(1) = imt
      lendims(2) = imt 
      lendims(3) = kl
      lendims(4) = kl
      lendims(5) = jl
      lendims(6) = jl
c
c------------------------------------------------------------------------
c     save data from this row to depth given by "kl"
c     Using NCSTORE routine
c------------------------------------------------------------------------
c
      ivar = 1
      do m=1,nt
        istart(1) = 1
        istart(2) = 1
        istart(3) = jrow-jssnap+1
        icount(1) = lendims(1)
        icount(2) = lendims(3)
        icount(3) = 1
        call ncstore (ivar, istart, icount, t(1,1,j,m,tau), recval)
        ivar = ivar + 1
      enddo
c
      istart(1) = 1
      istart(2) = 1
      istart(3) = jrow-jssnap+1
      icount(1) = lendims(2)
      icount(2) = lendims(4)
      icount(3) = 1 
      call ncstore (ivar, istart,icount,u(1,1,j,1,tau), recval)
      ivar = ivar + 1
c
      istart(1) = 1
      istart(2) = 1
      istart(3) = jrow-jssnap+1
      icount(1) = lendims(2)
      icount(2) = lendims(4)
      icount(3) = 1
      call ncstore (ivar, istart, icount,u(1,1,j,2,tau), recval)
      ivar = ivar + 1
c
      istart(1) = 1
      istart(2) = 1
      istart(3) = jrow-jssnap+1
      icount(1) = lendims(1)
      icount(2) = lendims(3)
      icount(3) = 1
      call ncstore (ivar, istart, icount, adv_vbt(1,1,j), recval)
      ivar = ivar + 1
# ifdef tcvmix
c
      istart(1) = 1
      istart(2) = 1
      istart(3) = jrow-jssnap+1
      icount(1) = lendims(1)
      icount(2) = lendims(3)
      icount(3) = 1
      call ncstore (ivar, istart, icount, q2(1,1,j), recval)
      ivar = ivar + 1
c
      istart(1) = 1
      istart(2) = 1
      istart(3) = jrow-jssnap+1
      icount(1) = lendims(1)
      icount(2) = lendims(3)
      icount(3) = 1
      call ncstore (ivar, istart,icount,vdc(1,1,j), recval)
      ivar = ivar + 1
# endif
c
      do m = 1,nt
        istart(1) = 1
        istart(2) = jrow-jssnap+1
        icount(1) = lendims(1)
        icount(2) = 1
        call ncstore (ivar, istart, icount,stf(1,j,m), recval)
        ivar = ivar + 1
      enddo
c
c
      do m = 1,2
        istart(1) = 1
        istart(2) = jrow-jssnap+1
        icount(1) = lendims(2)
        icount(2) = 1
        call ncstore (ivar, istart, icount, stf(1,j,m), recval)
        ivar = ivar + 1
      enddo
c
# if defined rigid_lid_surface_pressure || defined implicit_free_surface
c
c     ps is at "tau" since the external mode has not been updated yet.
c
      istart(1) = 1
      istart(2) = jrow-jssnap+1
      icount(1) = lendims(1)
      icount(2) = 1
      call ncstore (ivar, istart,icount,ps(1,jrow,1), recval)
      ivar = ivar + 1
c
# endif
# ifdef stream_function 
c
c     psi is at "tau" since the external mode has not been updated yet.
c
      istart(1) = 1
      istart(2) = jrow-jssnap+1
      icount(1) = lendims(1)
      icount(2) = 1
      call ncstore (ivar,istart, icount, psi(1,jrow,1), recval)
# endif
c
1001  continue
c
c-----------------------------------------------------------------------
c     Finished writing NetCDF snapshot file 
c-----------------------------------------------------------------------
c
      if (jrow .eq. jesnap) then
         call ncdone
      endif
c
# ifdef timing
      call toc ('diagnostic', 'snapshots')      
# endif
#endif
#endif
      return
      end


#ifdef netcdf
      function NON_BLANK(string)
c
c     defines length of "string" ignoring leading blanks
c
      character * (*) string
      integer NON_BLANK, i
c
      NON_BLANK = 0
      do i =  1,LEN (string)
        if (string(i:i) .NE. ' ') then
          NON_BLANK = i
          goto 999
        end if
      end do
999   continue
c      
      return
      end
#endif

