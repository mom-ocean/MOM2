#ifdef makesbc
# define driver_only
# include "iomngr.F"
# include "grids.F"
# include "size_check.F"
# include "util.F"
      program sbc
c
c-----------------------------------------------------------------------
c     sbc prepares surface boundary conditions using ...
c     Hellerman wind stress monthly climatology and 
c     Oort air temperature monthly climatology
c     for MOM
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
# include "size.h"
      character*60 name
# include "stdunits.h"
# include "coord.h"
c
c     intrph  = 1 when MOM grid resolution is less than Hellerman
c               2 when MOM grid resolution is greater than Hellerman
c     intrpo  = 1 when MOM grid resolution is less than Oort
c               2 when MOM grid resolution is greater than Oort
c
      data  name/' MOM boundary conditions '/
      data intrph /1/, intrpo /1/
      namelist /input/ name, intrph, intrpo
c
# ifdef cyclic
      write (stdout,*) '  ==> cyclic option set'
# else
      write (stdout,*) '  ==> cyclic option NOT set'
# endif
      read  (stdin,input)
      write (stdout,input)
c
      if (intrph .eq. 1) then
        print *,' Assuming MOM resolution is less than Hellerman'
      elseif (intrph .eq. 2) then
        print *,' Assuming MOM resolution is greater than Hellerman'
      else
        print *,' intrph must be 1 or 2'
	stop
      endif
c
      if (intrpo .eq. 1) then
        print *,' Assuming MOM resolution is less than Oort'
      elseif (intrpo .eq. 2) then
        print *,' Assuming MOM resolution is greater than Oort'
      else
        print *,' intrpo must be 1 or 2'
	stop
      endif
c
c-----------------------------------------------------------------------
c     initialize the units
c-----------------------------------------------------------------------
c
      call ioinit
c
c-----------------------------------------------------------------------
c     construct the grid.
c-----------------------------------------------------------------------
c
      write (stdout,'(//1x,a,//)') 'Construct the MOM grid'
c
      maxlen = max(imt,jmt,km)
      call gcoord (maxlen, imt2, jmt2, km2, dxtdeg, dytdeg, dxudeg
     &,            dyudeg, dzt, dzw, xt, xu, yt, yu, zt, zw)
c
      call size_check (imt2, jmt2, km2, 'sbc', 'stop')
c
c-----------------------------------------------------------------------
c     interpolate sol`s monthly wind stress to MOM "u,v" grid
c-----------------------------------------------------------------------
c
      call stress (intrph, xu, yu)
c
c-----------------------------------------------------------------------
c     interpolate oorts monthly air temp to MOM "t" grid
c-----------------------------------------------------------------------
c
      call airtem (intrpo, xt, yt)
c
      call showunits
      call showfiles
      stop
      end 


      subroutine stress (intrp, xu, yu)
c
c-----------------------------------------------------------------------
c     linearly interpolate Sol Hellerman`s monthly wind stress
c     climatology to the MOM grid
c-----------------------------------------------------------------------
c
# include "size.h"
      parameter (ib=182, jb=92, lenw=imt*jmt)
# include "stdunits.h"
      character*32 stamps
      character*80 iotext
      dimension xu(imt), yu(jmt)
      dimension stamps(13), dpm(13)
      dimension at(ib), cxt(ib)
      common/tds/ d(imt,jmt), a(ib,jb), b(ib,jb)
      common /tau/ cx(ib), cy(jb), clmtau(ib,jb,13,2)
      common /tau1/ work(lenw)
c
c-----------------------------------------------------------------------
c     read in Sol Hellerman`s monthly wind stress climatology
c     (the 13th record is the annual mean. k=(1,2) = (taux,tauy))
c-----------------------------------------------------------------------
c
      call getunit (nutauc, 'hellerman.tau', 'u s r ieee')
      do n=1,13
        read (nutauc)
        read (nutauc) stamps(n), dpm(n), lons, lats, nn
     &, cx, cy, (((clmtau(i,j,n,k),i=1,ib),j=1,jb),k=1,2)
      enddo        
c
c-----------------------------------------------------------------------
c     set up the units
c     lrec is in bytes for the cray (words*8)
c-----------------------------------------------------------------------
c
      call getunit (nutaux, 'taux.mom', 'u s r ieee')
      call getunit (nutauy, 'tauy.mom', 'u s r ieee')
c
      mmm    = 0
      do n=1,13
        do m=1,2
	  mmm = mmm + 1
c
c-----------------------------------------------------------------------
c         translate data grid so it entirely encloses the MOM "u,v"
c         grid. interpolate stress data to model "u,v" grid.
c         set last row (=jmt) to zero stress
c-----------------------------------------------------------------------
c
# ifdef cyclic
          call tranlon (clmtau(1,1,n,m), ib, ib-2, jb, at
     &,                 cx, xu(2), imt-2, cxt)
c
          if (intrp .eq. 1) then
            call ftc (clmtau(1,1,n,m), ib, jb, cxt, cy, d, imt, jmt, 2
     &,         imt-1, 1, jmt-1, xu, yu, mmm, work, lenw)
          else if (intrp .eq. 2) then
            call ctf (clmtau(1,1,n,m), ib, jb, cxt, cy, d, imt, jmt, 2
     &,         imt-1, 1, jmt-1, xu, yu, mmm, work, lenw)
          else
            write (stdout,*) ' Error: intrp must = 1 or 2'
	    stop
          endif
          do j=1,jmt
            d(1,j)   = d(imt-1,j)
            d(imt,j) = d(2,j)
          enddo
# else
          call tranlon (clmtau(1,1,n,m), ib, ib-2, jb, at
     &,                cx, xu, imt, cxt)
c
          if (intrp .eq. 1) then
            call ftc (clmtau(1,1,n,m), ib, jb, cxt, cy, d, imt, jmt, 1
     &,         imt, 1, jmt-1, xu, yu, mmm, work, lenw)
          else if (intrp .eq. 2) then
            call ctf (clmtau(1,1,n,m), ib, jb, cxt, cy, d, imt, jmt, 1
     &,         imt, 1, jmt-1, xu, yu, mmm, work, lenw)
          else
            write (stdout,*) ' Error: intrp must = 1 or 2'
	    stop
          endif
# endif
c
          if (n .eq. 1 .and. m .eq. 1) then
            write(stdout,'(///,1x, a38,/)') 
     &       'Hellerman wind stress grid definition:'
            write (stdout,'(//,1x,a30,/)') '            latitudes:'
            write (stdout,'(1x,10f10.5)') cy
            write (stdout,'(//,1x,a30,/)') '  original longitudes:'
            write (stdout,'(1x,10f10.5)') cx
            write (stdout,'(//,1x,a30,/)') 'translated longitudes:'
            write (stdout,'(1x,10f10.5)') cxt
          endif
c
          do i=1,imt
            d(i,jmt) = 0.0
          enddo
c
          if (m .eq. 1) then
         iotext='dim x(i),y(j),tx(i,j);read(nu)stamp,avg,i,j,mon,x,y,tx'
            print *, iotext
            write (nutaux) iotext
            write (nutaux) stamps(n), dpm(n), imt, jmt, n, xu, yu, d
            write (stdout,*) ' ==> taux done ',stamps(n),'  ', dpm(n)
     &,     ' day average'
            call scope (d, imt, imt, jmt, 'taux')
          else
         iotext='dim x(i),y(j),ty(i,j);read(nu)stamp,avg,i,j,mon,x,y,ty'
            print *, iotext
            write (nutauy) iotext
            write (nutauy) stamps(n), dpm(n), imt, jmt, n, xu, yu, d
            write (stdout,*) ' ==> tauy done ',stamps(n),'  ', dpm(n)
     &,     ' day average'
            call scope (d, imt, imt, jmt, 'tauy')
          endif
        enddo
      enddo
      return
      end


      subroutine airtem (intrp, xt, yt)
c
c-----------------------------------------------------------------------
c     interp Oort`s climatological air temperature to MOM "t" grid.
c-----------------------------------------------------------------------
c
# include "size.h"
      parameter (iw=73, jw=75, lenw=imt*jmt)
# include "stdunits.h"
      character*32 stamps
      character*80 iotext
      dimension xt(imt), yt(jmt)
      dimension stamps(13), dpm(13)
      common /airs/ cx(iw), cy(jw), clmair(iw,jw,13)
      common /air1/ work(lenw)
      dimension d(imt,jmt), airt(iw,jw)
      dimension airtp(iw), cxt(iw)
c
c-----------------------------------------------------------------------
c     read in Oort`s monthly air temperature climatology
c     (the 13th record is the annual mean)
c-----------------------------------------------------------------------
c
      call getunit (nuairc, 'oorts.air', 'u s r ieee')
      do n=1,13
        read (nuairc)
        read (nuairc) stamps(n), dpm(n), lons, lats, nn
     &, cx, cy, ((clmair(i,j,n),i=1,iw),j=1,jw)
      enddo        
c
c-----------------------------------------------------------------------
c     set up the units
c     lrec is in bytes for the cray (words*8)
c-----------------------------------------------------------------------
c
      call getunit (nuair, 'airtemp.mom', 'u s r ieee')
c
      do n=1,13
c
c-----------------------------------------------------------------------
c       translate data then interp to MOM "t" grid
c-----------------------------------------------------------------------
c
# ifdef cyclic
        call tranlon (clmair(1,1,n), iw, iw-1, jw, airtp
     &,               cx, xt(2), imt-2, cxt)
c
        if (intrp .eq. 1) then
          call ftc (clmair(1,1,n), iw, jw, cxt, cy, d, imt, jmt, 2
     &,         imt-1, 1, jmt, xt, yt, n, work, lenw)
        else if (intrp .eq. 2) then
          call ctf (clmair(1,1,n), iw, jw, cxt, cy, d, imt, jmt, 2
     &,         imt-1, 1, jmt, xt, yt, n, work, lenw)
        else
          write (stdout,*) ' Error: intrp must = 1 or 2'
	  stop
        endif
c
        do j=1,jmt
          d(1,j)   = d(imt-1,j)
          d(imt,j) = d(2,j)
        enddo
# else
        call tranlon (clmair(1,1,n), iw, iw-1, jw, airtp
     &,               cx, xt, imt, cxt)
        if (intrp .eq. 1) then
          call ftc (clmair(1,1,n), iw, jw, cxt, cy, d, imt, jmt, 1
     &,         imt, 1, jmt, xt, yt, n, work, lenw)
        else if (intrp .eq. 2) then
          call ctf (clmair(1,1,n), iw, jw, cxt, cy, d, imt, jmt, 1
     &,         imt, 1, jmt, xt, yt, n, work, lenw)
        else
          write (stdout,*) ' Error: intrp must = 1 or 2'
	  stop
        endif
# endif
        if (n .eq. 1) then
          write(stdout,'(///,1x, a30,/)') 
     &     'Oort air temp grid definition:'
          write (stdout,'(//,1x,a30,/)') '            latitudes:'
          write (stdout,'(1x,10f10.5)') cy
          write (stdout,'(//,1x,a30,/)') '  original longitudes:'
          write (stdout,'(1x,10f10.5)') cx
          write (stdout,'(//,1x,a30,/)') 'translated longitudes:'
          write (stdout,'(1x,10f10.5)') cxt
        endif
c
      iotext='dim x(i),y(j),at(i,j);read(nu)stamp,avg,i,j,mon,x,y,at'
        print *, iotext
        write (nuair) iotext
        write (nuair) stamps(n), dpm(n), imt, jmt, n, xt, yt, d
        write (stdout,*) ' ==> air temp done ',stamps(n), '  ',dpm(n)
     &, ' day average'
        call scope (d, imt, imt, jmt, 'air temp')
      enddo
      return
      end
#else
      subroutine sbc
      return
      end
#endif
