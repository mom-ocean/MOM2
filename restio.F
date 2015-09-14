      subroutine read_restart
c
c=======================================================================
c
c     read the restart file
c     
c     author:        a. rosati             e-mail=> ar@gfdl.gov
c                    r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#ifdef trajectories
# include "ptraj.h"
#endif
#include "coord.h"
#include "emode.h"
#include "grdvar.h"
#include "iounit.h"
#include "levind.h"
#include "mw.h"
#include "tmngr.h"
      character*32 ostamp
      dimension bufsl(imt,km,nvar)
c
c-----------------------------------------------------------------------
c     read restart file (unformatted)
c
c     record 1 on the restart file is the control block (containing
c              the time step and relative year)
c     record`s 2 & 3 are for stream function time levels
c     record`s 4 & 5 are for the d(stream function)/dt fields used
c              to construct a guess for the elliptic solver
c     record 6 is for "kmt" indicating the number of ocean "t" cells
c              from the surface to the ocean floor.
c     the next "jmt" records contain prognostic variables at "tau"
c     the last "jmt" records contain prognostic variables at "tau+1"
c-----------------------------------------------------------------------
c
      call tic ('ocean', 'reading restart')
c
      call getunit (iorest, 'restart.dta','u s r')
c
c-----------------------------------------------------------------------
c     read the time step and stamp from the restart file
c     timestep "itt" at this point corresponds to "tau+1" which was the
c     condition at the end of the previous timestep
c-----------------------------------------------------------------------
c
      read (iorest) stamp, iotext
      read (iorest) itt, irstdy, msrsdy, im, jm, kkm
      write (stdout,'(a,i10,1x,a)')
     &'===> Reading MOM restart from file restart.dta on ts=',itt, stamp
      if (im .ne. imt .or. jm .ne. jmt .or. kkm .ne. km) then
        write (stdout,*) '=>Error restio. im=',im,' jm=',jm,' km=',kkm
      endif
c
c-----------------------------------------------------------------------
c     verify that this restart file is the right one. The time step
c     and stamp must match those in the "archive.time" file which the
c     run script should restore upon restarting. (security measure)
c-----------------------------------------------------------------------
c
      ittarc = 0
      ostamp = 'no archive.temp file'
      call getunit (ioarch, 'archive.time','f s r')
      do n=1,100000
        read (ioarch,8888,end=100) ittarc, ostamp
	write (stdout,'(1x,a,1x,i10,1x,a32)')
     & ' read archive.time record for:', ittarc, ostamp
      enddo
100   continue
      call relunit (ioarch)
8888  format (1x,i10,1x,a32)
      if (itt .ne. ittarc .or. stamp .ne. ostamp) then
        if (ittarc .eq. 0) then
          write (stdout,*) '=>Error: could not find "archive.time" file'
     &,   ' to verify this is the correct restart.'
          write (stdout,*) ' restart itt=',itt,' restart stamp = ',stamp
	else
	  write (stdout,*) '=>Error: archive.time itt=',ittarc
     &,                    ' but restart has itt = ',itt	  
	  write (stdout,*) '         archive.time stamp=',ostamp
     &,                    ' but restart stamp = ',stamp
	endif
        write (stdout,*) '=>Error in restio.F'
        stop	  
      endif
c
c-----------------------------------------------------------------------
c     update pointers to tau-1, tau, & tau+1 data on disk based on itt
c     for latitude rows they point to latdisk(1) or latdisk(2)
c     for 2D fields they point to records on kflds
c-----------------------------------------------------------------------
c
      taum1disk = mod(itt+1,2) + 1
      taudisk   = mod(itt  ,2) + 1
      taup1disk = taum1disk
c
c-----------------------------------------------------------------------
c     update pointers to tau-1, tau, & tau+1 data in the MW based on itt
c-----------------------------------------------------------------------
c
      if (wide_open_mw) then
c
c       rotate time levels instead of moving data
c
        taum1 = mod(itt+0,3) - 1
        tau   = mod(itt+1,3) - 1
        taup1 = mod(itt+2,3) - 1
      endif
c
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
# ifdef restnosp
c
c-----------------------------------------------------------------------
c     when restarting a prognostic surface pressure run from a  
c     stream function restart...
c     initialize ubar and vbar from the stream function case
c     first read psi at "tau" then psi at "tau+1"
c     warning: ubar and vbar will be multiplied by "hr" below when
c     "hr" is available
c-----------------------------------------------------------------------
c
      do m=2,1,-1
        read (iorest) stamp, iotext
        read (iorest) iimt, ijmt, ptd
        if (iimt .ne. imt .or. ijmt .ne. jmt) then
          write (stdout,*) ' => Error restio #0,imt=',iimt,',jmt=',ijmt
     &,   ' m=',m
        endif
        do jrow=1,jmtm1
          do i=2,imtm1
            diag1   = ptd(i+1,jrow+1) - ptd(i  ,jrow)
            diag0   = ptd(i  ,jrow+1) - ptd(i+1,jrow)
	    if (m .eq. 1) then
              ubar(i,jrow,1) = -(diag1+diag0)*dyu2r(jrow)
              ubar(i,jrow,2) =  (diag1-diag0)*dxu2r(i)*csur(jrow)
	    else
              ubarm1(i,jrow,1) = -(diag1+diag0)*dyu2r(jrow)
              ubarm1(i,jrow,2) =  (diag1-diag0)*dxu2r(i)*csur(jrow)
	    endif
          enddo
        enddo
      enddo
      call setbcx (ubar(1,1,1), imt, jmt)
      call setbcx (ubarm1(1,1,1), imt, jmt)
      call setbcx (ubar(1,1,2), imt, jmt)
      call setbcx (ubarm1(1,1,2), imt, jmt)
      read (iorest) stamp, iotext
      read (iorest)
      read (iorest) stamp, iotext
      read (iorest)
c
c     set the surface pressure equal to zero
c
      do it=1,2
        do jrow=1,jmt
          do i=1,imt
	    ps(i,jrow,it) = c0
          enddo
        enddo
      enddo
# else
c
      read (iorest) stamp, iotext
      read (iorest) iimt, ijmt, ps
      if (iimt .ne. imt .or. ijmt .ne. jmt) then
          write (stdout,*) ' => Error restio #A,imt=',iimt,',jmt=',ijmt
      endif
c
c     first read "tau" then "tau+1"
c
      read (iorest) stamp, iotext
      read (iorest) iimt, ijmt, ((ubarm1(i,jrow,1),i=1,imt),jrow=1,jmt)
     &,                         ((ubarm1(i,jrow,2),i=1,imt),jrow=1,jmt)
      if (iimt .ne. imt .or. ijmt .ne. jmt) then
          write (stdout,*) ' => Error restio #B,imt=',iimt,',jmt=',ijmt
      endif
c
      read (iorest) stamp, iotext
      read (iorest) iimt, ijmt, ((ubar(i,jrow,1),i=1,imt),jrow=1,jmt)
     &,                         ((ubar(i,jrow,2),i=1,imt),jrow=1,jmt)
      if (iimt .ne. imt .or. ijmt .ne. jmt) then
          write (stdout,*) ' => Error restio #C,imt=',iimt,',jmt=',ijmt
      endif
# endif
c
c     initialize the kflds disk which is used for surface pressure. 
c
      do jrow=1,jmt
        do i=1,imt
          pguess(i,jrow) = c0
          ptd(i,jrow)    = ps(i,jrow,1)
        enddo
      enddo
      call oput (kflds, nwds, nkflds-1, ptd)
      call oput (kflds, nwds, nkflds, ptd)
c
#endif
c
#ifdef stream_function
c
c     first read "tau" then "tau+1" stream function
c
      read (iorest) stamp, iotext
      read (iorest) iimt, ijmt, ((psi(i,jrow,2),i=1,imt),jrow=1,jmt)
      if (iimt .ne. imt .or. ijmt .ne. jmt) then
        write (stdout,*) ' => Error restart #1, imt=',iimt,', jmt=',ijmt
      endif
c
      read (iorest) stamp, iotext
      read (iorest) iimt, ijmt, ((psi(i,jrow,1),i=1,imt),jrow=1,jmt) 
      if (iimt .ne. imt .or. ijmt .ne. jmt) then
        write (stdout,*) ' => Error restart #2, imt=',iimt,', jmt=',ijmt
      endif
c
c     read the two guess fields
c
      read (iorest) stamp, iotext
      read (iorest) iimt, ijmt, ptd
      if (iimt .ne. imt .or. ijmt .ne. jmt) then
        write (stdout,*) ' => Error restart #3, imt=',iimt,', jmt=',ijmt
      endif
      call oput (kflds, nwds, nkflds-1, ptd)
c
      read (iorest) stamp, iotext
      read (iorest) iimt, ijmt, ptd
      if (iimt .ne. imt .or. ijmt .ne. jmt) then
        write (stdout,*) ' => Error restart #4, imt=',iimt,', jmt=',ijmt
      endif
      call oput (kflds, nwds, nkflds, ptd)
#endif
c
c     read "kmt" which indicates the number of t" cells to ocean floor
c
      read (iorest) stamp, iotext
      read (iorest) iimt, ijmt, kmt
      if (iimt .ne. imt .or. ijmt .ne. jmt) then
        write (stdout,*) ' => Error restart #5, imt=',iimt,', jmt=',ijmt
      endif
c
c-----------------------------------------------------------------------
c     construct depth arrays associated with "u" cells
c-----------------------------------------------------------------------
c
      call depth_u (kmt, imt, jmt, zw, km, kmu, h, hr)
c
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
# ifdef restnosp
c
c     complete external mode calculation using "hr"
c
      do m=1,2
        do jrow=1,jmt
	  do i=1,imt
	    ubar(i,jrow,m)   = ubar(i,jrow,m)*hr(i,jrow)
	    ubarm1(i,jrow,m) = ubarm1(i,jrow,m)*hr(i,jrow)
	  enddo
	enddo
	call setbcx (ubar(1,1,m), imt, jmt)
	call setbcx (ubarm1(1,1,m), imt, jmt)
      enddo        
# endif
#endif
c
c-----------------------------------------------------------------------
c     read the "tau" latitude rows
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
c
        if (wide_open_mw) then
	  j = jrow
	else
	  j = jmw
	endif
c
        read (iorest) stamp, iotext
        read (iorest) ilensl, bufsl
        if (ilensl .ne. nslab) then
          write (stdout,*) ' => Error rest #6, lensl=',ilensl
     &,                         ',jrow=',jrow
        endif
#ifdef tcvmix
c       How are the tcvmix variables read in???????
        uu = 1.0/undef
	uu = uu**2
#endif
c
c       velocity is internal mode component only
c       note that tracers are stored first, then velocities on restart
c       for compatability with previous versions
c
	do k=1,km
	  do i=1,imt
	    u(i,k,j,1,tau) = bufsl(i,k,nvar-1)
	    u(i,k,j,2,tau) = bufsl(i,k,nvar)
	    do n=1,nvar-2
	      t(i,k,j,n,tau) = bufsl(i,k,n)
	    enddo
	  enddo
	enddo
c
c       initialize every latitude
c
        if (wide_open_mw) then
c
c         do nothing since "tau" data is in place in the MW
c
	else
          call putrow (latdisk(taudisk), nslab, jrow, u(1,1,j,1,tau)
     &,                                               t(1,1,j,1,tau))
	endif
      enddo
c
      if (wide_open_mw) then
c
c       Initialze 1st and last latitude row for tau-1 to prevent
c       use of uninitialized values on boundary row.
c
        do j=1,jmt,jmt-1
	  do k=1,km
	    do i=1,imt
	      u(i,k,j,1,taum1) = u(i,k,j,1,tau)
	      u(i,k,j,2,taum1) = u(i,k,j,2,tau)
	      do n=1,nvar-2
	        t(i,k,j,n,taum1) = t(i,k,j,n,tau)
	      enddo
	    enddo
	  enddo
	enddo
      endif
c
c-----------------------------------------------------------------------
c     read the "tau+1" latitude rows
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
c
        if (wide_open_mw) then
	  j = jrow
	else
	  j = jmw
	endif
c
        read (iorest) stamp, iotext
        read (iorest) ilensl, bufsl
        if (ilensl .ne. nslab) then
          write (stdout,*) ' => Error restart #7, lensl=',ilensl
     &,                    ',jrow=',jrow
        endif
#ifdef tcvmix
c       How are the tcvmix variables read in???????
        uu = 1.0/undef
	uu = uu**2
#endif
c
c       velocity is internal mode component only
c       note that tracers are stored first, then velocities on restart
c       for compatability with previous versions
c
	do k=1,km
	  do i=1,imt
	    u(i,k,j,1,taup1) = bufsl(i,k,nvar-1)
	    u(i,k,j,2,taup1) = bufsl(i,k,nvar)
	    do n=1,nvar-2
	      t(i,k,j,n,taup1) = bufsl(i,k,n)
	    enddo
	  enddo
	enddo
c
c       initialize every latitude
c
        if (wide_open_mw) then
c
c         do nothing since "tau+1" data is in place in the MW
c
	else
          call putrow (latdisk(taup1disk), nslab, jrow
     &,                u(1,1,j,1,taup1), t(1,1,j,1,taup1))
	endif
      enddo
c
c
#ifdef trajectories
      read (iorest, err=210) stamp, iotext
      read (iorest)
      read (iorest) stamp, iotext
# ifdef lyapunov
      read (iorest,err=205) npart, pxyz, pijk, em
# else
      read (iorest) npart, pxyz, pijk
# endif
      write (stdout, '(a,i8,a)') 
     &'===> also read trajectory restart for ',npart,' particles'
      if (npart .ne. nptraj) then
        write (stdout,*) ' => Particle trajectory restart error:'
        write (stdout,*) '    read ',npart,' particles but expected '
     &,                     nptraj 
        stop "=>rdrest"
      endif
      go to 220
# ifdef lyapunov
205   continue
      write (stdout,'(/a/a/)')
     & '==>Note: Did not find lyapunov part of trajectory restart...'
     &,'         Setting all particles to Initial Positions'   
# endif
210   continue
c
c     if no restart data is available. initialize the particles
c
      call ptraji
220   continue
#endif
      rewind iorest
c
      call relunit (iorest)
      call toc ('ocean', 'reading restart')
c
      return
      end



      subroutine write_restart
c
c=======================================================================
c
c     read the restart file
c     
c     author:        a. rosati             e-mail=> ar@gfdl.gov
c                    r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
#include "param.h"
#ifdef trajectories
# include "ptraj.h"
#endif
#include "coord.h"
#include "emode.h"
#include "grdvar.h"
#include "iounit.h"
#include "levind.h"
#include "mw.h"
#include "tmngr.h"
      character*32 ostamp
      dimension bufsl(imt,km,nvar)
#ifdef stream_function
      dimension ext(imt,2)
#endif
c
      call getunit (iorest, 'restart.dta','u s r')
c
      call tic ('ocean', 'writing restart')
c
      write (stdout,'(a,i8,1x,a)')
     & '==> Writing MOM restart to file restart.dta at ts=', itt, stamp
c
c-----------------------------------------------------------------------
c     write restart file using unformatted fortran i/o
c
c     record 1 on the restart file is the control block (containing
c              the time step and relative year)
c     record`s 2 & 3 are for stream function time levels
c     record`s 4 & 5 are for the d(stream function)/dt fields used
c              to construct a guess for the elliptic solver
c     record 6 is for "kmt" indicating the number of ocean "t" cells
c              from the surface to the ocean floor.
c     the next "jmt" records contain prognostic variables at "tau"
c     the last "jmt" records contain prognostic variables at "tau+1"
c-----------------------------------------------------------------------
c
c     note that timestep "itt" corresponds to "tau+1"
c
      iotext = ' read (iorest) itt, irstdy, msrsdy, imt, jmt, km'
      write (iorest) stamp, iotext, expnam
      write (iorest) itt, iday(imodeltime), msday(imodeltime)
     &, imt, jmt, km
c
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
c
      iotext = 
     &' read (iorest) imt,jmt, ((ps(i,j,1),i=1,imt),j=1,jmt), ps(,,2)'
      write (iorest) stamp, iotext, expnam
      write (iorest) imt, jmt, ps
c
c     (ubarm1,ubar) is at (tau,tau+1) here
c
      iotext = 
     &' read (iorest) imt,jmt,((ubarm1(i,j),i=1,imt),j=1,jmt), vbarm1()'
      write (iorest) stamp, iotext, expnam
      write (iorest) imt, jmt, ((ubarm1(i,jrow,1),i=1,imt),jrow=1,jmt)
     &,                        ((ubarm1(i,jrow,2),i=1,imt),jrow=1,jmt)
c
      iotext = 
     &' read (iorest) imt, jmt, ((ubar(i,j),i=1,imt),j=1,jmt), vbar()'
      write (iorest) stamp, iotext, expnam
      write (iorest) imt, jmt, ((ubar(i,jrow,1),i=1,imt),jrow=1,jmt)
     &,                        ((ubar(i,jrow,2),i=1,imt),jrow=1,jmt)
#endif
c
#ifdef stream_function
c
c     first do psi at "tau" then at "tau+1"
c
      iotext = ' read (iorest) imt, jmt, ((psi(i,j,2),i=1,imt),j=1,jmt)'
      write (iorest) stamp, iotext, expnam
      write (iorest) imt, jmt, ((psi(i,jrow,2),i=1,imt),jrow=1,jmt)
c
      iotext = ' read (iorest) imt, jmt, ((psi(i,j,1),i=1,imt),j=1,jmt)'
      write (iorest) stamp, iotext, expnam
      write (iorest) imt, jmt, ((psi(i,jrow,1),i=1,imt),jrow=1,jmt)
c
c     guess fields
c
      call oget (kflds, nwds, nkflds-1, ptd)
      iotext = ' read (iorest) imt, jmt, ((g1(i,j),i=1,imt),j=1,jmt)'
      write (iorest) stamp, iotext, expnam
      write (iorest) imt, jmt, ptd
c
      call oget (kflds, nwds, nkflds, ptd)
      iotext = ' read (iorest) imt, jmt, ((g2(i,j),i=1,imt),j=1,jmt)'
      write (iorest) stamp, iotext, expnam
      write (iorest) imt, jmt, ptd
#endif
c
      iotext = ' read (iorest) imt, jmt, ((kmt(i,j),i=1,imt),j=1,jmt)'
      write (iorest) stamp, iotext, expnam
      write (iorest) imt, jmt, kmt
c
c-----------------------------------------------------------------------
c     save the "tau" latitude rows
c-----------------------------------------------------------------------
c
      ilensl = imt*km*nvar
      do jrow=1,jmt
c
        if (wide_open_mw) then
	  j = jrow
c
c         remove external mode from "tau". since psi has been updated
c         psi(,,2) is at "tau" 
c
          if (jrow .lt. jmt) then
#ifdef stream_function
	    do i=2,imt-1
              diag1    = psi(i+1,jrow+1,2) - psi(i  ,jrow,2)
              diag0    = psi(i  ,jrow+1,2) - psi(i+1,jrow,2)
	      ext(i,1)  = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
	      ext(i,2)  =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)
     &                       *csur(jrow)
            enddo
            do n=1,2
	      do k=1,km
                do i=2,imt-1
		  if (k .le. kmu(i,jrow)) then
                    u(i,k,j,n,tau) = (u(i,k,j,n,tau) - ext(i,n))
		  endif
                enddo
	      enddo
	    enddo
#endif
#if defined rigid_lid_surface_pressure || defined implicit_free_surface
c
c           since external mode has been updated, ubarm1 is at "tau"
c
            do n=1,2
	      do k=1,km
                do i=2,imt-1
		  if (k .le. kmu(i,jrow)) then
                    u(i,k,j,n,tau) = (u(i,k,j,n,tau)-ubarm1(i,jrow,n))
		  endif
                enddo
              enddo
            enddo
#endif
            do n=1,2
	      call setbcx (u(1,1,j,n,tau), imt, km)
	    enddo
	  endif
c
	else
	  j = jmw
          call getrow (latdisk(taudisk), nslab, jrow, u(1,1,j,1,tau)
     &,                                               t(1,1,j,1,tau))
	endif
c
	do k=1,km
	  do i=1,imt
	    bufsl(i,k,nvar-1) = u(i,k,j,1,tau)
	    bufsl(i,k,nvar) = u(i,k,j,2,tau)
	    do n=1,nvar-2
	      bufsl(i,k,n) = t(i,k,j,n,tau)
	    enddo
	  enddo
	enddo
        write (iotext,'(a15,i4)') ' taudisk jrow =',jrow
        iotext(20:) = ' read(iorest) lensl, (bufsl(i),i=1,lensl)'
        write (iorest) stamp, iotext, expnam
        write (iorest) ilensl, bufsl
      enddo
c
c-----------------------------------------------------------------------
c     save the "tau+1" latitude rows
c-----------------------------------------------------------------------
c
      do jrow=1,jmt
c
        if (wide_open_mw) then
	  j = jrow
c
c         do nothing since data is already in "tau+1" position in MW
c         and velocity is already internal mode only
c
	else
	  j = jmw
          call getrow (latdisk(taup1disk), nslab, jrow
     &,                u(1,1,j,1,taup1), t(1,1,j,1,taup1))
	endif
c
	do k=1,km
	  do i=1,imt
	    bufsl(i,k,nvar-1) = u(i,k,j,1,taup1)
	    bufsl(i,k,nvar) = u(i,k,j,2,taup1)
	    do n=1,nvar-2
	      bufsl(i,k,n) = t(i,k,j,n,taup1)
	    enddo
	  enddo
	enddo
        write (iotext,'(a17,i4)') ' taup1disk jrow =',jrow
        iotext(22:) = ' read(iorest) lensl, (bufsl(i),i=1,lensl)'
        write (iorest) stamp, iotext, expnam
        write (iorest) ilensl, bufsl
      enddo
c
#ifdef trajectories
      reltim = relyr
      iotext = 'read (iorest) reltim'
      write (iorest) stamp, iotext, expnam
      write (iorest) reltim
# ifdef lyapunov
      iotext ='read (iorest) nptraj, pxyz, pijk, em'
      write (iorest) stamp, iotext, expnam
      write (iorest) nptraj, pxyz, pijk, em
# else
      iotext ='read (iorest) nptraj, pxyz, pijk'
      write (iorest) stamp, iotext, expnam
      write (iorest) nptraj, pxyz, pijk
# endif
      write (stdout, '(a,i8,a)')
     &'===> also writing trajectory restart for ',nptraj, 'particles'
#endif
c
c-----------------------------------------------------------------------
c     write time step and stamp to file archive.time for purpose of
c     verifying that the restart file "restart.dta" is the correct one
c     when beginning the next job. The run script should save this file
c     and restore it when restarting. (security measure)
c-----------------------------------------------------------------------
c     
      call getunit (ioarch, 'archive.time','f s a')
      write (ioarch,8888) itt, stamp
      call relunit (ioarch)
c
      call relunit (iorest)
      call toc ('ocean', 'writing restart')
c
      return
8888  format (1x,i10,1x,a32)
      end
