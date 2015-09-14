#ifdef held_larichev
      subroutine hlmixi (error, cifdef, ifdmax, nifdef)
c
      logical error
      character*(*) cifdef(ifdmax)
c
# include "size.h"
# include "accel.h"
# include "coord.h"
# include "hmixc.h"
# include "iounit.h"
# include "scalar.h"
# include "stdunits.h"
c
c=======================================================================
c     Initialization for ...
c     A Scaling Theory for Horizontally Homogeneous, Baroclinically
c     Unstable Flow on a Beta Plane (Submitted to JAS, March 1995). 
c    
c     input:
c       error  = logical to signal problems
c       cifdef = array of character strings for listing enabled "ifdefs"
c       ifdmax = size of "cifdef"
c       nifdef = current number of enabled "ifdefs" 
c
c     output:
c       hl_depth = integration depth
c       hl_back  = minimum diffusion coefficient
c       hl_max   = maximum diffusion coefficient
c       
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
c
      namelist /hlmix/ hl_depth, hl_back, hl_max
c
      write (stdout,'(/,20x,a,/)')
     & 'H L M I X    I N I T I A L I Z A T I O N'
c
c-----------------------------------------------------------------------
c     initialize variables (all mixing units are cm**2/sec.)
c-----------------------------------------------------------------------
c
      hl_depth = 500.0e2
      hl_back  = 1.e4
      hl_max   = 1.e9      
      c0       = 0.0
c
c-----------------------------------------------------------------------
c     provide for namelist over-ride of above settings + documentation
c-----------------------------------------------------------------------
c
      call getunit (io, 'namelist', 'fsr')
      read (io,hlmix,end=100)
100   continue
c
c-----------------------------------------------------------------------
c     set no-flux condition on density difference across bottom level
c     initialize inverse richardson number = zero at all levels
c-----------------------------------------------------------------------
c
      do j=1,jmw
        do i=1,imt
	  droz(i,km,j) = c0
        enddo
      enddo
      do j=1,jemw
        do k=1,km
          do i=1,imt
	    rich_inv(i,k,j) = c0
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     add character string to "ifdef option list" indicating that this
c     option is enabled
c-----------------------------------------------------------------------
c
      nifdef = nifdef + 1
      cifdef(nifdef) = 'held_larichev'
c
c-----------------------------------------------------------------------
c     check for problems
c-----------------------------------------------------------------------
c
# if !defined isopycmix
        write (stdout,'(/,(1x,a))')
     & '==> Error:"isopycmix" must be enabled since "held_larichev"'
     &,'            predicts mixing coeffs used by "isopycmix"         '        
        error = .true.
# endif
c
c     write out namelist values
c
      write (stdout,hlmix)
      call relunit (io)
      call getunit (iodoc, 'document.dta', 'f s a')
      write (iodoc, hlmix)
      call relunit (iodoc)
      return
      end


      subroutine hlmix (joff, js, je, is, ie)
c
c=======================================================================
c     Compute mixing coefficients based on...
c     A Scaling Theory for Horizontally Homogeneous, Baroclinically
c     Unstable Flow on a Beta Plane (Submitted to JAS, March 1995). 
c
c
c     inputs:
c
c      joff   = offset between rows in the MW and latitude rows
c               "joff" > 0 moves variables
c      js     = starting row for loading variables to calculate
c               coefficients. calculations start at jstrt=max(js-1,jsmw)
c      je     = ending row for loading variables to calculate
c               coefficients. calculations end at je-1
c      is     = starting index for calculating coefficients in the 
c               longitude direction
c      ie     = ending index for calculating coefficients in the 
c               longitude direction
c      grav   = gravity (cm/sec**2)
c      umask  = land/sea mask on "u" grid (land=0.0, sea=1.0)
c      tmask  = land/sea mask on "t" grid (land=0.0, sea=1.0)
c
c     outputs:
c
c      hl_coeff = held_larichev diffusion coefficient
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c=======================================================================
c
      integer tlev
# include "size.h"
# include "coord.h"
# include "grdvar.h"
# include "hmixc.h"
# include "mw.h"
# include "scalar.h"
# include "switch.h"
      dimension ro(imt,km,1:jmw)
      dimension sumz(imt,1:jemw)
c#define debug_hlmix
# ifdef debug_hlmix
      common /debughl/ hlcoeff(imt,jmt)
      if (joff .eq. 0) then
        do j=1,jmt
	  do i=1,imt
	    hlcoeff(i,j) = 0.0
	  enddo
	enddo
      endif
# endif
c
# ifdef timing
      call tic ('hmixc', 'hlmix')
# endif
c
c-----------------------------------------------------------------------
c     set local constants 
c-----------------------------------------------------------------------
c
      p25   = 0.25
      c0    = 0.0
      p5    = 0.5
      epsln = 1.e-20
      fx    = -p25/grav
      istrt = max(2,is)
      iend  = min(imt-1,ie)
      tlev  = taum1
      kstart= 1
      kend  = max (1,min(indp(hl_depth, zt, km),km-1))
c
c-----------------------------------------------------------------------
c     set "ro" (density) at j=1 for 1st memory window otherwise ... move 
c     variables from top two rows to bottom two rows to eliminate
c     redundant calculation
c-----------------------------------------------------------------------
c
      if (joff .eq. 0) then
        do k=1,km
	  do i=istrt-1,iend+1
	    ro(i,k,1) = c0
	  enddo
	enddo
      else
        call movehlmix (istrt-1, iend+1)
      endif
c
c-----------------------------------------------------------------------
c     compute density difference across bottom of "t" cells at tau-1
c     for rows js through je in the MW. Set density difference = zero
c     across bottom and in land areas 
c-----------------------------------------------------------------------
c
      do ks=1,2
        call statec (t(1,1,1,1,tlev), t(1,1,1,2,tlev), ro(1,1,jsmw)
     &,              max(js,jsmw), je, istrt-1, iend+1, ks)
        do j=js,je
          do k=ks,km-1,2
            do i=istrt-1,iend+1
              droz(i,k,j) = (ro(i,k,j) - ro(i,k+1,j))*tmask(i,k+1,j)
            enddo
          enddo
	enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute inverse richardson numbers on bottom of "u" cells
c-----------------------------------------------------------------------
c
      jsriu = max(js,jsmw)-1
      do j=jsriu,je-1
        do k=kstart,kend
          t1 = fx/dzw(k)
          do i=istrt,iend
            rich_inv(i,k,j) = t1*umask(i,k+1,j)*(
     &                    (u(i,k,j,1,tlev) - u(i,k+1,j,1,tlev))**2 +
     &                    (u(i,k,j,2,tlev) - u(i,k+1,j,2,tlev))**2) / 
     &                    (droz(i,k,j+1) + droz(i+1,k,j+1) +
     &                     droz(i,k,j)   + droz(i+1,k,j) + epsln)
          enddo
        enddo
	call setbcx (rich_inv(1,1,j), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     integrate inverse Richardson number vertically over "u" cells
c     assume Ri at ocean surface is the same as at bottom of first cell
c-----------------------------------------------------------------------
c
      do j=jsriu,je-1
        k = kstart
        do i=istrt,iend
          hl_u(i,j) = rich_inv(i,k,j)*dzw(k-1)
          sumz(i,j) = dzw(k-1) + epsln
        enddo
        do k=kstart,kend
          do i=istrt,iend
            hl_u(i,j) = hl_u(i,j) + rich_inv(i,k,j)*dzw(k)
            sumz(i,j) = sumz(i,j) + dzw(k)*umask(i,k+1,j)
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     constrain integrated inverse Richardson numbers to be >= epsln
c-----------------------------------------------------------------------
c
      do j=jsriu,je-1
        do i=istrt,iend
	  hl_u(i,j) = max(hl_u(i,j),epsln)
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     construct Diffusion coefficient = 1/(beta**2*T**3) on "u" cells 
c-----------------------------------------------------------------------
c
      do j=jsriu,je-1
        jrow = joff + j
        fsq = cori(jrow,1)**2
	beta_sq_r = 1.0/(2.0*omega*csu(jrow)/radius)**2
        do i=istrt,iend
          tm2 = fsq*hl_u(i,j)/sumz(i,j)
          hl_u(i,j) = beta_sq_r*tm2**1.5
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     Limit Diffusion coefficient to be between hl_max and hl_back
c-----------------------------------------------------------------------
c
      do j=jsriu,je-1
        do i=istrt,iend
	  hl_u(i,j) = max(min(hl_u(i,j),hl_max),hl_back)
        enddo
# ifdef cyclic
	hl_u(1,j)   = hl_u(imt-1,j)
	hl_u(imt,j) = hl_u(2,j)
# else
	hl_u(1,j)   = c0
	hl_u(imt,j) = c0
# endif
      enddo
c
c-----------------------------------------------------------------------
c     construct Diffusion coefficient on eastern face of "T" cells
c-----------------------------------------------------------------------
c
      do j=jsmw,je-1
        do i=istrt,iend
	  hl_e(i,j) = p5*(hl_u(i,j) + hl_u(i,j-1))
        enddo
# ifdef cyclic
        hl_e(1,j)   = hl_e(imt-1,j)
	hl_e(imt,j) = hl_e(2,j)
# else
        hl_e(1,j)   = c0
	hl_e(imt,j) = c0
# endif
      enddo
c
c-----------------------------------------------------------------------
c     construct Diffusion coefficient centered in "T" cells
c-----------------------------------------------------------------------
c
      do j=jsmw,je-1
        do i=istrt,iend
	  hl_b(i,j) = p25*(hl_u(i,j) + hl_u(i,j-1) + hl_u(i-1,j)
     &                   + hl_u(i-1,j-1))
        enddo
# ifdef cyclic
        hl_b(1,j)   = hl_b(imt-1,j)
	hl_b(imt,j) = hl_b(2,j)
# else
        hl_b(1,j)   = c0
	hl_b(imt,j) = c0
# endif
      enddo
c
c-----------------------------------------------------------------------
c     construct Diffusion coefficient on northern face of "T" cells
c-----------------------------------------------------------------------
c
      do j=1,je-1
        do i=istrt,iend
	  hl_n(i,j) = p5*(hl_u(i,j) + hl_u(i-1,j))
        enddo
# ifdef cyclic
        hl_n(1,j)   = hl_n(imt-1,j)
	hl_n(imt,j) = hl_n(2,j)
# else
        hl_n(1,j)   = c0
	hl_n(imt,j) = c0
# endif
      enddo
c
#ifdef debug_hlmix
      do j=jsriu,je-1
        do i=1,imt
	  hlcoeff(i,jrow) = hl_u(i,j)
	enddo
	jrow = j + joff
	if (jrow .ge. 12 .and. jrow .le. 20 .and. prxzts .and. eots)then
	  print *,' tp: jrow=',jrow
	  call scope (rich_inv(1,1,j), imt, imt, km, 'rich_inv')
	endif
	if (jrow .eq. jmt-1 .and. prxzts .and. eots) then
	  print *,'hlcoeff='
	  call matrix (hlcoeff, imt, 1, imt, -1, -jmt, 0.0)
	endif
      enddo
#endif
# ifdef timing
      call toc ('hmixc', 'hlmix')
# endif
      return
      end


      
      subroutine movehlmix (is, ie)
c
c=======================================================================
c     as the MW moves northward, move data from the last two rows
c     into the first two rows. (last 3 rows if using biharmonic option)
c     
c     author:  r.c.pacanowski   e-mail  rcp@gfdl.gov
c=======================================================================
c
# include "param.h"
# include "hmixc.h"
c
      nrows = jmw - ncrows
      do move=1,nrows
        jfrom = jmw - (nrows - move)
	jto   = move
c
c-----------------------------------------------------------------------
c       move quantities with rows dimensioned (1:jmw)
c-----------------------------------------------------------------------
c
        do k=1,km
	  do i=is,ie
	    droz(i,k,jto)  = droz(i,k,jfrom)
          enddo
        enddo
c
c-----------------------------------------------------------------------
c       move quantities with rows dimensioned (1:jemw)
c-----------------------------------------------------------------------
c
        if (jfrom .le. jemw) then
          do k=1,km
	    do i=is,ie
	      rich_inv(i,k,jto)  = rich_inv(i,k,jfrom)
            enddo
          enddo
	  do i=is,ie
	    hl_u(i,jto) = hl_u(i,jfrom)
	  enddo
        endif
c
c-----------------------------------------------------------------------
c       move quantities with rows dimensioned (jsmw:jmw)
c-----------------------------------------------------------------------
c
        if (jto .ge. jsmw) then
c
c         nothing to move
c
        endif
c
      enddo
      return
      end
#else
      subroutine hlmix
      return
      end
#endif
