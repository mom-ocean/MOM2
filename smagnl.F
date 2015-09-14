#ifdef smagnlmix
      subroutine smagnli (error, cifdef, ifdmax, nifdef, hmixset)
c
c=======================================================================
c     Initialization for the Smagorinsky nonlinear horizontal viscosity 
c     as described in Rosati & Miyakoda (jpo,vol 18,#11,1988)
c     see Smagorinsky 1963, Mon Wea Rev, 91, 99-164. 
c     Also see Deardorff 1973 J. Fluid Eng. Sep., 429-438.
c     
c     input:
c
c       error  = logical to signal problems
c       cifdef = array of character strings for listing enabled "ifdefs"
c       ifdmax = size of "cifdef"
c       nifdef = current number of enabled "ifdefs" 
c       hmixset= logical to determine if a vertical mixing scheme was
c                chosen
c     output:
c
c      author: R. C. Pacanowski (implemention for MOM 2)   rcp@gfdl.gov 
c              A. Rosati  (MOM 1 version)     e-mail: ar@gfdl.gov
c=======================================================================
c
      logical hmixset, error
      character*(*) cifdef(ifdmax)
c
# include "size.h"
# include "hmixc.h"
# include "iounit.h"
# include "switch.h"
# include "stdunits.h"
      namelist /smagnl/ diff_c_back
c
      write (stdout,'(/,20x,a,/)')
     & 'S M A G M I X    I N I T I A L I Z A T I O N'
c
c-----------------------------------------------------------------------
c     initialize variables (all mixing units are cm**2/sec.)
c-----------------------------------------------------------------------
c
      c0          = 0.0
      diff_c_back = c0
c
c-----------------------------------------------------------------------
c     provide for namelist over-ride of above settings + documentation
c-----------------------------------------------------------------------
c
      call getunit (io, 'namelist', 'fsr')
      read (io,smagnl,end=100)
100   continue
      write (stdout,smagnl)
      call relunit (io)
      call getunit (iodoc, 'document.dta', 'f s a')
      write (iodoc, smagnl)
      call relunit (iodoc)
c
c-----------------------------------------------------------------------
c     add character string to "ifdef option list" indicating that this
c     option is enabled
c-----------------------------------------------------------------------
c
      nifdef = nifdef + 1
      cifdef(nifdef) = 'smagnlmix      '
c
c-----------------------------------------------------------------------
c     check for problems
c-----------------------------------------------------------------------
c
      if (hmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "smagnlmix"  cannot be enabled because another   '
     &,'            horizontal mixing scheme has been enabled    '        
        error = .true.
      else
        hmixset = .true.
      endif
c
c-----------------------------------------------------------------------
c     initialize arrays
c-----------------------------------------------------------------------
c
      do j=1,jemw
        do i=1,imt
	  do k=1,km
	    strain(i,k,j,1) = c0
	    strain(i,k,j,2) = c0
	    am_lambda(i,k,j)= c0 
	    am_phi(i,k,j)   = c0 
	    visc_cnu(i,k,j) = c0
	    diff_cnt(i,k,j) = c0 
          enddo
        enddo
      enddo
      do j=jsmw,jemw
        do i=1,imt
	  do k=1,km
	    visc_ceu(i,k,j) = c0
	    diff_cet(i,k,j) = c0 
          enddo
        enddo
      enddo
      return
      end


      subroutine smagnlc (joff, js, je, is, ie)
c
# include "size.h"
# include "grdvar.h"
# include "hmixc.h"
# include "mw.h"
# include "switch.h"
# include "stdunits.h"
      dimension clam(imt)
c
c=======================================================================
c     compute tension and shearing rates of strain, total deformation
c     and diffusion coefficients for the Smagorinsky nonlinear 
c     viscosity as described in Rosati & Miyakoda (jpo,vol 18,#11,1988).
c     
c     input:
c
c      joff   = offset between row j in the MW and latitude jrow on disk
c               "joff" > 0 moves variables
c      js     = starting row for calculations
c      je     = ending row for calculations
c      is     = starting longitude index for calculations
c      ie     = ending longitude index for calculations
c
c     output:
c
c      strain  = tension (1) and shearing (2) rates of strain
c      defor  = total deformation
c
c=======================================================================
c
# ifdef timing
      call tic ('hmixc', 'non-linear')
# endif
c
c-----------------------------------------------------------------------
c     set local constants
c-----------------------------------------------------------------------
c
      p25    = 0.25
      p5     = 0.5
      c0     = 0.0
      c1     = 1.0
      c2     = 2.0
      sqrt2r = c1/sqrt(c2)
      c14    = 0.14
c
      istrt  = max(2,is)
      iend   = min(imt-1,ie)
      jstrt  = max(1,js-1)
      jend   = je-1
c
c-----------------------------------------------------------------------
c     move variables from top two rows to bottom two rows to eliminate
c     redundant computation after the first MW
c-----------------------------------------------------------------------
c
      if (joff .ne. 0) then
        call movesmagnl (istrt-1, iend+1)
      endif
c
c-----------------------------------------------------------------------
c     compute rates of strain: tension and shear on north face of 
c     "u" cell. Rosati & Miyakoda (jpo,vol 18,#11,1988) Eqn 2.24 & 2.25 
c-----------------------------------------------------------------------
c
      do j=jstrt,jend
        jrow = j + joff
        f1   = p5*cstr(jrow+1)
        f2   = dytr(jrow+1)*cst(jrow+1)
        do k=1,km
          do i=istrt,iend
            strain(i,k,j,1) = ((u(i+1,k,j,1,taum1)+u(i+1,k,j+1,1,taum1))
     &                      - (u(i-1,k,j,1,taum1)+u(i-1,k,j+1,1,taum1)))
     &                          *f1*dxu2r(i)
     &                      - (csur(jrow+1)*u(i,k,j+1,2,taum1) - 
     &                         csur(jrow  )*u(i,k,j,2,taum1))*f2
            strain(i,k,j,2) = ((u(i+1,k,j,2,taum1)+u(i+1,k,j+1,2,taum1))
# ifdef free_slip
     &                          *umask(i+1,k,j)*umask(i+1,k,j+1)
# endif
     &                      - (u(i-1,k,j,2,taum1)+u(i-1,k,j+1,2,taum1)))
# ifdef free_slip
     &                          *umask(i-1,k,j)*umask(i-1,k,j+1)
# endif
     &                          *f1*dxu2r(i)
     &                      + (csur(jrow+1)*u(i,k,j+1,1,taum1) - 
     &                         csur(jrow  )*u(i,k,j,1,taum1))*f2
# ifdef free_slip
     &                          *umask(i,k,j)*umask(i,k,j+1)
# endif
          enddo
        enddo
        call setbcx (strain(1,1,j,1), imt, km)
        call setbcx (strain(1,1,j,2), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     compute effective anisentropic wavenumber of diffusing turbulence
c     for effective 2d isentropic wavenumber set coeffs
c     to: (c14*sqrt(csu(jrow)*dxu(i)*dyu(jrow)))**2*sqrt2r
c     Rosati & Miyakoda(jpo,vol 18,#11,1988). Eqn: 2.28 & 2.29
c     (note: Eqn 2.28 should have m**-1)
c     compute total deformation + viscosity coefficients on north face 
c     of "u" cell
c-----------------------------------------------------------------------
c
      do j=jstrt,jend
        jrow     = joff + j
        cphi  = (c14*dyu(jrow))**2*sqrt2r
        do i=istrt,iend
          clam(i) = (c14*csu(jrow)*dxu(i))**2*sqrt2r 
        enddo
        do k=1,km
          do i=istrt,iend
            deform     = sqrt(c2*(strain(i,k,j,1)**2 
     &                                + strain(i,k,j,2)**2))
            am_lambda(i,k,j) = clam(i)*deform
            am_phi(i,k,j)    = cphi*deform
          enddo
        enddo
        call setbcx (am_lambda(1,1,j), imt, km)
        call setbcx (am_phi(1,1,j), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     set j index for calculating diffusive coefficient on north face
c     of cells. 
c-----------------------------------------------------------------------
c
      jsc = max(jsmw,js-1)
c
# ifdef matrix_sections
      if (prxzts .and. eots) then
c
c-----------------------------------------------------------------------
c       calculate the mixing coefficients for momentum and diffusion
c       Rosati & Miyakoda(jpo,vol 18,#11,1988)
c       Eqn: 2.26 & 2.27  mixing coeff for momentum
c       Both "visc_cnu" and "visc_ceu" are purely diagnostic here since
c       momentum flux "diff_fe" and "diff_fn" is calculated in terms of
c       shearing and tension rates of strain in subroutine "smagnlm"
c-----------------------------------------------------------------------
c
c       viscosity coeff on north face of "u" cells
c
        do j=jstrt,jend
          do k=1,km
            do i=istrt,iend
              visc_cnu(i,k,j) = am_phi(i,k,j)
            enddo
          enddo
	  call setbcx (visc_cnu(1,1,j), imt, km)
        enddo
c
c       viscosity coeff on east face of "u" cells
c 
        do j=jsc,jend
          do k=1,km
            do i=istrt-1,iend
              visc_ceu(i,k,j) = p25*(am_lambda(i,k,j)
     &                             + am_lambda(i+1,k,j)
     &                             + am_lambda(i,k,j-1)
     &                             + am_lambda(i+1,k,j-1))
	    enddo
	  enddo
	  call setbcx (visc_ceu(1,1,j), imt, km)
        enddo
      endif
# endif
c
c-----------------------------------------------------------------------
c     calculate the diffusion coefficients for "t" cells
c     Rosati & Miyakoda(jpo,vol 18,#11,1988)
c     Eqn: 2.35 & 2.36  mixing coeff for tracers
c-----------------------------------------------------------------------
c
c     diffusion coeff on east face of "t" cells
c 
      do j=jsc,jend
        do k=1,km
          do i=istrt-1,iend
            diff_cet(i,k,j) = am_lambda(i,k,j-1) + diff_c_back
	  enddo
	enddo
	call setbcx (diff_cet(1,1,j), imt, km)
      enddo
c
c     diffusion coeff on north face of "t" cells. Index jsc is used
c     because "j-1" is needed. The southern wall is taken care of by
c     masking in the diffusion operator.
c
      do j=jsc,jend
        jrow = j + joff
        do k=1,km
          do i=istrt,iend
	    diff_cnt(i,k,j) = p25*(am_phi(i,k,j) + am_phi(i-1,k,j)
     &                           + am_phi(i,k,j-1) + am_phi(i-1,k,j-1))
     &                           + diff_c_back                    
          enddo
        enddo
        call setbcx (diff_cnt(1,1,j), imt, km)
      enddo
c
# ifdef matrix_sections
      if (prxzts .and. eots) then
        call diagnl (joff, jsc, jend)
      endif
# endif
# ifdef trace_indices
      write (stdout,'(4x,6(a,i4))')
     & "=> In smagnlc: js=",js," je=",je," joff=",joff," jstrt=",jstrt
     &," jend=",jend," jsc=",jsc
# endif
# ifdef timing
      call toc ('hmixc', 'non-linear')
# endif
      return
      end




      subroutine smagnlm (joff, js, je, is, ie, n)
c
# include "size.h"
# include "grdvar.h"
# include "hmixc.h"
# include "mw.h"
# include "scalar.h"
# include "stdunits.h"
c
c=======================================================================
c     compute diffusive flux across north and east face of "u" cells
c     for velocity component "n"
c
c     input:
c
c      joff   = offset between row j in the MW and latitude jrow
c               "joff" > 0 moves variables
c      js     = starting row for calculations
c      je     = ending row for calculations
c      is     = starting longitude index for calculations
c      ie     = ending longitude index for calculations
c      n      = velocity component: 1 = "u", 2 = "v"
c        + quantities from subroutine "dform"
c
c     output:
c
c      diff_fe  = diffusive flux across eastern face of "u" cells
c      diff_fn  = diffusive flux across northern face of "u" cells
c
c=======================================================================
c
# ifdef timing
      call tic ('hmixc', 'non-linear')
# endif
c
c-----------------------------------------------------------------------
c     set local constants
c-----------------------------------------------------------------------
c
      p25    = 0.25
      p5     = 0.5
      c0     = 0.0
c
      istrt  = max(2,is)
      iend   = min(imt-1,ie)
      jstrt  = max(2,js)
      jend   = min(jemw,je)
c
c-----------------------------------------------------------------------
c     compute zonal flux components of the stress tensor on the eastern
c     face of the "u" cell
c     Rosati & Miyakoda(jpo,vol 18,#11,1988). Eqn: 2.18
c-----------------------------------------------------------------------
c
      do j=jstrt,jend
        do k=1,km
          do i=istrt-1,ie
            diff_fe(i,k,j) = p25*(
     &                am_lambda(i,k,j)*strain(i,k,j,n)
     &              + am_lambda(i,k,j-1)*strain(i,k,j-1,n)
     &              + am_lambda(i+1,k,j)*strain(i+1,k,j,n)
     &              + am_lambda(i+1,k,j-1)*strain(i+1,k,j-1,n))
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute meridional flux components of the stress tensor on the 
c     north face of the "u" cell
c     Rosati & Miyakoda(jpo,vol 18,#11,1988). Eqn: 2.19
c-----------------------------------------------------------------------
c
      if (n .eq. 1) then
c
c       northward flux term for zonal memoentum eqn is zero
c
        do j=jstrt-1,jend
          do k=1,km
            do i=istrt,iend
	      diff_fn(i,k,j) = c0
            enddo
          enddo
        enddo
c
c       compute second term (which is not a flux term) when working on 
c       the zonal momentum eqn.
c
        do j=jstrt,jend
          jrow   = joff + j
          rmsq   = rho0r*csur(jrow)**2*dyur(jrow)
          cstsq  = cst(jrow+1)**2
          cstsqm = cst(jrow)**2
          do k=1,km
            do i=istrt,iend
	      smag_metric(i,k,j) = rmsq*
     &                     (am_phi(i,k,j)*strain(i,k,j,2)*cstsq-
     &                      am_phi(i,k,j-1)*strain(i,k,j-1,2)*cstsqm)
            enddo
          enddo
        enddo
      elseif (n .eq. 2) then
c
c       northward flux for meridional momentum equation.
c
        do j=jstrt-1,jend
	  jrow = joff + j
          do k=1,km
            do i=istrt,iend
	      diff_fn(i,k,j) =-cst(jrow+1)*am_phi(i,k,j)*strain(i,k,j,1)
            enddo
          enddo
        enddo
c
c       compute third term (which is not a flux term) when working on 
c       the meridional component of the momentum eqn.
c
        do j=jstrt,jend
          jrow = joff + j
          f1   = csur(jrow)*sine(jrow)*p5/radius
          do k=1,km
            do i=istrt,iend
	      smag_metric(i,k,j) = f1*(am_lambda(i,k,j)*strain(i,k,j,1)
     &                          +  am_lambda(i,k,j-1)*strain(i,k,j-1,1))
            enddo
          enddo
        enddo
      endif
c
# ifdef timing
      call toc ('hmixc', 'non-linear')
# endif
# ifdef trace_indices
      write (stdout,'(4x,5(a,i4))')
     & "=> In smagnlm: jstrt=",jstrt," jend=",jend," joff=",joff
     &," jrows=",jstrt+joff," to ",jend+joff 
# endif
      return
      end




# ifdef matrix_sections
      subroutine diagnl (joff, js, je)
#  include "param.h"
#  include "coord.h"
#  include "cprnts.h"
#  include "hmixc.h"
#  include "iounit.h"
#  include "switch.h"
#  include "tmngr.h"
c
c=======================================================================
c     Show some diagnostics
c
c     input:
c
c      joff   = offset between row j in the MW and latitude jrow on disk
c               "joff" > 0 moves variables
c      js     = starting row for calculations
c      je     = ending row for calculations
c=======================================================================
c
#  ifdef timing
      call tic ('diagnostic', 'matrix sections')
#  endif
      do j=js,je
	jrow = j + joff
        reltim = relyr
        do jlat=1,nlatpr
          jj = indp (prlat(jlat), yt, jmt)
          if (jj .eq. jrow .and. prlat(jlat) .le. yt(jmt)) then
            is = indp (prslon(jlat), xt, imt)
            ie = indp (prelon(jlat), xt, imt)
            ks = indp (prsdpt(jlat), zt, km)
            ke = indp (predpt(jlat), zt, km)
            fx = 1.0e-2
c
c           write out the diffusion coeffs for tracers
c
            scl = 1.e7
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'diff_cet', itt, jrow 
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (diff_cet(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
	      call getunit (io, 'sections.dta','u s a ieee')
              write (stdout,*) ' => diff_cet ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((diff_cet(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, diff_cet(1,1,j), imt*km)
              call relunit (io)
            endif
c
            scl = 1.e7
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'diff_cnt', itt, jrow 
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (diff_cnt(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
	      call getunit (io, 'sections.dta','u s a ieee')
              write (stdout,*) ' => diff_cnt ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((diff_cnt(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, diff_cnt(1,1,j), imt*km)
              call relunit (io)
            endif
c
c           write out the viscosity coeffs for momentum
c
            scl = 1.e7
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'visc_ceu', itt, jrow 
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (visc_ceu(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
	      call getunit (io, 'sections.dta','u s a ieee')
              write (stdout,*) ' => visc_ceu ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((visc_ceu(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, visc_ceu(1,1,j), imt*km)
              call relunit (io)
            endif
c
            scl = 1.e7
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'visc_cnu', itt, jrow 
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (visc_cnu(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
	      call getunit (io, 'sections.dta','u s a ieee')
              write (stdout,*) ' => visc_cnu ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((visc_cnu(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, visc_cnu(1,1,j), imt*km)
              call relunit (io)
            endif
          endif
        enddo
      enddo
#  ifdef timing
      call toc ('diagnostic', 'matrix sections')
#  endif
      return
9100    format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     &,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     &,', scaling=',1pg10.3)
      end
# endif



      

      subroutine movesmagnl (is, ie)
c
c=======================================================================
c     as the MW moves northward, move data from the last two rows
c     into the first two rows.
c     
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
c          nothing to move
c-----------------------------------------------------------------------
c       move quantities with rows dimensioned (1:jemw)
c-----------------------------------------------------------------------
c
        if (jfrom .le. jemw) then
          do k=1,km
	    do i=is,ie
	      strain(i,k,jto,1)  = strain(i,k,jfrom,1)
	      strain(i,k,jto,2)  = strain(i,k,jfrom,2)
	      am_lambda(i,k,jto) = am_lambda(i,k,jfrom)
	      am_phi(i,k,jto)    = am_phi(i,k,jfrom)
	      visc_cnu(i,k,jto)  = visc_cnu(i,k,jfrom)
	      diff_cnt(i,k,jto)  = diff_cnt(i,k,jfrom)
            enddo
          enddo
# ifdef trace_indices
          write (stdout,'(4x,2(a,i4))')
     & "=> In movesmagnl: moving variables on row ",jfrom," to row ",jto
# endif
        endif
c
c-----------------------------------------------------------------------
c       move quantities with rows dimensioned (jsmw:jmw)
c-----------------------------------------------------------------------
c
        if (jto .ge. jsmw) then
c         nothing to move
        endif
c
      enddo
      return
      end
#else
      subroutine smagnl
      return
      end
#endif

