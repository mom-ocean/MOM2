#if defined isopycmix || defined gent_mcwilliams
      subroutine isopi (error, cifdef, ifdmax, nifdef, am, ah)
c
c=======================================================================
c
c               Initialization for isopycnal mixing scheme
c
c      -Disopycmix gives Redi/Cox version using the mixing tensor from
c       Gent_Mcwilliams but no explicit use of advective velocity to
c       parameterize effect of eddies on isopycnals 
c
c      -Disopycmix -Dgent_mcwilliams gives adds the advective velocity
c       parameterization of Gent_McWilliams
c
c     input:
c       error  = logical to signal problems
c       cifdef = array of character strings for listing enabled "ifdefs"
c       ifdmax = size of "cifdef"
c       nifdef = current number of enabled "ifdefs" 
c       slmxr  = 1/(max slope of isopycnals)
c       ahisop = isopycnal tracer diffusivity(cm**2/sec) 
c       athkdf = isopycnal thickness diffusivity (cm**2/sec)
# ifdef isopycmixspatialvar
c       dciso1 = isopycnal tracer diffusivity coeffs modified based
c              on the slopes of the isopycnal surfaces on the east face
c              of "T" cells.
c       dciso2 = isopycnal tracer diffusivity coeffs modified based
c              on the slopes of the isopycnal surfaces on the north face
c              of "T" cells.
c       dslope = half length of the interval in which "ahisop" changes
c              with a steep slope from about 0.9*"ahisop" to about
c              0.1*"ahisop"
c       slopec = slope at which "ahisop" is equal to half of its 
c              original  value
# endif
c       dptlim = depth limits for the reference pressure levels (in cm). 
c
c     output:
c       The above input is potentially reset via namelist
c
c
c     authors: 
c             R.C.Pacanowski  (implemented Gokhan's version for MOM 2)
c                              rcp@gfdl.gov
c             Gokhan Danabasoglu (pre MOM 2 version including
c                                 Gent/Mcwilliams transport velocities) 
c                                 gokhan@isis.cgd.ucar.edu
c             K. Dixon        (implemented Cox version in MOM 1)
c                              kd@gfdl.gov
c             M. Cox          (original isopycnal code in COX model)
c=======================================================================
c
      logical hmixset, error
      character*(*) cifdef(ifdmax)
c
# include "param.h"
# include "accel.h"
# include "coord.h"
# include "grdvar.h"
# include "iounit.h"
# include "isopyc.h"
# include "levind.h"
# include "mw.h"
# include "switch.h"
# include "vmixc.h"
c
      namelist /isopyc/ slmxr, ahisop, athkdf, dptlim, dslope, slopec
c
      write (stdout,'(/,20x,a,/,20x,a,/)')
     & 'I S O P Y C M I X    I N I T I A L I Z A T I O N'
# ifdef gent_mcwilliams
     &,'Gent-McWilliams version including thickness diffusion'  
# else
     &,'Redi/Cox version. (Gent-McWilliams without thickness diffusion)'  
# endif
c
c-----------------------------------------------------------------------
c     USER INPUT
c     initialize variables (all mixing units are cm**2/sec.)
c-----------------------------------------------------------------------
c
      slmxr  = 100.0
      ahisop = 2.e7
c
c     define the isopycnal thickness diffusion coefficient
c
      athkdf = 2.0e7
c
c     reference pressure level intervals are defined (see "isopyc.h").
c     "dptlim" must have "nrpl+1" elements (see "param.h"). also,
c     "dptlim" are in cm.
c
c
c     REMARK: the first and the last elements of "dptlim" must be the 
c             depth at the top (0cm) and the maximum bottom depth,
c             respectively. Also, the elements of "dptlim" must be in
c             increasing order.
c
      dptlim(1) = 0.0e2
      dptlim(2) = 1000.0e2
      dptlim(3) = 2000.0e2
      dptlim(4) = 3000.0e2
      dptlim(5) = 4000.0e2
      dptlim(6) = zw(km)
c
# ifdef isopycmixspatialvar
      dslope = 0.001
      slopec = 0.004
# endif
c
c-----------------------------------------------------------------------
c     provide for namelist over-ride of above settings + documentation
c-----------------------------------------------------------------------
c
      call getunit (io, 'namelist', 'fsr')
      read (io,isopyc,end=100)
100   continue
      write (stdout,isopyc)
      call relunit (io)
      call getunit (iodoc, 'document.dta', 'f s a')
      write (iodoc, isopyc)
      call relunit (iodoc)
c
c-----------------------------------------------------------------------
c     add character string to "ifdef option list" indicating that this
c     option is enabled
c-----------------------------------------------------------------------
c
      nifdef = nifdef + 1
      cifdef(nifdef) = 'isopycmix      '
c
c-----------------------------------------------------------------------
c     check for problems
c-----------------------------------------------------------------------
c
      if (nrpl+1 .ne. 6) then
        error = .true.
	write (stdout,*) '=>Error: "dptlim" is set for nrpl=5 in isopi'
     &, ' but nrpl=',nrpl,'. reset dptlim in isop0i or change nrpl'
      endif     
# if defined gent_mcwilliams && !defined isopycmix
        write (stdout,*)
     & '==> Error: "isopycmix" must be enabled with "gent_mcwilliams"'
        error = .true.
# endif
c
c-----------------------------------------------------------------------
c     print out tracer diffusion coefficients in isopycnal case
c-----------------------------------------------------------------------
c
# ifdef held_larichev
      write(stdout,9102) athkdf, ah, slmxr, am
      call getunit (iodoc, 'document.dta', 'f s a')
      write (iodoc,'(a,e14.7)') 'am=', am, 'ah=',ah
     &,     'athkdf=',athkdf, 'slmxr=',slmxr
      call relunit (iodoc)
c
9102  format(/' ahisop is from held_larichev along isopyncal surfaces '
     &,  '(cm**2/sec) '/' athkdf = ',e12.6,' isopycnal thickness '
     &,'diffusion (cm**2/sec) '/' ah = ',e12.6,' cm**2/sec for '
     &,'background horizontal tracer diffusion'/' slmxr = ',e12.6
     &,', to limit the slope used in computing mixing tensors',/
     &,' am = ',e12.6,' cm**2/sec for horizontal mixing of momentum'/)
# else
      write(stdout,9102) ahisop, athkdf, ah, slmxr, am
      call getunit (iodoc, 'document.dta', 'f s a')
      write (iodoc,'(a,e14.7)') 'am=', am, 'ah=',ah
     &,     'ahisop=',ahisop, 'athkdf=',athkdf, 'slmxr=',slmxr
      call relunit (iodoc)
c
9102  format(/' ahisop = ',e12.6,' along isopyncal tracer mixing '
     &,  '(cm**2/sec) '/' athkdf = ',e12.6,' isopycnal thickness '
     &,'diffusion (cm**2/sec) '/' ah = ',e12.6,' cm**2/sec for '
     &,'background horizontal tracer diffusion'/' slmxr = ',e12.6
     &,', to limit the slope used in computing mixing tensors',/
     &,' am = ',e12.6,' cm**2/sec for horizontal mixing of momentum'/)
# endif
c
c-----------------------------------------------------------------------
c     store the square root of the tracer timestep acceleration values
c     into variable "dtxsqr" for use in isopycnal mixing
c-----------------------------------------------------------------------
c
      do k=1,km
       dtxsqr(k) = sqrt(dtxcel(k))
      enddo
c
c-----------------------------------------------------------------------
c     determine the isopycnal reference pressure levels for the "t"
c     grid point levels, using the depths at the "t" grid points as the
c     reference depth (pressure)
c-----------------------------------------------------------------------
c
      do k=1,km
        do m=2,nrpl+1
          if (zt(k).gt.dptlim(m-1) .and. zt(k).le.dptlim(m)) then
            kisrpl(k) = m-1
            go to 101
          endif
        enddo
101     continue
        if (kisrpl(k) .lt. 1 .or. kisrpl(k) .gt. nrpl) then
          write (stdout,9100) kisrpl(k), k
          error = .true.
        endif
      enddo
9100  format (/,' =>Error: kisrpl is ',i3,' at k ',1x,i3,' in isopyc.F')
c
c-----------------------------------------------------------------------
c     the indices used in isopycnal mixing indicating the location of
c     the reference pressure levels in the 20-level table of polynomial
c     expansion variables are computed
c
c     REMARK: because the polynomial expansion coefficients are
c             functions of the reference potential temperature and
c             salinity profiles, at the reference pressure level
c             the corresponding potential temperature and salinity
c             values will be used.
c-----------------------------------------------------------------------
c
      do m=1,nrpl
        krplin(m) = 0
      enddo
c
      do m=2,nrpl+1
        dptmid = p5*(dptlim(m-1)+dptlim(m))
        if (dptmid .le. zt(1)) then
          krplin(m-1) = 1
        elseif (dptmid .gt. zt(km)) then
          krplin(m-1) = km
        elseif (dptmid.gt.zt(1) .and. dptmid.le.zt(km)) then
          do k=2,km
            if (zt(k) .ge. dptmid) then
              t1 = zt(k)-dptmid
              t2 = dptmid-zt(k-1)
              if (t1 .gt. t2) then
                krplin(m-1) = k-1
              else
                krplin(m-1) = k
              endif
              go to 102 
            endif
          enddo
102       continue
        endif
        if (krplin(m-1) .lt. 1 .or. krplin(m-1) .gt. km) then
          write (stderr,9110) krplin(m-1),m-1
          error = .true.
        endif
      enddo
c
      write (stdout,96) (kisrpl(k),k=1,km)
      write (stdout,97) (krplin(m),m=1,nrpl)
c
c-----------------------------------------------------------------------
c     the isopycnal diffusion coefficient may be a function of depth. in 
c     the default configuration, the isopycnal diffusion coefficient is 
c     a constant: "fzisop", which multiplies "ahisop", is set to unity. 
c     if "ahisop" varies in the vertical, "fzisop" should contain this 
c     variation. the value of "ahisop" should be adjusted accordingly.
c-----------------------------------------------------------------------
c
      do k=1,km
        fzisop(k) = c1
      enddo
c
      write (stdout,98)
      write (stdout,99) (fzisop(k),k=1,km)
c
  96  format (/,' isopycnal reference pressure levels (kisrpl) = ',
     &        20(1x,i4))
  97  format (/,' reference pressure level indices (krplin) = ',
     &        20(1x,i4))
  98  format (/,' vertical variation of "ahisop" (fzisop) = ')
  99  format (5(1x,e12.6))
9110  format (/,' =>Error: krplin is ',i3,' at m ',1x,i3)
c
c-----------------------------------------------------------------------
c     initialize arrays
c-----------------------------------------------------------------------
c
      do j=jsmw,jemw
        do i=1,imt
	  do k=1,km
            K1(i,k,j,3)   = c0
            K3(i,k,j,1)   = c0
            K3(i,k,j,2)   = c0
            K3(i,k,j,3)   = c0
# ifdef gent_mcwilliams
            adv_vetiso(i,k,j) = c0
#  ifdef isopycmixspatialvar
            dciso1(i,k,j)  = c0
#  endif
# endif
          enddo
        enddo
      enddo
c
      do j=1,jemw
        do i=1,imt
	  do k=1,km
            K2(i,k,j,3)    = c0
# ifdef gent_mcwilliams
            adv_vntiso(i,k,j) = c0
#  ifdef isopycmixspatialvar
            dciso2(i,k,j)  = c0
#  endif
# endif
          enddo
        enddo
      enddo
c
      do m=1,nrpl
        do j=1,jmw
          do i=1,imt
	    do k=1,km
	      rhoi(i,k,j,m) = c0
	    enddo
	  enddo
        enddo
      enddo
c
      do m=1,3
        do j=1,jemw
          do i=1,imt
	    do k=1,km+1
	      e(i,k,j,m) = c0
	    enddo
	  enddo
        enddo
      enddo
c
# ifdef gent_mcwilliams
      do j=jsmw,jemw
        do k=0,km
          do i=1,imt
            adv_vbtiso(i,k,j) = c0
            adv_fbiso(i,k,j)  = c0
	  enddo
        enddo
      enddo
# endif
      return
      end




      subroutine isopyc (joff, js, je, is, ie)
c
c=======================================================================
c
c     Compute the isopycnal mixing tensor components and the
c     isopycnal advection velocities which parameterize the effect
c     of eddies on the isopycnals.
c
c
c     Mixing tensor "K" is ...
c
c          | 1.0            K1(,,,2)          K1(,,,3) | 
c          |                                           |
c     K =  | K2(,,,1)        1.0              K2(,,,3) |
c          |                                           |
c          | K3(,,,1)       K3(,,,2)          K3(,,,3) |
c
c     where K1(,,,2) and K2(,,,1) are set to 0.0 (neglected)
c
c
c     input:
c       joff = offset relating row "j" in the MW to latitude "jrow"
c       js   = starting row within the MW for calculations
c       je   = ending row within the MW for calculations
c       is   = starting index longitude within the MW
c       ie   = ending index longitude within the MW
c
c     output:
c       rhoi = density at tau-1 referenced to pressure levels
c       K1   = tensor components (1,2), (1,3) centered on east face
c              of "T" cells
c       K2   = tensor components (2,1), (2,3) centered on north face
c              of "T" cells
c       K3   = tensor components (3,1), (3,2), (3,3) centered on
c              bottom face of "T" cells  
# ifdef gent_mcwilliams
c       adv_vetiso = isopycnal advective vel on east face of "T" cell
c       adv_vntiso = isopycnal advective vel on north face of "T" cell
c               (Note: this includes the cosine factor as in "adv_vnt")
c       adv_vbtiso = isopycnal advective vel on bottom face of "T" cell       
# endif
c
c=======================================================================
c
# include "param.h"
# include "accel.h"
# include "grdvar.h"
# include "iounit.h"
# include "isopyc.h"
# include "levind.h"
# include "mw.h"
# include "switch.h"
# include "state.h"
# include "dens.h"
c
# ifdef timing
      call tic ('isopycnal ', 'isopyc')
# endif
c
c-----------------------------------------------------------------------
c     set local constants 
c-----------------------------------------------------------------------
c
      istrt = max(2,is)
      iend  = min(imt-1,ie)
c
c-----------------------------------------------------------------------
c     after the first MW, move variables from the top two rows to the  
c     bottom two rows of the MW to eliminate redundant calculation
c-----------------------------------------------------------------------
c
      if (joff .ne. 0) then
        call moveiso (istrt-1, iend+1)
      endif
c
c-----------------------------------------------------------------------
c     compute normalized densities for each isopycnal reference pressure
c     level using a 3rd order polynomial fit to the equation of state.
c     for each isopycnal reference pressure level, the same reference
c     potential temperature, reference salinity and expansion coeff
c     values are used at all of the vertical levels.
c
c     Note: this density is used for the mixing tensor in both the 
c     Redi/Cox and Gent/McWilliams options
c
c     first MW     : calculate for rows 1..jmw
c     1 < MW < last: calculate for rows jsmw+1..jmw 
c     last MW      : calculate for rows jsmw+1..<=jmw
c-----------------------------------------------------------------------
c
      do m=1,nrpl
        tref = to(krplin(m))
        sref = so(krplin(m))
	do j=js, je
          do k=1,km
            do i=1,imt
	      tprime = t(i,k,j,1,taum1) - tref
	      sprime = t(i,k,j,2,taum1) - sref
              rhoi(i,k,j,m) = dens (tprime, sprime, krplin(m))
             enddo
          enddo
        enddo
      enddo
#ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     & "=> In isopyc: rhoi from js=",js," je=",je," joff=",joff
     &," jrows=",js+joff," to ",je+joff
#endif
c
c-----------------------------------------------------------------------
c     evaluate K2(,,3) centered on the northern face of "T" cells
c     first MW     : calculate for rows 1..jemw
c     1 < MW < last: calculate for rows jsmw..jemw 
c     last MW      : calculate for rows jsmw..<=jemw
c-----------------------------------------------------------------------
c
      call k2_3 (joff, max(js-1,1), je-1, is, ie)
c
c-----------------------------------------------------------------------
c     evaluate K1(,,3) centered on eastern face of "T" cells
c     1 <= MW < last: calculate for rows jsmw..jemw 
c     last MW       : calculate for rows jsmw..<=jemw
c-----------------------------------------------------------------------
c
      call k1_3 (joff, max(js-1,jsmw), je-1, is, ie)
c
c-----------------------------------------------------------------------
c     evaluate K3(,,1..3) centered on bottom face of "T" cells
c     1 <= MW < last: calculate for rows jsmw..jemw 
c     last MW       : calculate for rows jsmw..<=jemw
c-----------------------------------------------------------------------
c
      call k3_123 (joff, max(js-1,jsmw), je-1, is, ie)
c
# ifdef gent_mcwilliams
c
c-----------------------------------------------------------------------
c     compute isopycnal advective velocities for tracers
c     first MW     : calculate (viso)           for rows 1..jemw
c                              (uiso,wiso)      for rows 2..jemw
c     1 < MW < last: calculate (uiso,viso,wiso) for rows jsmw..jemw 
c     last MW      : calculate (uiso,viso,wiso) for rows jsmw..<=jemw
c-----------------------------------------------------------------------
c
      call isopyc_adv (joff, max(js-1,1), je-1, is, ie)
c
# endif
c
# ifdef timing
      call toc ('isopycnal ', 'isopyc')
# endif
      return
      end


      subroutine moveiso (is, ie) 
c
c=======================================================================
c     move the MW up (northward) by moving data from the last two rows
c     into the first two rows.
c
c     input:
c      is = starting longitude index in the MW
c      ie = ending longitude index in the MW
c=======================================================================
c
#include "param.h"
#include "hmixc.h"
#include "isopyc.h"
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
        do n=1,nrpl
          do k=1,km
	    do i=is,ie
	      rhoi(i,k,jto,n) = rhoi(i,k,jfrom,n)
	    enddo
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
	      K2(i,k,jto,3)       = K2(i,k,jfrom,3)
# ifdef gent_mcwilliams
	      adv_vntiso(i,k,jto)  = adv_vntiso(i,k,jfrom)
# endif
            enddo
          enddo
	endif
# ifdef trace_indices
        write (stdout,'(4x,2(a,i4))')
     & "=> In moveiso: moving variables on row ",jfrom," to row ",jto
# endif
      enddo
c
      return
      end




      subroutine k1_3 (joff, js, je, is, ie)
c
c=======================================================================
c     
c     compute "K1(,,3)" at the center of the eastern face of "T" cells
c     use "c1e10" to keep the exponents in range.
c
c=======================================================================
c
# include "param.h"
# include "accel.h"
# include "coord.h"
# include "grdvar.h"
# include "iounit.h"
# include "isopyc.h"
# include "levind.h"
# include "mw.h"
# include "switch.h"
c
c-----------------------------------------------------------------------
c     set local constants 
c-----------------------------------------------------------------------
c
      c1e10 = 1.0e10
      eps   = 1.0e-25
c
c-----------------------------------------------------------------------
c     d(rho_barx_barz)/dz centered on eastern face of "t" cells
c     Note: values involving ocean surface and ocean bottom are
c           estimated afterwards using a linear extrapolation
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=2,kmm1
          m = kisrpl(k)
          fxd = c1e10*p5*dzt2r(k)
          do i=1,imtm1
            e(i,k,j,3) = fxd*(rhoi(i  ,k-1,j,m) - rhoi(i  ,k+1,j,m)
     &                       +rhoi(i+1,k-1,j,m) - rhoi(i+1,k+1,j,m))
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     linearly extrapolate densities to ocean surface for calculation
c     of d(rho_barx_barz)/dz involving level 1.
c
c     REMARK: requires min(kmt(i,jrow)) = 2 cells in ocean.
c-----------------------------------------------------------------------
c
      k   = 1
      fxd = c1e10*dztr(k)
      fxe = dzw(k-1)+dzw(k)
      m   = kisrpl(k)
      do j=js,je
        do i=1,imtm1
          fxa        = p5*(rhoi(i,k+1,j,m) + rhoi(i+1,k+1,j,m))
          fxb        = p5*(rhoi(i,k,j,m) + rhoi(i+1,k,j,m))
          fxc        = dzwr(k)*(fxb*fxe - fxa*dzw(k-1))
          e(i,k,j,3) = fxd*(fxc - p5*(fxa+fxb))
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     linearly extrapolate densities to ocean bottom for calculation
c     of d(rho_barx_barz)/dz involving bottom level.
c-----------------------------------------------------------------------
c
      do j=js,je
        do i=1,imtm1
          e(i,km,j,3) = c0
        enddo
      enddo
c
      do j=js,je
        jrow = j + joff
        do i=1,imtm1
          k = min(kmt(i,jrow),kmt(i+1,jrow))
          if (k .ne. 0) then
            fxe        = dzw(k-1)+dzw(k)
            m          = kisrpl(k)
            fxa        = p5*(rhoi(i,k-1,j,m) + rhoi(i+1,k-1,j,m))
            fxb        = p5*(rhoi(i,k  ,j,m) + rhoi(i+1,k  ,j,m))
            fxc        = dzwr(k-1)*(fxb*fxe - fxa*dzw(k))
            e(i,k,j,3) = dztr(k)*c1e10*(p5*(fxa+fxb) - fxc)
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     "e(,,,1)" = d(rho)/dx centered on east face of "T" cells
c     "e(,,,2)" = d(rho_barx_bary)/dy centered on east face of "T" cells
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        do k=1,km
          m = kisrpl(k)
          do i=1,imtm1
            e(i,k,j,1) = tmask(i,k,j)*tmask(i+1,k,j)*cstr(jrow)*dxur(i)
     &                     *c1e10*(rhoi(i+1,k,j,m) - rhoi(i,k,j,m))
            e(i,k,j,2) = dyt4r(jrow)*c1e10*(
     &                     rhoi(i  ,k,j+1,m) - rhoi(i  ,k,j-1,m)
     &                   + rhoi(i+1,k,j+1,m) - rhoi(i+1,k,j-1,m))
           enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     if any one of the 4 neighboring corner grid points is a land point,
c     set "e(i,k,j,2)" to zero. note that "e(i,k,j,2)" will be used
c     only in the slope check.
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=1,km
          do i=1,imtm1
            olmask = tmask(i,k,j-1)*tmask(i,k,j+1)*tmask(i+1,k,j-1)
     &   	    *tmask(i+1,k,j+1)
            if (olmask .eq. c0)  e(i,k,j,2) = c0
           enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     impose zonal boundary conditions at "i"=1 and "imt"
c-----------------------------------------------------------------------
c
      do j=js,je
        call setbcx (e(1,1,j,1), imt, km)
        call setbcx (e(1,1,j,2), imt, km)
        call setbcx (e(1,1,j,3), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     compute "K1", using "slmxr" to limit vertical slope of isopycnal
c     to guard against numerical instabilities. 
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=1,km
# ifdef isopycmixspatialvar
          do i=1,imt
            fxa = c0
            fxb = sign(1.,e(i,k,j,3))/(abs(e(i,k,j,3))+eps)
            slope = fxb*sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)
            if (slope .le. c0 .and. slope .ge. (-c1/slmxr)) then
              fxa = p5*(c1+tanh((slope+slopec)/dslope))
            endif
            dciso1(i,k,j) = fxa
            K1(i,k,j,3)  = -fxb*e(i,k,j,1)*dciso1(i,k,j)
          enddo
# else
          do i=1,imt
            chkslp = -sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)*slmxr*dtxsqr(k)
            if (e(i,k,j,3) .gt. chkslp)  e(i,k,j,3) = chkslp
          enddo
          do i=1,imt
            K1(i,k,j,3) = (-e(i,k,j,1)*e(i,k,j,3)*fzisop(k))
     &                     /(e(i,k,j,3)**2+eps)
          enddo
# endif
        enddo
      enddo
#ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     & "=> In K1_3: js=",js," je=",je," joff=",joff
     &," jrows=",js+joff," to ",je+joff
#endif
      return
      end





      subroutine k2_3 (joff, js, je, is, ie)
c     
c=======================================================================
c     compute "K2(,,3)" at the center of the northern face of "T" cells
c     use "c1e10" to keep the exponents in range.
c=======================================================================
c
# include "param.h"
# include "accel.h"
# include "coord.h"
# include "grdvar.h"
# include "iounit.h"
# include "isopyc.h"
# include "levind.h"
# include "mw.h"
# include "switch.h"
c
c-----------------------------------------------------------------------
c     set local constants 
c-----------------------------------------------------------------------
c
      c1e10 = 1.0e10
      eps   = 1.0e-25
c
c-----------------------------------------------------------------------
c     d(rho_bary_barz)/dz centered on northern face of "T" cells
c     Note: values involving ocean surface and ocean bottom are
c           estimated afterwards using a linear extrapolation
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=2,kmm1
          m = kisrpl(k)
          fxd = c1e10*p5*dzt2r(k)
          do i=2,imtm1
            e(i,k,j,3) = fxd*(rhoi(i,k-1,j  ,m) - rhoi(i,k+1,j  ,m)
     &                       +rhoi(i,k-1,j+1,m) - rhoi(i,k+1,j+1,m))
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     linearly extrapolate densities to ocean surface for calculation
c     of d(rho_bary_barz)/dz involving level 1.
c
c     REMARK: requires min(kmt(i,jrow)) = 2 cells in ocean.
c-----------------------------------------------------------------------
c
      k   = 1
      fxd = c1e10*dztr(k) 
      fxe = dzw(k-1)+dzw(k)
      m   = kisrpl(k)
      do j=js,je
        do i=2,imtm1
          fxa        = p5*(rhoi(i,k+1,j,m) + rhoi(i,k+1,j+1,m))
          fxb        = p5*(rhoi(i,k  ,j,m) + rhoi(i,k  ,j+1,m))
          fxc        = dzwr(k)*(fxb*fxe - fxa*dzw(k-1))
          e(i,k,j,3) = fxd*(fxc - p5*(fxa+fxb))
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     linearly extrapolate densities to ocean bottom for calculation
c     of d(rho_bary_barz)/dz involving bottom level.
c-----------------------------------------------------------------------
c
      do j=js,je
        do i=2,imtm1
          e(i,km,j,3) = c0
        enddo
      enddo
c
      do j=js,je
        jrow = j + joff
        do i=2,imtm1
          k = min(kmt(i,jrow),kmt(i,jrow+1))
          if (k .ne. 0) then
            fxe        = dzw(k-1)+dzw(k)
            m          = kisrpl(k)
            fxa        = p5*(rhoi(i,k-1,j,m) + rhoi(i,k-1,j+1,m))
            fxb        = p5*(rhoi(i,k  ,j,m) + rhoi(i,k  ,j+1,m))
            fxc        = dzwr(k-1)*(fxb*fxe - fxa*dzw(k))
            e(i,k,j,3) = dztr(k)*c1e10*(p5*(fxa+fxb) - fxc)
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     "e(,,,1)" = d(rho_barx_bary)/dx centered on north face of "T" cells
c     "e(,,,2)" = d(rho)/dy on north face of "T" cells
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        do k=1,km
          m = kisrpl(k)
          do i=2,imtm1
            e(i,k,j,1) = csur(jrow)*dxt4r(i)*c1e10*(
     &                    rhoi(i+1,k,j+1,m) - rhoi(i-1,k,j+1,m)
     &                  + rhoi(i+1,k,j  ,m) - rhoi(i-1,k,j  ,m))
            e(i,k,j,2) = tmask(i,k,j)*tmask(i,k,j+1)*dyur(jrow)*c1e10
     &                  *(rhoi(i,k,j+1,m) - rhoi(i,k,j,m))   
           enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     if any one of the 4 neighboring corner grid points is a land point,
c     set "e(i,k,j,1)" to zero. note that "e(i,k,j,1)" will be used
c     only in the slope check.
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=1,km
          do i=2,imtm1
            olmask = tmask(i-1,k,j+1)*tmask(i+1,k,j+1)*tmask(i-1,k,j)
     &              *tmask(i+1,k,j)
          if (olmask .eq. c0)  e(i,k,j,1) = c0
           enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     impose zonal boundary conditions at "i"=1 and "imt"
c-----------------------------------------------------------------------
c
      do j=js,je
        call setbcx (e(1,1,j,1), imt, km)
        call setbcx (e(1,1,j,2), imt, km)
        call setbcx (e(1,1,j,3), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     compute "K2", using "slmxr" to limit vertical slope of isopycnal
c     to guard against numerical instabilities. 
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=1,km
# ifdef isopycmixspatialvar
          do i=1,imt
            fxa = c0
            fxb = sign(1.,e(i,k,j,3))/(abs(e(i,k,j,3))+eps)
            slope = fxb*sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)
            if (slope .le. c0 .and. slope .ge. (-c1/slmxr)) then
              fxa = p5*(c1+tanh((slope+slopec)/dslope))
            endif
            dciso2(i,k,j) = fxa
            K2(i,k,j,3)  = -fxb*e(i,k,j,2)*dciso2(i,k,j)
          enddo
# else
          do i=1,imt
            chkslp = -sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)*slmxr*dtxsqr(k)
            if (e(i,k,j,3) .gt. chkslp)  e(i,k,j,3) = chkslp
          enddo
          do i=1,imt
            K2(i,k,j,3) = (-e(i,k,j,2)*e(i,k,j,3)*fzisop(k))
     &                     /(e(i,k,j,3)**2+eps)
          enddo
# endif
        enddo
      enddo
#ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     & "=> In K2_3: js=",js," je=",je," joff=",joff
     &," jrows=",js+joff," to ",je+joff
#endif
      return
      end



      subroutine k3_123 (joff, js, je, is, ie)
c
c=======================================================================
c     compute K2(,,,1:3) at the center of the bottom face of "T" cells
c     use "c1e10" to keep the exponents in range.
c=======================================================================
c
# include "param.h"
# include "accel.h"
# include "grdvar.h"
# include "iounit.h"
# include "isopyc.h"
# include "levind.h"
# include "mw.h"
# include "switch.h"
# include "vmixc.h"
c
c-----------------------------------------------------------------------
c     set local constants 
c-----------------------------------------------------------------------
c
      c1e10 = 1.0e10
      eps   = 1.0e-25
c
      do j=js,je
        jrow = j + joff
        do k=2,km
          m = kisrpl(k)
          do i=2,imtm1
            e(i,k-1,j,1) = cstr(jrow)*dxt4r(i)*c1e10
     &               *(tmask(i-1,k-1,j)*tmask(i,k-1,j)*(rhoi(i,k-1,j,m)
     &               -rhoi(i-1,k-1,j,m))
     &               +tmask(i,k-1,j)*tmask(i+1,k-1,j)*(rhoi(i+1,k-1,j,m)
     &               -rhoi(i,k-1,j,m))
     &               +tmask(i-1,k,j)*tmask(i,k,j)*(rhoi(i,k,j,m)
     &               -rhoi(i-1,k,j,m))
     &               +tmask(i,k,j)*tmask(i+1,k,j)*(rhoi(i+1,k,j,m)
     &               -rhoi(i,k,j,m)))
            e(i,k-1,j,2) = dyt4r(jrow)*c1e10
     &               *(tmask(i,k-1,j-1)*tmask(i,k-1,j)*(rhoi(i,k-1,j,m)
     &               -rhoi(i,k-1,j-1,m))
     &               +tmask(i,k-1,j)*tmask(i,k-1,j+1)*(rhoi(i,k-1,j+1,m)
     &               -rhoi(i,k-1,j,m))
     &               +tmask(i,k,j-1)*tmask(i,k,j)*(rhoi(i,k,j,m)
     &               -rhoi(i,k,j-1,m))
     &               +tmask(i,k,j)*tmask(i,k,j+1)*(rhoi(i,k,j+1,m)
     &               -rhoi(i,k,j,m)))
            e(i,k-1,j,3) = dzwr(k-1)*tmask(i,k-1,j)*tmask(i,k,j)*c1e10
     &               *(rhoi(i,k-1,j,m) - rhoi(i,k,j,m))
           enddo
        enddo
        k = km
        do i=2,imtm1
	  e(i,k,j,1) = c0
	  e(i,k,j,2) = c0
	  e(i,k,j,3) = c0
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute "K3", using "slmxr" to limit vertical slope of isopycnal
c     to guard against numerical instabilities.  
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=1,km
# ifdef isopycmixspatialvar
          do i=2,imt-1
            fxa = c0
            fxb = sign(1.,e(i,k,j,3))/(abs(e(i,k,j,3))+eps)
            slope = fxb*sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)
            if (slope .le. c0 .and. slope .ge. (-c1/slmxr)) then
              fxa = p5*(c1+tanh((slope+slopec)/dslope))
            endif
            fxc = fxb*fxa
            K3(i,k,j,1) = -fxc*e(i,k,j,1)
            K3(i,k,j,2) = -fxc*e(i,k,j,2)
            K3(i,k,j,3) = fxb*fxb*(e(i,k,j,1)**2+e(i,k,j,2)**2)
     &                    *fxa
          enddo
# else
          do i=2,imt-1
            chkslp = -sqrt(e(i,k,j,1)**2+e(i,k,j,2)**2)*slmxr*dtxsqr(k)
            if (e(i,k,j,3) .gt. chkslp)  e(i,k,j,3) = chkslp
          enddo
          fxa = p5*(fzisop(min(km,k+1))+fzisop(k))
          do i=2,imt-1
            ahfctr = fxa/(e(i,k,j,3)**2+eps)
            K3(i,k,j,1) = -e(i,k,j,3)*e(i,k,j,1)*ahfctr
            K3(i,k,j,2) = -e(i,k,j,3)*e(i,k,j,2)*ahfctr
            K3(i,k,j,3) = (e(i,k,j,1)**2+e(i,k,j,2)**2)*ahfctr
          enddo
# endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     impose zonal boundary conditions at "i"=1 and "imt"
c-----------------------------------------------------------------------
c
      do j=js,je
        call setbcx (K3(1,1,j,1), imt, km)
        call setbcx (K3(1,1,j,2), imt, km)
        call setbcx (K3(1,1,j,3), imt, km)
      enddo
#ifdef trace_indices
      write (stdout,'(2x,5(a,i4))')
     & "=> In K3_123: js=",js," je=",je," joff=",joff
     &," jrows=",js+joff," to ",je+joff
#endif
      return
      end               


# ifdef gent_mcwilliams
      subroutine isopyc_adv (joff, js, je, is, ie)
c
c=======================================================================
c     compute isopycnal transport velocities.
c=======================================================================
c
# include "param.h"
# include "coord.h"
# include "grdvar.h"
# include "isopyc.h"
# include "levind.h"
# include "mw.h"
# include "switch.h"
c
c-----------------------------------------------------------------------
c     compute the meridional component of the isopycnal mixing velocity
c     at the center of the northern face of the "t" cells.
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        do k=2,kmm1
          fxa = -dzt2r(k)*athkdf*csu(jrow)
          do i=1,imt
            adv_vntiso(i,k,j) = fxa*tmask(i,k,j)*tmask(i,k,j+1)*(
     &                          K2(i,k-1,j,3) - K2(i,k+1,j,3))
          enddo
        enddo
      enddo
c
c     consider the top and bottom levels. "K2" is assumed to be zero
c     at the ocean top and bottom.
c
      k = 1
      fxa = -dzt2r(k)*athkdf
      do j=js,je
        jrow = j + joff
        do i=1,imt
          adv_vntiso(i,k,j) = -fxa*tmask(i,k,j)*tmask(i,k,j+1)*csu(jrow)
     &                        *(K2(i,k,j,3) + K2(i,k+1,j,3))
        enddo
      enddo
c
      do j=js,je
        do i=1,imt
          adv_vntiso(i,km,j) = c0
        enddo
      enddo
c
      do j=js,je
        jrow = j + joff
        do i=1,imt
          k = min(kmt(i,jrow),kmt(i,jrow+1))
          if (k .ne. 0) then
            adv_vntiso(i,k,j) = -dzt2r(k)*athkdf*csu(jrow)*tmask(i,k,j)
     &                   *tmask(i,k,j+1)*(K2(i,k,j,3) + K2(i,k-1,j,3))
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute the zonal component of the isopycnal mixing velocity
c     at the center of the eastern face of the "t" grid box.
c-----------------------------------------------------------------------
c
      jstrt = max(js,jsmw)
      do j=jstrt,je
        do k=2,kmm1
          fxa = -dzt2r(k)*athkdf
          do i=1,imtm1
            adv_vetiso(i,k,j) = fxa*tmask(i,k,j)*tmask(i+1,k,j)
     &                          *(K1(i,k-1,j,3) - K1(i,k+1,j,3))
          enddo
        enddo
      enddo
c
c     consider the top and bottom levels. "K1" is assumed to be zero
c     at the ocean top and bottom.
c
      k = 1
      fxa = -dzt2r(k)*athkdf
      do j=jstrt,je
        do i=1,imtm1
          adv_vetiso(i,k,j) = -fxa*tmask(i,k,j)*tmask(i+1,k,j)
     &                        *(K1(i,k,j,3)+K1(i,k+1,j,3))
        enddo
      enddo
c
      do j=jstrt,je
        do i=1,imtm1
          adv_vetiso(i,km,j) = c0
        enddo
      enddo
c
      do j=jstrt,je
        jrow = j + joff
        do i=1,imtm1
          k = min(kmt(i,jrow),kmt(i+1,jrow))
          if (k .ne. 0) then
            adv_vetiso(i,k,j) = -dzt2r(k)*athkdf*tmask(i,k,j)
     &                     *tmask(i+1,k,j)*(K1(i,k,j,3)+K1(i,k-1,j,3))
          endif
        enddo
      enddo
c
c     set the boundary conditions
c
      do j=jstrt,je
	call setbcx (adv_vetiso(1,1,j), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     compute the vertical component of the isopycnal mixing velocity
c     at the center of the bottom face of the "t" cells, using the
c     continuity equation for the isopycnal mixing velocities
c-----------------------------------------------------------------------
c
      do j=jstrt,je
        do i=1,imt
          adv_vbtiso(i,0,j) = c0
        enddo
      enddo
c
      do j=jstrt,je
        jrow = j + joff
        do k=1,kmm1
          do i=2,imt
            adv_vbtiso(i,k,j) = dzt(k)*cstr(jrow)*(
     &      (adv_vetiso(i,k,j) - adv_vetiso(i-1,k,j))*dxtr(i) + 
     &      (adv_vntiso(i,k,j) - adv_vntiso(i,k,j-1))*dytr(jrow)) 
          enddo
        enddo
      enddo
c
      do j=jstrt,je
        do k=1,kmm1
          do i=2,imt
            adv_vbtiso(i,k,j) = adv_vbtiso(i,k,j) + adv_vbtiso(i,k-1,j)
          enddo
        enddo
      enddo
c
      do j=jstrt,je
        jrow = j + joff
        do i=2,imt
          adv_vbtiso(i,kmt(i,jrow),j) = c0
        enddo
      enddo
c
c     set the boundary conditions
c
      do j=jstrt,je
	call setbcx (adv_vbtiso(1,0,j), imt, km+1)
      enddo
c
# ifdef trace_indices
      write (stdout,'(2x,3(a,i4),6(/4x,4(a,i4),a))')
     & "=> In isopyc_adv: js=",js," je=",je," joff=",joff
     &, " adv_vntiso calculated on jrows ",js+joff," through ",je+joff
     &, " (rows ",js," through ",je," in the MW)"
     &, " adv_vetiso calculated on jrows ",jstrt+joff," through ",je+joff
     &, " (rows ",jstrt," through ",je," in the MW)"
     &, " adv_vbtiso calculated on jrows ",jstrt+joff," through ",je+joff
     &, " (rows ",jstrt," through ",je," in the MW)"
# endif
c#  define debug_adv_vel_iso
#  ifdef debug_adv_vel_iso
c
c-----------------------------------------------------------------------
c     to inspect the divergence of "t" cells at point (ipt,jpt)
c     for all levels
c-----------------------------------------------------------------------
c
      sml  = -1.0
      ipt = 61
      jpt = 16
      do j=jstrt,je
        jrow = j + joff
        if (jrow .eq. jpt .and. eorun) then
          write (stdout,*) 'Showing isopycnal advection velocities and '
     &,   'divergence for jrow=',jpt, ' i=',ipt
	  do i=ipt,ipt
            do k=1,km
	      divgt =
     &        ((adv_vetiso(i,k,j) - adv_vetiso(i-1,k,j))
     &        *cstr(jrow)*dxtr(i)
     & +      (adv_vntiso(i,k,j) - adv_vntiso(i,k,j-1))
     &        *cstr(jrow)*dytr(jrow)
     & +      (adv_vbtiso(i,k-1,j) - adv_vbtiso(i,k,j))
     &        *dztr(k))*tmask(i,k,j)
              if (abs(divgt) .ge. sml) then
                write (stdout,98) i,k,divgt
     &,                         adv_vetiso(i,k,j), adv_vetiso(i-1,k,j)
     &,                         adv_vntiso(i,k,j), adv_vntiso(i,k,j-1)
     &,                         adv_vbtiso(i,k-1,j), adv_vbtiso(i,k,j)
                write (stdout,*) ' tmask=',tmask(i,k,j)
     &,                        ' kmt=',kmt(i,jrow),' kmu='
     &,                        kmu(i,jrow), ', jrow=',jrow
                write (stdout,*) ' '
	      endif
            enddo
	    write (stdout,*) 'adv_vbtiso(i,0,j) = ',adv_vbtiso(i,0,j)
	    write (stdout,'(///)')
	  enddo
        endif
      enddo
98      format (1x,'tp: i=',i3,', k=',i2,'divtiso=',7e14.7)
#  endif
      return
      end
# endif




      subroutine isoflux (joff, js, je, is, ie, mtrace)
c
c=======================================================================
c     isopycnal diffusive tracer fluxes are computed.
c=======================================================================
c
# include "param.h"
# include "coord.h"
# include "grdvar.h"
# include "hmixc.h"
# include "isopyc.h"
# include "levind.h"
# include "mw.h"
# include "vmixc.h"
#ifdef debug_tensor
# include "tmngr.h"
#endif
      dimension temp(imt,km,jmw)
c
# ifdef timing
      call tic ('isopycnal ', 'isoflux')
# endif
c
      m = mtrace
c
c-----------------------------------------------------------------------
c     first compute the vertical tracer flux "temp" at the northern
c     face of "t" cells.
c-----------------------------------------------------------------------
c
      do j=js-1,je
        do k=2,kmm1
          do i=1,imt
            temp(i,k,j) = p5*dzt2r(k)*(
     &                      t(i,k-1,j+1,m,taum1) - t(i,k+1,j+1,m,taum1)
     &                     +t(i,k-1,j  ,m,taum1) - t(i,k+1,j  ,m,taum1))
     &                     
           enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     now consider the top level, assuming that the surface tracer 
c     values are the same as the ones at "k"=1
c-----------------------------------------------------------------------
c
      k = 1
      do j=js-1,je
        do i=1,imt
          temp(i,k,j) = p5*dzt2r(k)*(
     &                     t(i,k,j+1,m,taum1) - t(i,k+1,j+1,m,taum1)
     &                    +t(i,k,j  ,m,taum1) - t(i,k+1,j  ,m,taum1))
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     finally, consider the bottom level. the extrapolative estimator
c     is used to compute the tracer values at the ocean bottom.
c-----------------------------------------------------------------------
c
      do j=js-1,je
        do i=1,imt
          temp(i,km,j) = c0
        enddo
      enddo
c
      do j=js-1,je
        jrow = j + joff
        do i=1,imt
          k = min(kmt(i,jrow),kmt(i,jrow+1))
          if (k .ne. 0) then
            fxe = dzw(k-1)+dzw(k)
            fxa = p5*(t(i,k-1,j+1,m,taum1) + t(i,k-1,j,m,taum1))
            fxb = p5*(t(i,k  ,j+1,m,taum1) + t(i,k  ,j,m,taum1))
            fxc = dzwr(k-1)*(fxb*fxe-fxa*dzw(k))
            temp(i,k,j) = dztr(k)*(p5*(fxa+fxb) - fxc)
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute of meridional tracer flux at northern face of "T" cells.
c     first calculate the effects of purely horizontal diffusion, using
c     the sum of the along isopycnal and background horizontal diffusion
c     coefficients as the total diffusion coefficient.     
c     add in the effects of the along isopycnal diffusion computed
c     using "K2" component of the tensor and apply land/sea masks
c-----------------------------------------------------------------------
c
      do j=js-1,je
        jrow = j + joff
        call setbcx (diff_fn(1,1,j), imt, km)
        do k=1,km
          do i=1,imt
#ifdef held_larichev
            ahisop = hl_n(i,j)
#endif
#ifdef isopycmixspatialvar
            fxa = ahisop*dciso2(i,k,j)
#else
            fxa = ahisop*fzisop(k)
#endif
            diff_fn(i,k,j) = (diff_fn(i,k,j) + fxa*dyur(jrow)*(
     &                        t(i,k,j+1,m,taum1)-t(i,k,j,m,taum1))
     &                       + ahisop*K2(i,k,j,3)*temp(i,k,j)
     &                       )*csu(jrow)*tmask(i,k,j+1)*tmask(i,k,j)
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute the vertical tracer flux "temp" at the eastern
c     face of "t" cells.
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=2,kmm1
          do i=1,imtm1
            temp(i,k,j) = p5*dzt2r(k)*(
     &                    t(i+1,k-1,j,m,taum1) - t(i+1,k+1,j,m,taum1)
     &                   +t(i  ,k-1,j,m,taum1) - t(i  ,k+1,j,m,taum1))
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     now consider the top level, assuming that the surface tracer 
c     values are the same as the ones at "k"=1
c-----------------------------------------------------------------------
c
      k = 1
      do j=js,je
        do i=1,imtm1
          temp(i,k,j) = p5*dzt2r(k)*(
     &                  t(i+1,k,j,m,taum1) -t(i+1,k+1,j,m,taum1)
     &                 +t(i  ,k,j,m,taum1) -t(i  ,k+1,j,m,taum1))
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     finally, consider the bottom level. the extrapolative estimator
c     is used to compute the tracer values at the ocean bottom.
c-----------------------------------------------------------------------
c
      do j=js,je
        do i=1,imtm1
          temp(i,km,j) = c0
        enddo
      enddo
c
      do j=js,je
        jrow = j + joff
        do i=1,imtm1
          k = min(kmt(i,jrow),kmt(i+1,jrow))
          if (k .ne. 0) then
            fxe          = dzw(k-1)+dzw(k)
            fxa          = p5*(t(i,k-1,j,m,taum1)+t(i+1,k-1,j,m,taum1))
            fxb          = p5*(t(i,k,j,m,taum1)+t(i+1,k,j,m,taum1))
            fxc          = dzwr(k-1)*(fxb*fxe - fxa*dzw(k))
            temp(i,k,j) = dztr(k)*(p5*(fxa+fxb)-fxc)
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     compute of zonal tracer flux at eastern face of "T" cells.
c     first calculate the effects of purely horizontal diffusion, using
c     the sum of the along isopycnal and background horizontal diffusion
c     coefficients as the total diffusion coefficient.     
c     add in the effects of the along isopycnal diffusion computed
c     using "K1" component of the tensor and apply land/sea masks      
c-----------------------------------------------------------------------
c      
      do j=js,je
        jrow = j + joff
        call setbcx (diff_fe(1,1,j), imt, km)
        do k=1,km
          do i=1,imtm1
#ifdef held_larichev
            ahisop = hl_e(i,j)
#endif
#ifdef isopycmixspatialvar
            fxa = ahisop*dciso1(i,k,j)
#else
            fxa = ahisop*fzisop(k)
#endif
            diff_fe(i,k,j) = (diff_fe(i,k,j) + fxa*cstr(jrow)*dxur(i)*(
     &                        t(i+1,k,j,m,taum1) - t(i,k,j,m,taum1))
     &                        + ahisop*K1(i,k,j,3)*temp(i,k,j)
     &                       )*tmask(i+1,k,j)*tmask(i,k,j)
          enddo
        enddo
        call setbcx (diff_fe(1,1,j), imt, km)
      enddo
c
c-----------------------------------------------------------------------
c     compute the vertical tracer flux "diff_fbiso" containing the K31
c     and K32 components which are to be solved explicitly. The K33
c     component will be treated semi-implicitly by diff_fb.
c-----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
        do k=2,km
          do i=2,imtm1
#ifdef held_larichev
            ahisop = hl_b(i,j)
#endif
            diff_fbiso(i,k-1,j) = tmask(i,k-1,j)*tmask(i,k,j)
     &                  *(ahisop*K3(i,k-1,j,1)*cstr(jrow)*dxt4r(i)*(
     &      tmask(i-1,k,j)*(t(i,k,j,m,taum1)-t(i-1,k,j,m,taum1))
     &     +tmask(i-1,k-1,j)*(t(i,k-1,j,m,taum1)-t(i-1,k-1,j,m,taum1))
     &     +tmask(i+1,k,j)*(t(i+1,k,j,m,taum1)-t(i,k,j,m,taum1))
     &     +tmask(i+1,k-1,j)*(t(i+1,k-1,j,m,taum1)-t(i,k-1,j,m,taum1)))
     &                  +ahisop*K3(i,k-1,j,2)*dyt4r(jrow)*(
     &      tmask(i,k,j-1)*(t(i,k,j,m,taum1)-t(i,k,j-1,m,taum1))
     &     +tmask(i,k-1,j-1)*(t(i,k-1,j,m,taum1)-t(i,k-1,j-1,m,taum1))
     &     +tmask(i,k,j+1)*(t(i,k,j+1,m,taum1)-t(i,k,j,m,taum1))
     &     +tmask(i,k-1,j+1)*(t(i,k-1,j+1,m,taum1)-t(i,k-1,j,m,taum1)))
     &                   )
          enddo
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     now consider the top and bottom levels. at the ocean surface,
c     "diff_fb" is set to zero to reflect the no tracer flux condition.
c     the same condition is also imposed at the ocean bottom.
c-----------------------------------------------------------------------
c
      do j=js,je
        do i=2,imtm1
          diff_fbiso(i,0,j) = c0
          diff_fbiso(i,km,j) = c0
        enddo
      enddo
c
c-----------------------------------------------------------------------
c     impose the boundary conditions at "i"=1 and "imt"
c-----------------------------------------------------------------------
c
      do j=js,je
        call setbcx (diff_fbiso(1,0,j), imt, km+1)
      enddo
c
# ifdef gent_mcwilliams
c
c-----------------------------------------------------------------------
c     compute advective tracer flux at the center of the bottom face of
c     the "T" cells
c-----------------------------------------------------------------------
c
      do j=js,je
        do k=1,kmm1
          do i=2,imt-1
            adv_fbiso(i,k,j) = adv_vbtiso(i,k,j)*
     &                         (t(i,k,j,m,taum1) + t(i,k+1,j,m,taum1))
          enddo
        enddo
      enddo
c
c     now consider the top and bottom boundaries
c
      do j=js,je
        do i=2,imt-1
          adv_fbiso(i,0,j)  = c0
          adv_fbiso(i,km,j) = c0
        enddo
      enddo
# endif
c
#ifdef trace_indices
      write (stdout,'(2x,9(a,i4))')
     & "=> In isoflux: n=",n," js=",js," je=",je," joff=",joff
     &," diff_fn from row ",js-1," to ",je
     &,", diff_fe from row ",js," to ",je
     &,", diff_fbiso from row ",js," to ",je
#endif
c
#ifdef debug_tensor
      i    = 61
      k    = 3
      j    = 2
      jrow = j + joff
      if (mtrace .eq. 1 .and. (jrow .eq. 30 .or. jrow .eq. 31)) then
        diagonal = 1.0
        xyterm   = 0.0
        write (stdout,*) 'Showing tensor at jrow=',jrow,' i=',i,' k=',k
        write (stdout,'(a,3(1x,e14.7))')
     &            "K1=",diagonal, xyterm, K1(i,k,j,3)
        write (stdout,'(a,3(1x,e14.7))')
     &            "K2=",xyterm, diagonal, K2(i,k,j,3)
        write (stdout,'(a,3(1x,e14.7))')
     &            "K3=",K3(i,k,j,1), K3(i,k,j,2), K3(i,k,j,3)
c
c       construct the d/dx flux term = 1*dt/dx + 0*dt/dy + K13*dt/dz
c       on the western face of the T cell
c
        a1 = Ahisop*(t(i,k,j,1,taum1)-t(i-1,k,j,1,taum1))/
     &               (cst(jrow)*dxu(i-1))
        tbarxz = 0.25*(t(i-1,k-1,j,1,taum1) + t(i,k-1,j,1,taum1)
     &               + t(i-1,k,j,1,taum1) + t(i,k,j,1,taum1))
        tbarxzp = 0.25*(t(i-1,k,j,1,taum1) + t(i,k,j,1,taum1)
     &               + t(i-1,k+1,j,1,taum1) + t(i,k+1,j,1,taum1))
        a3 = Ahisop*K1(i-1,k,j,3)*(tbarxz - tbarxzp)/dzt(k)
        write (stdout,'(4(1x,a,e14.7))')
     &            " Ahisop*dT/dx(i-1)=",a1,"Ahisop*K13(i-1)*dT/dz=",a3
     &,           " sum(i-1)=",a1+a3,  " K13(i-1)=",K1(i-1,k,j,3)
c
c       construct the d/dy flux term = 0*dt/dx + 1*dt/dy + K23*dt/dz
c       on the southern face of the T cell
c
c
        b2 = Ahisop*csu(jrow-1)*(t(i,k,j,1,taum1)-t(i,k,j-1,1,taum1))/
     &               (dyu(jrow-1))
        tbaryz = 0.25*(t(i,k-1,j-1,1,taum1) + t(i,k-1,j,1,taum1)
     &               + t(i,k,j-1,1,taum1) + t(i,k,j,1,taum1))
        tbaryzp = 0.25*(t(i,k,j-1,1,taum1) + t(i,k,j,1,taum1)
     &               + t(i,k+1,j-1,1,taum1) + t(i,k+1,j,1,taum1))
        b3 = Ahisop*csu(jrow-1)*K2(i,k,j-1,3)*(tbaryz - tbaryzp)/dzt(k)
        write (stdout,'(4(1x,a,e14.7))')
     &            " Ahisop*dT/dy(j-1)=",b2,"Ahisop*K23(j-1)*dT/dz=",b3
     &,           " sum(j-1)=",b2+b3,  " K23(j-1)=",K2(i,k,j-1,3)
c
c       construct the d/dz flux term = K31*dt/dx + K32*dt/dy + K33*dt/dz
c       on the top face of the T cell
c
        dtbarxz = 0.25*(t(i,k,j,1,taum1) - t(i-1,k,j,1,taum1)
     &                + t(i+1,k,j,1,taum1) - t(i,k,j,1,taum1)
     &                + t(i,k-1,j,1,taum1) - t(i-1,k-1,j,1,taum1)
     &                + t(i+1,k-1,j,1,taum1) - t(i,k-1,j,1,taum1))
     &                 / (cst(jrow)*dxt(i))
        t31 = Ahisop*K3(i,k-1,j,1)*dtbarxz
c
        dtbaryz = 0.25*(t(i,k,j,1,taum1) - t(i,k,j-1,1,taum1)
     &                + t(i,k,j+1,1,taum1) - t(i,k,j,1,taum1)
     &                + t(i,k-1,j,1,taum1) - t(i,k-1,j-1,1,taum1)
     &                + t(i,k-1,j+1,1,taum1) - t(i,k-1,j,1,taum1))
     &                 / dyt(jrow)
        t32 = Ahisop*K3(i,k-1,j,2)*dtbaryz
c
        dtdz = (t(i,k-1,j,1,taum1) - t(i,k,j,1,taum1))/dzw(k-1)
	t33  = ahisop*K3(i,k-1,j,3)*dtdz
        write (stdout,'(3(1x,a,e14.7),/2(1x,a,e14.7),/3(1x,a,e14.7))')
     &            " Ahisop*K31(k-1)*dT/dx=",t31
     &,           " Ahisop*K32(k-1)*dT/dy=",t32
     &,           " Ahisop*K33(k-1)*dT/dz=",t33
     &,           " sum12(k-1)=",t31+t32," sum123(k-1)=",t31+t32+t33
     &,           " K31(k-1)=",K3(i,k-1,j,1)
     &,           " K32(k-1)=",K3(i,k-1,j,2)
     &,           " K33(k-1)=",K3(i,k-1,j,3)
c
c       fluxes calculated from MOM
c
        write (stdout,*) 'itt=',itt,' isoflux:n=',mtrace
     &,    ': showing calculated fluxes at jrow=', jrow,' i=',i,' k=',k
        write (stdout,'(a,3(1x,e20.13))')
     &            "diff_fe(i-1)=",diff_fe(i-1,k,j)
        write (stdout,'(a,3(1x,e20.13))')
     &            "diff_fn=",diff_fn(i,k,j-1),diff_fn(i,k,j)
        write (stdout,'(a,3(1x,e20.13))')
     &          "diff_fbiso=",diff_fbiso(i,k-1,j)
        write (stdout,'(a,3(1x,e20.13))')
     &            "tau-1 t(1,2,3)=",(t(i,k,j,mtrace,taum1),j=1,jmw)
      endif
#endif
# ifdef timing
      call toc ('isopycnal ', 'isoflux')
# endif
      return
      end

#else
      subroutine isop0
      return
      end
#endif

