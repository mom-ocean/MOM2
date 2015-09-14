#ifdef time_varying_sbc_data
      subroutine setatm (dtatm)
c
c-----------------------------------------------------------------------
c     set up the monthly climatological data which has been prepared
c     using the run scripts in PREP_DATA
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
# include "param.h"
# include "calendar.h"
# include "csbc.h"
# include "ctdbc.h"
# include "tmngr.h"
# include "switch.h"
      dimension xxu(imt), yyu(jmt), obc_data(imt,jmt)
c
c     names of sequential unformatted input data (S.B.C.)
c
      namein(1) = 'taux.mom'
      namein(2) = 'tauy.mom'
      namein(3) = 'sst.mom'
      namein(4) = 'sss.mom'
c
c     names of direct access data  (S.B.C.)
c
      nameout(1) = 'taux.ud'
      nameout(2) = 'tauy.ud'
      nameout(3) = 'sst.ud'
      nameout(4) = 'sss.ud'
c
       write (stdout,'(/a/)')
     &   ' ==> S.B.C. for MOM are comming from time varying S.B.C. data'
c
c-----------------------------------------------------------------------
c     transfer the monthly climiatolgical surface boundary conditions
c     from the sequential unformatted files (prepared by PREP_DATA)
c     to direct access files.
c     disk resource factor of 13 is for 12 months + 1 annual mean
c-----------------------------------------------------------------------
c
c     define how dataset is to be used: (periodic or not) and number
c     of records for each dataset
c
      do n=1,ntdbc
        ntdrec(n) = 12
	period(n) = .true.
        if (ntdrec(n) .gt. maxbcr) then
          write (stdout,'(/,1x,a)') 'Error: "ntdrec" exceeds "maxbcr"'
          stop '=>setatm'
        endif
      enddo
c
c     lrec needs room for "period" + jmt + character*32 = 6 words
c     assuming 4 characters per word
c
      lrec = imt+1+1+4
c
      totrec = 13     
c              13 is for 12 months + 1 annual mean
c
# ifdef cray_ymp
      write (opt_tdbc,'(a,i10,a,i10,a)')
     &' u d words =',lrec,' sds maxsize words='
     &,lrec*jmt*int(totrec),' novfl'
# else
      write (opt_tdbc,'(a,i10)') ' u d words =',lrec
# endif
c
      do n=1,ntdbc
        write (stdout,'(/a,a,a,1pg10.3,a)')
     & ' Sequential access disk resource for file ', namein(n),' = '
     &,  lrec*jmt*totrec*1.e-6,' (MW)'
        call getunit (nu, namein(n), 'u s r ieee')
        write (stdout,'(/a,a,a,1pg10.3,a)')
     & ' Direct access disk resource for file ', nameout(n),' = '
     &,  lrec*jmt*ntdrec(n)*1.e-6,' (MW)'
        call getunit (iotdbc(n), nameout(n), opt_tdbc)
c
c       transfer the input data to direct access files for each month
c
        write (stdout,'(/a,a)') 'Initializing ',nameout(n)
        cksum = 0.0
        do m=1,ntdrec(n)
          if (n .eq. 1 .or. n .eq. 2) then
c
c           read taux & tauy
c
	    read (nu)
	    read (nu) dstamp(m,n), aprec(m,n), iimt, jjmt, nn 
     &,            xxu, yyu, obc_data
          else
c
c           read sst & sss
c
	    read (nu)
	    read (nu) dstamp(m,n), aprec(m,n), iimt, jjmt, k, zz, nn
     &,               xxu, yyu, obc_data
	  endif
c
c         adjust feburary length when using data with leap year calendar
c
          if (.not. eqyear .and. nint(aprec(m,n)) .eq. 28) then
	    aprec(m,n) = aprec(m,n) + 0.2425
	    write (stdout,*)
     &     '=>Warning: adding 0.2425 days to feb for leap year calendar'
     &,    ' for ',nameout(n) 
	  endif
c
          cksum = cksum + checksum (obc_data, imt, jmt)
c
          call write_obc_data (iotdbc(n), m, dstamp(m,n), aprec(m,n)
     &,                        obc_data) 
	enddo
c
        print *,
     & ' checksum for monthly climatology = ', cksum, ' on ', nameout(n)
c
c       build time for start of S.B.C. data records from time stamp at
c       end of record. month is one less than that given by the stamp
c
	call rdstmp (dstamp(1,n), isbcyear, isbcmonth, isbcday
     &,              isbchour, isbcmin, isbcsec)
        isbcmon   = mod(isbcmonth+ntdrec(n)-2,ntdrec(n)) + 1
c
c       initialize time interpolation arrays and check record integrity
c
        call timeinterpi (ntdrec(n), dstamp(1,n), aprec(1,n), tdrec(1,n)
     &,                  isbcstart(n), period(n))
c
        call relunit (nu)
      enddo
c
# if defined minimize_sbc_memory
c
      do n=1,ntdbc 
c
c       express model time in days after start of S.B.C. by adding time
c       of I.C. to current model time then subtract time at start of
c       S.B.C.. Note that "itemptime" was allocated in settmngr and is
c       only needed as a temporary.
c
        call addtime (initial, imodeltime, itemptime)
        call subtime (itemptime, isbcstart(n), itemptime)
        daysbc(n) = realdays(itemptime)
#ifdef debug_atmos
        print *,' tp:n=',n,' daysbc(n)=',daysbc(n)
#endif
c
c       initialize the pointers (indices) to memory buffers and choose
c       interpolation method
c
        iprevm(n) = 1
        inextm(n) = 2
	method    = 3
c
c       find the disk indices and interpolation factor
c
        call timeinterp (daysbc, n, tdrec(1,n), aprec(1,n), ntdrec(n)
     &,      period(n), method, inextd(n), iprevd(n), wprev(n)
     &,      rdtdbc(n), inextm(n), iprevm(n))
      enddo
# else
      write (stdout,'(/a,a,a/)') '==> Note: to save memory when using'
     &,' "time_varying_sbc_data", enable "minimize_sbc_memory"'
     &,' which will use disk instead'
c
c-----------------------------------------------------------------------
c     determine time interpolation factors and disk record indices for 
c     interpolating time dependent surface boundary conditions to the
c     model time. the S.B.C. will be updated at the mid points of each 
c     month. For generality, each S.B.C. is treated separately.
c     This (for example) allows for a mixture of monthly heatflux with
c     daily winds.
c     this test case assumes, however, that all S.B.C are monthly
c
c     note: to understand how the time interpolator works, use
c           run_timeinterp.
c-----------------------------------------------------------------------
c
      do n=1,ntdbc 
c
c       express model time in days after start of S.B.C. by adding time
c       of I.C. to current model time then subtract time at start of
c       S.B.C.. Note that "itemptime" was allocated in settmngr and is
c       only needed as a temporary.
c
        call addtime (initial, imodeltime, itemptime)
        call subtime (itemptime, isbcstart(n), itemptime)
        daysbc(n) = realdays(itemptime)
#ifdef debug_atmos
        print *,' tp:n=',n,' daysbc(n)=',daysbc(n)
#endif
c
c       initialize the pointers (indices) to memory buffers and choose
c       interpolation method
c
        iprevm(n) = 1
        inextm(n) = 2
	method    = 3
c
c       find the disk indices and interpolation factor
c
        call timeinterp (daysbc, n, tdrec(1,n), aprec(1,n), ntdrec(n)
     &,      period(n), method, inextd(n), iprevd(n), wprev(n)
     &,      rdtdbc(n), inextm(n), iprevm(n))
c
c       read data straddling (just before and after) current model time
c
        js = 1
	je = jmt
	joff = 0
        call read_obc_data (iotdbc(n), iprevd(n), rstamp1
     &,                     obc(1,1,n,iprevm(n)), js, je, joff)
        call read_obc_data (iotdbc(n), inextd(n), rstamp2
     &,                     obc(1,1,n,inextm(n)), js, je, joff)
c
        write (stdout,8800) n, iprevd(n), rstamp1, iprevm(n) 
     &,                        inextd(n), rstamp2, inextm(n)
     &,                        nameout(n), stamp
        write (stdout,'(/24x,a12,a,i2,a,i2/)')
     &  nameout(n),' is dataset index ',n
     &, ' for time interpolation using method #', method
c
        call relunit (iotdbc(n))
      enddo
8800  format (/1x,'==> Reading S.B.C. #',i2,' from rec ',i6, ' at ', a
     &,' into index ',i3,/24x,' and rec ',i6,' at ', a, ' into index '
     &,i3,/,24x, ' for ',a12, ' at MOM time ',a) 
# endif
c
c-----------------------------------------------------------------------
c     zero all fluxes in case a zero flux S.B.C. is desired
c     "m" relates how S.B.C. "n" is ordered within "sbcocn"
c-----------------------------------------------------------------------
c
      do n=1,numosbc
        m = mapsbc(n)
        do jrow=1,jmt
          do i=1,imt
            sbcocn(i,jrow,m) = 0.0
          enddo
        enddo
      enddo
c
# ifdef coupled
c
c     verify that "coupled" is not enabled
c
      write (stdout,*) '==> Error: "time_varying_sbc_data" and '
     &, '"coupled" are incompatable'
      stop '=>setatm'
# else
c
c-----------------------------------------------------------------------
c     since this is not an atmospheric model, the atmospheric 
c     time step "dtatm" does not have to be set. mom.F will set it.
c-----------------------------------------------------------------------
c
# endif
# ifdef simple_sbc
      write (stdout,'(/,a,a,/)')
     &'==> Error: "time_varying_sbc_data" and '
     &,' "simple_sbc" are incompatable'
      stop '=>setatm'
# endif
c
c-----------------------------------------------------------------------
c     initialize euler2 because it is tested within atmos but not set
c     until the first call to subroutine mom
c-----------------------------------------------------------------------
c
      euler2 = .false.
      return
      end

      subroutine write_obc_data (iou, m, stamp, period, obc)
# include "size.h"
# include "stdunits.h"
      dimension obc(imt,jmt), buf(imt)
      character*(*) stamp
      do j=1,jmt
        nrec = (m-1)*jmt + j
        do i=1,imt
	  buf(i) = obc(i,j)
        enddo
        write (iou, rec=nrec) stamp, period, j, buf
      enddo
      write (stdout,'(a,i2,1x,a,a,g14.7,a)') ' month=',m, stamp
     &,  '.  averaged over preceeding ', period, 'days'
      return
      end


      subroutine read_obc_data (iou, m, stamp, obc, js, je, joff)
# include "size.h"
      dimension obc(imt,js:je), buf(imt)
      character*(*) stamp
      do j=js,je
        jrow = j + joff
        nrec = (m-1)*jmt + jrow
        read (iou, rec=nrec) stamp, period, jj, buf
        do i=1,imt
	  obc(i,j) = buf(i)
        enddo
      enddo
      return
      end


#else
      subroutine nodatasetatm
      return
      end
#endif
