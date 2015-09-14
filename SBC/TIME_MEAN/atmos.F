#ifdef time_mean_sbc_data
      subroutine atmos
c
c-----------------------------------------------------------------------
c     this "atmos" model is an EXAMPLE of how to supply MOM with surface
c     boundary conditions which are a function of latitude and
c     longitude but NOT of time.
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
# include "param.h"
# include "csbc.h"
      character*32 stamp
      dimension xxt(imt), yyt(jmt), xxu(imt), yyu(jmt)
      common /catmos/ sst(imt,jmt), sss(imt,jmt)
c
c
c-----------------------------------------------------------------------
c     Work on all rows assuming S.B.C. arrays are (imt,jmt)
c-----------------------------------------------------------------------
c
      joff = 0
      js   = 1
      je   = jmt
c
c-----------------------------------------------------------------------
c     set indices which relate how S.B.C. "n" is arranged within arrays
c     eg: m = mapsbc(n) where sbcocn(,,m) is S.B.C. "n"
c-----------------------------------------------------------------------
c
      itaux = mapsbc(1)
      itauy = mapsbc(2)
      ihflx = mapsbc(3)
      ipme  = mapsbc(4)
c     ishort= mapsbc(5)
      isst  = mapsbc(6)
      isss  = mapsbc(7)
c
      if (afirst) then
c
c-----------------------------------------------------------------------
c       on the first call to "atmos", read the time mean data "taux" 
c       and "tauy" along with the time mean SST and SSS for restoring
c       surface tracers. It is assumed that this data has been prepared
c       with the run scripts in PREP_DATA
c-----------------------------------------------------------------------
c
        call getunit (nutaux, 'taux.mom', 'u s r ieee')
        call getunit (nutauy, 'tauy.mom', 'u s r ieee')
        call getunit (nusst, 'sst.mom', 'u s r ieee')
        call getunit (nusss, 'sss.mom', 'u s r ieee')
c
c-----------------------------------------------------------------------
c       skip the monthly values and read the annual means
c-----------------------------------------------------------------------
c
        do n=1,12
          read (nutaux)
          read (nutaux)
c
          read (nutauy)
          read (nutauy)
c
	  read (nusst)
	  read (nusst)
c
	  read (nusss)
	  read (nusss)
        enddo
c
c       annual mean surface wind stress components (13th record)
c
        read (nutaux)
        read (nutaux) stamp, dpm, iimt, jjmt, nn, xxu, yyu
     &,             ((sbcocn(i,j,itaux),i=1,imt),j=1,jmt)
        write (stdout,'(/a,g14.7,a/)')
     & ' initilized annual mean taux. annual period=',dpm, 'days'
        read (nutauy)
        read (nutauy) stamp, dpm, iimt, ijmt, nn, xxu, yyu
     &,             ((sbcocn(i,j,itauy),i=1,imt),j=1,jmt)
        write (stdout,'(/a,g14.7,a/)')
     & ' initilized annual mean tauy. annual period=',dpm, 'days'
c
c       annual mean SST and SSS  (13th record)
c
        read (nusst)
        read (nusst) stamp, dpm, iimt, ijmt, k, zp, nn, xxt, yyt, sst
        write (stdout,'(/a,g14.7,a/)')
     & ' initilized annual mean SST. annual period=',dpm, 'days'
        read (nusss)
        read (nusss) stamp, dpm, iimt, ijmt, k, zp, nn, xxt, yyt, sss
        write (stdout,'(/a,g14.7,a/)')
     & ' initilized annual mean SSS. annual period=',dpm, 'days'
c
        call print_checksum (sbcocn(1,1,itaux), imt, jmt
     &,                    ' checksum for annual mean taux =')
        call print_checksum (sbcocn(1,1,itauy), imt, jmt
     &,                    ' checksum for annual mean tauy =')
        call print_checksum (sst, imt, jmt
     &,                    ' checksum for annual mean sst =')
        call print_checksum (sss, imt, jmt
     &,                    ' checksum for annual mean sss =')
c
        call relunit (nutaux)
        call relunit (nutauy)
        call relunit (nusst)
        call relunit (nusss)
      endif
c
c----------------------------------------------------------------------
c     Taux and Tauy don`t change with time, but heatflux and saltflux
c     may change if restoring SST and SSS to time mean data
c----------------------------------------------------------------------
c
      do j=js,je
        jrow = j + joff
# ifdef restorst
c
c       restore surface tracers to prescribed data using a newtonian
c       damping term with a restoring time scale of "dampts" (days)
c       (converted to a surface tracer flux by factor "damp").
c
c
	if (isst .ne. 0 .and. ihflx .ne. 0) then
          damp = dampdz(1)*secday/dampts(1)
	  do i=2,imtm1
	    sbcocn(i,jrow,ihflx) = damp*
     &                             (sst(i,jrow) - sbcocn(i,jrow,isst))
          enddo
        endif
	if (isss .ne. 0 .and. ipme .ne. 0) then
          damp = dampdz(2)*secday/dampts(2)
	  do i=2,imtm1
	    sbcocn(i,jrow,ipme) = damp*
     &                            (sss(i,jrow) - sbcocn(i,jrow,isss))
          enddo
        endif
# else
c
c       supply zero heat flux and salt flux
c
        if (ihflx .ne. 0) then
          do i=1,imt
            sbcocn(i,jrow,ihflx) = c0
          enddo
        endif
        if (ipme .ne. 0) then
          do i=1,imt
            sbcocn(i,jrow,ipme) = c0
          enddo
        endif
# endif
      enddo
c
      return
      end
#else
      subroutine nodataatmos
      return
      end
#endif
