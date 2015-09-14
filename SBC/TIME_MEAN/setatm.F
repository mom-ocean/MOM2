#ifdef time_mean_sbc_data
      subroutine setatm (dtatm)
c
c-----------------------------------------------------------------------
c     set up all quantities that need to be defined only once for this
c     simple atmosphere.
c     
c     coded by:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c-----------------------------------------------------------------------
c
# include "stdunits.h"
c
c     The atmospheric time step "dtatm" is not needed since this is
c     not an atmospheric model. it will be defined in mom.F 
c
      write (stdout,'(/,a,/)')
     &'==> S.B.C. for MOM are comming from annual mean S.B.C. data'
c
# ifdef simple_sbc
      write (stdout,'(/,a,a,/)')
     &'==> Error: "time_mean_sbc_data" and '
     &,' "simple_sbc" are incompatable'
      stop '=>setatm'
# endif
# ifdef coupled
c
c     verify that "coupled" is not enabled.
c
      write (stdout,*) '==> Error: "time_mean_sbc_data" and '
     &, '"coupled" are incompatable'
       stop '=>setatm'
# endif
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
