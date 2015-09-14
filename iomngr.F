c=======================================================================
c
c     iomngr module:  see subroutine print_iomngr_description below
c
c     The following iomngr define options may be selected:
c
c     stand-alone test of iomngr
c #define test_iomngr
c     print description of iomngr features and their use
c #define print_iomngr_description
c     generate "io_history" file showing i/o file use
c #define io_history
c     print record of calls to getunit and relunit, etc.
c #define verbose_iomngr
c     [cray_ymp only:] print "assign" system calls made by getunit
c #define show_assigns
c     print cray_ymp debugging information about file attributes
c #define debug_iomngr
c     print more detailed debugging information
c #define debug_iomngr2

c=======================================================================

#if defined verbose_iomngr && defined cray_ymp
# define show_assigns
#endif

#ifdef print_iomngr_description
      subroutine io_description

c=======================================================================
c
c     author c. h. goldberg          e-mail: chg@gfdl.gov
c
c                                    revised May 17, 1995
c     
c     Please inform him of either problems with or desired extensions to
c     the functionality offered by this package.
c
c     Additions are very easy to make, but I would like to maintain some
c     consistency in the syntax of specifying options.
c
c=======================================================================
c

#include "iomngr.h"

      write (stdout, '(a)') ' '
     &,'      I/O Unit Manager enabled'
     &,' '
     &,'  To use i/o unit manager, first initialize with'
     &,'      call ioinit'
     &,'  Then files should be opened with'
     &,'      call getunit (iounit, filename, optionlist)'
     &,'  where "iounit" is an integer variable that will be assigned'
     &,'  a unit number not currently in use, "filename" is a'
     &,'  character constant or variable giving the [unix] name of the'
     &,'  file to be opened, and "optionlist" is a character constant'
     &,'  or variable choosing options from the following list:'
     &,' '
     &,'      s   sequential'
     &,'      d   direct [recl] [:] [=] size_spec'
     &,'      f   formatted'
     &,'      u   unformatted'
     &,'      r   rewind'
     &,'      a   append'
     &,'          ieee'
     &,'      ro  read [only]'
     &,'      wo  write [only]'
     &,'      rw  readwrite' 
#ifdef cray_ymp
     &,'      w   word [buffer] [:] size_spec   <used with wread>'
     &,'          sdsalloc [max] [:] size_spec  <used with ssread>'
     &,' '
     &,'       Cray ymp memory and ssd options'
     &,' '
     &,'          mr [max] [:] size_spec'
     &,'          sds [max] [:] size_spec'
     &,'          ovfl'
     &,'          novfl'
     &,'    scr   scratch'
     &,'          save'
#endif
     &,' '
     &,'       Size_spec options'
     &,' '
#ifdef cray_ymp
     &,'          blocks [=] digits'
     &,'          words [=] digits        <64-bit words>'
#else
# ifdef sgi
     &,'          words [=] digits        <32-bit words>'
# else
     &,'          words [=] digits        <32-bit words>'
# endif
#endif
     &,'          bytes [=] digits'
     &,'          digits'
     &,' '
     &,'       all size_spec options except digits are automatically'
     &,'       converted to proper machine-dependent units.'
     &,' '
     &,'       Blocking options'
     &,' '
     &,'          fortran'
     &,'          f77              <synonym for ''fortran''>'
     &,'          unblocked'
     &,'          stream           <synonym for ''unblocked''>'
     &,'          null             <synonym for ''unblocked''>'
#ifdef cray_ymp
     &,'          cos'
     &,'          cache'
#endif
     &,' '
     &,'  Either upper or lowercase may be used, and the short forms'
     &,'  in the left column may be mixed with the more explanatory'
     &,'  long forms.  Spaces, commas, colons, and equal signs may'
     &,'  appear between options.' 
      write (stdout, '(a)') ' '
     &,'    For example:' 
     &,'       call getunit (idiag, ''diagnostics'', '
     &,'      &           ''unformatted, sequential, ieee, append'')'
     &,'       call getunit (islab, ''slabs'', ''ud1000'')'
     &,'       call getunit (init, ''init.data'', ''fsr'')'
     &,'       write (string,''(a,i4)'') '
     &,'      &     ''unformatted, direct, recl words='', imt'
     &,'       call getunit (islice, ''slice'', string)'
     &,' '
     &,'  Files should be closed with'
     &,'      call relunit (iounit)'
     &,'  to recycle their unit numbers.'
      write (stdout, '(a)') ' '
     &,'  If io_history is defined, a file "iohistory" is created in'
     &,'  which all calls to getunit and relunit are recorded.'
     &,' '
     &,'  The subroutine call'
     &,'      call showunits'
     &,'  lists all i/o units currently in use and their filenames.'
     &,'  '
     &,'  The subroutine call'
     &,'      call showfiles'
     &,'  lists in chronological order all files that have been opened'
     &,'  '
     &,'  The subroutine call'
     &,'      call release_all'
     &,'  releases and closes all files that are still open.  To avoid'
     &,'  all potential conflicts, one must always call release_all'
     &,'  before executing a section of code that may open units'
     &,'  without using the io_manager.  The io_manager is resumed'
     &,'  using the subroutine call'
     &,'      call iomngr_resume'
#ifdef cray_ymp
     &,'  '
     &,'  On the cray ymp, appends to pre-existing ieee unformatted'
     &,'  sequential files will not be completed unless this call is'
     &,'  made.  Without this call, such files will have only the new' 
     &,'  appends, and the pre-existing file will be renamed with the'
     &,'  suffix ''.old''.'
#endif
     &,'  '

      return
      end
#endif

#ifdef test_iomngr

      program driver
c
c=======================================================================
c
c     This driver exercises various features of the io-manager.
c     
c
c     To run the driver, look at and execute the "run_iomngr" script
c
c     author:         c. h. goldberg    e-mail=> chg@gfdl.gov
c=======================================================================
c
#include "stdunits.h"
      character*160 options
      character*15 reclen, direct_file_size
      character*80 iotext
      parameter (imt=100,jmt=50)
      dimension a(imt,jmt), b(imt,jmt)
      logical init, restrt
      namelist /contrl/ init, days, ntasks, restrt
      namelist /tsteps/ dtts, dtuv, dtsf
c
      call ioinit
c
      init = .false.
      days = 2.0
      ntasks = 5
      restrt = .true.
      dtts = 1.0
      dtuv = 2.0
      dtsf = 3.0
c
      write (stdout, contrl)
      write (stdout, tsteps)
c
#ifdef cray_ymp
c     testing formatted sequential append memory resident file
c
      print '(/,a,/)', 'Testing option: fsa mr save'
      call getunit (iou, 'iotest.doc', 'fsa mr save')
      write (iou,contrl)
      call relunit (iou)
      call getunit (iou, 'iotest.doc', 
     &              'formatted sequential append mr save')
      write (iou,tsteps)
      write (iou,'(1x,a,i10)') 'imt=',imt,' jmt=',jmt
      call relunit (iou)
#endif
      print '(/,a,/)', 'Testing option: formatted sequential append'
      call getunit (iou, 'iotest.doc', 'fsa ')
      write (iou,contrl)
      call relunit (iou)
      call getunit (iou, 'iotest.doc', 
     &              'formatted sequential append')
      write (iou,tsteps)
      write (iou,'(1x,a,i10)') 'imt=',imt,' jmt=',jmt
      call relunit (iou)
c
c     testing unformatted sequential append IEEE byte stream file
c
      print '(/,a,/)', 
     &  'Testing unformatted sequential append IEEE byte stream file'
      call getunit (iou, 'iotest.udoc', 'usa ieee stream')
      call showunits
      iotext = 'read(io) imt, jmt'
      write (iou) iotext
      write (iou) imt, jmt
      call relunit (iou)
      call showunits
c
      call getunit (iou, 'iotest.udoc', 'usr ieee stream')
      call showunits
      read (iou) iotext
      read (iou) iimt, jjmt
      print '(a)',iotext
      print '(a,i4,a,i4)','imt=',iimt,' jmt=',jjmt
      call relunit (iou)
      call showunits
c
c     testing unformatted sequential append IEEE f77 blocked file
c
      print '(/,a,/)', 
     &  'Testing unformatted sequential append IEEE f77 blocked file'
      call getunit (iou, 'iotest.vdoc', 'usa ieee ')
      call showunits
      iotext = '(1) read(io) imt, jmt'
      write (iou) iotext
      write (iou) imt, jmt
      call relunit (iou)
      call showunits
c
      call getunit (iou, 'iotest.vdoc', 'usa ieee ')
      call showunits
      iotext = '(2) read(io) imt, jmt'
      write (iou) iotext
      write (iou) imt, jmt
      call relunit (iou)
      call showunits
c
c     testing isolation of io_manager allocated files from code that
c     has unit number hard wired.
c
      call release_all
      call showunits
c
c     open three hard-wired units without using the io_manager
c
      print '(/,a,/)', 'writing formatted to hard-wired units 1 to 3.'
c
      open (unit=1, file='iotest.fort1')
      write (1, '(a)') 'formatted stuff to unit 1'
c
      open (unit=2, file='iotest.fort2')
      write (2, '(a)') 'formatted stuff to unit 2'
c
      open (unit=3, file='iotest.fort3')
      write (3, '(a)') 'formatted stuff to unit 3'
c
c     close unit 2, but neglect to close units 1 and 3
c
      print '(/,a,/)', 'closing unit 2, but not units 1 and 3'
      close (2)
c
c     resume io_manager
c
      call iomngr_resume
c
      print '(/,a,/)', 'resuming io_manager tests: usa ieee'
      call getunit (iou, 'iotest.vdoc', 'usa ieee ')
      call showunits
      iotext = '(3) read(io) imt, jmt'
      write (iou) iotext
      write (iou) imt, jmt
      call relunit (iou)
      call showunits
c
      call getunit (iou, 'iotest.vdoc', 'usr ieee ')
      call showunits
      print '(/,a,/)', 'reading back file iotest.vdoc'
      do i=1,1000
        read (iou, end=99) iotext
        read (iou, end=99) iimt, jjmt
        print '(a,i3)', '  RECORD', i
        print '(a)',iotext
        print '(tr2,a,i4,a,i4)','imt=',iimt,' jmt=',jjmt
      end do
   99 continue
      call relunit (iou)
      call showunits
c
c     testing unformatted sequential append IEEE f77 blocked file
c     buffer size is chosen for efficiency
c
      print '(/,a,/)', 'test1: usa ieee buffer=393216'
      call getunit (iou, 'iotest.buf', 'usa ieee buffer=393216')
      write (iou) imt, jmt
      call showunits
      call relunit (iou)
      call showunits
c
      call getunit (iou, 'iotest.buf', 'usa ieee buffer=393216')
      call showunits
      write (iou) imt, jmt
      call relunit (iou)
      call showunits
c
      call getunit (iou, 'iotest.buf', 'usr ieee buffer=393216')
      call showunits
      read (iou) i, j
      read (iou) k, l
      print *,'imt=',i, 'jmt=',j,'imt=',k, 'jmt=',l
      call relunit (iou)
c 
      call showunits     
      call showfiles
c
c     direct access
c
c     note: there are imt=100 values in each record.  getunit automatically
c           converts 'direct recl:words=100' to the correct record length
c           for each machine: 800 bytes on cray_ymp and 100 words on sgi.
c
      print '(/,a,/)', 'test2: direct access'
      do j=1,jmt
        do i=1,imt
          a(i,j) = 100*i + j
        enddo
      enddo
c
      write (reclen,'(i10)') imt
      write (direct_file_size,'(i10)') imt*jmt
#ifdef cray_ymp
      options = 'unformatted direct recl:words=' // reclen
     &       // 'sds maxsize words=' // direct_file_size //' novfl'
#else
      options = 'unformatted direct recl:words=' // reclen
#endif
      call getunit (iou, 'iotest.direct', options)
      do j=1,jmt
        write (iou,rec=j) (a(i,j),i=1,imt)
	if (j .eq. jmt/2) then
	  call relunit (iou)
          call getunit (iou, 'iotest.direct', options)
	endif
      enddo
      do j=jmt,1,-1
        read (iou,rec=j) (b(i,j),i=1,imt)
      enddo
      kerror = 0
      do j=1,jmt
        do i=1,imt
	  if (a(i,j) .ne. b(i,j)) then
	    kerror = kerror + 1
	  endif
	enddo
      enddo
      print '(/a,i2)','Direct access: number of reread errors =',kerror
c      
      call showunits
      call showfiles
      call release_all
      stop
      end



      function LEN_TRIM (string)
c
c     defines length of "string" ignoring trailing blanks
c
      character * (*) string
      integer LEN_TRIM, i
c
      LEN_TRIM = 0
      do i = LEN (string), 1, -1
        if (string(i:i) .NE. ' ') then
          LEN_TRIM = i
          goto 999
        end if
      end do
999   continue
c      
      return
      end

#endif

      subroutine ioinit
c
c=======================================================================
c     Ioinit initializes logical array "inuse" to show which i/o units
c     are available to the subroutine getunitnumber.  It also 
c     initializes the array fname of file names associated with each
c     i/o unit, and opens a file "iohistory", if enabled.
c
c     author c. h. goldberg          e-mail: chg@gfdl.gov
c=======================================================================
c
#include "iomngr.h"
      parameter (last_reserved_unit=0)
c
#if !defined cray_ymp && !defined sgi
      print '(a)'
     &, 'WARNING:  This i/o manager has been tested extensively only'
     &, '          on the Cray YMP and Silicon Graphics Personal Iris'
     &, '          and Indigo computers.  If you are using one of these'
     &, '          you must enable the -Dcray_ymp or -Dsgi define'
     &, '          options for proper execution.'
     &, ' '
     &, '          If you are not using one of these computers, there'
     &, '          are generic choices made in the i/o manager. However'
     &, '          certain aspects of Fortran I/O such as size units'
     &, '          for record lengths in direct access are not'
     &, '          standardized, and others such as INQUIRE and APPEND'
     &, '          do not become standard until Fortran 90, although'
     &, '          they appear as possibly not standard-conforming'
     &, '          extensions in many Fortran compilers.'
     &, ' '
     &, '          Also, there are options in MOM to use high'
     &, '          performance features of vendor-specific I/O on'
     &, '          computers available to the developers.'
     &, ' '
     &, '          If there are no I/O problems running MOM on your'
     &, '          computer, please add its name to the cpp test'
     &, '          that surpresses this message and error stop.'  
     &, '          We suggest running the stand-alone driver for the'  
     &, '          i/o manager and other test programs designed to'  
     &, '          exercise the specific i/o options you are using.'  
     &, '          Pay particular attention to record and file sizes'  
     &, '          and to record blocking on your system.' 
      stop '==>iomngr' 
#endif
      do i=0,maxfilenames
        hide_file(i) = .false.
        iunit(i) = 0
        fname(i) = ' '
      end do
c
      fname(0) = ' '
      fname(1) = 'reserved'
      nfiles = 1
c
c     reserve units 1..last_reserved_unit (making them not available)
c             units last_reserved_unit+1..maxunit=99 are available
c
c     unnecessary assignment to avoid WARNING on sgi f77 compiler
      last_reserved = last_reserved_unit
c
      do i=1,last_reserved
        inuse(i) = .true.
        ifile(i) = 1
        scratch_file(i) = .false.
        unix_name(i) = .true.
      end do
      do i=last_reserved_unit+1,maxunit
        inuse(i) = .false.
        ifile(i) = 0
        scratch_file(i) = .false.
        unix_name(i) = .true.
      end do
c
c     reserve units 5 and 6, whether or not they are stdin/stdout
c
      inuse(5) = .true.
      inuse(6) = .true.
      ifile(5) = 1
      ifile(6) = 1
c
c     reserve standard i/o files:  stdin, stdout, stderr
c
      call link_unit (stdin, 'stdin')
      call link_unit (stdout, 'stdout')
      if (stderr .ne. stdout) then
        call link_unit (stderr, 'stderr')
      else
        fname(ifile(stdout)) = 'stdout/stderr'
      end if
c
#ifdef print_iomngr_description
      call io_description
#endif
#ifdef io_history
      call getunit (iohist, 'iohistory', 'sequential formatted rewind')
#endif
#ifdef verbose_iomngr
      call showunits
#endif
#ifdef cray_ymp
c     A blanket "assign -R" in initialization is a mistake because it 
c     wipes out assignments in the run script, such as assign to sds or
c     assign an actual file name.  There are few i/o/ options not 
c     handled by getunit, but still the backwards compatibility with
c     shell script assigns is maintained by using incremental assigns
c     in this module.
#endif
c
      return
      end


      subroutine iomngr_resume
c
c=======================================================================
c     iomngr_resume resumes normal operation of the iomngr after a
c     call release_all has closed all files cleanly to permit executing
c     a section of code that does not use the iomngr to allocate
c     units.  It reopens the file "iohistory", if enabled.
c
c     author c. h. goldberg          e-mail: chg@gfdl.gov
c=======================================================================
c
#include "iomngr.h"
c
#ifdef verbose_iomngr
      print '(/,a,/)','Resuming io_manager'
#endif
#ifdef io_history
      call getunit (iohist, 'iohistory', 'sequential formatted append')
#endif
#ifdef verbose_iomngr
      call showunits
#endif
c
      return
      end


      subroutine getunitnumber (iounit)
c
c=======================================================================
c     getunitnumber assigns a Fortran i/o unit number that does not
c     conflict with any i/o unit currently in use.
c
c     comments:
c
c     see relunit which releases a Fortran i/o unit number.
c
c     author c. h. goldberg          e-mail: chg@gfdl.gov
c=======================================================================
c
#include "iomngr.h"
      character * 2 aiounit
      logical       exists,named
      character * 64 name
      character * 64 string

      do i = 1,maxunit
#ifdef debug_iomngr2
        print '(a,i2,a,l1)','in getunitnumber,inuse(',i,') = ',inuse(i)
#endif
        if (.NOT. inuse(i)) then
          inquire(i,EXIST=exists,NAMED=named,NAME=name)
#ifdef debug_iomngr2
          write (stdout,*) 'unit ',i, ' exists = ', exists
          write (stdout,*) 'unit ',i, ' named = ', named
          if (named) write (stdout,*) 'NAME=',name
#endif
c
c         A unit opened, even without the knowledge of Getunit will
c         still be NAMED.  Getunit skips over it.
c
          if (.NOT. named) then
c
c           Cray ymp: (5/94) (Unicos 7.x, cf77)
c           A unit given a file name in an "assign -a filename"
c           statement in the Unix script, but not yet opened by a
c           Fortran open statement will inquire as ".NOT. NAMED",
c           but in an attempted open of the unit with a different name,
c           IOSTAT will return error value 1062.  Such units should
c           not be used.
c           Note:  this test does no harm in sgi Fortran (5/94).
c                  It has not been tested with sgi equivalent of
c                  "assign -a filename", if there is one.
c
            open (i,FILE='dummy',IOSTAT=iostat)
c
#ifdef debug_iomngr2
            write (stdout,*) 'iostat=',iostat
            inquire(i,EXIST=exists,NAMED=named,NAME=name)
            write (stdout,*) 'iostat=',iostat
            write (stdout,*) 'NAME=',name
#endif
            if (iostat .EQ. 0) then
              close (i)
              goto 101
            end if
          end if
        end if
      end do
c
c     "Normal" termination of the do loop indicates all units are in use
c
      write (stdout, '(a4,)')
     &    'ERROR:  No i/o units available -- Execution terminated'
      call showunits
      stop '==>Getunit'
c
  101 continue
      iounit = i
      inuse(i) = .true.
c
#ifdef verbose_iomngr
c     write (stdout,*)
c     write (stdout,'(a,i3)') 'acquiring i/o unit',i
#endif

      return
      end


      subroutine getunit (iounit, filename, optionlist)
c
c=======================================================================
c     getunit gets an unused Fortran unit number, assigns it to the
c     dummy argument iounit, and opens it with FILE=filename and options 
c     given in optionlist.  Enable cpp directive "print_io_description"
c     to see the list of options.
c
c     author c. h. goldberg          e-mail: chg@gfdl.gov
c=======================================================================
c
#include "iomngr.h"
c
      character * (*) filename, optionlist
      character * 160 options
      character * 2 aiounit
      character * 160 string
      character * 15 ssequential, sdirect
      character * 15 sformatted, sunformatted
      character * 15 srewind, sappend
      character * 15 sread, swrite, sreadwrite, sunknown, sold, snew
      character * 15 sstatus, saccess, sform, sposition, saction
      character * 15 sieee, srec_len, sdummy, sbyte, sword, sblock
      character * 15 sbuffer, sblocking, psblocking, snumeric, snum
      character * 15 ssave, sovfl
      character * 15 sffio
      character * 40 slayer
      dimension slayer(6)
      logical scan_option, scan_number, scan_size, exists, named
      integer sdsalloc, errors
      logical skip_open

      data ssequential  /'SEQUENTIAL'/
      data sdirect      /'DIRECT'/
      data sformatted   /'FORMATTED'/
      data sunformatted /'UNFORMATTED'/
      data srewind      /'REWIND'/
      data sappend      /'APPEND'/
      data sread        /'READ'/
      data swrite       /'WRITE'/
      data sreadwrite   /'READWRITE'/
      data sunknown     /'UNKNOWN'/
      data sold         /'OLD'/
      data snew         /'NEW'/
c
c     set default file attributes
c
      saccess     = ssequential
      sform       = ' '
      sposition   = srewind
      saction     = sreadwrite
      sstatus     = sunknown
      sblocking   = ' '
      psblocking  = ' '
      sieee       = ' '
      snumeric    = ' '
      ssave       = 'save'
      sovfl       = 'novfl'
      irec_len    = 0
      iblock      = 0
      iword       = 0
      ibyte       = 0
      inum        = 0
      ibuffer     = 0
      maxsize     = 0
      layer       = 0
      sffio       = ' '
c
c     convert options list to lowercase
c
      options = optionlist
      call tolower (options)
#ifdef show_assigns
      print '(a,a)', 'options = ',options(1:len_trim(options))
#endif
c
c     decode options:
c
      i = 1
      errors = 0
100   continue
      if (i .gt. len_trim(options)) go to 101
c
c       scan options by full names
c
        if     (scan_option(options, i, 
     &          'sequential', saccess, ssequential)) then
        else if(scan_option(options, i, 
     &          'direct', saccess, sdirect)) then
          if(scan_option(options, i, 'recl', sdummy, sdummy)) then
          end if
          if(scan_size(options, i, irec_len, iblock, iword, ibyte)) then
          end if
	  call setreclen (irec_len,srec_len,sform,iblock,iword,ibyte)
	  if (irec_len .eq. 0) then
            print '(a,a,/,a,a)',
     &         '=>ERROR in getunit options list:',
     &         ' ''direct'' requires a record length',
     &         ' options=',options
            errors = errors + 1
          end if
        else if(scan_option(options, i, 'recl', sdummy, sdummy)) then
          if(scan_size(options, i, irec_len, iblock, iword, ibyte)) then
          end if
          call setreclen (irec_len,srec_len,sform,iblock,iword,ibyte)
	  if (irec_len .eq. 0) then
            print '(a,a,/,a,a)',
     &         '=>ERROR in getunit options list:',
     &         ' ''recl='' requires a record length',
     &         ' options=',options
            errors = errors + 1
          end if
        else if(scan_option(options, i, 
     &          'formatted', sform, sformatted)) then
        else if(scan_option(options, i, 
     &          'unformatted', sform, sunformatted)) then
        else if(scan_option(options, i, 
     &          'readwrite', saction, sreadwrite)) then
        else if(scan_option(options, i, 
     &          'read', saction, sread)) then
	  if (scan_option(options, i, 'only', sdummy, sdummy)) then
	  end if
          print '(a,a)', 'WARNING:  ACTION=''read'' not supported.  ',  
     &                            'ACTION=''READWRITE'' is used.'
        else if(scan_option(options, i, 
     &          'write', saction, swrite)) then
	  if (scan_option(options, i, 'only', sdummy, sdummy)) then
	  end if
          print '(a,a)', 'WARNING:  ACTION=''write'' not supported.  ',  
     &                            'ACTION=''READWRITE'' is used.'
        else if(scan_option(options, i, 
     &          'rewind', sposition, srewind)) then
        else if(scan_option(options, i, 
     &          'append', sposition, sappend)) then
        else if(scan_option(options, i, 'ieee', sieee, 'ieee')) then
c
c       scan for blocking options
c
        else if(scan_option(options, i, 'fortran', 
     &          sblocking, 'f77')) then
        else if(scan_option(options, i, 'f77', 
     &          sblocking, 'f77')) then
        else if(scan_option(options, i, 'stream', 
     &          sblocking, 'null')) then
        else if(scan_option(options, i, 'unblocked', 
     &          sblocking, 'null')) then
        else if(scan_option(options, i, 'null', 
     &          sblocking, 'null')) then
#ifdef cray_ymp
        else if(scan_option(options, i, 'cos', 
     &          sblocking, 'cos')) then
	  call new_layer(layer, sblocking, psblocking, slayer, sffio, 
     &                   inum, ibuffer, maxsize, iblock, iword, ibyte,
     &                   ssave, sovfl)
          if(scan_option(options, i, 'buffer', sdummy, sdummy)) then
	  end if
          if (scan_size(options, i, ibuffer, iblock, iword, ibyte)) then
	    call set2block (ibuffer, iblock, iword, ibyte)
            inum = ibuffer
          end if
#endif
c
c       the 'buffer' option is parsed, but ignored on many systems
c       f77 blocking is assumed
c
        else if(scan_option(options, i, 'buffer', sdummy, sdummy)) then
c
c         force new_layer to recognize f77 blocking
c
	  sblocking = 'f77'
          if (scan_size(options, i, ibuffer, iblock, iword, ibyte)) then
          end if
c         why assume that buffers are in bytes?   f77 buffer units=bytes.
	  call set2byte (ibuffer, iblock, iword, ibyte)
            inum = ibuffer
	  if (ibuffer .eq. 0) then
            print '(a,a,/,a,a)',
     &         '=>ERROR in getunit options list:',
     &         ' ''buffer='' requires a size',
     &         ' options=',options
            errors = errors + 1
          end if
#ifdef cray_ymp
c
c       scan for additional cray ymp access methods
c
        else if(scan_option(options, i, 'word', saccess, 'word')) then
          if(scan_option(options, i, 'buffer', sdummy, sdummy)) then
          end if
          if (scan_size(options, i, ibuffer, iblock, iword, ibyte)) then
            inum = ibuffer
          end if
          call set2block (ibuffer, iblock, iword, ibyte)
	  if (ibuffer .eq. 0) then
            print '(a,a,/,a,a)',
     &         '=>ERROR in getunit options list:',
     &         ' ''word='' requires a buffer size',
     &         ' options=',options
            errors = errors + 1
          end if
        else if(scan_option(options, i, 'sdsalloc', 
     &                         saccess, 'sdsalloc')) then
	  if(scan_option(options, i, 'maxsize', sdummy, sdummy)) then
	  end if
	  if(scan_option(options, i, 'max', sdummy, sdummy)) then
	  end if
          if (scan_size(options, i, maxsize, iblock, iword, ibyte)) then
            inum = maxsize
          end if
          call set2block (maxsize, iblock, iword, ibyte)
	  if (maxsize .eq. 0) then
            print '(a,a,/,a,a)',
     &         '=>ERROR in getunit options list:',
     &         ' ''sdsalloc'' requires a maximum size in blocks',
     &         ' options=',options
            errors = errors + 1
          end if
c
c       scan for additional cray ymp layers
c
        else if(scan_option(options, i, 'mr', sffio, 'mr')) then
	  call new_layer(layer, sblocking, psblocking, slayer, sffio, 
     &                   inum, ibuffer, maxsize, iblock, iword, ibyte,
     &                   ssave, sovfl)
	  if(scan_option(options, i, 'maxsize', sdummy, sdummy)) then
	  end if
	  if(scan_option(options, i, 'max', sdummy, sdummy)) then
	  end if
          if (scan_size(options, i, maxsize, iblock, iword, ibyte)) then
            inum = maxsize
          end if
          call set2block (maxsize, iblock, iword, ibyte)
        else if(scan_option(options, i, 'sds', sffio, 'sds')) then
	  call new_layer(layer, sblocking, psblocking, slayer, sffio, 
     &                   inum, ibuffer, maxsize, iblock, iword, ibyte,
     &                   ssave, sovfl)
	  if(scan_option(options, i, 'maxsize', sdummy, sdummy)) then
	  end if
	  if(scan_option(options, i, 'max', sdummy, sdummy)) then
	  end if
          if (scan_size(options, i, maxsize, iblock, iword, ibyte)) then
            inum = maxsize
          end if
          call set2block (maxsize, iblock, iword, ibyte)
        else if(scan_option(options, i, 'scratch', ssave, 'scr')) then
        else if(scan_option(options, i, 'scr', ssave, 'scr')) then
        else if(scan_option(options, i, 'save', ssave, 'save')) then
        else if(scan_option(options, i, 'ovfl', sovfl, 'ovfl')) then
        else if(scan_option(options, i, 'novfl', sovfl, 'novfl')) then
        else if(scan_option(options, i, 'cache', sffio, 'cache')) then
	  call new_layer(layer, sblocking, psblocking, slayer, sffio, 
     &                   inum, ibuffer, maxsize, iblock, iword, ibyte,
     &                   ssave, sovfl)
          if(scan_option(options, i, 'buffer', sdummy, sdummy)) then
          end if
          if (scan_size(options, i, ibuffer, iblock, iword, ibyte)) then
          end if
          call set2block (ibuffer, iblock, iword, ibyte)
          if(scan_option(options, i, 'pages', sdummy, sdummy)) then
          end if
          if(scan_number(options, i, ipage, snum)) then
            inum = ipage
          end if
	  if (ibuffer .eq. 0 .or. ipage .eq. 0) then
            print '(a,a,/,a,a)',
     &         '=>ERROR in getunit options list:'
     &,        ' ''cache'' requires a buffer size and a number of pages',
     &,        ' options=',options
            errors = errors + 1
          end if
#endif

c
c       scan options by abbreviated names
c
        else if(scan_option(options, i, 's', saccess, ssequential)) then
        else if(scan_option(options, i, 'd', saccess, sdirect)) then
          if(scan_option(options, i, 'recl', sdummy, sdummy)) then
          end if
          if(scan_size(options, i, irec_len, iblock, iword, ibyte)) then
          end if
	  call setreclen (irec_len,srec_len,sform,iblock,iword,ibyte)
	  if (irec_len .eq. 0) then
            print '(a,a,/,a,a)',
     &         '=>ERROR in getunit options list:',
     &         ' ''direct'' requires a record length',
     &         ' options=',options
            errors = errors + 1
          end if
        else if(scan_option(options, i, 'f', sform, sformatted)) then
        else if(scan_option(options, i, 'u', sform, sunformatted)) then
        else if(scan_option(options, i, 'rw', saction, sreadwrite)) then
        else if(scan_option(options, i, 'ro', saction, sread)) then
        else if(scan_option(options, i, 'wo', saction, swrite)) then
        else if(scan_option(options, i, 'r', sposition, srewind)) then
        else if(scan_option(options, i, 'a', sposition, sappend)) then
#ifdef cray_ymp
        else if(scan_option(options, i, 'w', saccess, 'word')) then
          if(scan_option(options, i, 'buffer', sdummy, sdummy)) then
          end if
          if (scan_size(options, i, ibuffer, iblock, iword, ibyte)) then
            inum = ibuffer
          end if
          call set2block (ibuffer, iblock, iword, ibyte)
	  if (ibuffer .eq. 0) then
            print '(a,a,/,a,a)',
     &         '=>ERROR in getunit options list:',
     &         ' ''word='' requires a buffer size',
     &         ' options=',options
            errors = errors + 1
          end if
#endif
c
c       scan for size options with no introductory keyword
c       the introductory keyword must be several options back
c
        else if(scan_size(options, i, inum, iblock, iword, ibyte)) then
        else
          print '(4a)', '=>ERROR:  Getunit unable to parse options=',
     &                    options(1:i-1), '|',
     &                    options(i:len_trim(options))
          errors = errors + 1
          i = i + 1
        end if
#ifdef cray_ymp
	if (sblocking .ne. psblocking) then
	  call new_layer(layer, sblocking, psblocking, slayer, sffio, 
     &                   inum, ibuffer, maxsize, iblock, iword, ibyte,
     &                   ssave, sovfl)
	end if
#endif
      goto 100
c
c     end of scan loop
c     a do loop is not be used because i is changed within the loop
c
101   continue
#ifdef cray_ymp
      if (layer .gt. 0 .or. sieee .eq. 'ieee') then
        sffio = 'END OPTIONS'
        call new_layer(layer, sblocking, psblocking, slayer, sffio,
     &                 inum, ibuffer, maxsize, iblock, iword, ibyte,
     &                 ssave, sovfl)
      end if
#endif
      if (errors .gt. 0) then
        print '(/a/a,i3,a)'
     &,   'STOP in getunit'
     &,   'There were',errors,' errors scanning the list of I/O options'
        stop
      end if
c
c     in case neither "formatted" or "unformatted" is specified, choose
c     the proper default for sequential and direct access files.
c
      if (sform .eq. ' ') then
        if (saccess .eq. ssequential) then
	  sform = sformatted
	else if (saccess .eq. sdirect) then
	  sform = sunformatted
	end if
      end if
c-----------------------------------------------------------------------
c     See if the file has previously been opened in this run
c-----------------------------------------------------------------------
c
#ifdef debug_iomngr
      print '(/,a,a,/)', 'Searching for a previously opened file named '
     &,                filename
      print '(a,tr2,a,/)', 'unit', 'file name'
#endif
c
      do i=1,nfiles 
#ifdef debug_iomngr
        print '(i3,tr2,a30,a)',iunit(i),fname(i)
#endif
        if (fname(i) .eq. filename) then
          sstatus = sold
          if (iunit(i) .ne. 0) then
c
c           This happens when a call relunit was executed for this file,
c           but due to system bugs, it was necessary to merely hide the
c           file while keeping it connected to a unit to permit later 
c           appending.
c
            iounit = iunit(i)
            ifile(iounit) = i
            inuse(iounit) = .true.
            if (sposition .eq. sappend) then
              skip_open = .true.
            else
c
c             A previously hidden file that is now being opened with
c             POSITION=REWIND.  It should first be closed to permit 
c             reading from it.
c
#ifdef verbose_iomngr
              write (stdout,'(a,i3,a,a)')
     &            'closing   i/o unit', iounit, '  ',
     &            fname (i)(1:LEN_TRIM(fname(i)))
#endif
              call close_file (iounit, filename)
              skip_open = .false.
c
#ifdef cray_ymp
c             This might also be a good time to do the concatenation
c             trick for unformatted ieee files which have pre-existing
c             parts, so that the reopen with POSITION=REWIND will reach
c             the beginning of the .old part.
c
              call cray_ymp_ieee_concat_kludge(fname(i))
c
#endif
            end if
          else
c
c           This file is now closed, but had previously been opened.
c           Get a new unit number and reconnect it.
c
            call getunitnumber (iounit)
            ifile(iounit) = i
            iunit(i) = iounit 
c
c           STATUS=UNKNOWN may cause rewind on some computers.
c           This is incompatible with append
c      
            if (sposition .eq. sappend) then
              sstatus = sold
            end if
            skip_open = .false.
#ifdef cray_ymp
c
c           Bug fix:
c           if a previously accessed ieee append file is now closed
c           instead of hidden, it must be copied to filename.old
c           and a new copy started.  Later, when release_all is called,
c           filename.old and filename will be concatenated to produce
c           the correct output file.
c
            if (saccess .eq. ssequential .and. sform .eq. sunformatted
     &          .and. sieee .eq. 'ieee') then
# ifdef debug_iomngr
              print '(tr1,a,i3,tr2,a,a)',
     &     'a. hiding unit ', iounit, 'filename = ',fname(ifile(iounit))
# endif
              hide_file(ifile(iounit)) = .true.
              sposition = srewind
              string =  'mv '// filename(1:len_trim(filename)) // ' '
     &                  // filename(1:len_trim(filename)) // '.old' 
# ifdef debug_iomngr
c  tp:
	      call assign('assign -V')
c  tp:
# endif
# ifdef verbose_iomngr
              print '(a,a)',
     &         'RENAMING PRE-EXISTING APPEND FILE--SHELL CMD: ', string
# endif
              istatus = ishell (string)
              sstatus = snew
            end if
#endif
          end if
          go to 201
        end if
      end do
c
c         
c     This file has not previously been opened in this run.
c
c     Get a free unit number and bind it to the file
c     unless the request is 'sdsalloc' which does not
c     use a true i/o unit number.
c
      if (saccess .ne. 'sdsalloc') then
        call getunitnumber (iounit)
        call link_unit (iounit, filename)
        if (sposition .eq. sappend) then
c
c         If file pre-exists, append to it, but if no such file
c         exists, open a new file and rewind it.
c
          inquire (FILE=fname(nfiles), EXIST=exists, NAMED=named)
#ifdef debug_iomngr2
          write (stdout,*) 'file ',filename(1:len_trim(filename)),
     &          ' exists = ', exists
          write (stdout,*) 'file ',filename(1:len_trim(filename)),
     &          ' named = ', named
#endif
          if (.not. exists .or. .not. named) then
#ifdef debug_iomngr
            print '(a,a,a)',
     &            'file ',fname(nfiles)(1:len_trim(fname(nfiles))),
     &            ' does not exist--status changed to NEW' 
#endif
            sstatus = snew
            sposition = srewind
          else
#ifdef debug_iomngr
            print '(a,a,a)',
     &            'file ',fname(nfiles)(1:len_trim(fname(nfiles))),
     &            ' already exists--status changed to OLD' 
#endif
	    sstatus = sold
          end if

#ifdef cray_ymp
c
c         Bug fix:  The cray operating system does not support
c         appending to sequential access, unformatted, ieee files.
c         This option is handled as follows:  (1) it is marked as
c         a file that must be hidden [and not closed] when the
c         program asks to release it.  (2) if it is to be reopened
c         with append [by getunit] it is simply un-hidden. (3) if
c         it is to be reopened with rewind, it is closed and then
c         reopened [to permit switching from write to read]. (4) If
c         such a file exists before the first time it is opened by
c         getunit in this run, the old file is renamed with a suffix
c         '.old' and a new file is opened.  At the conclusion of a
c         run by a call to 'release_all', or when reopened with
c         POSITION=REWIND, the new file is closed and the old and new
c         files are concatenated together and renamed with the correct
c         user-specified name.  In all cases, the result at end-of-run
c         appears as though append were supported.
c
          if (saccess .eq. ssequential .and. sform .eq. sunformatted
     &        .and. sieee .eq. 'ieee') then
# ifdef debug_iomngr
            print '(tr1,a,i3,tr2,a,a)',
     &     'b. hiding unit ', iounit, 'filename = ',fname(ifile(iounit))
# endif
            hide_file(ifile(iounit)) = .true.
            sposition = srewind
            if (exists) then
              string =  'mv '// filename(1:len_trim(filename)) // ' '
     &                  // filename(1:len_trim(filename)) // '.old' 
# ifdef debug_iomngr
c  tp:
	      call assign('assign -V')
c  tp:
# endif
# ifdef verbose_iomngr
              print '(a,a)',
     &         'RENAMING PRE-EXISTING APPEND FILE--SHELL CMD: ', string
# endif
              istatus = ishell (string)
              sstatus = snew
            end if
          end if
#endif
        end if
      end if
      skip_open = .false.
c
201   continue 
#ifdef cray_ymp
c
c      bug fix
         if (saccess .eq. ssequential .and. sform .eq. sunformatted
     &       .and. sieee .eq. 'ieee') then
#ifdef debug_iomngr
           print '(tr1,a,i3,tr2,a,a)',
     &   'c. hiding unit ', iounit, 'filename = ',fname(ifile(iounit))
#endif
          hide_file(ifile(iounit)) = .true.
        end if
#endif
c 
#ifdef cray_ymp
c
c     generate Cray "ffio" assign environment
c
c     ieee is a non-native format on cray ymp
c     cray blocking and buffer size options are set here
c
      if (sieee .eq. 'ieee') then
        snumeric = '-N ieee'
      end if
      if (layer .gt. 0 .or. sieee .eq. 'ieee') then
        string = 'assign -I -F'
        do l=1,layer-1
          if (l .eq. 1) then
            string = string(1:len_trim(string)) // ' ' //
     &             slayer(l)(1:len_trim(slayer(l)))
          else
            string = string(1:len_trim(string)) // ',' //
     &             slayer(l)(1:len_trim(slayer(l)))
          end if
#ifdef debug_iomngr
          print *,l,' ',string
#endif
        end do
        if (sieee .eq. 'ieee') then
	  if (layer .eq. 1) then
            string = string(1:len_trim(string)) // ' ' // 'f77'
	  end if
          string = string(1:len_trim(string)) // ' ' // snumeric
        end if
	if (saccess .eq. 'word') then
          call i2a(iounit,aiounit)
          string = string(1:len_trim(string)) // ' u:' // aiounit
	else
c
c         tp:  try assigning a unit as well as a filename.
c         tp:  this seems to fail.  assign may not permit it.
c
          call i2a(iounit,aiounit)
          string = string(1:len_trim(string)) // ' f:' //
     &           fname(ifile(iounit))
ccccc          string = string(1:len_trim(string)) // ' u:' // aiounit
        end if
#if defined debug_iomngr || defined show_assigns
        print *, '       '
        print *, string
#endif
        call assign (string, ierr)
        if (ierr .ne. 0) then
          print *,'bug in iomngr: ierr=', ierr
          call assign ('assign -V')
        end if
      end if
#endif
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
      
      write (aiounit,'(i2)') iounit
#ifdef debug_iomngr
c
#if defined verbose_iomngr || defined io_history
      string =         'acquiring i/o unit ' // aiounit //
     &            '  ' // filename(1:LEN_TRIM(filename)) //
     &            '  ' // saccess(1:LEN_TRIM(saccess)) //
     &            '  ' // sieee(1:LEN_TRIM(sieee)) //
     &            '  ' // sform(1:LEN_TRIM(sform)) //
     &            '  ' // sstatus(1:LEN_TRIM(sstatus)) //
     &            '  ' // sposition(1:LEN_TRIM(sposition))
      if (saccess .eq. sdirect) then
        string = string(1:len_trim(string)) // ' RECL=' // srec_len
      end if
#endif
c
#ifdef io_history
c     write (iohist,'(a)') string(1:len_trim(string))
#endif
#ifdef verbose_iomngr
      write (stdout,'(a)') string(1:len_trim(string))
#endif
#endif
c
c-----------------------------------------------------------------------
c
#ifdef sgi
c     fix for SGI Fortran which does not support POSITION=
c
      if (sposition .eq. sappend .and. saccess .eq. ssequential) then
        saccess = sappend
        sposition = ' '
      end if
#endif
c
c     set implied file attributes corresponding to 'stream' options.
c
#ifdef sgi
c
c     make 'stream' mean: sequential access, system
c
      if (sblocking .eq. 'null') then
        saccess = ssequential
        sform = 'SYSTEM'
      end if
#endif
#ifdef cray_ymp
c
c     make 'stream' mean: unblocked, sequential access, unformatted
c
      if (sblocking .eq. 'null') then
	saccess = ssequential
	sform = sunformatted
      end if
#endif
c-----------------------------------------------------------------------
c
c     select open statement compatible with options and computer platform
c
c-----------------------------------------------------------------------
      if (skip_open) then
        goto 301
      end if
#ifdef debug_iomngr
        print *,' ...executing open: ',
     &          sform, saccess, sstatus, sposition, 'recl=', irec_len
# ifdef cray_ymp
        call assign ('assign -V')
# endif
#endif
      if (saccess .eq. ssequential .or. saccess .eq. sappend) then
	open (iounit, 
     &        FILE=filename, 
     &        FORM=sform, 
     &        ACCESS=saccess, 
     &        STATUS=sstatus
#ifndef sgi
     &,       POSITION=sposition
#endif
     &        )
        unix_name(iounit) = .true.
      else if (saccess .eq. sdirect) then
        open (iounit, 
     &        FILE=filename, 
     &        FORM=sform, 
     &        ACCESS=saccess,
     &        RECL=irec_len, 
     &        STATUS=sstatus
     &        )
        unix_name(iounit) = .true.
#ifdef cray_ymp
      else if (saccess .eq. 'word') then
	call set2block (ibuffer, iblock, iword, ibyte)
	if (ibuffer .eq. 0) then
          print '(a,a,/,a,a)',
     &       '=>ERROR in getunit options list:',
     &       ' ''word'' requires a buffer size in words or blocks',
     &       ' options=',options
        end if
        call wopen (iounit, ibuffer, 1)
        unix_name(iounit) = .false.
      else if (saccess .eq. 'sdsalloc') then
c       allocate blocks on sds and return block number of first block
c       in place of unit number in argument iounit
        call set2block (maxsize, iblock, iword, ibyte)
	if (maxsize .eq. 0) then
          print '(a,a,/,a,a)',
     &     '=>ERROR in getunit options list:',
     &     ' ''sdsalloc'' requires a max size ',
     &     ' options=',options
        end if
        iounit = sdsalloc (maxsize, ierr)
        if (ierr .ne. 0) then
          print '(a,a,/,a,i8,a)',
     &      '=>ERROR:  sdsalloc failed for file ', filename,
     &      '          cannot allocate ',maxsize, ' blocks.'
          stop
        unix_name(iounit) = .false.
        end if
#endif
      end if
c
301   continue
#if defined verbose_iomngr || defined io_history
      string =         'acquiring i/o unit ' // aiounit //
     &            '  ' // filename(1:LEN_TRIM(filename)) //
     &            '  ' // saccess(1:LEN_TRIM(saccess)) //
     &            '  ' // sieee(1:LEN_TRIM(sieee)) //
     &            '  ' // sform(1:LEN_TRIM(sform)) //
     &            '  ' // sstatus(1:LEN_TRIM(sstatus)) //
     &            '  ' // sposition(1:LEN_TRIM(sposition))
      if (saccess .eq. sdirect) then
        string = string(1:len_trim(string)) // ' RECL=' // srec_len
      end if
#endif
c
#ifdef io_history
      write (iohist,'(a)') string(1:len_trim(string))
#endif
#ifdef verbose_iomngr
      write (stdout,'(a)') string(1:len_trim(string))
#endif
c
      return
      end



      subroutine relunit (iounit)
c
c=======================================================================
c     Relunit releases a Fortran i/o unit number currently in use and
c     closes the file.
c     see Getunit which gets a Fortran i/o unit number.
c
c     author c. h. goldberg          e-mail: chg@gfdl.gov
c=======================================================================
c
#include "iomngr.h"
      character * 2 aiounit
c
#ifdef cray_ymp
c       We only pretend to release and close ieee, append files as a
c       work-around for a cray bug in which ieee unformatted sequential
c       files cannot be properly opened with POSITION=APPEND.  These
c       files are flagged by setting hide_file=.true..
c
c       Points to check carefully include "sdsalloc" unit numbers, which
c       are really starting block numbers in sds, and "word" files,
c       which have a unit number, but not a file name in the actual 
c       wopen statement.
c
#endif
c
      ifn = ifile(iounit)
#ifdef verbose_iomngr
      write (stdout,'(a,i3,a,a)')
     &            'releasing i/o unit', iounit, '  ',
     &            fname(ifn)(1:LEN_TRIM(fname(ifn)))
#endif
#ifdef io_history
      write (iohist,'(a,i3,a,a)')
     &            'releasing i/o unit', iounit, '  ',
     &            fname(ifn)(1:LEN_TRIM(fname(ifn)))
#endif
c
      if (iounit .ge. 1 .and. iounit .le. maxunit) then
        inuse(iounit) = .false.
        if (.not. hide_file(ifn) .and. unix_name(iounit)) then
#ifdef verbose_iomngr
          write (stdout,'(a,i3,a,a)')
     &       'closing   i/o unit', iounit, '  ',
     &       fname(ifn)(1:LEN_TRIM(fname(ifn)))
#endif
          call close_file (iounit, fname(ifn))
          iunit(ifn) = 0
        end if
        ifile(iounit) = 0
      else
        print *,'=>ERROR: iounit is bad. iounit=',iounit
        stop 'relunit'
      endif

      return
      end



      subroutine release_all

c=======================================================================
c
c     Releases all user i/o units currently in use and their file names.
c     All files except stdin, stdout, and stderr are closed.
c     Hidden files are closed at this time.
c 
c     Ater a call to release_all, it is safe to execute a sectiion of
c     code that does not use the iomngr to allocate units.
#ifdef cray_ymp
c
c     Appends to pre-existing ieee unformatted sequential files take
c     place at this time using the shell command 'cat'.
c     Unicos "word" and "sdsalloc" files can never be closed because
c     they have no unix names and the iomngr cannont guarantee to use
c     the same unit number to reopen them.
#endif
c
c=======================================================================

#include "iomngr.h"
      logical exists
      character *160 string, oldfname
c
#ifdef verbose_iomngr
      write(stdout, '(/,a,/)') 'Releasing and closing all files'
      write(stdout, '(a5,tr1,a,tr1,a,tr1,a)')
     &              'unit', 'hidden', 'file name'
#endif
c
      do ifn=1,nfiles
        if (iunit(ifn) .ne. 0) then
          if ((inuse(iunit(ifn)) .or. hide_file(ifn)) .and.
     &        unix_name(iunit(ifn)) .and.
     &        fname(ifn) .ne. 'reserved' .and.
     &        fname(ifn) .ne. 'stdin' .and.
     &        fname(ifn) .ne. 'stdout' .and.
     &        fname(ifn) .ne. 'stderr' .and.
     &        fname(ifn) .ne. 'stdout/stderr') then
#ifdef verbose_iomngr
            write(stdout, '(tr2,i2,tr3,l2,tr4,a)') 
     &            iunit(ifn), hide_file(ifn),
     &              fname(ifn)(1:len_trim(fname(ifn)))
#endif
            call close_file (iunit(ifn), fname(ifn))
#ifdef cray_ymp
c
c           Check to see if any hidden files are the new part of a
c           pre-existing ieee sequential unformatted append file.
c
            if (hide_file(ifn)) then
              call cray_ymp_ieee_concat_kludge(fname(ifn))
            end if
#endif
c
c           mark file as released
c
            inuse(iunit(ifn)) = .false.
            iunit(ifn) = 0
          end if
        end if
      end do


      return
      end


      subroutine close_file (iounit, filename)
c=======================================================================
c
c       close Fortran unit iounit
c
#ifdef cray_ymp
c       release all assigned file attributes for this unit number
c
#endif
c=======================================================================
c
c
#ifdef cray_ymp
# include "iomngr.h"
#endif
        character * 64 filename, string1, string2
        character * 2 aiounit
c
        close (iounit)
#ifdef cray_ymp
        string1 = 'assign -R'
        call i2a(iounit,aiounit)
        string1 = string1(1:len_trim(string1)) // ' u:' // aiounit
        call assign (string1, ierr)
        string2 = 'assign -R f:' // filename(1:len_trim(filename))
        call assign (string2, ierr)
# if defined debug_iomngr || defined show_assigns
        print *, string1
        print *, string2
# endif
#endif
        return
        end





#ifdef cray_ymp
c=======================================================================
c
      subroutine cray_ymp_ieee_concat_kludge(filename)
c
c=======================================================================
c
c     the only way to append to an unformatted sequential ieee file
c     on the current cray ymp operating system is to use the unix
c     command 'cat' to concatenate the preexixting portion of the file
c     with the new portion.
c
c     this subroutine inquires about 'filename.old', and if it is 
c     found, concatenates 'filename.old' and 'filename', leaving the
c     result in 'filename'.
c
c=======================================================================
c
      character * (*) filename
      character * 128 oldfname, string
      logical exists
c
      oldfname = filename(1:len_trim(filename)) // '.old'
      inquire (FILE=oldfname, EXIST=exists)
      if (exists) then
        string = 'cat ' // filename(1:len_trim(filename)) //
     &           ' >> ' // oldfname(1:len_trim(oldfname)) 
# ifdef verbose_iomngr
        print '(a,a)', 'SHELL CMD: ', string
# endif
        istatus = ishell (string)
        string = 'mv ' // oldfname(1:len_trim(oldfname)) //
     &           ' ' // filename(1:len_trim(filename))
# ifdef verbose_iomngr
        print '(a,a)', 'SHELL CMD: ', string
# endif
        istatus = ishell (string)
      end if
c
      return
      end
#endif

      subroutine showunits

c=======================================================================
c
c     lists all i/o units currently in use and their file names.
c
c=======================================================================

#include "iomngr.h"

      write(stdout, '(/,a,/)') 'i/o units in use'
      write(stdout, '(a5,tr1,a,tr1,a)') ' unit', 'hidden', 'file name'
      do i=1,maxunit
        if (inuse(i)) then
          write(stdout, '(i4,tr3,l2,tr4,a)') 
     &              i, hide_file(ifile(i)),
     &              fname(ifile(i))(1:len_trim(fname(ifile(i))))
        end if
      end do

#ifdef io_history
      write(iohist, '(/,a,/)') 'i/o units in use'
      write(iohist, '(a4,tr1,a)') 'unit', 'file name'
      do i=1,maxunit
        if (inuse(i)) then
          write(iohist, '(i4,tr3,l2,tr4,a)') 
     &              i, hide_file(ifile(i)),
     &              fname(ifile(i))(1:len_trim(fname(ifile(i))))
        end if
      end do
#endif
      return
      end
c=======================================================================


      subroutine showfiles

c=======================================================================
c
c     lists all file names that have been used [in chronological order]
c
c=======================================================================

#include "iomngr.h"

      write(stdout, '(/,a,/)') '  Files Used in This Run'
      write(stdout, '(a5,tr1,a,tr1,a)') '     ', 'hidden', 'file name'
      do i=1,nfiles
        write(stdout, '(i4,tr3,l2,tr4,a)') 
     &            i, hide_file(i),
     &              fname(i)(1:len_trim(fname(i)))
      end do

#ifdef io_history
      write(iohist, '(/,a,/)') '  Files Used in This Run'
      write(iohist, '(a5,tr1,a,tr1,a)') '     ', 'hidden', 'file name'
      do i=1,nfiles
          write(iohist, '(i4,tr3,l2,tr2,a)') 
     &              i, hide_file(i),
     &              fname(i)(1:len_trim(fname(i)))
      end do
c     call relunit (iohist)
#endif
      return
      end
c=======================================================================

      subroutine link_unit (iounit, filename)

c=======================================================================
c     links an i/o unit number to a new file name
c-----------------------------------------------------------------------

#include "iomngr.h"
      character * (*) filename

      nfiles = nfiles + 1
      inuse(iounit) = .true.
      fname(nfiles) = filename
      ifile(iounit) = nfiles
      iunit(nfiles) = iounit
c
      return
      end

c=======================================================================

      subroutine setreclen (irec_len,srec_len,sform,iblock,iword,ibyte)

c-----------------------------------------------------------------------
c     calculates record length for direct access files
c-----------------------------------------------------------------------
      character *(*) srec_len, sform
      character *15 sformatted
      data sformatted   /'FORMATTED'/
c
#ifdef cray_ymp
      call set2byte (irec_len, iblock, iword, ibyte)
#else
# ifdef sgi
      if (sform .eq. sformatted) then
        call set2byte (irec_len, iblock, iword, ibyte)
      else
        call set2word (irec_len, iblock, iword, ibyte)
      end if
# else
      if (sform .eq. sformatted) then
        call set2byte (irec_len, iblock, iword, ibyte)
      else
        call set2word (irec_len, iblock, iword, ibyte)
      end if
# endif
#endif
      call i2a(irec_len, srec_len)
c
      return
      end

c=======================================================================

      subroutine set2byte (inum, iblock, iword, ibyte)

c-----------------------------------------------------------------------
c     if first argument is zero, sets first argument to value in bytes
c-----------------------------------------------------------------------
      if (inum .eq. 0) then
	if (ibyte .ne. 0) then
	  inum = ibyte
	else if (iword .ne. 0) then
#ifdef cray_ymp
	  inum = 8 * iword
#else
# ifdef sgi
          inum = 4 * iword
# else
	  inum = 4 * iword
# endif
#endif
#ifdef cray_ymp
	else if (iblock .ne. 0) then
	  inum = 512 * 8 * iblock
#endif
	end if
      end if

      return
      end

c=======================================================================

      subroutine set2word (inum, iblock, iword, ibyte)

c-----------------------------------------------------------------------
c     if first argument is zero, sets first argument to value in bytes
c-----------------------------------------------------------------------
      if (inum .eq. 0) then
	if (ibyte .ne. 0) then
#ifdef cray_ymp
	  inum = 1 + (ibyte - 1) / 8
#else
# ifdef sgi
          inum = 1 + (ibyte - 1) / 4
# else
          inum = 1 + (ibyte - 1) / 4
# endif
#endif
	else if (iword .ne. 0) then
	  inum = iword
#ifdef cray_ymp
	else if (iblock .ne. 0) then
	  inum = 512 * iblock
#endif
	end if
      end if

      return
      end

#ifdef cray_ymp
c=======================================================================
      subroutine set2block (inum, iblock, iword, ibyte)

c-----------------------------------------------------------------------
c     if first argument is zero, sets first argument to value in bytes
c-----------------------------------------------------------------------
      if (inum .eq. 0) then
	if (ibyte .ne. 0) then
	  inum = 1 + (ibyte - 1) / 4096
	else if (iword .ne. 0) then
	  inum = 1 + (iword - 1) / 512
	else if (iblock .ne. 0) then
	  inum = iblock
	end if
      end if

      return
      end

#endif
c=======================================================================
      subroutine tolower (string)

c-----------------------------------------------------------------------
c     converts all alphabetic characters in string to lowercase
c-----------------------------------------------------------------------

      character*(*) string
      character*1 c
      parameter (lcshift=32)

      do i=1,len(string)
        c = string(i:i)
        if ('A' .le. c .and. c .le. 'Z') then
          string(i:i) = char(ichar(c) + lcshift)
        end if
      end do

      return
      end


c=======================================================================

      subroutine i2a (i, a)

      character * (*) a
      character * 15 string

      write (string, '(i15)') i
      do j=1,15
        if (string(j:j) .ne. ' ') go to 101
      end do
101   continue
      a = string(j:)

      return
      end

c=======================================================================

      function scan_option (options, i, spattern, svar, soption)

c-----------------------------------------------------------------------
c     scans the options string for spattern starting at position i.
c     if found, svar is set to soption.
c-----------------------------------------------------------------------

      logical scan_option
      character*(*) options, spattern, svar, soption
      character*1 c

c     skip delimiters

100   continue
      if (i .le. len(options)) then
	c = options(i:i)
        if (c.eq.'=' .or. c.eq.' ' .or. c.eq.',' .or. c.eq.':') then
	  i = i + 1 
          goto 100
          end if
        end if
c 
c     scan for option 
c
      scan_option = .false.
      j = i + len_trim(spattern) - 1
      if (j .le. len(options)) then
	if (options(i:j) .eq. spattern) then
	  scan_option = .true.
          svar = soption
#ifdef debug_iomngr
	  print '(tr1,a,a,tr1,a,i3)',
     &     'scan option: found option-->',spattern
#endif
	  i = i + len_trim(spattern)
        end if
      end if

      return
      end

c=======================================================================
c
      function scan_size(options, i, inum, iblock, iword, ibyte)
c
c-----------------------------------------------------------------------
c     scans the options string for a size specification starting at
c     position i.
c     if found, one of inum, iblock, iword, ibyte is nonzero.
c-----------------------------------------------------------------------
      logical scan_size, scan_number, scan_option
      character*(*) options
      character*15 snum, sblock, sword, sbyte, sdummy

c     scan for size options
c
      inum = 0
      iblock = 0
      iword = 0
      ibyte = 0
      scan_size = .false.
c
      if(scan_number(options, i, inum, snum)) then
        scan_size = .true.
#ifdef cray_ymp
      else if(scan_option(options, i, 'blocks', sdummy, sdummy)) then
        if (scan_number(options, i, iblock, sblock)) then
          scan_size = .true.
        else
          print '(a,a,/,a,a)',
     &       '=>ERROR in getunit options list:',
     &       ' ''blocks='' requires a size in 512-word blocks',
     &       ' options=',options
        end if
#endif
      else if(scan_option(options, i, 'words', sdummy, sdummy)) then
        if (scan_number(options, i, iword, sword)) then
          scan_size = .true.
        else
          print '(a,a,/,a,a)',
     &       '=>ERROR in getunit options list:',
     &       ' ''words='' requires a size in words',
     &       ' options=',options
        end if
      else if(scan_option(options, i, 'bytes', sdummy, sdummy)) then
        if (scan_number(options, i, ibyte, sbyte)) then
          scan_size = .true.
        else
          print '(a,a,/,a,a)',
     &       '=>ERROR in getunit options list:',
     &       ' ''bytes='' requires a size in bytes',
     &       ' options=',options
        end if
      end if
c
      return
      end

c=======================================================================
c
      function scan_number (options, i, n, sn)

c-----------------------------------------------------------------------
c     scans the options string for a number starting at position i.
c     if found, n is the number as an integer and sn as a string.
c-----------------------------------------------------------------------

      logical scan_number
      character*(*) options, sn
      character*1 c


c     skip delimiters

100   continue
      if (i .le. len(options)) then
        c = options(i:i)
        if (c.eq.'=' .or. c.eq.' ' .or. c.eq.',' .or. c.eq.':') then
          i = i + 1
          goto 100
        end if
      end if
c
c     scan for digits
c
      if (i .le. len(options)) then
        c = options(i:i)
        if ('0' .le. c .and. c .le. '9') then
          sn = c
          n = ichar(c) - ichar('0')
          scan_number = .true.
200       continue
          i = i + 1
          if (i .le. len(options)) then
            c = options(i:i)
            if ('0' .le. c .and. c .le. '9') then
              sn = sn(1:len_trim(sn)) // c
              n = 10 * n + ichar(c) - ichar('0')
              goto 200
            end if
          end if
#ifdef debug_iomngr
	  print '(tr1,a,a,tr1,a,i3)', 
     &     'scan number: found option-->',sn(1:len_trim(sn))
#endif
        else
          scan_number = .false.
        end if
      end if

      return
      end

#ifdef cray_ymp
c=======================================================================
c

      subroutine new_layer(layer, sblocking, psblocking, slayer, sffio,
     &                     inum, ibuffer, maxsize, iblock, iword, ibyte,
     &                     ssave, sovfl)
c
c-----------------------------------------------------------------------
c     completes character string for options on previous ffio layer
c-----------------------------------------------------------------------
c
      character *(*) sblocking, psblocking, slayer, sffio, ssave, sovfl
      dimension slayer(6)
      character *15 snum, sbuffer
      character *40 string

      layer = layer + 1
#ifdef debug_iomngr
      if (layer .gt. 1) then
        print *,'    lower layer-->',slayer(layer-1)
      end if
# ifdef debug_iomngr2
      print *,'    sffio-->',sffio
      print *,'    sblocking-->',sblocking
      print *,'    psblocking-->',psblocking
      print *,'    inum=',inum,' ibuffer=',ibuffer,' maxsize=',maxsize
      print *,'    iblock=',iblock,' iword=',iword,' ibyte=',ibyte
# endif
#endif
      if (sffio .ne. ' ') then
	slayer(layer) = sffio
	sffio = ' '
      else if (sblocking .ne. psblocking) then
        if (psblocking .eq. ' ') then
	  slayer(layer) = sblocking
          psblocking = sblocking
        else 
          print '(a,a)',
     &          'ERROR:  only one blocking option is permitted',
     &          ' for a file'
          print '(a,a,a,a)', 'options = ...',
     &          psblocking(1:len_trim(psblocking)), '...',
     &          sblocking(1:len_trim(sblocking))
          stop '==>iomngr'
        end if
      else
	print '(a,a)', 
     &        'ERROR:  only blocking or ffio layer options should',
     &        ' invoke subroutine new_layer'
      end if
      if (layer .ne. 1) then
	if (slayer(layer-1) .eq. 'cos') then
	  call set2block (inum, iblock, iword, ibyte)
	  if (ibuffer .ne. 0) then
            call i2a (ibuffer, snum)
	    string = slayer(layer-1)(1:len_trim(slayer(layer-1))) //
     &               ':' // snum(1:len_trim(snum))
	    slayer(layer-1) = string
          end if
	else if (slayer(layer-1) .eq. 'f77' .or.
     &           slayer(layer-1) .eq. 'null') then
          if (ibuffer .ne. 0) inum = ibuffer
	  call set2byte (inum, iblock, iword, ibyte)
	  if (inum .ne. 0) then
            call i2a (inum, snum)
	    string = slayer(layer-1)(1:len_trim(slayer(layer-1))) //
     &               '.::' // snum(1:len_trim(snum))
          else
            string = slayer(layer-1)
          end if
          slayer(layer-1) = string
	else if (slayer(layer-1) .eq. 'sds' .or.
     &           slayer(layer-1) .eq. 'mr') then
	  call set2block (maxsize, iblock, iword, ibyte)
	  if (maxsize .ne. 0) then
	    call i2a (maxsize, snum)
	  else
	    snum = ' '
	  end if
	  string = slayer(layer-1)(1:len_trim(slayer(layer-1))) //
     &             '.' // ssave(1:len_trim(ssave)) // '.' // 
     &             sovfl(1:len_trim(sovfl)) // '::' //
     &             snum(1:len_trim(snum)) // ':'
	  slayer(layer-1) = string
	else if (slayer(layer-1) .eq. 'cache') then
          call set2block (ibuffer, iblock, iword, ibyte)
          call i2a (ibuffer, sbuffer)
          call i2a (inum, snum)
          string = slayer(layer-1)(1:len_trim(slayer(layer-1))) //
     &             ':' // sbuffer(1:len_trim(sbuffer)) // 
     &             ':' //snum(1:len_trim(snum))
	  slayer(layer-1) = string
	end if
#ifdef debug_iomngr
          print '(tr1,a,i2,a,a)', 
     &     '              found layer',layer-1,'-->', slayer(layer-1)
#endif
      end if
c
      return
      end

#endif

#ifdef cray_ymp
      function io_blocks (nword)
c
c     number of cray_ymp i/o blocks needed to store nword words
c
      io_blocks = 1 + (nword - 1)/512

      return
      end
#endif

