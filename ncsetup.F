      SUBROUTINE NCSETUP ( FNAME, LCLOBBER, GTITLE, LGSPVAL,GSPVAL,
     *                     IRECDIMTYP, CDATIM,
     *                     IYY,IMM,IDD,IHH,IMIN,ISEC,
     *                     IYR_REFMIN, IMO_REFMIN,  IDA_REFMIN, 
     *                     IHR_REFMIN, IMIN_REFMIN, ISEC_REFMIN,
     *                     NDIMS, LENDIMS, CDIMNAM, CDIMLNAM,
     *                            CDIMUNITS, DIMVALS, INDIM1,
     *                     NVARS, CVARNAM, CVARLNAM, CVARUNITS, 
     *                            NVDIMS, IDIMINDX, INDIM2,
     *                            VALIDR, LSPVAL, SPVAL )
#ifdef netcdf
c
c            these netcdf routines were borrowed from 
c  author:   john sheldon     e-mail=> jps@gfdl.gov
c            and have not been extensively tested with MOM
c
C-----------------------------------------------------------------------------
C NOTE: for convenience, the length of character strings is indicated numerically  
C       below;  actually, their length is set through PARAMTERS which could
C       conceivably change in the future, so be careful.  The code here checks
C       that the strings are the proper length and STOPS if something doesn't 
C       match.  Use of "Flint" is strongly recommended.
C
C INPUT:
C =====
C   FNAME        : name of netCDF file to be created  (C*128)
C   LCLOBBER     : if TRUE, it is OK to overwrite any file already by that name
C                  if FALSE, and a file already exists by that name, STOP.
C   GTITLE       : a "Global" attribute string  (C*256)
C   LGSPVAL      : logical indicating that a GLOBAL SPVAL should be set
C   GSPVAL       : a value which indicates a missing value
C   IRECDIMTYP   : scheme to use for the "record" dimension, and date/time stamps:
C                  * 1) global attribute CDATIM (string) will be used as 'date';
C                  * 2) date/time will be constructed from GLOBAL attributes
C                        IYY,IMM,IDD,IHH,IMIN,ISEC (someday, also "forecast time");
C                  * 3) record dimension is numeric "time", with or without units;
C                         units may be in either straightforward form (eg, hours) 
C                         or "...since..." form; NOTE: not yet DOUBLE-PRECISION!
C                    4) coordinate variable "record" will contain strings
C                         containing date/time info suitable for labeling a plot
C                           (need to experiment to see if this will work)
C                    5) GLOBAL attribute "record" will contain the names of 
C                         variables ("reftime,valtime") which can be used to
C                         construct both the valid date/time and the forecast time
C                         (this routine will define the vars "reftime" and "valtime"
C                          to be DOUBLE, assuming they will use the "___ since..."
C                          approach)
C                    6) if a variable references the 'record' dimension, and it
C                         has attached to it an attribute "record", the attribute
C                         "record" will contain the names of variables (typically
C                          "reftime,valtime") which can be used to construct both 
C                          the valid date/time and the forecast time
C                ( * = currently supported )
C
C   CDATIM       : a character string containing a desired date/time stamp which 
C                    applies to the whole file (see option #1 above); (C*256)
C                        NOTE: if null string is sent in, CDATIM will not be
C                              written to the file, avoiding later confusion
C
C   IYY,IMM,IDD,IHH,IMIN,ISEC : year, month, day, hour, minute, second for entire
C                                 file, from which a date/time label could be 
C                                 constructed (see option 2 above)
C  (RV_TIME_UNITS: T.B.D. - not yet used)
C   I{*}_REFMIN  : a BASE reference date/time; all subsequent derived dates/times
C                     (using "valtime" and "reftime") will use this in stating
C                    'time' as "seconds since....", for example (options 4 and 5)
C   NDIMS         : total number of coordinate axes used in all the data arrays
C                     to be stored in this file
C   LENDIMS(*)   : a array NDIMS long containing the length of all the axes;
C                    a "record" dimension, if any, is indicated by LENDIMS=0 for
C                    that dimension
C   CDIMNAM(*)   : a character*128 array of length NDIMS containing the names 
C                    of each of the axes (FORTRAN-like names, no spaces or wierd 
C                   characters))
C   CDIMLNAM(*)  : a character*256 array of length NDIMS containing the  
C                    "long-names" (ie, descriptive strings) of each of the axes 
C                    (arbitrary strings OK)
C   CDIMUNITS(*) : a character*256 array of length NDIMS containing the units 
C                    of each of the axes (arbitrary strings OK, but more useful if
C                   in a form compatible with UDUNITS, a Unidata units-handling
C                   library));  if you are defining a "record dimension" 
C                    for use with the options for IRECDIMTYP involving
C                   the "record" dimension which references "reftime"
C                   and "valtime", this should be the units that 
C                   the "....since..." string begins with; eg, set 
C                   CDIMUNITS(4)='hours' if dimension #4 is the record dimension and
C                   "reftime" and valtime" will be measured in "hours since...."
C   DIMVALS(*,*) : array dimensioned (INDIM1,NDIMS) containing the coordinates
C                    of the points along each dimension; 
C                        eg, DIMVALS( 1:LENDIMS(i) ,i ) are the coordinates along
C                                    i-th dimension
C                    (eventually, we'll have to figure out how to get the time
C                      dimension to be DOUBLE PRECISION, a la NUWG convention;
C                      probably add new args)
C   INDIM1       : inner dimension of DIMVALS
C
C   NVARS        : total number of data variables to be stored in this file (MAX 512)
C                        IMPORTANT: DO NOT INCLUDE THE AXES IN THE "VARIABLE" COUNT...!
C                                  (even tho they are, in netCDF, treated as vars)
C   CVARNAM(*)   : character*128 array of length NVARS containing the names of 
C                     variables (FORTRAN-like name, ie, alphanumeric and '_' OK)
C   CVARLNAM(*)  : character*256 array of length NVARS containing the long-names of 
C                     the variables (arbitrary string)
C   CVARUNITS(*) : character*256 array of length NVARS containing the units for 
C                     the variables (arbitrary string, BUT best to follow UDUNITS
C                    conventions)
C   NVDIMS(*)    : array of length NVARS containing the number of dimensions 
C                     for each variable
C   IDIMINDX(*,*): array of size (INDIM2,NVARS) containing the indices of the
C                     axes above which apply to each variable;
C                      eg, if IDIMINDX(2,4)=5, this means that the second dimension
C                          of the 4th data variable is represented by axis #5  
C                          defined above
C   INDIM2       : inner dimension of array IDIMINDX
C   VALIDR(2,*)  : array of size (2,NVARS) containing the valid range (min,max)
C                      of data values for each variable
C   LSPVAL(*)    : array of size NVARS containing indicating whether there is a
C                     "SPVAL" for each variable
C   SPVAL(*)     : array of size NVARS containing a value which indicates a
C                      missing value (this overrides the GLOBAL GSPVAL)
C
C
C LIMITATIONS:
C ===========
C     1) Data to be stored must be single-precision floating point;
C     2) coordinate axis values are SINGLE precision; someday, when better 
C         support for the "record"/"time" dimension is installed, "time" 
C         will have to be DOUBLE precision;
C
C===============================================================================

C..#include "/net/jps/unidata.ucar.edu/pub/netcdf/include/netcdf.inc"
#include "netcdf.h"
      PARAMETER (MAXLATT=256)

C...dummy args...
	CHARACTER*(MAXNCNAM)  FNAME
	CHARACTER*(MAXLATT)   GTITLE,CDATIM

	LOGICAL        LCLOBBER, LGSPVAL, LSPVAL(NVARS)
	DIMENSION      LENDIMS(NDIMS)
	REAL           DIMVALS(INDIM1,NDIMS)
	
	CHARACTER*(MAXNCNAM)  CDIMNAM(NDIMS),   CVARNAM(NVARS)
	CHARACTER*(MAXLATT)  CDIMLNAM(NDIMS),  CVARLNAM(NVARS)
	CHARACTER*(MAXLATT) CDIMUNITS(NDIMS), CVARUNITS(NVARS)

	DIMENSION  VALIDR(2,NVARS), SPVAL(NVARS)
	DIMENSION  NVDIMS(NVARS), IDIMINDX(INDIM2,NVARS)

C...local stuff....
	CHARACTER*64  TIME_UNITS                    ! for future use
	DIMENSION IDIMS(4), ISTART(4), ICOUNT(4)
	LOGICAL LEXIST
	CHARACTER*16 CTMP16
	CHARACTER*(MAXNCNAM) CTMPNAM
	CHARACTER*(MAXLATT)  CTMPATT
	
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C #include "ncwrapper.h"
	PARAMETER (MAXNDIMS=MAXNCDIM, MAXNVARS=MAXNCVAR)
	INTEGER CDFID
	COMMON /NCIDSTUFF/ CDFID, NDIMS1, ID_DIM(MAXNDIMS), 
     *                            NVARS1, ID_VAR(MAXNVARS),
     *			          IRECDIM
C     *                           ,ID_RTIME, ID_VTIME       ! for future use
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	
	INTEGER CLEN
	EXTERNAL CLEN


C================================================================
C23456789_C23456789_C23456789_C23456789_C23456789_C23456789_C23456789_12
C================================================================

C...First, check for things which will muck everything up...
C  ...max # dims for netCDF...
      IF (NDIMS.GT.MAXNCDIM) THEN
          PRINT *,'In NCSETUP, NDIMS too large (=',NDIMS,'), max=',
     *			MAXNCDIM
          STOP '04111028'
      ENDIF
C  ..."indefinite" values (check SOME); if user uses indef initialization, 
C  .....we can force/detect problems here by doing simple math....
      IF (LGSPVAL) TMP=GSPVAL/2.
      TMP=IYY/2
      TMP=IDD/2
      TMP=IHH/2
      TMP=IMIN/2
      TMP=ISEC/2
      TMP=NDIMS/2
      TMP=INDIM1/2
      DO L2=1,NDIMS
         TMP=LENDIMS(L2)/2
         DO L1=1,LENDIMS(L2)
	    TMP=DIMVALS(L1,L2)/2.
	 ENDDO
      ENDDO
      TMP=NVARS/2
      TMP=INDIM2/2
      DO L2=1,NVARS
         TMP=NVDIMS(L2)/2
         DO L1=1,NVDIMS(L2)
	    TMP=IDIMINDX(L1,L2)/2
	 ENDDO
      ENDDO
      IF (TMP.EQ.0) PRINT *,' '
C...strings not the correct length.......
      NSTOP=0
      IF (LEN(FNAME).NE.MAXNCNAM)    NSTOP=NSTOP+1
      IF (LEN(GTITLE).NE.MAXLATT)    NSTOP=NSTOP+1
      IF (LEN(CDATIM).NE.MAXLATT)    NSTOP=NSTOP+1
      DO I=1,NDIMS
         IF (LEN(CDIMNAM(I)).NE.MAXNCNAM)  NSTOP=NSTOP+1
         IF (LEN(CDIMLNAM(I)).NE.MAXLATT)  NSTOP=NSTOP+1
         IF (LEN(CDIMUNITS(I)).NE.MAXLATT) NSTOP=NSTOP+1
      ENDDO
      DO I=1,NVARS
         IF (LEN(CVARNAM(I)).NE.MAXNCNAM)  NSTOP=NSTOP+1
         IF (LEN(CVARLNAM(I)).NE.MAXLATT)  NSTOP=NSTOP+1
         IF (LEN(CVARUNITS(I)).NE.MAXLATT) NSTOP=NSTOP+1
      ENDDO
      IF (NSTOP.GT.0) THEN
          PRINT *,'In NCSETUP, ',NSTOP,' strings are wrong length.'
	  STOP '04121016'
      ENDIF
C  ...invalid characters in variable names...
      DO I=1,NDIMS
          CTMPNAM=CDIMNAM(I)
          LL=CLEN(CTMPNAM)
	  NG=0
	  DO L=1,LL
             IC=ICHAR(CTMPNAM(L:L))
	     IF (IC.LT.48)                NG=NG+1
	     IF (IC.GE.58 .AND. IC.LE.64) NG=NG+1
	     IF (IC.GE.91 .AND. IC.LE.94) NG=NG+1
	     IF (IC.EQ.96)                NG=NG+1
	     IF (IC.GT.122)               NG=NG+1
	  ENDDO
	  IF (NG.GT.0) THEN
             PRINT *,' '
	     PRINT *,'[1m Dimension name #',I,' contains ',NG,
     *               ' invalid characters:'
             PRINT *,'"',CTMPNAM(1:LL),'" [21m'
	     STOP '03172348'
	  ENDIF
      ENDDO
      DO I=1,NVARS
          CTMPNAM=CVARNAM(I)
          LL=CLEN(CTMPNAM)
	  NG=0
	  DO L=1,LL
             IC=ICHAR(CTMPNAM(L:L))
	     IF (IC.LT.48)                NG=NG+1
	     IF (IC.GE.58 .AND. IC.LE.64) NG=NG+1
	     IF (IC.GE.91 .AND. IC.LE.94) NG=NG+1
	     IF (IC.EQ.96)                NG=NG+1
	     IF (IC.GT.122)               NG=NG+1
	  ENDDO
	  IF (NG.GT.0) THEN
             PRINT *,' '
	     PRINT *,'[1m Variable name #',I,' contains ',NG,
     *               ' invalid characters:'
             PRINT *,'"',CTMPNAM(1:LL),'" [21m'
	     STOP '03172349'
	  ENDIF
      ENDDO


C...set it so that any illegal netCDF action will force an abort....
       CALL NCPOPT (NCVERBOS+NCFATAL) 
C      CALL NCPOPT (NCVERBOS) 


C...Open the output file, checking if LCLOBBER=F...
	print *,'Entering netCDF DEFINE mode...'
	IF (LCLOBBER) THEN
	    LL=CLEN(FNAME)
	    PRINT *,'Opening file [1m',FNAME(1:LL),'[21m...'
	    CDFID = NCCRE (FNAME, NCCLOB, IRET)	   ! write over any existing file...
            IF (IRET.NE.0) STOP '03211315'
	ELSE
            INQUIRE(FILE=FNAME,EXIST=LEXIST)
	    IF (LEXIST) THEN
	        PRINT *,'File already exists:'
		PRINT *,'[1M ',FNAME,' [21M'
		STOP '12221300'
	    ELSE
	        CDFID = NCCRE (FNAME, NCNOCLOB, IRET)
                IF (IRET.NE.0) STOP '03211316'
	    ENDIF
	ENDIF

C...save the values of NDIMS and NVARS for later use....
        NDIMS1=NDIMS
	NVARS1=NVARS
	PRINT *,'File will have ',NDIMS1,' dimensions and ',
     *				  NVARS1,' variables.'


C..Assign global attributes to this file.....
	print *,'  ...global attributes..'
	print *,'...title'
	CTMP16='title           '
	CALL NCAPTC(CDFID, NCGLOBAL, CTMP16, NCCHAR, CLEN(GTITLE),
     *						        GTITLE, IRET)
        IF (IRET.NE.0) STOP '03181835'

C  ...if the user sends in a CDATIM string which is not null, write it out...
	LL=CLEN(CDATIM)
	IF (LL.GT.0) THEN
	   print *,'...date_time'

C   ...check for illegal characters...
             LL=LEN(CDATIM)
	     NG=0
	     DO L=1,LL
                IC=ICHAR(CDATIM(L:L))
	        IF (IC.LT.32)                NG=NG+1
	        IF (IC.GE.33 .AND. IC.LE.47) NG=NG+1
	        IF (IC.GE.58 .AND. IC.LE.64) NG=NG+1
	        IF (IC.GE.91 .AND. IC.LE.94) NG=NG+1
	        IF (IC.EQ.96)                NG=NG+1
	        IF (IC.GT.122)               NG=NG+1
	     ENDDO
	     IF (NG.GT.0) THEN
                PRINT *,' '
	        PRINT *,'[1m CDATIM contains ',NG,
     *                  ' invalid characters:'
                PRINT *,'"',CDATIM,'" [21m'
	        STOP '03191156'
	     ENDIF

	   CTMP16='date_time       '
	   CALL NCAPTC(CDFID, NCGLOBAL,CTMP16,NCCHAR,CLEN(CDATIM),
     *                                                   CDATIM,IRET)
           IF (IRET.NE.0) STOP '03181836'
     	ENDIF

C  ...does the user want to specify a "Special-value" to apply GLOBAL-ly?...
	IF (LGSPVAL) THEN
	   print *,'missing-value'
	   CTMP16='missing-value   '
	   CALL NCAPT(CDFID, NCGLOBAL, CTMP16, NCFLOAT,
     *                                                1,GSPVAL,IRET)
           IF (IRET.NE.0) STOP '03181837'
     	ENDIF

C  ...in any case, put out a _FillValue; for now, we just use the default
C  ....netCDF fill value NCFILL, but it will be informative to have this
C  ....show up as a global attribute in, eg, an ncdump...
C? **this approach aborted - actual fill value not available in FORTRAN!
c           IOLDMODE=NCSFIL(CDFID,NCFILL,IRET)
c	   CTMP16='_FillValue      '
c	   CALL NCAPT(CDFID, NCGLOBAL, CTMP16, NCFLOAT,
c     *                                              1,FILL_FLOAT,IRET)


C  ...assume the user sent in something for IYY,...etc., even if zeros.........	
	print *,'...year=',IYY
	   CTMP16='year            '
	CALL NCAPT(CDFID, NCGLOBAL, CTMP16, NCLONG, 1, IYY, IRET)
        IF (IRET.NE.0) STOP '03181838'
	print *,'...month=',IMM
	   CTMP16='month           '
	CALL NCAPT(CDFID, NCGLOBAL, CTMP16, NCLONG, 1, IMM, IRET)
        IF (IRET.NE.0) STOP '03181839'
	print *,'...day=',IDD
	   CTMP16='day             '
	   IONE=1
C	print *,'   args: CDFID,NCGLOBAL,CTMP16,NCLONG,IONE,IDD,IRET'
C	print *,'       ',CDFID,NCGLOBAL,CTMP16,NCLONG,IONE,IDD,IRET
	CALL NCAPT(CDFID, NCGLOBAL, CTMP16, NCLONG, IONE, IDD, IRET)
c23456789_c23456789_c23456789_c23456789_c23456789_c23456789_c23456789_12
C	print *,' on rtn: CDFID,NCGLOBAL,CTMP16,NCLONG,IONE,IDD,IRET'
C	print *,'       ',CDFID,NCGLOBAL,CTMP16,NCLONG,IONE,IDD,IRET
        IF (IRET.NE.0) STOP '03181840'
	print *,'...hour=',IHH
	   CTMP16='hour            '
	CALL NCAPT(CDFID, NCGLOBAL, CTMP16, NCLONG, 1, IHH, IRET)
        IF (IRET.NE.0) STOP '03181841'
	print *,'...minute=',IMIN
	   CTMP16='minute          '
	CALL NCAPT(CDFID, NCGLOBAL, CTMP16, NCLONG, 1, IMIN, IRET)
        IF (IRET.NE.0) STOP '03181842'
	print *,'...second=',ISEC
	   CTMP16='second          '
	CALL NCAPT(CDFID, NCGLOBAL, CTMP16, NCLONG, 1, ISEC, IRET)
        IF (IRET.NE.0) STOP '03181843'
	
	


C..Tell netCDF about the various Dimensions we have and get a dimension ID for each...
C...-these ID's get associated with any data variables using that dimension......
C...NOTE: it is important that the LABELS here match the label used when storing
C..........the actual X/Y/Z coordinate values...........
C......For this file, then:
C........1) label each dimension,
C........2) tell how big it is,
C........3) get an "ID" for it.
C.....If a dimen. length=0, ===> this is the "record" (UNLIMITED) dimension...
	print *,'   ....defining dimensions...'
	IRECDIM=-1
	DO I=1,NDIMS
	   IF (LENDIMS(I).EQ.0) THEN
	    IF (IRECDIM.EQ.-1) THEN
	     IRECDIM=I
             print *, 'Calling NCDDEF (LEN=0) ==>',cdfid, cdimnam(I)
     *,               ncunlim, iret
	     ID_DIM(I)=NCDDEF(CDFID,CDIMNAM(I),NCUNLIM,IRET)
            ELSE
	     PRINT *,'NCSETUP Error: record dimension already defined.'
	     PRINT *,'LENDIMS(*):',(LENDIMS(II),II=1,NDIMS)
	     STOP '12201105'
	    ENDIF
	   ELSE
              print *, 'Calling NCDDEF ==>',cdfid, cdimnam(I)
     *,               lendims(I), iret
	      ID_DIM(I)=NCDDEF(CDFID, CDIMNAM(I), LENDIMS(I), IRET)
              print *, 'ID_DIM= ',ID_DIM(I)
	   ENDIF
	ENDDO
	IF (IRECDIM.NE.-1) PRINT *,' IRECDIM defined, = ',IRECDIM


C..Next, tell netCDF what Variables we have and get ID's for them...
C  ....
C  ....The first "variables" to define are the axes themselves....
C  ........1) Label the axis array,
C  ........2) tell what kind of number is used as coord. pts (ie, integer, float,...)
C  ........3) this is a 1 dimensional "variable",
C  ........4) use a dimension-ID to point to which "Dimension" this "variable" uses...
	print *,'  ...Describing dimension variables...'
	DO I=1,NDIMS
	   ID_VAR(I) = NCVDEF(CDFID, CDIMNAM(I), NCFLOAT, 1, 
     *						 ID_DIM(I), IRET)
           IF (IRET.NE.0) STOP '03172317'
	ENDDO

C  ......Assign "attributes" to these "variables"...
	print *,'  ...attributes for dimensions....'
	DO I=1,NDIMS
	   CTMP16='long_name       '
	   CTMPATT=CDIMLNAM(I)
	   CALL NCAPTC(CDFID, ID_VAR(I), CTMP16, NCCHAR, 
     *				 MAXLATT, CTMPATT, IRET)
	   CTMP16='units           '
	   CTMPATT=CDIMUNITS(I)
	   CALL NCAPTC(CDFID, ID_VAR(I), CTMP16, NCCHAR, 
     *				 MAXLATT, CTMPATT, IRET)
	ENDDO


      IF (IRECDIMTYP.EQ.3 .OR. IRECDIMTYP.GE.5) THEN
C     ...check for legal IRECDIM...
	 IF (IRECDIM.LT.1 .OR. IRECDIM.GT.NDIMS) THEN
	    PRINT *,'In NCSETUP, IRECDIM not legal, =',IRECDIM
	    STOP '12201124'
	 ENDIF
	 IF (IRECDIMTYP.GE.5) THEN
C       ...construct a string of the form "_____ since ..."...
	    CTMPATT=CDIMUNITS(IRECDIM)
            LL=CLEN(CTMPATT)
	    WRITE (TIME_UNITS,9009) CDIMUNITS(IRECDIM)(1:LL),
     *                           IYR_REFMIN, IMO_REFMIN,  IDA_REFMIN, 
     *                           IHR_REFMIN, IMIN_REFMIN, ISEC_REFMIN
 9009	    FORMAT ( A,' since (',I4.4,'-',I2.2,'-',I2.2,' ',
     *				I2.2,':',I2.2,':',I2.2,'.0)' )
         ELSE
	    TIME_UNITS=CDIMUNITS(IRECDIM)(1:64)
	 ENDIF
         PRINT *,'In NCSETUP, TIME_UNITS ='
	 PRINT *,'     ',TIME_UNITS
      ENDIF


      IF (IRECDIMTYP.GE.5) THEN
      
         print *,' From NCSETUP: IRECDIMTYP=5 not yet supported.'
	 stop '12301408'
      
C   ...Now, the "Quasi-Dimensional" variables;
C   .....If we adopt the NUWG convention regarding the time dimension of a data
C   ......set, the UNLIMITED dimension is "record" (ie, 1,2,3,...), each of which 
C   ......contains a 3-D snapshot of the data.  Ancillary variables "reftime" and
C   ......"valtime" utilize the "record" dimension; ie, every time a 3-D snapshot
C   ......is stored, a values of "reftime" and "valtime" are stored. Variable
C   ......"reftime" is a reference time (ie, an IC time) and "valtime" is a
C   ......validation time, allowing many forecasts from several initial times
C   ......to be stored in the same data set. Both are measured, for example, in
C   ......hours from a base time. The "unit" attributes of "reftime" and "valtime"
C   ......are expressed such that the UDUNITS package can determine the acutal
C   ......year, month, day,... values...............
c1	 print *,'  ...Describing time-reference variables..'
c1	 ID_RTIME = NCVDEF (CDFID, 'reftime' , NCDOUBLE, 1, IDREC, IRET)
c1	 ID_VTIME = NCVDEF (CDFID, 'valtime' , NCDOUBLE, 1, IDREC, IRET)
c1C   ...Assign "attributes" to these "variables"...
c1	 print *,'  ...attributes for time-reference dimensions....'
c1	 CTMP16='long_name       '
c1	 CALL NCAPTC(CDFID, ID_RTIME, CTMP16, NCCHAR, 14,
c1     *						'reference time', IRET)
c1	 CTMP16='long_name       '
c1	 CALL NCAPTC(CDFID, ID_VTIME, CTMP16, NCCHAR, 10,
c1     *						'valid time', IRET)
c1	 CTMP16='units           '
c1	 CALL NCAPTC(CDFID, ID_RTIME, CTMP16, NCCHAR, CLEN(TIME_UNITS),
c1     *					              TIME_UNITS , IRET)
c1	 CTMP16='units           '
c1	 CALL NCAPTC(CDFID, ID_VTIME, CTMP16, NCCHAR, CLEN(TIME_UNITS),
c1     *					              TIME_UNITS , IRET)
      ENDIF



C...Now describe the DATA variables...
      print *,'  ....Describing data variables....'
      DO IV=1,NVARS
      
        LL=CLEN(CVARNAM(IV))
	print *,' '
        print *,'Var #',IV,'=',CVARNAM(IV)(1:LL)
	print *,'  #dims=',NVDIMS(IV)
      
         II = IV+NDIMS

C      ...First, load IDIMS with the ID's of the pertinent dimensions for this variable..
         ND=NVDIMS(IV)
	 DO J=1,ND
	    IDIMS(J)=ID_DIM(IDIMINDX(J,IV))
	    
c	    print *,'  J,IV=',J,IV
c	    print *,'  IDIMINDX(J,IV)=',IDIMINDX(J,IV)
	 print *,
     $     '     Index of dim #',J,' (ie, IDIMS(',J,') )=',IDIMS(J)
	    
	 ENDDO

C     ...Now, tell netCDF what "name" and dimensions to associate with this array..
C       ...1) the label is CVARNAM(IV) (eg, "Z"),
C       ...2) the data is floating point,
C       ...3) there are NDIMS dimensions,
C       ...4) array IDIMS contains the ID's of the dimensions;
c        print *,'  ...',CVARNAM(IV)
c	 print *,'IDIMS:',idims
	 
c	 print *,'CDFID, CVARNAM(IV), NCFLOAT, ND, IDIMS:'
c	 print *, CDFID
c	 print *, CVARNAM(IV)
c	 print *, NCFLOAT
c	 print *, ND
c	 print *, IDIMS
	 
	 print *,'  CVARLNAM(IV)=',CVARLNAM(IV)
	 print *,'  CVARUNITS(IV)=',CVARUNITS(IV)
	 print *,'  VALIDR(1,IV)=',VALIDR(1,IV)
	 print *,'  VALIDR(2,IV)=',VALIDR(2,IV)
	 print *,'  LSPVAL(IV)=',LSPVAL(IV)
	 
	 ID_VAR(II) = NCVDEF(CDFID, CVARNAM(IV), NCFLOAT,ND,IDIMS, IRET)
	 IF (IRET.NE.0) STOP '03172332'
	
C     ...Assign whatever attributes are appropriate for this data....
	 CTMP16='long_name       '
	 CTMPATT=CVARLNAM(IV)
	 CALL NCAPTC(CDFID, ID_VAR(II), CTMP16, NCCHAR, 
     *				MAXLATT, CTMPATT, IRET)
	 CTMP16='units           '
	 CTMPATT=CVARUNITS(IV)
	 CALL NCAPTC(CDFID, ID_VAR(II), CTMP16, NCCHAR,  
     *				MAXLATT, CTMPATT, IRET)
	 CTMP16='valid_range     '
	 CALL NCAPT(CDFID, ID_VAR(II), CTMP16, NCFLOAT, 2,  
     *						    VALIDR(1,IV), IRET)
	 IF (LSPVAL(IV)) CALL NCAPT(CDFID, ID_VAR(II), 'missing_value', 
     *					   NCFLOAT, 1, SPVAL(IV), IRET)

      ENDDO


C...Leave netCDF DEFINE mode......
     	print *,'Leaving DEFINE mode....',CDFID
	CALL NCENDF(CDFID,IRET)



C  ..Put out the non-record dimensions....
      print *,'Putting non-record dimension values out....'
      DO I=1,NDIMS
         IF (I.NE.IRECDIM) THEN
            ISTART(1)=1
            ICOUNT(1)=LENDIMS(I)
            CALL NCVPT(CDFID, ID_DIM(I), ISTART,ICOUNT, 
     *					DIMVALS(1,I), IRET)
         ENDIF
      ENDDO	

      PRINT *,'-----------------'
      PRINT *,'NCSETUP complete.'
      PRINT *,'-----------------'

      RETURN
      END




 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER FUNCTION CLEN (S) 
C...Translated by Pacific-Sierra Research VAST-90 1.02V   17:06:41  10/15/93  
C...Switches: -rr -rk -np -nq -nc -nb RENUMB=1000:100,FORMAT=9000:5,CO
C...Switches: ONCHR=*,INDAL=3
C**************************************************************                 
C                                                             *                 
C SUBROUTINES AND FUNCTIONS CALLED FROM THIS ROUTINE          *                 
C                                                             *                 
C BEEP                                                        *                 
C                                                             *                 
C**************************************************************                 
C Function to find length of string, without trailing blanks;
C  Looks thru entire string (up to 256 chars) and finds last non-blank char.
      CHARACTER*(*) S 
      CHARACTER*256 S1 
 
c     CLEN=LEN(S)
c     RETURN
 
      ILEN = LEN(S) 
      IF (ILEN.GT.256) THEN 
         WRITE (15,*) 'String length>256 in CLEN;', 
     *      ' Checking only first 256 chars.' 
c         CALL BEEP 
      ENDIF 
      ILEN = MIN(ILEN,256) 
      S1(1:ILEN) = S(1:ILEN) 
 
c.    write (6,*) 'In CLEN,S =*',S,'*'
c.    write (6,*) '    LEN(S)=',len(s)
c.    write (6,*) 'In CLEN,S1=*',S1,'*'
c.    write (6,*) '    LEN(S1)=',len(s1)
C     ILEN=LEN(S1)
c.    CLEN=0
      CLEN = ILEN 
C.    DO 10 I=1,ILEN
      DO I = ILEN,1,-1 
c       write (6,*) 'In loop 10, I, S1(I:I)=',I,'  "', S1(I:I),'"'
c         IF (S1(I:I).NE.' ' .and. S1(I:I).NE.'') THEN 
c         IF (S1(I:I).NE.' ') THEN 
         IF (ichar(S1(I:I)).NE.32 .and. ichar(S1(I:I)).NE.0) THEN 
            CLEN = I 
            GO TO 1000 
C       write(6,*)'I, S1(I:I)=',I,'  "',S1(I:I),'"'
         ENDIF 
      END DO 
      CLEN = 0
 
 1000 CONTINUE 
      RETURN  
      END  

C========================================================================
C========================================================================
C========================================================================
      SUBROUTINE NCSTORE(IVAR,ISTART,ICOUNT,ARR,RECVAL)

C
C This routine performs the actual write to the netCDF file, as set up by
C  subroutine NCSETUP.
C
C  INPUT:
C  -----
C     IVAR     : index of variable to be stored, consistent with how it was
C                  defined in NCSETUP call;
C     ISTART   : array of starting indices in the disk-stored array at which
C                  to begin storing data;
C     ICOUNT   : array containing the number of points to store along each axis
C     ARR      : array of data to be stored;
C     RECVAL   : "record" axis coordinate values for incoming data, if
C                  applicable (used only if data variable uses "record"
C                  dimension)
C
C  NOTE: array ARR can be an arbitrary "hyperslab" of data representing the 
C         entire array, a subvolue of the entire array, or all or part of an
C         orthogonal slice through the data volume.  For example, if the data
C         variable were defined to NCSETUP as being dimensioned (30,20,10), 
C         calling NCSTORE with ISTART(*)=2,2,5 and ICOUNT(*)=28,18,1 would 
C         store the fifth slice along the Z-axis, except for the edge rows 
C         and columns.
C     
C LIMITATIONS:
C -----------
C     1. data and, if applicable, "record" axis coordinate value(s) are assumed 
C          to be real, single-precision;
C     2. for now, the "record" axis coordinate values are output to the coordinate
C          variable in the same order as they are received by ncstore; no attempt
C          is made to "re-order" the data along the "record" axis.
C---------------------------------------------------------------------------------
C...dummy args....
      DIMENSION ISTART(*),ICOUNT(*),ARR(*),RECVAL(*)


C..#include "/net/jps/unidata.ucar.edu/pub/netcdf/include/netcdf.inc"
#include "netcdf.h"

      CHARACTER*(MAXNCNAM)  CVARNAM,CDIMNAM
      DIMENSION  IVDIMS(MAXNCDIM)


c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C #include "ncwrapper.h"
	PARAMETER (MAXNDIMS=MAXNCDIM, MAXNVARS=MAXNCVAR)
	INTEGER CDFID
        COMMON /NCIDSTUFF/ CDFID, NDIMS1, ID_DIM(MAXNDIMS), 
     *                            NVARS1, ID_VAR(MAXNVARS),
     *			          IRECDIM
C     *                           ,ID_RTIME, ID_VTIME       ! for future use
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C==========================================================================

      IF (CDFID.LT.0) THEN
         PRINT *,'In NCSTORE: netCDF file not open yet (CDFID<0).'
	 STOP '12301412'
      ENDIF

      IF (IVAR.GT.NVARS1) THEN
         PRINT *,'NCSTORE Error: IVAR=',IVAR,' > max. # of variables'
	 STOP '12281739'
      ENDIF
      
      CALL NCVINQ (CDFID, ID_VAR(NDIMS1+IVAR), CVARNAM, IVTYPE, 
     *               NVDIMS, IVDIMS, NVATTS, IRET)

C...check to make sure user is not trying to store more than allowed...
C...IVRI will be the dim. number which refers to the record index, if any...
      IVRI=0
c      write(*,*) 'ICOUNT =',  (icount(i),i=1,nvdims)
      DO L=1,NVDIMS
          CALL NCDINQ(CDFID,IVDIMS(L),CDIMNAM,IDIMSIZ, IRET)
	  ILAST=ISTART(L)+ICOUNT(L)-1
          IF (ISTART(L).LE.0) THEN
	     PRINT *,'NCSTORE Error: bad ISTART=',ISTART(L),
     *		    ' (ie, less than 1) '
             STOP '12281814'
	  ENDIF
          IF (ICOUNT(L).LE.0) THEN
	     PRINT *,'NCSTORE Error: bad ICOUNT=',ICOUNT(L),
     *		    ' (ie, less than zero) '
             STOP '12281812'
	  ENDIF
          IF (IVDIMS(L).NE.IRECDIM .AND. ISTART(L).GT.IDIMSIZ) THEN
	     PRINT *,'NCSTORE Error: ISTART=',ISTART(L),
     *		   ' values along dim that is only ',IDIMSIZ,' long.'
             STOP '12281815'
	  ENDIF
          IF (IVDIMS(L).NE.IRECDIM .AND. ILAST.GT.IDIMSIZ) THEN
	     PRINT *,'NCSTORE Error: trying to store ',ILAST,
     *		   'values along dim that is only ',IDIMSIZ,' long.'
             STOP '12281813'
	  ENDIF
C       ...save the index associated with the record coordinate if
C       ....data uses the "record" dimension...
          IF (IVDIMS(L).EQ.IRECDIM .AND. IVRI.EQ.0) IVRI=L
      ENDDO

C...if the data variable utilizes the "record" dimension, store the coordinate
C....values along the record axis which are associated with this data...
       IF (IVRI.NE.0) CALL NCVPT(CDFID, ID_VAR(IRECDIM), 
     *			      ISTART(IVRI),ICOUNT(IVRI),RECVAL, IRET)

C...store the data.....
      CALL NCVPT(CDFID, ID_VAR(NDIMS1+IVAR), ISTART,ICOUNT,ARR, IRET)

      RETURN
      END


C========================================================================
C========================================================================
C========================================================================
      SUBROUTINE NCDONE
C..#include "/net/jps/unidata.ucar.edu/pub/netcdf/include/netcdf.inc"
#include "netcdf.h"
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C #include "ncwrapper.h"
	PARAMETER (MAXNDIMS=MAXNCDIM, MAXNVARS=MAXNCVAR)
	INTEGER CDFID
	COMMON /NCIDSTUFF/ CDFID, NDIMS1, ID_DIM(MAXNDIMS), 
     *                            NVARS1, ID_VAR(MAXNVARS),
     *			          IRECDIM
C     *                           ,ID_RTIME, ID_VTIME       ! for future use
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C==========================================================================

      IF (CDFID.LT.0) THEN
         PRINT *,'In NCDONE: netCDF file not open yet (CDFID<0).'
	 STOP '12301423'
      ENDIF

      CALL NCCLOS(CDFID,IRET)
      CDFID=-1
      NDIMS1=0
      NVARS1=0
       
      END
#else
      return
      end
#endif
