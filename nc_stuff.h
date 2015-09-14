C
c
c            these netcdf routines were borrowed from 
c  author:   john sheldon     e-mail=> jps@gfdl.gov
c            and have not been extensively tested with MOM
c
C   "nc_stuff.h" - May 1994 - JPS
C
C  This include file is intended to help get organized when 
C   developing code for reading/writing netCDF files.  Include it
C   in your code to help make sure all the definitions are set
C   correctly.  Your code must set the values of MXDIMIN, MXDIMOUT, 
C   NDIMS, and NVARS before this include file (defined below).
C
C  For details re: the variables below, set READNC or NCSETUP source 
C......code.  If you use these variables in your code, you'll be able 
C......to interface to READNC/NCSETUP, etc. easily.
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                                                                        !+
C...NETCDF-RELATED STUFF...                                              !+
C   --------------------                                                 !+
C.....Once the vaules of parameters MXDIMIN, MXDIMOUT, NDIMS, and NVARS  !+
C......have been set to accomodate the specific application, the code    !+
C......in the rest of the group delineated by "C++++++" lines does not   !+
C......require any further changes.                                      !+
C         MXDIMIN  : max # of pts along longest INPUT netCDF coord axis  !+
C         MXDIMOUT : max # of pts along longest OUTPUT netCDF coord axis !+
C         NDIMS    : # of different dimensions (axes) to be defined      !+
C                      to OUTPUT file                                    !+
C         NVARS    : # of data variables (NOT coordinate variables)      !+
C                      (a la netCDF parlance) to be sent to OUTPUT file  !+
C                                                                        !+
C -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - !+
C>>...Stuff for both reading and writing netCDF files...                 !+
C ...most of the netCDF definitions...                                   !+
#include "netcdf.h"                                                    !+
C  ...define a maximum length for a character attribute                  !+
      PARAMETER (MAXLATT=256)                                            !+
C                                                                        !+
C -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - !+
C>>...Stuff for READING a netCDF file (via READNC_1 or its derivatives)  !+
C......Use these variables in your code and you'll be able to interface  !+
C......to READNC easily...(see documentation for READNC)...              !+
      INTEGER    CDFID                                                   !+
      LOGICAL    LGETDATA, LSPVA1, LFILLVAL                              !+
C                                                                        !+
      CHARACTER*(MAXNCNAM) FNAMEIN                                        !+
      CHARACTER*(MAXNCNAM) CVARNAM1                                      !+
      CHARACTER*(MAXNCNAM) CXAXNAM, CYAXNAM, CZAXNAM, CTAXNAM            !+
      CHARACTER*(MAXLATT)  CVARLNAM1,CVUNITS                             !+
      CHARACTER*(MAXLATT)  CXUNITS, CYUNITS, CZUNITS, CTUNITS            !+
      CHARACTER*(MAXLATT)  CXAXLNAM,CYAXLNAM,CZAXLNAM,CTAXLNAM           !+
      CHARACTER*(MAXLATT)  GTITLE1                                       !+
C                                                                        !+
      DIMENSION  XCOORV(MXDIMIN), YCOORV(MXDIMIN),                       !+
     *           ZCOORV(MXDIMIN), TCOORV(MXDIMIN)                        !+
C                                                                        !+
C -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - !+
C>>...Stuff for WRITING netCDF file (via NCSETUP/NCSTORE/NCDONE calls).. !+
C......Use these variables in your code and you'll be able to interface  !+
C......to NCSETUP, etc. easily...(see documentation for NCSETUP)...      !+
      CHARACTER*(MAXNCNAM)  FNAME                                      !+
      DIMENSION     LENDIMS(NDIMS)                                       !+
      DIMENSION     DIMVALS(MXDIMOUT,NDIMS)                              !+
      LOGICAL       LCLOBBER, LGSPVAL, LSPVAL(NVARS)                     !+
C                                                                        !+
      CHARACTER*(MAXNCNAM) CDIMNAM(NDIMS),   CVARNAM(NVARS)              !+
      CHARACTER*(MAXLATT)  CDIMLNAM(NDIMS),  CVARLNAM(NVARS)             !+
      CHARACTER*(MAXLATT)  CDIMUNITS(NDIMS), CVARUNITS(NVARS)            !+
      CHARACTER*(MAXLATT)  GTITLE,CDATIM                                 !+ 
C                                                                        !+
      REAL       GSPVAL                                                  !+
      DIMENSION  VALIDR(2,NVARS), SPVAL(NVARS)                          !+
      DIMENSION  NVDIMS(NVARS)                                           !+
C  ....use of NDIMS as inner dim of IDIMINDX is overestimate, but safe.. !+
C  ....it really only needs to be as big as the max # of dims per var... !+
      DIMENSION  IDIMINDX(NDIMS,NVARS)                                   !+
      DIMENSION  ISTART(NDIMS),ICOUNT(NDIMS)                             !+
C  ....handling the record dimesion is sticky; you can output as many    !+
C  .....slices in the record dimension as you want, but typically it's   !+
C  .....one-at-a-time; here, we allow for 200 at a time; mod if needed;  !+
      DIMENSION  RECVAL(200)                                             !+
C                                                                        !+
C                                                                        !+
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

