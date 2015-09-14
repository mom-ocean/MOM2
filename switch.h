
c====================== include file "switch.h" ========================
c
c     all time dependent decisions are made by time manager "tmngr.F"
c     and communicated elsewhere to the model via logical switches.
c
c     inputs: (defaulted in "blkdta.F", optionally reset via namelist)
c
c
c     runlen  = integration period (see rununits). note "runlen" should
c               be an integral number of density time steps. if not,
c               then "runlen" is automatically adjusted to insure this.
c               fractional days are supported but not fractional months
c               or years.
c     rununits= units of "runlen". may be "days", "months", or "years".
c               tmngr will convert "runlen" which is in "rununits"
c               to "rundays" in units of days.
c
c     segtim  = the integration time "runlen" is broken into a number of
c               segments each of length "segtim" days. updated surface
c               boundary conditions are applied to MOM every "segtim" 
c               days. this is useful when coupling to atmospheric models
c               in which case both models exchange surface boundary
c               conditions every "segtim" days where "segtim"
c               is 1/(coupling frequency). without an atmospheric model,
c               when getting surface boundary conditions from data, 
c               "segtim" is set to the time step (in days) by mom.F. in
c               either case, "runlen" (in days) should be an integral
c               number of "segtim".     
c
c     nmix    = number of time steps between mixing timesteps. used
c               to damp timestep splitting due to centered leapfrog.
c
c     init    = (true,false)  indicates that this run is a
c               (start from initial conditions, restart)
c
c     restrt  = (true,false) = (do,don`t) write a restart at the end
c               of the run
c
c     eb      = (true,false) configures for the use of a 
c               (euler backward,forward) type mixing timestep
c
c     diagnostic switches (inputs) follow:
c
c     note: switches are used to control the interval between doing
c           diagnostics. units for all switches are in days. 
c           setting a switch < 0.0 disables whatever the switch is 
c           controlling. setting it = 0.0 causes the diagnostic to be 
c           done every time step, and setting it > 0.0 causes the 
c           diagnostic to be done repeatedly on the specified interval. 
c
c     cmixint = number of days between writing estimated mixing coeffs
c               on faces of T cells and U cells
c
c     exconvint = number of days between writing temperature rate of
c                 change due to explicit convection
c
c     glenint =  number of days between global energetics integrals.
c
c     trmbint =  number of days between momentum and tracer term
c                balances (global and regional).
c
c     itrmb   = (true,false) = (do,don`t) write regional mask info for
c               the term balance diagnostic. Typically set true
c               at the beginning of a run; otherwise false since it is
c               not necessary to keep writing a time independent field
c               particularly when it may be a significant part of the
c               time dependent part of the diagnostic.
c
c     gyreint =  number of days between calculation of tracer northward
c                transport.
c
c     vmsfint =  number of days between calculation of vertical and
c                meridional stream function.
c
c     prxzint =  number of days between printouts of x-z data.
c
c     extint  =  number of days between printouts of external mode.
c
c     dspint  =  number of days between surface pressure calculation.
c                Note: only when "stream_function" is enabled.
c
c     tavgint = number of days between regional tracer averages (under
c               horizontal regions).
c
c     itavg   = (true,false) = (do,don`t) write regional mask info for
c               the tracer average diagnostic. Typically set true
c               at the beginning of a run; otherwise false since it is
c               not necessary to keep writing a time independent field
c               particularly when it may be a significant part of the
c               time dependent part of the diagnostic.
c
c     tmbint  = number of days over which tracer equation in averaged
c               in depth and longitude to determine the meridional 
c               balance among storage, divergence, dissipation and 
c               forcing.
c
c     itmb    = (true,false) = (do,don`t) write "msktmb" for tracer
c               the meridional balance diagnostic. Typically set true
c               at the beginning of a run; otherwise false since it is
c               not necessary to keep writing a time independent field
c               particularly when it may be a significant part of the
c               time dependent part of the diagnostic.
c
c
c     tsiint  = number of days between printing of time step integrals.
c
c     stabint = number of days between sampling for various stability
c               criteria.
c
c     snapint = number of days between writing instantaneous samples
c               of data from the MOM grid. note: only when "snapshots"
c               is enabled. see "iounit.h" for more details.
c               all data is at time level "tau" so it is one time step
c               before the date shown in "stamp".
c
c     timavgint= averaging period (days) for writing time mean data from
c               the "averaging" grid (only when "time_averages" is
c               enabled). if "timavgint" is not an integral number of
c               density time steps,"timavgint" is automatically adjusted
c               to insure this. if the number of days to integrate is
c               not an integral number of "timavgint" then the last
c               averaging period will be less than "timavgint" days.this 
c               may lead to one more averaging period than expected. 
c               see "iounit.h" for more details.
c
c     xbtint  = averaging period (days) for writing XBT data (only when
c               "xbts" is enabled). if "xbtint" is not an integral
c               number of density time steps, "xbtint" is automatically
c               adjusted to insure this. if the number of days to 
c               integrate is not an integral number of "xbtint" then the
c               last averaging period will be less than "xbtint" days.
c               this may lead to one more averaging period than  
c               expected. see "iounit.h" for more details.
c
c     zmbcint = number of days between calculation of zonal mean
c               surface boundary conditions (and related  quantities)
c
c     trajint = number of days between writing particle trajectories.
c               particle trajectories are always written on the last
c               time step of the run to allow for restarting.
c
c
c
c
c
c     outputs: (from tmngr.F)
c
c     rundays = integration time in days (from "runlen")
c
c     the following are logical counterparts to the above switches are 
c     set within "tmngr" every time step. logical switches control all
c     decisions about when to do things in MOM.
c
c     cmixts  = (false,true) = (don`t, do) do write estimated mixing
c               coefficients on this time step.
c               based on "cmixint".
c
c     exconvts  = (false,true) = (don`t, do) do write temperature change
c               due to explicit convection on this time step.
c               based on "exconvint".
c
c     glents  = (false,true) = (don`t, do) do calculation of global
c               energy integrals on this time step. based on "glenint".
c
c     trmbts  = (false,true) = (don`t, do) do calculation of momentum &
c               tracer term balance on this timestep. based on "trmbint"
c
c     gyrets  = (false,true) = (don`t, do) do calculation of tracer
c               northward transport on this timestep. based on "gyreint"
c
c     vmsfts  = (false,true) = (don`t, do) do calculation of vertical
c               and meridional stream function on this time step.
c               based on "vmsfint"
c
c     prxzts  = (false,true) = (don`t, do) do printouts of x-z data
c               on this time step. based on "prxzint"
c
c     extts  = (false,true) = (don`t, do) do printout of external mode
c               on this time step. based on "extint"
c
c     dspts  = (false,true) = (don`t, do) do calculation of diagnostic 
c              surface pressure on this time step. based on "dspint"
c               
c
c     stabts  = (false,true) = (don`t, do) test for stability on this
c               time step. based on "stabint"
c
c     tavgts  = (false,true) = (don`t do) do tracer averages on this
c               time step. based on "tavgint"
c
c     tmbts   = (false,true) = (don`t, do) write out tracer meridional .
c               balance on this time step. based on "tmbint"
c
c     tsits   = (false,true) = (don`t, do) print time step integrals
c               on this time step. based on "tsiint"
c
c     zmbcts  = (false,true) = (don`t, do) print zonal mean boundary
c               conditions on this time step.  based on "zmbcint"
c
c     trajts  = (false,true) = (don`t, do) print particle trajectories
c               on this time step. based on "trajint"
c
c     snapts  = (false,true) = (don`t, do) save a snapshot of the data
c               on this time step. based on "snapint"
c
c     timats  = (false,true) = (don`t, do) write time mean data
c               on this time step. based on "timavgint"
c                                          
c     xbtts   = (false,true) = (don`t, do) write averaged XBT data on
c               this time step based on "xbtint"
c
c
c     leapfrog= (false,true) on a (mixing, normal leapfrog) time step
c                based on "nmix"
c                                          
c     euler1  = true on the 1st pass of an euler backward time step
c               otherwise false. (applies when "eb" = true)
c     euler2  = true on the 2nd pass of an euler backward time step
c               otherwise false. (applies when "eb" = true)
c     forward = true on a forward time step. otherwise false
c		(applies when "eb" = false)		  
c                                          
c
c     the following logical switches are based on the model time step. 
c    
c     first   = (true,false) =  when it`s (the first, not the first)
c                               time step of a run
c     eots    = end of a time step. always true except for first
c               pass of an euler backward time step
c     eorun   = last time step of a run. always false except during the
c               last time step of the run.
c
c     eoday   = true when within 1/2 time step of the end of a day
c               else ... false
c     eoweek  = true when within 1/2 time step of the end of a 7 day
c               week (referenced to the start of a year) else ...false
c     eo2wks  = true when within 1/2 time step of the end of two weeks
c               (referenced to the start of a year) else ... false
c     midmon  = true when within 1/2 time step of the middle of a month
c               else ... false
c     eomon   = true when within 1/2 time step of the end of a month
c               else ... false
c     eoyear  = true when within 1/2 time step of the end of a year
c               else ... false
c     osegs   = true on the 1st time step of an ocean segment in mom.F
c               otherwise false.
c     osege  =  true on the last time step of an ocean segment in mom.F
c               otherwise false.
c
c
c
c     author:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c
      logical eb, leapfrog, euler1, euler2, forward, eots
      logical init, first, restrt
      logical itavg, itmb, itrmb
      character *8 rununits
c
      common /switcc/ rununits
      common /switcr/ runlen, rundays
      common /switci/ nmix
      common /switcl/ eb, leapfrog, euler1, euler2, forward, eots
      common /switcl/ init, first, restrt
      common /switcl/ itavg, itmb, itrmb
c
      logical          eoday, eoweek, eo2wks
      common /switci/ ieoday,ieoweek,ieo2wks
      common /switcl/  eoday, eoweek, eo2wks
c
      logical          eomon, midmon, eoyear, eorun
      common /switci/ ieomon,imidmon,ieoyear,ieorun
      common /switcl/  eomon, midmon, eoyear, eorun
c
c-----------------------------------------------------------------------
c     each interval switch needs three variables in common
c         1) a real interval          (e.g.  glenint  )   
c         2) an integer alarm index   (e.g. iglenint  )
c         3) a logical switch         (e.g.  glents)
c
c     user must set the real variable, tmngr takes care of integer and
c     logical counterparts
c-----------------------------------------------------------------------
c
      logical           tavgts,   tmbts,   glents,    trmbts,   prxzts
      common /switcr/   tavgint,  tmbint,  glenint,   trmbint,  prxzint
      common /switci/ iitavgint,iitmbint, iglenint, iitrmbint, iprxzint
      common /switcl/   tavgts,   tmbts,   glents,    trmbts,   prxzts
c
      logical          extts,   trajts,   exconvts,   cmixts
      common /switcr/  extint,  trajint,  exconvint,  cmixint
      common /switci/ iextint, itrajint, iexconvint, icmixint
      common /switcl/  extts,   trajts,   exconvts,   cmixts
c
      logical          vmsfts,   gyrets,   dspts,   tsits,   snapts
      common /switcr/  vmsfint,  gyreint,  dspint,  tsiint,  snapint
      common /switci/ ivmsfint, igyreint, idspint, itsiint, isnapint
      common /switcl/  vmsfts,   gyrets,   dspts,   tsits,   snapts
c
      logical          timats,      xbtts,   stabts,   zmbcts
      common /switcr/  timavgint,   xbtint,  stabint,  zmbcint
      common /switci/ itimavgint, iixbtint, istabint, izmbcint
      common /switcl/  timats,      xbtts,   stabts,   zmbcts
c
      logical          osegs,  osege
      common /switcr/          segtim
      common /switci/         iosege
      common /switcl/  osegs,  osege
