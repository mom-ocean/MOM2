c====================== include file "aparam.h" ========================
c
c     parameter file which sets the atmospheric characteristics needed
c     for coupling to MOM
c
c     ima   = number of grid points in the longitudinal direction
c              for the atmospheric S.B.C. grid
c     jma   = number of grid points in the latitudinal direction
c              for the atmospheric S.B.C. grid
c     imap2 = ima + two extra points for the cyclic condition needed
c             for interpolating S.B.C. to MOM
c     
c     coded by:      r. c. pacanowski      e-mail=> rcp@gfdl.gov
c
c
      parameter (ima=36, jma=18, imap2=ima+2)
