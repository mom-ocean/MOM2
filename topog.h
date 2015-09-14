c====================== include file "topog.h" =========================
c                                                                       
c                      topography changes             
c                                                                       
c     kmt_opt     = a list of possible modifications to kmt to fix
c                   one flagged problem. 
c     kmt_changes = list of collected changes to kmt
c     nchanges    = number of changes to kmt in kmt_changes
c     auto_kmt_changes = tells whether any changes have been made to the
c                        kmt field as a result define options or
c                        interactive user actions.
c     max_opt     = max number of options for kmt_opt
c     len_opt     = max number of kmt points changed per option
c     io_del_kmt  = io unit for writing delta.kmt.###x###x###.h file
c     n_del_kmt   = number of kmt changes in delta.kmt.###x###x###.h
c
      parameter (max_opt=3, len_opt=10, max_change=100)
      logical auto_kmt_changes
c
      common /kmtchg/ kmt_opt(max_opt, len_opt, 4)
      common /kmtchg/ kmt_changes(max_change, 4)
      common /kmtchg/ nchanges, io_del_kmt, n_del_kmt
      common /kmtchg/ auto_kmt_changes
