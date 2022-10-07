!     $Id: ice_indexes.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $

      module ice_indexes

      implicit none

      integer (kind=4), parameter ::
     &   nsnd         = 23  ! number of fields sent to coupler
     &,  nrcv         = 21  ! number of fields sent from coupler

      ! buffer location indexes
      integer (kind=4), parameter ::
        !  ocean states
     &   kbufr_sst     =  1 ! sea surface temp         (K)
     &,  kbufr_sss     =  2 ! sea surface salinity     (o/oo)
     &,  kbufr_uocn    =  3 ! ocean current            (m/s)
     &,  kbufr_vocn    =  4 ! 
        !  atm states
     &,  kbufr_zlvl    =  5 ! atm level height         (m)
     &,  kbufr_uatm    =  6 ! wind                     (m/s)
     &,  kbufr_vatm    =  7 ! 
     &,  kbufr_potT    =  8 ! atm potential T          (K)
     &,  kbufr_Tair    =  9 ! atm temperature          (K) 
     &,  kbufr_Qa      = 10 ! atm specific humidity    (kg/kg)
     &,  kbufr_rhoa    = 11 ! atm air density          (kg/m^3)
        !  more ocean states
     &,  kbufr_tiltx   = 12 ! sea surface slope  
     &,  kbufr_tilty   = 13 !
        !  ocean -> ice flux
     &,  kbufr_fw      = 14 ! ptntl to form/melt ice   (W/m^2)
        !  atm -> ice fluxes
     &,  kbufr_swdidr  = 15 !          near IR,   drct
     &,  kbufr_swdvdr  = 16 ! sw down, vsbl,      drct (W/m^2)
     &,  kbufr_swdidf  = 17 !          near IR,   dffs
     &,  kbufr_swdvdf  = 18 !          vsbl,      dffs
     &,  kbufr_flw     = 19 ! longwave down            (W/m^2)
     &,  kbufr_rain    = 20 ! precip, rain             (kg/m^2 s)
     &,  kbufr_snow    = 21 ! precip, snow             (kg/m^2 s)
        !  ice states
     &,  kbufs_ifrc    =  1 ! ice fraction
     &,  kbufs_tsfc    =  2 ! surface temperature      (K)
     &,  kbufs_alb1    =  3 ! albedo, visible, direct
     &,  kbufs_alb2    =  4 !         near IR, direct
     &,  kbufs_alb3    =  5 !         visible, diffuse
     &,  kbufs_alb4    =  6 !         near IR, diffuse
        !  ice -> atm fluxes
     &,  kbufs_tauxa   =  7 ! wind stress              (N/m^2)
     &,  kbufs_tauya   =  8
     &,  kbufs_lat     =  9 ! latent heat flux         (W/m^2)
     &,  kbufs_sens    = 10 ! sensible heat flux       (W/m^2)
     &,  kbufs_lwup    = 11 ! outgoing longwave radiation (W/m^2)
     &,  kbufs_evap    = 12 ! evaporated water         (kg/m^2 s)
        !  2m atm reference temperature
     &,  kbufs_Tref    = 13 ! (K)
        !  ice -> ocean fluxes
     &,  kbufs_netsw   = 14 ! penetrating shortwave -> ocean (W/m^2)
     &,  kbufs_melth   = 15 ! net ocean heat used      (W/m^2)
     &,  kbufs_meltw   = 16 ! water flux -> ocean      (kg/m^2 s)
     &,  kbufs_salt    = 17 ! salt flux -> ocean       (kg/m^2 s)
     &,  kbufs_tauxo   = 18 ! ice/ocean stress         (N/m^2)
     &,  kbufs_tauyo   = 19
     &,  kbufs_dtauo   = 21 ! ice/ocean stress         (N/m^2)
     &,  kbufs_swabs   = 20 ! absorbed shortwave
     &,  kbufs_we      = 22 ! water eqvivalent (m)
     &,  kbufs_index   = 23 ! a counter

c=======================================================================

      end module ice_indexes

c=======================================================================
