!     $Id: atm_indexes.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $

      module atm_indexes

      implicit none

      integer (kind=4), parameter ::
     &   nsnd         = 16  ! number of fields sent to coupler
     &,  nrcv         =  8  ! number of fields sent from coupler

      ! buffer location indexes
      integer (kind=4), parameter ::
     &    kbufr_tsurf   =  1 ! surface temperature (K)
     &,   kbufr_cu      =  2 ! bulk coeff for wind stress
     &,   kbufr_ct      =  3 ! bulk coeff for heat flux
     &,   kbufr_ce      =  4 ! bulk coeff for heat flux
     &,   kbufr_albed   =  5 ! surface albedo
     &,   kbufr_dicec   =  6 ! surface albedo
     &,   kbufr_co2     =  7 ! CO2 flux
     &,   kbufr_index   =  8 ! pointer to grid nodes
        !  atm states
      integer (kind=4), parameter ::
     &    kbufs_zlvl    = 1  !  level at which velocity is measured (m)
     &,   kbufs_ztlvl   = 2  !  level at which temperature and humidity are measured (m)
     &,   kbufs_uatm    = 3  !  U comp of wind velocity (m/s)
     &,   kbufs_vatm    = 4  !  V comp of wind velocity (m/s)
     &,   kbufs_potT    = 5  !  potential temperature (K)
     &,   kbufs_Tair    = 6  !  air temperature (K)
     &,   kbufs_Qa      = 7  !  specific humidity (kg/kg)
     &,   kbufs_rhoa    = 8  !  air density (kg/m3)
     &,   kbufs_Fsw     = 9  !  downward shortwave radiation (W/m2): visual-radiative,
     &,   kbufs_Flw     = 10 !  downward longwave radiation (W/m2)
     &,   kbufs_Frain   = 11 !  rainfall rate (kg/(s*m2))
     &,   kbufs_Fsnow   = 12 !  snowfall rate (kg/(s*m2))
     &,   kbufs_cld     = 13 !  total cloudiness (unit)
     &,   kbufs_small_riv = 14  !  river runoffs to be regarded as rainfall (kg/(s*m2))
     &,   kbufs_co2     = 15 !  CO2 concentration
     &,   kbufs_index   = 16 ! pointer to grid nodes

c=======================================================================

      end module atm_indexes

c=======================================================================
