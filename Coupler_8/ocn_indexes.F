!     $Id: ocn_indexes.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $

      module ocn_indexes

      implicit none

      integer (kind=4), parameter ::
     &   nsnd         =  9  ! number of fields sent to coupler
     &,  nrcv         = 21  ! number of fields sent from coupler

      ! buffer location indexes
      integer (kind=4), parameter ::
        !  ocean states
     &   kbufr_sst     =  1 ! sea surface temp         (K)
     &,  kbufr_sss     =  2 ! sea surface salinity     (o/oo)
     &,  kbufr_uocn    =  3 ! ocean current            (m/s)
     &,  kbufr_vocn    =  4 ! 
     &,  kbufr_tiltx   =  5 ! sea surface slope  
     &,  kbufr_tilty   =  6 !
     &,  kbufr_frzmlt  =  7 !
     &,  kbufr_co2     =  8 ! CO2 concentration
     &,  kbufr_index   =  9 !
        !  to send to ocean
     &,  kbufs_swave   =  1 ! penetrative short wave radiation
     &,  kbufs_netheat =  2 ! net heat flux
     &,  kbufs_prec    =  3 ! precipitation flux
     &,  kbufs_evap    =  4 ! evaporation flux
     &,  kbufs_taux    =  5 ! surface stress components
     &,  kbufs_tauy    =  6 ! 
     &,  kbufs_dtau    =  7 ! surface stress components
     &,  kbufs_ifrc    =  8 ! ice fraction
     &,  kbufs_iswave  =  9 ! through ice short wave radiation
     &,  kbufs_inetheat= 10 ! net heat flux under ice cover
     &,  kbufs_itaux   = 11 ! surface stress components under ice
     &,  kbufs_itauy   = 12 ! 
     &,  kbufs_idtau   = 13 ! surface stress components under ice
     &,  kbufs_melt    = 14 ! melting-freezing water flux
     &,  kbufs_dnetheat= 15 ! net heat flux increament 
     &,  kbufs_utidal  = 16 ! u component of tidal velocity m/s
     &,  kbufs_vtidal  = 17 ! v component of tidal velocity m/s
     &,  kbufs_sriv    = 18 ! small river as precipitation flux
     &,  kbufs_co2     = 19 ! CO2 flux
     &,  kbufs_we      = 20 ! water equivalent m
     &,  kbufs_index   = 21 ! a counter
c     &,  kbufs_salt    = 14 ! melting-freezing water flux

c=======================================================================

      end module ocn_indexes

c=======================================================================
