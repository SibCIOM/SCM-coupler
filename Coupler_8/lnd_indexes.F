!     $Id: lnd_indexes.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $

      module lnd_indexes

      implicit none

      integer (kind=4), parameter ::
     &   nsnd         =  4  ! number of fields sent to coupler
     &,  nrcv         =  2  ! number of fields sent from coupler

      ! buffer location indexes
      integer (kind=4), parameter ::
     &    kbufr_Work    =  1 ! reserved space for future
     &,   kbufr_index   =  2 ! pointer to grid nodes
        !  lnd states
      integer (kind=4), parameter ::
     &    kbufs_NatCO2   = 1  !  Natural CO2 emission kg/m2/s
     &,   kbufs_AntCO2   = 2  !  Antropogenic CO2 emission kg/m2/s
     &,   kbufs_FosCO2   = 3  !  Fossil fuel CO2 emission kg/m2/s
     &,   kbufs_index    = 4  !  pointer to grid nodes

c=======================================================================

      end module lnd_indexes

c=======================================================================
