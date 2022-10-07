!     $Id: tides_data.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $

      module tides_data
      Real*8, allocatable, dimension(:,:) ::
     &   etide,utide,vtide,xtide,ytide

      Real*8 td_next

      end module tides_data

