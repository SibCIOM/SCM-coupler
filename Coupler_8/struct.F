!     $Id: struct.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $

	module struct
      Type state
	   integer*4 idate, sec, imt, jmt, nadv
         real*8 next
	End type
      Type grid
	   real*8, pointer, dimension(:,:) :: lon, lat, tarea
	   real*8, pointer, dimension(:,:) :: work
         integer*4, pointer :: mask(:,:)
	End type
      Type from_atm
	   real*8, pointer, dimension(:,:) ::
     &  zlvl, ztlvl, uatm, vatm, potT, Tair, Qa, rhoa, Fsw,
     &  Flw, Frain, Fsnow, cld, small_riv, cco2
	End type
      Type to_atm
        real*8, pointer, dimension(:,:) ::
     &  tsurf, cu, ct, ce, albed, dicec, Fco2
      End type

      Type from_ocn
	   real*8, pointer, dimension(:,:) ::
     &  sst, sss, uocn, vocn, tiltx, tilty, frzmlt, cco2
	End type
      Type to_ocn
	   real*8, pointer, dimension(:,:) ::
     &  Fswthru, Fhnet, Fprec, Fevap, Fmelt, tauxo, tauyo, dtauo,
     &  Fswthru_ice, Fhnet_ice, tauxoi, tauyoi, dtauoi, Frac,
     &  Fwe, dFhnet, utidal, vtidal, Frivprec, Fsalt, Fco2
	End type

      Type from_ice
	   real*8, pointer, dimension(:,:) ::
     &  ailohi, Tsrf, alvdr, alidr, alvdf, alidf, Flatent, Fsensible,
     &  Flwout, evap, tauxa, tauya, Fswthru, Fhnet, Fresh, Fsalt, tauxo,
     &  tauyo, Tref, sabs, dtauo, Fwe
	End type
      Type to_ice
	   real*8, pointer, dimension(:,:) ::
     *   zlvl,ztlvl,uatm,vatm,potT,Tair,Qa,rhoa,swvdr,swvdf,swidr,swidf,
     *   Flw,Frain,Fsnow,sst,sss,uocn,vocn,strtltx,strtlty,frzmlt
	End type

      Type from_lnd
	   real*8, pointer, dimension(:,:) ::
     &  NatCO2,AntCO2,FosCO2
	End type
      Type to_lnd
        real*8, pointer, dimension(:,:) ::
     &  Work
      End type

	type AO_flux
         real*8, pointer, dimension(:,:) ::
     &   taux,  ! wind stress, zonal
     &   tauy,  ! wind stress, meridional
     &   lat,   ! latent          heat flux
     &   sen,   ! sensible        heat flux
     &   dsen,  ! sensible        heat flux increament
     &   lwup,  ! upward longwave heat flux
     &   dlwup, ! upward longwave heat flux increament
     &   evap,  ! evaporation    water flux
     &   fco2,  ! CO2 flux
     &   tref   ! 2m reference temperature
	End type
	Type AO
         real*8, pointer, dimension(:,:) ::
     &   z,       ! bottom atm level height for velocity
     &   zt,      ! bottom atm level height for temperature and humidity
     &   u,       ! bottom atm level wind, zonal
     &   v,       ! bottom atm level wind, meridional
     &   ptem,    ! bottom atm level potential temperature
     &   shum,    ! bottom atm level specific humidity
     &   dens,    ! bottom atm level air density
     &   tbot,    ! bottom atm level temperature
     &   lwdn,    ! downward longwave radiation
     &   swdn,    ! downward shortwave radiation
     &   prec,    ! precipitation rate
     &   aco2,    ! CO2 concentration in atm
     &   sst,     ! sea surface temperature
     &   oco2     ! CO2 concentration in ocn
	End type
      Integer, parameter ::
     &   i2o=1,
     &   i2a=2,
     &   o2i=3,
     &   a2i=4,
     &   o2a=5,
     &   a2o=6,
     &   i2ov=7,
     &   i2av=8,
     &   o2iv=9,
     &   a2iv=10,
     &   o2av=11,
     &   a2ov=12,
     &   a2l=13,
     &   l2a=14,
     &   o2l=15,
     &   l2o=16,
     &   i2l=17,
     &   l2i=18,
     &   a2lv=19,
     &   l2av=20,
     &   o2lv=21,
     &   l2ov=22,
     &   i2lv=23,
     &   l2iv=24

      Real*8, parameter ::
     &   vdr = 0.28,
     &   vdf = 0.24,
     &   idr = 0.31,
     &   idf = 0.17

      end module struct
