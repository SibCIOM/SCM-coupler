!     $Id: prepare.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $

      Subroutine prepare_ice_support(   ToI,
     &               Ice, Igr, Ist,
     &               Ocn, Ogr, Ost,
     &               Lnd, Lgr, Lst,
     &               Atm, Agr, Ast   )

      use struct
      use conversion
      use timer
      use tides_data

      Type (state) Ist,Ost,Lst,Ast
      Type (grid) Igr,Ogr,Lgr,Agr
      Type (from_ice) Ice
      Type (from_ocn) Ocn
      Type (from_lnd) Lnd
      Type (from_atm) Atm
      Type (to_ice) ToI
      Integer i,j
      Real*8, allocatable, dimension(:,:) ::   cld
      Real*8, parameter ::
     &   stebol = 567.0e-10          !  W/m^2/K^4
     &,  icemiss = 0.95              ! emissivity of snow and ice
     &,  ocemiss = 0.985             ! emissivity of snow and ice
     &,  radm = 111000.0  ! *180.0/3.14159265358979 ! m/radian
	logical use_ciaf,tidal,use_frac,use_frcp,use_xie,use_xie4ncep,
     &     use_datahum, atm_shift, old_tides,aomip_rules,prognostic_run,
     &     use_tref,use_small_rivers,var_albo,use_coare,only_land_co2
        integer ciaf_vers,cld_lw_vers,cld_sw_vers,atm_nadv,td_nadv,
     &   first_riv_year,last_riv_year
        real*8 solar_avg,ice_lw_alb
        character*80 file_small_riv
	common/pars/cld_power,use_ciaf,ciaf_vers,tidal,
     &  use_frac,use_frcp,use_tref,
     &  cld_lw_vers,cld_sw_vers,use_xie,use_xie4ncep,
     &  use_datahum,old_tides,aomip_rules,prognostic_run,solar_avg,
     &  use_small_rivers,file_small_riv,first_riv_year,last_riv_year,
     &  use_coare,only_land_co2
        common/rad/var_albo,ice_lw_alb
      logical firsttime
      data firsttime/.true./

      if(.not.allocated(cld)) allocate(cld(Ist%imt,Ist%jmt))

      if(firsttime)then
         firsttime=.false.
c         Open(898,file='rad.txt',form='formatted')
      endif

      if(tidal)then
         nnn=oldtimer
         call clock_starttimer(9)
c         allocate(etide(Ist%imt,Ist%jmt),utide(Ist%imt,Ist%jmt),
c     *      vtide(Ist%imt,Ist%jmt),xtide(Ist%imt,Ist%jmt),
c     *      ytide(Ist%imt,Ist%jmt))
         if(old_tides)
     *      call tides(Ist%idate,Ist%sec,Ist%imt,Ist%jmt,
     *      Igr%mask,Igr%lon,Igr%lat,etide,utide,vtide)
         do i=1,Ist%imt
            do j=1,Ist%jmt
               if(Igr%mask(i,j).gt.0)then
                  i0=max(i-1,1)
                  if(Igr%mask(i0,j).le.0)i0=i
                  i1=min(i+1,Ist%imt)
                  if(Igr%mask(i1,j).le.0)i1=i
                  j0=max(j-1,1)
                  if(Igr%mask(i,j0).le.0)j0=j
                  j1=min(j+1,Ist%jmt)
                  if(Igr%mask(i,j1).le.0)j1=j
                  dx=radm*sqrt(
     *            ((Igr%lon(i1,j)-Igr%lon(i0,j))*cos(Igr%lat(i,j)))**2+
     *            (Igr%lat(i1,j)-Igr%lat(i0,j))**2)
                  dy=radm*sqrt(
     *            ((Igr%lon(i,j1)-Igr%lon(i,j0))*cos(Igr%lat(i,j)))**2+
     *            (Igr%lat(i,j1)-Igr%lat(i,j0))**2)
                  xtide(i,j)=(etide(i1,j)-etide(i0,j))/dx
                  ytide(i,j)=(etide(i,j1)-etide(i,j0))/dy
               endif
            enddo
         enddo
         call clock_starttimer(nnn)
      endif

c      write(nu_diag,*)'a2i'
      call convert_scalar(a2i, Atm%zlvl,   Agr, Ast, ToI%zlvl,    Igr, Ist)
c      write(nu_diag,*)'a2iv'
      call convert_vector(a2iv,Atm%uatm,Atm%vatm, Agr, Ast, ToI%uatm, ToI%vatm, Igr, Ist)
      call convert_scalar(a2i, Atm%potT,   Agr, Ast, ToI%potT,    Igr, Ist)
      call convert_scalar(a2i, Atm%Tair,   Agr, Ast, ToI%Tair,    Igr, Ist)
      call convert_scalar(a2i, Atm%Qa,     Agr, Ast, ToI%Qa,      Igr, Ist)
      call convert_scalar(a2i, Atm%rhoa,   Agr, Ast, ToI%rhoa,    Igr, Ist)
      call convert_scalar(a2i, Atm%Fsw,    Agr, Ast, ToI%swvdr,   Igr, Ist)
      call convert_scalar(a2i, Atm%Frain,  Agr, Ast, ToI%Frain,   Igr, Ist)
      call convert_scalar(a2i, Atm%Fsnow,  Agr, Ast, ToI%Fsnow,   Igr, Ist)
      call convert_scalar(a2i, Atm%Flw,    Agr, Ast, ToI%Flw,     Igr, Ist)
c      write(nu_diag,*)'o2i'
      call convert_scalar(o2i, Ocn%sst,    Ogr, Ost, ToI%sst,     Igr, Ist)

c   Longwave downward radiation calc

c      use_frac=.true.
c      use_frcp=.true.
c      use_max=0.0
      call convert_scalar(a2i, Atm%cld,    Agr, Ast, cld,         Igr, Ist)
c      TatMin=1000.0
c      TatMax=-1000.0
      TrefMax=maxval(Ice%Tref)
      Do i=1,Ist%imt
         Do j=1,Ist%jmt
            frcp=1.0
            frcp0=1.0
            if(use_frac)frcp=Ice%ailohi(i,j)
            if(use_frcp)frcp0=Ice%ailohi(i,j)
            qa=ToI%Qa(i,j)
            if(qa.lt.0.0)qa=0.0
            ToI%Qa(i,j)=qa
            ea=1000.0*qa/(0.622+0.378*qa)
c            ToI%Flw(i,j)=ToI%Flw(i,j)
            if(aomip_rules)then
                if(.false.)then
c            if(use_max.lt.frcp)use_max=frcp
                   sstk=Ice%Tsrf(i,j)*frcp+ToI%sst(i,j)*(1.0-frcp)
                   sstk4=Ice%Tsrf(i,j)**4*frcp+ToI%sst(i,j)**4*(1.0-frcp)
                   sstk40=Ice%Tsrf(i,j)**4
                   sstk3=Ice%Tsrf(i,j)**3*frcp+ToI%sst(i,j)**3*(1.0-frcp)
                   sstk3i=Ice%Tsrf(i,j)**3*frcp*(Ice%Tsrf(i,j)-ToI%Tair(i,j))
                   sstk3o=ToI%sst(i,j)**3*(1.0-frcp)*
     *                (ToI%sst(i,j)-ToI%Tair(i,j))
                   ssti=Ice%Tsrf(i,j)
                   emiss=icemiss*frcp+ocemiss*(1.0-frcp)

c            ToI%Flw(i,j)=stebol*(
c     *            sstk40-emiss/icemiss*sstk4* 
c     *            (0.39-0.05*sqrt(ea))*(1.0-0.8*cld(i,j)) -
c     *            4.0*sstk3*(sstk-ToI%Tair(i,j))
c     *         )*frcp0
                   if(cld_lw_vers.eq.1)then
                      ToI%Flw(i,j)=stebol*(
     *                   sstk40-emiss/icemiss*sstk4* 
     *                   (0.39-0.05*sqrt(ea))*(1.0-0.8*cld(i,j)) -
     *                   4.0*(sstk3i+sstk3o)
     *                )
                   elseif(cld_lw_vers.eq.2)then
                      ToI%Flw(i,j)=stebol*(
     *                   sstk40-emiss/icemiss*sstk4* 
     *                   (0.39-0.05*sqrt(ea))*(1.0-0.6*cld(i,j)**2) -
     *                   4.0*(sstk3i+sstk3o)
     *                )
                   endif
                else
                   Tat=ToI%Tair(i,j)
                   if(use_tref.and.TrefMax.ne.0.0)Tat=Ice%Tref(i,j)
c                if(Ice%ailohi(i,j).gt.0.0.and.Igr%mask(i,j).ne.0.and.
c     *             TatMin.gt.Tat)TatMin=Tat
c                if(Ice%ailohi(i,j).gt.0.0.and.Igr%mask(i,j).ne.0.and.
c     *             TatMax.lt.Tat)TatMax=Tat
                   sstk3i=Ice%Tsrf(i,j)**3*(Ice%Tsrf(i,j)-Tat)
                   sstk3o=ToI%sst(i,j)**3*
     *                (ToI%sst(i,j)-Tat)
                   sstk4i=Ice%Tsrf(i,j)**4
                   sstk4o=ToI%sst(i,j)**4
                   if(cld_lw_vers.eq.1)then
                      Flw_o=ocemiss*stebol*(
     *                   sstk4o*(1.0-(0.39-0.05*sqrt(ea))*(1.0-0.8*Cld(i,j)))-
     *                   4.0*sstk3i
     *                )
                      Flw_i=icemiss*stebol*(
     *                   sstk4i*(1.0-(0.39-0.05*sqrt(ea))*(1.0-0.8*Cld(i,j)))-
     *                   4.0*sstk3i
     *                )
                   elseif(cld_lw_vers.eq.2)then
                      Flw_o=ocemiss*stebol*(
     *                   sstk4o*(1.0-(0.39-0.05*sqrt(ea))*(1.0-0.6*Cld(i,j)**2))-
     *                   4.0*sstk3i
     *                )
                      Flw_i=icemiss*stebol*(
     *                   sstk4i*(1.0-(0.39-0.05*sqrt(ea))*(1.0-0.6*Cld(i,j)**2))-
     *                   4.0*sstk3i
     *                )
                   endif
                   ToI%Flw(i,j)=(frcp*Flw_i+(1.0-frcp)*Flw_o)/icemiss
                endif
            endif
            ToI%Flw(i,j)=ice_lw_alb*ToI%Flw(i,j)*frcp0
            ToI%swvdf(i,j)=ToI%swvdr(i,j)*frcp0*vdf
            ToI%swidr(i,j)=ToI%swvdr(i,j)*frcp0*idr
            ToI%swidf(i,j)=ToI%swvdr(i,j)*frcp0*idf
            ToI%swvdr(i,j)=ToI%swvdr(i,j)*frcp0*vdr
c            if(int(i/20)*20.eq.i.and.int(j/20)*20.eq.j.and.Igr%mask(i,j).ne.0)then
c               Write(898,'(2i4,2f8.2,f7.3,3f7.2,2f7.3)')
c     *            i,j,ToI%swvdr(i,j)+ToI%swvdf(i,j)+ToI%swidr(i,j)+ToI%swidf(i,j),
c     *            ToI%Flw(i,j),frcp,ToI%Tair(i,j),Ice%Tsrf(i,j),ToI%sst(i,j),Cld(i,j),sqrt(ea)
c            endif
       Enddo
      Enddo
c      Print *,'Prepare ice: Tatmin=',TatMin,', TatMax=',TatMax
c      Print *, 'IIIICCCCEEEE compactness',use_max
      call convert_scalar(o2i, Ocn%sss,    Ogr, Ost, ToI%sss,     Igr, Ist)
c      write(nu_diag,*)'o2iv'
      call convert_vector(o2iv, Ocn%uocn, Ocn%vocn, Ogr, Ost, ToI%uocn, ToI%vocn, Igr, Ist)
c      call convert(o2iv, Ocn%uocn,   Ogr, Ost, ToI%uocn,    Igr, Ist)
c      call convert(o2iv, Ocn%vocn,   Ogr, Ost, ToI%vocn,    Igr, Ist)
      call convert_vector(o2iv, Ocn%tiltx, Ocn%tilty, Ogr, Ost, 
     *   ToI%strtltx, ToI%strtlty, Igr, Ist)
      call convert_scalar(o2i, Ocn%frzmlt, Ogr, Ost, ToI%frzmlt,  Igr, Ist)
      if(tidal)then
         ToI%strtltx=ToI%strtltx+xtide
         ToI%strtlty=ToI%strtlty+ytide
         ToI%uocn=ToI%uocn+utide
         ToI%vocn=ToI%vocn+vtide
      endif

      Return
      End

      Subroutine prepare_ocean_support(   ToO,
     &               Ice, Igr, Ist,
     &               Ocn, Ogr, Ost,
     &               Lnd, Lgr, Lst,
     &               Atm, Agr, Ast   )

      use struct
      use conversion
      use tides_data
      use timer

      Type (state) Ist,Ost,Lst,Ast
      Type (grid) Igr,Ogr,Lgr,Agr
      Type (from_ice) Ice
      Type (from_ocn) Ocn
      Type (from_lnd) Lnd
      Type (from_atm) Atm
      Type (to_ocn) ToO
      Type (AO_flux) F
      Type (AO) A
      Real*8, allocatable, target, dimension(:,:) :: frc, sw, swo, swi, swco2
      Real*8, pointer, dimension(:,:) :: frcp, swp, swop, swip, swco2p
      Real*8, allocatable, target, dimension(:,:) ::
     &   taux,  ! wind stress, zonal
     &   tauy,  ! wind stress, meridional
     &   lat,   ! latent          heat flux
     &   sen,   ! sensible        heat flux
     &   lwup,  ! upward longwave heat flux
     &   dsen,  ! sensible        heat flux
     &   dlwup, ! upward longwave heat flux
     &   evap,  ! evaporation    water flux
     &   tref,   ! 2m reference temperature
     &   lwdn,   ! downward longwave heat flux
     &   Fco2,   ! CO2 flux
     &   Qa,
     &   Ta,
     &   Cld,
     &   Tsrf_c
      Real*8, allocatable, target, dimension(:,:) ::
     &   z,       ! bottom atm level height
     &   zt,      ! bottom atm level height for temp
     &   u,       ! bottom atm level wind, zonal
     &   v,       ! bottom atm level wind, meridional
     &   ptem,    ! bottom atm level potential temperature
     &   shum,    ! bottom atm level specific humidity
     &   dens,    ! bottom atm level air density
     &   tbot,    ! bottom atm level temperature
     &   prec,    ! precipitation rate
     &   lwbot,   ! dnwrd longwave radiation
     &   swbot,   ! dnwrd shortwave radiation
     &   a_co2,   ! CO2 concentration in atm
     &   sst,     ! surface temperature
     &   o_co2    ! CO2 concentration in ocn
      Real*8 albo ! ocn albedo
      Real*8, parameter ::
     &   stebol = 567.0e-10          !  W/m^2/K^4
     &,  icemiss = 0.95              ! emissivity of snow and ice
     &,  ocemiss = 0.985             ! emissivity of snow and ice
      logical, save :: initflag
      integer, save :: dddate,cuntr,idate0
      data initflag/.true./,dddate/0/,idate0/0/
      integer nu_diag
      common/iodiag/nu_diag
      character*80 datapath,respath,rpath
      character*1 writemode
      logical use_mask,opennc,var_albo
      integer, save :: nrec
      common/progn_run/datapath,respath,use_mask,writemode
	logical use_ciaf,tidal,use_frac,use_frcp,use_xie,use_xie4ncep,
     &     use_datahum, atm_shift, old_tides,aomip_rules,prognostic_run,
     &     use_tref,use_small_rivers,use_coare,only_land_co2
        integer ciaf_vers,cld_lw_vers,cld_sw_vers,atm_nadv,td_nadv,
     &   first_riv_year,last_riv_year
        real*8 solar_avg, ice_lw_alb
        character*80 file_small_riv
	common/pars/cld_power,use_ciaf,ciaf_vers,tidal,
     &  use_frac,use_frcp,use_tref,
     &  cld_lw_vers,cld_sw_vers,use_xie,use_xie4ncep,
     &  use_datahum,old_tides,aomip_rules,prognostic_run,solar_avg,
     &  use_small_rivers,file_small_riv,first_riv_year,last_riv_year,
     &  use_coare,only_land_co2
        integer, save :: monlen(12)
        common/rad/var_albo,ice_lw_alb
	data monlen/0,31,59,90,120,151,181,212,243,273,304,334/
      
      albo=0.069
c      var_albo=.true.

      imt=Ost%imt
      jmt=Ost%jmt
c      if(initflag)then
      if(.not.allocated(frc))
     &   allocate(frc(imt,jmt),lat(imt,jmt),sen(imt,jmt),lwup(imt,jmt),
     &   evap(imt,jmt),taux(imt,jmt),tauy(imt,jmt),tref(imt,jmt),
     &   Qa(imt,jmt),Ta(imt,jmt),Cld(imt,jmt),swo(imt,jmt),
     &   sw(Ast%imt,Ast%jmt),swi(Ist%imt,Ist%jmt),Tsrf_c(imt,jmt),
     &   dsen(imt,jmt),dlwup(imt,jmt),lwdn(imt,jmt),Fco2(imt,jmt))
      if(.not.allocated(swco2).and..not.only_land_co2)
     &   allocate(swco2(imt,jmt))
c      endif
      frcp=>frc
      F%taux=>taux
      F%tauy=>tauy
      F%lat=>lat
      F%sen=>sen
      F%lwup=>lwup
      F%dsen=>dsen
      F%dlwup=>dlwup
      F%evap=>evap
      F%tref=>tref
      F%fco2=>Fco2
      swp=>sw
      swop=>swo
      swip=>swi
      if(.not.allocated(z)) allocate(z(imt,jmt),u(imt,jmt),v(imt,jmt),ptem(imt,jmt),
     &  shum(imt,jmt),dens(imt,jmt),tbot(imt,jmt),sst(imt,jmt),prec(imt,jmt),
     &  lwbot(imt,jmt),swbot(imt,jmt),zt(imt,jmt),a_co2(imt,jmt),o_co2(imt,jmt))

      A%z=>z
      A%zt=>zt
      A%u=>u
      A%v=>v
      A%ptem=>ptem
      A%shum=>shum
      A%dens=>dens
      A%tbot=>tbot
      A%prec=>prec
      A%lwdn=>lwbot
      A%swdn=>swbot
      A%aco2=>a_co2
      A%sst=>sst
      A%oco2=>o_co2



c      write(nu_diag,*)'a2o'
      call convert_scalar(a2o, Atm%zlvl, Agr, Ast, A%z,    Ogr, Ost)
      call convert_scalar(a2o, Atm%ztlvl, Agr, Ast, A%zt,    Ogr, Ost)
c      write(nu_diag,*)'a2o'
      call convert_vector(a2ov, Atm%uatm, Atm%vatm, Agr, Ast, A%u, A%v, Ogr, Ost)
c      call convert(a2ov, Atm%vatm, Agr, Ast, A%v,    Ogr, Ost)
      call convert_scalar(a2o, Atm%potT, Agr, Ast, A%ptem, Ogr, Ost)
      call convert_scalar(a2o, Atm%Tair, Agr, Ast, A%tbot, Ogr, Ost)
      call convert_scalar(a2o, Atm%Qa,   Agr, Ast, A%shum, Ogr, Ost)
      call convert_scalar(a2o, Atm%rhoa, Agr, Ast, A%dens, Ogr, Ost)
      call convert_scalar(a2o, Atm%Frain, Agr, Ast, A%prec, Ogr, Ost)
      call convert_scalar(a2o, Atm%Flw, Agr, Ast, A%lwdn, Ogr, Ost)
      call convert_scalar(a2o, Atm%Fsw, Agr, Ast, A%swdn, Ogr, Ost)
      call convert_scalar(a2o, Atm%cco2, Agr, Ast, A%aco2, Ogr, Ost)
      sst=Ocn%sst
      o_co2=Ocn%cco2

      call clock_print_time_table('calculate fluxes ')
      call flux_ao(imt, jmt, A, F, Ogr)

      do i=1,imt
      do j=1,jmt
      if (isnan(sen(i,j)))then
      write(nu_diag,*)'i,j = ',i,j
      write(nu_diag,*)'z:',z(i,j)
      write(nu_diag,*)'zt:',zt(i,j)
      write(nu_diag,*)'u:',u(i,j)
      write(nu_diag,*)'v:',v(i,j)
      write(nu_diag,*)'ptem:',ptem(i,j)
      write(nu_diag,*)'shum:',shum(i,j)
      write(nu_diag,*)'rho:',dens(i,j)
      write(nu_diag,*)'tbot:',tbot(i,j)
      write(nu_diag,*)'prec:',prec(i,j)
      write(nu_diag,*)'lwbot:',lwbot(i,j)
      write(nu_diag,*)'swbot:',swbot(i,j)
      write(nu_diag,*)'sst:',sst(i,j)
      write(nu_diag,*)'taux:',taux(i,j)
      write(nu_diag,*)'tauy:',tauy(i,j)
      write(nu_diag,*)'lat:',lat(i,j)
      write(nu_diag,*)'sen:',sen(i,j)
      write(nu_diag,*)'lwup:',lwup(i,j)
      write(nu_diag,*)'dsen:',dsen(i,j)
      write(nu_diag,*)'dlwup:',dlwup(i,j)
      write(nu_diag,*)'evap:',evap(i,j)
      write(nu_diag,*)'tref:',tref(i,j)
      stop
      endif
      enddo
      enddo
      
      call clock_print_time_table('end fluxes       ')
      
c      call correct_wind(Ost%idate, imt, jmt,
c     *   Ogr%lon, Ogr%lat, Ogr%mask, F%taux, F%tauy)

c      call recorrect_wind(Ost%idate, imt, jmt,
c     *   Ogr%lon, Ogr%lat, Ogr%mask, F%taux, F%tauy)

c      write(nu_diag,*)'i2o'
      call convert_scalar(i2o, Ice%ailohi, Igr, Ist, frcp, Ogr, Ost)

c      sw=sw*(1.0-albo)
      ToO%Fswthru=A%swdn
      call convert_scalar(i2o, Ice%Fswthru, Igr, Ist, ToO%Fswthru_ice, Ogr, Ost)

      if(var_albo)then
         year=int(Ost%idate/10000)
         month=int((Ost%idate-year*10000)/100)
         day=Ost%idate-year*10000-month*100
         crnt=(365*year+monlen(month)+day-1)+Ost%sec/86400.0
         n=min(int(crnt-int(crnt/365.0)*365.0)+1,365)
         ToD=360.0*(crnt-int(crnt))
         pp=atan2(1.0,1.0)*4.0/180.0
         DoY=dble(n-1)*360.0/365.0
         delta=23.44*cos((172-DoY)*pp)
         Do i=1,imt
            Do j=1,jmt
               HoD=Ogr%lon(i,j)-180.0+ToD
               cosZ=sin(delta*pp)*sin(Ogr%lat(i,j)*pp)+
     *              cos(delta*pp)*cos(Ogr%lat(i,j)*pp)*cos(HoD*pp)
               if(cosZ.lt.0.0) cosZ=0.0
c               cosZ=90.0-acos(cosZ)/pp
c               albo=0.04+0.211*exp(-cosZ/10.0)
               albo=AMIN1(0.05/(cosZ+0.15),0.15)
               ToO%Fswthru(i,j)=(1.0-albo)*ToO%Fswthru(i,j)
            Enddo
         Enddo
      else
         Do i=1,imt
            Do j=1,jmt
               ToO%Fswthru(i,j)=(1.0-albo)*ToO%Fswthru(i,j)
            Enddo
         Enddo
      endif

c   Longwave downward radiation calc

      call convert_scalar(a2o, Atm%cld,  Agr, Ast, Cld,   Ogr, Ost)
      call convert_scalar(a2o, Atm%Tair,  Agr, Ast, Ta,   Ogr, Ost)
      call convert_scalar(a2o, Atm%Qa,  Agr, Ast, Qa,   Ogr, Ost)
      call convert_scalar(a2o, Atm%Flw,  Agr, Ast, lwup,   Ogr, Ost)
c      Print *,'Ta'
c      Print *,(Ta(i,60),i=165,175)
c      Print *,'To'
c      Print *,(Ocn%sst(i,60),i=165,175)
c      Print *,'Qa'
c      Print *,(Qa(i,60),i=165,175)
c      Print *,'Mask'
c      Print *,(Ogr%mask(i,60),i=165,175)
      lwdn=lwup
      dlwup = 0.0
      call convert_scalar(i2o, Ice%Tsrf, Igr, Ist, Tsrf_c, Ogr, Ost)
c      TatMin=1000.0
c      TatMax=-1000.0
      if(aomip_rules)then
         Do i=1,Ost%imt
            Do j=1,Ost%jmt
               qe=Qa(i,j)
               if(qe.lt.0.0)qe=0.0
               ea=1000.0*qe/(0.622+0.378*qe)
               frcp0=1.0
               if(use_frac)frcp0=1.0-frcp(i,j)
               if(.false.)then
                  emiss=icemiss*(1.0-frcp0)+ocemiss*frcp0
                  sstk=Tsrf_c(i,j)*(1.0-frcp0)+
     *               Ocn%sst(i,j)*frcp0
                  sstk4=Tsrf_c(i,j)**4*(1.0-frcp0)+
     *               Ocn%sst(i,j)**4*frcp0
                  sstk40=Ocn%sst(i,j)**4
                  sstk3=Tsrf_c(i,j)**3*(1.0-frcp0)+
     *               Ocn%sst(i,j)**3*frcp0
                  sstk3i=Tsrf_c(i,j)**3*(1.0-frcp0)*(Tsrf_c(i,j)-Ta(i,j))
                  sstk3o=Ocn%sst(i,j)**3*frcp0*(Ocn%sst(i,j)-Ta(i,j))

                  if(cld_lw_vers.eq.1)then
                     lwup(i,j)=-emiss*stebol*(
     *                  sstk4*(0.39-0.05*sqrt(ea))*(1.0-0.8*Cld(i,j))+
     *                  4.0*(sstk3i+sstk3o)
     *               )
                  elseif(cld_lw_vers.eq.2)then
                     lwup(i,j)=-emiss*stebol*(
     *                  sstk4*(0.39-0.05*sqrt(ea))*(1.0-0.6*Cld(i,j)**2)+
     *                  4.0*(sstk3i+sstk3o)
     *               )
                  endif
c            lwup(i,j)=-emiss*stebol*(
c     *         sstk4*(0.39-0.05*sqrt(ea))*(1.0-0.8*Cld(i,j))+
c     *         4.0*sstk3*(sstk-Ta(i,j))
c     *      )
               else
                  Tat=Ta(i,j)
                  if(use_tref)Tat=tref(i,j)
c                if(Ogr%mask(i,j).ne.0.and.TatMin.gt.Tat)TatMin=Tat
c                if(Ogr%mask(i,j).ne.0.and.TatMax.lt.Tat)TatMax=Tat
                  sstk3i=Tsrf_c(i,j)**3*(Tsrf_c(i,j)-Tat)
                  sstk3o=-Ocn%sst(i,j)**3*Tat
                  dsstk3o=Ocn%sst(i,j)**3
                  sstk4i=Tsrf_c(i,j)**4
                  sstk4o=-3.0*Ocn%sst(i,j)**4
                  dsstk4o=4.0*Ocn%sst(i,j)**3
                  if(cld_lw_vers.eq.1)then
                     Flw_o=ocemiss*stebol*(
     *                  sstk4o*(1.0-(0.39-0.05*sqrt(ea))*(1.0-0.8*Cld(i,j)))-
     *                  4.0*sstk3o
     *               )
                     dFlw_o=ocemiss*stebol*(
     *                  dsstk4o*(1.0-(0.39-0.05*sqrt(ea))*(1.0-0.8*Cld(i,j)))-
     *                  4.0*dsstk3o
     *               )
                     Flw_i=icemiss*stebol*(
     *                  sstk4i*(1.0-(0.39-0.05*sqrt(ea))*(1.0-0.8*Cld(i,j)))-
     *                  4.0*sstk3o
     *               )
                     dFlw_i=icemiss*stebol*(
     *                  -4.0*dsstk3o
     *               )
                  elseif(cld_lw_vers.eq.2)then
                     Flw_o=ocemiss*stebol*(
     *                  sstk4o*(1.0-(0.39-0.05*sqrt(ea))*(1.0-0.6*Cld(i,j)**2))-
     *                  4.0*sstk3o
     *               )
                     dFlw_o=ocemiss*stebol*(
     *                  dsstk4o*(1.0-(0.39-0.05*sqrt(ea))*(1.0-0.6*Cld(i,j)**2))-
     *                  4.0*dsstk3o
     *               )
                     Flw_i=icemiss*stebol*(
     *                  sstk4i*(1.0-(0.39-0.05*sqrt(ea))*(1.0-0.6*Cld(i,j)**2))-
     *                  4.0*sstk3o
     *               )
                     dFlw_i=icemiss*stebol*(
     *                  -4.0*dsstk3o
     *               )
                  endif
                  lwup(i,j)=frcp0*Flw_o+(1.0-frcp0)*Flw_i-ocemiss*stebol*sstk4o
                  dlwup(i,j)=frcp0*dFlw_o+(1.0-frcp0)*dFlw_i-ocemiss*stebol*dsstk4o
               endif
            Enddo
         Enddo
      else
         Do i=1,Ost%imt
            Do j=1,Ost%jmt
               if(Ogr%mask(i,j).ne.0)then
                  lwup(i,j)=lwup(i,j)+3.0*ocemiss*stebol*Ocn%sst(i,j)**4
                  dlwup(i,j)=dlwup(i,j)-4.0*ocemiss*stebol*Ocn%sst(i,j)**3
               endif
            enddo
         enddo
      endif

c      Print *,'LW'
c      Print *,(lwup(i,60),i=165,175)
      Do i=1,Ost%imt
         Do j=1,Ost%jmt
            if(Ogr%mask(i,j).ne.0)then
               lwdn(i,j)=(1.0-frcp(i,j))*(
     *            lwup(i,j)+dlwup(i,j)*Ocn%sst(i,j)
     *            )
            endif
         enddo
      enddo
c      Print *,'Prepare ocn: Tatmin=',TatMin,', TatMax=',TatMax
      
      call convert_scalar(i2o, Ice%Fwe,     Igr, Ist, ToO%Fwe,         Ogr, Ost)
      call convert_scalar(i2o, Ice%Fhnet,   Igr, Ist, ToO%Fhnet_ice,   Ogr, Ost)
      Do i=1,imt
         Do j=1,jmt
            ToO%Fhnet(i,j)=sen(i,j)+lat(i,j)+lwup(i,j)
            ToO%dFhnet(i,j)=dsen(i,j)+dlwup(i,j)
         Enddo
      Enddo

      Do i=1,Ast%imt
         Do j=1,Ast%jmt
            sw(i,j)=Atm%Frain(i,j)+Atm%Fsnow(i,j)
         Enddo
      Enddo
      Do i=1,Ist%imt
         Do j=1,Ist%jmt
c            Dice=Ice%Fsalt(i,j)*Ice%sss(i,j)/4.0
            swi(i,j)=Ice%Fresh(i,j)
c     *           -Dice+Dice*(Ice%sss(i,j)-4.0)/Ice%sss(i,j)
         Enddo
      Enddo
        call convert_scalar(a2o, swp, Agr, Ast, swop, Ogr, Ost)
      if(.not.only_land_co2)call convert_scalar(l2o, Lnd%NatCO2, Lgr, Lst, swco2, Ogr, Ost)
      Do i=1,imt
         Do j=1,jmt
            ToO%Fprec(i,j)=swo(i,j)
            ToO%Fevap(i,j)=F%evap(i,j)
	    if(only_land_co2)then
               ToO%Fco2(i,j)=-F%fco2(i,j)
	    else
	       ToO%Fco2(i,j)=-swco2(i,j)
	    endif
         Enddo
      Enddo
      if(.not.only_land_co2)then
         call convert_scalar(l2o, Lnd%AntCO2, Lgr, Lst, swco2, Ogr, Ost)
         Do i=1,imt
            Do j=1,jmt
	       ToO%Fco2(i,j)=ToO%Fco2(i,j)-swco2(i,j)
            Enddo
         Enddo
!         call convert_scalar(l2o, Lnd%FosCO2, Lgr, Lst, swco2, Ogr, Ost)
!         Do i=1,imt
!            Do j=1,jmt
!	       ToO%Fco2(i,j)=ToO%Fco2(i,j)-swco2(i,j)
!            Enddo
!         Enddo
      endif
      if(only_land_co2)write(nu_diag,*)'Only land co2 flux is considered'

        call convert_scalar(i2o, swip, Igr, Ist, swop, Ogr, Ost)
      Do i=1,imt
         Do j=1,jmt
              ToO%Fmelt(i,j)=swo(i,j)
         Enddo
      Enddo

      ToO%Frivprec=0.0
      if(use_small_rivers)then
         iyear0=int(Ost%idate/10000)
         iyear=min(max(iyear0,first_riv_year),last_riv_year)
            Open(23,file=file_small_riv,form='unformatted',
     *         access='direct',recl=imt*jmt*8)
            kmo=int((Ost%idate-10000*iyear0)/100)
            nrec=(iyear-first_riv_year)*12+kmo
            write(nu_diag,*)'Rivers from file: ',nrec,kmo
            write(nu_diag,*)iyear,iyear0,first_riv_year,Ost%idate
            Read(23,rec=nrec)ToO%Frivprec
            Close(23)
      endif

        call convert_scalar(a2o, Atm%small_riv, Agr, Ast, swop, Ogr, Ost)
      Do i=1,imt
         Do j=1,jmt
            ToO%Frivprec(i,j)=ToO%Frivprec(i,j)+swo(i,j)
         Enddo
      Enddo

c      write(nu_diag,*)'i2ov'
        call convert_vector(i2ov, Ice%tauxo, Ice%tauyo, Igr, Ist, 
     *     ToO%tauxoi, ToO%tauyoi, Ogr, Ost)
        call convert_scalar(i2o, Ice%dtauo, Igr, Ist, ToO%dtauoi, Ogr, Ost)
c        call convert(i2o, frcp,      Igr, Ist, ToO%Frac,   Ogr, Ost)
        if (tidal) then
           call convert_vector(i2ov, utide, vtide, Igr, Ist, 
     *        ToO%utidal, ToO%vtidal, Ogr, Ost)
        else
           ToO%utidal = 0.0
           ToO%vtidal = 0.0
        endif
      Do i=1,imt
         Do j=1,jmt
              ToO%tauxo(i,j)=F%taux(i,j)
              ToO%tauyo(i,j)=F%tauy(i,j)
              ToO%dtauo(i,j)=0.0
c              ToO%tauxoi(i,j)=0.1*F%taux(i,j)
c              ToO%tauyoi(i,j)=0.1*F%tauy(i,j)
c              ToO%dtauoi(i,j)=0.0
              ToO%Frac(i,j)=frcp(i,j)
         Enddo
      Enddo
      

      initflag=.false.

c      call cdf2Dvar('Fwe',ToO%Fwe,Ogr%mask,imt,jmt,.true.,Ost%idate)

      if(Ost%idate.ne.idate0)then
         opennc=idate0.eq.0
         idate0=Ost%idate
         if(writemode.eq.'w')then
            mh=imt
            nh=jmt
            rpath=trim(datapath)//trim(respath)
            Write(nu_diag,*)'Write data to '//trim(rpath)
            if(opennc)then
               call open_w2w(trim(rpath)//'Flw0.nc','Flw',mh,nh)
               opennc=.false.
               nrec=1
            endif
            call write_w(trim(rpath)//'Flw0.nc',nrec,idate0,'Flw',lwdn,Ogr%mask,mh,nh)
            nrec=nrec+1
         endif
      endif

c      deallocate(frc,lat,sen,lwup,evap,taux,tauy,tref,Qa,Ta,Cld,swo,sw,swi,Tsrf_c,
c     &   dsen,dlwup,lwdn)
c      deallocate(z,u,v,ptem,shum,dens,tbot,sst)

      Return
      End

      Subroutine prepare_atmos_support(   ToA,
     &               Ice, Igr, Ist,
     &               Ocn, Ogr, Ost,
     &               Lnd, Lgr, Lst,
     &               Atm, Agr, Ast   )

      use struct
      use conversion
      use tides_data
c      use timer

      Type (state) Ist,Ost,Lst,Ast
      Type (grid) Igr,Ogr,Lgr,Agr
      Type (from_ice) Ice
      Type (from_ocn) Ocn
      Type (from_lnd) Lnd
      Type (from_atm) Atm
      Type (to_atm) ToA
      Type (AO_flux) F
      Type (AO) A

      Real*8, allocatable, dimension(:,:) :: omask

      Real*8, allocatable, target, dimension(:,:) ::
     &   taux,  ! wind stress, zonal
     &   tauy,  ! wind stress, meridional
     &   lat,   ! latent          heat flux
     &   sen,   ! sensible        heat flux
     &   lwup,  ! upward longwave heat flux
     &   dsen,  ! sensible        heat flux
     &   dlwup, ! upward longwave heat flux
     &   evap,  ! evaporation    water flux
     &   tref,  ! 2m reference temperature
     &   Fco2,
     &   ifrc,  ! ice fraction
     &   asen,
     &   aevap,
     &   ataux,
     &   atauy

      Real*8, allocatable, target, dimension(:,:) ::
     &   z,       ! bottom atm level height
     &   zt,      ! bottom atm level height
     &   u,       ! bottom atm level wind, zonal
     &   v,       ! bottom atm level wind, meridional
     &   ptem,    ! bottom atm level potential temperature
     &   shum,    ! bottom atm level specific humidity
     &   dens,    ! bottom atm level air density
     &   tbot,    ! bottom atm level temperature
     &   prec,
     &   swbot,
     &   lwbot,
     &   a_co2,
     &   sst,     ! bottom atm level temperature
     &   o_co2,   !
     &   ln_natco2,   !
     &   ln_antco2,   !
     &   ln_fosco2    !

      Real*8 albo ! ocn albedo
      Real*8, allocatable :: albi(:,:)

      Real*8, parameter ::
     &   stebol = 567.0e-10          !  W/m^2/K^4
     &,  icemiss = 0.95              ! emissivity of snow and ice
     &,  ocemiss = 0.985             ! emissivity of snow and ice

      integer nu_diag
      common/iodiag/nu_diag

      logical use_ciaf,tidal,use_frac,use_frcp,use_xie,use_xie4ncep,
     &   use_datahum, atm_shift, old_tides,aomip_rules,prognostic_run,
     &   use_tref,use_small_rivers,use_coare,only_land_co2
      integer ciaf_vers,cld_lw_vers,cld_sw_vers,atm_nadv,td_nadv,
     &   first_riv_year,last_riv_year
        real*8 solar_avg,ice_lw_alb
        character*80 file_small_riv
	common/pars/cld_power,use_ciaf,ciaf_vers,tidal,
     &  use_frac,use_frcp,use_tref,
     &  cld_lw_vers,cld_sw_vers,use_xie,use_xie4ncep,
     &  use_datahum,old_tides,aomip_rules,prognostic_run,solar_avg,
     &  use_small_rivers,file_small_riv,first_riv_year,last_riv_year,
     &  use_coare,only_land_co2
        common/rad/var_albo,ice_lw_alb
      real*8, parameter ::
     &  qqqice       = 11637800.0       ! for qsat over ice
     &, TTTice       = 5897.8           ! for qsat over ice
        integer, save :: monlen(12)
	data monlen/0,31,59,90,120,151,181,212,243,273,304,334/

      logical var_albo
      
      esat(Tk) = 0.98*exp((0.7859+0.03477*Tk)/
     *   (1.0+0.00412*Tk)*2.302585093)
      qsat(Tk) = 0.622*esat(Tk-273.16)/(1000.0-0.378*esat(Tk-273.16))

c     Prepare  mask, tsurf, sens, lat, evap, taux, tauy, ifrc, swave, lwave

      albo=0.069
c      var_albo = .true.
      ueps=0.0000001
      teps=0.0000001
      qeps=0.000000001

      imt=Ast%imt
      jmt=Ast%jmt

      if(.not.allocated(lat))
     &   allocate(lat(imt,jmt),sen(imt,jmt),lwup(imt,jmt),
     &   evap(imt,jmt),taux(imt,jmt),tauy(imt,jmt),tref(imt,jmt),
     &   dsen(imt,jmt),dlwup(imt,jmt),ifrc(imt,jmt),asen(imt,jmt),aevap(imt,jmt),
     &   ataux(imt,jmt),atauy(imt,jmt),Fco2(imt,jmt))

      F%taux=>taux
      F%tauy=>tauy
      F%lat=>lat
      F%sen=>sen
      F%lwup=>lwup
      F%dsen=>dsen
      F%dlwup=>dlwup
      F%evap=>evap
      F%tref=>tref
      F%fco2=>Fco2

      if(.not.allocated(z))
     &  allocate(z(imt,jmt),zt(imt,jmt),u(imt,jmt),v(imt,jmt),ptem(imt,jmt),
     &  shum(imt,jmt),dens(imt,jmt),tbot(imt,jmt),sst(imt,jmt),prec(imt,jmt),
     &  swbot(imt,jmt),lwbot(imt,jmt),a_co2(imt,jmt),o_co2(imt,jmt),
     &  ln_natco2(imt,jmt),ln_antco2(imt,jmt),ln_fosco2(imt,jmt))

      A%z=>z
      A%zt=>zt
      A%u=>u
      A%v=>v
      A%ptem=>ptem
      A%shum=>shum
      A%dens=>dens
      A%tbot=>tbot
      A%prec=>prec
      A%lwdn=>lwbot
      A%swdn=>swbot
      A%aco2=>a_co2
      A%sst=>sst
      A%oco2=>o_co2

      z=Atm%zlvl
      zt=Atm%ztlvl
      u=Atm%uatm
      v=Atm%vatm
      ptem=Atm%potT
      tbot=Atm%Tair
      shum=Atm%Qa
      dens=Atm%rhoa
      prec=Atm%Frain
      lwbot=Atm%Flw
      a_co2=Atm%cco2
      Do i=1,Ast%imt
         Do j=1,Ast%jmt
            swbot(i,j)=Atm%Fsw(i,j)
         Enddo
      Enddo

      call convert_scalar(i2a, Ice%ailohi, Igr, Ist, ifrc, Agr, Ast)
      call convert_scalar(i2a, Ice%Tsrf, Igr, Ist, ToA%tsurf, Agr, Ast)
      call convert_scalar(o2a, Ocn%sst, Ogr, Ost, sst, Agr, Ast)
      call convert_scalar(o2a, Ocn%cco2, Ogr, Ost, o_co2, Agr, Ast)

      call flux_ao(imt, jmt, A, F, Agr)
      
      allocate(albi(imt,jmt))

      if(var_albo)then
         year=int(Ast%idate/10000)
         month=int((Ast%idate-year*10000)/100)
         day=Ast%idate-year*10000-month*100
         crnt=(365*year+monlen(month)+day-1)+Ast%sec/86400.0
         n=min(int(crnt-int(crnt/365.0)*365.0)+1,365)
         ToD=360.0*(crnt-int(crnt))
         pp=atan2(1.0,1.0)*4.0/180.0
         DoY=dble(n-1)*360.0/365.0
         delta=23.44*cos((172-DoY)*pp)
         Do i=1,imt
            Do j=1,jmt
               HoD=Agr%lon(i,j)-180.0+ToD
               cosZ=sin(delta*pp)*sin(Agr%lat(i,j)*pp)+
     *              cos(delta*pp)*cos(Agr%lat(i,j)*pp)*cos(HoD*pp)
               if(cosZ.lt.0.0) cosZ=0.0
c               cosZ=90.0-acos(cosZ)/pp
c               albo=0.04+0.211*exp(-cosZ/10.0)
               albo=AMIN1(0.05/(cosZ+0.15),0.15)
               ToA%albed(i,j)=albo*(1.-ifrc(i,j))
            Enddo
         Enddo
      else
         Do i=1,imt
            Do j=1,jmt
               ToA%albed(i,j)=albo*(1.-ifrc(i,j))
            Enddo
         Enddo
      endif

      call convert_scalar(i2a, Ice%alvdr,Igr,Ist, albi,Agr,Ast)
      Do i=1,imt
         Do j=1,jmt
            ToA%albed(i,j)=ToA%albed(i,j)+albi(i,j)*ifrc(i,j)*vdr
         Enddo
      Enddo

      call convert_scalar(i2a, Ice%alvdf,Igr,Ist, albi,Agr,Ast)
      Do i=1,imt
         Do j=1,jmt
            ToA%albed(i,j)=ToA%albed(i,j)+albi(i,j)*ifrc(i,j)*vdf
         Enddo
      Enddo
      call convert_scalar(i2a, Ice%alidr,Igr,Ist, albi,Agr,Ast)
      Do i=1,imt
         Do j=1,jmt
            ToA%albed(i,j)=ToA%albed(i,j)+albi(i,j)*ifrc(i,j)*idr
         Enddo
      Enddo
      call convert_scalar(i2a, Ice%alidf,Igr,Ist, albi,Agr,Ast)
      call convert_scalar(l2a, Lnd%NatCO2,Lgr,Lst, ln_natco2,Agr,Ast)
      call convert_scalar(l2a, Lnd%AntCO2,Lgr,Lst, ln_antco2,Agr,Ast)
!      call convert_scalar(l2a, Lnd%FosCO2,Lgr,Lst, ln_fosco2,Agr,Ast)

      r=1.0     ! rate to be changed
      land=2    ! land=2, ocean=1
      Do i=1,imt
         Do j=1,jmt
            ToA%albed(i,j)=ToA%albed(i,j)+albi(i,j)*ifrc(i,j)*idf
            if(Agr%mask(i,j).eq.land)then
!               ToA%fco2(i,j)=ln_natco2(i,j)+ln_fosco2(i,j)+
!     *            r*(ln_antco2(i,j)-ln_fosco2(i,j))
               ToA%fco2(i,j)=-ln_natco2(i,j)-ln_antco2(i,j)
            else
               if(only_land_co2)then
                  ToA%fco2(i,j)=-Fco2(i,j)
               else
                  ToA%fco2(i,j)=-ln_natco2(i,j)-ln_antco2(i,j)
               endif
            endif
         Enddo
      Enddo

      deallocate(albi)

      call convert_scalar(i2a, Ice%Fsensible,Igr, Ist, asen,  Agr, Ast)
      call convert_scalar(i2a, Ice%Flatent,Igr, Ist, aevap,  Agr, Ast)
      call convert_vector(i2av,Ice%tauxa,Ice%tauya,Igr,Ist,ataux,atauy,Agr,Ast)

      Do i=1,imt
         Do j=1,jmt
            frci=ifrc(i,j)
            frco=1.0-frci

c            asen(i,j) = frci*asen(i,j) + frco*(sen(i,j)+dsen(i,j)*sst(i,j))
c            aevap(i,j) = frci*aevap(i,j) + frco*evap(i,j)
            ataux(i,j) = frci*ataux(i,j) + frco*taux(i,j)
            atauy(i,j) = frci*atauy(i,j) + frco*tauy(i,j)
            ustar=max(u(i,j)**2+v(i,j)**2,ueps)
            ToA%cu(i,j)=sqrt((ataux(i,j)**2+atauy(i,j)**2)/ustar)/dens(i,j)

c            qsatur = qqqIce * exp(-TTTIce/ToA%tsurf(i,j))*frci+qsat(sst(i,j))*frco

c            cp=1005.0 + 1870.0*qsatur
c            ToA%ct(i,j)=asen(i,j)/dens(i,j)/cp/(ptem(i,j)-ToA%tsurf(i,j))
c            ToA%ce(i,j)=aevap(i,j)/dens(i,j)/(shum(i,j)-qsatur)

            qsatur = qsat(sst(i,j))
            cp=1005.0 + 1870.0*qsatur
            cto=0.0012
            ceo=0.0015
            if(use_coare)then
               if(dabs(tbot(i,j)-sst(i,j)).gt.teps)
     *            cto=(sen(i,j)+dsen(i,j)*sst(i,j))/dens(i,j)/cp/(tbot(i,j)-sst(i,j))
            else
               if(dabs(ptem(i,j)-sst(i,j)).gt.teps)
     *            cto=(sen(i,j)+dsen(i,j)*sst(i,j))/dens(i,j)/cp/(ptem(i,j)-sst(i,j))
            endif
            if(dabs(shum(i,j)-qsatur).gt.qeps)
     *         ceo=evap(i,j)/dens(i,j)/(shum(i,j)-qsatur)

            cti=0.0012
            cei=0.0015
            if (dabs(ToA%tsurf(i,j)).gt.teps)then
               qsatur = qqqIce * exp(-TTTIce/ToA%tsurf(i,j))
               cp=1005.0 + 1870.0*qsatur
               if(dabs(ptem(i,j)-ToA%tsurf(i,j)).gt.teps)
     *            cti=asen(i,j)/dens(i,j)/cp/(ptem(i,j)-ToA%tsurf(i,j))
               if(dabs(shum(i,j)-qsatur).gt.qeps)
     *            cei=aevap(i,j)/2.5e+6/dens(i,j)/(shum(i,j)-qsatur)
            endif

            ToA%ct(i,j)=frci*cti+frco*cto
c            ToA%ct(i,j)=cto
c            if(cto<0)then
c               Print *,'!!!! ',i,j,sen(i,j),dsen(i,j),sst(i,j)
c               Print *,'!!!! ',use_coare,tbot(i,j),ptem(i,j)
c               Print *,'!!!! ',sen(i,j)+dsen(i,j)*sst(i,j)
c               Print *,'!!!! ',dens(i,j),tbot(i,j)-sst(i,j)
c            endif
            ToA%ce(i,j)=frci*cei+frco*ceo

            ToA%tsurf(i,j)=frci*ToA%tsurf(i,j)+frco*sst(i,j)
            ToA%dicec(i,j)=frci

         Enddo
      Enddo

c      Do i=1,imt
c         Do j=1,jmt
c            if(ifrc(i,j).gt.0.0)then
c               ToA%tsurf(i,j)=ToA%tsurf(i,j)*ifrc(i,j)+
c     +            sst(i,j)*(1.0-ifrc(i,j))
c            else
c               ToA%tsurf(i,j)=sst(i,j)
c            endif
c         Enddo
c      Enddo

c      deallocate(lat,sen,lwup,evap,taux,tauy,tref,dsen,dlwup)
c      deallocate(z,zt,u,v,ptem,shum,dens,tbot,sst)

      Return
      End

      Subroutine prepare_land_support(   ToA,
     &               Ice, Igr, Ist,
     &               Ocn, Ogr, Ost,
     &               Lnd, Lgr, Lst,
     &               Atm, Agr, Ast   )

      use struct

      Type (state) Ist,Ost,Lst,Ast
      Type (grid) Igr,Ogr,Lgr,Agr
      Type (from_ice) Ice
      Type (from_ocn) Ocn
      Type (from_lnd) Lnd
      Type (from_atm) Atm
      Type (to_lnd) ToL

      ToL%Work = 0.0

      Return
      End

