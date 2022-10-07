!     $Id: coupler.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $

      program coupler

      use struct
      use timer
      use tides_data
      use cpl_indexes
      use conversion
      use msgtypes

      integer*4, parameter ::
     &   nbuff = 100

      integer cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs(20)
      common/allids/cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs

      integer ::
     &  cpl_task
     &, my_task           ! master task for ice
     &, MPI_COMM_MY
      common/tasks/my_task, cpl_task, MPI_COMM_MY

	integer*4, dimension(nbuff) ::
     &   buff

      integer*4 ierr,stop_now,ctrl
      common/stop_flag/stop_now

	Type (state) Ast,Ist,Ost,Lst

      integer nu_diag
      common/iodiag/nu_diag

      Type (grid) Agr,Igr,Ogr,Lgr
	Type (from_atm) A
	Type (from_ocn) O
	Type (from_ice) I
	Type (from_lnd) FromL
	Type (to_ocn) ToO
	Type (to_ice) ToI
	Type (to_atm) ToA
	Type (to_lnd) ToL

	!  Atmosphere state
      integer*4 a_idate,a_sec,a_imt,a_jmt,a_nadv
      common/atm_state/a_idate,a_sec,a_imt,a_jmt,a_nadv
	real*8, allocatable, target ::
     & a_lon(:,:), a_lat(:,:), a_work(:,:), a_tarea(:,:)
      integer*4, allocatable, target ::
     & a_mask(:,:)
	!  From atmosphere
	real*8, allocatable, target, dimension(:,:) ::
     &  a_zlvl, a_ztlvl, a_uatm, a_vatm, a_potT, a_Tair, a_Qa, a_rhoa, 
     &  a_Fsw, a_Flw, a_Frain, a_Fsnow, a_cld,
     &  a_small_riv, a_co2
	!  To atmosphere
	real*8, allocatable, target, dimension(:,:) ::
     &  a_tsurf, a_Cu, a_Ct, a_Ce, a_albed, a_dicec, a_fco2

        real*8 atm_mpt

	!  Ice state
      integer*4 i_idate,i_sec,i_imt,i_jmt,i_nadv
      common/ice_state/i_idate,i_sec,i_imt,i_jmt,i_nadv
	real*8, allocatable, target ::
     & i_lon(:,:), i_lat(:,:), i_work(:,:), i_tarea(:,:)
      integer*4, allocatable, target ::
     & i_mask(:,:)
	!  From ice
	real*8, allocatable, target, dimension(:,:) ::
     &  ailohi, Tsrf, alvdr, alidr, alvdf, alidf, Flatent, Fsensible,
     &  Flwout, evap, tauxa, tauya, Fswthru, Fhnet, Fresh, Fsalt, tauxo,
     &  tauyo, Tref, sabs, dtauo, Fwe
	!  To ice
	real*8, allocatable, target, dimension(:,:) ::
     *   zlvl,ztlvl,uatm,vatm,potT,Tair,Qa,rhoa
     *  ,swvdr,swvdf,swidr,swidf,Flw,Frain,Fsnow,sst,sss,uocn,vocn
     *  ,strtltx,strtlty,frzmlt

        real*8 ice_mpt

	!  Ocean state
      integer*4 o_idate,o_sec,o_imt,o_jmt,o_nadv
      common/ocn_state/o_idate,o_sec,o_imt,o_jmt,o_nadv
	real*8, allocatable, target ::
     & o_lon(:,:), o_lat(:,:), o_work(:,:), o_tarea(:,:)
      integer*4, allocatable, target ::
     & o_mask(:,:)
	! From ocean
	real*8, allocatable, target, dimension(:,:) ::
     &  o_sst, o_sss, o_uocn, o_vocn, o_strtltx, o_strtlty, o_frzmlt, o_co2
	! To ocean
	real*8, allocatable, target, dimension(:,:) ::
     &  o_Fswthru, o_Fhnet, o_Fprec, o_Fevap, o_tauxo, o_tauyo, 
     &  o_dtauo, o_Fmelt, o_Fsalt,
     &  o_Fswthru_ice, o_Fhnet_ice, o_tauxoi, o_tauyoi, 
     &  o_utidal, o_vtidal, o_sriv,
     &  o_dtauoi, o_Frac, o_Fwe, o_dFhnet, o_fco2

        real*8 ocn_mpt

	!  Land state
      integer*4 ln_idate,ln_sec,ln_imt,ln_jmt,ln_nadv
      common/lnd_state/ln_idate,ln_sec,ln_imt,ln_jmt,ln_nadv
	real*8, allocatable, target ::
     & ln_lon(:,:), ln_lat(:,:), ln_wrk(:,:), ln_tarea(:,:)
      integer*4, allocatable, target ::
     & ln_mask(:,:)
	!  From land
	real*8, allocatable, target, dimension(:,:) ::
     &  ln_NatCO2, ln_AntCO2, ln_FosCO2
	!  To land
	real*8, allocatable, target, dimension(:,:) ::
     &  ln_Work

        real*8 lnd_mpt

	Logical running
	real*8 next,current,finish
	integer cur_date,cur_sec,old_date
        character*80 wap_file,progn_atm_path,res_atm_path
        character*21 timestring
        character*1 progn_write_mode

	logical use_ciaf,tidal,use_frac,use_frcp,use_xie,use_xie4ncep,
     &     use_datahum, old_tides,aomip_rules,prognostic_run,
     &     progn_use_mask,use_tref,use_rivers,var_ocnalb
        integer ciaf_vers,cld_lw_vers,cld_sw_vers,atm_nadv,td_nadv,nrivers
        real*8 solar_avg, ice_lw_alb
        real*8, allocatable :: rivcoords(:),rtrans(:)
        integer first_riv_year,last_riv_year
        logical use_small_rivers,use_coare, only_land_co2
        character*80 file_small_riv
        common/progn_run/progn_atm_path,res_atm_path,progn_use_mask,progn_write_mode
        logical rec_ocn,rec_atm,rec_ice,rec_lnd,
     &     rd_rec_ocn,rd_rec_atm,rd_rec_ice,rd_rec_lnd
        common/recstate/rec_ocn,rec_atm,rec_ice,rec_lnd,
     &     rd_rec_ocn,rd_rec_atm,rd_rec_ice,rd_rec_lnd
        common/rivers/use_rivers,nrivers
        common/rad/var_ocnalb,ice_lw_alb
        common/xtiming/cur_date,cur_sec
	common/pars/cld_power,use_ciaf,ciaf_vers,tidal,
     &  use_frac,use_frcp,use_tref,
     &  cld_lw_vers,cld_sw_vers,use_xie,use_xie4ncep,
     &  use_datahum,old_tides,aomip_rules,prognostic_run,solar_avg,
     &  use_small_rivers,file_small_riv,first_riv_year,last_riv_year,
     &  use_coare,only_land_co2
        character*80 input_path, tides_file
        common/input_data/input_path, tides_file
        namelist/cpl_nml/wap_file,use_ciaf,ciaf_vers,tidal,use_frac,
     &     use_frcp,use_tref,cld_lw_vers,cld_sw_vers,use_xie,use_xie4ncep,
     &     use_datahum,atm_nadv,td_nadv,old_tides,cld_power,
     &     aomip_rules,prognostic_run,progn_atm_path,res_atm_path,progn_use_mask,
     &     progn_write_mode,input_path,tides_file,solar_avg,use_rivers,
     &     use_small_rivers,file_small_riv,first_riv_year,last_riv_year,
     &     ncv_a2o,a2o_type,ncv_a2i,a2i_type,ncv_o2a,o2a_type,ncv_o2i,o2i_type,
     &     ncv_i2a,i2a_type,ncv_i2o,i2o_type,ncv_l2a,l2a_type,ncv_l2o,l2o_type,
     &     use_coare,setdate,startdate,finishdate,
     &     var_ocnalb,ice_lw_alb,rec_ocn,rec_ice,rec_atm,rec_lnd,
     &     rd_rec_ocn,rd_rec_atm,rd_rec_ice,rd_rec_lnd,only_land_co2

	integer*4 startyear,monlen(12),year,month,day
        integer*4 startdate,finishdate
        logical setdate
	data monlen/0,31,59,90,120,151,181,212,243,273,304,334/,old_date/0/,ds/86400.0/

      wap_file='../Output/cpl/state.html' ! html file for current state summery
      use_ciaf=.true.                     ! use CORE data
      ciaf_vers=1                         ! use 1 - CORE v1.0, 2 - CORE v2.0
      tidal=.false.                       ! use tidal forcing
      var_ocnalb = .false.                ! .true. if ocean albedo depends on angle of incidence,
                                          !    otherwise it is constant
      ice_lw_alb = 1.0                    ! this is justifying parameter regulating a portion of
                                          ! lw radiation accepted by the ice surface: ideally it is 1.0
      use_frac=.false.                    ! use combined lwdn equal for ocn and ice
      use_frcp=.false.                    ! multiply atm-ice rad fluxes by frcp
      use_xie=.true.                      ! use Xie & Arkin precipitation for Arctic
      use_xie4ncep=.true.                 ! use Xie & Arkin precipitation instead of NCEP/NCAR
      use_datahum=.false.                 ! use humidity data instead of rhum=0.9
      cld_lw_vers=1                       ! LW clouds: 1 - AOMIP, 2 - A.Gill
      cld_sw_vers=1                       ! SW clouds: 1 - AOMIP, 2 - Rosati&Miyakoda
      atm_nadv=8                          ! times a day to update atm state
      td_nadv=16                          ! times a day to update tidal state
      old_tides=.true.                    ! use old tide calculation mode (slow)
      cld_power=1.0                       ! use power for cloud data
      aomip_rules=.true.                  ! use AOMIP rules to calculate dw radiation
      prognostic_run=.false.              ! read 2000-2100 atmospheric data
      progn_atm_path=''                   ! directory path where prognostic atm files are placed
      res_atm_path=''                     ! storage subdirectory for interpolated atm forcing
      progn_use_mask=.true.               ! use area mask from atm data in interpolation routine
      progn_write_mode='w'                ! n -- just interpolate, w - store interpolated, 
                                          ! r - just read interpolated
      solar_avg=24.0                      ! average period for solar radiation in hours
      use_tref=.false.                    ! use t at ref z=2m instead of Tair at z0
      use_rivers=.false.
      use_small_rivers=.false.
      file_small_riv=''
      first_riv_year=0
      last_riv_year=0
      use_coare = .false.                 ! use COARE fluxes and drag coeffs instead of defaults
      setdate = .false.
      startdate = 19480101
      finishdate = 20160101
      
      only_land_co2 = .false.             ! allow CO2 fluxes from land only, i.e. interaction
                                          ! between atm and ocn, otherwise fluxes are prescribed
					  ! everywhere from lnd module

      rec_ocn = .false.
      rec_ice = .false.
      rec_atm = .false.
      rec_lnd = .false.
      rd_rec_ocn = .false.
      rd_rec_ice = .false.
      rd_rec_atm = .false.
      rd_rec_lnd = .false.

      ncv_a2o = 1
      a2o_type = cnv_exp
      ncv_a2i = 1
      a2i_type = cnv_exp
      ncv_a2l = 1
      a2l_type = cnv_exp
      ncv_o2a = 1
      o2a_type = cnv_exp
      ncv_o2i = 1
      o2i_type = cnv_exp
      ncv_i2a = 1
      i2a_type = cnv_exp
      ncv_i2o = 1
      i2o_type = cnv_exp
      ncv_l2a = 1
      l2a_type = cnv_exp
      ncv_l2o = 1
      l2o_type = cnv_exp

      pp=4.0*datan2(1.0,1.0)/180.0

      call init_idxs

      input_path='../../MassStore/'
      tides_file='tides/tides.dat'

      Open(20,file='../Input/cpl/parms.nml',status='old')
      Read(20,nml=cpl_nml)
      Close(20)
      
      ! Check input parameters consistensy
      !   If data are not recorded they cannot be read
      if(rec_ocn)rd_rec_ocn=.false.
      if(rec_atm)rd_rec_atm=.false.
      if(rec_ice)rd_rec_ice=.false.
      if(rec_lnd)rd_rec_lnd=.false.
            
      call clock_ini
      
      nu_diag=48
      if(nu_diag.ne.6)Open(unit=nu_diag,file='cpl.log',form='formatted')
      buff=0
      Write(nu_diag,*)'Prognostic Run: ',prognostic_run
      Write(nu_diag,*)'Use Rivers:     ',use_rivers
      Write(nu_diag,*)'Tides are from file ',trim(input_path)//trim(tides_file)

      call setup_mpi('cpl')    ! initial setup for message passing
      !================================================================
      ! Start components
      !================================================================

        if(setdate)then
           buff(1)=startdate
           buff(2)=finishdate
           startyear=int(startdate/10000)
           year=int(finishdate/10000)-startyear
           month=int((finishdate-(year+startyear)*10000)/100)
           day=finishdate-(year+startyear)*10000-month*100
           finish=(365*year+monlen(month)+day-1)*ds+Ist%sec
        else
           finish=1.0e+20
        endif

        Write(nu_diag,*) 'Start components'
	call send_i(buff, nbuff, msgtype_c2ii, ice_id, ierr)

	call send_i(buff, nbuff, msgtype_c2ai, atm_id, ierr)

	call send_i(buff, nbuff, msgtype_c2li, lnd_id, ierr)

        if(use_rivers)buff(5)=1
	call send_i(buff, nbuff, msgtype_c2oi, ocn_id, ierr)

	Write(nu_diag,*) 'Starting buffer is sent'

      !================================================================
      ! Receive description buffer and grid from ice
      !================================================================

	Write(nu_diag,*) 'Waiting ice model configuration'
	call recv_i(buff, nbuff, msgtype_i2ci, ice_id, ierr)
	Write(nu_diag,*) 'Got ice model configuration'
      if(buff(3).eq.1) stop_now=1   !stop now flag
      Ist%idate    = buff( 4)    ! initial date (coded: yyyymmdd)
      Ist%sec      = buff( 5)    ! elapsed seconds into date
      Ist%imt      = buff( 7)    ! grid size in x-direction
      Ist%jmt      = buff( 8)    ! grid size in y-direction
      Ist%nadv     = buff( 9)    ! ice steps per day

      solar_avg = max(solar_avg,24.0/dble(Ist%nadv))
c      print *,'Ice Solar_avg = ',solar_avg
	
      Write(nu_diag,'(a29,5(/,a29,i8))')
     * 'Received ice information:    ',
     * '    initial date (yyyymmdd): ', Ist%idate,
     * '            elapsed seconds: ', Ist%sec,
     * '   grid size in x-direction: ', Ist%imt,
     * '   grid size in y-direction: ', Ist%jmt,
     * '          ice steps per day: ', Ist%nadv
	
      
      startyear=int(Ist%idate/10000)
	year=int(Ist%idate/10000)-startyear
	month=int((Ist%idate-(year+startyear)*10000)/100)
	day=Ist%idate-(year+startyear)*10000-month*100
      current=(365*year+monlen(month)+day-1)*ds+Ist%sec
	Ist%next=current+ds/Ist%nadv
c	Write(nu_diag,'(a29,f15.7)')
c     * '       next date (yyyymmdd): ', Ist%next
	
c      if(stop_now == 1) call exit_coupler
	
      imt=Ist%imt
	jmt=Ist%jmt
      allocate(i_lon(imt,jmt),i_lat(imt,jmt),i_work(imt,jmt),
     &  i_tarea(imt,jmt),i_mask(imt,jmt))
	Igr%lon=>i_lon
	Igr%lat=>i_lat
	Igr%work=>i_work
	Igr%tarea=>i_tarea
	Igr%mask=>i_mask
	
      Write(nu_diag,*)'Waiting Ice grid'
      call get_grid(Igr, Ist, msgtype_i2ci, ice_id, ierr)
      Igr%lon = Igr%lon/pp
      Igr%lat = Igr%lat/pp
      Igr%tarea = Igr%tarea/pp/pp*111000.0*111000.0
      Write(nu_diag,*)'Got Ice grid'
      Write(nu_diag,*)'Ice Latitude : ',minval(i_lat),maxval(i_lat)
      Write(nu_diag,*)'Ice Longitude: ',minval(i_lon),maxval(i_lon)
      Write(nu_diag,*)'Ice Area     : ',minval(i_tarea)/1.0e+10,maxval(i_tarea)/1.0e+10

      !================================================================
	! Allocate Ice arrays and pointers
      !================================================================

	allocate(ailohi(imt,jmt),Tsrf(imt,jmt),alvdr(imt,jmt),
     &  alidr(imt,jmt),alvdf(imt,jmt),alidf(imt,jmt),Flatent(imt,jmt),
     &  Fsensible(imt,jmt),Flwout(imt,jmt),evap(imt,jmt),
     &  tauxa(imt,jmt),tauya(imt,jmt),Fswthru(imt,jmt),Fhnet(imt,jmt),
     &  Fresh(imt,jmt),Fsalt(imt,jmt),tauxo(imt,jmt),tauyo(imt,jmt),
     &  dtauo(imt,jmt),ztlvl(imt,jmt),Fwe(imt,jmt),
     &  Tref(imt,jmt), sabs(imt,jmt), zlvl(imt,jmt),uatm(imt,jmt),
     &  vatm(imt,jmt),potT(imt,jmt),Tair(imt,jmt),Qa(imt,jmt),
     &  rhoa(imt,jmt),swvdr(imt,jmt),swvdf(imt,jmt),swidr(imt,jmt),
     &  swidf(imt,jmt),Flw(imt,jmt),Frain(imt,jmt),Fsnow(imt,jmt),
     &  sst(imt,jmt),sss(imt,jmt),uocn(imt,jmt),vocn(imt,jmt),
     &  strtltx(imt,jmt),strtlty(imt,jmt),frzmlt(imt,jmt))
      I%ailohi=>ailohi
      I%Tsrf=>Tsrf
      I%Fwe=>Fwe
      I%alvdr=>alvdr
      I%alidr=>alidr
      I%alvdf=>alvdf
      I%alidf=>alidf
      I%Flatent=>Flatent
      I%Fsensible=>Fsensible
      I%Flwout=>Flwout
      I%evap=>evap
      I%tauxa=>tauxa
      I%tauya=>tauya
      I%Fswthru=>Fswthru
      I%Fhnet=>Fhnet
      I%Fresh=>Fresh
      I%Fsalt=>Fsalt
      I%tauxo=>tauxo
      I%tauyo=>tauyo
      I%dtauo=>dtauo
      I%Tref=>Tref
      I%sabs=>sabs
      ToI%zlvl=>zlvl
      ToI%ztlvl=>ztlvl
      ToI%uatm=>uatm
      ToI%vatm=>vatm
      ToI%potT=>potT
      ToI%Tair=>Tair
      ToI%Qa=>Qa
      ToI%rhoa=>rhoa
      ToI%swvdr=>swvdr
      ToI%swvdf=>swvdf
      ToI%swidr=>swidr
      ToI%swidf=>swidf
      ToI%Flw=>Flw
      ToI%Frain=>Frain
      ToI%Fsnow=>Fsnow
      ToI%sst=>sst
      ToI%sss=>sss
      ToI%uocn=>uocn
      ToI%vocn=>vocn
      ToI%strtltx=>strtltx
      ToI%strtlty=>strtlty
      ToI%frzmlt=>frzmlt
      

      !================================================================
      ! Receive description buffer and grid from ocean
      !================================================================

	Write(nu_diag,*) 'Waiting ocean model configuration'
	call recv_i(buff, nbuff, msgtype_o2ci, ocn_id, ierr)
	Write(nu_diag,*) 'Got ocean model configuration'
      if(buff(3).eq.1) stop_now=1   !stop now flag
      Ost%idate    = buff( 4)    ! initial date (coded: yyyymmdd)
      Ost%sec      = buff( 5)    ! elapsed seconds into date
      Ost%imt      = buff( 7)    ! grid size in x-direction
      Ost%jmt      = buff( 8)    ! grid size in y-direction
      Ost%nadv     = buff( 9)    ! ocean steps per day
      if(use_rivers)nrivers=buff(10)

      solar_avg = max(solar_avg,24.0/dble(Ost%nadv))
c      print *,'Ocn Solar_avg = ',solar_avg
	Write(nu_diag,'(a29,5(/,a29,i8))')
     * 'Received ocean information:    ',
     * '    initial date (yyyymmdd): ', Ost%idate,
     * '            elapsed seconds: ', Ost%sec,
     * '   grid size in x-direction: ', Ost%imt,
     * '   grid size in y-direction: ', Ost%jmt,
     * '        ocean steps per day: ', Ost%nadv

	year=int(Ost%idate/10000)-startyear
	month=int((Ost%idate-(year+startyear)*10000)/100)
	day=Ost%idate-(year+startyear)*10000-month*100
      if(current.ne.(365*year+monlen(month)+day-1)*ds+Ost%sec)then
	if(current.eq.0) then
         current=(365*year+monlen(month)+day-1)*ds+Ost%sec
	else
	   Write(nu_diag,*)'start time mismatch'
         call exit_coupler
	endif
	endif
	Ost%next=current+ds/Ost%nadv
c	Write(nu_diag,'(a29,f15.7)')
c     * '       next date (yyyymmdd): ', Ost%next

c	if(stop_now == 1) call exit_coupler

      imt=Ost%imt
	jmt=Ost%jmt
	allocate(o_lon(imt,jmt),o_lat(imt,jmt),o_work(imt,jmt),
     &  o_tarea(imt,jmt),o_mask(imt,jmt))
	Ogr%lon=>o_lon
	Ogr%lat=>o_lat
	Ogr%work=>o_work
	Ogr%tarea=>o_tarea
	Ogr%mask=>o_mask

      Write(nu_diag,*)'Waiting Ocn grid'
      call get_grid(Ogr, Ost, msgtype_o2ci, ocn_id, ierr)
      Ogr%lon = Ogr%lon/pp
      Ogr%lat = Ogr%lat/pp
      Write(nu_diag,*)'Got Ocn grid'
      Write(nu_diag,*)'Ocn Latitude : ',minval(o_lat),maxval(o_lat)
      Write(nu_diag,*)'Ocn Longitude: ',minval(o_lon),maxval(o_lon)
      Write(nu_diag,*)'Ocn Area     : ',minval(o_tarea)/1.0e+10,maxval(o_tarea)/1.0e+10


      !================================================================
	! Allocate Ocean arrays and pointers
      !================================================================

	allocate(o_sst(imt,jmt),o_sss(imt,jmt),o_uocn(imt,jmt),
     &  o_vocn(imt,jmt),o_strtltx(imt,jmt),o_strtlty(imt,jmt),
     &  o_frzmlt(imt,jmt),o_Fswthru(imt,jmt),o_Fhnet(imt,jmt),
     &  o_Fprec(imt,jmt),o_Fevap(imt,jmt),o_tauxo(imt,jmt),
     &  o_tauyo(imt,jmt),o_dtauo(imt,jmt),o_Fmelt(imt,jmt),
     &  o_Fswthru_ice(imt,jmt), o_Fhnet_ice(imt,jmt), 
     &  o_tauxoi(imt,jmt), o_tauyoi(imt,jmt),o_Fsalt(imt,jmt),
     &  o_utidal(imt,jmt), o_vtidal(imt,jmt), o_sriv(imt,jmt),
     &  o_dtauoi(imt,jmt), o_Frac(imt,jmt), o_dFhnet(imt,jmt),
     &  o_co2(imt,jmt), o_fco2(imt,jmt), o_Fwe(imt,jmt))
      O%sst=>o_sst
      O%sss=>o_sss
      O%uocn=>o_uocn
      O%vocn=>o_vocn
      O%tiltx=>o_strtltx
      O%tilty=>o_strtlty
      O%frzmlt=>o_frzmlt
      O%cco2=>o_co2
      ToO%Fswthru=>o_Fswthru
      ToO%Fhnet=>o_Fhnet
      ToO%Fprec=>o_Fprec
      ToO%Fevap=>o_Fevap
      ToO%Fmelt=>o_Fmelt
      ToO%Fsalt=>o_Fsalt
      ToO%tauxo=>o_tauxo
      ToO%tauyo=>o_tauyo
      ToO%dtauo=>o_dtauo
      ToO%Fswthru_ice=>o_Fswthru_ice
      ToO%Fhnet_ice=>o_Fhnet_ice
      ToO%tauxoi=>o_tauxoi
      ToO%tauyoi=>o_tauyoi
      ToO%dtauoi=>o_dtauoi
      ToO%Frac=>o_Frac
      ToO%dFhnet=>o_dFhnet
      ToO%utidal=>o_utidal
      ToO%vtidal=>o_vtidal
      ToO%Frivprec=>o_sriv
      ToO%Fco2=>o_fco2
      ToO%Fwe=>o_Fwe

      !================================================================
      ! Receive description buffer and grid from atmosphere
      !================================================================

	Write(nu_diag,*) 'Waiting atmosphere model configuration'
	call recv_i(buff, nbuff, msgtype_a2ci, atm_id, ierr)
	Write(nu_diag,*) 'Got atmosphere model configuration'
      if(buff(3).eq.1) stop_now=1   !stop now flag
      Ast%idate    = buff( 4)    ! initial date (coded: yyyymmdd)
      Ast%sec      = buff( 5)    ! elapsed seconds into date
      Ast%imt      = buff( 7)    ! grid size in x-direction
      Ast%jmt      = buff( 8)    ! grid size in y-direction
      Ast%nadv     = buff( 9)    ! atm steps per day

      solar_avg = max(solar_avg,24.0/dble(Ast%nadv))
	Write(nu_diag,'(a29,5(/,a29,i8))')
     * 'Received atm information   : ',
     * '    initial date (yyyymmdd): ', Ast%idate,
     * '            elapsed seconds: ', Ast%sec,
     * '   grid size in x-direction: ', Ast%imt,
     * '   grid size in y-direction: ', Ast%jmt,
     * '   atmosphere steps per day: ', Ast%nadv

      year=int(Ast%idate/10000)-startyear
      month=int((Ast%idate-(year+startyear)*10000)/100)
      day=Ast%idate-(year+startyear)*10000-month*100

      if(current.ne.(365*year+monlen(month)+day-1)*ds+Ast%sec)then
	if(current.eq.0) then
           current=(365*year+monlen(month)+day-1)*ds+Ast%sec
	else
	   Write(nu_diag,*)'start time mismatch'
           call exit_coupler
	endif
      endif
      Ast%next=current+ds/Ast%nadv

c      if(stop_now == 1) call exit_coupler

      imt=Ast%imt
      jmt=Ast%jmt
      allocate(a_lon(imt,jmt),a_lat(imt,jmt),a_work(imt,jmt),
     &  a_tarea(imt,jmt),a_mask(imt,jmt))
      Agr%lon=>a_lon
      Agr%lat=>a_lat
      Agr%work=>a_work
      Agr%tarea=>a_tarea
      Agr%mask=>a_mask

      Write(nu_diag,*)'Waiting Atm grid'
      call get_grid(Agr, Ast, msgtype_a2ci, atm_id, ierr)
c      Agr%mask(:,:)=Agr%mask(:,:)+1

      Write(nu_diag,*)'Got Atm grid'
      Write(nu_diag,*)'Atm Latitude : ',minval(a_lat),maxval(a_lat)
      Write(nu_diag,*)'Atm Longitude: ',minval(a_lon),maxval(a_lon)
      Write(nu_diag,*)'Atm Area     : ',minval(a_tarea)/1.0e+10,maxval(a_tarea)/1.0e+10

      !================================================================
      ! Allocate Atmosphere arrays and pointers
      !================================================================

	allocate(a_zlvl(imt,jmt),a_uatm(imt,jmt),a_vatm(imt,jmt),
     &  a_potT(imt,jmt),a_Tair(imt,jmt),a_Qa(imt,jmt),a_rhoa(imt,jmt),
     &  a_Fsw(imt,jmt),a_Flw(imt,jmt),a_Frain(imt,jmt),
     &  a_Fsnow(imt,jmt),a_cld(imt,jmt),a_small_riv(imt,jmt),
     &  a_ztlvl(imt,jmt),a_tsurf(imt,jmt),
     &  a_Cu(imt,jmt),a_Ct(imt,jmt),a_Ce(imt,jmt),a_albed(imt,jmt),
     &  a_dicec(imt,jmt),a_co2(imt,jmt),a_fco2(imt,jmt))

      A%zlvl=>a_zlvl
      A%ztlvl=>a_ztlvl
      A%uatm=>a_uatm
      A%vatm=>a_vatm
      A%potT=>a_potT
      A%Tair=>a_Tair
      A%Qa=>a_Qa
      A%rhoa=>a_rhoa
      A%Fsw=>a_Fsw
      A%Flw=>a_Flw
      A%Frain=>a_Frain
      A%Fsnow=>a_Fsnow
      A%cld=>a_cld
      A%small_riv=>a_small_riv
      A%cco2=>a_co2

      ToA%tsurf=>a_tsurf
      ToA%Cu=>a_Cu
      ToA%Ct=>a_Ct
      ToA%Ce=>a_Ce
      ToA%albed=>a_albed
      ToA%dicec=>a_dicec
      ToA%Fco2=>a_fco2

      !================================================================
      ! Receive description buffer and grid from land
      !================================================================

	Write(nu_diag,*) 'Waiting land model configuration'
	call recv_i(buff, nbuff, msgtype_l2ci, lnd_id, ierr)
	Write(nu_diag,*) 'Got land model configuration'
      if(buff(3).eq.1) stop_now=1   !stop now flag
      Lst%idate    = buff( 4)    ! initial date (coded: yyyymmdd)
      Lst%sec      = buff( 5)    ! elapsed seconds into date
      Lst%imt      = buff( 7)    ! grid size in x-direction
      Lst%jmt      = buff( 8)    ! grid size in y-direction
      Lst%nadv     = buff( 9)    ! land steps per day

	Write(nu_diag,'(a29,5(/,a29,i8))')
     * 'Received lnd information   : ',
     * '    initial date (yyyymmdd): ', Lst%idate,
     * '            elapsed seconds: ', Lst%sec,
     * '   grid size in x-direction: ', Lst%imt,
     * '   grid size in y-direction: ', Lst%jmt,
     * '   atmosphere steps per day: ', Lst%nadv

      year=int(Lst%idate/10000)-startyear
      month=int((Lst%idate-(year+startyear)*10000)/100)
      day=Lst%idate-(year+startyear)*10000-month*100

      if(current.ne.(365*year+monlen(month)+day-1)*ds+Lst%sec)then
	if(current.eq.0) then
           current=(365*year+monlen(month)+day-1)*ds+Lst%sec
	else
	   Write(nu_diag,*)'start time mismatch'
           call exit_coupler
	endif
      endif
      Lst%next=current+ds/Lst%nadv

c      if(stop_now == 1) call exit_coupler

      imt=Lst%imt
      jmt=Lst%jmt
      allocate(ln_lon(imt,jmt),ln_lat(imt,jmt),ln_wrk(imt,jmt),
     &  ln_tarea(imt,jmt),ln_mask(imt,jmt))
      Lgr%lon=>ln_lon
      Lgr%lat=>ln_lat
      Lgr%work=>ln_wrk
      Lgr%tarea=>ln_tarea
      Lgr%mask=>ln_mask

      Write(nu_diag,*)'Waiting Lnd grid'
      call get_grid(Lgr, Lst, msgtype_l2ci, lnd_id, ierr)

      Write(nu_diag,*)'Got Lnd grid'
      Write(nu_diag,*)'Lnd Latitude : ',minval(ln_lat),maxval(ln_lat)
      Write(nu_diag,*)'Lnd Longitude: ',minval(ln_lon),maxval(ln_lon)
      Write(nu_diag,*)'Lnd Area     : ',minval(ln_tarea)/1.0e+10,maxval(ln_tarea)/1.0e+10

      !================================================================
      ! Allocate Land arrays and pointers
      !================================================================

	allocate(ln_NatCO2(imt,jmt),ln_AntCO2(imt,jmt),ln_FosCO2(imt,jmt),
     &  ln_Work(imt,jmt))

      FromL%NatCO2=>ln_NatCO2
      FromL%AntCO2=>ln_AntCO2
      FromL%FosCO2=>ln_FosCO2

      ToL%Work=>ln_Work

      Write(nu_diag,*) 'Start receiving initial states'
      !================================================================
      ! Get initial Atm-Ocn-Ice-Lnd state
      !================================================================

      Write(nu_diag,*)'Waiting Ocn state'
      call get_ocean_state(O, Ost,nrivers)
      Write(nu_diag,*) 'Ocean state received'
      if(use_rivers)then
         Write(nu_diag,*) 'External river transport is on. Number of rivers: ',nrivers
         allocate(rivcoords(2*nrivers),rtrans(nrivers))
         Write(nu_diag,*) 'Waiting river coordinates'
         call recv_d(rivcoords, 2*nrivers, msgtype_o2ci, ocn_id, ierr)
         Write(nu_diag,*) 'Got river coordinates'
         do k=1,nrivers
            write(nu_diag,*)rivcoords(2*k-1),rivcoords(2*k)
         enddo
      endif
      year=int(Ost%idate/10000)-startyear
      month=int((Ost%idate-(year+startyear)*10000)/100)
      day=Ost%idate-(year+startyear)*10000-month*100
      current=(365*year+monlen(month)+day-1)*ds+Ost%sec
      Ost%next=current+ds/Ost%nadv
	Write(nu_diag,'(a15,f15.7)')
     * 'Ocean next date: ', Ost%next

      Write(nu_diag,*)'Sending Ocn state to Ice'
      call send_ocean_init_to_ice( ToI,
     &               I, Igr, Ist,
     &               O, Ogr, Ost)
      Write(nu_diag,*)'Done'

      Write(nu_diag,*)'Waiting Ice state'
      call get_ice_state(I, Ist)
      Write(nu_diag,*) 'Ice state received'
      year=int(Ist%idate/10000)-startyear
      month=int((Ist%idate-(year+startyear)*10000)/100)
      day=Ist%idate-(year+startyear)*10000-month*100
      current=(365*year+monlen(month)+day-1)*ds+Ist%sec
      Ist%next=current+ds/Ist%nadv
      Write(nu_diag,'(a15,f15.7)')
     * 'Ice next date: ', Ist%next

      Write(nu_diag,*)'Waiting Lnd state'
      call get_lnd_state(FromL, Lst)
      Write(nu_diag,*) 'Lnd state received'
      year=int(Lst%idate/10000)-startyear
      month=int((Lst%idate-(year+startyear)*10000)/100)
      day=Lst%idate-(year+startyear)*10000-month*100
      current=(365*year+monlen(month)+day-1)*ds+Lst%sec
      Lst%next=current+ds/Lst%nadv
      Write(nu_diag,'(a15,f15.7)')
     * 'Lnd next date: ', Lst%next

      Write(nu_diag,*)'Waiting Atm state'
      call get_atm_state(A, Ast)
      Write(nu_diag,*) 'Atm state received'
      year=int(Ast%idate/10000)-startyear
      month=int((Ast%idate-(year+startyear)*10000)/100)
      day=Ast%idate-(year+startyear)*10000-month*100
      current=(365*year+monlen(month)+day-1)*ds+Ast%sec
      Ast%next=current+ds/Ast%nadv
      Write(nu_diag,'(a15,f15.7)')
     * 'Atm next date: ', Ast%next
      Write(nu_diag,*)Ast%idate,Ast%nadv

      solar_avg = max(solar_avg,24.0/dble(Ast%nadv))
c      print *,'Atm Solar_avg = ',solar_avg

c      Write(nu_diag,*)'Waiting Atm state'
c      call get_atm_state(current, A, Agr, Ast)
      Write(nu_diag,*) 'Atmosphere state received'
      
      if(use_rivers)then
         call get_rivers(rivcoords,rtrans,nrivers,Ist%idate,Ist%sec)
      endif

      if(tidal)then
      allocate(etide(Ist%imt,Ist%jmt),utide(Ist%imt,Ist%jmt),
     *      vtide(Ist%imt,Ist%jmt),xtide(Ist%imt,Ist%jmt),
     *      ytide(Ist%imt,Ist%jmt))
      call tides(Ist%idate,Ist%sec,Ist%imt,Ist%jmt,
     *      Igr%mask,Igr%lon,Igr%lat,etide,utide,vtide,
     *      trim(input_path)//trim(tides_file))
      td_next = current+ds/td_nadv
      endif

c            call clock_print_time_table('prepare ice ')
c	call prepare_ice_support(   ToI,
c     &               I, Igr, Ist,
c     &               O, Ogr, Ost,
c     &               FromL, Lgr, Lst,
c     &               A, Agr, Ast   )

            call clock_print_time_table('prepare ocn ')
	call prepare_ocean_support(   ToO,
     &               I, Igr, Ist,
     &               O, Ogr, Ost,
     &               FromL, Lgr, Lst,
     &               A, Agr, Ast   )

            call clock_print_time_table('prepare lnd ')
	call prepare_land_support(   ToL,
     &               I, Igr, Ist,
     &               O, Ogr, Ost,
     &               FromL, Lgr, Lst,
     &               A, Agr, Ast   )

            call clock_print_time_table('prepare atm ')
	call prepare_atmos_support(   ToA,
     &               I, Igr, Ist,
     &               O, Ogr, Ost,
     &               FromL, Lgr, Lst,
     &               A, Agr, Ast   )

c      Write(nu_diag,*)'Sending data to Ice'
c      call send_ice_data(ToI, Igr, Ist)
c      Write(nu_diag,*)'Done'
      Write(nu_diag,*)'Sending data to Ocn'
      call send_ocean_data(current, ToO, Ogr, Ost)
      Write(nu_diag,*)'Done'
      if(use_rivers)then
         write(nu_diag,*)'sending river transports ',nrivers
         call send_d(rtrans, nrivers, msgtype_c2o, ocn_id, ierr)
         Write(nu_diag,*)'Done'
      endif
      Write(nu_diag,*)'Sending data to Lnd'
      call send_land_data(ToL, Lgr, Lst)
      Write(nu_diag,*)'Done'
      Write(nu_diag,*)'Sending data to Atm'
      call send_atmos_data(ToA, Agr, Ast)
      Write(nu_diag,*)'Done'

	Open(20,file='../Input/cpl/control',form='formatted')
	Read(20,*)ctrl
	Close(20)
	if(ctrl.eq.0)stop_now=1

        atm_mpt = 0.
        ocn_mpt = 0.
        ice_mpt = 0.
        lnd_mpt = 0.
c============================================================================================
	Write(nu_diag,*) 'Start running'
        running=.true.

      Do while(running)

	next=min(Ost%next,Ist%next,Ast%next)
        if(next.gt.finish)call exit_coupler
c        Write(nu_diag,*)'ice ',Ist%next-next,',  ocn ',Ost%next-next,', atm ',Ast%next-next
        
      if(next.ge.Ist%next)then
         cur_date = Ist%idate
         cur_sec  = Ist%sec  
      endif
      if(next.ge.Ost%next)then
         cur_date = Ost%idate
         cur_sec  = Ost%sec  
      endif
      if(next.ge.Ast%next)then
         cur_date = Ast%idate
         cur_sec = Ast%sec
      endif
      if(next.ge.Lst%next)then
         cur_date = Lst%idate
         cur_sec = Lst%sec
      endif

      if(old_date.ne.cur_date)then
         if(reftime.ne.0)call clock_print_time_table('finish ')
         call clock_start_time_table
         Write(nu_diag,*)'Date: ',cur_date
         old_date=cur_date
      endif

      if(tidal.and.(.not.old_tides))then
         if(next.ge.Ist%next.and.next.ge.td_next)then
            call clock_starttimer(9)
            call clock_print_time_table('updating tides ')
            ice_day=int(Ist%next/ds)+1
            ice_sec=int(Ist%next-(ice_day-1)*ds)
            ice_year=startyear+int(ice_day/365)
            ice_day=ice_day-(ice_year-startyear)*365
            do ice_month=12,1,-1
               if(ice_day.gt.monlen(ice_month))exit
            enddo
            ice_day=ice_day-monlen(ice_month)
            ice_date=10000*ice_year+100*ice_month+ice_day
c            Write(nu_diag,*)'Waiting tides update'
            call tides(ice_date,ice_sec,Ist%imt,Ist%jmt,
     *         Igr%mask,Igr%lon,Igr%lat,etide,utide,vtide,
     *         trim(input_path)//trim(tides_file))
c            Write(nu_diag,*) 'Tides are updated'
            td_next = td_next+ds/td_nadv
         endif
      endif

      call clock_starttimer(8)
      call clock_print_time_table('printing out ')
	Open(21,file='../Output/cpl/current.date',form='formatted')
        write(timestring,
     * "(I2.2,'/',I2.2,'/',I4.4,'   ',I2.2,':',I2.2,':',I2.2)")
     * int(Ist%idate-100*int(Ist%idate/100)),
     * int((Ist%idate-10000*int(Ist%idate/10000))/100),
     * int(Ist%idate/10000), int(Ist%sec/3600),
     * int((Ist%sec-3600*int(Ist%sec/3600))/60),
     * int(Ist%sec-60*int(Ist%sec/60))
	write(21,"(a)")timestring
        write(21,'(a18,f12.3)')'Next ice time is :',Ist%next
	write(21,
     * "(I2.2,'/',I2.2,'/',I4.4,'   ',I2.2,':',I2.2,':',I2.2)")
     * int(Ost%idate-100*int(Ost%idate/100)),
     * int((Ost%idate-10000*int(Ost%idate/10000))/100),
     * int(Ost%idate/10000), int(Ost%sec/3600),
     * int((Ost%sec-3600*int(Ost%sec/3600))/60),
     * int(Ost%sec-60*int(Ost%sec/60))
        write(21,'(a18,f12.3)')'Next ocn time is :',Ost%next
        write(21,'(a18,f12.3)')'Next cpl time is :',next
        Close(21)
	Open(21,file=wap_file,form='formatted')
        Write(21,*)'<html>'
        Write(21,*)'<head>'
        Write(21,'(a)')'  <title>Ice-Ocean State: '//timestring//' </title>'
        Write(21,*)'  <meta http-equiv="refresh" content="30" >'
        Write(21,*)'</head>'
        Write(21,*)'  <body>'
        Write(21,*)'<h4>Current time:</h4>'
	Write(21,'(a)')'<h4>'//timestring//'</h4>'
        ntot=0; Stot=0.0; Ttot=0.0
        do m=1,Ast%imt
           do n=1,Ast%jmt
              if(Agr%mask(m,n).ne.0.and.Agr%lat(m,n).gt.87.5)then
                 Stot=Stot+A%potT(m,n)
                 Ttot=Ttot+A%Tair(m,n)
                 ntot=ntot+1
              endif
           enddo
        enddo
        Write(21,'(a10,f7.3,a2,f7.3,a6)')
     *     '<h4>Tair: ',Stot/dble(ntot)-273.16,' (',
     *     Ttot/dble(ntot)-273.16,')</h4>'
        ntot=0; Stot=0.0
        do m=1,Ist%imt
           do n=1,Ist%jmt
              if(Igr%mask(m,n).ne.0.and.Igr%lat(m,n).gt.89.0)then
                 Stot=Stot+I%Tsrf(m,n)
                 ntot=ntot+1
              endif
           enddo
        enddo
        Write(21,'(a10,f7.3,a5)')
     *     '<h4>Tsrf: ',Stot/dble(ntot)-273.16,'</h4>'
        ntot=0; Stot=0.0; Ttot=0.0
        do m=1,Ost%imt
           do n=1,Ost%jmt
              if(Ogr%mask(m,n).ne.0.and.Ogr%lat(m,n).gt.89.0)then
                 Stot=Stot+O%sss(m,n)
                 Ttot=Ttot+O%sst(m,n)
                 ntot=ntot+1
              endif
           enddo
        enddo
        Write(21,'(a10,f7.3,a5)')
     *     '<h4>SST:  ',Ttot/dble(ntot)-273.16,'</h4>'
        Write(21,'(a10,f7.3,a5)')
     *     '<h4>SSS:  ',Stot/dble(ntot),'</h4>'
        Write(21,*)'</body>'
        Write(21,*)'</html>'
        Close(21)
	Open(21,file='../Output/cpl/timers',form='formatted')
	write(21,
     * "('Calculating the atm state               : ',F6.2)")ptimers(1)
	write(21,
     * "('Wait for comm and receive the ice state : ',F6.2)")ptimers(2)
	write(21,
     * "('Wait for comm and receive the ocn state : ',F6.2)")ptimers(3)
	write(21,
     * "('Prepare ice support                     : ',F6.2)")ptimers(4)
	write(21,
     * "('Prepare ocn support                     : ',F6.2)")ptimers(5)
	write(21,
     * "('Prepare tides for ice                   : ',F6.2)")ptimers(9)
	write(21,
     * "('Wait for comm and send ice support      : ',F6.2)")ptimers(6)
	write(21,
     * "('Wait for comm and send ocn support      : ',F6.2)")ptimers(7)
	write(21,
     * "('Print out text files                    : ',F6.2)")ptimers(8)
        Close(21)

      if(next.ge.Ist%next)then
         call clock_starttimer(4)
         call clock_print_time_table('prepare ice data ')

         rtime0=dble(newtime-reftime)/dble(timescale)
         call prepare_ice_support(   ToI,
     &               I, Igr, Ist,
     &               O, Ogr, Ost,
     &               FromL, Lgr, Lst,
     &               A, Agr, Ast   )
         call clock_starttimer(6)
c         Write(nu_diag,*)'Sending data to Ice'
         call clock_print_time_table('sending ice data ')
         rtime1=dble(newtime-reftime)/dble(timescale)
c         Write(nu_diag,*)'ICE Preparation time: ',rtime1-rtime0
         ice_mpt=(ice_mpt*dble(Ist%nadv*30-1)+(rtime1-rtime0))/dble(Ist%nadv*30)

         call send_ice_data(ToI, Igr, Ist)
c         Write(nu_diag,*)'Done'
         cur_date=Ist%idate
         cur_sec=Ist%sec
      endif
      if(next.ge.Ost%next)then
         call clock_starttimer(3)
         call clock_print_time_table('getting ocean state ')
c         Write(nu_diag,*)'Waiting Ocn state'
         call get_ocean_state(O, Ost,nriv0)
c         Write(nu_diag,*) 'Ocean state received'
         cur_date=Ost%idate
         cur_sec=Ost%sec
      endif
      if(next.ge.Ast%next)then
         call clock_print_time_table('getting atmosphere state ')
c         Write(nu_diag,*)'Waiting Atm state'
         call get_atm_state(A, Ast)
c         Write(nu_diag,*) 'Atm state received'
         cur_date=Ast%idate
         cur_sec=Ast%sec
      endif
      if(next.ge.Lst%next)then
         call clock_print_time_table('getting land state ')
c         Write(nu_diag,*)'Waiting Lnd state'
         call get_lnd_state(FromL, Lst)
c         Write(nu_diag,*) 'Lnd state received'
         cur_date=Lst%idate
         cur_sec=Lst%sec
      endif
      current=next
      if(use_rivers)then
         call clock_print_time_table('getting rivers ')
         call get_rivers(rivcoord,rtrans,nrivers,cur_date,cur_sec)
      endif
c      print *,'Got from ocn ',O%uocn(105,157)
      
	Open(20,file='../Input/cpl/control',form='formatted')
	Read(20,*)ctrl
	Close(20)
	if(ctrl.eq.0)stop_now=1

c      print *,'*************',next,Ist%next,Ost%next
      if(next.ge.Ist%next)then
         call clock_starttimer(2)
         call clock_print_time_table('getting ice state ')
c         Write(nu_diag,*)'Waiting Ice state'
         call get_ice_state(I, Ist)
c         Write(nu_diag,*) 'Ice state received'
         Ist%next=current+ds/Ist%nadv
         if(Ist%next.gt.finish.or.stop_now.eq.1) Ist%next=-999
         if(Ist%next.lt.-100) Write(nu_diag,*)'Stop signal is to be send to ice'
      endif
      if(next.ge.Ost%next)then
         call clock_starttimer(5)
         call clock_print_time_table('prepare ocn data ')

         rtime0=dble(newtime-reftime)/dble(timescale)
         call prepare_ocean_support(   ToO,
     &               I, Igr, Ist,
     &               O, Ogr, Ost,
     &               FromL, Lgr, Lst,
     &               A, Agr, Ast   )

         call clock_starttimer(7)
c         Write(nu_diag,*)'Sending data to Ocn'
         call clock_print_time_table('sending ocn data ')
         rtime1=dble(newtime-reftime)/dble(timescale)
c         Write(nu_diag,*)'OCN Preparation time: ',rtime1-rtime0
         ocn_mpt=(ocn_mpt*dble(Ost%nadv*30-1)+(rtime1-rtime0))/dble(Ost%nadv*30)

         call send_ocean_data(current, ToO, Ogr, Ost)
c         Write(nu_diag,*)'Done'
         if(use_rivers)then
            call clock_print_time_table('sending river transport ')
c            write(nu_diag,*)'sending river transports ',nrivers
            call send_d(rtrans, nrivers, msgtype_c2o, ocn_id, ierr)
c            Write(nu_diag,*)'Done'
         endif
         Ost%next=current+ds/Ost%nadv
         if(Ost%next.gt.finish.or.stop_now.eq.1) Ost%next=-999
         if(Ost%next.lt.-100) Write(nu_diag,*)'Stop signal is to be send to ocn'
      endif
      if(next.ge.Ast%next)then
         call clock_print_time_table('prepare atm data ')

         rtime0=dble(newtime-reftime)/dble(timescale)
         call prepare_atmos_support(   ToA,
     &               I, Igr, Ist,
     &               O, Ogr, Ost,
     &               FromL, Lgr, Lst,
     &               A, Agr, Ast   )
c         Write(nu_diag,*)'Sending data to Atm'
         call clock_print_time_table('sending atm data ')
         rtime1=dble(newtime-reftime)/dble(timescale)
c         Write(nu_diag,*)'ATM Preparation time: ',rtime1-rtime0
         atm_mpt=(atm_mpt*dble(Ast%nadv*30-1)+(rtime1-rtime0))/dble(Ast%nadv*30)

         call send_atmos_data(ToA, Agr, Ast)
c         Write(nu_diag,*)'Done'
         Ast%next=current+ds/Ast%nadv
         if(Ast%next.gt.finish.or.stop_now.eq.1) Ast%next=-999
         if(Ast%next.lt.-100) Write(nu_diag,*)'Stop signal is to be send to atm'
      endif

      if(next.ge.Lst%next)then
         call clock_print_time_table('prepare lnd data ')

         rtime0=dble(newtime-reftime)/dble(timescale)
         call prepare_land_support(   ToA,
     &               I, Igr, Ist,
     &               O, Ogr, Ost,
     &               A, Agr, Ast   )
c         Write(nu_diag,*)'Sending data to Lnd'
         call clock_print_time_table('sending Lnd data ')
         rtime1=dble(newtime-reftime)/dble(timescale)
c         Write(nu_diag,*)'LND Preparation time: ',rtime1-rtime0
         lnd_mpt=(lnd_mpt*dble(Lst%nadv*30-1)+(rtime1-rtime0))/dble(Lst%nadv*30)

         call send_land_data(ToL, Lgr, Lst)
c         Write(nu_diag,*)'Done'
         Lst%next=current+ds/Lst%nadv
         if(Lst%next.gt.finish.or.stop_now.eq.1) Lst%next=-999
         if(Lst%next.lt.-100) Write(nu_diag,*)'Stop signal is to be send to lnd'
      endif

c      Write(nu_diag,*)'Mean ATM Preparation time: ',atm_mpt
c      Write(nu_diag,*)'Mean OCN Preparation time: ',ocn_mpt
c      Write(nu_diag,*)'Mean ICE Preparation time: ',ice_mpt
c      Write(nu_diag,*)'Mean LND Preparation time: ',lnd_mpt

c      Write(nu_diag,*)'Ast ',Ast%next,', Ost ',Ost%next,', Ist ',Ist%next
      running=(Ast%next.ge.0 .or. Ist%next.ge.0 .or. Ost%next.ge.0 .or. Lst%next.ge.0 .or.
     * Ast%next.lt.-100 .or. Ist%next.lt.-100 .or. Ost%next.lt.-100 .or. Lst%next.lt.-100)

      enddo

      if(rec_ocn.or.rd_rec_ocn) close(127)
      if(rec_atm.or.rd_rec_atm) close(129)
      if(rec_ice.or.rd_rec_ice) close(131)
      if(rec_lnd.or.rd_rec_lnd) close(133)

	call exit_coupler

      end


      subroutine get_ice_state(V, S)

      use struct
      use ice_indexes
      use msgtypes
	Type (state) S
	Type (from_ice) V

      integer cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs(20)
      common/allids/cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs
      integer ::
     &  cpl_task
     &, my_task           ! master task for ice
     &, MPI_COMM_MY
      common/tasks/my_task, cpl_task, MPI_COMM_MY
      integer nu_diag
      common/iodiag/nu_diag
      integer*4, parameter ::
     &   nbuff = 100
	integer*4, dimension(nbuff) ::
     &   buff
      integer*4 ierr,stop_now,nsend,ice_recno
      common/stop_flag/stop_now
	real*8, allocatable :: sbuff(:,:),ssbuff(:,:)
        logical rec_ocn,rec_atm,rec_ice,rec_lnd,
     &     rd_rec_ocn,rd_rec_atm,rd_rec_ice,rd_rec_lnd
        common/recstate/rec_ocn,rec_atm,rec_ice,rec_lnd,
     &     rd_rec_ocn,rd_rec_atm,rd_rec_ice,rd_rec_lnd

      if(S%next.lt.0) return

c      Write(nu_diag,*) 'Waiting ice model configuration'
	call recv_i(buff, nbuff, msgtype_i2c, ice_id, ierr)
      if(buff(3).eq.1) stop_now=1   !stop now flag
      S%idate = buff(4)          ! date (coded: yyyymmdd)
      S%sec   = buff(5)            ! elapsed seconds into date
c      S%nadv  = buff(9)          ! ice comm pairs per day
      nsend   = buff(10)          ! size of real buffer to send
c	Write(nu_diag,'(a29,4(/,a29,i8))')
c     * 'Received ice information:    ',
c     * '            date (yyyymmdd): ', S%idate,
c     * '            elapsed seconds: ', S%sec,
c     * 'size of real buffer to send: ', nsend
c	if(stop_now == 1) call exit_coupler
c      allocate(sbuff(nsend,nsnd))

      if(rec_ice .or. rd_rec_ice)then
         ice_recno=dayofyear(S%idate)
         ice_step=S%sec*S%nadv/86400
         ice_recno=(ice_recno-1)*S%nadv+ice_step+1
         Write(nu_diag,*) 'Receive ice model data for record ',ice_recno
c         Write(nu_diag,*) 'Record parameters: ',nsnd,nsend,nsnd*nsend
         allocate(ssbuff(nsnd,S%imt*S%jmt))
      endif

      allocate(sbuff(nsnd,nsend))

c      Write(nu_diag,*) 'Receive ice model data'
	call recv_d(sbuff, nsnd*nsend, msgtype_i2c, ice_id, ierr)

      do j=1,S%jmt
      do i=1,S%imt
         V%ailohi      (i,j) = 0.0
         V%Tsrf        (i,j) = 0.0
         V%alvdr       (i,j) = 0.0
         V%alidr       (i,j) = 0.0
         V%alvdf       (i,j) = 0.0
         V%alidf       (i,j) = 0.0
         V%Flatent     (i,j) = 0.0
         V%Fsensible   (i,j) = 0.0
         V%Flwout      (i,j) = 0.0
         V%evap        (i,j) = 0.0
         V%tauxa       (i,j) = 0.0
         V%tauya       (i,j) = 0.0
         V%Fswthru     (i,j) = 0.0
         V%Fhnet       (i,j) = 0.0
         V%Fresh       (i,j) = 0.0
         V%Fsalt       (i,j) = 0.0
         V%tauxo       (i,j) = 0.0
         V%tauyo       (i,j) = 0.0
         V%dtauo       (i,j) = 0.0
         V%Tref        (i,j) = 0.0
         V%sabs        (i,j) = 0.0
         V%Fwe         (i,j) = 0.0
      enddo
      enddo
      if(rec_ice)ssbuff=0.0

      if(.not.rd_rec_ice)then
      
      do n2=1,nsend
        i = mod(nint(sbuff(nsnd,n2))-1,S%imt)+1
        j = (nint(sbuff(nsnd,n2))-1)/S%imt+1
          V%ailohi      (i,j)= sbuff(kbufs_ifrc, n2)  ! none
          V%Tsrf        (i,j)= sbuff(kbufs_tsfc, n2)  ! K
          V%alvdr       (i,j)= sbuff(kbufs_alb1, n2)  ! none
          V%alidr       (i,j)= sbuff(kbufs_alb2, n2)  ! none
          V%alvdf       (i,j)= sbuff(kbufs_alb3, n2)  ! none
          V%alidf       (i,j)= sbuff(kbufs_alb4, n2)  ! none
          V%Flatent     (i,j)= sbuff(kbufs_lat , n2)  ! W/m^2
          V%Fsensible   (i,j)= sbuff(kbufs_sens, n2)  ! W/m^2
          V%Flwout      (i,j)= sbuff(kbufs_lwup, n2)  ! W/m^2
          V%evap        (i,j)= sbuff(kbufs_evap, n2)  ! kg/m^2 s
          V%tauxa       (i,j)= sbuff(kbufs_tauxa,n2)  ! N/m^2
          V%tauya       (i,j)= sbuff(kbufs_tauya,n2)  ! N/m^2
          V%Fswthru     (i,j)= sbuff(kbufs_netsw,n2)  ! W/m^2
          V%Fhnet       (i,j)= sbuff(kbufs_melth,n2)  ! W/m^2
          V%Fresh       (i,j)= sbuff(kbufs_meltw,n2)  ! kg/m^2 s
          V%Fsalt       (i,j)= sbuff(kbufs_salt ,n2)  ! kg/m^2 s
          V%tauxo       (i,j)= sbuff(kbufs_tauxo,n2)  ! N/m^2
          V%tauyo       (i,j)= sbuff(kbufs_tauyo,n2)  ! N/m^2
          V%dtauo       (i,j)= sbuff(kbufs_dtauo,n2) ! N*s/m^3
          V%Tref        (i,j)= sbuff(kbufs_Tref, n2)  ! K
          V%sabs        (i,j)= sbuff(kbufs_swabs,n2)  ! W/m^2
          V%Fwe         (i,j)= sbuff(kbufs_we,   n2)  ! m
          if(rec_ice)then
             n1=(j-1)*S%imt+i
             ssbuff(:,n1)=sbuff(:,n2)
          endif
c        i = mod(nint(sbuff(n2,nsnd))-1,S%imt)+1
c        j = (nint(sbuff(n2,nsnd))-1)/S%imt+1
c        Write(nu_diag,*)i,j,S%imt
c          V%ailohi      (i,j)= sbuff(n2,kbufs_ifrc)  ! none
c          V%Tsrf        (i,j)= sbuff(n2,kbufs_tsfc)  ! K
c          V%alvdr       (i,j)= sbuff(n2,kbufs_alb1)  ! none
c          V%alidr       (i,j)= sbuff(n2,kbufs_alb2)  ! none
c          V%alvdf       (i,j)= sbuff(n2,kbufs_alb3)  ! none
c          V%alidf       (i,j)= sbuff(n2,kbufs_alb4)  ! none
c          V%Flatent     (i,j)= sbuff(n2,kbufs_lat )  ! W/m^2
c          V%Fsensible   (i,j)= sbuff(n2,kbufs_sens)  ! W/m^2
c          V%Flwout      (i,j)= sbuff(n2,kbufs_lwup)  ! W/m^2
c          V%evap        (i,j)= sbuff(n2,kbufs_evap)  ! kg/m^2 s
c          V%tauxa       (i,j)= sbuff(n2,kbufs_tauxa)  ! N/m^2
c          V%tauya       (i,j)= sbuff(n2,kbufs_tauya)  ! N/m^2
c          V%Fswthru     (i,j)= sbuff(n2,kbufs_netsw)  ! W/m^2
c          V%Fhnet       (i,j)= sbuff(n2,kbufs_melth)  ! W/m^2
c          V%Fresh       (i,j)= sbuff(n2,kbufs_meltw)  ! kg/m^2 s
c          V%Fsalt       (i,j)= sbuff(n2,kbufs_salt )  ! kg/m^2 s
c          V%tauxo       (i,j)= sbuff(n2,kbufs_tauxo)  ! N/m^2
c          V%tauyo       (i,j)= sbuff(n2,kbufs_tauyo)  ! N/m^2
c          V%Tref        (i,j)= sbuff(n2,kbufs_Tref)  ! K
c          V%sabs        (i,j)= sbuff(n2,kbufs_swabs)  ! W/m^2
      enddo

      if(rec_ice) call write_ice_record(ssbuff,nsnd*S%imt*S%jmt, ice_recno)
      
      else

      call read_ice_record(ssbuff,nsnd*S%imt*S%jmt, ice_recno)

      do n2=1,nsend
        i = mod(nint(sbuff(nsnd,n2))-1,S%imt)+1
        j = (nint(sbuff(nsnd,n2))-1)/S%imt+1
        n1=(j-1)*S%imt+i
          V%ailohi      (i,j)= ssbuff(kbufs_ifrc, n1)  ! none
          V%Tsrf        (i,j)= ssbuff(kbufs_tsfc, n1)  ! K
          V%alvdr       (i,j)= ssbuff(kbufs_alb1, n1)  ! none
          V%alidr       (i,j)= ssbuff(kbufs_alb2, n1)  ! none
          V%alvdf       (i,j)= ssbuff(kbufs_alb3, n1)  ! none
          V%alidf       (i,j)= ssbuff(kbufs_alb4, n1)  ! none
          V%Flatent     (i,j)= ssbuff(kbufs_lat , n1)  ! W/m^2
          V%Fsensible   (i,j)= ssbuff(kbufs_sens, n1)  ! W/m^2
          V%Flwout      (i,j)= ssbuff(kbufs_lwup, n1)  ! W/m^2
          V%evap        (i,j)= ssbuff(kbufs_evap, n1)  ! kg/m^2 s
          V%tauxa       (i,j)= ssbuff(kbufs_tauxa,n1)  ! N/m^2
          V%tauya       (i,j)= ssbuff(kbufs_tauya,n1)  ! N/m^2
          V%Fswthru     (i,j)= ssbuff(kbufs_netsw,n1)  ! W/m^2
          V%Fhnet       (i,j)= ssbuff(kbufs_melth,n1)  ! W/m^2
          V%Fresh       (i,j)= ssbuff(kbufs_meltw,n1)  ! kg/m^2 s
          V%Fsalt       (i,j)= ssbuff(kbufs_salt ,n1)  ! kg/m^2 s
          V%tauxo       (i,j)= ssbuff(kbufs_tauxo,n1)  ! N/m^2
          V%tauyo       (i,j)= ssbuff(kbufs_tauyo,n1)  ! N/m^2
          V%dtauo       (i,j)= ssbuff(kbufs_dtauo,n1) ! N*s/m^3
          V%Tref        (i,j)= ssbuff(kbufs_Tref, n1)  ! K
          V%sabs        (i,j)= ssbuff(kbufs_swabs,n1)  ! W/m^2
          V%Fwe         (i,j)= ssbuff(kbufs_we,   n2)  ! m

      enddo

      endif

      write(nu_diag,*) 'Just received from ice'
      i=168
      j=299
c      write(nu_diag,*) 'atm states and fluxes'
c      write(nu_diag,*) '            frac  = ',V%ailohi      (i,j)
c      write(nu_diag,*) '            Tsrf  = ',V%Tsrf        (i,j)
c      write(nu_diag,*) '            avdr  = ',V%alvdr       (i,j)
c      write(nu_diag,*) '            aidr  = ',V%alidr       (i,j)
c      write(nu_diag,*) '            avdf  = ',V%alvdf       (i,j)
c      write(nu_diag,*) '            aidf  = ',V%alidf       (i,j)
c      write(nu_diag,*) '            flat  = ',V%Flatent     (i,j)
c      write(nu_diag,*) '            fsens = ',V%Fsensible   (i,j)
c      write(nu_diag,*) '            flw   = ',V%Flwout      (i,j)
c      write(nu_diag,*) '            evap  = ',V%evap        (i,j)
c      write(nu_diag,*) '            taux  = ',V%tauxa       (i,j)
c      write(nu_diag,*) '            tauy  = ',V%tauya       (i,j)
      write(nu_diag,*) '            swthru= ',V%Fswthru     (i,j)
      write(nu_diag,*) '            Fhnet = ',V%Fhnet       (i,j)
c      write(nu_diag,*) '            Fresh = ',V%Fresh       (i,j)
c      write(nu_diag,*) '            fsalt = ',V%Fsalt       (i,j)
c      write(nu_diag,*) '            taxo  = ',V%tauxo       (i,j)
c      write(nu_diag,*) '            tayo  = ',V%tauyo       (i,j)
c      write(nu_diag,*) '            dtau  = ',V%dtauo       (i,j)
c      write(nu_diag,*) '            Tref  = ',V%Tref        (i,j)
c      write(nu_diag,*) '            sabs  = ',V%sabs        (i,j)
c      write(nu_diag,*) '            we    = ',V%we        (i,j)

      deallocate(sbuff)
      if(rec_ice .or. rd_rec_ice) deallocate(ssbuff)

      return
	end

      subroutine send_ice_data(ToI, Igr, Ist)

	use ice_indexes
	use struct
        use msgtypes

	Type (state) Ist
	Type (grid) Igr
	Type (to_ice) ToI

      integer cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs(20)
      common/allids/cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs
      integer ::
     &  cpl_task
     &, my_task           ! master task for ice
     &, MPI_COMM_MY
      common/tasks/my_task, cpl_task, MPI_COMM_MY
      integer nu_diag
      common/iodiag/nu_diag
      integer*4, parameter ::
     &   nbuff = 100
	integer*4, dimension(nbuff) ::
     &   buff
	integer*4 nrecv, nrs, nrn
	real*8, allocatable :: sbuff(:,:)

      if(abs(Ist%next+99).lt.1) then
         Write(nu_diag,*)'Sending to ice, but stop signal is already sent'
         return
      endif

c      Write(nu_diag,*)'Ice data are about to be sent',Ist%imt,Ist%jmt

      buff=0
      imt=Ist%imt
      jmt=Ist%jmt
        nrecv=imt*jmt
	  nrs=0
	  nrn=0
        buff(13) = nrecv
        buff(14) = nrs
        buff(15) = nrn

      if(Ist%next.lt.0) buff(3)=1

      call SEND_I(buff, nbuff,msgtype_c2i, ice_id, ierr)

      if(buff(3).eq.1) then
         Ist%next=-99
         Write(nu_diag,*)'Stop signal to ice!!!!!!!!'
c        return
      endif

      allocate(sbuff(nrecv,nrcv))

      do j = 1,jmt
      do i = 1,imt
          i2 = (j-1)*imt + i
        !  atm states                                   
          sbuff(i2,kbufr_zlvl)   = ToI%zlvl (i,j)       ! m
          sbuff(i2,kbufr_uatm)   = ToI%uatm (i,j)       ! m/s
          sbuff(i2,kbufr_vatm)   = ToI%vatm (i,j)       ! m/s
          sbuff(i2,kbufr_potT)   = ToI%potT (i,j)       ! K
          sbuff(i2,kbufr_Tair)   = ToI%Tair (i,j)       ! K
          sbuff(i2,kbufr_Qa)     = ToI%Qa   (i,j)       ! kg/kg
          sbuff(i2,kbufr_rhoa)   = ToI%rhoa (i,j)       ! kg/m^3
        !  atm -> ice fluxes
          sbuff(i2,kbufr_swdvdr) = ToI%swvdr(i,j)       ! W/m^2
          sbuff(i2,kbufr_swdvdf) = ToI%swvdf(i,j)       ! W/m^2
          sbuff(i2,kbufr_swdidr) = ToI%swidr(i,j)       ! W/m^2
          sbuff(i2,kbufr_swdidf) = ToI%swidf(i,j)       ! W/m^2
          sbuff(i2,kbufr_flw)    = ToI%Flw  (i,j)       ! W/m^2
          sbuff(i2,kbufr_rain)   = ToI%Frain(i,j)       ! kg/m^2/s
          sbuff(i2,kbufr_snow)   = ToI%Fsnow(i,j)       ! kg/m^2/s liquid
        !  ocean states
          sbuff(i2,kbufr_sst)    = ToI%sst  (i,j)       ! K
          sbuff(i2,kbufr_sss)    = ToI%sss  (i,j)       ! psu
          sbuff(i2,kbufr_uocn)   = ToI%uocn (i,j)       ! m/s
          sbuff(i2,kbufr_vocn)   = ToI%vocn (i,j)       ! m/s
          sbuff(i2,kbufr_tiltx)  = ToI%strtltx(i,j)     ! m/m
          sbuff(i2,kbufr_tilty)  = ToI%strtlty(i,j)     ! m/m
        !  ocean -> ice flux
          sbuff(i2,kbufr_fw)     = ToI%frzmlt(i,j)      ! W/m^2
      enddo
      enddo
      write(nu_diag,*) 'About to send to ice'
      i=168
      j=299
c      write(nu_diag,*) 'atm states and fluxes'
c      write(nu_diag,*) '            uatm  = ',minval(ToI%uatm,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%uatm,MASK=(Igr%mask.ne.0))
c      write(nu_diag,*) '            vatm  = ',minval(ToI%vatm,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%vatm,MASK=(Igr%mask.ne.0))
c      write(nu_diag,*) '            potT  = ',minval(ToI%potT,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%potT,MASK=(Igr%mask.ne.0))
c      write(nu_diag,*) '            Tair  = ',minval(ToI%Tair,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%Tair,MASK=(Igr%mask.ne.0))
c      write(nu_diag,*) '            Qa    = ',minval(ToI%Qa,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%Qa,MASK=(Igr%mask.ne.0))
c      write(nu_diag,*) '            rhoa  = ',minval(ToI%rhoa,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%rhoa,MASK=(Igr%mask.ne.0))
c      write(nu_diag,*) '            sw    = ',minval(ToI%swvdr/0.28,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%swvdr/0.28,MASK=(Igr%mask.ne.0))
      write(nu_diag,*) '            sw = ',ToI%swvdr(i,j)+ToI%swvdf(i,j)+ToI%swidr(i,j)+ToI%swidf(i,j)
c      write(nu_diag,*) '            lw    = ',minval(ToI%Flw,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%Flw,MASK=(Igr%mask.ne.0))
c      write(nu_diag,*) '            frain = ',minval(ToI%Frain,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%Frain,MASK=(Igr%mask.ne.0))
c      write(nu_diag,*) '            fsnow = ',minval(ToI%Fsnow,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%Fsnow,MASK=(Igr%mask.ne.0))
c      write(nu_diag,*) ' '
c      write(nu_diag,*) 'ocn states and fluxes'
c      write(nu_diag,*) '            frzmlt  = ',minval(ToI%frzmlt,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%frzmlt,MASK=(Igr%mask.ne.0))
c      write(nu_diag,*) '            sst     = ',minval(ToI%sst,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%sst,MASK=(Igr%mask.ne.0))
c      write(nu_diag,*) '            sss     = ',minval(ToI%sss,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%sss,MASK=(Igr%mask.ne.0))
c      write(nu_diag,*) '            uocn    = ',minval(ToI%uocn,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%uocn,MASK=(Igr%mask.ne.0))
c      write(nu_diag,*) '            vocn    = ',minval(ToI%vocn,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%vocn,MASK=(Igr%mask.ne.0))
c      write(nu_diag,*) '            strtltx = ',minval(ToI%strtltx,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%strtltx,MASK=(Igr%mask.ne.0))
c      write(nu_diag,*) '            strtlty = ',minval(ToI%strtlty,MASK=(Igr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToI%strtlty,MASK=(Igr%mask.ne.0))

      call SEND_D(sbuff, nrecv*nrcv, msgtype_c2i, ice_id, ierr)

      deallocate(sbuff)
      
c      Write(nu_diag,*)'Ice data are sent successfully'
      
      return
	end

      subroutine send_ocean_init_to_ice(ToI,
     &               Ice, Igr, Ist,
     &               Ocn, Ogr, Ost)

	use ice_indexes
	use struct
	use conversion
        use msgtypes

	Type (state) Ist,Ost
	Type (grid) Igr,Ogr
	Type (from_ice) Ice
	Type (from_ocn) Ocn
	Type (to_ice) ToI

      integer cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs(20)
      common/allids/cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs
      integer ::
     &  cpl_task
     &, my_task           ! master task for ice
     &, MPI_COMM_MY
      common/tasks/my_task, cpl_task, MPI_COMM_MY
      integer nu_diag
      common/iodiag/nu_diag

      Write(nu_diag,*)'Ocean initial data are about to be sent to ice',Ist%imt,Ist%jmt

      write(nu_diag,*)'i2o',maxval(Ocn%sst),minval(Ocn%sst)
      call convert_scalar(o2i, Ocn%sst,    Ogr, Ost, ToI%sst,     Igr, Ist)
      call convert_scalar(o2i, Ocn%sss,    Ogr, Ost, ToI%sss,     Igr, Ist)
      write(nu_diag,*)'i2ov'
      call convert_vector(o2iv, Ocn%uocn, Ocn%vocn, Ogr, Ost, 
     *   ToI%uocn, ToI%vocn, Igr, Ist)

      write(nu_diag,*)'i2o',maxval(ToI%sst),minval(ToI%sst)
      call SEND_D(ToI%sst, Ist%imt*Ist%jmt, msgtype_c2i, ice_id, ierr)
      call SEND_D(ToI%sss, Ist%imt*Ist%jmt, msgtype_c2i, ice_id, ierr)
      call SEND_D(ToI%uocn, Ist%imt*Ist%jmt, msgtype_c2i, ice_id, ierr)
      call SEND_D(ToI%vocn, Ist%imt*Ist%jmt, msgtype_c2i, ice_id, ierr)
     
      call convert_scalar(o2i, Ogr%work,    Ogr, Ost, Igr%work,     Igr, Ist)
      call SEND_D(Igr%work, Ist%imt*Ist%jmt, msgtype_c2i, ice_id, ierr)

      Write(nu_diag,*)'Ocean initial data are sent successfully to ice'
      
      return
	end

      subroutine get_ocean_state(V, S, nrivers)

	use ocn_indexes
	use struct
        use msgtypes

	Type (state) S
	Type (from_ocn) V

      integer cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs(20)
      common/allids/cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs
      integer ::
     &  cpl_task
     &, my_task           ! master task for ice
     &, MPI_COMM_MY
      common/tasks/my_task, cpl_task, MPI_COMM_MY
      integer nu_diag
      common/iodiag/nu_diag
      integer*4, parameter ::
     &   nbuff = 100
	integer*4, dimension(nbuff) ::
     &   buff
      integer*4 ierr,stop_now,nsend,ocn_recno
      common/stop_flag/stop_now
	real*8, allocatable :: sbuff(:,:)
        logical rec_ocn,rec_atm,rec_ice,rec_lnd,
     &     rd_rec_ocn,rd_rec_atm,rd_rec_ice,rd_rec_lnd
        common/recstate/rec_ocn,rec_atm,rec_ice,rec_lnd,
     &     rd_rec_ocn,rd_rec_atm,rd_rec_ice,rd_rec_lnd

c      if(S%next.lt.0) return

c      Write(nu_diag,*) 'Receive ocean model configuration'
	call recv_i(buff, nbuff, msgtype_o2c, ocn_id, ierr)
      if(buff(3).eq.1) stop_now=1   !stop now flag
      S%idate = buff( 4)        ! date (coded: yyyymmdd)
      S%sec   = buff( 5)        ! elapsed seconds into date
c      S%nadv  = buff( 9)        ! ice comm pairs per day
      nsend   = buff(10)        ! size of real buffer to send
      nrivers = buff(11)
c	Write(nu_diag,'(a29,4(/,a29,i8))')
c     * 'Received ocean information:    ',
c     * '            date (yyyymmdd): ', S%idate,
c     * '            elapsed seconds: ', S%sec,
c     * '          ocn steps per day: ', S%nadv,
c     * 'size of real buffer to send: ', nsend,
c     * '                   stop now: ', stop_now
       Write(nu_diag,*)'received from ocn nrivers = ',nrivers
	
      if(rec_ocn .or. rd_rec_ocn)then
         ocn_recno=dayofyear(S%idate)
         ocn_step=S%sec*S%nadv/86400
         ocn_recno=(ocn_recno-1)*S%nadv+ocn_step+1
         Write(nu_diag,*) 'Receive ocn model data for record ',ocn_recno
      endif

c      if(stop_now == 1) call exit_coupler
      allocate(sbuff(nsnd,nsend))

c      Write(nu_diag,*) 'Receive ocean model data'
c      Write(nu_diag,*) '               dimension: ',S%imt,S%jmt
	call recv_d(sbuff, nsnd*nsend, msgtype_o2c, ocn_id, ierr)
c      Write(nu_diag,*) 'Received'

      do j=1,S%jmt
      do i=1,S%imt
         V%sst     (i,j) = 0.0
         V%sss     (i,j) = 0.0
         V%uocn    (i,j) = 0.0
         V%vocn    (i,j) = 0.0
         V%tiltx   (i,j) = 0.0
         V%tilty   (i,j) = 0.0
         V%frzmlt  (i,j) = 0.0
         V%cco2    (i,j) = 0.0
      enddo
      enddo
c      Write(nu_diag,*) 'Initialized'

      if(.not.rd_rec_ocn)then

         if(rec_ocn) call write_ocn_record(sbuff,nsnd*nsend, ocn_recno)

      else
      
         call read_ocn_record(sbuff,nsnd*nsend, ocn_recno)
         
      endif

      do n2=1,nsend
        i = mod(nint(sbuff(nsnd,n2))-1,S%imt)+1
        j = (nint(sbuff(nsnd,n2))-1)/S%imt+1
c      Write(nu_diag,*) '     No',n2,':',sbuff(nsnd,n2),'--(',i,j,')'
         V%sst    (i,j) = sbuff(kbufr_sst,    n2)  ! C
         V%sss    (i,j) = sbuff(kbufr_sss,    n2)  ! psu
         V%uocn   (i,j) = sbuff(kbufr_uocn,   n2)  ! m/s
         V%vocn   (i,j) = sbuff(kbufr_vocn,   n2)  ! m/s
         V%tiltx  (i,j) = sbuff(kbufr_tiltx,  n2)  ! m/m
         V%tilty  (i,j) = sbuff(kbufr_tilty,  n2)  ! m/m
         V%frzmlt (i,j) = sbuff(kbufr_frzmlt, n2)  ! W/m^2
         V%cco2   (i,j) = sbuff(kbufr_co2,    n2)  !
c      Write(nu_diag,*) '      ... done!'
      enddo
c      write(nu_diag,*) 'Just received from ocn'
c      i=60
c      j=2
c      write(nu_diag,*) '            SST   = ',V%sst      (i,j)
c      write(nu_diag,*) '            SSS   = ',V%sss      (i,j)
c      write(nu_diag,*) '            U     = ',V%uocn     (i,j)
c      write(nu_diag,*) '            V     = ',V%vocn     (i,j)
c      write(nu_diag,*) '            TiltX = ',V%tiltx    (i,j)
c      write(nu_diag,*) '            TiltY = ',V%tilty    (i,j)
c      write(nu_diag,*) '            Frz   = ',V%frzmlt   (i,j)
c      write(nu_diag,*) V%sst(80,:)

c      write(nu_diag,*) '            SST   = ',minval(V%sst),' ',maxval(V%sst)
c      write(nu_diag,*) '            SSS   = ',minval(V%sss),' ',maxval(V%sss)
c      write(nu_diag,*) '            U     = ',minval(V%uocn),' ',maxval(V%uocn)
c      write(nu_diag,*) '            V     = ',minval(V%vocn),' ',maxval(V%vocn)
c      write(nu_diag,*) '            TiltX = ',minval(V%tiltx),' ',maxval(V%tiltx)
c      write(nu_diag,*) '            TiltY = ',minval(V%tilty),' ',maxval(V%tilty)
c      write(nu_diag,*) '            Frz   = ',minval(V%frzmlt),' ',maxval(V%frzmlt)
c      write(nu_diag,*) '            CO2   = ',minval(V%cco2),' ',maxval(V%cco2)

      deallocate(sbuff)

      return
	end

      subroutine send_ocean_data(current, ToO, Ogr, Ost)

	use ocn_indexes
	use struct
        use msgtypes

	Type (state) Ost
	Type (grid) Ogr
	Type (to_ocn) ToO

      integer cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs(20)
      common/allids/cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs
      integer ::
     &  cpl_task
     &, my_task           ! master task for ice
     &, MPI_COMM_MY
      common/tasks/my_task, cpl_task, MPI_COMM_MY
      integer nu_diag, nold
      character*80 fname
      save nold
      common/iodiag/nu_diag
      integer*4, parameter ::
     &   nbuff = 100
	integer*4, dimension(nbuff) ::
     &   buff
	integer*4 nrecv
	real*8, allocatable :: sbuff(:,:)
      data nold/-1000/

      if(abs(Ost%next+99).lt.1) then
         Write(nu_diag,*)'Sending to ocn, but stop signal is already sent'
         return
      endif

      rn=current/86400.0-int(current/86400.0/365.0)*365.0
      n=int(rn*24.0*60.0)
      Write(nu_diag,*)'Minutes: ',n

      imt=Ost%imt
	jmt=Ost%jmt
        nrecv=imt*jmt
        buff(13) = nrecv
        buff(14) = n

      if(Ost%next.lt.0) buff(3)=1

      call SEND_I(buff, nbuff,msgtype_c2o, ocn_id, ierr)

      if(buff(3).eq.1) then
         Ost%next=-99
         Write(nu_diag,*)'Stop signal to ocn!!!!!!!!'
c         return
      endif

      allocate(sbuff(nrecv,nrcv))

      do j = 1,jmt
      do i = 1,imt
          i2 = (j-1)*imt + i
          sbuff(i2,kbufs_swave)  = ToO%Fswthru (i,j)
          sbuff(i2,kbufs_netheat)= ToO%Fhnet (i,j)
          sbuff(i2,kbufs_prec)   = ToO%Fprec (i,j)
          sbuff(i2,kbufs_evap)   = ToO%Fevap (i,j)
          sbuff(i2,kbufs_melt)   = ToO%Fmelt (i,j)
c          sbuff(i2,kbufs_salt)   = ToO%Fsalt (i,j)
          sbuff(i2,kbufs_taux)   = ToO%tauxo (i,j)
          sbuff(i2,kbufs_tauy)   = ToO%tauyo (i,j)
          sbuff(i2,kbufs_dtau)   = ToO%dtauo (i,j)
          sbuff(i2,kbufs_iswave)  = ToO%Fswthru_ice (i,j)
          sbuff(i2,kbufs_inetheat)= ToO%Fhnet_ice (i,j)
          sbuff(i2,kbufs_itaux)   = ToO%tauxoi (i,j)
          sbuff(i2,kbufs_itauy)   = ToO%tauyoi (i,j)
          sbuff(i2,kbufs_idtau)   = ToO%dtauoi (i,j)
          sbuff(i2,kbufs_ifrc)   = ToO%Frac (i,j)
          sbuff(i2,kbufs_dnetheat) = ToO%dFhnet(i,j)
          sbuff(i2,kbufs_utidal)   = ToO%utidal (i,j)
          sbuff(i2,kbufs_vtidal)   = ToO%vtidal (i,j)
          sbuff(i2,kbufs_sriv)   = ToO%Frivprec (i,j)
          sbuff(i2,kbufs_co2)   = ToO%Fco2 (i,j)
          sbuff(i2,kbufs_we)   = ToO%Fwe (i,j)
          sbuff(i2,kbufs_index)  = i*1000+j
      enddo
      enddo

c      write(nu_diag,*) 'About to send to ocn'
cc      i=60
cc      j=2
c      write(nu_diag,*) '            SW rad  = ',minval(ToO%Fswthru,MASK=(Ogr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToO%Fswthru,MASK=(Ogr%mask.ne.0))
c      write(nu_diag,*) '            Net Ht  = ',minval(ToO%Fhnet,MASK=(Ogr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToO%Fhnet,MASK=(Ogr%mask.ne.0))
c      write(nu_diag,*) '            Fprec   = ',minval(ToO%Fprec,MASK=(Ogr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToO%Fprec,MASK=(Ogr%mask.ne.0))
c      write(nu_diag,*) '            Fevap   = ',minval(ToO%Fevap,MASK=(Ogr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToO%Fevap,MASK=(Ogr%mask.ne.0))
c      write(nu_diag,*) '            Fmelt   = ',minval(ToO%Fmelt,MASK=(Ogr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToO%Fmelt,MASK=(Ogr%mask.ne.0))
c      write(nu_diag,*) '            TauX    = ',minval(ToO%tauxo,MASK=(Ogr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToO%tauxo,MASK=(Ogr%mask.ne.0))
c      write(nu_diag,*) '            TauY    = ',minval(ToO%tauyo,MASK=(Ogr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToO%tauyo,MASK=(Ogr%mask.ne.0))
c      write(nu_diag,*) '            dTau    = ',minval(ToO%dtauo,MASK=(Ogr%mask.ne.0)),' <---> ',
c     *                                        maxval(ToO%dtauo,MASK=(Ogr%mask.ne.0))
cc      write(nu_diag,*) ToO%Fhnet(80,:)

      call SEND_D(sbuff, nrecv*nrcv, msgtype_c2o, ocn_id, ierr)

      deallocate(sbuff)

c      call cdf2Dvar('Fswthru',ToO%Fswthru,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('Fhnet',ToO%Fhnet,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('Fprec',ToO%Fprec,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('Fevap',ToO%Fevap,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('Fmelt',ToO%Fmelt,Ogr%mask,imt,jmt,.true.,Ost%idate)
cc      call cdf2Dvar('Fsalt',ToO%Fsalt,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('tauxo',ToO%tauxo,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('tauyo',ToO%tauyo,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('dtauo',ToO%dtauo,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('Fswthru_ice',ToO%Fswthru_ice,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('Fhnet_ice',ToO%Fhnet_ice,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('tauxoi',ToO%tauxoi,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('tauyoi',ToO%tauyoi,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('dtauoi',ToO%dtauoi,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('Frac',ToO%Frac,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('dFhnet',ToO%dFhnet,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('utidal',ToO%utidal,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      call cdf2Dvar('vtidal',ToO%vtidal,Ogr%mask,imt,jmt,.true.,Ost%idate)
c      if(nold.ne.n.and.nold.ne.-1000) stop
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c        Write data tranfered to the ocean
      n=int(current/86400.0-int(current/86400.0/365.0)*365.0)+1
      if(n.gt.365)n=365

      if(n.ne.nold)then
      nold=n

      if(.false.)then
      fname='../Output/cpl/ocnDDD.bin'
      ii=index(fname,'DDD');
      write(fname(ii:ii+2),'(i3.3)') n

      open(23,file=fname,form='unformatted',access='direct',recl=8*imt*jmt)
      Write(23,rec=1) ToO%Fswthru
      Write(23,rec=2) ToO%Fhnet
      Write(23,rec=3) ToO%Fprec
      Write(23,rec=4) ToO%Fevap
      Write(23,rec=5) ToO%tauxo
      Write(23,rec=6) ToO%tauyo
      Write(23,rec=7) ToO%dtauo
      Write(23,rec=8) ToO%Fmelt
c      Print *, 'cpl ',ToO%tauxo(78,159)
      close(23)
      endif
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      endif

      return
	end

      Subroutine get_grid(G, S, key, id, ierr)
	
	use struct
	Type (grid) G
	Type (state) S
	Integer key, id
	Integer*4 ierr
      integer nu_diag
      common/iodiag/nu_diag

	call recv_d(G%lon, S%imt*S%jmt, key, id, ierr)
	Write(nu_diag,*) 'Grid longitudes are received'
c        Write(nu_diag,'(i3,i3,f12.3,/)') ((i,j,G%lon(i,j),i=1,S%imt),j=1,S%jmt)
	call recv_d(G%lat, S%imt*S%jmt, key, id, ierr)
	Write(nu_diag,*) 'Grid latitudes are received'
        Write(nu_diag,*) (G%lat(S%imt/2,j),j=1,S%jmt)
c	call recv_d(G%work, 4*S%imt*S%jmt, key, id, ierr)
c	Write(nu_diag,*) 'Grid array 1 is received'
c	call recv_d(G%work, 4*S%imt*S%jmt, key, id, ierr)
c	Write(nu_diag,*) 'Grid array 2 is received'
	call recv_d(G%tarea, S%imt*S%jmt, key, id, ierr)
	Write(nu_diag,*) 'Grid cell area array is received'
	call recv_i(G%mask, S%imt*S%jmt, key, id, ierr)
	Write(nu_diag,*) 'Grid mask is received'
	call recv_d(G%work, S%imt*S%jmt, key, id, ierr)
c	Write(nu_diag,*) 'Grid array 1 is received'

        do i=1,S%imt
           do j=1,S%jmt
              if(G%mask(i,j).ne.0) G%mask(i,j)=1
           enddo
        enddo

      Return
	End

      subroutine get_atm_state(V, S)

      use struct
      use atm_indexes
      use msgtypes

      Type (state) S
      Type (from_atm) V

      integer cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs(20)
      common/allids/cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs
      integer ::
     &  cpl_task
     &, my_task           ! master task for ice
     &, MPI_COMM_MY
      common/tasks/my_task, cpl_task, MPI_COMM_MY
      integer nu_diag
      common/iodiag/nu_diag
      integer*4, parameter ::
     &   nbuff = 100
	integer*4, dimension(nbuff) ::
     &   buff
      integer*4 ierr,stop_now,nsend,cunt,atm_recno
      common/stop_flag/stop_now
	real*8, allocatable :: sbuff(:,:)
        integer, allocatable :: msk(:,:)
        logical rec_ocn,rec_atm,rec_ice,rec_lnd,
     &     rd_rec_ocn,rd_rec_atm,rd_rec_ice,rd_rec_lnd
        common/recstate/rec_ocn,rec_atm,rec_ice,rec_lnd,
     &     rd_rec_ocn,rd_rec_atm,rd_rec_ice,rd_rec_lnd
      logical init
      data cunt/5/,init/.true./

c      if(S%next.lt.0) return

c      Write(nu_diag,*) 'Waiting atm model configuration'
	call recv_i(buff, nbuff, msgtype_a2c, atm_id, ierr)
      if(buff(3).eq.1) stop_now=1   !stop now flag
      S%idate = buff(4)          ! date (coded: yyyymmdd)
      S%sec   = buff(5)            ! elapsed seconds into date
c      S%nadv  = buff(9)          ! atm comm pairs per day
      nsend   = buff(10)          ! size of real buffer to send
c	Write(nu_diag,'(a29,4(/,a29,i8))')
c     * 'Received atm information:    ',
c     * '            date (yyyymmdd): ', S%idate,
c     * '            elapsed seconds: ', S%sec,
c     * 'size of real buffer to send: ', nsend
c	if(stop_now == 1) call exit_coupler
c      allocate(sbuff(nsend,nsnd))

      if(rec_atm .or. rd_rec_atm)then
         atm_recno=dayofyear(S%idate)
         atm_step=S%sec*S%nadv/86400
         atm_recno=(atm_recno-1)*S%nadv+atm_step+1
         Write(nu_diag,*) 'Receive atm model data for record ',atm_recno
      endif

      allocate(sbuff(nsnd,nsend))

c      Write(nu_diag,*) 'Receive atm model data'
	call recv_d(sbuff, nsnd*nsend, msgtype_a2c, atm_id, ierr)

      do j=1,S%jmt
      do i=1,S%imt
         V%zlvl        (i,j) = 0.0
         V%ztlvl       (i,j) = 0.0
         V%uatm        (i,j) = 0.0
         V%vatm        (i,j) = 0.0
         V%potT        (i,j) = 0.0
         V%Tair        (i,j) = 0.0
         V%Qa          (i,j) = 0.0
         V%rhoa        (i,j) = 0.0
         V%Fsw         (i,j) = 0.0
         V%Flw         (i,j) = 0.0
         V%Frain       (i,j) = 0.0
         V%Fsnow       (i,j) = 0.0
         V%cld         (i,j) = 0.0
         V%small_riv   (i,j) = 0.0
         V%cco2        (i,j) = 0.0
      enddo
      enddo
      
      if(.not.rd_rec_atm)then

         if(rec_atm) call write_atm_record(sbuff,nsnd*nsend, atm_recno)
      
      else
      
         call read_atm_record(sbuff,nsnd*nsend, atm_recno)
         
      endif

      do n2=1,nsend
        i = mod(nint(sbuff(nsnd,n2))-1,S%imt)+1
        j = (nint(sbuff(nsnd,n2))-1)/S%imt+1
          V%zlvl        (i,j)= sbuff(kbufs_zlvl,  n2)  ! m
          V%ztlvl       (i,j)= sbuff(kbufs_ztlvl, n2)  ! m
          V%uatm        (i,j)= sbuff(kbufs_uatm,  n2)  ! m/s
          V%vatm        (i,j)= sbuff(kbufs_vatm,  n2)  ! m/s
          V%potT        (i,j)= sbuff(kbufs_potT,  n2)  ! K
          V%Tair        (i,j)= sbuff(kbufs_Tair,  n2)  ! K
          V%Qa          (i,j)= sbuff(kbufs_Qa,    n2)  ! kg/kg
          V%rhoa        (i,j)= sbuff(kbufs_rhoa,  n2)  ! kg/m^3
          V%Fsw         (i,j)= sbuff(kbufs_Fsw,   n2)  ! W/m^2
          V%Flw         (i,j)= sbuff(kbufs_Flw,   n2)  ! W/m^2
          V%Frain       (i,j)= sbuff(kbufs_Frain, n2)  ! kg/m^2 s
          V%Fsnow       (i,j)= sbuff(kbufs_Fsnow, n2)  ! kg/m^2 s
          V%cld         (i,j)= sbuff(kbufs_cld,   n2)  ! none
          V%small_riv   (i,j)= sbuff(kbufs_small_riv,n2)  ! kg/m^2 s
          V%cco2        (i,j)= sbuff(kbufs_co2,   n2)  ! 
      enddo

c      allocate(msk(S%imt,S%jmt))
c      msk=1
c
c      call cdf2Dvar('zlvl',V%zlvl,msk,S%imt,S%jmt,init,S%idate)
c      call cdf2Dvar('uatm',V%uatm,msk,S%imt,S%jmt,init,S%idate)
c      call cdf2Dvar('vatm',V%vatm,msk,S%imt,S%jmt,init,S%idate)
c      call cdf2Dvar('potT',V%potT,msk,S%imt,S%jmt,init,S%idate)
c      call cdf2Dvar('Tair',V%Tair,msk,S%imt,S%jmt,init,S%idate)
c      call cdf2Dvar('Qa',V%Qa,msk,S%imt,S%jmt,init,S%idate)
c      call cdf2Dvar('rhoa',V%rhoa,msk,S%imt,S%jmt,init,S%idate)
c      call cdf2Dvar('Fsw',V%Fsw,msk,S%imt,S%jmt,init,S%idate)
c      call cdf2Dvar('Flw',V%Flw,msk,S%imt,S%jmt,init,S%idate)
c      call cdf2Dvar('Frain',V%Frain,msk,S%imt,S%jmt,init,S%idate)
c      call cdf2Dvar('Fsnow',V%Fsnow,msk,S%imt,S%jmt,init,S%idate)
c      call cdf2Dvar('cld',V%cld,msk,S%imt,S%jmt,init,S%idate)
c      call cdf2Dvar('runoff',V%small_riv,msk,S%imt,S%jmt,init,S%idate)
      init=.false.
cc      if(maxval(V%Flw).eq.0.0) calculate Flw
c      if(cunt.le.0)stop
c      cunt=cunt-1
c
c      deallocate(msk)

      deallocate(sbuff)

      return
      end

      subroutine send_atmos_data(ToA, Agr, Ast)

	use atm_indexes
	use struct
        use msgtypes

	Type (state) Ast
	Type (grid) Agr
	Type (to_atm) ToA

      integer cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs(20)
      common/allids/cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs
      integer ::
     &  cpl_task
     &, my_task           ! master task for ice
     &, MPI_COMM_MY
      common/tasks/my_task, cpl_task, MPI_COMM_MY
      integer nu_diag
      common/iodiag/nu_diag
      integer*4, parameter ::
     &   nbuff = 100
	integer*4, dimension(nbuff) ::
     &   buff
	integer*4 nrecv, nrs, nrn
	real*8, allocatable :: sbuff(:,:)
        integer, allocatable :: msk(:,:)
      logical init
      data init/.true./

      if(abs(Ast%next+99).lt.1) then
         Write(nu_diag,*)'Sending to atm, but stop signal is already sent'
         return
      endif

c      Write(nu_diag,*)'Atm data are about to be sent',Ast%imt,Ast%jmt

      buff=0
      imt=Ast%imt
      jmt=Ast%jmt
      nrecv=imt*jmt
	  nrs=0
	  nrn=0
      buff(13) = nrecv
      buff(14) = nrs
      buff(15) = nrn

      if(Ast%next.lt.0) buff(3)=1

      call SEND_I(buff, nbuff,msgtype_c2a, atm_id, ierr)
c      Write(nu_diag,*)'Integer buffer is sent'

      if(buff(3).eq.1) then
         Ast%next=-99
         Write(nu_diag,*)'Stop signal to atm!!!!!!!!'
c         return
      endif

      allocate(sbuff(nrecv,nrcv))

      do j = 1,jmt
      do i = 1,imt
          i2 = (j-1)*imt + i
          sbuff(i2,kbufr_tsurf)  = ToA%tsurf (i,j)
          sbuff(i2,kbufr_cu)  = ToA%Cu (i,j)
          sbuff(i2,kbufr_ct)  = ToA%Ct (i,j)
          sbuff(i2,kbufr_ce)   = ToA%Ce  (i,j)
          sbuff(i2,kbufr_albed)   = ToA%albed  (i,j)
          sbuff(i2,kbufr_dicec)   = ToA%dicec  (i,j)
          sbuff(i2,kbufr_co2)    = ToA%Fco2  (i,j)
          sbuff(i2,kbufr_index)  = i*1000+j
      enddo
      enddo

      allocate(msk(imt,jmt))
      msk=1
c      Write(nu_diag,*)'Data to atmosphere'
c      Write(nu_diag,*)imt,jmt,Ast%idate
c      call cdf2Dvar('tsurf',ToA%tsurf,msk,imt,jmt,init,Ast%idate)
c      Write(nu_diag,*)'tsurf is done'
c      call cdf2Dvar('Cu',ToA%Cu,msk,imt,jmt,init,Ast%idate)
c      Write(nu_diag,*)'Cu is done'
c      call cdf2Dvar('Ct',ToA%Ct,msk,imt,jmt,init,Ast%idate)
c      Write(nu_diag,*)'Ct is done'
c      call cdf2Dvar('Ce',ToA%Ce,msk,imt,jmt,init,Ast%idate)
c      Write(nu_diag,*)'Ce is done'
c      call cdf2Dvar('albedo',ToA%albed,msk,imt,jmt,init,Ast%idate)
c      Write(nu_diag,*)'albed is done'
c      call cdf2Dvar('dicec',ToA%dicec,msk,imt,jmt,init,Ast%idate)
c      Write(nu_diag,*)'dicec is done'
c      call cdf2Dvar('co2',ToA%Fco2,msk,imt,jmt,init,Ast%idate)
      init=.false.
      deallocate(msk)

      call SEND_D(sbuff, nrecv*nrcv, msgtype_c2a, atm_id, ierr)

      deallocate(sbuff)
c      stop
      
c      Write(nu_diag,*)'Atm data are sent successfully'
      
      return
	end

      subroutine get_lnd_state(V, S)

      use struct
      use lnd_indexes
      use msgtypes

      Type (state) S
      Type (from_lnd) V

      integer cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs(20)
      common/allids/cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs

      integer ::
     &  cpl_task
     &, my_task           ! master task for ice
     &, MPI_COMM_MY
      common/tasks/my_task, cpl_task, MPI_COMM_MY

      integer nu_diag
      common/iodiag/nu_diag

      integer*4, parameter ::
     &   nbuff = 100
	integer*4, dimension(nbuff) ::
     &   buff

      integer*4 ierr,stop_now,nsend,cunt,lnd_recno
      common/stop_flag/stop_now

	real*8, allocatable :: sbuff(:,:)
        integer, allocatable :: msk(:,:)
        logical rec_ocn,rec_atm,rec_ice,rec_lnd,
     &     rd_rec_ocn,rd_rec_atm,rd_rec_ice,rd_rec_lnd
        common/recstate/rec_ocn,rec_atm,rec_ice,rec_lnd,
     &     rd_rec_ocn,rd_rec_atm,rd_rec_ice,rd_rec_lnd
      logical init
      data cunt/5/,init/.true./

c      if(S%next.lt.0) return

      Write(nu_diag,*) 'Waiting lnd model configuration'
	call recv_i(buff, nbuff, msgtype_l2c, lnd_id, ierr)
      if(buff(3).eq.1) stop_now=1   !stop now flag
      S%idate = buff(4)             ! date (coded: yyyymmdd)
      S%sec   = buff(5)             ! elapsed seconds into date
c      S%nadv  = buff(9)            ! lnd comm pairs per day
      nsend   = buff(10)            ! size of real buffer to send
	Write(nu_diag,'(a29,4(/,a29,i8))')
     * 'Received lnd information:    ',
     * '            date (yyyymmdd): ', S%idate,
     * '            elapsed seconds: ', S%sec,
     * 'size of real buffer to send: ', nsend
c	if(stop_now == 1) call exit_coupler
c      allocate(sbuff(nsend,nsnd))

      if(rec_lnd .or. rd_rec_lnd)then
         lnd_recno=dayofyear(S%idate)
         lnd_step=S%sec*S%nadv/86400
         lnd_recno=(lnd_recno-1)*S%nadv+lnd_step+1
         Write(nu_diag,*) 'Receive lnd model data for record ',lnd_recno
      endif

      allocate(sbuff(nsnd,nsend))

      Write(nu_diag,*) 'Receive lnd model data',nsnd,nsend
	call recv_d(sbuff, nsnd*nsend, msgtype_l2c, lnd_id, ierr)
      Write(nu_diag,*) 'Lnd model data are received'

      do j=1,S%jmt
      do i=1,S%imt
         V%NatCO2      (i,j) = 0.0
         V%AntCO2      (i,j) = 0.0
         V%FosCO2      (i,j) = 0.0
      enddo
      enddo
      
      if(.not.rd_rec_lnd)then

         if(rec_lnd) call write_lnd_record(sbuff,nsnd*nsend, lnd_recno)
      
      else
      
         call read_lnd_record(sbuff,nsnd*nsend, lnd_recno)
         
      endif

      do n2=1,nsend
        i = mod(nint(sbuff(nsnd,n2))-1,S%imt)+1
        j = (nint(sbuff(nsnd,n2))-1)/S%imt+1
          V%NatCO2        (i,j)= sbuff(kbufs_NatCO2, n2)  ! kg/m2/s
          V%AntCO2        (i,j)= sbuff(kbufs_AntCO2, n2)  ! kg/m2/s
          V%FosCO2        (i,j)= sbuff(kbufs_FosCO2, n2)  ! kg/m2/s
      enddo

      allocate(msk(S%imt,S%jmt))
      msk=1

c      call cdf2Dvar('NatCO2',V%NatCO2,msk,S%imt,S%jmt,init,S%idate)
c      call cdf2Dvar('AntCO2',V%AntCO2,msk,S%imt,S%jmt,init,S%idate)
c      call cdf2Dvar('FosCO2',V%FosCO2,msk,S%imt,S%jmt,init,S%idate)
      init=.false.
c      if(cunt.le.0)stop
c      cunt=cunt-1
c
      deallocate(msk)

      deallocate(sbuff)

      return
      end

      subroutine send_land_data(ToL, Lgr, Lst)

	use lnd_indexes
	use struct
        use msgtypes

	Type (state) Lst
	Type (grid) Lgr
	Type (to_lnd) ToL

      integer cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs(20)
      common/allids/cpl_id, ice_id, ocn_id, atm_id, lnd_id, regs

      integer ::
     &  cpl_task
     &, my_task           ! master task for ice
     &, MPI_COMM_MY
      common/tasks/my_task, cpl_task, MPI_COMM_MY

      integer nu_diag
      common/iodiag/nu_diag

      integer*4, parameter ::
     &   nbuff = 100
	integer*4, dimension(nbuff) ::
     &   buff

	integer*4 nrecv, nrs, nrn
	real*8, allocatable :: sbuff(:,:)
        integer, allocatable :: msk(:,:)
      logical init
      data init/.true./

      if(abs(Lst%next+99).lt.1) then
         Write(nu_diag,*)'Sending to lnd, but stop signal is already sent'
         return
      endif

c      Write(nu_diag,*)'Lnd data are about to be sent',Lst%imt,Lst%jmt

      buff=0
      imt=Lst%imt
      jmt=Lst%jmt
      nrecv=imt*jmt
	  nrs=0
	  nrn=0
      buff(13) = nrecv
      buff(14) = nrs
      buff(15) = nrn

      if(Lst%next.lt.0) buff(3)=1

      call SEND_I(buff, nbuff, msgtype_c2l, lnd_id, ierr)
c      Write(nu_diag,*)'Integer buffer is sent'

      if(buff(3).eq.1) then
         Lst%next=-99
         Write(nu_diag,*)'Stop signal to lnd!!!!!!!!'
c         return
      endif

      allocate(sbuff(nrecv,nrcv))

      do j = 1,jmt
      do i = 1,imt
          i2 = (j-1)*imt + i
          sbuff(i2,kbufr_Work)  = ToL%Work (i,j)
          sbuff(i2,kbufr_index)  = i*1000+j
      enddo
      enddo

c      allocate(msk(imt,jmt))
c      msk=1
c      Write(nu_diag,*)'Data to land'
c      Write(nu_diag,*)imt,jmt,Lst%idate
c      call cdf2Dvar('Work',ToL%Work,msk,imt,jmt,init,Lst%idate)
      init=.false.
c      deallocate(msk)

      call SEND_D(sbuff, nrecv*nrcv, msgtype_c2l, lnd_id, ierr)

      deallocate(sbuff)
c      stop
      
c      Write(nu_diag,*)'Lnd data are sent successfully'
      
      return
	end

      integer*4 function dayofyear(idate)
      integer*4 idate, monlen(12)
      save monlen
      data monlen/0,31,59,90,120,151,181,212,243,273,304,334/
      n=idate-int(idate/10000)*10000
      m=int(n/100)
      n=n-m*100
      dayofyear=monlen(m)+n
      return
      end

      subroutine write_ice_record(sbuff,length,recno)
      real*8 sbuff(length)
      integer recno
      logical, save :: firsttime
      data firsttime/.true./
      if (firsttime) then
         firsttime=.false.
         Open(131,file='../Output/cpl/ice_rec.bin',form='binary',access='direct',recl=length*8)
         Print *,'Ice record length is equal to ',length
      endif
      write(131,rec=recno)sbuff
      Print *,'Number of recorded ice units is ',length
      return
      end

      subroutine write_ocn_record(sbuff,length,recno)
      real*8 sbuff(length)
      integer recno
      logical, save :: firsttime
      data firsttime/.true./
      if (firsttime) then
         firsttime=.false.
         Open(127,file='../Output/cpl/ocn_rec.bin',form='binary',access='direct',recl=length*8)
      endif
      write(127,rec=recno)sbuff
      return
      end

      subroutine write_atm_record(sbuff,length,recno)
      real*8 sbuff(length)
      integer recno
      logical, save :: firsttime
      data firsttime/.true./
      if (firsttime) then
         firsttime=.false.
         Open(129,file='../Output/cpl/atm_rec.bin',form='binary',access='direct',recl=length*8)
      endif
      write(129,rec=recno)sbuff
      return
      end

      subroutine write_lnd_record(sbuff,length,recno)
      real*8 sbuff(length)
      integer recno
      logical, save :: firsttime
      data firsttime/.true./
      if (firsttime) then
         firsttime=.false.
         Open(133,file='../Output/cpl/lnd_rec.bin',form='binary',access='direct',recl=length*8)
      endif
      write(133,rec=recno)sbuff
      return
      end

      subroutine read_ice_record(sbuff,length,recno)
      real*8 sbuff(length)
      integer recno
      logical, save :: firsttime
      data firsttime/.true./
      if (firsttime) then
         firsttime=.false.
         Open(131,file='../Output/cpl/ice_rec.bin',form='binary',access='direct',recl=length*8)
         Print *,'Ice record length is equal to ',length
      endif
      read(131,rec=recno)sbuff
      Print *,'Number of imported ice units is ',length
      return
      end

      subroutine read_ocn_record(sbuff,length,recno)
      real*8 sbuff(length)
      integer recno
      logical, save :: firsttime
      data firsttime/.true./
      if (firsttime) then
         firsttime=.false.
         Open(127,file='../Output/cpl/ocn_rec.bin',form='binary',access='direct',recl=length*8)
      endif
      read(127,rec=recno)sbuff
      return
      end

      subroutine read_atm_record(sbuff,length,recno)
      real*8 sbuff(length)
      integer recno
      logical, save :: firsttime
      data firsttime/.true./
      if (firsttime) then
         firsttime=.false.
         Open(129,file='../Output/cpl/atm_rec.bin',form='binary',access='direct',recl=length*8)
      endif
      read(129,rec=recno)sbuff
      return
      end

      subroutine read_lnd_record(sbuff,length,recno)
      real*8 sbuff(length)
      integer recno
      logical, save :: firsttime
      data firsttime/.true./
      if (firsttime) then
         firsttime=.false.
         Open(133,file='../Output/cpl/lnd_rec.bin',form='binary',access='direct',recl=length*8)
      endif
      read(133,rec=recno)sbuff
      return
      end
