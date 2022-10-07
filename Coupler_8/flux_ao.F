c=======================================================================
!     $Id: flux_ao.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $
c=======================================================================

      SUBROUTINE flux_ao(imt, jmt, A, F, Ogr)

      use struct
      use module_coare
      
	Type (grid) Ogr
	Type (AO) A
	Type (AO_flux) F

      integer i,j                  ! loop index
      integer nu_diag,cur_date,cur_sec
      common/iodiag/nu_diag
      common/xtiming/cur_date,cur_sec

c-----------------------------------------------------------------------
c PURPOSE:
c   compute certain atm/ocn fluxes (on the ocn grid)
c
c INPUT:
c   Sa  : atm state variables 
c   So  : ocn state variables 
c
c OUTPUT:
c   Faoc: atm/ocn fluxes
c
c NOTE: 
c   net heat flux = net solar + lw up + lw down + sensible + latent + Q
c
c ASSUMPTIONS:
c   all fields are on ocn grid
c
c-----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE(j)
      integer dindex
      real*8 cld_power, solar_avg
      logical use_ciaf, tidal, use_frac, use_frcp, use_tref, use_xie,
     *   use_xie4ncep, use_datahum, old_tides, aomip_rules, prognostic_run,
     *   use_small_rivers, use_coare
      integer ciaf_vers, cld_lw_vers, cld_sw_vers, first_riv_year,last_riv_year
      character*80 file_small_riv
      common/pars/cld_power,use_ciaf,ciaf_vers,tidal,
     &  use_frac,use_frcp,use_tref,
     &  cld_lw_vers,cld_sw_vers,use_xie,use_xie4ncep,
     &  use_datahum,old_tides,aomip_rules,prognostic_run,solar_avg,
     &  use_small_rivers,file_small_riv,first_riv_year,last_riv_year,
     &  use_coare
!      real*8, allocatable :: wo(:,:)

      dindex=1
      F%sen = 0.0
      F%dsen = 0.0
      F%lat = 0.0
      F%lwup = 0.0
      F%dlwup = 0.0
      F%evap = 0.0
      F%taux = 0.0
      F%tauy = 0.0
      F%tref = 0.0
      F%fco2 = 0.0

      if(.not.use_coare)then
      DO j=1,jmt
         call srfflx_ao(imt,
     &               Ogr%mask (1,j),
     $               A%z  (1,j),
     $               A%u  (1,j),
     $               A%v  (1,j),
     $               A%ptem (1,j),
     $               A%shum (1,j),
     $               A%dens (1,j),
     $               A%tbot (1,j),
     $               A%sst(1,j),
     $               F%sen (1,j), 
     $               F%dsen (1,j), 
     $               F%lat (1,j),
     $               F%lwup(1,j), 
     $               F%dlwup(1,j), 
     $               F%evap(1,j),
     $               F%taux(1,j),
     $               F%tauy(1,j),
     $               F%tref(1,j) )
      END DO
c      Write(nu_diag,*) 'u=',A%u(87,130),A%u(88,130),A%u(89,130)
c      Write(nu_diag,*) 'v=',A%v(87,130),A%v(88,130),A%v(89,130)
c      Write(nu_diag,*) 'taux=',F%taux(87,130),F%taux(88,130),F%taux(89,130)
c      Write(nu_diag,*) 'tauy=',F%tauy(87,130),F%tauy(88,130),F%tauy(89,130)
c      read *
      else
!     xtime (hours GMT)
         xtime=dble(cur_sec)/3600.0
         do j=1,jmt
            do i=1,imt
               if(Ogr%mask(i,j).ne.0)then
!     u (wind speed relative to the sea surface) m/s  i.e. add CURRENT VECTOR if available
               u=sqrt(A%u(i,j)**2+A%v(i,j)**2)
!     ts (sea temperature)  deg. C
               ts=A%sst(i,j)-273.16
!     t (air temperature) deg. C
               t=A%tbot(i,j)-273.16
!     q (specific humidity) g/kg or RH (as decimal) - code detects which then works internally in kg/kg
               q=A%shum(i,j)*1000.0
!     rs (shortwave radiation) W/m2
               rs=A%swdn(i,j)
!     rl (downwelling longwave) W/m2
               rl=A%lwdn(i,j)
!     rain (average rainrate in timestep) mm/hour
               rain=A%prec(i,j)*3600.0/1.022
!     xlat (latitude) degrees  [latitude north +ve, latitude south -ve]
               xlat=Ogr%lat(i,j)
!     xlon (longitude) degrees [Longitude east +ve, longitude west -ve]
               xlon=Ogr%lon(i,j)
               if(xlon.gt.180.0)xlon=xlon-360.0
               if(xlon.lt.-180.0)xlon=xlon+360.0
!     hwt alternative sea temperature (6m depth)
               htw=ts
!     p (pressure) mb; use 1008mb if data unavailable
               p=1008.0
!     zi (boundary-layer depth) m; use 600m if data unavailable
               zi=600.0
!     zu (wind measurement height) m   
               zu=A%z(i,j)
!     zt (T measurement height) m  
               zt=A%zt(i,j)
!     zq (q measurement height) m  
               zq=zt
!     zus (wind standard height) m - usually 10m
               zus=10.0
!     zts (T and q standard height) m - 
               zts=zus
!     ts_depth (depth of ts measurement) - positive m - for warm layer correction 
               ts_depth=2.5
!     jcool (=1 for cool skin calculation;  =0 if SST measured by IR radiometer)
               jcool=1
!     jwarm (=1 for warm layer calculation; =0 if SST measured by IR radiometer)
               jwarm=1
!     jwave (wave state options 0=Charnock,1=Oost et al,2=Taylor and Yelland)
               jwave=0
               call bulk_flux(dindex)
               F%sen (i,j)=-(hf+rf)
               F%dsen (i,j)=0.0
               F%lat (i,j)=-ef
               F%lwup(i,j)=-(rnl+0.97*rl)
               F%dlwup(i,j)=0.0
               F%evap(i,j)=-ef/xlv
               F%taux(i,j)=0.0
               F%tauy(i,j)=0.0
               if(u.gt.0.0)then
                  F%taux(i,j)=(tau+taur)*A%u(i,j)/u
                  F%tauy(i,j)=(tau+taur)*A%v(i,j)/u
               endif
               F%tref(i,j)=tref
               endif
            enddo
         enddo

      endif

!      allocate(wo(imt,jmt))
!      do i=1,imt
!         do j=1,jmt
!!            wo(i,j)=F%lwup(i,j)+F%dlwup(i,j)*A%sst(i,j)
!            wo(i,j)=F%lat(i,j)
!         enddo
!      enddo
!      call open_w2w('lat1.nc','lath',imt,jmt)
!      call write_w('lat1.nc',1,19480101,'lath',wo,Ogr%mask,imt,jmt)
!      deallocate(wo)
!
!      stop

      SCHM=450.
      SCH0=660.
      AK00=0.83
      do i=1,imt
         do j=1,jmt
            po=PCO2(A%oco2(i,j),A%sst(i,j))
            AKK=(-5.0+2.4*amax1(wind,4.0))/sqrt(SCHM/SCH0)                  ! AKK IN CM/HOUR, WIND IN M/S
!               C08=12.E-6*(PATM-POC)/29.*AK00*AKK/3600.*1.2E-3             ! G/(CM**2*S)
            F%fco2(i,j)=-12.E-5*(1.0e+5*A%aco2(i,j)/A%dens (i,j)-po)/29.*AK00*AKK/3600.*1.2E-3   ! KG/(M**2*S)
         enddo
      enddo

      RETURN
      END

c=======================================================================


      FUNCTION PCO2(CC,TT)
      IMPLICIT NONE
      REAL C,CC,CCORB,CCO2,HCO3,CO3,AS,PCO2
      REAL TT

      AS=2.22

      C=CC/0.012 ! CC IN KG/M**3, C IN MOLES/M**3

      CCORB=0.000575+0.000006*(TT-5.0)
      CCO2=0.035+0.0019*(TT-5.0)

      HCO3=(C-SQRT(C*C-AS*(2.*C-AS)*(1.-4.*CCORB)))/
     /(1.-4.*CCORB)

      CO3=(AS-HCO3)/2.

      PCO2=280.*CCO2*HCO3**2/CO3

      RETURN
      END

      subroutine srfflx_ao(imax, mask, zbot  ,ubot  ,vbot  ,thbot ,
     $                     qbot  ,rbot  ,tbot  ,ts    ,sen   ,dsen,
     $                     lat   ,lwup  ,dlwup ,evap  ,taux  ,tauy,
     $                     tref  )

      implicit none

      !--- input arguments --------------------------------
      integer imax        ! array dimensions
	integer*4 mask(imax)
      real    zbot (imax) ! atm level height  (m)
      real    ubot (imax) ! atm u wind        (m/s)
      real    vbot (imax) ! atm v wind        (m/s)
      real    thbot(imax) ! atm potential T   (K)
      real    qbot (imax) ! atm specific humidity (kg/kg)
      real    rbot (imax) ! atm air density   (kg/m^3)
      real    tbot (imax) ! atm T             (K) 
      real    ts   (imax) ! surface temperature

      !--- output arguments -------------------------------
      real    sen  (imax) ! sensible         heat flux  (W/m^2)
      real    dsen (imax) ! sensible         heat flux  (W/m^2)
      real    lat  (imax) ! latent           heat flux  (W/m^2)
      real    lwup (imax) ! long-wave upward heat flux  (W/m^2)
      real    dlwup(imax) ! long-wave upward heat flux  (W/m^2)
      real    evap (imax) ! evaporative water flux ((kg/s)/m^2)
      real    taux (imax) ! x surface stress (N)
      real    tauy (imax) ! y surface stress (N)
      real    tref (imax) ! 2m reference height temperature
 
      !--- constants --------------------------------------
      real    cpair       ! Specific heat of dry air
      real    cpwv        ! Specific heat of water vapor
      real    cpvir       ! Defined as cpwv/cpair - 1.
      real    gravit      ! Acceleration of gravity
      real    stebol      ! Stefan-Boltzmann's constant
      real    xkar        ! Von Karman constant
      real    umin        ! minimum wind speed (m/s)
      real    zref        ! reference height (m)
      real    zvir        ! rh2o/rair - 1.0
      real    latvap      ! latent heat of evaporation

      parameter( cpair  =  1.005e3        )
      parameter( cpwv   =  1.81e3         )
      parameter( cpvir  =  cpwv/cpair - 1.)
      parameter( gravit =  9.80616        )
      parameter( stebol = 56.7e-9         )
      parameter( xkar   =  0.4            )
      parameter( umin   =  1.0            )
      parameter( zref   = 10.0            )
      parameter( zvir   =  0.606          )
      parameter( latvap =  2.5e6          )

      !--- local variables --------------------------------
      integer i      ! vector index
      real    vmag   ! surface wind magnitude   (m/s)
      real    thvbot ! virtual temperature      (K)
      real    ssq    ! sea surface humidity     (kg/kg)
      real    delt   ! potential T difference   (K)
      real    delq   ! humidity difference      (kg/kg)
      real    stable ! stability factor
      real    rdn    ! sqrt of neutral exchange coefficient (momentum) 
      real    rhn    ! sqrt of neutral exchange coefficient (heat)     
      real    ren    ! sqrt of neutral exchange coefficient (water)    
      real    rd     ! sqrt of exchange coefficient (momentum)         
      real    rh     ! sqrt of exchange coefficient (heat)             
      real    re     ! sqrt of exchange coefficient (water)            
      real    ustar  ! ustar             
      real    qstar  ! qstar             
      real    tstar  ! tstar             
      real    hol    ! H (at zbot) over L
      real    xsq    ! temporary variable                         
      real    xqq    ! temporary variable                         
      real    psimh  ! stability function at zbot (momentum)      
      real    psixh  ! stability function at zbot (heat and water)
      real    alz    ! ln(zbot/z10)
      real    u10n   ! 10m neutral wind 
      real    tau    ! stress at zbot
      real    cp     ! specific heat of moist air

      real    bn     ! exchange coef funct for interpolation
      real    bh     ! exchange coef funct for interpolation
      real    fac    ! interpolation factor
      real    ln0    ! log factor for interpolation
      real    ln3    ! log factor for interpolation
      real    ztref  ! reference height for air temperature
      parameter (ztref=2.0)

      !--- local functions --------------------------------
      real    Tk     ! temperature (K)
      real    esat   ! the saturation vapor pressure (dPa)
      real    qsat   ! the saturation humididty of air (kg/m^3)
      real    Umps   ! wind velocity (m/s)
      real    cdn    ! function: neutral drag coeff at 10m
      real    psimhu ! unstable part of psimh
      real    psixhu ! unstable part of psimx
      real    xd     ! dummy argument  
      integer nu_diag
      common/iodiag/nu_diag
 
c      qsat(Tk)   = 640380. / exp(5107.4/Tk)
      esat(Tk) = 0.98*exp((0.7859+0.03477*Tk)/
     *   (1.0+0.00412*Tk)*2.302585093)
      qsat(Tk) = 0.622*esat(Tk-273.16)/(1000.0-0.378*esat(Tk-273.16))
 
      cdn(Umps)  = 0.0027 / Umps + .000142 + .0000764 * Umps
 
      psimhu(xd) = log((1.+xd*(2.+xd))*(1.+xd*xd)/8.)
     $             - 2.*atan(xd) + 1.571
 
      psixhu(xd) = 2. * log((1. + xd*xd)/2.)
 
c-----------------------------------------------------------------------
c PURPOSE:
c   computes certain atm/ocn surface fluxes
c
c NOTE: 
c   o all fluxes are positive downward
c   o net heat flux = net sw + lw up + lw down + sen + lat
c   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
c   o wind speeds should all be above a minimum speed (say 1.0 m/s)
c 
c ASSUME:
c   o Neutral 10m drag coeff: cdn = .0027/U10 + .000142 + .0000764 U10
c   o Neutral 10m stanton number: ctn = .0327 sqrt(cdn), unstable
c                                 ctn = .0180 sqrt(cdn), stable
c   o Neutral 10m dalton number:  cen = .0346 sqrt(cdn)
c   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
c-----------------------------------------------------------------------
 
      do i=1,imax

        !--- compute some needed quantities ---
	if(mask(i).ne.0)then
        vmag   = max(umin, sqrt(ubot(i)**2+vbot(i)**2))
        thvbot = thbot(i) * (1. + zvir * qbot(i)) ! virtual temp (K)
c        ssq    = 0.98 * qsat(ts(i)) / rbot(i)     ! sea surf hum (kg/kg)
        ssq    = qsat(ts(i))                      ! sea surf hum (kg/kg)
        delt   = thbot(i) - ts(i)                 ! pot temp diff (K)
        delq   = qbot(i) - ssq                    ! spec hum dif (kg/kg)
c        if(delq.gt.0.0)delq=0.0
        alz    = log(zbot(i)/zref) 
        cp     = cpair*(1. + cpvir*ssq) 
	
        !------------------------------------------------------------
        !     first estimate of Z/L and ustar, tstar and qstar
        !------------------------------------------------------------

        !--- neutral coefficients, z/L = 0.0 ---
        stable = 0.5 + sign(0.5 , delt)
        rdn    = sqrt(cdn(vmag))
        rhn    = (1.-stable) * 0.0327 + stable * 0.018 
        ren    = 0.0346 

        !--- ustar, tstar, qstar ---
        ustar = rdn * vmag
        tstar = rhn * delt  
        qstar = ren * delq  

        !--- compute stability & evaluate all stability functions ---
        hol    = xkar * gravit * zbot(i)
     $           * (tstar/thvbot+qstar/(1./zvir+qbot(i))) / ustar**2
        hol    = sign( min(abs(hol),10.), hol )
        stable = 0.5 + sign(0.5 , hol)
        xsq    = max(sqrt(abs(1. - 16.*hol)) , 1.)
        xqq    = sqrt(xsq)
        psimh  = -5.*hol*stable + (1.-stable)*psimhu(xqq)
        psixh  = -5.*hol*stable + (1.-stable)*psixhu(xqq)

        !--- shift wind speed using old coefficient ---
        rd   = rdn / (1.+rdn/xkar*(alz-psimh))
        u10n = vmag * rd / rdn 

        !--- update transfer coeffs at 10m and neutral stability ---
        rdn = sqrt(cdn(u10n))
        ren = 0.0346
        rhn = (1.-stable)*0.0327 + stable * 0.018 

        !--- shift all coeffs to measurement height and stability ---
        rd = rdn / (1.+rdn/xkar*(alz-psimh)) 
        rh = rhn / (1.+rhn/xkar*(alz-psixh)) 
        re = ren / (1.+ren/xkar*(alz-psixh)) 

        !--- update ustar, tstar, qstar using updated, shifted coeffs --
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 

        !------------------------------------------------------------
        !     iterate to converge on Z/L, ustar, tstar and qstar
        !------------------------------------------------------------

        !--- compute stability & evaluate all stability functions ---
        hol    = xkar * gravit * zbot(i)
     $           * (tstar/thvbot+qstar/(1./zvir+qbot(i))) / ustar**2
        hol    = sign( min(abs(hol),10.), hol )
        stable = 0.5 + sign(0.5 , hol)
        xsq    = max(sqrt(abs(1. - 16.*hol)) , 1.)
        xqq    = sqrt(xsq)
        psimh  = -5.*hol*stable + (1.-stable)*psimhu(xqq)
        psixh  = -5.*hol*stable + (1.-stable)*psixhu(xqq)

        !--- shift wind speed using old coeffs ---
        rd   = rdn / (1.+rdn/xkar*(alz-psimh))
        u10n = vmag * rd / rdn 

        !--- update transfer coeffs at 10m and neutral stability ---
        rdn = sqrt(cdn(u10n))
        ren = 0.0346
        rhn = (1.-stable)*0.0327 + stable * 0.018 

        !--- shift all coeffs to measurement height and stability ---
        rd = rdn / (1.+rdn/xkar*(alz-psimh)) 
        rh = rhn / (1.+rhn/xkar*(alz-psixh)) 
        re = ren / (1.+ren/xkar*(alz-psixh)) 

        !--- update ustar, tstar, qstar using updated, shifted coeffs --
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 

        !------------------------------------------------------------
        !                    compute the fluxes
        !------------------------------------------------------------

        tau = rbot(i) * ustar * ustar 
      
        !--- momentum flux ---
        taux(i) = tau * ubot(i) / vmag 
        tauy(i) = tau * vbot(i) / vmag 
       
        !--- heat flux ---
        sen (i) =      cp * tau * rh * thbot(i) / ustar 
        dsen (i) =     - cp * tau * rh / ustar 
        lat (i) =  latvap * tau * qstar / ustar
        dlwup(i) = -4.0*stebol * ts(i)**3
        lwup(i) = 3.0*stebol * ts(i)**4
     
        !--- water flux ---
        evap(i) = lat(i)/latvap 

        !------------------------------------------------------------
        !         compute 2m reference height temperature
        !------------------------------------------------------------

        !--- Compute function of exchange coefficients. Assume that 
        !--- cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore 
        !--- 1/sqrt(cn(i))=1/rdn and sqrt(cm(i))/ch(i)=1/rh 
         bn = xkar/rdn
         bh = xkar/rh

        !--- Interpolation factor for stable and unstable cases
         ln0 = log(1.0 + (ztref/zbot(i))*(exp(bn) - 1.0))
         ln3 = log(1.0 + (ztref/zbot(i))*(exp(bn - bh) - 1.0))
         fac = (ln0 - ztref/zbot(i)*(bn - bh))/bh * stable
     $       + (ln0 - ln3)/bh * (1.-stable)
         fac = min(max(fac,0.),1.)

        !--- Actual interpolation
         tref(i) = ts(i) + (tbot(i) - ts(i))*fac
	else
      sen(i)=0.0
      dsen(i)=0.0
      lat(i)=0.0
      lwup(i)=0.0
      dlwup(i)=0.0
      evap(i)=0.0
      taux(i)=0.0
      tauy(i)=0.0
      tref(i)=0.0
      endif

      enddo 

      return 
      end
 
c=======================================================================
