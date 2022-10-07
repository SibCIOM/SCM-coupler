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
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      parameter( latvap =  2.5e6          )  ! for evaporation from water surface
c      parameter( latvap =  3.337e5 + 2.5e6)  ! for evaporation from ice surface
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        rdn    = sqrt(cdn(vmag))           ! for open water surface
c        rdn    = xkar/log(zref/0.0005)     ! for ice surface
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
        sen (i) =      cp * tau * tstar / ustar 
        lat (i) =  latvap * tau * qstar / ustar
        lwup(i) = - stebol * ts(i)**4
     
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
         lat(i)=0.0
         lwup(i)=0.0
         evap(i)=0.0
         taux(i)=0.0
         tauy(i)=0.0
         tref(i)=0.0
        endif

      enddo 

      return 
      end
 
c=======================================================================
