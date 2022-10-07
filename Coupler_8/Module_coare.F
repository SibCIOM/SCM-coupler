MODULE module_coare

!......................................................
!   input data - introduced in calling program "fluxes" and passed to subroutines in COMMON
!     xtime (COARE convention yyyymnddhhmmss.ss GMT)
!     u (wind speed relative to the sea surface) m/s  i.e. add CURRENT VECTOR if available
!     ts (sea temperature)  deg. C
!     t (air temperature) deg. C
!     q (specific humidity) g/kg or RH (as decimal) - code detects which then works internally in kg/kg
!     rs (shortwave radiation) W/m2
!     rl (downwelling longwave) W/m2
!     rain (average rainrate in timestep) mm/hour
!     xlat (latitude) degrees  [latitude north +ve, latitude south -ve]
!     xlon (longitude) degrees [Longitude east +ve, longitude west -ve]
!     hwt alternative sea temperature (6m depth)
!     p (pressure) mb; use 1008mb if data unavailable
!     zi (boundary-layer depth) m; use 600m if data unavailable
!     zu (wind measurement height) m   
!     zt (T measurement height) m  
!     zq (q measurement height) m  
!     zus (wind standard height) m - usually 10m
!     zts (T and q standard height) m - 
!     ts_depth (depth of ts measurement) - positive m - for warm layer correction 
!     jcool (=1 for cool skin calculation;  =0 if SST measured by IR radiometer)
!     jwarm (=1 for warm layer calculation; =0 if SST measured by IR radiometer)
!     jwave (wave state options 0=Charnock,1=Oost et al,2=Taylor and Yelland)
!....................................................................
!   calculated inputs: qs (humidity at ocean interface) and tsw (ts with warm layer added)
!....................................................................
!   outputs available - write statements in subroutine bulk_flux
!     hf W/m**2     (Sensible heat flux) - turbulent part only
!     ef W/m**2     (Latent heat flux) - turbulent part only
!     rf W/m**2     (Heat flux due to rainfall )
!     tau N/m**2   (wind stress) - turbulent part only
!     taur - momentum flux due to rain N/m**2
!     usr m/s     (M-O velocity scaling parameter u* = friction velocity)
!     tsr C       (M-O temperature scaling parameter t*)
!     qsr kg/kg   (M-O humidity scaling parameter q*)
!     Cdn - neutral drag coefficient
!     Chn - neutral transfer coefficient for heat
!     Cen - neutral transfer coefficient for moisture
!     RR - Roughness Reynolds number
!     RT - Roughness Reynolds number for temperature
!     RQ - Roughness Reynolds number for moisture
!     zL - stability parameter height/L where L is the Obukhov length
!     zo - velocity roughness length
!     zot - temperature roughness length
!     zoq - humidity roughness length
!     dt_wrm - total warm layer temperature difference C
!     tk_pwp - thickness of warm layer
!     dsea - (=dt_wrm*ts_depth/tk_pwp) warming above sea temperature measurement
!     dter - cool skin temperature difference C 
!     tkt - cool skin thickness mm (printed out *1000)
!     sst - skin temperature C (sst = ts - dter + dsea)
!     Wg - gustiness factor m/s
!     Wbar - Webb mean vertical velocity m/s
!     u_zs - wind velocity at standard height
!     t_zs - air temperature ditto
!     q_zs,qa_zs,e_zs,rh_zs - humidity ditto in a variety of units
!
IMPLICIT NONE
SAVE
!      COMMON/indata/xtime,u,ts,t,q,rs,rl,rain,xlat,xlon,qs,tsw !qs added
!      COMMON/warm/qcol_ac,tau_ac,jamset,jday1,fxp,tk_pwp
!      COMMON/const/al,be,cpa,cpw,grav,xlv,rhoa,rhow,rgas,toK,von
!      COMMON/heights/zu,zt,zq,ts_depth,zs,p,zi
!      COMMON/options/Jwarm,Jcool,Jwave
!      COMMON/rad/rns,rnl
!      COMMON/ASLout/usr,tsr,qsr,zo,zot,zoq,zL,RR,RT,RQ,RI,dter,dqer,tkt,Du,Wg
!
! F90 equivalence of common blocks
real (KIND=8) :: xtime,u,ts,t,q,rs,rl,rain,xlat,xlon,qs,qa,tsw
!real (KIND=8) :: qcol_ac,tau_ac
REAL (KIND=8) :: al,be,cpa,cpw,cpv,grav,xlv,rhoa,rhow,rgas,toK,von
real (KIND=8) :: rns,rnl
real (KIND=8) :: usr,tsr,qsr,zo,zot,zoq,zL,RR,RT,RQ,RI
real (KIND=8) :: dter, dqer,tkt,Du,Wg
real (KIND=8) :: hwave,twave
real (KIND=8) :: tau,taur,hf,ef,rf,tref
real (KIND=8) :: tau_old,hf_old,ef_old,rf_old,time_old
!some defaults
real (KIND=8) :: zu=15.               !height of wind measurement Moana Wave COARE
real (KIND=8) :: zt=15.               !height of air temp MW COARE
real (KIND=8) :: zq=15.               !height of humidity MW COARE
real (KIND=8) :: ts_depth=0.05        !Chris Fairall's floating sensor (SeaSnake)
!  default values for standard height (normally 10m) pressure and mixed layer height
real (KIND=8) :: zs=10.
real (KIND=8) :: p=1008.
real (KIND=8) :: zi=600.
integer       :: jwarm = 1,jcool = 1,jwave=0
real (KIND=8) :: qcol_ac=0.
real (KIND=8) :: tau_ac=0.
real (KIND=8) :: fxp=0.5
real (KIND=8) :: tk_pwp=19.0
integer       :: jamset=0, jday1=1

CONTAINS
!*****************************************************************
      subroutine bulk_flux(dindex)
!
implicit none
      REAL (KIND=8) :: u_zs,t_zs,q_zs,e_zs,esat_zs,qa_zs
      REAL (KIND=8) :: CD,CE,CH
      REAL (KIND=8) :: sst,ea,es
! various constants
!      REAL (KIND=8) :: visa,visw,von
! rain heat variables
      REAL (KIND=8) :: alfac,dqs_dt,dwat,dtmp
! warm layer variables
      REAL (KIND=8) :: intime,sol,sol_time,ctd1,ctd2,rich
      REAL (KIND=8) :: dsea,qjoule,qr_out,dtime
      REAL (KIND=8) :: psiu,psit
      REAL (KIND=8) :: Cdn,Cen,Chn
      REAL (KIND=8) :: Wbar,q_pwp,R,rh_zs,dt_wrm
      integer yy,mn,dd,hh,mm,ss,iter1,dindex
      character*19 chtime   !only works with 19
!      
!
      if(dindex.eq.1) then   !first line of data
       tsw=ts               !henceforth redefined after warm/cool calculations
       dter=0.3*jcool       !or in "if" block below when jwarm=0.
       sst=ts-dter          !for initial Rnl
      endif                        
!
! Constants and coefficients (Stull 1988 p640). 
      Rgas=287.1                    !J/kg/K     gas const. dry air
      toK=273.16                    ! Celsius to Kelvin
      xlv=(2.501-0.00237*ts)*1e+6   !J/kg latent heat of vaporization at ts (3C warming=0.3%)
      Cpa=1004.67                   !J/kg/K specific heat of dry air (Businger 1982)
      Cpv=Cpa*(1+0.00084*q)         !Moist air - currently not used (Businger 1982)
      rhoa=P*100./(Rgas*(t+toK)*(1.+.00061*q)) !kg/m3  Moist air density NB q still g/kg
      call gravity(xlat,grav)       ! gravity equatorial value 9.72
      al=2.1e-5*(ts+3.2)**0.79      !water thermal expansion coefft.
      be=0.026                      !salinity expansion coefft.
      cpw=4000.                     !J/kg/K specific heat water
      rhow=1022.                    !kg/m3  density water
      von=0.4                       !von Karman's "constant"
!
!   compute net radiation, updated in flux loop of ASL
      Rnl= 0.97*(5.67e-8*(sst+toK)**4-rl) !Net longwave (up = +). Typically 3C warming=15W/m2
      Rns=0.945*rs                        !Net shortwave (into water)
!
!     START Warm Layer - check switch
!            
warm_layer: if(Jwarm.eq.1) then  !jump in warm layer calculation
   intime=xtime !eg 13.35
!
!   and then to local solar time in seconds
!
   sol_time=mod(xlon/15+intime+24,24.0)*3600   !eg 85517.
   sol=sol_time/3600.
!
     index_1: if(dindex.ne.1) then        !first line of data. Set time_old and compute fluxes in ASL
     sol_t : if(sol_time.lt.time_old) then   !reset all variables at local midnight
                  jday1=0                      !reset after last 6am test of first day
                  jamset=0                     !indicates heat integration not yet started
                  tau_ac=0.0
                  qcol_ac=0.0
                  dt_wrm=0.0  
! initial guess at warm layer parameters expected in early morning
! fxp=0.5 implies a shallow heating layer to start the integration;
! tk_pwp=19.0 implies warm layer thickness is a maximum from the day
! before and is not meant to match this timestep''s fxp.
                  fxp=0.5
                  tk_pwp=19.0
                  tsw=ts
             else if(sol_time.gt.21600..and.jday1.eq.1)then   !6 am too late to start on first day
               dt_wrm=0.
               tsw=ts
             else                       !compute warm layer. Rnl and "_old"s from previous timestep
               rich=.65                                    !critical Rich. No.
               ctd1=sqrt(2*rich*cpw/(al*grav*rhow))        !u*^2 integrated so
               ctd2=sqrt(2*al*grav/(rich*rhow))/(cpw**1.5) !has /rhow in both
               dtime=sol_time-time_old                     !time step for integrals
               qr_out=rnl+hf_old+ef_old+rf_old             !total cooling at surface
               q_pwp=fxp*rns-qr_out                        !total heat absorption in warm layer
     int_thres: if(.not.(q_pwp.lt.50.and.jamset.eq.0))then       !integration threshold
                  jamset=1                                    !indicates integration has started
                  tau_ac=tau_ac+max(.002,tau_old)*dtime       !momentum integral
                    if(qcol_ac+q_pwp*dtime.gt.0)then           
                      DO iter1=1,5                           !iterate warm layer thickness
                        fxp=1.-(0.28*0.014*(1-exp(-tk_pwp/0.014))   &
                           +0.27*0.357*(1-exp(-tk_pwp/0.357))      &
                           +.45*12.82*(1-exp(-tk_pwp/12.82)))/tk_pwp !Soloviev solar absorb. prof
                        qjoule=(fxp*rns-qr_out)*dtime
                        if((qcol_ac+qjoule.gt.0.0))then
                          tk_pwp=min(19.0,ctd1*tau_ac/sqrt(qcol_ac+qjoule))  !warm layer thickness
                        end if
                      END DO
                    else                                        !warm layer wiped out
                      fxp=.75
                      tk_pwp=19
                      qjoule=(fxp*rns-qr_out)*dtime             
                    endif                                       
                qcol_ac=qcol_ac+qjoule                      !heat integral
                  if(qcol_ac.gt.0) then                       !sign check on qcol_ac
                    dt_wrm=ctd2*(qcol_ac)**1.5/tau_ac         !pwp model warming
                  else
                    dt_wrm=0.
                  endif         
! original place of sol_time        endif sol_time
                  if(tk_pwp.lt.ts_depth) then           !sensor deeper than pwp layer
                     dsea=dt_wrm                         !all warming must be added to ts
                  else                                  !warming deeper than sensor
                    dsea=dt_wrm*ts_depth/tk_pwp         !assume linear temperature profile
                  endif
                tsw=ts+dsea                           !add warming above sensor to ts
                ELSE   !int_thres
                  tsw=ts
                END IF int_thres
             END IF sol_t
           END IF index_1
           time_old=sol_time                     !all in local solar time
         ELSEIF (jwarm==0)then
           tsw=ts
         END IF warm_layer

      call humidity(t,p,ea)         !Teten's returns sat. air qa in mb
      if(q.lt.2.) then              !checks whether humidity in g/Kg or RH      
         R=q
         ea=ea*R                    !convert from RH using vapour pressure      
         q=.62197*(ea/(p-0.378*ea)) !Spec. humidity kg/kg
      else
         q=q/1000.                  !g/kg to kg/kg
      endif
      qa=.62197*(ea/(p-0.378*ea))   !convert from mb to spec. humidity  kg/kg
      call humidity(tsw,p,es)       !sea qs returned in mb
      es=es*0.98                    !reduced for salinity Kraus 1972 p. 46
      qs=.62197*(es/(p-0.378*es))   !convert from mb to spec. humidity  kg/kg
!
      call ASL(dindex)
!
! compute surface fluxes and other parameters
       sst=tsw-dter*jcool             !final skin temperature this timestep
       tau=rhoa*usr*usr*u/Du          !stress N/m2
       hf=-cpa*rhoa*usr*tsr           !sensible W/m2
       ef=-xlv*rhoa*usr*qsr           !latent W/m2
! compute heat flux due to rainfall
       dwat=2.11e-5*((t+toK)/toK)**1.94                    !water vapour diffusivity
       dtmp=(1.+3.309e-3*t-1.44e-6*t*t)*0.02411/(rhoa*cpa) !heat diffusivity
       dqs_dt=qa*xlv/(rgas*(t+toK)**2)                     !Clausius-Clapeyron
       alfac= 1/(1+0.622*(dqs_dt*xlv*dwat)/(cpa*dtmp))     !wet bulb factor
       rf= rain*alfac*cpw*((sst-t)+(qs-q-dqer)*xlv/cpa)/3600.
! compute momentum flux due to rainfall
       taur=0.85*rain/3600*u  
! Webb correction to latent heat flux already in ef via zoq/rr function so return Wbar
       Wbar=-1.61*usr*qsr/(1+1.61*q)-usr*tsr/(t+toK)
! save fluxes for next timestep warm layer integrals
       tau_old=tau 
       ef_old=ef
       hf_old=hf
       rf_old=rf
! compute transfer coefficients
       CD=(USR/Du)**2
       CH=USR*TSR/(Du*(T-sst+.0098*zt)) 
       CE=USR*QSR/(Du*(Q-QS+dqer))                                      
! compute neutral transfer coefficients and met variables at standard height
       Cdn=(0.4/log(zs/zo))**2
       Chn=0.4*0.4/(log(zs/zo)*log(zs/zot))
       Cen=0.4*0.4/(log(zs/zo)*log(zs/zoq))
! adjust met. variables to standard height
       call sub_psiu(zL,psiu)
       call sub_psit(zL,psit)
       u_zs=usr/von*(log(zs/zo)-psiu)
       t_zs=sst+tsr/von*(log(zs/zot)-psit)
       q_zs=(qs-dqer)+qsr/von*(log(zs/zoq)-psit) !kg/kg
       qa_zs=1000.*q_zs                               !g/kg
       e_zs=q_zs*p/(0.62197+0.378*q_zs)               !mb
       call humidity(t_zs,p,esat_zs)                  !mb
       rh_zs=e_zs/esat_zs                             !RH as decimal
       tref=t_zs+tok
!
! output fluxes
!     
!      write(4,200)dindex,xtime,hf,ef,sst,tau,Wbar,rf,rain, &
!      dter,dt_wrm,tk_pwp,tkt*1000.,Wg
!200   format(i6,',',f18.0,3(',',f8.2),2(',',f9.5),7(',',f8.2))
!
      end subroutine bulk_flux
!
! ------------------------------------------------------------------
      subroutine ASL( dindex)
!
! TO EVALUATE SURFACE FLUXES, SURFACE ROUGHNESS AND STABILITY OF
! THE ATMOSPHERIC SURFACE LAYER FROM BULK PARAMETERS BASED ON
! LIU ET AL. (79) JAS 36 1722-1735 
!    
implicit none
      REAL (KIND=8) :: L10,L
! constants
      REAL (KIND=8) :: visa
      REAL (KIND=8) :: visw,charn,psiu,psit
! cool skin quantities
      REAL (KIND=8) :: wetc,bigc,Bf,tcw
      REAL (KIND=8) :: hsb,hlb,alq,qcol,qout,dels,xlamx,dq,dt,ta
! Grachev and Fairall variables
      REAL (KIND=8) :: u10,zo10,cd10,ch10,ct10,zot10,cd,ct,cc,ribcu,ribu
      REAL (KIND=8) :: cwave,lwave,twopi,pst,psu,beta,zetu
INTEGER :: iter,nits,dindex
!
!
! Factors
      Beta=1.2     !Given as 1.25 in Fairall et al.(1996)
      twopi=3.14159*2.
! 
! Additional constants needed for cool skin
      visw=1.e-6                   !m2/s kinematic viscosity water
      tcw=0.6                      !W/m/K   Thermal conductivity water
      bigc=16.*grav*cpw*(rhow*visw)**3/(tcw*tcw*rhoa*rhoa)
      wetc=0.622*xlv*qs/(rgas*(tsw+toK)**2) !correction for dq;slope of sat. vap.
      visa=1.326e-5*(1+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t)   !m2/s
          !Kinematic viscosity of dry air - Andreas (1989) CRREL Rep. 89-11
! 
! Wave parameters
      cwave=grav*twave/twopi
      lwave=cwave*twave
!
! Initial guesses
      dter=0.3*jcool              !cool skin Dt
      dqer=wetc*dter              !cool skin Dq
      zo=0.0001
      Wg=0.5                      !Gustiness factor initial guess
      tkt= 0.001*jcool                  !Cool skin thickness first guess
!
! Air-sea differences - includes warm layer in Dt and Dq
      Du=(u**2.+Wg**2.)**.5       !include gustiness in wind spd. difference
      Dt=tsw-t-0.0098*zt          !potential temperature difference.
      Dq=qs-q                     
!
! **************** neutral coefficients ******************
!
      u10=Du*log(10/zo)/log(zu/zo)
      usr=0.035*u10
      zo10=0.011*usr*usr/grav+0.11*visa/usr
      Cd10=(von/log(10/zo10))**2
      Ch10=0.00115
      Ct10=Ch10/sqrt(Cd10)
      zot10=10/exp(von/Ct10)
      Cd=(von/log(zu/zo10))**2
!      
! ************* Grachev and Fairall (JAM, 1997) **********
!
      Ct=von/log(zt/zot10)         ! Temperature transfer coefficient
      CC=von*Ct/Cd                  ! z/L vs Rib linear coefficient
      Ribcu=-zu/(zi*0.004*Beta**3)  ! Saturation or plateau Rib 
      ta=t+toK
      Ribu=-grav*zu*((Dt-dter)+0.61*ta*Dq)/(ta*Du**2)
      if (Ribu.lt.0.) then
          zetu=CC*Ribu/(1+Ribu/Ribcu)   ! Unstable G and F
      else
          zetu=CC*Ribu*(1+27/9*Ribu/CC) ! Stable
      endif
      L10=zu/zetu                       ! MO length
      if (zetu.gt.50) then
        nits=1
      else
        nits=3   ! number of iterations
      endif
!
! First guess M-O stability dependent scaling params.(u*,t*,q*) to estimate zo and z/L
!
      call sub_psiu(zu/L10,psiu)
      usr= Du*von/(log(zu/zo10)-psiu)
      call sub_psit(zu/L10,psit)
      tsr=-(Dt-dter)*von/(log(zt/zot10)-psit)
      call sub_psit(zq/L10,psit)
      qsr=-(Dq-dqer)*von/(log(zq/zot10)-psit)
!      
      charn=0.011     !then modify Charnock for high wind speeds Chris' data
      if(Du.gt.10) charn=0.011+(0.018-0.011)*(Du-10)/(18-10)
      if(Du.gt.18) charn=0.018
!      
! **** Iterate across u*(t*,q*),zo(zot,zoq) and z/L including cool skin ****
!
  DO iter=1,nits
       if(Jwave.eq.0) then
        zo=charn*usr*usr/grav + 0.11*visa/usr    !after Smith 1988
       else if(Jwave.eq.1) then
        zo=(50./twopi)*lwave*(usr/cwave)**4.5+0.11*visa/usr !Oost et al.
       else if(Jwave.eq.2) then
        zo=1200.*hwave*(hwave/lwave)**4.5+0.11*visa/usr !Taylor and Yelland 
       endif 
      rr=zo*usr/visa
!
! *** zoq and zot fitted to results from several ETL cruises ************
!
      zoq=min(1.15e-4,5.5e-5/rr**0.6)
      zot=zoq
!
      call zeta(t,q,usr,tsr,qsr,zu,zL,tok,grav,von)
      L=zu/zL
      call sub_psiu(zu/L,psiu)
      dqer=wetc*dter*jcool
      usr=Du*von/(log(zu/zo)-psiu)
      call sub_psit(zt/L,psit)
      tsr=-(Dt-dter)*von/(log(zt/zot)-psit)
      call sub_psit(zq/L,psit)
      qsr=-(Dq-dqer)*von/(log(zq/zoq)-psit)
      Bf=-grav/ta*usr*(tsr+0.61*ta*qsr)
       if (Bf.gt.0) then
          Wg=Beta*(Bf*zi)**.333
       else
          Wg=0.2
       endif
         Du=sqrt(u**2.+Wg**2.)        !include gustiness in wind spd.
!      
      rnl= 0.97*(5.67e-8*(tsw-dter+toK)**4-rl)  !Recompute net longwave; cool skin=-2W/m2
!
!   Cool skin
!
 jcool_if: if(Jcool.ne.0)then
           hsb=-rhoa*cpa*usr*tsr
           hlb=-rhoa*xlv*usr*qsr
           qout=rnl+hsb+hlb
           dels=rns*(.065+11*tkt-6.6e-5/tkt*(1-exp(-tkt/8.0e-4))) !Eq.16 Ohlmann 
           qcol=qout-dels
           alq=Al*qcol+be*hlb*cpw/xlv                      !Eq. 7 Buoy flux water
           if(alq.gt.0.) then                              !originally (qcol.gt.0)
             xlamx=6/(1+(bigc*alq/usr**4)**.75)**.333      !Eq 13 Saunders coeff.
             tkt=xlamx*visw/(sqrt(rhoa/rhow)*usr)          !Eq.11 Sublayer thickness
           else
             xlamx=6.                                      !prevent excessive warm skins
             tkt=min(.01,xlamx*visw/(sqrt(rhoa/rhow)*usr)) !Limit tkt
           endif
         dter=qcol*tkt/tcw                                 ! Eq.12 Cool skin
         dqer=wetc*dter
       ENDIF jcool_if
  END DO
!
!      idum=index          ! avoids warning on compilation
      end subroutine ASL
!
!------------------------------------------------------------------
      subroutine humidity(T,P,esat)                                 
!
! Tetens' formula for saturation vp Buck(1981) JAM 20, 1527-1532 
      implicit none     
      REAL (KIND=8) :: T,P,esat
!     
      esat = (1.0007+3.46e-6*P)*6.1121*exp(17.502*T/(240.97+T)) !mb
      end subroutine humidity
!
!------------------------------------------------------------------
      subroutine sub_psiu(zL,psiu)
!
! psiu and psit evaluate stability function for wind speed and scalars
! matching Kansas and free convection forms with weighting f
! convective form follows Fairall et al (1996) with profile constants
! from Grachev et al (2000) BLM
! stable form from Beljaars and Holtslag (1991)
!
      implicit none
      REAL (KIND=8) :: zL,x,y,psik,psic,f,psiu,c
      if(zL.lt.0) then
       x=(1-15.*zL)**.25                        !Kansas unstable
       psik=2.*log((1.+x)/2.)+log((1.+x*x)/2.)-2.*atan(x)+2.*atan(1.)
       y=(1.-10.15*zL)**.3333                   !Convective
       psic=1.5*log((1.+y+y*y)/3.)-sqrt(3.)*atan((1.+2.*y)/sqrt(3.)) &
           +4.*atan(1.)/sqrt(3.)
       f=zL*zL/(1.+zL*zL)
       psiu=(1.-f)*psik+f*psic
      else
       c=min(50.,0.35*zL)                       !Stable
       psiu=-((1.+1.*zL)**1.+.6667*(zL-14.28)/exp(c)+8.525)
      endif
      return
      end subroutine sub_psiu

!--------------------------------------------------------------  
      subroutine sub_psit(zL,psit)
      implicit none
      REAL (KIND=8) :: zL,x,y,psik,psic,f,psit,c
      if(zL.lt.0) then
       x=(1-15.*zL)**.5                          !Kansas unstable
       psik=2.*log((1.+x)/2.)
       y=(1.-34.15*zL)**.3333                    !Convective
       psic=1.5*log((1.+y+y*y)/3.)-sqrt(3.)*atan((1.+2.*y)/sqrt(3.))&
           +4.*atan(1.)/sqrt(3.)
       f=zL*zL/(1.+zL*zL)
       psit=(1.-f)*psik+f*psic
      else
       c=min(50.,0.35*zL)                        !Stable
       psit=-((1.+2.*zL/3.)**1.5+.6667*(zL-14.28)/exp(c)+8.525)
      endif
       return
      end subroutine sub_psit
          
!-------------------------------------------------------------
      subroutine zeta(t,q,usr,tsr,qsr,z,zL,tok,grav,von)
!
! TO EVALUATE OBUKHOVS STABILITY PARAMETER z/L FROM AVERAGE
! TEMP T IN DEG C, AVERAGE HUMIDITY Q IN GM/GM, HEIGHT IN M,
! AND FRICTIONAL VEL,TEMP.,HUM. IN MKS UNITS SEE LIU ET AL.(JAM 1979,36,1722-1735)
! which actually contains slight errors in these equations (EFB)
!     
      IMPLICIT NONE
      REAL (KIND=8) :: t,q,ob,tvsr,tv,ta,sgn
      REAL (KIND=8) :: usr,tsr,qsr,z,zL
      REAL (KIND=8) :: grav,toK,von
      ta=t+toK
      tv=ta*(1.+0.61*q)
      tvsr=tsr*(1.+0.61*q)+0.61*ta*qsr    
      sgn=sign(1.,tvsr)               !added this to avoid program
      if(abs(tvsr) .lt. 1.e-3) then   !failure when TVSR is very small
         tvsr=sgn*tvsr
      endif
      ob=tv*usr*usr/(grav*von*tvsr) 
      zL=z/ob                     
      end subroutine zeta
       
!-------------------------------------------------------------------------
      Subroutine gravity(lat,g)
!       calculates g as a function of latitude using the 1980 IUGG formula
!         
!       Bulletin Geodesique, Vol 62, No 3, 1988 (Geodesist's Handbook)
!       p 356, 1980 Gravity Formula (IUGG, H. Moritz)
!       units are in m/sec^2 and have a relative precision of 1 part
!       in 10^10 (0.1 microGal)
!       code by M. Zumberge.
!
!       check values are:
!
!        g = 9.780326772 at latitude  0.0
!        g = 9.806199203 at latitude 45.0
!        g = 9.832186368 at latitude 90.0
!
      implicit none
      REAL (KIND=8) :: gamma, c1, c2, c3, c4, phi, lat, g
      gamma = 9.7803267715
      c1 = 0.0052790414
      c2 = 0.0000232718
      c3 = 0.0000001262
      c4 = 0.0000000007
      phi = lat * 3.14159265358979 / 180.0
      g = gamma * (1.0 &
      + c1 * ((sin(phi))**2) &
      + c2 * ((sin(phi))**4) &
      + c3 * ((sin(phi))**6) &
      + c4 * ((sin(phi))**8))
!
      end subroutine gravity
 
END MODULE module_coare 
 
 
