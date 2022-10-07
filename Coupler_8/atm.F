!     $Id: atm.F 5 2015-03-31 07:04:05Z plat $
!     $HeadURL$

      subroutine init_atm(Atm, Ast, Ocn, Ost)

      use struct
      Type (grid) Atm, Ocn
      Type (state) Ast, Ost
      integer nu_diag
      common/iodiag/ nu_diag
      
      integer*4 i,j

      if(Ast%imt.ne.Ost%imt.or.Ast%jmt.ne.Ost%jmt)then
         write(nu_diag,*) 'atmosphere grid is different from ocn'
         write(nu_diag,*) 'make atmosphere module more smart'
c         call exit_coupler
      endif

      Do j=1,Ast%jmt
      Do i=1,Ast%imt
         Atm%lon(i,j)    = Ocn%lon(i,j)
         Atm%lat(i,j)    = Ocn%lat(i,j)
         Atm%tarea(i,j)  = Ocn%tarea(i,j)
         Atm%mask(i,j)   = Ocn%mask(i,j)
      enddo
      enddo
      
      return
      end

      subroutine get_atm_old_state(current, V, G, S)
      use struct
      Type (state) S
      Type (grid) G
      Type (from_atm) V

      real*8, allocatable, save, dimension(:,:,:) :: prec,cld,newprec
      real*8, allocatable, save, dimension(:,:) :: pdata,udata,vdata,tdata,
     *   humdata,prdata,rhum
      real*8 curr, tm, pp
      character*80 fname,file_small_riv
      logical, save :: firstcall,useciaf,oldstuff
      integer, save :: nold
      integer nu_diag, n
      logical, external :: isnan8
      logical use_ciaf,tidal,use_frac,use_frcp,use_xie,use_xie4ncep,
     &     use_datahum, atm_shift, old_tides,aomip_rules,prognostic_run,
     &     use_tref,use_small_rivers
      integer ciaf_vers,cld_lw_vers,cld_sw_vers,atm_nadv,td_nadv,
     &   first_riv_year,last_riv_year
      real*8 solar_avg
      common/pars/cld_power,use_ciaf,ciaf_vers,tidal,
     &  use_frac,use_frcp,use_tref,
     &  cld_lw_vers,cld_sw_vers,use_xie,use_xie4ncep,
     &  use_datahum,old_tides,aomip_rules,prognostic_run,solar_avg,
     &  use_small_rivers,file_small_riv,first_riv_year,last_riv_year
      logical use_corrected,firsttime,opennc,use_mask
      integer ciafvers
      integer, save :: idate0,nrec
      character*80 datapath,rpath,respath
      character*1 writemode
      common/progn_run/datapath,respath,use_mask,writemode
      common/iodiag/ nu_diag
      data firstcall/.true./,nold/0/,kmold/0/,firsttime/.true./,nsolar_old/0/,
     *   idates/0/,opennc/.true./
      

      use_corrected=.false.

      pp=atan2(1.0,1.0)*4.0/180.0

      iyear=int(S%idate/10000)

c    Comment next 10 lines and also one line in the end to make normal run
c      if(iyear.gt.1970)then
c         idate_old=S%idate
c         if(iyear.lt.1980)then
c            iyear_new=1970-(iyear-1970)
c         else
c            iyear_new=iyear-20
c         endif
c         S%idate=iyear_new*10000+(idate_old-10000*iyear)
c         iyear=iyear_new
c      endif

      if(.not.prognostic_run)then

      ciafvers=-1
      if(use_ciaf)then
         if(ciaf_vers.eq.0)
     *      useciaf=(use_ciaf.and.iyear.ge.1948.and.iyear.le.2500)
         if(ciaf_vers.eq.1)
     *      useciaf=(use_ciaf.and.iyear.ge.1958.and.iyear.le.2004)
         if(ciaf_vers.eq.2)
     *      useciaf=(use_ciaf.and.iyear.ge.1948.and.iyear.le.2007)
         if(useciaf) ciafvers=ciaf_vers
      endif
c      useciaf=.false.
      
c      Write(nu_diag,*) 'In Get_Atm_State'

      imt=S%imt
      jmt=S%jmt
      if(firstcall)then
         firstcall=.false.
         allocate(prec(imt,jmt,12),cld(imt,jmt,12),newprec(imt,jmt,12),
     *      pdata(imt,jmt),udata(imt,jmt),vdata(imt,jmt),tdata(imt,jmt),
     *      humdata(imt,jmt),prdata(imt,jmt),rhum(imt,jmt))

         call clouds(S%imt,S%jmt,G%lon,G%lat,G%mask,cld,cld_power)

         call precip(S%idate,S%imt,S%jmt,G%lon,G%lat,G%mask,
     *      prec,use_xie,use_xie4ncep,ciafvers,oldstuff)

         if(use_corrected)then
         call newprecip(S%imt,S%jmt,G%lon,G%lat,G%mask,newprec)
         do k=1,12
c            Wmax=0.
c            Wmin=0.
            do j=1,jmt
               do i=1,imt
                  if(G%lat(i,j)/pp.ge.60.0)then
                     if(G%lat(i,j)/pp.gt.65.0)then
                        prec(i,j,k)=newprec(i,j,k)
                     else
                        prec(i,j,k)=(
     *                     (65.0-G%lat(i,j)/pp)*prec(i,j,k)+
     *                     (G%lat(i,j)/pp-60.0)*newprec(i,j,k)
     *                      )/5.0
                     endif
                  endif
c                  if(Wmax.lt.prec(i,j,k))Wmax=prec(i,j,k)
c                  if(Wmin.gt.prec(i,j,k))Wmin=prec(i,j,k)
               enddo
            enddo
c            Write(nu_diag,'(a2,i2,a8,e12.4,a8,e12.4)')
c     *         'pr',k,': Max = ',Wmax,', Min = ',Wmin
         enddo
         endif
         
c       write(nu_diag,*)'First time initialization is finished'
      else
         call precip(S%idate,S%imt,S%jmt,G%lon,G%lat,G%mask,
     *      prec,use_xie,use_xie4ncep,ciafvers,oldstuff)
         if((.not.oldstuff).and.(use_corrected))then
            do k=1,12
c               Wmax=0.
c               Wmin=0.
               do j=1,jmt
                  do i=1,imt
                     if(G%lat(i,j)/pp.ge.60.0)then
                        if(G%lat(i,j)/pp.gt.65.0)then
                           prec(i,j,k)=newprec(i,j,k)
                        else
                           prec(i,j,k)=(
     *                        (65.0-G%lat(i,j)/pp)*prec(i,j,k)+
     *                        (G%lat(i,j)/pp-60.0)*newprec(i,j,k)
     *                        )/5.0
                        endif
                     endif
c                     if(Wmax.lt.prec(i,j,k))Wmax=prec(i,j,k)
c                     if(Wmin.gt.prec(i,j,k))Wmin=prec(i,j,k)
                  enddo
               enddo
c              Write(nu_diag,'(a2,i2,a8,e12.4,a8,e12.4)')
c     *          'pr',k,': Max = ',Wmax,', Min = ',Wmin
            enddo
         endif
      endif

      crnt=current/86400.0
      n=min(int(crnt-int(crnt/365.0)*365.0)+1,365)
      timesolar=(crnt-int(crnt/365.0)*365.0-dble(n-1))*24.0    ! time of the day in hours
      nsolar=int(timesolar/solar_avg)+1                        ! number of the solar average period during a day
      nsolar_tot=nsolar+int((n-1)*24.0/solar_avg)              ! number of the solar average period during a year

      curr=crnt
      rhofresh=1000.0
      
      kmo=int((S%idate-10000*iyear)/100)
c      if(n.ne.nold)then
      if(kmo.ne.kmold)then
         tm=(curr-int(curr/365.0)*365.0)*360.0/365.0-15.0
         if(tm.lt.0)then
            k1=12
            k2=1
            r1=-tm/30.0
         else
            k1=int(tm/30.0)+1
            k2=k1+1
            kpr1=k1
	    if(k2.gt.12)k2=k2-12
            r1=(30.0*k1-tm)/30.0
         endif
         r2=1.0-r1
         do j=1,jmt
            do i=1,imt
               if(G%mask(i,j).eq.1)then
c                  V%Fsnow(i,j)=(r1*prec(i,j,kpr1)+r2*prec(i,j,k2))*
c     *               rhofresh
c                  V%Frain(i,j)=0.0
c                  V%cld(i,j)=r1*cld(i,j,k1)+r2*cld(i,j,k2)
                  V%Fsnow(i,j)=prec(i,j,kmo)*rhofresh
                  V%Frain(i,j)=0.0
                  V%cld(i,j)=cld(i,j,kmo)
c                  V%cld(i,j)=V%cld(i,j)+
c     *               2.0*V%cld(i,j)*(1.0-V%cld(i,j))/3.0
               else
                  V%Fsnow(i,j)=0.0
                  V%Frain(i,j)=0.0
                  V%cld(i,j)=0.0
               endif
            enddo
         enddo
      end if

      if(useciaf)then

         if(firsttime)then
         print *,'CIAF data are in use: v.',ciaf_vers
         firsttime=.false.
         endif
         call ciaf_pressure(S%idate,S%sec,S%imt,S%jmt,G%lon,G%lat,
     *      G%mask,pdata,ciaf_vers)
         call ciafwind(S%idate,S%sec,S%imt,S%jmt,G%lon,G%lat,G%mask,
     *      udata,vdata,ciaf_vers)
         call ciaf_airtemp(S%idate,S%sec,S%imt,S%jmt,G%lon,G%lat,
     *      G%mask,tdata,ciaf_vers)
         if(use_datahum)then
         call ciaf_airhum(S%idate,S%sec,S%imt,S%jmt,G%lon,G%lat,
     *      G%mask,humdata,ciaf_vers)
         endif

      else

         if(.not.firsttime)then
         print *,'CIAF data are not used'
         firsttime=.true.
         endif
         if(n.ne.nold)then
           call pressure(S%idate,S%imt,S%jmt,G%lon,G%lat,G%mask,pdata)
           call wind(S%idate,S%imt,S%jmt,G%lon,G%lat,G%mask,udata,vdata)
            call airtemp(S%idate,S%imt,S%jmt,G%lon,G%lat,G%mask,tdata)
            if(use_datahum)then
            call rhumid(S%idate,S%imt,S%jmt,G%lon,G%lat,G%mask,rhum)
            endif
         endif

      endif
         
c         Write(nu_diag,*)'finish data reading'

      if(n.ne.nold.or.useciaf)then

         do j=1,jmt
            j0=max(j-1,1)
            j1=min(j+1,jmt)
            do i=1,imt
               i0=max(i-1,1)
               i1=min(i+1,imt)
               if(G%mask(i,j).eq.1)then
                  V%uatm(i,j)=udata(i,j)
                  V%vatm(i,j)=vdata(i,j)
c
c      The following commented out lines are to work with 2m temperature
c
c                  V%potT(i,j)=tdata(i,j)*(1.e+5/pdata(i,j))**0.29
c                  V%Tair(i,j)=tdata(i,j)*(0.995e+5/pdata(i,j))**0.29
c
c      The following lines are to work with 995 potential temperature
c
                  if(useciaf)then
                     V%zlvl(i,j)=10.0
                     V%Tair(i,j)=tdata(i,j)
                     V%potT(i,j)=tdata(i,j)*(1000.0/999.0)**0.29
                  else
                     V%zlvl(i,j)=50.0
                     V%potT(i,j)=tdata(i,j)
                     V%Tair(i,j)=tdata(i,j)*(995.0/1000.0)**0.29
                  endif
                  Tk=V%Tair(i,j)-273.16
                  esat = 0.98*exp((0.7859+0.03477*Tk)/
     *               (1.0+0.00412*Tk)*2.302585093)
                  qsat = 0.622*esat/(1000.0-0.378*esat)
                  if(use_datahum)then
                     if(useciaf)then
                        V%Qa(i,j)=humdata(i,j)
                     else
                        V%Qa(i,j)=rhum(i,j)*qsat
                     endif
                  else
                     V%Qa(i,j)=0.9*qsat
                  endif
                  
                  QaV=V%Qa(i,j)
                  if(QaV.le.0.0)QaV=0.0

                  V%rhoa(i,j)=pdata(i,j)/(287.04*V%Tair(i,j))

c   Note: AOMIP Long Wave radiative flux is modifyed in prepare.F

               endif
            enddo
         enddo
      endif

      if(kmo.ne.kmold .or. n.ne.nold .or. useciaf)then

         do j=1,jmt
            do i=1,imt
               if(G%mask(i,j).eq.1)then
                  Rate=V%Frain(i,j)+V%Fsnow(i,j)
                  if(V%Tair(i,j).ge.273.16)then
                     V%Frain(i,j)=Rate
                     V%Fsnow(i,j)=0.0
                  else
                     V%Frain(i,j)=0.0
                     V%Fsnow(i,j)=Rate
                  endif
               endif
            enddo
         enddo
      endif
c   AOMIP Short Wave radiative flux

      if(nsolar_tot.ne.nsolar_old)then
         SS=0.0
         nSS=0
         do j=1,jmt
            do i=1,imt
               if(G%mask(i,j).eq.1)then
                  DoY=dble(n-1)*360.0/365.0
                  delta=23.44*cos((172-DoY)*pp)
                  Qz=0.0
                  QaV=V%Qa(i,j)
                  do ltime=1,int(solar_avg*60)
                     HoD=(12.0+G%lon(i,j)/15.0-dble(ltime)/60.0+solar_avg*dble(nsolar-1))*15.0
                     cosZ=sin(delta*pp)*sin(G%lat(i,j)*pp)+
     *                  cos(delta*pp)*cos(G%lat(i,j)*pp)*cos(HoD*pp)
                     if(cosZ.lt.0.0) cosZ=0.0
                     Qz=Qz+1376.0*cosZ**2/((cosZ+2.7)*QaV/(0.622+0.387*QaV)+
     *                  1.085*cosZ+0.1)
                  enddo
                  if(cld_sw_vers.eq.1)then
                     sw=(1.0-0.6*V%cld(i,j)**3)*Qz/(solar_avg*60.0)
                  elseif(cld_sw_vers.eq.2)then
                     betta=90.0-abs(G%lat(i,j)-delta)
                     if(betta.lt.0.0)betta=0.0
                     shading=1.0-0.62*V%cld(i,j)+0.0019*betta
                     if(shading.gt.1.0)shading=1.0
                     sw=shading*Qz/(solar_avg*60.0)
                  endif
                  V%Fsw(i,j)=sw
                  if(G%lon(i,j).gt.108.78.and.G%lon(i,j).lt.142.44.and.
     *               G%lat(i,j).gt.68.53.and.G%lat(i,j).lt.77.5)then
                     SS=SS+sw
                     nSS=nSS+1
                  endif
               endif
            enddo
         enddo
         if(nSS.gt.0)then
c            Print *,'Atm Solar: Day',DoY
c            Print *,'Atm Solar: Time',solar_avg*dble(nsolar-1)
c            Print *,'Atm Solar: Atm Time',dble(S%sec)/3600.0,' next',24.0*(S%next-int(S%next))
            Print '(a31,i3,a5)','Atm Solar: Mean Radiation Flux ',
     *         int(SS/nSS+0.5),' W/m2'
         endif
      endif

c      Write(nu_diag,*)'Pole Ta: ', V%Tair(88,136)-273.16,tdata(88,136)-273.16
c      Write(nu_diag,*)'Pole Qa: ', V%Qa(88,136)
c      Write(nu_diag,*)'Pole SW: ', V%swidr(88,136)/0.31
c      Write(nu_diag,*)'Eq.  Ta: ', V%Tair(82,21)-273.16,tdata(82,21)-273.16
c      Write(nu_diag,*)'Eq.  Qa: ', V%Qa(82,21)
c      Write(nu_diag,*)'Eq.  SW: ', V%swidr(82,21)/0.31
c      if(n.ne.nold.or.useciaf)then
c          Write(nu_diag,*)'Pole   Ta: ', V%Tair(88,136)-273.16
c          Write(nu_diag,*)'Pole   Th: ', V%potT(88,136)-273.16
c          Write(nu_diag,*)'Pole   Qa: ', V%Qa(88,136)
c          Write(nu_diag,*)'Pole  Rho: ', V%rhoa(88,136)
c          Write(nu_diag,*)'Pole Snow: ', V%Fsnow(88,136)
c         Write(nu_diag,*)'Pole Test'
c         Write(nu_diag,*)'zlvl'
c         Write(nu_diag,*)(V%zlvl(i,136),i=87,90)
c         Write(nu_diag,*)'uatm'
c         Write(nu_diag,*)(V%uatm(i,136),i=87,90)
c         Write(nu_diag,*)'vatm'
c         Write(nu_diag,*)(V%vatm(i,136),i=87,90)
c         Write(nu_diag,*)'potT'
c         Write(nu_diag,*)(V%potT(i,136),i=87,90)
c         Write(nu_diag,*)'Tair'
c         Write(nu_diag,*)(V%Tair(i,136),i=87,90)
c         Write(nu_diag,*)'Qa'
c         Write(nu_diag,*)(V%Qa(i,136),i=87,90)
c         Write(nu_diag,*)'rhoa'
c         Write(nu_diag,*)(V%rhoa(i,136),i=87,90)
c         Write(nu_diag,*)'sw'
c         Write(nu_diag,*)(V%Fsw(i,136),i=87,90)
c         Write(nu_diag,*)'Fsnow'
c         Write(nu_diag,*)(V%Fsnow(i,136),i=87,90)
c         Write(nu_diag,*)'Frain'
c         Write(nu_diag,*)(V%Frain(i,136),i=87,90)
c         Write(nu_diag,*)'cld'
c         Write(nu_diag,*)(V%cld(i,136),i=87,90)
c      endif
      
      if(n.ne.nold)then
         if(.false.)then
         fname='../Output/cpl/atmDDD.bin'
         ii=index(fname,'DDD');
         write(fname(ii:ii+2),'(i3.3)') n

         open(23,file=fname,form='unformatted',access='direct',
     *      recl=8*imt*jmt)
         Write(23,rec=1) V%zlvl
         Write(23,rec=2) pdata
         Write(23,rec=3) V%uatm
         Write(23,rec=4) V%vatm
         Write(23,rec=5) V%potT
         Write(23,rec=6) V%Tair
         Write(23,rec=7) V%Qa
         Write(23,rec=8) V%rhoa
         Write(23,rec=9) V%Fsw
         Write(23,rec=10) V%Flw
         Write(23,rec=11) V%Fsnow
         Write(23,rec=12) V%Frain
         close(23)
         endif
      endif
      if(S%idate.ne.idate0)then
         idate0=S%idate
         if(writemode.eq.'w')then
         mh=imt
         nh=jmt
         rpath=trim(datapath)//trim(respath)
         Write(nu_diag,*)'Write data to '//trim(rpath)
         if(opennc)then
            call open_w2w(trim(rpath)//'zlvl0.nc','zlvl',mh,nh)
            call open_w2w(trim(rpath)//'uatm0.nc','uatm',mh,nh)
            call open_w2w(trim(rpath)//'vatm0.nc','vatm',mh,nh)
            call open_w2w(trim(rpath)//'potT0.nc','potT',mh,nh)
            call open_w2w(trim(rpath)//'Tair0.nc','Tair',mh,nh)
            call open_w2w(trim(rpath)//'rhoa0.nc','rhoa',mh,nh)
            call open_w2w(trim(rpath)//'Qa0.nc','Qa',mh,nh)
            call open_w2w(trim(rpath)//'Fsw0.nc','Fsw',mh,nh)
            call open_w2w(trim(rpath)//'cld0.nc','cld',mh,nh)
            call open_w2w(trim(rpath)//'Fsnow0.nc','Fsnow',mh,nh)
            call open_w2w(trim(rpath)//'Frain0.nc','Frain',mh,nh)
            opennc=.false.
            nrec=1
         endif
         call write_w(trim(rpath)//'zlvl0.nc',nrec,idate0,'zlvl',V%zlvl,G%mask,mh,nh)
         call write_w(trim(rpath)//'uatm0.nc',nrec,idate0,'uatm',V%uatm,G%mask,mh,nh)
         call write_w(trim(rpath)//'vatm0.nc',nrec,idate0,'vatm',V%vatm,G%mask,mh,nh)
         call write_w(trim(rpath)//'potT0.nc',nrec,idate0,'potT',V%potT,G%mask,mh,nh)
         call write_w(trim(rpath)//'Tair0.nc',nrec,idate0,'Tair',V%Tair,G%mask,mh,nh)
         call write_w(trim(rpath)//'rhoa0.nc',nrec,idate0,'rhoa',V%rhoa,G%mask,mh,nh)
         call write_w(trim(rpath)//'Qa0.nc',nrec,idate0,'Qa',V%Qa,G%mask,mh,nh)
         call write_w(trim(rpath)//'Fsw0.nc',nrec,idate0,'Fsw',V%Fsw,G%mask,mh,nh)
         call write_w(trim(rpath)//'Fsnow0.nc',nrec,idate0,'Fsnow',V%Fsnow,G%mask,mh,nh)
         call write_w(trim(rpath)//'Frain0.nc',nrec,idate0,'Frain',V%Frain,G%mask,mh,nh)
         call write_w(trim(rpath)//'cld0.nc',nrec,idate0,'cld',V%cld,G%mask,mh,nh)
         nrec=nrec+1
         endif
      endif
      
      nold=n
      nsolar_old=nsolar_tot
      kmold=kmo

c        Here is that line that is to be commented to make a normal run
c      S%idate=idate_old
      
      else
         call read_atm_data(S%idate,V,G,S,aomip_rules,cld_sw_vers)
      endif

      return
      end

      logical function isnan8(r8)
      real*8 r8,q8
      integer*4 i8(2)
      equivalence (i8,q8)
      q8=r8
      isnan8=i8(2).eq.Z'fff80000'
      return
      end

      logical function isnan4(r4)
      real*4 r4,q4
      integer*4 i4
      equivalence (i4,q4)
      q4=r4
      isnan4=i4.eq.Z'fff80000'
      return
      end
