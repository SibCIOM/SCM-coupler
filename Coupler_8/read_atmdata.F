!     $Id: read_atmdata.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $

      Subroutine read_atm_data(idate,V,G,S,aomip_rules,cldv)
      use struct
      Type (state) S
      Type (grid) G
      Type (from_atm) V
      Type (from_atm) dV
      integer idate

      real*8, allocatable :: lon(:),lat(:),pdat(:,:),rmask(:,:),
     *   upolar(:,:),vpolar(:,:),udat(:,:),vdat(:,:)
      real*8, allocatable :: rbox00(:,:),rbox10(:,:),
     *   rbox01(:,:),rbox11(:,:)
      integer, allocatable :: i0box(:,:),j0box(:,:), 
     *   i1box(:,:),j1box(:,:), nodetype(:,:)
      real*8, allocatable, target, dimension(:,:) :: 
     *         dVcld,dVFlw,dVFsw,dVFrain,dVFsnow,
     *         dVpotT,dVTair,dVQa,dVrhoa,dVuatm,dVvatm,dVzlvl
      integer cldv
      logical aomip_rules,firsttime
      character*7 frmt
      integer
     *        id_uas
     *,       id_vas
     *,       id_tas
     *,       id_rsds
     *,       id_rlds
     *,       id_psl
     *,       id_pr
     *,       id_huss
     *,       id_clt
     *,       nrec0
     *,       nrec
     *,       nold
      integer
     *        id_zlvl
     *,       id_uatm
     *,       id_vatm
     *,       id_potT
     *,       id_Tair
     *,       id_rhoa
     *,       id_Qa
     *,       id_Fsw
     *,       id_Flw
     *,       id_Fsnow
     *,       id_Frain
     *,       id_cld
     
      integer monlen(12),nyr0,nyr1,ly,nyear0
      logical use_mask,use_huss,use_increment
      character*80 datapath,rpath,fname,respath
      character*1 writemode
      common/progn_run/datapath,respath,use_mask,writemode

      save

      data firsttime/.true./,id_uas/0/,id_vas/0/,id_tas/0/
     *,id_rsds/0/,id_rlds/0/,id_psl/0/,id_pr/0/,id_huss/0/
     *,id_clt/0/,nrec0/0/,nrec/0/,nold/0/,nfi0/-1/,nfi/0/,idate0/0/
     *,nyear0/0/
c     *datapath/'./../../MassStore/UKMO_HadCM3/'/,
c     *datapath/'./../../MassStore/CNRM_CM5/'/,
c     *datapath/'./../../MassStore/GFDL_CM3/'/,
     *monlen/0,31,59,90,120,151,181,212,243,273,304,334/

      use_increment = .true.

      nyear=int(idate/10000)
      month=int((idate-nyear*10000)/100)
      nday=idate-nyear*10000-month*100
c      use_mask = .false.
c      writemode = 'w'    ! 'n' - read source and interpolate,
c                         ! 'w' - read source, interpolate and write workfile
c                         ! 'r' - just read workfile
      rpath=trim(datapath)//trim(respath)
      mh=S%imt
      nh=S%jmt
      pp=atan2(1.0,1.0)*4.0/180.0

c      nfi = int((1.0*nyear-2006.0)/5.0)
c      nfi = int((1.0*nyear-2006.0)/10.0)
      if(nyear.ne.nyear0)then
      nfi=0
      Open(1002,file=trim(datapath)//'year_content',form='formatted')
      Write(*,*)'Reading a file "'//trim(datapath)//'year_content"'
      do while(.true.)
         Read(1002,*,err=1,end=1)nfir,nyr0,nyr1,ly
c         Write(*,*)nfir,nyr0,nyr1,ly
         if(nyr0.lt.3000)then                ! the range of years is used
            if(nyear.ge.nyr0.and.nyear.le.nyr1)then
               nfi=nfir
               goto 1
            endif
         else                                ! the range of dates is used
            if(idate.ge.nyr0.and.idate.le.nyr1)then
               nfi=nfir
               goto 1
            endif
         endif
      enddo
    1 Close(1002)
      if(nfi.eq.0)Stop 'Cannot read atm files for a specified year'
      nyear0=nyear
      nfi=nfi-1
      else
      nfi=nfi0
      endif

      if(nfi0.ne.nfi)then

         Open(1001,file='nc_op_cl',form='formatted',access='append')
         Write(1001,'(a6,i2,1h/,i2,1h/,i4)')'Date: ',nday,month,nyear

      endif

      if(firsttime)then
         Open(1002,file=trim(datapath)//'file_content',form='formatted')
         nfir=0
         use_huss=.false.
         do while(.true.)
            Read(1002,*,err=2,end=2)fname
            if(trim(fname).eq.'sftlf')nfir=1
            if(trim(fname).eq.'huss')use_huss=.true.
         enddo
    2    Close(1002)
         if(nfir.eq.0)use_mask=.false.
         firsttime=.false.
         if(use_increment)then
            allocate(dVcld(mh,nh),dVFlw(mh,nh),dVFsw(mh,nh),dVFrain(mh,nh),dVFsnow(mh,nh),
     *         dVpotT(mh,nh),dVTair(mh,nh),dVQa(mh,nh),dVrhoa(mh,nh),
     *         dVuatm(mh,nh),dVvatm(mh,nh),dVzlvl(mh,nh))
            dVcld   = 0.0
            dVFlw   = 0.0
            dVFsw   = 0.0
            dVFrain = 0.0
            dVFsnow = 0.0
            dVpotT  = 0.0
            dVTair  = 0.0
            dVQa    = 0.0
            dVrhoa  = 0.0
            dVuatm  = 0.0
            dVvatm  = 0.0
            dVzlvl  = 0.0

            dV%cld   => dVcld
            dV%Flw   => dVFlw
            dV%Fsw   => dVFsw
            dV%Frain => dVFrain
            dV%Fsnow => dVFsnow
            dV%potT  => dVpotT
            dV%Tair  => dVTair
            dV%Qa    => dVQa
            dV%rhoa  => dVrhoa
            dV%uatm  => dVuatm
            dV%vatm  => dVvatm
            dV%zlvl  => dVzlvl

            call open_nc(id_zlvl, trim(rpath)//'zlvl.nc')
            call open_nc(id_uatm, trim(rpath)//'uatm.nc')
            call open_nc(id_vatm, trim(rpath)//'vatm.nc')
            call open_nc(id_potT, trim(rpath)//'potT.nc')
            call open_nc(id_Tair, trim(rpath)//'Tair.nc')
            call open_nc(id_rhoa, trim(rpath)//'rhoa.nc')
            call open_nc(id_Qa,   trim(rpath)//'Qa.nc')
            call open_nc(id_Fsw,  trim(rpath)//'Fsw.nc')
            call open_nc(id_Flw,  trim(rpath)//'Flw.nc')
            call open_nc(id_Fsnow,trim(rpath)//'Fsnow.nc')
            call open_nc(id_Frain,trim(rpath)//'Frain.nc')
            call open_nc(id_cld,  trim(rpath)//'cld.nc')

            call read_nc(id_zlvl, 1,'zlvl', dV%zlvl)
            call read_nc(id_uatm, 1,'uatm', dV%uatm)
            call read_nc(id_vatm, 1,'vatm', dV%vatm)
            call read_nc(id_potT, 1,'potT', dV%potT)
            call read_nc(id_Tair, 1,'Tair', dV%Tair)
            call read_nc(id_rhoa, 1,'rhoa', dV%rhoa)
            call read_nc(id_Qa,   1,'Qa',   dV%Qa)
            call read_nc(id_Fsw,  1,'Fsw',  dV%Fsw)
            call read_nc(id_Flw,  1,'Flw',  dV%Flw)
            call read_nc(id_Fsnow,1,'Fsnow',dV%Fsnow)
            call read_nc(id_Frain,1,'Frain',dV%Frain)
            call read_nc(id_cld,  1,'cld',  dV%cld)

            call close_nc(id_zlvl)
            call close_nc(id_uatm)
            call close_nc(id_vatm)
            call close_nc(id_potT)
            call close_nc(id_Tair)
            call close_nc(id_rhoa)
            call close_nc(id_Qa)
            call close_nc(id_Fsw)
            call close_nc(id_Flw)
            call close_nc(id_Fsnow)
            call close_nc(id_Frain)
            call close_nc(id_cld)
         endif
      endif

      if(nfi0.ne.nfi)then

         if(nfi0.ge.0)then

            if (writemode.ne.'r')then

               call close_nc(id_uas)

               call close_nc(id_vas)

               call close_nc(id_tas)

               call close_nc(id_rsds)

               call close_nc(id_rlds)

               call close_nc(id_psl)

               call close_nc(id_pr)

               if(use_huss)call close_nc(id_huss)

               call close_nc(id_clt)

            endif

            if (writemode.eq.'r')then

               call close_nc(id_zlvl)

               call close_nc(id_uatm)

               call close_nc(id_vatm)

               call close_nc(id_potT)

               call close_nc(id_Tair)

               call close_nc(id_rhoa)

               call close_nc(id_Qa)

               call close_nc(id_Fsw)

               call close_nc(id_Flw)

               call close_nc(id_Fsnow)

               call close_nc(id_Frain)

               call close_nc(id_cld)

            endif


         else

            if (writemode.ne.'r')then

               call getsize_nc(trim(datapath)//'clt1.nc',imt,jmt)

               allocate(lon(imt),lat(jmt),pdat(imt,jmt),rmask(imt,jmt),
     *            udat(imt,jmt),vdat(imt,jmt))

               call getlonlat_nc(trim(datapath)//'clt1.nc',lon,lat)

               rmask=1.0
               if(use_mask)then
                  call getmask_nc(trim(datapath)//'sftlf.nc','sftlf',
     *               imt,jmt,rmask)
               endif

               allocate(rbox00(mh,nh),rbox10(mh,nh),
     *            rbox01(mh,nh),rbox11(mh,nh),i0box(mh,nh),
     *            j0box(mh,nh),i1box(mh,nh),j1box(mh,nh),
     *            nodetype(mh,nh),upolar(mh,nh),vpolar(mh,nh))

               call prepare_coeffs( rbox00, rbox10, rbox01, rbox11, 
     *            i0box, j0box,
     *            i1box, j1box, nodetype,
     *            imt, jmt, lon, lat, rmask,
     *            mh,  nh,  G%lon,  G%lat, G%mask)

            endif

         endif

         if(writemode.ne.'r')then

            frmt='(a3,i1)'
            if(nfi+1.gt.9)frmt='(a3,i2)'
            
            Write(fname,frmt)'uas',nfi+1
            call open_nc(id_uas,trim(datapath)//trim(fname)//'.nc')

            Write(fname,frmt)'vas',nfi+1
            call open_nc(id_vas,trim(datapath)//trim(fname)//'.nc')

            Write(fname,frmt)'tas',nfi+1
            call open_nc(id_tas,trim(datapath)//trim(fname)//'.nc')

            frmt='(a4,i1)'
            if(nfi+1.gt.9)frmt='(a4,i2)'
            
            Write(fname,frmt)'rsds',nfi+1
            call open_nc(id_rsds,trim(datapath)//trim(fname)//'.nc')

            Write(fname,frmt)'rlds',nfi+1
            call open_nc(id_rlds,trim(datapath)//trim(fname)//'.nc')

            frmt='(a3,i1)'
            if(nfi+1.gt.9)frmt='(a3,i2)'
            
            Write(fname,frmt)'psl',nfi+1
            call open_nc(id_psl,trim(datapath)//trim(fname)//'.nc')

            frmt='(a2,i1)'
            if(nfi+1.gt.9)frmt='(a2,i2)'
            
            Write(fname,frmt)'pr',nfi+1
            call open_nc(id_pr,trim(datapath)//trim(fname)//'.nc')

            if(use_huss)then

               frmt='(a4,i1)'
               if(nfi+1.gt.9)frmt='(a4,i2)'
            
               Write(fname,frmt)'huss',nfi+1
               call open_nc(id_huss,trim(datapath)//trim(fname)//'.nc')

            endif

            frmt='(a3,i1)'
            if(nfi+1.gt.9)frmt='(a3,i2)'
            
            Write(fname,frmt)'clt',nfi+1
            call open_nc(id_clt,trim(datapath)//trim(fname)//'.nc')

         endif

         if(writemode.eq.'w')then

            frmt='(a4,i1)'
            if(nfi+1.gt.9)frmt='(a4,i2)'
            
            Write(fname,frmt)'zlvl',nfi+1
            call open_w2w(trim(rpath)//trim(fname)//'.nc',
     *         'zlvl',mh,nh)

            Write(fname,frmt)'uatm',nfi+1
            call open_w2w(trim(rpath)//trim(fname)//'.nc',
     *         'uatm',mh,nh)

            Write(fname,frmt)'vatm',nfi+1
            call open_w2w(trim(rpath)//trim(fname)//'.nc',
     *         'vatm',mh,nh)

            Write(fname,frmt)'potT',nfi+1
            call open_w2w(trim(rpath)//trim(fname)//'.nc',
     *         'potT',mh,nh)

            Write(fname,frmt)'Tair',nfi+1
            call open_w2w(trim(rpath)//trim(fname)//'.nc',
     *         'Tair',mh,nh)

            Write(fname,frmt)'rhoa',nfi+1
            call open_w2w(trim(rpath)//trim(fname)//'.nc',
     *         'rhoa',mh,nh)

            frmt='(a2,i1)'
            if(nfi+1.gt.9)frmt='(a2,i2)'
            
            Write(fname,frmt)'Qa',nfi+1
            call open_w2w(trim(rpath)//trim(fname)//'.nc',
     *         'Qa',mh,nh)

            frmt='(a3,i1)'
            if(nfi+1.gt.9)frmt='(a3,i2)'
            
            Write(fname,frmt)'Fsw',nfi+1
            call open_w2w(trim(rpath)//trim(fname)//'.nc',
     *         'Fsw',mh,nh)

            Write(fname,frmt)'Flw',nfi+1
            call open_w2w(trim(rpath)//trim(fname)//'.nc',
     *         'Flw',mh,nh)

            Write(fname,frmt)'cld',nfi+1
            call open_w2w(trim(rpath)//trim(fname)//'.nc',
     *         'cld',mh,nh)

            frmt='(a5,i1)'
            if(nfi+1.gt.9)frmt='(a5,i2)'
            
            Write(fname,frmt)'Fsnow',nfi+1
            call open_w2w(trim(rpath)//trim(fname)//'.nc',
     *         'Fsnow',mh,nh)

            Write(fname,frmt)'Frain',nfi+1
            call open_w2w(trim(rpath)//trim(fname)//'.nc',
     *         'Frain',mh,nh)

         endif

         if(writemode.eq.'r')then

            frmt='(a2,i1)'
            if(nfi+1.gt.9)frmt='(a2,i2)'
            
            Write(fname,frmt)'zlvl',nfi+1
            call open_nc(id_zlvl,trim(datapath)//trim(fname)//'.nc')

            Write(fname,frmt)'uatm',nfi+1
            call open_nc(id_uatm,trim(datapath)//trim(fname)//'.nc')

            Write(fname,frmt)'vatm',nfi+1
            call open_nc(id_vatm,trim(datapath)//trim(fname)//'.nc')

            Write(fname,frmt)'potT',nfi+1
            call open_nc(id_potT,trim(datapath)//trim(fname)//'.nc')

            Write(fname,frmt)'Tair',nfi+1
            call open_nc(id_Tair,trim(datapath)//trim(fname)//'.nc')

            Write(fname,frmt)'rhoa',nfi+1
            call open_nc(id_rhoa,trim(datapath)//trim(fname)//'.nc')

            frmt='(a2,i1)'
            if(nfi+1.gt.9)frmt='(a2,i2)'
            
            Write(fname,frmt)'Qa',nfi+1
            call open_nc(id_Qa,trim(datapath)//trim(fname)//'.nc')

            frmt='(a3,i1)'
            if(nfi+1.gt.9)frmt='(a3,i2)'
            
            Write(fname,frmt)'Fsw',nfi+1
            call open_nc(id_Fsw,trim(datapath)//trim(fname)//'.nc')

            Write(fname,frmt)'Flw',nfi+1
            call open_nc(id_Flw,trim(datapath)//trim(fname)//'.nc')

            frmt='(a5,i1)'
            if(nfi+1.gt.9)frmt='(a5,i2)'
            
            Write(fname,frmt)'Fsnow',nfi+1
            call open_nc(id_Fsnow,trim(datapath)//trim(fname)//'.nc')

            Write(fname,frmt)'Frain',nfi+1
            call open_nc(id_Frain,trim(datapath)//trim(fname)//'.nc')

            frmt='(a3,i1)'
            if(nfi+1.gt.9)frmt='(a3,i2)'
            
            Write(fname,frmt)'cld',nfi+1
            call open_nc(id_cld,trim(datapath)//trim(fname)//'.nc')

         endif

         Close(1001)

         nfi0=nfi

         nrec = 0

      endif

      if(idate.ne.idate0)then

c         idate0=idate

         nrec=nrec+1
         if( ly.eq.1  .and.
     *       4*int(0.25*nyear).eq.nyear  .and.
     *       month.eq.3  .and.
     *       nday.eq.1  )nrec=nrec+1

         if (writemode.ne.'r') then

            V%zlvl = 10.0

            call read_nc(id_uas,nrec,'uas',udat)
c            pdat = pdat * 1.414213562
            call projection(imt,jmt,udat,mh,nh,G%mask,V%uatm, 
     *         rbox00,rbox10,rbox01,rbox11,i0box,j0box,i1box,j1box,nodetype)

            call read_nc(id_vas,nrec,'vas',vdat)
c            pdat = pdat * 1.414213562
            call projection(imt,jmt,vdat,mh,nh,G%mask,V%vatm, 
     *         rbox00,rbox10,rbox01,rbox11,i0box,j0box,i1box,j1box,nodetype)

c-------------------------------------------------------------------------------
c     The following is to fix polar vector problem
c-------------------------------------------------------------------------------
            do i=1,imt
               cs=cos(lon(i)*pp)
               sn=sin(lon(i)*pp)
               do j=1,jmt
                  uwrk=udat(i,j)*cs-vdat(i,j)*sn
                  vdat(i,j)=udat(i,j)*sn+vdat(i,j)*cs
                  udat(i,j)=uwrk
               enddo
            enddo
            call projection(imt,jmt,udat,mh,nh,G%mask,upolar, 
     *         rbox00,rbox10,rbox01,rbox11,i0box,j0box,i1box,j1box,nodetype)
            call projection(imt,jmt,vdat,mh,nh,G%mask,vpolar, 
     *         rbox00,rbox10,rbox01,rbox11,i0box,j0box,i1box,j1box,nodetype)
            do i=1,mh
               do j=1,nh
                  if(G%mask(i,j).ne.0.and.abs(G%lat(i,j)).gt.60.0)then
                     cs=cos(G%lon(i,j)*pp)
                     sn=sin(G%lon(i,j)*pp)
                     V%uatm(i,j)=upolar(i,j)*cs+vpolar(i,j)*sn
                     V%vatm(i,j)=-upolar(i,j)*sn+vpolar(i,j)*cs
                  endif
               enddo
            enddo
c-------------------------------------------------------------------------------

            call rotate_uv(mh,nh,G%lon,G%lat,V%uatm,V%vatm,G%mask)

            call read_nc(id_psl,nrec,'psl',pdat)
            call projection(imt,jmt,pdat,mh,nh,G%mask,V%potT, 
     *         rbox00,rbox10,rbox01,rbox11,i0box,j0box,i1box,j1box,nodetype)

            call read_nc(id_tas,nrec,'tas',pdat)
            call projection(imt,jmt,pdat,mh,nh,G%mask,V%Tair, 
     *         rbox00,rbox10,rbox01,rbox11,i0box,j0box,i1box,j1box,nodetype)

            if(use_huss)then
               call read_nc(id_huss,nrec,'huss',pdat)
               call projection(imt,jmt,pdat,mh,nh,G%mask,V%Qa, 
     *            rbox00,rbox10,rbox01,rbox11,i0box,j0box,i1box,j1box,nodetype)
            endif

            do i=1,mh
               do j=1,nh
                  if(G%mask(i,j).ne.0)then
                     V%rhoa(i,j)=V%potT(i,j)/(287.04*V%Tair(i,j))
                     V%potT(i,j)=V%Tair(i,j)*(1.0e+5/V%potT(i,j))**0.29
                     if(.not.use_huss)then
                        Tk=V%Tair(i,j)-273.16
                        esat = 0.98*exp((0.7859+0.03477*Tk)/
     *                     (1.0+0.00412*Tk)*2.302585093)
                        qsat = 0.622*esat/(1000.0-0.378*esat)
                        V%Qa(i,j)=0.9*qsat
                     endif
                  endif
               enddo
            enddo

            if(.not.aomip_rules)then
         
               call read_nc(id_rsds,nrec,'rsds',pdat)
               call projection(imt,jmt,pdat,mh,nh,G%mask,V%Fsw, 
     *            rbox00,rbox10,rbox01,rbox11,i0box,j0box,i1box,j1box,nodetype)

               call read_nc(id_rlds,nrec,'rlds',pdat)
               call projection(imt,jmt,pdat,mh,nh,G%mask,V%Flw, 
     *            rbox00,rbox10,rbox01,rbox11,i0box,j0box,i1box,j1box,nodetype)
            
            endif

            call read_nc(id_pr,nrec,'pr',pdat)
            call projection(imt,jmt,pdat,mh,nh,G%mask,V%Fsnow, 
     *         rbox00,rbox10,rbox01,rbox11,i0box,j0box,i1box,j1box,nodetype)

            do i=1,mh
               do j=1,nh
                  if(G%mask(i,j).ne.0)then
                     V%Frain(i,j)=0.0
                     if(V%Tair(i,j).gt.273.15)then
                        V%Frain(i,j)=V%Fsnow(i,j)
                        V%Fsnow(i,j)=0.0
                     endif
                  endif
               enddo
            enddo

            call read_nc(id_clt,nrec,'clt',pdat)
            pdat=pdat/100.0
            call projection(imt,jmt,pdat,mh,nh,G%mask,V%cld, 
     *         rbox00,rbox10,rbox01,rbox11,i0box,j0box,i1box,j1box,nodetype)

            if(use_increment)then
               DoY=dble(n-1)*360/365
               delta=pp*23.44*cos((172-DoY)*pp)
               do i=1,mh
               do j=1,nh
                  if(G%mask(i,j).ne.0)then
                     V%uatm(i,j)=V%uatm(i,j)+dV%uatm(i,j)
                     V%vatm(i,j)=V%vatm(i,j)+dV%vatm(i,j)

                     V%Tair(i,j)=V%Tair(i,j)+dV%Tair(i,j)

                     V%rhoa(i,j)=V%rhoa(i,j)+dV%rhoa(i,j)

                     V%potT(i,j)=V%potT(i,j)+dV%potT(i,j)

                     V%Qa(i,j)=V%Qa(i,j)+dV%Qa(i,j)
                     if(V%Qa(i,j).lt.0.0) V%Qa(i,j)=0.0
                        Tk=V%Tair(i,j)-273.16
                        esat = 0.98*exp((0.7859+0.03477*Tk)/
     *                     (1.0+0.00412*Tk)*2.302585093)
                        qsat = 0.622*esat/(1000.0-0.378*esat)
                     if(V%Qa(i,j).gt.qsat) V%Qa(i,j)=qsat

                     V%cld(i,j)=V%cld(i,j)+dV%cld(i,j)
                     if(V%cld(i,j).lt.0.0) V%cld(i,j)=0.0
                     if(V%cld(i,j).gt.1.0) V%cld(i,j)=1.0

                     if(V%Fsw(i,j).ne.0.0) V%Fsw(i,j)=V%Fsw(i,j)+dV%Fsw(i,j)
                     if(V%Fsw(i,j).lt.0.0) V%Fsw(i,j)=0.0
                     Qz=0.0
                     do ltime=1,24*60
                        HoD=(12.0-dble(ltime)/60.0)*pp*15.0
                        cosZ=sin(delta)*sin(G%lat(i,j)*pp)+
     *                     cos(delta)*cos(G%lat(i,j)*pp)*cos(HoD)
                        if(cosZ.lt.0.0) cosZ=0.0
                        Qz=Qz+1376.0*cosZ**2/(1.085*cosZ+0.1)
                     enddo
                     if(cldv.eq.2)then
                       betta=90.0-abs(G%lat(i,j)*pp-delta)/pp
                       if(betta.lt.0.0)betta=0.0
                       shading=1.0+0.0019*betta
                       if(shading.gt.1.0)shading=1.0
                       Qz=shading*Qz
                     endif
                     sw=Qz/dble(24*60)
                     if(V%Fsw(i,j).gt. sw) V%Fsw(i,j)=sw

                     if(V%Frain(i,j).ne.0.0) V%Frain(i,j)=V%Frain(i,j)+dV%Frain(i,j)
                     if(V%Fsnow(i,j).ne.0.0) V%Fsnow(i,j)=V%Fsnow(i,j)+dV%Fsnow(i,j)
                     prec=V%Fsnow(i,j)+V%Frain(i,j)
                     if(prec.lt.0.0) then
                        V%Fsnow(i,j)=0.0
                        V%Frain(i,j)=0.0
                     elseif(V%Fsnow(i,j).lt.0 .or. V%Frain(i,j).lt.0.0) then
                        if(V%Tair(i,j).gt.273.15)then
                           V%Fsnow(i,j)=0.0
                           V%Frain(i,j)=prec
                        else
                           V%Fsnow(i,j)=prec
                           V%Frain(i,j)=0.0
                        endif
                     endif
                     V%Flw(i,j)=V%Flw(i,j)+dV%Flw(i,j)
                     if(V%Flw(i,j).lt. 0.0) V%Flw(i,j)=0.0
                  endif
               enddo
               enddo
            endif
         endif

         Print *,'X wind:            ',minval(V%uatm,MASK=G%mask.ne.0),
     *                                 maxval(V%uatm,MASK=G%mask.ne.0)
         Print *,'Y wind:            ',minval(V%vatm,MASK=G%mask.ne.0),
     *                                 maxval(V%vatm,MASK=G%mask.ne.0)
         Print *,'potential temp.:   ',minval(V%potT,MASK=G%mask.ne.0),
     *                                 maxval(V%potT,MASK=G%mask.ne.0)
         Print *,'air temperature:   ',minval(V%Tair,MASK=G%mask.ne.0),
     *                                 maxval(V%Tair,MASK=G%mask.ne.0)
         Print *,'air humidity:      ',minval(V%Qa,MASK=G%mask.ne.0),
     *                                 maxval(V%Qa,MASK=G%mask.ne.0)
         Print *,'solar radiation:   ',minval(V%Fsw,MASK=G%mask.ne.0),
     *                                 maxval(V%Fsw,MASK=G%mask.ne.0)
         Print *,'lw downward rad:   ',minval(V%Flw,MASK=G%mask.ne.0),
     *                                 maxval(V%Flw,MASK=G%mask.ne.0)
         Print *,'precipitation:     ',minval(V%Frain+V%Fsnow,MASK=G%mask.ne.0),
     *                                 maxval(V%Frain+V%Fsnow,MASK=G%mask.ne.0)
         Print *,'cloudiness:        ',minval(V%cld,MASK=G%mask.ne.0),
     *                                 maxval(V%cld,MASK=G%mask.ne.0)

c         if(writemode.eq.'w')then
c
c         endif

      endif
      
      
      n=monlen(month)+nday
      if(nold.ne.n.and.aomip_rules)then
         nold=n
         DoY=dble(n-1)*360/365
         delta=pp*23.44*cos((172-DoY)*pp)
         do j=1,nh
            do i=1,mh
               if(G%mask(i,j).ne.0)then

c   AOMIP Short Wave radiative flux

                  QaV=V%Qa(i,j)
                  if(QaV.le.0.0)QaV=0.0

                  Qz=0.0

                  do ltime=1,24*60
                     HoD=(12.0-dble(ltime)/60.0)*pp*15.0
                     cosZ=sin(delta)*sin(G%lat(i,j)*pp)+
     *                  cos(delta)*cos(G%lat(i,j)*pp)*cos(HoD)
                     if(cosZ.lt.0.0) cosZ=0.0
                     Qz=Qz+1376.0*cosZ**2/((cosZ+2.7)*QaV/(0.622+0.387*QaV)+
     *                  1.085*cosZ+0.1)
                  enddo

                  if(cldv.eq.1)then
                     sw=(1.0-0.6*V%cld(i,j)**3)*Qz/dble(24*60)
                  elseif(cldv.eq.2)then
                     betta=90.0-abs(G%lat(i,j)*pp-delta)/pp
                     if(betta.lt.0.0)betta=0.0
                     shading=1.0-0.62*V%cld(i,j)+0.0019*betta
                     if(shading.gt.1.0)shading=1.0
                     sw=shading*Qz/dble(24*60)
                  endif
                  V%Fsw(i,j)=sw
c                  Print *,i,j,sw

               endif
            enddo
         enddo
      endif

      if(idate.ne.idate0)then

         idate0=idate

c         nrec=nrec+1

         if (writemode.eq.'w') then
            Write(fname,'(a4,i1)')'zlvl',nfi+1
            call write_w(trim(rpath)//trim(fname)//'.nc',nrec,idate,'zlvl',V%zlvl,G%mask,mh,nh)
            Write(fname,'(a4,i1)')'uatm',nfi+1
            call write_w(trim(rpath)//trim(fname)//'.nc',nrec,idate,'uatm',V%uatm,G%mask,mh,nh)
            Write(fname,'(a4,i1)')'vatm',nfi+1
            call write_w(trim(rpath)//trim(fname)//'.nc',nrec,idate,'vatm',V%vatm,G%mask,mh,nh)
            Write(fname,'(a4,i1)')'potT',nfi+1
            call write_w(trim(rpath)//trim(fname)//'.nc',nrec,idate,'potT',V%potT,G%mask,mh,nh)
            Write(fname,'(a4,i1)')'Tair',nfi+1
            call write_w(trim(rpath)//trim(fname)//'.nc',nrec,idate,'Tair',V%Tair,G%mask,mh,nh)
            Write(fname,'(a4,i1)')'rhoa',nfi+1
            call write_w(trim(rpath)//trim(fname)//'.nc',nrec,idate,'rhoa',V%rhoa,G%mask,mh,nh)
            Write(fname,'(a2,i1)')'Qa',nfi+1
            call write_w(trim(rpath)//trim(fname)//'.nc',nrec,idate,'Qa',V%Qa,G%mask,mh,nh)
            Write(fname,'(a3,i1)')'Fsw',nfi+1
            call write_w(trim(rpath)//trim(fname)//'.nc',nrec,idate,'Fsw',V%Fsw,G%mask,mh,nh)
            Write(fname,'(a3,i1)')'Flw',nfi+1
            call write_w(trim(rpath)//trim(fname)//'.nc',nrec,idate,'Flw',V%Flw,G%mask,mh,nh)
            Write(fname,'(a5,i1)')'Fsnow',nfi+1
            call write_w(trim(rpath)//trim(fname)//'.nc',nrec,idate,'Fsnow',V%Fsnow,G%mask,mh,nh)
            Write(fname,'(a5,i1)')'Frain',nfi+1
            call write_w(trim(rpath)//trim(fname)//'.nc',nrec,idate,'Frain',V%Frain,G%mask,mh,nh)
            Write(fname,'(a3,i1)')'cld',nfi+1
            call write_w(trim(rpath)//trim(fname)//'.nc',nrec,idate,'cld',V%cld,G%mask,mh,nh)
         endif

         if (writemode.eq.'r') then
               call read_nc(id_zlvl,nrec,'zlvl',V%zlvl)
               call read_nc(id_uatm,nrec,'uatm',V%uatm)
               call read_nc(id_vatm,nrec,'vatm',V%vatm)
               call read_nc(id_potT,nrec,'potT',V%potT)
               call read_nc(id_Tair,nrec,'Tair',V%Tair)
               call read_nc(id_rhoa,nrec,'rhoa',V%rhoa)
               call read_nc(id_Qa,nrec,'Qa',V%Qa)
               call read_nc(id_Fsw,nrec,'Fsw',V%Fsw)
               call read_nc(id_Flw,nrec,'Flw',V%Flw)
               call read_nc(id_Fsnow,nrec,'Fsnow',V%Fsnow)
               call read_nc(id_Frain,nrec,'Frain',V%Frain)
               call read_nc(id_cld,nrec,'cld',V%cld)
         endif

      endif

      return
      end

      Subroutine read_nc(nid, nrec, varname, array)

      integer, intent(in) :: nid
      integer, intent(in) :: nrec
      character*(*), intent(in) :: varname
      real*8, intent(out) :: array(*)
      
      include 'netcdf.inc'

      integer ierr, ncid0, nDimensions, nVariables,
     *   nAttributes, unlimitedDimId, leng, dimids(4),
     *   xtype, ndims, natts, vtype, nidold, vid, 
     *   ip, jp, kp

      real*4, allocatable :: r4vals(:,:)
      real*4 add_offset,scale_factor,missing,fillval
      integer start(3),counts(3)

      ncid0 = nid

      ierr=nf_inq_dimid(ncid0,'lon',vid)
      if(ierr.ne.NF_NOERR)then
         ierr=nf_inq_dimid(ncid0,'i',vid)
         if(ierr.ne.NF_NOERR)stop 'dont know how to open lon'
      endif
      ierr = nf_inq_dimlen(ncid0, vid, ip)

      ierr=nf_inq_dimid(ncid0,'lat',vid)
      if(ierr.ne.NF_NOERR)then
         ierr=nf_inq_dimid(ncid0,'j',vid)
         if(ierr.ne.NF_NOERR)stop 'dont know how to open lat'
      endif
      ierr = nf_inq_dimlen(ncid0, vid, jp)

      ierr=nf_inq_dimid(ncid0,'time',vid)
      if(ierr.ne.NF_NOERR)then
         kp=1
      else
         ierr = nf_inq_dimlen(ncid0, vid, kp)
      endif


      ierr=nf_inq_varid(ncid0,trim(varname),vid)
      kk=nrec
      if(kk.lt.1)then
       Write(*,*) 'read_nc: recno: ',kk,' < 1     ncid=',nid
       kk=1
      endif
      if(kk.gt.kp)then
       Write(*,*) 'read_nc: recno: ',kk,' > ',kp,'     ncid=',nid
       kk=kp
      endif

      allocate(r4vals(ip,jp))

      start(1)=1
      start(2)=1
      start(3)=kk
      counts(1)=ip
      counts(2)=jp
      counts(3)=1
      ierr=nf_get_vara_real(ncid0,vid,start,counts,r4vals)

      ierr=nf_get_att_real(ncid0,vid,'missing_value',missing)
      if(ierr.ne.NF_NOERR)missing=1.e+20
      ierr=nf_get_att_real(ncid0,vid,'_FillValue',fillval)
      if(ierr.ne.NF_NOERR)fillval=1.e+20

      do i=1,ip
         do j=1,jp
            if(r4vals(i,j).eq.missing.or.r4vals(i,j).eq.fillval)then
               array(i+ip*(j-1))=1.0e+10
            else
               array(i+ip*(j-1))=dble(r4vals(i,j))
            endif
         enddo
      enddo

      do i=1,ip
         do j=1,jp
            if(array(i+ip*(j-1)).gt.0.9e+10)then
               do k=1,10
                  sums=0.0
                  kns=0
                  do ii=max(i-k,1),min(i+k,ip)
                     do jj=max(j-k,1),min(j+k,jp)
                        if (array(ii+ip*(jj-1)).le.0.9e+10)then
                           sums=sums+array(ii+ip*(jj-1))
                           kns=kns+1
                        endif
                     enddo
                  enddo
                  if(kns.gt.0)exit
               enddo
               if(kns.eq.0)then
                  print *,'ncid=',ncid0,',  i=',i,',  j=',j
                  stop 'this radius in read_nc is small'
               endif
               array(i+ip*(j-1))=sums/kns
            endif
         enddo
      enddo

      deallocate(r4vals)

      Return
      End

      Subroutine getsize_nc(filename, imt, jmt)

      character*(*), intent(in) :: filename
      integer, intent(out) :: imt, jmt
      
      include 'netcdf.inc'

      integer ierr, ncid0, nDimensions, nVariables,
     *   nAttributes, unlimitedDimId, leng, dimids(4),
     *   xtype, ndims, natts, vtype, nidold, vid, 
     *   ip, jp, kp

c         NF_OPEN ! open existing netCDF dataset

      ierr = nf_open(trim(filename), NF_NOWRITE, ncid0)

      if(ierr.ne.NF_NOERR)then
         Write(*,*)'Error: file '//trim(filename)//
     *      ' is not present or damaged'
         Stop
      endif

      ierr=nf_inq_dimid(ncid0,'lon',vid)
      ierr = nf_inq_dimlen(ncid0, vid, imt)

      ierr=nf_inq_dimid(ncid0,'lat',vid)
      ierr = nf_inq_dimlen(ncid0, vid, jmt)

      ierr = nf_close(ncid0)

      Return
      End

      Subroutine getlonlat_nc(filename, lon, lat)

      character*(*), intent(in) :: filename
      real*8, intent(out) :: lon(*), lat(*)
      
      include 'netcdf.inc'

      integer ierr, ncid0, nDimensions, nVariables,
     *   nAttributes, unlimitedDimId, leng, dimids(4),
     *   xtype, ndims, natts, vtype, nidold, vid, 
     *   ip, jp, kp

      integer start(3),counts(3)

c         NF_OPEN ! open existing netCDF dataset

      ierr = nf_open(trim(filename), NF_NOWRITE, ncid0)

      if(ierr.ne.NF_NOERR)then
         Write(*,*) nf_strerror(ierr)
         Write(*,*)'Error: file '//trim(filename)//
     *      ' is not present or damaged'
         Stop
      endif

      ierr=nf_inq_dimid(ncid0,'lon',vid)
      ierr = nf_inq_dimlen(ncid0, vid, ip)

      ierr=nf_inq_varid(ncid0,'lon',vid)
      start(1)=1
      counts(1)=ip
      ierr = nf_get_vara_double(ncid0,vid,start, 
     *            counts,lon)

      ierr=nf_inq_dimid(ncid0,'lat',vid)
      ierr = nf_inq_dimlen(ncid0, vid, jp)
        
      ierr=nf_inq_varid(ncid0,'lat',vid)
      start(1)=1
      counts(1)=jp
      ierr = nf_get_vara_double(ncid0,vid,start, 
     *            counts,lat)

      ierr = nf_close(ncid0)

      Return
      End

      Subroutine getmask_nc(filename, varname, imt, jmt, rmask)

      integer imt,jmt
      character*(*), intent(in) :: filename, varname
      real*8, dimension(imt,jmt),intent(out) :: rmask
      
      include 'netcdf.inc'

      integer ierr, ncid0, nDimensions, nVariables,
     *   nAttributes, unlimitedDimId, leng, dimids(4),
     *   xtype, ndims, natts, vtype, nidold, vid, 
     *   ip, jp, kp

      integer start(3),counts(3)

c         NF_OPEN ! open existing netCDF dataset

      ierr = nf_open(trim(filename), NF_NOWRITE, ncid0)

      if(ierr.ne.NF_NOERR)then
         Write(*,*)'Error: file '//trim(filename)//
     *      ' is not present or damaged'
         Stop
      endif

      ierr=nf_inq_varid(ncid0,trim(varname),vid)
      start(1)=1
      start(2)=1
      counts(1)=imt
      counts(2)=jmt
      ierr = nf_get_vara_real(ncid0,vid,start, 
     *            counts,rmask)

      ierr = nf_close(ncid0)
      
      rmask=(100.0-rmask)/100.0

      Return
      End

      Subroutine open_nc(nid, filename)

      integer, intent(inout) :: nid
      character*(*), intent(in) :: filename
      
      include 'netcdf.inc'

      integer ierr, ncid0, nDimensions, nVariables,
     *   nAttributes, unlimitedDimId, leng, dimids(4),
     *   xtype, ndims, natts, vtype, nidold, vid, 
     *   ip, jp, kp

c         NF_OPEN ! open existing netCDF dataset

      ierr = nf_open(trim(filename), NF_NOWRITE, ncid0)

      if(ierr.ne.NF_NOERR)then
         Write(*,*)'Error: file '//trim(filename)//
     *      ' is not present or damaged'
         Stop
      else
         Write(1001,'(a)')'File '//trim(filename)//' is open'
         Write(1001,'(10x,a7,i8)')'ncid = ',ncid0
      endif

      nid=ncid0

      Return
      End

      Subroutine open_w2w(filename, varname, im, jm)

c      integer, intent(inout) :: nid
      character*(*), intent(in) :: filename, varname
      integer, intent(in) ::  im, jm
      
      include 'netcdf.inc'

      integer ierr, ncid0, nDimensions, nVariables,
     *   nAttributes, unlimitedDimId, leng, dimids(4),
     *   xtype, ndims, natts, vtype, nidold, vid, 
     *   ip, jp, kp
      real*4 , allocatable :: ix(:), iy(:)

c         NF_CREATE ! open new netCDF dataset

      ierr = nf_create(trim(filename), NF_CLOBBER, ncid0)

      if(ierr.ne.NF_NOERR)then
         Write(*,*)'Error: file '//trim(filename)//
     *      ' cannot be created'
         Stop
      else
         Write(1001,'(a)')'File '//trim(filename)//' is created'
         Write(1001,*)'ncid = ',ncid0
      endif

      ierr = nf_def_dim(ncid0, 'lon', im, dimids(1))
      ierr = nf_def_dim(ncid0, 'lat', jm, dimids(2))
      ierr = nf_def_dim(ncid0, 'time', NF_UNLIMITED, dimids(3))
      ierr = nf_def_var(ncid0, 'lon', NF_INT, 1,
     *                         dimids(1), vid)
      ierr = nf_def_var(ncid0, 'lat', NF_INT, 1,
     *                         dimids(2), vid)
      ierr = nf_def_var(ncid0, 'time', NF_INT, 1,
     *                         dimids(3), vid)

      ierr = nf_def_var(ncid0, trim(varname), NF_FLOAT, 3,
     *                         dimids, vid)
      ierr = nf_put_att_real(ncid0, vid, 'missing_value', NF_FLOAT,
     *                         1, 2.0e+10)
      ierr = nf_put_att_real(ncid0, vid, '_FillValue', NF_FLOAT,
     *                         1, NF_FILL_FLOAT)
      ierr = nf_enddef(ncid0)
      ierr = nf_close(ncid0)
      
c      ierr = nf_open(trim(filename), NF_WRITE, ncid0)
c      ierr=nf_inq_varid(ncid0,'lon',vid)
c      allocate(ix(im))
c      Do i=1,im
c         ix(i)=1.0*i
c      Enddo
c      ierr=nf_put_vara_real(ncid0,vid,1,im,ix)

c      ierr=nf_inq_varid(ncid0,'lat',vid)
c      allocate(iy(jm))
c      Do i=1,jm
c         iy(i)=1.0*i
c      Enddo
c      ierr=nf_put_vara_real(ncid0,vid,1,jm,iy)
      
c      deallocate(ix,iy)
      
c      ierr=nf_close(ncid0)

      Return
      End

      Subroutine write_w(filename, nrec, idate, varname, array, mask, mh, nh)

c      integer, intent(inout) :: nid
      character*(*), intent(in) :: filename, varname
      integer, intent(in) ::  nrec, mh, nh
      real*8, intent(in) :: array(mh,nh)
      integer, intent(in) :: mask(mh,nh)
      
      include 'netcdf.inc'

      integer ierr, ncid0, nDimensions, nVariables,
     *   nAttributes, unlimitedDimId, leng, dimids(4),
     *   xtype, ndims, natts, vtype, nidold, vid, 
     *   ip, jp, kp

      real*4, allocatable :: r4vals(:,:)
      integer start(3),counts(3)

      ierr = nf_open(trim(filename), NF_WRITE, ncid0)
c      ncid0 = nid

      ip=mh
      jp=nh
      ierr=nf_inq_varid(ncid0,'time',vid)
      ierr=nf_put_vara_int(ncid0,vid,nrec,1,idate)
      
      ierr=nf_inq_varid(ncid0,trim(varname),vid)
      if(ierr.ne.NF_NOERR)then
         Write(*,*)'Error: file '//trim(filename)//
     *      '. No variable '//trim(varname)
      endif
      kk=nrec

      allocate(r4vals(ip,jp))

      Do i=1,ip
         Do j=1,jp
            if(mask(i,j).ne.0)then
               r4vals(i,j)=SNGL(array(i,j))
            else
               r4vals(i,j)=NF_FILL_FLOAT
            endif
         enddo
      enddo

      start(1)=1
      start(2)=1
      start(3)=kk
      counts(1)=ip
      counts(2)=jp
      counts(3)=1
      ierr=nf_put_vara_real(ncid0,vid,start,counts,r4vals)
      if(ierr.ne.NF_NOERR)then
         Write(*,*)'Error: file '//trim(filename)//
     *      '.  Variable '//trim(varname)
      endif

      deallocate(r4vals)
      
      ierr = nf_close(ncid0)

      Return
      End

      Subroutine close_nc(nid)

      integer, intent(inout) :: nid
c      character*(*), intent(in) :: filename
      
      include 'netcdf.inc'

      integer ierr

c         NF_CLOSE ! close netCDF dataset

      ierr = nf_close(nid)

      if(ierr.ne.NF_NOERR)then
         Write(*,*)'Error: file ',nid,' cannot be closed'
         Stop
      else
         Write(1001,'(a5,i8,a10)')'File ',nid,' is closed'
      endif

      Return
      End

      subroutine prepare_coeffs(
     *   rbox00, rbox10, rbox01, rbox11, i0box, j0box,
     *   i1box, j1box, nodetype,
     *   idat, jdat, Xdat, Ydat, Rmask,
     *   imt,  jmt,  X,    Y,    mask)
      
      include 'netcdf.inc'
      
      integer*4, intent(in) :: idat, jdat
      real*8, intent(in) :: Xdat(idat), Ydat(jdat), Rmask(idat,jdat)
      integer*4, intent(in) :: imt, jmt, mask(imt,jmt)
      real*8, intent(in) :: X(imt,jmt), Y(imt,jmt)

      real*8, dimension(imt,jmt), intent(out) :: 
     *   rbox00,rbox10,rbox01,rbox11
      integer, dimension(imt,jmt), intent(out) :: i0box,j0box, 
     *   i1box,j1box, nodetype

      real*8,parameter :: pi=3.14159265/180.0
      character*7 mode

      common/iodiag/nu_diag
      
      mode='levitus'
      eps=0.001

      do i=1,imt
         do j=1,jmt
            xx=X(i,j)
            yy=Y(i,j)

            if(xx.lt.0)xx=xx+360.0

            i0=0
		  do ii=1,idat-1
               if(Xdat(ii).le.xx.and.Xdat(ii+1).ge.xx)then
                  i0=ii
                  i1=ii+1
                  rx0=(xx-Xdat(i0))+eps
                  rx1=(Xdat(i1)-xx)+eps
                  exit
               endif
            enddo

            if(i0.eq.0)then
               if(Xdat(idat).le.xx.and.(Xdat(1)+360).ge.xx)then
                  i0=idat
                  i1=1
                  rx0=(xx-Xdat(i0))+eps
                  rx1=(Xdat(i1)+360-xx)+eps
               elseif(Xdat(1).ge.xx.and.(Xdat(idat)-360).le.xx)then
                  i0=idat
                  i1=1
                  rx0=(xx-Xdat(i0)+360)+eps
                  rx1=(Xdat(i1)-xx)+eps
               else
                  Write(nu_diag,*) 'Grid from NetCDF is not covering'
                  Write(nu_diag,*) xx,Xdat(1),Xdat(idat)
                  call exit_coupler
               endif
            endif

            j0=0
            if(yy.lt.Ydat(1))then
               j0=1
               j1=1
               nodetype(i,j)=1
               ry0=(yy+90.0)+eps
               ry1=(Ydat(1)-yy)+eps
            else if(yy.gt.Ydat(jdat))then
               j0=jdat
               j1=jdat
               nodetype(i,j)=1
               ry0=(yy-Ydat(jdat))+eps
               ry1=(90.0-yy)+eps
            else
               do jj=1,jdat-1
                  if(Ydat(jj).le.yy.and.Ydat(jj+1).ge.yy)then
                     j0=jj
                     j1=jj+1
                     ry0=(yy-Ydat(j0))+eps
                     ry1=(Ydat(j1)-yy)+eps
                     exit
                  endif
               enddo
               nodetype(i,j)=0
            endif
            
            if(mode.eq.'kriging')then
               w00=1.0/(rx0**2+ry0**2)
               w10=1.0/(rx1**2+ry0**2)
               w01=1.0/(rx0**2+ry1**2)
               w11=1.0/(rx1**2+ry1**2)
            elseif(mode.eq.'levitus')then
               rad2=2.25*((rx1+rx0)**2+(ry1+ry0)**2)/4.0
               w00=exp(-(rx0**2+ry0**2)/rad2)
               w10=exp(-(rx1**2+ry0**2)/rad2)
               w01=exp(-(rx0**2+ry1**2)/rad2)
               w11=exp(-(rx1**2+ry1**2)/rad2)
            else
               write(nu_diag,*)
     * 'Interpolation mode '//mode//' is not defined'
               Stop
            endif
            
            wt=
     *         (Rmask(i0,j0)+eps)*w00+
     *         (Rmask(i1,j0)+eps)*w10+
     *         (Rmask(i0,j1)+eps)*w01+
     *         (Rmask(i1,j1)+eps)*w11
     
            rbox00(i,j)=(Rmask(i0,j0)+eps)*w00/wt
            rbox10(i,j)=(Rmask(i1,j0)+eps)*w10/wt
            rbox01(i,j)=(Rmask(i0,j1)+eps)*w01/wt
            rbox11(i,j)=(Rmask(i1,j1)+eps)*w11/wt

            i0box(i,j)=i0
            i1box(i,j)=i1
            j0box(i,j)=j0
            j1box(i,j)=j1
         enddo
      enddo
	   
      return
      end


      subroutine projection(idat, jdat, PB, imt,jmt, mask, PP, 
     *   rbox00, rbox10, rbox01, rbox11, 
     *   i0box, j0box, i1box, j1box, nodetype)
      
      integer*4, intent(in) :: idat, jdat
      real*8, dimension(idat,jdat), intent(in) :: PB
      integer*4, intent(in) :: imt, jmt
      integer*4, dimension(imt,jmt), intent(in) :: mask
      real*8, dimension(imt,jmt), intent(out) :: PP
      real*8, dimension(imt,jmt), intent(in) :: 
     *   rbox00,rbox10,rbox01,rbox11
      integer, dimension(imt,jmt), intent(in) :: i0box,j0box, 
     *   i1box,j1box, nodetype

      logical flat
      integer nu_diag
      common/iodiag/nu_diag

      flat=.false.

      PMIN=1e+10
      PMAX=-1e+10

      do j=1,jmt
         do i=1,imt
            if(mask(i,j).ne.0)then
               i0=i0box(i,j)
               i1=i1box(i,j)
               j0=j0box(i,j)
               j1=j1box(i,j)

               if(nodetype(i,j).eq.0)then
                  PP(i,j)=
     *               rbox00(i,j)*PB(i0,j0)+rbox10(i,j)*PB(i1,j0)+
     *               rbox01(i,j)*PB(i0,j1)+rbox11(i,j)*PB(i1,j1)
c                  Print *,PP(i,j),rbox00(i,j),rbox10(i,j),rbox01(i,j),
c     *            rbox11(i,j),PB(i0,j0),PB(i1,j0),PB(i0,j1),PB(i1,j1)
               else
                  PP(i,j)=0.0

                  rr=1.0/idat
                  PPOLE=0.0
                  do ii=1,idat
                     PPOLE=PPOLE+rr*PB(ii,j0)
                  enddo
                  if(flat)then
                     PP(i,j)=PPOLE
                  else
                     PP(i,j) = 
     *                  rbox00(i,j)*PB(i0,j0)+rbox10(i,j)*PB(i1,j0)+
     *                  (rbox01(i,j)+rbox11(i,j))*PPOLE
                  endif
c                  Print *,PP(i,j),rbox00(i,j),rbox10(i,j),rbox01(i,j),
c     *            rbox11(i,j),PB(i0,j0),PB(i1,j0),PPOLE
               endif
	         if(PMAX.lt.PP(i,j))PMAX=PP(i,j)
	         if(PMIN.gt.PP(i,j))PMIN=PP(i,j)
            endif
         enddo
      enddo
c      Stop

c      Write(nu_diag,*)'Max: ', PMAX
c      Write(nu_diag,*)'Min: ', PMIN
      
      return
      end

      Subroutine rotate_uv(imt,jmt,X,Y,U,V,mask)
      real*8, dimension(imt,jmt) :: X,Y,U,V
      integer, dimension(imt,jmt) :: mask

      real*8, allocatable :: csn(:,:),snn(:,:)
      logical firsttime
      data firsttime/.true./

      save csn, snn

      pp=3.14159265/180.0

      if (firsttime) then

         firsttime=.false.
         allocate(csn(imt,jmt),snn(imt,jmt))

         do i=1,imt
            do j=1,jmt
               if(mask(i,j).ne.0)then

                  if(i.ne.1.and.i.ne.imt)then
                     xx1=(90.0-Y(i+1,j))*sin(X(i+1,j)*pp)
                     yy1=-(90-Y(i+1,j))*cos(X(i+1,j)*pp)
                     xx2=(90-Y(i-1,j))*sin(X(i-1,j)*pp)
                     yy2=-(90-Y(i-1,j))*cos(X(i-1,j)*pp)
                  else
                     if(i.eq.1)then
                        xx1=(90-Y(i+1,j))*sin(X(i+1,j)*pp)
                        yy1=-(90-Y(i+1,j))*cos(X(i+1,j)*pp)
                        xx2=(90-Y(i,j))*sn
                        yy2=-(90-Y(i,j))*cs
                     else
                        xx1=(90-Y(i,j))*sn
                        yy1=-(90-Y(i,j))*cs
                        xx2=(90-Y(i-1,j))*sin(X(i-1,j)*pp)
                        yy2=-(90-Y(i-1,j))*cos(X(i-1,j)*pp)
                     endif
                  endif
                  dx=xx1-xx2
                  dy=yy1-yy2

                  alph=atan2(dy,dx)
                  csn(i,j)=cos(alph)
                  snn(i,j)=sin(alph)

               endif

            enddo
         enddo
      endif

      do i=1,imt
         do j=1,jmt
            if(mask(i,j).ne.0)then

               sn=sin(X(i,j)*pp)
               cs=cos(X(i,j)*pp)

               uu0=U(i,j)*cs-V(i,j)*sn
               vv0=U(i,j)*sn+V(i,j)*cs
               U(i,j)=uu0*csn(i,j)+vv0*snn(i,j)
               V(i,j)=-uu0*snn(i,j)+vv0*csn(i,j)

            endif
         enddo
      enddo

      return
      end

