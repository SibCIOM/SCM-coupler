!     $Id: rivers.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $

      subroutine get_rivers(rivcoords,rtrans,nrivers,idate,isec)

      include 'netcdf.inc'

      real*8 rivcoords(*),rtrans(*)
      integer ierr,start(3),counts(3)

      logical use_mask
      character*80 datapath,respath,fname,filename,file_small_riv,i_name,j_name
      character*7 frmt
      character*1 writemode,freq
      common/progn_run/datapath,respath,use_mask,writemode
      logical use_ciaf,tidal,use_frac,use_frcp,use_xie,use_xie4ncep,
     &   use_datahum,old_tides,aomip_rules,prognostic_run,use_tref,
     &   use_small_rivers
      integer ciaf_vers,cld_lw_vers,cld_sw_vers,first_riv_year,
     &   last_riv_year,nyear0
      real*8 solar_avg,cld_power
      common/pars/cld_power,use_ciaf,ciaf_vers,tidal,
     &  use_frac,use_frcp,use_tref,
     &  cld_lw_vers,cld_sw_vers,use_xie,use_xie4ncep,
     &  use_datahum,old_tides,aomip_rules,prognostic_run,solar_avg,
     &  use_small_rivers,file_small_riv,first_riv_year,last_riv_year

      integer, save :: nfi0,nre0,id_nc,vid,mon(12),imt,jmt,nyr0,nyr1
      integer, save, allocatable :: river_mask(:,:)
      logical, save :: firsttime, mask_is_empty
      real*8, allocatable, save :: dat(:,:),area(:,:)
      real*8, allocatable :: lon(:,:),lat(:,:),lonv(:,:,:),latv(:,:,:)
      real*8 vec1(3),vec2(3),vec0(3)
      real*4 missing,fillval
      real*4, allocatable :: r4vals(:,:)
      data nfi0/-1/,nre0/-1/,firsttime/.true./,nyear0/0/,
     &  mon/0,31,59,90,120,151,181,212,243,273,304,334/

      if(prognostic_run)then
         nyear=floor(dble(idate)/10000.0)
         nmon=idate-10000*nyear
         nmon=floor(dble(nmon)/100.0)
         nday=idate-10000*nyear-100*nmon
         if(firsttime)then
            Open(1002,file=trim(datapath)//'file_content',form='formatted')
            nfir=0
            do while(.true.)
               Read(1002,*,err=2,end=2)fname
               if(trim(fname).eq.'friver')nfir=1
            enddo
    2       Close(1002)
            if(nfir.eq.0)Stop 'Cannot use friver, it doesnot exist'
         endif
         if(nyear.ne.nyear0)then
         nfi=0
         Open(1002,file=trim(datapath)//'friv_content',form='formatted')
         Read(1002,*)freq
         Read(1002,*)i_name
         Read(1002,*)j_name
         do while(.true.)
            Read(1002,*,err=1,end=1)nfir,nyr0,nyr1
            if(nyear.ge.nyr0.and.nyear.le.nyr1)then
               nfi=nfir
               goto 1
            endif
         enddo
    1    Close(1002)
         if(nfi.eq.0)Stop 'Cannot read friver files for a specified year'
         nyear0=nyear
         else
         nfi=nfi0
         endif

c         nfi=floor(0.1*dble(nyear-2006))+1
         if(nfi.ne.nfi0)then

            if(nfi0.ge.0)then

               call close_nc(id_nc)

            endif

            nfi0=nfi

            frmt='(a6,i1)'
            if(nfi.gt.9)frmt='(a6,i2)'

            Write(fname,frmt)'friver',nfi
            filename=trim(datapath)//trim(fname)//'.nc'
               

c         NF_OPEN ! open existing netCDF dataset

            ierr = nf_open(trim(filename), NF_NOWRITE, id_nc)

            if(ierr.ne.NF_NOERR)then
               Write(*,*)nf_strerror(ierr)
               Write(*,*)'Error: file "'//trim(filename)//
     *            '" is not present or damaged'
               Stop
            else
               Write(*,'(a)')'File '//trim(filename)//' is open'
               Write(*,'(10x,a7,i8)')'ncid = ',id_nc
            endif

            nre0=-1

            if(firsttime)then

               firsttime=.false.

               ierr=nf_inq_dimid(id_nc,trim(i_name),vid)
               if(ierr.ne.NF_NOERR)then
                  Write(*,*)'Error: dimid i does not exist'
                  Stop
               endif
               ierr = nf_inq_dimlen(id_nc, vid, imt)

               ierr=nf_inq_dimid(id_nc,trim(j_name),vid)
               if(ierr.ne.NF_NOERR)then
                  Write(*,*)'Error: dimid j does not exist'
                  Stop
               endif
               ierr = nf_inq_dimlen(id_nc, vid, jmt)

               Write(*,*)'Size: ',imt,jmt
c               Write(*,*)'Allocating arrays'
               allocate(lon(imt,jmt),lat(imt,jmt),lonv(4,imt,jmt),latv(4,imt,jmt),dat(imt,jmt),
     *            river_mask(imt,jmt),area(imt,jmt))

               start=1
               counts(1)=imt
               counts(2)=jmt

c               Write(*,*)'Fetching longitude'
               ierr=nf_inq_varid(id_nc,'lon',vid)
               ierr = nf_get_vara_double(id_nc,vid,start,counts,lon)

c               Write(*,*)'Fetching latitude'
               ierr=nf_inq_varid(id_nc,'lat',vid)
               ierr = nf_get_vara_double(id_nc,vid,start,counts,lat)

               counts(1)=4
               counts(2)=imt
               counts(3)=jmt

c               Write(*,*)'Fetching longitude vertices'
               ierr=nf_inq_varid(id_nc,'lon_vertices',vid)
               ierr = nf_get_vara_double(id_nc,vid,start,counts,lonv)

c               Write(*,*)'Fetching latitude vertices'
               ierr=nf_inq_varid(id_nc,'lat_vertices',vid)
               ierr = nf_get_vara_double(id_nc,vid,start,counts,latv)
               
               Rearth=6371000.0   ! Earth radius in m
c               Write(*,*)'Calculating areas'
               area_mean=0.0
               pi=4.0*datan2(1.0,1.0)/180.0
               do i=1,imt
                  do j=1,jmt
                  ! Calculate 3D coordinates of vertices as if the Earth were of unit radius
c                     Write(*,*)'Lonv: ', lonv(:,i,j)
c                     Write(*,*)'Latv: ', latv(:,i,j)
                     xglo=cos(latv(2,i,j)*pi)*cos(lonv(2,i,j)*pi)
                     yglo=cos(latv(2,i,j)*pi)*sin(lonv(2,i,j)*pi)
                     zglo=sin(latv(2,i,j)*pi)
                     xgl1=cos(latv(1,i,j)*pi)*cos(lonv(1,i,j)*pi)
                     ygl1=cos(latv(1,i,j)*pi)*sin(lonv(1,i,j)*pi)
                     zgl1=sin(latv(1,i,j)*pi)
                     xgl2=cos(latv(3,i,j)*pi)*cos(lonv(3,i,j)*pi)
                     ygl2=cos(latv(3,i,j)*pi)*sin(lonv(3,i,j)*pi)
                     zgl2=sin(latv(3,i,j)*pi)
                  ! Calculate upper-left triangle area
                     vec1(1)=xgl1-xglo
                     vec1(2)=ygl1-yglo
                     vec1(3)=zgl1-zglo
c                     Write(*,*)'Vec1: ',vec1
                     vec2(1)=xgl2-xglo
                     vec2(2)=ygl2-yglo
                     vec2(3)=zgl2-zglo
c                     Write(*,*)'Vec2: ',vec2
                     vec0(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
                     vec0(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
                     vec0(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)
c                     Write(*,*)'Vec0: ',vec0
                     area(i,j)=0.0
                     Sarea=0.0
                     do k=1,3
                        Sarea=Sarea+vec0(k)**2
                     enddo
                     area(i,j)=0.5*sqrt(Sarea)
                     xglo=cos(latv(4,i,j)*pi)*cos(lonv(4,i,j)*pi)
                     yglo=cos(latv(4,i,j)*pi)*sin(lonv(4,i,j)*pi)
                     zglo=sin(latv(4,i,j)*pi)
                  ! Calculate lover-right triangle area
                     vec1(1)=xgl1-xglo
                     vec1(2)=ygl1-yglo
                     vec1(3)=zgl1-zglo
c                     Write(*,*)'Vec1: ',vec1
                     vec2(1)=xgl2-xglo
                     vec2(2)=ygl2-yglo
                     vec2(3)=zgl2-zglo
c                     Write(*,*)'Vec2: ',vec2
                     vec0(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
                     vec0(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
                     vec0(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)
c                     Write(*,*)'Vec0: ',vec0
                     Sarea=0.0
                     do k=1,3
                        Sarea=Sarea+vec0(k)**2
                     enddo
                     area(i,j)=area(i,j)+0.5*sqrt(Sarea)
                  ! Substitute the real Earth's radius
                     area(i,j)=area(i,j)*Rearth**2
                     area_mean=area_mean+area(i,j)
                  enddo
               enddo
               Write(*,*)'Mean area: ',area_mean/dble(imt*jmt)
c              We don't need vertices coordinates any more
c               Write(*,*)'deallocating'
               deallocate(lonv,latv)
               
               river_mask=0
               mask_is_empty=.true.

            endif

         endif

         if(freq.eq.'d')then           ! daily records without leap years
            nre=(nyear-nyr0)*365+mon(nmon)+nday
         elseif(freq.eq.'v')then       ! daily records with leap years
            nre=(nyear-nyr0)*365+mon(nmon)+nday
            do k=nyr0,nyear
               if(4*int(0.25*k).eq.k.and.nmon.gt.2)nre=nre+1
            enddo
         elseif(freq.eq.'m')then       ! monthly records
            nre=(nyear-nyr0)*12+nmon
         endif

         if(nre.ne.nre0)then

c               Write(*,*)'allocating r4'
               allocate(r4vals(imt,jmt))

               start=1
               start(3)=nre
               counts(1)=imt
               counts(2)=jmt
               counts(3)=1
               nre0=nre

               Write(*,*)'Fetching friver: nre= ',nre
               ierr=nf_inq_varid(id_nc,'friver',vid)
               ierr=nf_get_vara_real(id_nc,vid,start,counts,r4vals)
               Write(*,*)'Min= ',minval(r4vals),'; Max= ',maxval(r4vals)

c               Write(*,*)'Fetching missing_value'
               ierr=nf_get_att_real(id_nc,vid,'missing_value',missing)
               if(ierr.ne.NF_NOERR)missing=1.e+20
c               Write(*,*)'Fetching _FillValue'
               ierr=nf_get_att_real(id_nc,vid,'_FillValue',fillval)
               if(ierr.ne.NF_NOERR)fillval=1.e+20
c               Write(*,*)'Fetching scale_factor'
c               ierr=nf_get_att_real(id_nc,vid,'scale_factor',sc)
c               if(ierr.ne.NF_NOERR)sc=1.0
c               Write(*,*)'Fetching add_offset'
c               ierr=nf_get_att_real(id_nc,vid,'add_offset',ad)
c               if(ierr.ne.NF_NOERR)ad=0.0
c               Write(*,*)'Sc= ',sc,'; Ad= ',ad

               if(mask_is_empty)then
c               Write(*,*)'Filling area mask: 1-st step'
c         Mask out missing values and fillvals
                  do i=1,imt
                     do j=1,jmt
                        if(lat(i,j).lt.-20)then
                           river_mask(i,j)=-1
                        endif
                     enddo
                  enddo
c               Write(*,*)'Filling area mask: 2-d step'
c         Mask out all nodes with latitudes < -20S
                  do i=1,imt
                     do j=1,jmt
                        if(r4vals(i,j).eq.missing.or.r4vals(i,j).eq.fillval)then
                           river_mask(i,j)=-1
                        endif
                     enddo
                  enddo
c               Write(*,*)'Filling area mask: 3-d step'
c         Mask out all nodes beyond Bering Strait
                  do i=1,imt
                     do j=1,jmt
                        if(lat(i,j).lt.65.776.and.
     *                     ((lon(i,j).gt.-174.7.and.lon(i,j).lt.-160.3).or.
     *                     (lon(i,j).gt.185.3.and.lon(i,j).lt.199.7)))then
                           river_mask(i,j)=-1
                        endif
                     enddo
                  enddo
c               Write(*,*)'Filling area mask: 4-th step'
c         Mask out all nodes beyond Gibraltar
                  do i=1,imt
                     do j=1,jmt
                        if(lat(i,j).lt.40.0.and.lat(i,j).gt.31.and.
     *                     ((lon(i,j).gt.-5.6.and.lon(i,j).lt.0.0).or.
     *                     (lon(i,j).gt.354.4)))then
                           river_mask(i,j)=-1
                        endif
                     enddo
                  enddo
c               Write(*,*)'Filling area mask: 5.1-th step'
c         Mark square in the middle of North Atlantic
                  do i=1,imt
                     do j=1,jmt
                        if(lat(i,j).lt.41.0.and.lat(i,j).gt.19.0.and.
     *                     ((lon(i,j).gt.-58.0.and.lon(i,j).lt.-25.0).or.
     *                     (lon(i,j).gt.302.0.and.lon(i,j).lt.335.0)))then
                           if(river_mask(i,j).ne.-1)river_mask(i,j)=999
                        endif
                     enddo
                  enddo
c               Write(*,*)'Filling area mask: 5.2-th step'
c         Mark square in the middle of Arctic
                  do i=1,imt
                     do j=1,jmt
                        if(lat(i,j).lt.83.3.and.lat(i,j).gt.72.6.and.
     *                     lon(i,j).gt.154.0.and.lon(i,j).lt.180.0)then
                           if(river_mask(i,j).ne.-1)river_mask(i,j)=999
                        endif
                     enddo
                  enddo
c               Write(*,*)'Filling area mask: 6-th step'
c         Mark all nodes which are not masked out and which are close to marked nodes
                  k=999
                  do while(k.gt.0)
                     k=0
                     do i=1,imt
                        im=max(1,i-1)
                        ip=min(imt,i+1)
                        do j=1,jmt
                           jm=max(1,j-1)
                           jp=min(jmt,j+1)
                           if(river_mask(i,j).eq.0.and.(
     *                        river_mask(im,j).eq.999.or.river_mask(i,jm).eq.999.or.
     *                        river_mask(ip,j).eq.999.or.river_mask(i,jp).eq.999))then
                              river_mask(i,j)=999
                              k=k+1
                           endif
                        enddo
                     enddo
                  enddo
c               Write(*,*)'Filling area mask: 7-th step'
c          Mask out all rested nodes                  
                  do i=1,imt
                     do j=1,jmt
                        if(river_mask(i,j).ne.999)river_mask(i,j)=0
                     enddo
                  enddo
c          Search nearest river mouth to each masked nodes
                  pi=4.0*datan2(1.0,1.0)/180.0
c                  Write(*,*)'Lats: ',lat(imt/2,:)
c                  Write(*,*)'Rivs: ',(rivcoords(2*k),k=1,nrivers)
c                  Write(*,*)'Lons: ',lon(:,jmt/2)
c                  Write(*,*)'Rivs: ',(rivcoords(2*k-1),k=1,nrivers)
c                  Write(*,*)'Searching nearest mouth'
                  do i=1,imt
                     do j=1,jmt
                        if(river_mask(i,j).ne.0)then
                           xglo=cos(lat(i,j)*pi)*cos(lon(i,j)*pi)
                           yglo=cos(lat(i,j)*pi)*sin(lon(i,j)*pi)
                           zglo=sin(lat(i,j)*pi)
                           rmin=1e+20
                           kk=0
                           do k=1,nrivers
                              ulon=rivcoords(2*k-1)
                              ulat=rivcoords(2*k)
                              xglu=cos(ulat*pi)*cos(ulon*pi)
                              yglu=cos(ulat*pi)*sin(ulon*pi)
                              zglu=sin(ulat*pi)
                              rcur=(xglu-xglo)**2+(yglu-yglo)**2+(zglu-zglo)**2
                              if(rcur.lt.rmin)then
                                 kk=k
                                 rmin=rcur
                              endif
                           enddo
                           if(kk.eq.0)then
                              Print *,'get_river: something wrong with node ',i,j
                              Print *,'           cannot find nearest river mouth'
                              Stop
                           endif
                           river_mask(i,j)=kk
                        endif
                     enddo
                  enddo
                  mask_is_empty=.false.
c              We don't need coordinates any more
c                  Write(*,*)'deallocating lat lon'
                  open(1002,file='mask',form='formatted')
                  write(1002,'(362(292i2,/))')((river_mask(i,j),j=1,jmt),i=1,imt)
                  close(1002)
                  deallocate(lat,lon)
               endif

               do i=1,imt
                  do j=1,jmt
                     if(r4vals(i,j).eq.missing.or.r4vals(i,j).eq.fillval)then
                        dat(i,j)=0.0
                     else
                        dat(i,j)=dble(r4vals(i,j))
                     endif
                  enddo
               enddo
               
c         Calculate river transport at by searching nearest nodes
               Write(*,*)'Calculating river transports'
               rho_fresh=1000.0
               rtrans(1:nrivers)=0.0
               do i=1,imt
                  do j=1,jmt
                     if(river_mask(i,j).ne.0)then
                        k=river_mask(i,j)
                        rtrans(k)=rtrans(k)+dat(i,j)*area(i,j)/rho_fresh
c                        Write(*,*)k,dat(i,j),area(i,j),rtrans(k)
                     endif
                  enddo
               enddo

c               Write(*,*)'deallocating r4'
               deallocate(r4vals)

         endif

      else

         rtrans(1:nrivers)=0.0

      endif
c      Write(*,*)'Done get_river'
c      Write(*,*)'Rtrans: ',rtrans(16),rtrans(10)
c      stop
      return
      end
