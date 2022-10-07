!     $Id: ncepdata.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $

      subroutine pressure(idate, imt, jmt, X, Y, mask, PBOT)
      
      include 'netcdf.inc'
      
      integer*4, intent(in) :: idate, imt, jmt, mask(imt,jmt)
      real*8, intent(in) :: X(imt,jmt), Y(imt,jmt)
      real*8, intent(out) :: PBOT(imt,jmt)
      real*8, allocatable, dimension(:,:) ::
     &  rxbox,rybox,PP
      integer, allocatable, dimension(:,:) ::
     &  i0box,j0box, i1box,j1box, nodetype
      integer dindex(12)
      integer i,j,n,ncid,vid,did,beg(3),cnt(3),ierr,nold
      integer dday, dmonth, dyear, dyearold
      character*80 attstr
      character*80 fname
      character, allocatable :: str(:)
      integer,parameter :: IDAT=144
      integer,parameter :: JDAT=73
      real*8,parameter :: XSTEPDAT=2.5
      real*8,parameter :: YSTEPDAT=2.5
      real*8,parameter :: pi=3.14159265/180.0
      real*8 ::  Xdat(IDAT), Ydat(JDAT)
      real*4 ::  Xdat0(IDAT), Ydat0(JDAT)
      Integer*2,dimension(IDAT,JDAT) ::  PB0
      Real*4,dimension(IDAT,JDAT) ::  PB4
      Real*4 attrib
      Real*8 offset,scale
      Real*8,dimension(IDAT,JDAT) ::  PB
      character*80 input_path
      common/input_data/input_path
      character*4 cyear

      save

      logical firstcall, newyear, flat
      integer nu_diag
      common/iodiag/nu_diag
      save firstcall
      data firstcall/.true./,nold/-1/,dyearold/-1/,dindex
     */0,31,59,90,120,151,181,212,243,273,304,334/

      flat=.false.
      dyear=int(idate/10000)
      newyear=dyearold.ne.dyear
      
      dmonth=int((idate-10000*dyear)/100)
      dday=int(idate-10000*dyear-100*dmonth)
      
      n=dindex(dmonth)+dday
c      if(4*int(dyear/4).eq.4.and.dmonth.gt.2)n=n+1
      
c      Write(nu_diag,*)'In Pressure'
      
	if(firstcall)then
	
      allocate(rxbox(imt,jmt),rybox(imt,jmt),i0box(imt,jmt),
     *   j0box(imt,jmt),i1box(imt,jmt),j1box(imt,jmt),nodetype(imt,jmt),
     *   PP(imt,jmt))

         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'slp/slp.'//cyear//'.nc'
         ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'Pressure Error: reading NetCDF file '''//fname//''''
            call exit_coupler
         endif

      beg(1)=1
      beg(2)=1
      beg(3)=1
      cnt(1)=IDAT
      cnt(2)=1
      cnt(3)=1
      
      ierr=nf_inq_varid(ncid,'lon',vid)
      ierr=nf_get_vara_real(ncid,vid,beg,cnt,Xdat0(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'Pressure Error: Error reading lon'
         call exit_coupler
      endif

      cnt(1)=JDAT

      ierr=nf_inq_varid(ncid,'lat',vid)
      ierr=nf_get_vara_real(ncid,vid,beg,cnt,Ydat0(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'Pressure Error: Error reading lat'
         call exit_coupler
      endif
      ierr=nf_close(ncid)

	   do i=1,IDAT
	      Xdat(i)=dble(Xdat0(i))
	   enddo
	   do j=1,JDAT
	      Ydat(j)=dble(Ydat0(JDAT-j+1))
           enddo

	   do i=1,imt
	      do j=1,jmt
	         xx=X(i,j)
		 if(xx.lt.0)xx=xx+360.0
		 yy=Y(i,j)
		 mm=mask(i,j)
		    i0=0
		    do ii=1,IDAT-1
		       if(Xdat(ii).le.xx.and.Xdat(ii+1).ge.xx)then
		          i0=ii
			  i1=ii+1
			  rxbox(i,j)=(xx-Xdat(i0))/(Xdat(i1)-Xdat(i0))
			  exit
		       endif
		    enddo
		    if(i0.eq.0)then
		       if(Xdat(IDAT).le.xx.and.(Xdat(1)+360).ge.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0))/
     *  		    (Xdat(i1)+360-Xdat(i0))
                   elseif(Xdat(1).ge.xx.and.(Xdat(IDAT)-360).le.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0)+360)/
     *  		    (Xdat(i1)+360-Xdat(i0))
		       else
                          Write(nu_diag,*) 'Pressure: Grid from NetCDF is not covering'
			  Write(nu_diag,*) xx,Xdat(1),Xdat(IDAT)
                          call exit_coupler
		       endif
		    endif
		    j0=0
		    if(yy.le.Ydat(2))then
		       j0=2
		       j1=2
		       nodetype(i,j)=1
                       rybox(i,j)=(yy+90.0)/(Ydat(j1)+90.0)
		    else if(yy.ge.Ydat(JDAT-1))then
		       j0=JDAT-1
		       j1=JDAT-1
		       nodetype(i,j)=1
                       rybox(i,j)=(yy-Ydat(j0))/(90.0-Ydat(j0))
		    else
		    do jj=1,JDAT-1
		       if(Ydat(jj).le.yy.and.Ydat(jj+1).ge.yy)then
		          j0=jj
			  j1=jj+1
			  rybox(i,j)=(yy-Ydat(j0))/(Ydat(j1)-Ydat(j0))
			  exit
		       endif
		    enddo
		    nodetype(i,j)=0
		    endif
		    i0box(i,j)=i0
		    i1box(i,j)=i1
		    j0box(i,j)=j0
		    j1box(i,j)=j1
	      enddo
	   enddo
	   
        endif
	
      if(n.ne.nold)then
      
      if(newyear)then

         if(dyearold.ge.0)then      
	   ierr=nf_close(ncid)
	 endif
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'slp/slp.'//cyear//'.nc'
c	 Write(nu_diag,*)'Reading file: ',fname
	   ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
	   if(ierr.NE.NF_NOERR)then
              Write(nu_diag,*) 'Pressure Error: reading NetCDF file'
              call exit_coupler
c           else
c              Write(nu_diag,'(a,i3)')'Open CDF ', ncid
           endif
	   dyearold=dyear
      endif

c      Write(nu_diag,'(a,i3)')'Read record no.',n,' from CDF ', ncid

      beg(1)=1
      beg(2)=1
      beg(3)=n
      cnt(1)=IDAT
      cnt(2)=JDAT
      cnt(3)=1
      
      ierr=nf_inq_varid(ncid,'slp',vid)
      if(dyear.lt.2014)then
         ierr=nf_get_att_real(ncid,vid,'add_offset',attrib)
         offset=dble(attrib)
         ierr=nf_get_att_real(ncid,vid,'scale_factor',attrib)
         scale=dble(attrib)
c         Write(nu_diag,*)'Offset=',offset
c         Write(nu_diag,*)'Factor=',scale
         ierr=nf_get_vara_int2(ncid,vid,beg,cnt,PB0(1,1))

         Do i=1,IDAT
            Do j=1,JDAT
               PB(i,j)=dble(PB0(i,JDAT-j+1))
            Enddo
         Enddo
      else
         offset=0.0
         scale=1.0
         ierr=nf_get_vara_real(ncid,vid,beg,cnt,PB4(1,1))

         Do i=1,IDAT
            Do j=1,JDAT
               PB(i,j)=dble(PB4(i,JDAT-j+1))
            Enddo
         Enddo
      endif

c      Write(nu_diag,'(a,f7.2,a)')'Read Control: P(115,31)= ',
c     *   (PB(115,31)*scale+offset)/100.0,' hPa'

      PMIN=1e+10
      PMAX=-1e+10
      do j=1,jmt
         do i=1,imt
               i0=i0box(i,j)
	       i1=i1box(i,j)
	       j0=j0box(i,j)
	       j1=j1box(i,j)
	       rx0=rxbox(i,j)
	       rx1=1.0-rx0
	       ry0=rybox(i,j)
	       ry1=1.0-ry0
              if(nodetype(i,j).eq.0)then
               PP(i,j)= scale*(
     *            ry1*(rx1*PB(i0,j0)+rx0*PB(i1,j0))+
     *            ry0*(rx1*PB(i0,j1)+rx0*PB(i1,j1))) +
     *            offset
              else
               PP(i,j)=0.0
	       rr=1.0/IDAT
               if(flat)then
	       do ii=1,IDAT
                  PP(i,j)=PP(i,j)+rr*PB(ii,j0)
	       enddo
	       PP(i,j)=scale*PP(i,j)+offset
               else
c        Components at the pole as mean values
                  PPOLE=0.0
                  do ii=1,IDAT
                     uur=PB(ii,j0)*scale+offset
                     PPOLE=PPOLE+rr*uur
                  enddo

                  uu0=PB(i0,j0)*scale+offset
                  uu1=PB(i1,j0)*scale+offset

                  uu=rx1*uu0+rx0*uu1

                  PP(i,j)=ry1*uu+ry0*PPOLE
               endif
            endif
	    if(PMAX.lt.PP(i,j))PMAX=PP(i,j)
	    if(PMIN.gt.PP(i,j))PMIN=PP(i,j)
         enddo
      enddo

c      Write(nu_diag,*)'Max pressure: ', PMAX/100.0
c      Write(nu_diag,*)'Min pressure: ', PMIN/100.0
      
      endif
      
      nold=n
      
      do j=1,jmt
         do i=1,imt
	    PBOT(i,j)=PP(i,j)
         enddo
      enddo

      
      firstcall=.false.

      return
      end

      subroutine wind(idate, imt, jmt, X, Y, mask, UBOT, VBOT)
      
      include 'netcdf.inc'
      
      integer*4, intent(in) :: idate, imt, jmt, mask(imt,jmt)
      real*8, intent(in) :: X(imt,jmt), Y(imt,jmt)
      real*8, intent(out) :: UBOT(imt,jmt), VBOT(imt,jmt)
      real*8, allocatable, dimension(:,:) ::
     &  rxbox,rybox,UP,VP
      integer, allocatable, dimension(:,:) ::
     &  i0box,j0box, i1box,j1box, nodetype
      integer dindex(12)
      integer i,j,n,ncidu,ncidv,vid,did,beg(3),cnt(3),ierr,nold
      integer dday, dmonth, dyear, dyearold
      character*80 attstr
      character*80 fname
      character, allocatable :: str(:)
      integer,parameter :: IDAT=144
      integer,parameter :: JDAT=73
      real*8,parameter :: XSTEPDAT=2.5
      real*8,parameter :: YSTEPDAT=2.5
      real*8,parameter :: pi=3.14159265/180.0
      real*8 ::  Xdat(IDAT), Ydat(JDAT)
      real*4 ::  Xdat0(IDAT), Ydat0(JDAT)
      Integer*2,dimension(IDAT,JDAT) ::  WORK
      Real*4,dimension(IDAT,JDAT) ::  WORK4
      Real*4 attrib
      Real*8 uoffset,uscale,voffset,vscale
      Real*8,dimension(IDAT,JDAT) ::  UB, VB
      character*80 input_path
      common/input_data/input_path
      character*4 cyear

      save

      logical firstcall, newyear, flat, first
      integer nu_diag
      common/iodiag/nu_diag
      save firstcall
      data firstcall/.true./,nold/-1/,dyearold/-1/,dindex
     */0,31,59,90,120,151,181,212,243,273,304,334/

      flat=.false.
      dyear=int(idate/10000)
      newyear=dyearold.ne.dyear
      
      dmonth=int((idate-10000*dyear)/100)
      dday=int(idate-10000*dyear-100*dmonth)
      
      n=dindex(dmonth)+dday
c      if(4*int(dyear/4).eq.4.and.dmonth.gt.2)n=n+1
      
c      Write(nu_diag,*)'In Wind'
      
	if(firstcall)then
	
      allocate(rxbox(imt,jmt),rybox(imt,jmt),i0box(imt,jmt),
     *   j0box(imt,jmt),i1box(imt,jmt),j1box(imt,jmt),nodetype(imt,jmt),
     *   UP(imt,jmt),VP(imt,jmt))

         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'wnd/uwnd.sig995.'//cyear//'.nc'
         ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'U-Wind Error: reading NetCDF file'
            call exit_coupler
         endif

      beg(1)=1
      beg(2)=1
      beg(3)=1
      cnt(1)=IDAT
      cnt(2)=1
      cnt(3)=1
      
      ierr=nf_inq_varid(ncid,'lon',vid)
      ierr=nf_get_vara_real(ncid,vid,beg,cnt,Xdat0(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'U-Wind Error: Error reading lon'
         call exit_coupler
      endif

      cnt(1)=JDAT

      ierr=nf_inq_varid(ncid,'lat',vid)
      ierr=nf_get_vara_real(ncid,vid,beg,cnt,Ydat0(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'U-Wind Error: Error reading lat'
         call exit_coupler
      endif
      ierr=nf_close(ncid)

	   do i=1,IDAT
	      Xdat(i)=Xdat0(i)
	   enddo
	   do j=1,JDAT
	      Ydat(j)=Ydat0(JDAT-j+1)
           enddo
	   do i=1,imt
	      do j=1,jmt
	         xx=X(i,j)
		 if(xx.lt.0)xx=xx+360.0
		 yy=Y(i,j)
		 mm=mask(i,j)
		    i0=0
		    do ii=1,IDAT-1
		       if(Xdat(ii).le.xx.and.Xdat(ii+1).ge.xx)then
		          i0=ii
			  i1=ii+1
			  rxbox(i,j)=(xx-Xdat(i0))/(Xdat(i1)-Xdat(i0))
			  exit
		       endif
		    enddo
		    if(i0.eq.0)then
		       if(Xdat(IDAT).le.xx.and.(Xdat(1)+360).ge.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0))/
     *  		    (Xdat(i1)+360-Xdat(i0))
                   elseif(Xdat(1).ge.xx.and.(Xdat(IDAT)-360).le.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0)+360)/
     *  		    (Xdat(i1)+360-Xdat(i0))
		       else
                          Write(nu_diag,*) 'Wind: Grid from NetCDF is not covering'
			  Write(nu_diag,*) xx,Xdat(1),Xdat(IDAT)
                          call exit_coupler
		       endif
		    endif
		    j0=0
		    if(yy.le.Ydat(2))then
		       j0=2
		       j1=2
		       nodetype(i,j)=1
                       rybox(i,j)=(yy+90.0)/(Ydat(j1)+90.0)
		    else if(yy.ge.Ydat(JDAT-1))then
		       j0=JDAT-1
		       j1=JDAT-1
		       nodetype(i,j)=1
                       rybox(i,j)=(yy-Ydat(j0))/(90.0-Ydat(j0))
		    else
		    do jj=1,JDAT-1
		       if(Ydat(jj).le.yy.and.Ydat(jj+1).ge.yy)then
		          j0=jj
			  j1=jj+1
			  rybox(i,j)=(yy-Ydat(j0))/(Ydat(j1)-Ydat(j0))
			  exit
		       endif
		    enddo
		    nodetype(i,j)=0
		    endif
		    i0box(i,j)=i0
		    i1box(i,j)=i1
		    j0box(i,j)=j0
		    j1box(i,j)=j1
	      enddo
	   enddo
	   
        endif
	
      if(n.ne.nold)then
      
      if(newyear)then

         if(dyearold.ge.0)then      
	   ierr=nf_close(ncidu)
	 endif
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'wnd/uwnd.sig995.'//cyear//'.nc'
c	 Write(nu_diag,*)'Reading file: ',fname
	   ierr=nf_open(trim(fname),NF_NOWRITE,ncidu)
	   if(ierr.NE.NF_NOERR)then
              Write(nu_diag,*) 'U-Wind Error: reading NetCDF file'
              call exit_coupler
c           else
c              Write(nu_diag,'(a,i3)')'Open CDF ', ncidu
           endif
      endif

c      Write(nu_diag,'(a,i3)')'Read record no.',n,' from CDF ', ncidu

      beg(1)=1
      beg(2)=1
      beg(3)=n
      cnt(1)=IDAT
      cnt(2)=JDAT
      cnt(3)=1
      
      ierr=nf_inq_varid(ncidu,'uwnd',vid)
      if(dyear.lt.2014)then
         ierr=nf_get_att_real(ncidu,vid,'add_offset',attrib)
         uoffset=dble(attrib)
         ierr=nf_get_att_real(ncidu,vid,'scale_factor',attrib)
         uscale=dble(attrib)
c         Write(nu_diag,*)'Offset=',uoffset
c         Write(nu_diag,*)'Factor=',uscale
         ierr=nf_get_vara_int2(ncidu,vid,beg,cnt,WORK(1,1))
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'U-Wind Error: Error reading uwnd'
            call exit_coupler
         endif
         Do i=1,IDAT
            Do j=1,JDAT
               UB(i,j)=dble(WORK(i,JDAT-j+1))
            Enddo
         Enddo
c         UB=uscale*UB+uoffset
      else
         uoffset=0.0
         uscale=1.0
         ierr=nf_get_vara_real(ncidu,vid,beg,cnt,WORK4(1,1))
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'U-Wind Error: Error reading uwnd'
            call exit_coupler
         endif
         Do i=1,IDAT
            Do j=1,JDAT
               UB(i,j)=dble(WORK4(i,JDAT-j+1))
            Enddo
         Enddo
c         UB=uscale*UB+uoffset
      endif

      if(newyear)then

         if(dyearold.ge.0)then      
	   ierr=nf_close(ncidv)
	 endif
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'wnd/vwnd.sig995.'//cyear//'.nc'
c	 Write(nu_diag,*)'Reading file: ',fname
	   ierr=nf_open(trim(fname),NF_NOWRITE,ncidv)
	   if(ierr.NE.NF_NOERR)then
              Write(nu_diag,*) 'V-Wind Error: reading NetCDF file'
              call exit_coupler
c           else
c              Write(nu_diag,'(a,i3)')'Open CDF ', ncidv
           endif
	   dyearold=dyear
      endif

c      Write(nu_diag,'(a,i3)')'Read record no.',n,' from CDF ', ncidv

      beg(1)=1
      beg(2)=1
      beg(3)=n
      cnt(1)=IDAT
      cnt(2)=JDAT
      cnt(3)=1
      
      ierr=nf_inq_varid(ncidv,'vwnd',vid)
      if(dyear.lt.2014)then
         ierr=nf_get_att_real(ncidu,vid,'add_offset',attrib)
         voffset=dble(attrib)
         ierr=nf_get_att_real(ncidu,vid,'scale_factor',attrib)
         vscale=dble(attrib)
c         Write(nu_diag,*)'Offset=',voffset
c         Write(nu_diag,*)'Factor=',vscale
         ierr=nf_get_vara_int2(ncidv,vid,beg,cnt,WORK(1,1))
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'V-Wind Error: Error reading vwnd'
            call exit_coupler
         endif
         Do i=1,IDAT
            Do j=1,JDAT
               VB(i,j)=dble(WORK(i,JDAT-j+1))
            Enddo
         Enddo
c         VB=vscale*VB+voffset
      else
         voffset=0.0
         vscale=1.0
         ierr=nf_get_vara_real(ncidv,vid,beg,cnt,WORK4(1,1))
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'V-Wind Error: Error reading vwnd'
            call exit_coupler
         endif
         Do i=1,IDAT
            Do j=1,JDAT
               VB(i,j)=dble(WORK4(i,JDAT-j+1))
            Enddo
         Enddo
c         VB=vscale*VB+voffset
      endif

c      Write(nu_diag,'(a,f6.2,a,f6.2,a)')'Read Control: U(115,31)= ',
c     *   UB(115,31)*uscale+uoffset,' m/s,  V(115,31)= ',
c     *   VB(115,31)*vscale+voffset,' m/s'

      PMIN=1e+10
      PMAX=-1e+10
      pp=atan2(1.0,1.0)*4.0/180.0
c      if(i.eq.80.and.j.eq.148)Write(nu_diag,*)0,pp

      first=.true.
      do j=1,jmt
         do i=1,imt
            if(nodetype(i,j).eq.0)then
               i0=i0box(i,j)
               i1=i1box(i,j)
               j0=j0box(i,j)
               j1=j1box(i,j)
               rx0=rxbox(i,j)
               rx1=1.0-rx0
               ry0=rybox(i,j)
               ry1=1.0-ry0
               uu= 
     *            ry1*(rx1*UB(i0,j0)+rx0*UB(i1,j0))+
     *            ry0*(rx1*UB(i0,j1)+rx0*UB(i1,j1))
               uu=uu*uscale+uoffset

               vv=
     *            ry1*(rx1*VB(i0,j0)+rx0*VB(i1,j0))+
     *            ry0*(rx1*VB(i0,j1)+rx0*VB(i1,j1))
               vv=vv*vscale+voffset

               if(Y(i,j).le.60.0)then
                  UP(i,j)=uu
                  VP(i,j)=vv
	       else
	          sn=sin(X(i,j)*pp)
		  cs=cos(X(i,j)*pp)
                  if(i.ne.1.and.i.ne.imt)then
                     xx1=(90.0-Y(i+1,j))*sin(X(i+1,j)*pp)
                     yy1=-(90.0-Y(i+1,j))*cos(X(i+1,j)*pp)
                     xx2=(90.0-Y(i-1,j))*sin(X(i-1,j)*pp)
                     yy2=-(90.0-Y(i-1,j))*cos(X(i-1,j)*pp)
                  else
                     if(i.eq.1)then
                        xx1=(90.0-Y(i+1,j))*sin(X(i+1,j)*pp)
                        yy1=-(90.0-Y(i+1,j))*cos(X(i+1,j)*pp)
                        xx2=(90.0-Y(i,j))*sn
                        yy2=-(90.0-Y(i,j))*cs
                     else
                        xx1=(90.0-Y(i,j))*sn
                        yy1=-(90.0-Y(i,j))*cs
                        xx2=(90.0-Y(i-1,j))*sin(X(i-1,j)*pp)
                        yy2=-(90.0-Y(i-1,j))*cos(X(i-1,j)*pp)
                     endif
                  endif
                  dx=xx1-xx2
                  dy=yy1-yy2
                  alph=atan2(dy,dx)
		  uu0=uu*cs-vv*sn
		  vv0=uu*sn+vv*cs
                  UP(i,j)=uu0*cos(alph)+vv0*sin(alph)
                  VP(i,j)=-uu0*sin(alph)+vv0*cos(alph)
	       endif
            else
               UP(i,j)=0.0
               VP(i,j)=0.0
               sn=sin(X(i,j)*pp)
               cs=cos(X(i,j)*pp)
               if(i.ne.1.and.i.ne.imt)then
                  xx1=(90.0-Y(i+1,j))*sin(X(i+1,j)*pp)
                  yy1=-(90.0-Y(i+1,j))*cos(X(i+1,j)*pp)
                  xx2=(90.0-Y(i-1,j))*sin(X(i-1,j)*pp)
                  yy2=-(90.0-Y(i-1,j))*cos(X(i-1,j)*pp)
               else
                  if(i.eq.1)then
                     xx1=(90.0-Y(i+1,j))*sin(X(i+1,j)*pp)
                     yy1=-(90.0-Y(i+1,j))*cos(X(i+1,j)*pp)
                     xx2=(90.0-Y(i,j))*sn
                     yy2=-(90.0-Y(i,j))*cs
                  else
                     xx1=(90.0-Y(i,j))*sn
                     yy1=-(90.0-Y(i,j))*cs
                     xx2=(90.0-Y(i-1,j))*sin(X(i-1,j)*pp)
                     yy2=-(90.0-Y(i-1,j))*cos(X(i-1,j)*pp)
                  endif
               endif
               dx=xx1-xx2
               dy=yy1-yy2
               alph=atan2(dy,dx)
	       uu=0.0
	       vv=0.0
               j0=j0box(i,j)
               j1=j1box(i,j)
               rr=1.0/dble(IDAT)
               if(flat)then
                  do ii=1,IDAT
                     uur=UB(ii,j0)*uscale+uoffset
                     vvr=VB(ii,j0)*vscale+voffset
                     uu=uu+
     *rr*(uur*cos(pp*Xdat(ii))-vvr*sin(pp*Xdat(ii)))
                     vv=vv+
     *rr*(uur*sin(pp*Xdat(ii))+vvr*cos(pp*Xdat(ii)))
                  enddo
               else

c        Components at the pole as mean values
                  if(first)then
                     first=.false.
                     UPOLE=0.0
                     VPOLE=0.0
                     do ii=1,IDAT
                        uur=UB(ii,j0)*uscale+uoffset
                        vvr=VB(ii,j0)*vscale+voffset
                        UPOLE=UPOLE+
     *                     rr*(uur*cos(pp*Xdat(ii))-vvr*sin(pp*Xdat(ii)))
                        VPOLE=VPOLE+
     *                     rr*(uur*sin(pp*Xdat(ii))+vvr*cos(pp*Xdat(ii)))
                     enddo
                  endif

                  i0=i0box(i,j)
                  i1=i1box(i,j)

c                  Write(nu_diag,*) i,j
c                  Write(nu_diag,*) i0,j0,i1,j1

                  rx0=rxbox(i,j)
                  rx1=1.0-rx0
                  ry0=rybox(i,j)
                  ry1=1.0-ry0

c                  Write(nu_diag,*) 'rx0=',rx0,
c     *                         ',   rx1=',rx1

c        Rotate data velocity vector at i0,j0 to make it cartesian
                  uu0=UB(i0,j0)*uscale+uoffset
                  vv0=VB(i0,j0)*vscale+voffset
                  uu=uu0*cos(pp*Xdat(i0))-vv0*sin(pp*Xdat(i0))
                  vv=uu0*sin(pp*Xdat(i0))+vv0*cos(pp*Xdat(i0))
                  uu0=uu
                  vv0=vv
c                  Write(nu_diag,*) 'u0=',uu0,',   v0=',vv0

c        Rotate data velocity vector at i1,j0 to make it cartesian
                  uu1=UB(i1,j0)*uscale+uoffset
                  vv1=VB(i1,j0)*vscale+voffset
                  uu=uu1*cos(pp*Xdat(i1))-vv1*sin(pp*Xdat(i1))
                  vv=uu1*sin(pp*Xdat(i1))+vv1*cos(pp*Xdat(i1))
                  uu1=uu
                  vv1=vv
c                  Write(nu_diag,*) 'u1=',uu1,',   v1=',vv1

c        Interpolate along longitude between i0,j0 and i1,j0
                  uu=rx1*uu0+rx0*uu1
                  vv=rx1*vv0+rx0*vv1

c        Apply offset and scaling factor
c                  uu=uu*uscale+uoffset
c                  vv=vv*vscale+voffset

c                  Write(nu_diag,*) 'u=',uu,',   v=',vv

c        Interpolate along latitude between j0 and pole,
c        considering pole value to be (UPOLE,VPOLE)
                  uu=ry1*uu+ry0*UPOLE
                  vv=ry1*vv+ry0*VPOLE

c                  Write(nu_diag,*) 'ry0=',ry0,
c     *                         ',   ry1=',ry1
c                  Write(nu_diag,*) 'U=',uu,',   V=',vv

               endif
               UP(i,j)=uu*cos(alph)+vv*sin(alph)
               VP(i,j)=-uu*sin(alph)+vv*cos(alph)
            endif
            PPU=sqrt(UP(i,j)**2+VP(i,j)**2)
	    if(PMAX.lt.PPU)PMAX=PPU
	    if(PMIN.gt.PPU)PMIN=PPU
         enddo
      enddo

c      Write(nu_diag,*)'Max wind: ', PMAX
c      Write(nu_diag,*)'Min wind: ', PMIN
      
      endif
      
      nold=n
      
      do j=1,jmt
         do i=1,imt
	    UBOT(i,j)=UP(i,j)
	    VBOT(i,j)=VP(i,j)
         enddo
      enddo

      
      firstcall=.false.

      return
      end

      subroutine airtemp(idate, imt, jmt, X, Y, mask, PBOT)
      
      include 'netcdf.inc'
      
      integer*4, intent(in) :: idate, imt, jmt, mask(imt,jmt)
      real*8, intent(in) :: X(imt,jmt), Y(imt,jmt)
      real*8, intent(out) :: PBOT(imt,jmt)
      real*8, allocatable, dimension(:,:) ::
     &  rxbox,rybox,PP
      integer, allocatable, dimension(:,:) ::
     &  i0box,j0box, i1box,j1box, nodetype
      integer dindex(12)
      integer i,j,n,ncid,vid,did,beg(3),cnt(3),ierr,nold
      integer dday, dmonth, dyear, dyearold
      character*80 attstr
      character*80 fname
      character, allocatable :: str(:)
c      integer,parameter :: IDAT=192
c      integer,parameter :: JDAT=94
      integer,parameter :: IDAT=144
      integer,parameter :: JDAT=73
      real*8,parameter :: XSTEPDAT=360.0/IDAT
      real*8,parameter :: YSTEPDAT=2*(90.0/(JDAT+1))
      real*8,parameter :: pi=3.14159265/180.0
      real*8 ::  Xdat(IDAT), Ydat(JDAT)
      real*4 ::  Xdat0(IDAT), Ydat0(JDAT)
      Integer*2,dimension(IDAT,JDAT) ::  PB0
      Real*4,dimension(IDAT,JDAT) ::  PB4
      Real*4 attrib
      Real*8 offset,scale
      Real*8,dimension(IDAT,JDAT) ::  PB
      character*80 input_path
      common/input_data/input_path
      character*4 cyear

      save

      logical firstcall, newyear, flat
      integer nu_diag
      common/iodiag/nu_diag
      save firstcall
      data firstcall/.true./,nold/-1/,dyearold/-1/,dindex
     */0,31,59,90,120,151,181,212,243,273,304,334/

      dyear=int(idate/10000)
      newyear=dyearold.ne.dyear
      
      flat=.false.
      
      dmonth=int((idate-10000*dyear)/100)
      dday=int(idate-10000*dyear-100*dmonth)
      
      n=dindex(dmonth)+dday
c      if(4*int(dyear/4).eq.4.and.dmonth.gt.2)n=n+1
      
c      Write(nu_diag,*)'In AirTemp'
      
	if(firstcall)then
	
      allocate(rxbox(imt,jmt),rybox(imt,jmt),i0box(imt,jmt),
     *   j0box(imt,jmt),i1box(imt,jmt),j1box(imt,jmt),nodetype(imt,jmt),
     *   PP(imt,jmt))

         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'pottmp995/pottmp.sig995.'//cyear//'.nc'
         ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'AirTemp Error: reading NetCDF file'
            call exit_coupler
         endif

      beg(1)=1
      beg(2)=1
      beg(3)=1
      cnt(1)=IDAT
      cnt(2)=1
      cnt(3)=1
      
      ierr=nf_inq_varid(ncid,'lon',vid)
      ierr=nf_get_vara_real(ncid,vid,beg,cnt,Xdat0(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'AirTemp Error: Error reading lon'
         call exit_coupler
      endif

      cnt(1)=JDAT

      ierr=nf_inq_varid(ncid,'lat',vid)
      ierr=nf_get_vara_real(ncid,vid,beg,cnt,Ydat0(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'AirTemp Error: Error reading lat'
         call exit_coupler
      endif
      ierr=nf_close(ncid)

	   do i=1,IDAT
	      Xdat(i)=dble(Xdat0(i))
	   enddo
	   do j=1,JDAT
	      Ydat(j)=dble(Ydat0(JDAT-j+1))
           enddo

	   do i=1,imt
	      do j=1,jmt
	         xx=X(i,j)
		 if(xx.lt.0)xx=xx+360.0
		 yy=Y(i,j)
		 mm=mask(i,j)
		    i0=0
		    do ii=1,IDAT-1
		       if(Xdat(ii).le.xx.and.Xdat(ii+1).ge.xx)then
		          i0=ii
			  i1=ii+1
			  rxbox(i,j)=(xx-Xdat(i0))/(Xdat(i1)-Xdat(i0))
			  exit
		       endif
		    enddo
		    if(i0.eq.0)then
		       if(Xdat(IDAT).le.xx.and.(Xdat(1)+360).ge.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0))/
     *  		    (Xdat(i1)+360-Xdat(i0))
                   elseif(Xdat(1).ge.xx.and.(Xdat(IDAT)-360).le.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0)+360)/
     *  		    (Xdat(i1)+360-Xdat(i0))
		       else
                          Write(nu_diag,*) 'Airtemp: Grid from NetCDF is not covering'
			  Write(nu_diag,*) xx,Xdat(1),Xdat(IDAT)
                          call exit_coupler
		       endif
		    endif
		    j0=0
		    if(yy.le.Ydat(2))then
		       j0=2
		       j1=2
		       nodetype(i,j)=1
                       rybox(i,j)=(yy+90.0)/(Ydat(j1)+90.0)
		    else if(yy.ge.Ydat(JDAT-1))then
		       j0=JDAT-1
		       j1=JDAT-1
		       nodetype(i,j)=1
                       rybox(i,j)=(yy-Ydat(j0))/(90.0-Ydat(j0))
		    else
		    do jj=1,JDAT-1
		       if(Ydat(jj).le.yy.and.Ydat(jj+1).ge.yy)then
		          j0=jj
			  j1=jj+1
			  rybox(i,j)=(yy-Ydat(j0))/(Ydat(j1)-Ydat(j0))
			  exit
		       endif
		    enddo
		    nodetype(i,j)=0
		    endif
		    i0box(i,j)=i0
		    i1box(i,j)=i1
		    j0box(i,j)=j0
		    j1box(i,j)=j1
	      enddo
	   enddo
	   
        endif
	
      if(n.ne.nold)then
      
      if(newyear)then

         if(dyearold.ge.0)then      
	   ierr=nf_close(ncid)
	 endif
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'pottmp995/pottmp.sig995.'//cyear//'.nc'
c	 Write(nu_diag,*)'Reading file: ',fname
	   ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
	   if(ierr.NE.NF_NOERR)then
              Write(nu_diag,*) 'AirTemp Error: reading NetCDF file'
              call exit_coupler
c           else
c              Write(nu_diag,'(a,i3)')'Open CDF ', ncid
           endif
	   dyearold=dyear
      endif

c      Write(nu_diag,'(a,i3)')'Read record no.',n,' from CDF ', ncid

      beg(1)=1
      beg(2)=1
      beg(3)=n
      cnt(1)=IDAT
      cnt(2)=JDAT
      cnt(3)=1
      
c      ierr=nf_inq_varid(ncid,'air',vid)
      ierr=nf_inq_varid(ncid,'pottmp',vid)
      if(dyear.lt.2014)then
         ierr=nf_get_att_real(ncid,vid,'add_offset',attrib)
         offset=dble(attrib)
         ierr=nf_get_att_real(ncid,vid,'scale_factor',attrib)
         scale=dble(attrib)
c         Write(nu_diag,*)'Offset=',offset
c         Write(nu_diag,*)'Factor=',scale
         ierr=nf_get_vara_int2(ncid,vid,beg,cnt,PB0(1,1))
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'AirTemp Error: Error reading air'
            call exit_coupler
         endif
      
         Do i=1,IDAT
            Do j=1,JDAT
               PB(i,j)=dble(PB0(i,JDAT-j+1))
            Enddo
         Enddo
      else
         offset=0.0
         scale=1.0
         ierr=nf_get_vara_real(ncid,vid,beg,cnt,PB4(1,1))
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'AirTemp Error: Error reading air'
            call exit_coupler
         endif
      
         Do i=1,IDAT
            Do j=1,JDAT
               PB(i,j)=dble(PB4(i,JDAT-j+1))
            Enddo
         Enddo
      endif

c      Write(nu_diag,'(a,f6.2,a)')'Read Control: TA(115,31)= ',
c     *   PB(115,31)*scale+offset,' K'

      PMIN=1e+10
      PMAX=-1e+10
      do j=1,jmt
         do i=1,imt
             i0=i0box(i,j)
             i1=i1box(i,j)
             j0=j0box(i,j)
             j1=j1box(i,j)
             rx0=rxbox(i,j)
             rx1=1.0-rx0
             ry0=rybox(i,j)
             ry1=1.0-ry0
	       if(nodetype(i,j).eq.0)then
               PP(i,j)= scale*(
     *            ry1*(rx1*PB(i0,j0)+rx0*PB(i1,j0))+
     *            ry0*(rx1*PB(i0,j1)+rx0*PB(i1,j1))) +
     *            offset
              else
               PP(i,j)=0.0
	       rr=1.0/IDAT
               if(flat)then
	       do ii=1,IDAT
                  PP(i,j)=PP(i,j)+rr*PB(ii,j0)
	       enddo
	       PP(i,j)=scale*PP(i,j)+offset
               else
c        Components at the pole as mean values
                  PPOLE=0.0
                  do ii=1,IDAT
                     uur=PB(ii,j0)*scale+offset
                     PPOLE=PPOLE+rr*uur
                  enddo
c      Write(nu_diag,*)i,j,'Tpole= ', PPOLE

                  uu0=PB(i0,j0)*scale+offset
                  uu1=PB(i1,j0)*scale+offset

                  uu=rx1*uu0+rx0*uu1

                  PP(i,j)=ry1*uu+ry0*PPOLE
               endif
            endif
	    if(PMAX.lt.PP(i,j))PMAX=PP(i,j)
	    if(PMIN.gt.PP(i,j))PMIN=PP(i,j)
         enddo
      enddo

c      Write(nu_diag,*)'Max temperature: ', PMAX
c      Write(nu_diag,*)'Min temperature: ', PMIN
c      Write(nu_diag,*)'Pole temperature: ', PP(88,136)-273.16
      
      endif
      
      nold=n
      
      do j=1,jmt
         do i=1,imt
	    PBOT(i,j)=PP(i,j)
         enddo
      enddo

      
      firstcall=.false.

      return
      end

      real*8 function dist(x1,y1,x0,y0)
      dist=6371000.0*acos(
     *   cos(y1)*cos(x1)*cos(y0)*cos(x0) +
     *   cos(y1)*sin(x1)*cos(y0)*sin(x0) +
     *   sin(y1)*sin(y0)
     *)
      return
      end

      subroutine correct_wind(idate, imt, jmt, X, Y, mask, UBOT, VBOT)
      
      include 'netcdf.inc'
      
      integer*4, intent(in) :: idate, imt, jmt, mask(imt,jmt)
      real*8, intent(in) :: X(imt,jmt), Y(imt,jmt)
      real*8, intent(out) :: UBOT(imt,jmt), VBOT(imt,jmt)
      real*8, allocatable, dimension(:,:) ::
     &  UP,VP
      real*8, allocatable, dimension(:,:,:) :: rnode
      integer, allocatable, dimension(:,:,:) ::
     &  inode,jnode
      integer dindex(12),monlen(12),ind(2)
      integer i,j,n,ncidu,ncidv,vid,did,beg(3),cnt(3),ierr,nold
      integer dday, dmonth, dyear, dyearold
      character*80 attstr
      character*80 fname
      character, allocatable :: str(:)
      integer,parameter :: IDAT=360
      integer,parameter :: NNO=8
      integer,parameter :: JDAT=170
      real*8,parameter :: XSTEPDAT=1.0
      real*8,parameter :: YSTEPDAT=1.0
      real*8,parameter :: pi=3.14159265/180.0
      real*8 ::  Xdat(IDAT,JDAT), Ydat(IDAT,JDAT)
      Real*4,dimension(IDAT,JDAT) ::  UB0, VB0, UB1, VB1, RRad
      Real*4 uoffset,voffset
      real*4 uu0,vv0
      real*8 uu,vv,Rru,Rrv
      character*80 input_path
      common/input_data/input_path

      save

      logical firstcall, newyear, leepyear
      integer nu_diag
      common/iodiag/nu_diag
      save firstcall
      data firstcall/.true./,nold/-1/,dyearold/-1/,dindex
     */0,31,59,90,120,151,181,212,243,273,304,334/,monlen
     */31,28,31,30,31,30,31,31,30,31,30,31/

      dyear=int(idate/10000)
c      leepyear=(dyear.eq.4*(int(float(dyear)/4.)))
      leepyear=.false.
      
      dmonth=int((idate-10000*dyear)/100)
      
      dday=int(idate-10000*dyear-100*dmonth)
      
      n=dindex(dmonth)+dday
      
      if (dday.gt.15) then
         n0=dmonth
         n1=dmonth+1
         if(n1.gt.12)n1=n1-12
         leep=0
         if(n0.eq.2.and.leepyear)leep=1
         dr=float(monlen(n0)+leep)
         ratio=float(dday-15)/dr
      else
         n1=dmonth
         n0=dmonth-1
         if(n0.lt.1)n0=n0+12
         leep=0
         if(n0.eq.2.and.leepyear)leep=1
         dr=float(monlen(n0)+leep)
         ratio=float(monlen(n0)+leep+dday-15)/dr
      endif
      
c      if(4*int(dyear/4).eq.4.and.dmonth.gt.2)n=n+1
      
c      Write(nu_diag,*)'In Correct_Wind'
      
	if(firstcall)then
	
      allocate(inode(imt,jmt,nno),jnode(imt,jmt,nno),rnode(imt,jmt,nno),
     *   UP(imt,jmt),VP(imt,jmt))

      do i=1,IDAT
      do j=1,JDAT
	 Xdat(i,j)=-179.5+(i-1)*XSTEPDAT
         Ydat(i,j)=-84.5+(j-1)*YSTEPDAT
      enddo
      enddo
      
      RRad0=1.5*max(XSTEPDAT,YSTEPDAT)*2
	   
	   do i=1,imt
c              Write(nu_diag,'(a,i3,a,i3)')'Layer ',i,' of ',imt
	      do j=1,jmt
	         xx=X(i,j)
		 yy=Y(i,j)
                 if(mask(i,j).gt.0)then
                    RRad=sqrt((xx-Xdat)**2+(yy-Ydat)**2)
                    RR0=0.0
                    do l=1,nno
                       ind=minloc(RRad,RRad.gt.RR0)
                       inode(i,j,l)=ind(1)
                       jnode(i,j,l)=ind(2)
                       RR0=RRad(ind(1),ind(2))
                       RRad(ind(1),ind(2))=RRad(ind(1),ind(2))-1.0
                       rnode(i,j,l)=exp(-(RR0/RRad0)**2)
                    enddo
                 endif
	      enddo
	   enddo

         fname=trim(input_path)//'wndstress/tauu.nc'
c	 Write(nu_diag,*)'Reading file: ',fname
	   ierr=nf_open(trim(fname),NF_NOWRITE,ncidu)
	   if(ierr.NE.NF_NOERR)then
              Write(nu_diag,*) 'U-Wind Stress Error: reading NetCDF file'
              call exit_coupler
c           else
c              Write(nu_diag,'(a,i3)')'Open CDF ', ncidu
           endif

         fname=trim(input_path)//'wndstress/tauv.nc'
c	 Write(nu_diag,*)'Reading file: ',fname
	   ierr=nf_open(trim(fname),NF_NOWRITE,ncidv)
	   if(ierr.NE.NF_NOERR)then
              Write(nu_diag,*) 'V-Wind Stress Error: reading NetCDF file'
              call exit_coupler
c           else
c              Write(nu_diag,'(a,i3)')'Open CDF ', ncidv
           endif

      endif

      if(n.ne.nold)then
      
c      Write(nu_diag,'(a,i3,a,i3,a,i3)')'Read record no.',n0,' and ',n1,
c     *   ' from CDF ', ncidu
c      Write(nu_diag,'(a,f5.3)')'Ratio coeff',ratio

      beg(1)=1
      beg(2)=1
      beg(3)=n0
      cnt(1)=IDAT
      cnt(2)=JDAT
      cnt(3)=1
      
      ierr=nf_inq_varid(ncidu,'tauu',vid)
      ierr=nf_get_att_real(ncidu,vid,'missing_value',uoffset)
      ierr=nf_get_vara_real(ncidu,vid,beg,cnt,UB0(1,1))
      beg(3)=n1
      ierr=nf_get_vara_real(ncidu,vid,beg,cnt,UB1(1,1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'U-Wind Error: Error reading uwnd'
         call exit_coupler
      endif

      beg(1)=1
      beg(2)=1
      beg(3)=n0
      cnt(1)=IDAT
      cnt(2)=JDAT
      cnt(3)=1
      
      ierr=nf_inq_varid(ncidv,'tauv',vid)
      ierr=nf_get_att_real(ncidv,vid,'missing_value',voffset)
      ierr=nf_get_vara_real(ncidv,vid,beg,cnt,VB0(1,1))
      beg(3)=n1
      ierr=nf_get_vara_real(ncidv,vid,beg,cnt,VB1(1,1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'V-Wind Error: Error reading vwnd'
         call exit_coupler
      endif

c      Write(nu_diag,'(a,f6.2,a,f6.2,a)')'Read Control: taux(141,134)= ',
c     *   UB0(141,134),' N/m2,  tauy(141,134)= ',
c     *   VB0(141,134),' N/m2'
     
      do i=1,IDAT
         do j=1,JDAT
            if((UB0(i,j).ne.uoffset).and.(UB1(i,j).ne.uoffset))then
               UB0(i,j)=UB0(i,j)*(1.0-ratio)+UB1(i,j)*ratio
            else
               UB0(i,j)=uoffset
            endif
            if((VB0(i,j).ne.voffset).and.(VB1(i,j).ne.voffset))then
               VB0(i,j)=VB0(i,j)*(1.0-ratio)+VB1(i,j)*ratio
            else
               VB0(i,j)=voffset
            endif
         enddo
      enddo

      PMIN=1e+10
      PMAX=-1e+10
      pp=atan2(1.0,1.0)*4.0/180.0
c      pp=1.0

      do j=1,jmt
         do i=1,imt
            if(mask(i,j).gt.0)then
               uu=0.0
               vv=0.0
               Rru=0.0
               Rrv=0.0
               Do l=1,nno
                  uu0=UB0(inode(i,j,l),jnode(i,j,l))
                  vv0=VB0(inode(i,j,l),jnode(i,j,l))
                  if(uu0.ne.uoffset)then
                     uu=uu+dble(rnode(i,j,l)*uu0)
                     Rru=Rru+dble(rnode(i,j,l))
                  endif
                  if(vv0.ne.voffset)then
                     vv=vv+dble(rnode(i,j,l)*vv0)
                     Rrv=Rrv+dble(rnode(i,j,l))
                  endif
               enddo
               if(Rru*Rrv.ne.0.0)then
                  uu=uu/Rru
                  vv=vv/Rrv

               if(Y(i,j).le.60.0)then
                  UP(i,j)=uu
                  VP(i,j)=vv
	       else if(Y(i,j).lt.65.0)then
	          sn=sin(X(i,j)*pp)
		  cs=cos(X(i,j)*pp)
                  if(i.ne.1.and.i.ne.imt)then
                     xx1=(90-Y(i+1,j))*sin(X(i+1,j)*pp)
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
		  Rru=uu*cs-vv*sn
		  Rrv=uu*sn+vv*cs
                  uu=Rru*cos(alph)+Rrv*sin(alph)
                  vv=-Rru*sin(alph)+Rrv*cos(alph)
                  UP(i,j)=(uu*(65.-Y(i,j))+UBOT(i,j)*(Y(i,j)-60.))/5.
                  VP(i,j)=(vv*(65.-Y(i,j))+VBOT(i,j)*(Y(i,j)-60.))/5.
               else
                  UP(i,j)=UBOT(i,j)
                  VP(i,j)=VBOT(i,j)
	       endif
               endif
            endif
            PPU=sqrt(UP(i,j)**2+VP(i,j)**2)
	    if(PMAX.lt.PPU)PMAX=PPU
	    if(PMIN.gt.PPU)PMIN=PPU
         enddo
      enddo

c      Write(nu_diag,*)'Max wind stress: ', PMAX
c      Write(nu_diag,*)'Min wind stress: ', PMIN
      
      endif
      
      nold=n
      
      do j=1,jmt
         do i=1,imt
	    UBOT(i,j)=UP(i,j)
	    VBOT(i,j)=VP(i,j)
         enddo
      enddo

      
      firstcall=.false.

      return
      end

      subroutine recorrect_wind(idate, imt, jmt, X, Y, mask, UBOT, VBOT)
      
      include 'netcdf.inc'
      
      integer*4, intent(in) :: idate, imt, jmt, mask(imt,jmt)
      real*8, intent(in) :: X(imt,jmt), Y(imt,jmt)
      real*8, intent(out) :: UBOT(imt,jmt), VBOT(imt,jmt)
      real*8, allocatable, dimension(:,:) ::
     &  UP,VP
      integer dindex(12),monlen(12),ind(2)
      integer i,j,n,ncidu,ncidv,vid,did,beg(3),cnt(3),ierr,nold
      integer dday, dmonth, dyear, dyearold
      character*80 attstr
      character*80 fname
      character, allocatable :: str(:)
      integer,parameter :: IDAT=360
      integer,parameter :: NNO=8
      integer,parameter :: JDAT=170
      real*8,parameter :: XSTEPDAT=1.0
      real*8,parameter :: YSTEPDAT=1.0
      real*8,parameter :: pi=3.14159265/180.0
      real*8 ::  Xdat(IDAT,JDAT), Ydat(IDAT,JDAT)
      Real*8, allocatable, dimension(:,:) ::  UB0, VB0, UB1, VB1
      Real*4 uoffset,voffset
      real*4 uu0,vv0
      real*8 uu,vv,Rru,Rrv
      character*80 input_path
      common/input_data/input_path

      save

      logical firstcall, newyear, leepyear
      logical, external :: isnan8
      integer nu_diag
      common/iodiag/nu_diag
      save firstcall
      data firstcall/.true./,nold/-1/,dyearold/-1/,dindex
     */0,31,59,90,120,151,181,212,243,273,304,334/,monlen
     */31,28,31,30,31,30,31,31,30,31,30,31/

      dyear=int(idate/10000)
c      leepyear=(dyear.eq.4*(int(float(dyear)/4.)))
      leepyear=.false.
            
      dmonth=int((idate-10000*dyear)/100)
      
      dday=int(idate-10000*dyear-100*dmonth)
      
      n=dindex(dmonth)+dday
      
      if (dday.gt.15) then
         n0=dmonth
         n1=dmonth+1
         if(n1.gt.12)n1=n1-12
         leep=0
         if(n0.eq.2.and.leepyear)leep=1
         dr=float(monlen(n0)+leep)
         ratio=float(dday-15)/dr
      else
         n1=dmonth
         n0=dmonth-1
         if(n0.lt.1)n0=n0+12
         leep=0
         if(n0.eq.2.and.leepyear)leep=1
         dr=float(monlen(n0)+leep)
         ratio=float(monlen(n0)+leep+dday-15)/dr
      endif
c      n0=dmonth
c      n1=n0+1
c      if(n1.gt.12)n1=n1-12            
c      dr=float(monlen(n0))      
c      ratio=float(dday-1)/dr
c      if(4*int(dyear/4).eq.4.and.dmonth.gt.2)n=n+1
      
c      Write(nu_diag,*)'In ReCorrect_Wind'
      
	if(firstcall)then
	
      allocate(UP(imt,jmt),VP(imt,jmt),UB0(imt,jmt),VB0(imt,jmt),
     * UB1(imt,jmt),VB1(imt,jmt))

         fname=trim(input_path)//'wndstress/tax.dat'
c	 Write(nu_diag,*)'Reading file: ',fname
#if COMPAQ
         open (41,file=trim(fname),form='unformatted',access='direct',
     *      recl=imt*jmt*2)
#else
         open (41,file=trim(fname),form='unformatted',access='direct',
     *      recl=imt*jmt*8)
#endif

         fname=trim(input_path)//'wndstress/tay.dat'
c	 Write(nu_diag,*)'Reading file: ',fname
#if COMPAQ
         open (42,file=trim(fname),form='unformatted',access='direct',
     *      recl=imt*jmt*2)
#else
         open (42,file=trim(fname),form='unformatted',access='direct',
     *      recl=imt*jmt*8)
#endif

      endif

      if(n.ne.nold)then
      
c      Write(nu_diag,'(a,i3,a,i3,a,i3)')'Read record no.',n0,' and ',n1,
c     *   ' from CDF ', ncidu
c      Write(nu_diag,'(a,f5.3)')'Ratio coeff',ratio

      Read(41,rec=n0) UB0
      Read(41,rec=n1) UB1
      
      Read(42,rec=n0) VB0
      Read(42,rec=n1) VB1

      do i=1,imt
         do j=1,jmt
            if(isnan8(UB0(i,j)))UB0(i,j)=0.0
            if(isnan8(VB0(i,j)))VB0(i,j)=0.0
            if(isnan8(UB1(i,j)))UB1(i,j)=0.0
            if(isnan8(VB1(i,j)))VB1(i,j)=0.0
            UB0(i,j)=UB0(i,j)*(1.0-ratio)+UB1(i,j)*ratio
            VB0(i,j)=VB0(i,j)*(1.0-ratio)+VB1(i,j)*ratio
         enddo
      enddo

c      Write(nu_diag,'(a,f6.2,a,f6.2,a)')'Read Control: taux(141,134)= ',
c     *   UB0(141,134),' N/m2,  tauy(141,134)= ',
c     *   VB0(141,134),' N/m2'
     


      PMIN=1e+10
      PMAX=-1e+10
      do j=1,jmt
         do i=1,imt
            if(mask(i,j).gt.0)then
               uu=UB0(i,j)*0.1
               vv=VB0(i,j)*0.1
               if(Y(i,j).le.60.0)then
                  UP(i,j)=uu
                  VP(i,j)=vv
	       else if(Y(i,j).lt.65.0)then
                  UP(i,j)=(uu*(65.-Y(i,j))+UBOT(i,j)*(Y(i,j)-60.))/5.
                  VP(i,j)=(vv*(65.-Y(i,j))+VBOT(i,j)*(Y(i,j)-60.))/5.
               else
                  UP(i,j)=UBOT(i,j)
                  VP(i,j)=VBOT(i,j)
               endif
            endif
            PPU=sqrt(UP(i,j)**2+VP(i,j)**2)
	    if(PMAX.lt.PPU)PMAX=PPU
	    if(PMIN.gt.PPU)PMIN=PPU
         enddo
      enddo

c      Write(nu_diag,*)'Max wind stress: ', PMAX
c      Write(nu_diag,*)'Min wind stress: ', PMIN
      
      endif
      
      nold=n
      
      do j=1,jmt
         do i=1,imt
	    UBOT(i,j)=UP(i,j)
	    VBOT(i,j)=VP(i,j)
         enddo
      enddo

      
      firstcall=.false.

      return
      end

      subroutine precip(idate, imt, jmt, X, Y, mask, PBOT, 
     *   use_xie,use_xie4ncep,ciafvers, oldstuff)
      
      include 'netcdf.inc'
      
      integer*4, intent(in) :: imt, jmt, mask(imt,jmt), idate, ciafvers
      logical use_xie,use_xie4ncep
      real*8, intent(in) :: X(imt,jmt), Y(imt,jmt)
      real*8, intent(out) :: PBOT(imt,jmt,12)
      real*8, allocatable, dimension(:,:) ::
     &  rxbox,rybox,PP
      integer, allocatable, dimension(:,:) ::
     &  i0box,j0box, i1box,j1box, nodetype
      integer i,j,n,ncid,vid,did,beg(3),cnt(3),ierr,nold
      character*80 attstr
      character*80 fname
      character, allocatable :: str(:)
      integer,parameter :: IDAT=192
      integer,parameter :: JDAT=94
      real*8,parameter :: pi=3.14159265/180.0
      real*8 ::  Xdat(IDAT), Ydat(JDAT)
      real*4 ::  Xdat0(IDAT), Ydat0(JDAT)
      real*4,dimension(IDAT,JDAT) ::  PB, PB0
      Real*4 offset,scale
      character*80 input_path
      common/input_data/input_path
      character*4 cyear

      logical flat, newyear, oldstuff, firstcall
      integer nu_diag, dyear, dyearold
      common/iodiag/nu_diag
      save dyearold
      data dyearold/-1/,firstcall/.true./

      flat=.false.
      dyear=int(idate/10000)
      newyear=dyearold.ne.dyear
      oldstuff=.true.
      
      if((ciafvers.ge.0.and..not.newyear).or.
     *   (.not.firstcall.and.use_xie))return
      
      oldstuff=.false.
      dyearold=dyear
      

c      Write(nu_diag,*)'In Precip'
      
      allocate(rxbox(imt,jmt),rybox(imt,jmt),i0box(imt,jmt),
     *   j0box(imt,jmt),i1box(imt,jmt),j1box(imt,jmt),nodetype(imt,jmt),
     *   PP(imt,jmt))

      if(.not.use_xie)then
      if(ciafvers.ge.0)then
         if(ciafvers.eq.1)then
            if(dyear.lt.1958)dyear=1958
            Write(cyear,'(i4)')dyear
            fname=trim(input_path)//'CIAF_1.0/ncar_precip.'//
     *         cyear//'.2004_08_03.nc'
            if(dyear.gt.2000)
     *         fname=trim(input_path)//'CIAF_1.0/ncar_precip.'//
     *            cyear//'.2006_04_06.nc'
         elseif(ciafvers.eq.2)then
            if(dyear.lt.1948)dyear=1948
            Write(cyear,'(i4)')dyear
            fname=trim(input_path)//'CIAF_2.0/ncar_precip.'//
     *         cyear//'.05APR2010.nc'
         elseif(ciafvers.eq.0)then
            fname=trim(input_path)//'CIAF_2.0/ncar_precip.'//
     *         '15JUNE2009.nc'
         endif
      else
         if(use_xie4ncep)then
            fname=trim(input_path)//'prec/prec.xa.clim.nc'
         else
            fname=trim(input_path)//'prec/prate.mon.mean.nc'
         endif
      endif
      else
         fname=trim(input_path)//'prec/prec.xa.clim.nc'
      endif
c         Write(nu_diag,*) 'Opening NetCDF file '//trim(fname)
         ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'Precipitation Error: reading NetCDF file'
            call exit_coupler
         endif

      beg(1)=1
      beg(2)=1
      beg(3)=1
      cnt(1)=IDAT
      cnt(2)=1
      cnt(3)=1
      
      if(ciafvers.ge.0.and..not.use_xie)then
         ierr=nf_inq_varid(ncid,'LON',vid)
         ierr=nf_get_vara_double(ncid,vid,beg,cnt,Xdat(1))
      else
         ierr=nf_inq_varid(ncid,'lon',vid)
         ierr=nf_get_vara_real(ncid,vid,beg,cnt,Xdat0(1))
         do i=1,IDAT
	    Xdat(i)=dble(Xdat0(i))
         enddo
      endif
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'Precipitation Error: Error reading lon'
         call exit_coupler
      endif

      cnt(1)=JDAT

      if(ciafvers.ge.0.and..not.use_xie)then
         if(ciafvers.eq.1)ierr=nf_inq_varid(ncid,'LAT1',vid)
         if(ciafvers.eq.2.or.ciafvers.eq.0)ierr=nf_inq_varid(ncid,'LAT',vid)
         ierr=nf_get_vara_double(ncid,vid,beg,cnt,Ydat(1))
      else
         ierr=nf_inq_varid(ncid,'lat',vid)
         ierr=nf_get_vara_real(ncid,vid,beg,cnt,Ydat0(1))
         do j=1,JDAT
	    Ydat(j)=dble(Ydat0(j))
         enddo
      endif
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'Precipitation Error: Error reading lat'
         call exit_coupler
      endif
      
      if(Ydat(1).gt.Ydat(JDAT))then
         do j=1,JDAT/2
            YY=Ydat(j)
            Ydat(j)=Ydat(JDAT-j+1)
            Ydat(JDAT-j+1)=YY
         enddo
      endif

	   do i=1,imt
	      do j=1,jmt
	         xx=X(i,j)
		 if(xx.lt.0)xx=xx+360.0
		 yy=Y(i,j)
		 mm=mask(i,j)
		    i0=0
		    do ii=1,IDAT-1
		       if(Xdat(ii).le.xx.and.Xdat(ii+1).ge.xx)then
		          i0=ii
			  i1=ii+1
			  rxbox(i,j)=(xx-Xdat(i0))/(Xdat(i1)-Xdat(i0))
			  exit
		       endif
		    enddo
		    if(i0.eq.0)then
		       if(Xdat(IDAT).le.xx.and.(Xdat(1)+360).ge.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0))/
     *  		    (Xdat(i1)+360-Xdat(i0))
                   elseif(Xdat(1).ge.xx.and.(Xdat(IDAT)-360).le.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0)+360)/
     *  		    (Xdat(i1)+360-Xdat(i0))
		       else
                          Write(nu_diag,*) 'Precip: Grid from NetCDF is not covering'
			  Write(nu_diag,*) xx,Xdat(1),Xdat(IDAT)
                          call exit_coupler
		       endif
		    endif
		    j0=0
		    if(yy.le.Ydat(2))then
		       j0=2
		       j1=2
		       nodetype(i,j)=1
                       rybox(i,j)=(yy+90.0)/(Ydat(j1)+90.0)
		    else if(yy.ge.Ydat(JDAT-1))then
		       j0=JDAT-1
		       j1=JDAT-1
		       nodetype(i,j)=1
                       rybox(i,j)=(yy-Ydat(j0))/(90.0-Ydat(j0))
		    else
		    do jj=1,JDAT-1
		       if(Ydat(jj).le.yy.and.Ydat(jj+1).ge.yy)then
		          j0=jj
			  j1=jj+1
			  rybox(i,j)=(yy-Ydat(j0))/(Ydat(j1)-Ydat(j0))
			  exit
		       endif
		    enddo
		    nodetype(i,j)=0
		    endif
		    i0box(i,j)=i0
		    i1box(i,j)=i1
		    j0box(i,j)=j0
		    j1box(i,j)=j1
	      enddo
	   enddo

      do n=1,12	

c      Write(nu_diag,'(a,i3)')'Read record no.',n,' from CDF ', ncid

      beg(1)=1
      beg(2)=1
      beg(3)=n
      cnt(1)=IDAT
      cnt(2)=JDAT
      cnt(3)=1
      
      scale=1.0
      offset=0.0
      if(.not.use_xie)then
      if(ciafvers.ge.0)then
         ierr=nf_inq_varid(ncid,'RAIN',vid)
         ierr=nf_get_vara_real(ncid,vid,beg,cnt,PB(1,1))
         ierr=nf_inq_varid(ncid,'SNOW',vid)
         ierr=nf_get_vara_real(ncid,vid,beg,cnt,PB0(1,1))
         Do i=1,IDAT
            Do j=1,JDAT
               PB(i,j)=(PB(i,j)+PB0(i,j))/1000.0  ! m/s from mm/s
            Enddo
         Enddo
      else
         if(use_xie4ncep)then
            Write(nu_diag,'(a20,i4)')'Read prec -- recno: ',beg(3)
            ierr=nf_inq_varid(ncid,'prc',vid)
            ierr=nf_get_vara_real(ncid,vid,beg,cnt,PB0(1,1))
            Do i=1,IDAT
               Do j=1,JDAT
                  PB(i,j)=PB0(i,j)/1000.0/86400.0  ! m/s from mm/day
               Enddo
            Enddo
         else
            beg(3)=n+(dyear-1948)*12
            Write(nu_diag,'(a21,i4)')'Read prate -- recno: ',beg(3)
            ierr=nf_inq_varid(ncid,'prate',vid)
            ierr=nf_get_vara_real(ncid,vid,beg,cnt,PB0(1,1))
            Do i=1,IDAT
               Do j=1,JDAT
                  PB(i,j)=PB0(i,JDAT-j+1)/1000.0  ! m/s from kg/m^2/s
               Enddo
            Enddo
         endif
      endif
      else
         Write(nu_diag,'(a19,i4)')'Read Xie -- recno: ',beg(3)
         ierr=nf_inq_varid(ncid,'prc',vid)
         ierr=nf_get_vara_real(ncid,vid,beg,cnt,PB0(1,1))
         Do i=1,IDAT
            Do j=1,JDAT
               PB(i,j)=PB0(i,j)/1000.0/86400.0  ! m/s from mm/day
            Enddo
         Enddo
      endif


c      Write(nu_diag,'(a,f7.2,a)')'Read Control: P(115,31)= ',
c     *   PB(115,31)*8640000.0,' cm/day'

      PMIN=1e+10
      PMAX=-1e+10
      do j=1,jmt
         do i=1,imt
               i0=i0box(i,j)
	       i1=i1box(i,j)
	       j0=j0box(i,j)
	       j1=j1box(i,j)
	       rx0=rxbox(i,j)
	       rx1=1.0-rx0
	       ry0=rybox(i,j)
	       ry1=1.0-ry0
              if(nodetype(i,j).eq.0)then
               PP(i,j)= 
     *            ry1*(rx1*dble(PB(i0,j0))+rx0*dble(PB(i1,j0)))+
     *            ry0*(rx1*dble(PB(i0,j1))+rx0*dble(PB(i1,j1)))
              else
               PP(i,j)=0.0
	       rr=1.0/IDAT
               if(flat)then
	       do ii=1,IDAT
                  PP(i,j)=PP(i,j)+rr*dble(PB(ii,j0))
	       enddo
               else
c        Components at the pole as mean values
                  PPOLE=0.0
                  do ii=1,IDAT
                     uur=dble(PB(ii,j0))*scale+offset
                     PPOLE=PPOLE+rr*uur
                  enddo

                  uu0=dble(PB(i0,j0))
                  uu1=dble(PB(i1,j0))

                  uu=rx1*uu0+rx0*uu1

                  PP(i,j)=ry1*uu+ry0*PPOLE
               endif
            endif
	    if(PMAX.lt.PP(i,j))PMAX=PP(i,j)
	    if(PMIN.gt.PP(i,j))PMIN=PP(i,j)
         enddo
      enddo

c      Write(nu_diag,*)'Max precipitation: ', PMAX*100.0*86400.0
c      Write(nu_diag,*)'Min precipitation: ', PMIN*100.0*86400.0
      
      do j=1,jmt
         do i=1,imt
	    PBOT(i,j,n)=PP(i,j)
         enddo
      enddo

      enddo
      
      ierr=nf_close(ncid)

      return
      end

      subroutine newprecip(imt, jmt, X, Y, mask, PBOT)
      
      integer*4, intent(in) :: imt, jmt, mask(imt,jmt)
      real*8, intent(in) :: X(imt,jmt), Y(imt,jmt)
      real*8 PBOT(imt,jmt,12)
      integer i,j,n,k
      character*80 fname
      integer,parameter :: IDAT=7921
      real*8,parameter :: pi=3.14159265/180.0
      real*8 ::  Xdat(IDAT), Ydat(IDAT)
      real*4 ::  Xdat0(IDAT), Ydat0(IDAT)
      real*4,dimension(IDAT,13) ::  PBA
      real*4,dimension(IDAT) ::  PB, PB0
      real*8, allocatable, dimension(:,:) ::
     &  XX,YY
      real*8 Rad,Val(12),Mes(12)
      character*80 input_path
      common/input_data/input_path

      integer nu_diag
      common/iodiag/nu_diag

c      Write(nu_diag,*)'In New Precip'
      
      allocate(XX(imt,jmt),YY(imt,jmt))

      fname=trim(input_path)//'prec/precip60.out'
      Open(791,file=trim(fname),form='formatted',status='old')

      do i=1,IDAT
         read(791,*) zgrid,Xdat0(i),Ydat0(i),(PBA(i,m),m=1,13)
      end do
      
      Close(791)

      do i=1,IDAT
	 call transcoord(Xdat(i),Ydat(i),dble(Xdat0(i)),dble(Ydat0(i)))
      enddo
      do j=1,jmt
         do i=1,imt
            call transcoord(XX(i,j),YY(i,j),X(i,j),Y(i,j))
         enddo
      enddo

      Rad=2.0
      PMIN=1e+10
      PMAX=-1e+10
      do i=1,imt
         do j=1,jmt
            Val=0.0
            Mes=0.0
            if(mask(i,j).gt.0)then
               do n=1,IDAT
                  nn=n
                  if(n.ne.nn)Write(nu_diag,*)i,j,n
                  RRad=sqrt((XX(i,j)-Xdat(n))**2+(YY(i,j)-Ydat(n))**2)
                  w=Weight2(RRad,Rad)
                  do k=1,12
                     if(PBA(n,k).lt.9999.0)then
                        Val(k)=Val(k)+w*PBA(n,k)
                        Mes(k)=Mes(k)+w
                     endif
                  enddo
               enddo
               do k=1,12
                  PBOT(i,j,k)=Val(k)/Mes(k)
                  PBOT(i,j,k)=PBOT(i,j,k)/1000.0/86400.0/30.0 ! mm/month to m/s
                  if(PMAX.lt.PBOT(i,j,k))PMAX=PBOT(i,j,k)
                  if(PMIN.gt.PBOT(i,j,k))PMIN=PBOT(i,j,k)
               enddo
            else
               do k=1,12
                  PBOT(i,j,k)=0.0
               enddo
            endif

         enddo
      enddo

c      Write(nu_diag,*)'Max new precipitation: ', PMAX*100.0*86400.0
c      Write(nu_diag,*)'Min new precipitation: ', PMIN*100.0*86400.0

      return
      end

      subroutine transcoord(xx,yy,x,y)
c      lam = atan2(
c     1   -( -cos(tet)*cos(tet0)*sin(lam)*cos(lam0) + 
c     2       cos(tet)*cos(tet0)*cos(lam)*sin(lam0) + 
c     3       sin(tet)*sin(tet0) ),
c     9   cos(tet)*( cos(lam)*cos(lam0) + sin(lam)*sin(lam0) ) )
c      z = r*(sin(tet)*cos(tet0) + sin(lam-lam0)*cos(tet)*sin(tet0))
      real*8,parameter :: pi=3.14159265/180.0
      real*8 xx,yy,x,y,lam,lam0,tet,lamp,z,uu,r
      
      lam=x*pi
      tet=y*pi
      r=6.7d+6;
      lam0=-42.0*pi
   
      lamp = datan2( dsin(tet), dcos(lam-lam0)*dcos(tet) )
      z = -r*dsin(lam-lam0)*dcos(tet);
   
      uu=dsqrt(r**2-z**2)
      xx=-z/uu/pi
      yy=lamp/pi
      return
      end
      
      real*8 function Weight2(R,Rad)
      Weight2=dexp(-(R/2.0/Rad)**2)
      return
      end
      
      subroutine clouds(imt, jmt, X, Y, mask, PBOT, power)
      
      include 'netcdf.inc'
      
      integer*4, intent(in) :: imt, jmt, mask(imt,jmt)
      real*8, intent(in) :: X(imt,jmt), Y(imt,jmt), power
      real*8, intent(out) :: PBOT(imt,jmt,12)
      real*8, allocatable, dimension(:,:) ::
     &  rxbox,rybox,PP
      integer, allocatable, dimension(:,:) ::
     &  i0box,j0box, i1box,j1box, nodetype
      integer i,j,n,ncid,vid,did,beg(3),cnt(3),ierr,nold
      character*80 attstr
      character*80 fname
      character, allocatable :: str(:)
      integer,parameter :: IDAT=320
      integer,parameter :: JDAT=160
      real*8,parameter :: pi=3.14159265/180.0
      real*8 ::  Xdat(IDAT), Ydat(JDAT)
      real*4 ::  Xdat0(IDAT), Ydat0(JDAT)
      real*4,dimension(IDAT,JDAT) ::  PB, PB0
      Real*4 offset,scale
      character*80 input_path
      common/input_data/input_path

      logical flat, first
      integer nu_diag
      common/iodiag/nu_diag

      flat=.false.

c      Write(nu_diag,*)'In Clouds'
      
      allocate(rxbox(imt,jmt),rybox(imt,jmt),i0box(imt,jmt),
     *   j0box(imt,jmt),i1box(imt,jmt),j1box(imt,jmt),nodetype(imt,jmt),
     *   PP(imt,jmt))

         fname=trim(input_path)//'cloud/tsm_org.nc'
         ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'Clouds Error: reading NetCDF file'
            Write(nu_diag,*) fname
            call exit_coupler
         endif

      beg(1)=1
      beg(2)=1
      beg(3)=1
      cnt(1)=IDAT
      cnt(2)=1
      cnt(3)=1
      
      ierr=nf_inq_varid(ncid,'longitude_geo_t',vid)
      ierr=nf_get_vara_real(ncid,vid,beg,cnt,Xdat0(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'Clouds Error: Error reading lon'
         call exit_coupler
      endif

      cnt(1)=JDAT

      ierr=nf_inq_varid(ncid,'latitude_geo_t',vid)
      ierr=nf_get_vara_real(ncid,vid,beg,cnt,Ydat0(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'Clouds Error: Error reading lat'
         call exit_coupler
      endif

      do i=1,IDAT
	 Xdat(i)=dble(Xdat0(i))
      enddo
      do j=1,JDAT
	 Ydat(j)=dble(Ydat0(j))
      enddo

	   do i=1,imt
	      do j=1,jmt
	         xx=X(i,j)
		 if(xx.lt.0)xx=xx+360.0
		 yy=Y(i,j)
		 mm=mask(i,j)
		    i0=0
		    do ii=1,IDAT-1
		       if(Xdat(ii).le.xx.and.Xdat(ii+1).ge.xx)then
		          i0=ii
			  i1=ii+1
			  rxbox(i,j)=(xx-Xdat(i0))/(Xdat(i1)-Xdat(i0))
			  exit
		       endif
		    enddo
		    if(i0.eq.0)then
		       if(Xdat(IDAT).le.xx.and.(Xdat(1)+360).ge.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0))/
     *  		    (Xdat(i1)+360-Xdat(i0))
                   elseif(Xdat(1).ge.xx.and.(Xdat(IDAT)-360).le.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0)+360)/
     *  		    (Xdat(i1)+360-Xdat(i0))
		       else
                          Write(nu_diag,*) 'Clouds: Grid from NetCDF is not covering'
                          Write(nu_diag,*) i,j
			  Write(nu_diag,*) xx,Xdat(1),Xdat(IDAT)
                          call exit_coupler
		       endif
		    endif
		    j0=0
		    if(yy.le.Ydat(2))then
		       j0=2
		       j1=2
		       nodetype(i,j)=1
                       rybox(i,j)=(yy+90.0)/(Ydat(j1)+90.0)
		    else if(yy.ge.Ydat(JDAT-1))then
		       j0=JDAT-1
		       j1=JDAT-1
		       nodetype(i,j)=1
                       rybox(i,j)=(yy-Ydat(j0))/(90.0-Ydat(j0))
		    else
		    do jj=1,JDAT-1
		       if(Ydat(jj).le.yy.and.Ydat(jj+1).ge.yy)then
		          j0=jj
			  j1=jj+1
			  rybox(i,j)=(yy-Ydat(j0))/(Ydat(j1)-Ydat(j0))
			  exit
		       endif
		    enddo
		    nodetype(i,j)=0
		    endif
		    i0box(i,j)=i0
		    i1box(i,j)=i1
		    j0box(i,j)=j0
		    j1box(i,j)=j1
	      enddo
	   enddo

      do n=1,12	

c      Write(nu_diag,'(a,i3)')'Read record no.',n,' from CDF ', ncid

      beg(1)=1
      beg(2)=1
      beg(3)=n
      cnt(1)=IDAT
      cnt(2)=JDAT
      cnt(3)=1
      
      ierr=nf_inq_varid(ncid,'cloud_cover',vid)
      ierr=nf_get_vara_real(ncid,vid,beg,cnt,PB0(1,1))

      Do i=1,IDAT
         Do j=1,JDAT
            PB(i,j)=PB0(i,j)
         Enddo
      Enddo

c      Write(nu_diag,'(a,f7.2,a)')'Read Control: C(115,31)= ',
c     *   PB(115,31)*100,' %'

      PMIN=1e+10
      PMAX=-1e+10
      first=.true.
      do j=1,jmt
         do i=1,imt
               i0=i0box(i,j)
	       i1=i1box(i,j)
	       j0=j0box(i,j)
	       j1=j1box(i,j)
	       rx0=rxbox(i,j)
	       rx1=1.0-rx0
	       ry0=rybox(i,j)
	       ry1=1.0-ry0
              if(nodetype(i,j).eq.0)then
               PP(i,j)= 
     *            ry1*(rx1*dble(PB(i0,j0))+rx0*dble(PB(i1,j0)))+
     *            ry0*(rx1*dble(PB(i0,j1))+rx0*dble(PB(i1,j1)))
              else
               PP(i,j)=0.0
	       rr=1.0/IDAT
               if(flat)then
	       do ii=1,IDAT
                  PP(i,j)=PP(i,j)+rr*dble(PB(ii,j0))
	       enddo
               else
c        Components at the pole as mean values
                  if(first)then
                  first=.false.
                  PPOLE=0.0
                  do ii=1,IDAT
                     uur=dble(PB(ii,j0))
                     PPOLE=PPOLE+rr*uur
                  enddo
c                  Write(nu_diag,*)'CL_POLE=',PPOLE
                  endif

                  uu0=dble(PB(i0,j0))
                  uu1=dble(PB(i1,j0))

                  uu=rx1*uu0+rx0*uu1

                  PP(i,j)=ry1*uu+ry0*PPOLE
               endif
            endif
	    if(PMAX.lt.PP(i,j))PMAX=PP(i,j)
	    if(PMIN.gt.PP(i,j))PMIN=PP(i,j)
         enddo
      enddo

c      Write(nu_diag,*)'Max precipitation: ', PMAX*100.0
c      Write(nu_diag,*)'Min precipitation: ', PMIN*100.0
      
      do j=1,jmt
         do i=1,imt
	    PBOT(i,j,n)=PP(i,j)**(power)
         enddo
      enddo

      enddo

      ierr=nf_close(ncid)

      return
      end

      subroutine ciafwind(idate, sec, imt, jmt, X, Y, mask, UBOT, VBOT,
     *                    vers)
      
      include 'netcdf.inc'
      
      integer*4, intent(in) :: idate, imt, jmt, mask(imt,jmt), sec, vers
      real*8, intent(in) :: X(imt,jmt), Y(imt,jmt)
      real*8, intent(out) :: UBOT(imt,jmt), VBOT(imt,jmt)
      real*8, allocatable, dimension(:,:) ::
     &  rxbox,rybox,UP,VP
      integer, allocatable, dimension(:,:) ::
     &  i0box,j0box, i1box,j1box, nodetype
      integer dindex(12)
      integer i,j,n,ncidu,ncidv,vid,did,beg(3),cnt(3),ierr,nold
      integer dday, dmonth, dyear, dyearold, dyear0
      character*80 attstr
      character*80 fname
      character, allocatable :: str(:)
      integer,parameter :: IDAT=192
      integer,parameter :: JDAT=94
      real*8,parameter :: XSTEPDAT=1.875
      real*8,parameter :: YSTEPDAT=1.888888888888889
      real*8,parameter :: pi=3.14159265/180.0
      real*8 ::  Xdat(IDAT), Ydat(JDAT)
      real*8 ::  Xdat0(IDAT), Ydat0(JDAT)
      Real*4,dimension(IDAT,JDAT) ::  UB, VB, WORK
      Real*4 uoffset,uscale,voffset,vscale
      character*80 input_path
      common/input_data/input_path
      character*4 cyear

      save

      logical firstcall, newyear, flat, first
      integer nu_diag
      common/iodiag/nu_diag
      save firstcall
      data firstcall/.true./,nold/-1/,dyearold/-1/,dindex
     */0,31,59,90,120,151,181,212,243,273,304,334/

      flat=.false.
      dyear0=int(idate/10000)
c  Data are started from 1958, so all previous years are considered as 1958
      dyear=max(dyear0,1958)
      newyear=dyearold.ne.dyear
      
      dmonth=int((idate-10000*dyear0)/100)
      dday=int(idate-10000*dyear0-100*dmonth)
      
      n=(dindex(dmonth)+dday-1)*4+int((1.0*sec+3*3600)/(6.0*3600.0))
c      if(4*int(dyear/4).eq.4.and.dmonth.gt.2)n=n+1
      
c      Write(nu_diag,*)'In CIafWind'
c      Write(nu_diag,*)'time: idate=',idate,' sec=',sec
c      Write(nu_diag,*)'record no.',n
      n=max(n,1)
      n=min(n,1460)
c      Write(nu_diag,*)'record no.',n
      
	if(firstcall)then
	
      allocate(rxbox(imt,jmt),rybox(imt,jmt),i0box(imt,jmt),
     *   j0box(imt,jmt),i1box(imt,jmt),j1box(imt,jmt),nodetype(imt,jmt),
     *   UP(imt,jmt),VP(imt,jmt))

      if(vers.eq.1)then
         if(dyear.lt.1958)dyear=1958
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_1.0/u_10.'//
     *      cyear//'.2004_08_03.nc'
         if(dyear.gt.2000)
     *      fname=trim(input_path)//'CIAF_1.0/u_10.'//
     *      cyear//'.2006_04_06.nc'
      elseif(vers.eq.2)then
         if(dyear.lt.1948)dyear=1948
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_2.0/u_10.'//
     *      cyear//'.05APR2010.nc'
      elseif(vers.eq.0)then
         fname=trim(input_path)//'CIAF_2.0/u_10.'//
     *      '15JUNE2009.nc'
      endif
         ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'U-Wind Error: reading NetCDF file'
            call exit_coupler
         endif

      beg(1)=1
      beg(2)=1
      beg(3)=1
      cnt(1)=IDAT
      cnt(2)=1
      cnt(3)=1
      
      if(vers.eq.1)ierr=nf_inq_varid(ncid,'LON1',vid)
      if(vers.eq.2.or.vers.eq.0)ierr=nf_inq_varid(ncid,'LON',vid)
      ierr=nf_get_vara_double(ncid,vid,beg,cnt,Xdat0(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'U-Wind Error: Error reading lon'
         call exit_coupler
      endif

      cnt(1)=JDAT

      if(vers.eq.1)ierr=nf_inq_varid(ncid,'LAT1',vid)
      if(vers.eq.2.or.vers.eq.0)ierr=nf_inq_varid(ncid,'LAT',vid)
      ierr=nf_get_vara_double(ncid,vid,beg,cnt,Ydat0(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'U-Wind Error: Error reading lat'
         call exit_coupler
      endif
      ierr=nf_close(ncid)

	   do i=1,IDAT
	      Xdat(i)=Xdat0(i)
	   enddo
	   do j=1,JDAT
c	      Ydat(j)=Ydat0(JDAT-j+1)
	      Ydat(j)=Ydat0(j)
           enddo
	   do i=1,imt
	      do j=1,jmt
	         xx=X(i,j)
		 if(xx.lt.0)xx=xx+360.0
		 yy=Y(i,j)
		 mm=mask(i,j)
		    i0=0
		    do ii=1,IDAT-1
		       if(Xdat(ii).le.xx.and.Xdat(ii+1).ge.xx)then
		          i0=ii
			  i1=ii+1
			  rxbox(i,j)=(xx-Xdat(i0))/(Xdat(i1)-Xdat(i0))
			  exit
		       endif
		    enddo
		    if(i0.eq.0)then
		       if(Xdat(IDAT).le.xx.and.(Xdat(1)+360).ge.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0))/
     *  		    (Xdat(i1)+360-Xdat(i0))
                   elseif(Xdat(1).ge.xx.and.(Xdat(IDAT)-360).le.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0)+360)/
     *  		    (Xdat(i1)+360-Xdat(i0))
		       else
                          Write(nu_diag,*) 'Ciafwind: Grid from NetCDF is not covering'
			  Write(nu_diag,*) xx,Xdat(1),Xdat(IDAT)
                          call exit_coupler
		       endif
		    endif
		    j0=0
		    if(yy.le.Ydat(2))then
		       j0=2
		       j1=2
		       nodetype(i,j)=1
                       rybox(i,j)=(yy+90.0)/(Ydat(j1)+90.0)
		    else if(yy.ge.Ydat(JDAT-1))then
		       j0=JDAT-1
		       j1=JDAT-1
		       nodetype(i,j)=1
                       rybox(i,j)=(yy-Ydat(j0))/(90.0-Ydat(j0))
		    else
		    do jj=1,JDAT-1
		       if(Ydat(jj).le.yy.and.Ydat(jj+1).ge.yy)then
		          j0=jj
			  j1=jj+1
			  rybox(i,j)=(yy-Ydat(j0))/(Ydat(j1)-Ydat(j0))
			  exit
		       endif
		    enddo
		    nodetype(i,j)=0
		    endif
		    i0box(i,j)=i0
		    i1box(i,j)=i1
		    j0box(i,j)=j0
		    j1box(i,j)=j1
	      enddo
	   enddo
	   
        endif
	
      if(n.ne.nold)then
      
      if(newyear)then

         if(dyearold.ge.0)then      
	   ierr=nf_close(ncidu)
	 endif
      if(vers.eq.1)then
         if(dyear.lt.1958)dyear=1958
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_1.0/u_10.'//
     *      cyear//'.2004_08_03.nc'
         if(dyear.gt.2000)
     *      fname=trim(input_path)//'CIAF_1.0/u_10.'//
     *      cyear//'.2006_04_06.nc'
      elseif(vers.eq.2)then
         if(dyear.lt.1948)dyear=1948
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_2.0/u_10.'//
     *      cyear//'.05APR2010.nc'
      elseif(vers.eq.0)then
         fname=trim(input_path)//'CIAF_2.0/u_10.'//
     *      '15JUNE2009.nc'
      endif
c	 Write(nu_diag,*)'Reading file: ',fname
         ierr=nf_open(trim(fname),NF_NOWRITE,ncidu)
	   if(ierr.NE.NF_NOERR)then
              Write(nu_diag,*) 'U-Wind Error: reading NetCDF file'
              call exit_coupler
c           else
c              Write(nu_diag,'(a,i3)')'Open CDF ', ncidu
           endif
      endif

c      Write(nu_diag,'(a,i3)')'Read record no.',n,' from CDF ', ncidu

      beg(1)=1
      beg(2)=1
      beg(3)=n
      cnt(1)=IDAT
      cnt(2)=JDAT
      cnt(3)=1
      
      ierr=nf_inq_varid(ncidu,'U_10_MOD',vid)
      uoffset=0.0
      uscale=1.0
      ierr=nf_get_vara_real(ncidu,vid,beg,cnt,WORK(1,1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'U-Wind Error: Error reading uwnd'
         call exit_coupler
      endif
      Do i=1,IDAT
         Do j=1,JDAT
            UB(i,j)=WORK(i,j)
         Enddo
      Enddo
c      UB=uscale*UB+uoffset

      if(newyear)then

         if(dyearold.ge.0)then      
	   ierr=nf_close(ncidv)
	 endif
      if(vers.eq.1)then
         if(dyear.lt.1958)dyear=1958
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_1.0/v_10.'//
     *      cyear//'.2004_08_03.nc'
         if(dyear.gt.2000)
     *      fname=trim(input_path)//'CIAF_1.0/v_10.'//
     *      cyear//'.2006_04_06.nc'
      elseif(vers.eq.2)then
         if(dyear.lt.1948)dyear=1948
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_2.0/v_10.'//
     *      cyear//'.05APR2010.nc'
      elseif(vers.eq.0)then
         fname=trim(input_path)//'CIAF_2.0/v_10.'//
     *      '15JUNE2009.nc'
      endif
	   ierr=nf_open(trim(fname),NF_NOWRITE,ncidv)
	   if(ierr.NE.NF_NOERR)then
              Write(nu_diag,*) 'V-Wind Error: reading NetCDF file'
              call exit_coupler
c           else
c              Write(nu_diag,'(a,i3)')'Open CDF ', ncidv
           endif
	   dyearold=dyear
      endif

c      Write(nu_diag,'(a,i3)')'Read record no.',n,' from CDF ', ncidv

      beg(1)=1
      beg(2)=1
      beg(3)=n
      cnt(1)=IDAT
      cnt(2)=JDAT
      cnt(3)=1
      
      ierr=nf_inq_varid(ncidv,'V_10_MOD',vid)
      voffset=0.0
      vscale=1.0
      ierr=nf_get_vara_real(ncidv,vid,beg,cnt,WORK(1,1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'V-Wind Error: Error reading vwnd'
         call exit_coupler
      endif
      Do i=1,IDAT
         Do j=1,JDAT
            VB(i,j)=WORK(i,j)
         Enddo
      Enddo
c      VB=vscale*VB+voffset

c      Write(nu_diag,'(a,f6.2,a,f6.2,a)')'Read Control: U(115,31)= ',
c     *   UB(115,31)*uscale+uoffset,' m/s,  V(115,31)= ',
c     *   VB(115,31)*vscale+voffset,' m/s'

      PMIN=1e+10
      PMAX=-1e+10
      pp=atan2(1.0,1.0)*4.0/180.0
c      if(i.eq.80.and.j.eq.148)Write(nu_diag,*)0,pp

      first=.true.
      do j=1,jmt
         do i=1,imt
            if(nodetype(i,j).eq.0)then
               i0=i0box(i,j)
               i1=i1box(i,j)
               j0=j0box(i,j)
               j1=j1box(i,j)
               rx0=rxbox(i,j)
               rx1=1.0-rx0
               ry0=rybox(i,j)
               ry1=1.0-ry0
               uu= 
     *            ry1*(rx1*dble(UB(i0,j0))+rx0*dble(UB(i1,j0)))+
     *            ry0*(rx1*dble(UB(i0,j1))+rx0*dble(UB(i1,j1)))
               uu=uu*uscale+uoffset

               vv=
     *            ry1*(rx1*dble(VB(i0,j0))+rx0*dble(VB(i1,j0)))+
     *            ry0*(rx1*dble(VB(i0,j1))+rx0*dble(VB(i1,j1)))
               vv=vv*vscale+voffset

               if(Y(i,j).le.60.0)then
                  UP(i,j)=uu
                  VP(i,j)=vv
	       else
	          sn=sin(X(i,j)*pp)
		  cs=cos(X(i,j)*pp)
                  if(i.ne.1.and.i.ne.imt)then
                     xx1=(90.0-Y(i+1,j))*sin(X(i+1,j)*pp)
                     yy1=-(90.0-Y(i+1,j))*cos(X(i+1,j)*pp)
                     xx2=(90.0-Y(i-1,j))*sin(X(i-1,j)*pp)
                     yy2=-(90.0-Y(i-1,j))*cos(X(i-1,j)*pp)
                  else
                     if(i.eq.1)then
                        xx1=(90.0-Y(i+1,j))*sin(X(i+1,j)*pp)
                        yy1=-(90.0-Y(i+1,j))*cos(X(i+1,j)*pp)
                        xx2=(90.0-Y(i,j))*sn
                        yy2=-(90.0-Y(i,j))*cs
                     else
                        xx1=(90.0-Y(i,j))*sn
                        yy1=-(90.0-Y(i,j))*cs
                        xx2=(90.0-Y(i-1,j))*sin(X(i-1,j)*pp)
                        yy2=-(90.0-Y(i-1,j))*cos(X(i-1,j)*pp)
                     endif
                  endif
                  dx=xx1-xx2
                  dy=yy1-yy2
                  alph=atan2(dy,dx)
		  uu0=uu*cs-vv*sn
		  vv0=uu*sn+vv*cs
                  UP(i,j)=uu0*cos(alph)+vv0*sin(alph)
                  VP(i,j)=-uu0*sin(alph)+vv0*cos(alph)
	       endif
              else
               UP(i,j)=0.0
               VP(i,j)=0.0
               sn=sin(X(i,j))
               cs=cos(X(i,j))
               if(i.ne.1.and.i.ne.imt)then
                  xx1=(90.0-Y(i+1,j))*sin(X(i+1,j)*pp)
                  yy1=-(90.0-Y(i+1,j))*cos(X(i+1,j)*pp)
                  xx2=(90.0-Y(i-1,j))*sin(X(i-1,j)*pp)
                  yy2=-(90.0-Y(i-1,j))*cos(X(i-1,j)*pp)
               else
                  if(i.eq.1)then
                     xx1=(90.0-Y(i+1,j))*sin(X(i+1,j)*pp)
                     yy1=-(90.0-Y(i+1,j))*cos(X(i+1,j)*pp)
                     xx2=(90.0-Y(i,j))*sn
                     yy2=-(90.0-Y(i,j))*cs
                  else
                     xx1=(90.0-Y(i,j))*sn
                     yy1=-(90.0-Y(i,j))*cs
                     xx2=(90.0-Y(i-1,j))*sin(X(i-1,j)*pp)
                     yy2=-(90.0-Y(i-1,j))*cos(X(i-1,j)*pp)
                  endif
               endif
               dx=xx1-xx2
               dy=yy1-yy2
               alph=atan2(dy,dx)
	       uu=0.0
	       vv=0.0
               j0=j0box(i,j)
               j1=j1box(i,j)
               rr=1.0/float(IDAT)
               if(flat)then
                  do ii=1,IDAT
                     uur=dble(UB(ii,j0))*uscale+uoffset
                     vvr=dble(VB(ii,j0))*vscale+voffset
                     uu=uu+
     *rr*(uur*cos(pp*Xdat(ii))-vvr*sin(pp*Xdat(ii)))
                     vv=vv+
     *rr*(uur*sin(pp*Xdat(ii))+vvr*cos(pp*Xdat(ii)))
                  enddo
               else

c        Components at the pole as mean values
                  if(first)then
                  first=.false.
                  UPOLE=0.0
                  VPOLE=0.0
                  do ii=1,IDAT
                     uur=dble(UB(ii,j0))*uscale+uoffset
                     vvr=dble(VB(ii,j0))*vscale+voffset
                     UPOLE=UPOLE+
     *rr*(uur*cos(pp*Xdat(ii))-vvr*sin(pp*Xdat(ii)))
                     VPOLE=VPOLE+
     *rr*(uur*sin(pp*Xdat(ii))+vvr*cos(pp*Xdat(ii)))
c                     Write(nu_diag,*) 'u=',
c     *                   uur*cos(pp*Xdat(ii))-vvr*sin(pp*Xdat(ii)),
c     *                   ',   v=',
c     *                   uur*sin(pp*Xdat(ii))+vvr*cos(pp*Xdat(ii))
                  enddo
c                  Write(nu_diag,*) 'UPOLE=',upole,',   VPOLE=',vpole
                  endif

                  i0=i0box(i,j)
                  i1=i1box(i,j)

c                  Write(nu_diag,*) i,j
c                  Write(nu_diag,*) i0,j0,i1,j1

                  rx0=rxbox(i,j)
                  rx1=1.0-rx0
                  ry0=rybox(i,j)
                  ry1=1.0-ry0

c                  Write(nu_diag,*) 'rx0=',rx0,
c     *                         ',   rx1=',rx1

c        Rotate data velocity vector at i0,j0 to make it cartesian
                  uu0=dble(UB(i0,j0))*uscale+uoffset
                  vv0=dble(VB(i0,j0))*vscale+voffset
                  uu=uu0*cos(pp*Xdat(i0))-vv0*sin(pp*Xdat(i0))
                  vv=uu0*sin(pp*Xdat(i0))+vv0*cos(pp*Xdat(i0))
                  uu0=uu
                  vv0=vv
c                  Write(nu_diag,*) 'u0=',uu0,',   v0=',vv0

c        Rotate data velocity vector at i1,j0 to make it cartesian
                  uu1=dble(UB(i1,j0))*uscale+uoffset
                  vv1=dble(VB(i1,j0))*vscale+voffset
                  uu=uu1*cos(pp*Xdat(i1))-vv1*sin(pp*Xdat(i1))
                  vv=uu1*sin(pp*Xdat(i1))+vv1*cos(pp*Xdat(i1))
                  uu1=uu
                  vv1=vv
c                  Write(nu_diag,*) 'u1=',uu1,',   v1=',vv1

c        Interpolate along longitude between i0,j0 and i1,j0
                  uu=rx1*uu0+rx0*uu1
                  vv=rx1*vv0+rx0*vv1

c        Apply offset and scaling factor
c                  uu=uu*uscale+uoffset
c                  vv=vv*vscale+voffset

c                  Write(nu_diag,*) 'u=',uu,',   v=',vv

c        Interpolate along latitude between j0 and pole,
c        considering pole value to be (UPOLE,VPOLE)
                  uu=ry1*uu+ry0*UPOLE
                  vv=ry1*vv+ry0*VPOLE

c                  Write(nu_diag,*) 'ry0=',ry0,
c     *                         ',   ry1=',ry1
c                  Write(nu_diag,*) 'U=',uu,',   V=',vv

               endif
               UP(i,j)=uu*cos(alph)+vv*sin(alph)
               VP(i,j)=-uu*sin(alph)+vv*cos(alph)
            endif
            PPU=sqrt(UP(i,j)**2+VP(i,j)**2)
	    if(PMAX.lt.PPU)PMAX=PPU
	    if(PMIN.gt.PPU)PMIN=PPU
         enddo
      enddo

c      Write(nu_diag,*)'Max wind: ', PMAX
c      Write(nu_diag,*)'Min wind: ', PMIN
      
      endif
      
      nold=n
      
      do j=1,jmt
         do i=1,imt
	    UBOT(i,j)=UP(i,j)
	    VBOT(i,j)=VP(i,j)
         enddo
      enddo

      
      firstcall=.false.

      return
      end

      subroutine ciaf_airtemp(idate, sec, imt, jmt, X, Y, mask, PBOT,
     *                        vers)
      
      include 'netcdf.inc'
      
      integer*4, intent(in) :: idate, imt, jmt, mask(imt,jmt), sec, vers
      real*8, intent(in) :: X(imt,jmt), Y(imt,jmt)
      real*8, intent(out) :: PBOT(imt,jmt)
      real*8, allocatable, dimension(:,:) ::
     &  rxbox,rybox,PP
      integer, allocatable, dimension(:,:) ::
     &  i0box,j0box, i1box,j1box, nodetype
      integer dindex(12)
      integer i,j,n,ncid,vid,did,beg(3),cnt(3),ierr,nold
      integer dday, dmonth, dyear, dyearold, dyear0
      character*80 attstr
      character*80 fname
      character, allocatable :: str(:)
      integer,parameter :: IDAT=192
      integer,parameter :: JDAT=94
      real*8,parameter :: XSTEPDAT=360.0/IDAT
      real*8,parameter :: YSTEPDAT=2*(90.0/(JDAT+1))
      real*8,parameter :: pi=3.14159265/180.0
      real*8 ::  Xdat(IDAT), Ydat(JDAT)
      Real*4,dimension(IDAT,JDAT) ::  PB
      Real*4 offset,scale
      character*80 input_path
      common/input_data/input_path
      character*4 cyear

      save

      logical firstcall, newyear, flat
      integer nu_diag
      common/iodiag/nu_diag
      save firstcall
      data firstcall/.true./,nold/-1/,dyearold/-1/,dindex
     */0,31,59,90,120,151,181,212,243,273,304,334/

      flat=.false.
      dyear0=int(idate/10000)
c  Data are started from 1958, so all previous years are considered as 1958
      dyear=max(dyear0,1958)
      newyear=dyearold.ne.dyear
      
      dmonth=int((idate-10000*dyear0)/100)
      dday=int(idate-10000*dyear0-100*dmonth)
      
      n=(dindex(dmonth)+dday-1)*4+int((1.0*sec+3*3600)/(6.0*3600.0))
      n=max(n,1)
      n=min(n,1460)
c      if(4*int(dyear/4).eq.4.and.dmonth.gt.2)n=n+1
      
c      Write(nu_diag,*)'In CIAF AirTemp'
      
	if(firstcall)then
	
      allocate(rxbox(imt,jmt),rybox(imt,jmt),i0box(imt,jmt),
     *   j0box(imt,jmt),i1box(imt,jmt),j1box(imt,jmt),nodetype(imt,jmt),
     *   PP(imt,jmt))

      if(vers.eq.1)then
         if(dyear.lt.1958)dyear=1958
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_1.0/t_10.'//
     *      cyear//'.2004_08_03.nc'
         if(dyear.gt.2000)
     *      fname=trim(input_path)//'CIAF_1.0/t_10.'//
     *      cyear//'.2006_04_06.nc'
      elseif(vers.eq.2)then
         if(dyear.lt.1948)dyear=1948
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_2.0/t_10.'//
     *      cyear//'.05APR2010.nc'
      elseif(vers.eq.0)then
         fname=trim(input_path)//'CIAF_2.0/t_10.'//
     *      '15JUNE2009.nc'
      endif
         ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'CIAF AirTemp Error: reading NetCDF file'
            call exit_coupler
         endif

      beg(1)=1
      beg(2)=1
      beg(3)=1
      cnt(1)=IDAT
      cnt(2)=1
      cnt(3)=1
      
      ierr=nf_inq_varid(ncid,'LON',vid)
      ierr=nf_get_vara_double(ncid,vid,beg,cnt,Xdat(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'AirTemp Error: Error reading lon'
         call exit_coupler
      endif

      cnt(1)=JDAT

      ierr=nf_inq_varid(ncid,'LAT',vid)
      ierr=nf_get_vara_double(ncid,vid,beg,cnt,Ydat(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'AirTemp Error: Error reading lat'
         call exit_coupler
      endif
      ierr=nf_close(ncid)

	   do i=1,imt
	      do j=1,jmt
	         xx=X(i,j)
		 if(xx.lt.0)xx=xx+360.0
		 yy=Y(i,j)
		 mm=mask(i,j)
		    i0=0
		    do ii=1,IDAT-1
		       if(Xdat(ii).le.xx.and.Xdat(ii+1).ge.xx)then
		          i0=ii
			  i1=ii+1
			  rxbox(i,j)=(xx-Xdat(i0))/(Xdat(i1)-Xdat(i0))
			  exit
		       endif
		    enddo
		    if(i0.eq.0)then
		       if(Xdat(IDAT).le.xx.and.(Xdat(1)+360).ge.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0))/
     *  		    (Xdat(i1)+360-Xdat(i0))
                   elseif(Xdat(1).ge.xx.and.(Xdat(IDAT)-360).le.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0)+360)/
     *  		    (Xdat(i1)+360-Xdat(i0))
		       else
                          Write(nu_diag,*) 'Ciaf_airtemp: Grid from NetCDF is not covering'
			  Write(nu_diag,*) xx,Xdat(1),Xdat(IDAT)
                          call exit_coupler
		       endif
		    endif
		    j0=0
		    if(yy.le.Ydat(2))then
		       j0=2
		       j1=2
		       nodetype(i,j)=1
                       rybox(i,j)=(yy+90.0)/(Ydat(j1)+90.0)
		    else if(yy.ge.Ydat(JDAT-1))then
		       j0=JDAT-1
		       j1=JDAT-1
		       nodetype(i,j)=1
                       rybox(i,j)=(yy-Ydat(j0))/(90.0-Ydat(j0))
		    else
		    do jj=1,JDAT-1
		       if(Ydat(jj).le.yy.and.Ydat(jj+1).ge.yy)then
		          j0=jj
			  j1=jj+1
			  rybox(i,j)=(yy-Ydat(j0))/(Ydat(j1)-Ydat(j0))
			  exit
		       endif
		    enddo
		    nodetype(i,j)=0
		    endif
		    i0box(i,j)=i0
		    i1box(i,j)=i1
		    j0box(i,j)=j0
		    j1box(i,j)=j1
	      enddo
	   enddo
	   
        endif
	
      if(n.ne.nold)then
      
      if(newyear)then

         if(dyearold.ge.0)then      
	   ierr=nf_close(ncid)
	 endif
      if(vers.eq.1)then
         if(dyear.lt.1958)dyear=1958
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_1.0/t_10.'//
     *      cyear//'.2004_08_03.nc'
         if(dyear.gt.2000)
     *      fname=trim(input_path)//'CIAF_1.0/t_10.'//
     *      cyear//'.2006_04_06.nc'
      elseif(vers.eq.2)then
         if(dyear.lt.1948)dyear=1948
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_2.0/t_10.'//
     *      cyear//'.05APR2010.nc'
      elseif(vers.eq.0)then
         fname=trim(input_path)//'CIAF_2.0/t_10.'//
     *      '15JUNE2009.nc'
      endif
	   ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
	   if(ierr.NE.NF_NOERR)then
              Write(nu_diag,*) 'AirTemp Error: reading NetCDF file'
              call exit_coupler
c           else
c              Write(nu_diag,'(a,i3)')'Open CDF ', ncid
           endif
	   dyearold=dyear
      endif

c      Write(nu_diag,'(a,i3)')'Read record no.',n,' from CDF ', ncid

      beg(1)=1
      beg(2)=1
      beg(3)=n
      cnt(1)=IDAT
      cnt(2)=JDAT
      cnt(3)=1
      
      ierr=nf_inq_varid(ncid,'T_10_MOD',vid)
      offset=0.0
      scale=1.0
      ierr=nf_get_vara_real(ncid,vid,beg,cnt,PB(1,1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'AirTemp Error: Error reading air'
         call exit_coupler
      endif
      
c      Write(nu_diag,'(a,f6.2,a)')'Read Control: TA(115,31)= ',
c     *   PB(115,31)*scale+offset,' K'

      PMIN=1e+10
      PMAX=-1e+10
      do j=1,jmt
         do i=1,imt
             i0=i0box(i,j)
             i1=i1box(i,j)
             j0=j0box(i,j)
             j1=j1box(i,j)
             rx0=rxbox(i,j)
             rx1=1.0-rx0
             ry0=rybox(i,j)
             ry1=1.0-ry0
	       if(nodetype(i,j).eq.0)then
               PP(i,j)= scale*(
     *            ry1*(rx1*dble(PB(i0,j0))+rx0*dble(PB(i1,j0)))+
     *            ry0*(rx1*dble(PB(i0,j1))+rx0*dble(PB(i1,j1)))) +
     *            offset
              else
               PP(i,j)=0.0
	       rr=1.0/IDAT
               if(flat)then
	       do ii=1,IDAT
                  PP(i,j)=PP(i,j)+rr*dble(PB(ii,j0))
	       enddo
	       PP(i,j)=scale*PP(i,j)+offset
               else
c        Components at the pole as mean values
                  PPOLE=0.0
                  do ii=1,IDAT
                     uur=dble(PB(ii,j0))*scale+offset
                     PPOLE=PPOLE+rr*uur
                  enddo

                  uu0=dble(PB(i0,j0))*scale+offset
                  uu1=dble(PB(i1,j0))*scale+offset

                  uu=rx1*uu0+rx0*uu1

                  PP(i,j)=ry1*uu+ry0*PPOLE
               endif
            endif
	    if(PMAX.lt.PP(i,j))PMAX=PP(i,j)
	    if(PMIN.gt.PP(i,j))PMIN=PP(i,j)
         enddo
      enddo

c      Write(nu_diag,*)'Max temperature: ', PMAX
c      Write(nu_diag,*)'Min temperature: ', PMIN
c      Write(nu_diag,*)'Pole temperature: ', PP(88,136)-273.16
      
      else

c      Write(nu_diag,*)'Old values not updated.'
      
      endif
      
      nold=n
      
      do j=1,jmt
         do i=1,imt
	    PBOT(i,j)=PP(i,j)
         enddo
      enddo

      
      firstcall=.false.

      return
      end

      subroutine ciaf_airhum(idate, sec, imt, jmt, X, Y, mask, PBOT, 
     *                       vers)
      
      include 'netcdf.inc'
      
      integer*4, intent(in) :: idate, imt, jmt, mask(imt,jmt), sec, vers
      real*8, intent(in) :: X(imt,jmt), Y(imt,jmt)
      real*8, intent(out) :: PBOT(imt,jmt)
      real*8, allocatable, dimension(:,:) ::
     &  rxbox,rybox,PP
      integer, allocatable, dimension(:,:) ::
     &  i0box,j0box, i1box,j1box, nodetype
      integer dindex(12)
      integer i,j,n,ncid,vid,did,beg(3),cnt(3),ierr,nold
      integer dday, dmonth, dyear, dyearold, dyear0
      character*80 attstr
      character*80 fname
      character, allocatable :: str(:)
      integer,parameter :: IDAT=192
      integer,parameter :: JDAT=94
      real*8,parameter :: XSTEPDAT=360.0/IDAT
      real*8,parameter :: YSTEPDAT=2*(90.0/(JDAT+1))
      real*8,parameter :: pi=3.14159265/180.0
      real*8 ::  Xdat(IDAT), Ydat(JDAT)
      Real*4,dimension(IDAT,JDAT) ::  PB
      Real*4 offset,scale
      character*80 input_path
      common/input_data/input_path
      character*4 cyear

      save

      logical firstcall, newyear, flat
      integer nu_diag
      common/iodiag/nu_diag
      save firstcall
      data firstcall/.true./,nold/-1/,dyearold/-1/,dindex
     */0,31,59,90,120,151,181,212,243,273,304,334/

      flat=.false.
      dyear0=int(idate/10000)
c  Data are started from 1958, so all previous years are considered as 1958
      dyear=max(dyear0,1958)
      newyear=dyearold.ne.dyear
      
      dmonth=int((idate-10000*dyear0)/100)
      dday=int(idate-10000*dyear0-100*dmonth)
      
      n=(dindex(dmonth)+dday-1)*4+int((1.0*sec+3*3600)/(6.0*3600.0))
      n=max(n,1)
      n=min(n,1460)
c      if(4*int(dyear/4).eq.4.and.dmonth.gt.2)n=n+1
      
c      Write(nu_diag,*)'In CIAF AirHum'
      
	if(firstcall)then
	
      allocate(rxbox(imt,jmt),rybox(imt,jmt),i0box(imt,jmt),
     *   j0box(imt,jmt),i1box(imt,jmt),j1box(imt,jmt),nodetype(imt,jmt),
     *   PP(imt,jmt))

      if(vers.eq.1)then
c         fname='../../MassStore/sst/air.2m.gauss.YYYY.nc'
         if(dyear.lt.1958)dyear=1958
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_1.0/q_10.'//
     *      cyear//'.2004_08_03.nc'
         if(dyear.gt.2000)
     *      fname=trim(input_path)//'CIAF_1.0/q_10.'//
     *      cyear//'.2006_04_06.nc'
      elseif(vers.eq.2)then
         if(dyear.lt.1948)dyear=1948
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_2.0/q_10.'//
     *      cyear//'.05APR2010.nc'
      elseif(vers.eq.0)then
         fname=trim(input_path)//'CIAF_2.0/q_10.'//
     *      '15JUNE2009.nc'
      endif
         ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'CIAF AirHum Error: reading NetCDF file'
            call exit_coupler
         endif

      beg(1)=1
      beg(2)=1
      beg(3)=1
      cnt(1)=IDAT
      cnt(2)=1
      cnt(3)=1
      
      ierr=nf_inq_varid(ncid,'LON',vid)
      ierr=nf_get_vara_double(ncid,vid,beg,cnt,Xdat(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'AirHum Error: Error reading lon'
         call exit_coupler
      endif

      cnt(1)=JDAT

      ierr=nf_inq_varid(ncid,'LAT',vid)
      ierr=nf_get_vara_double(ncid,vid,beg,cnt,Ydat(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'AirHum Error: Error reading lat'
         call exit_coupler
      endif
      ierr=nf_close(ncid)

	   do i=1,imt
	      do j=1,jmt
	         xx=X(i,j)
		 if(xx.lt.0)xx=xx+360.0
		 yy=Y(i,j)
		 mm=mask(i,j)
		    i0=0
		    do ii=1,IDAT-1
		       if(Xdat(ii).le.xx.and.Xdat(ii+1).ge.xx)then
		          i0=ii
			  i1=ii+1
			  rxbox(i,j)=(xx-Xdat(i0))/(Xdat(i1)-Xdat(i0))
			  exit
		       endif
		    enddo
		    if(i0.eq.0)then
		       if(Xdat(IDAT).le.xx.and.(Xdat(1)+360).ge.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0))/
     *  		    (Xdat(i1)+360-Xdat(i0))
                   elseif(Xdat(1).ge.xx.and.(Xdat(IDAT)-360).le.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0)+360)/
     *  		    (Xdat(i1)+360-Xdat(i0))
		       else
                          Write(nu_diag,*) 'Ciaf_airhum: Grid from NetCDF is not covering'
			  Write(nu_diag,*) xx,Xdat(1),Xdat(IDAT)
                          call exit_coupler
		       endif
		    endif
		    j0=0
		    if(yy.le.Ydat(2))then
		       j0=2
		       j1=2
		       nodetype(i,j)=1
                       rybox(i,j)=(yy+90.0)/(Ydat(j1)+90.0)
		    else if(yy.ge.Ydat(JDAT-1))then
		       j0=JDAT-1
		       j1=JDAT-1
		       nodetype(i,j)=1
                       rybox(i,j)=(yy-Ydat(j0))/(90.0-Ydat(j0))
		    else
		    do jj=1,JDAT-1
		       if(Ydat(jj).le.yy.and.Ydat(jj+1).ge.yy)then
		          j0=jj
			  j1=jj+1
			  rybox(i,j)=(yy-Ydat(j0))/(Ydat(j1)-Ydat(j0))
			  exit
		       endif
		    enddo
		    nodetype(i,j)=0
		    endif
		    i0box(i,j)=i0
		    i1box(i,j)=i1
		    j0box(i,j)=j0
		    j1box(i,j)=j1
	      enddo
	   enddo
	   
        endif
	
      if(n.ne.nold)then
      
      if(newyear)then

         if(dyearold.ge.0)then      
	   ierr=nf_close(ncid)
	 endif
      if(vers.eq.1)then
         if(dyear.lt.1958)dyear=1958
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_1.0/q_10.'//
     *      cyear//'.2004_08_03.nc'
         if(dyear.gt.2000)
     *      fname=trim(input_path)//'CIAF_1.0/q_10.'//
     *      cyear//'.2006_04_06.nc'
      elseif(vers.eq.2)then
         if(dyear.lt.1948)dyear=1948
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_2.0/q_10.'//
     *      cyear//'.05APR2010.nc'
      elseif(vers.eq.0)then
         fname=trim(input_path)//'CIAF_2.0/q_10.'//
     *      '15JUNE2009.nc'
      endif
c	 Write(nu_diag,*)'Reading file: ',fname
	   ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
	   if(ierr.NE.NF_NOERR)then
              Write(nu_diag,*) 'AirHum Error: reading NetCDF file'
              call exit_coupler
c           else
c              Write(nu_diag,'(a,i3)')'Open CDF ', ncid
           endif
	   dyearold=dyear
      endif

c      Write(nu_diag,'(a,i4)')'Read record no.',n,' from CDF ', ncid

      beg(1)=1
      beg(2)=1
      beg(3)=n
      cnt(1)=IDAT
      cnt(2)=JDAT
      cnt(3)=1
      
      ierr=nf_inq_varid(ncid,'Q_10_MOD',vid)
      offset=0.0
      scale=1.0
      ierr=nf_get_vara_real(ncid,vid,beg,cnt,PB(1,1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'AirHum Error: Error reading air'
         call exit_coupler
      endif
      
c      Write(nu_diag,'(a,f6.2,a)')'Read Control: QA(115,31)= ',
c     *   PB(115,31)*scale+offset,' K'

      PMIN=1e+10
      PMAX=-1e+10
      do j=1,jmt
         do i=1,imt
             i0=i0box(i,j)
             i1=i1box(i,j)
             j0=j0box(i,j)
             j1=j1box(i,j)
             rx0=rxbox(i,j)
             rx1=1.0-rx0
             ry0=rybox(i,j)
             ry1=1.0-ry0
	       if(nodetype(i,j).eq.0)then
               PP(i,j)= scale*(
     *            ry1*(rx1*dble(PB(i0,j0))+rx0*dble(PB(i1,j0)))+
     *            ry0*(rx1*dble(PB(i0,j1))+rx0*dble(PB(i1,j1)))) +
     *            offset
              else
               PP(i,j)=0.0
	       rr=1.0/IDAT
               if(flat)then
	       do ii=1,IDAT
                  PP(i,j)=PP(i,j)+rr*dble(PB(ii,j0))
	       enddo
	       PP(i,j)=scale*PP(i,j)+offset
               else
c        Components at the pole as mean values
                  PPOLE=0.0
                  do ii=1,IDAT
                     uur=dble(PB(ii,j0))*scale+offset
                     PPOLE=PPOLE+rr*uur
                  enddo

                  uu0=dble(PB(i0,j0))*scale+offset
                  uu1=dble(PB(i1,j0))*scale+offset

                  uu=rx1*uu0+rx0*uu1

                  PP(i,j)=ry1*uu+ry0*PPOLE
               endif
            endif
	    if(PMAX.lt.PP(i,j))PMAX=PP(i,j)
	    if(PMIN.gt.PP(i,j))PMIN=PP(i,j)
         enddo
      enddo

c      Write(nu_diag,*)'Max humidity: ', PMAX
c      Write(nu_diag,*)'Min humidity: ', PMIN
c      Write(nu_diag,*)'Pole humidity: ', PP(88,136)
      
      else

c      Write(nu_diag,*)'Old values not updated.'
      
      endif
      
      nold=n
      
      do j=1,jmt
         do i=1,imt
	    PBOT(i,j)=PP(i,j)
         enddo
      enddo

      
      firstcall=.false.

      return
      end

      subroutine rhumid(idate, imt, jmt, X, Y, mask, PBOT)
      
      include 'netcdf.inc'
      
      integer*4, intent(in) :: idate, imt, jmt, mask(imt,jmt)
      real*8, intent(in) :: X(imt,jmt), Y(imt,jmt)
      real*8, intent(out) :: PBOT(imt,jmt)
      real*8, allocatable, dimension(:,:) ::
     &  rxbox,rybox,PP
      integer, allocatable, dimension(:,:) ::
     &  i0box,j0box, i1box,j1box, nodetype
      integer dindex(12)
      integer i,j,n,ncid,vid,did,beg(3),cnt(3),ierr,nold
      integer dday, dmonth, dyear, dyearold
      character*80 attstr
      character*80 fname
      character, allocatable :: str(:)
      integer,parameter :: IDAT=144
      integer,parameter :: JDAT=73
      real*8,parameter :: XSTEPDAT=2.5
      real*8,parameter :: YSTEPDAT=2.5
      real*8,parameter :: pi=3.14159265/180.0
      real*8 ::  Xdat(IDAT), Ydat(JDAT)
      real*4 ::  Xdat0(IDAT), Ydat0(JDAT)
      Integer*2,dimension(IDAT,JDAT) ::  PB0
      Real*4,dimension(IDAT,JDAT) ::  PB4
      Real*4 attrib
      Real*8,dimension(IDAT,JDAT) ::  PB
      Real*8 offset,scale
      character*80 input_path
      common/input_data/input_path
      character*4 cyear

      save

      logical firstcall, newyear, flat
      integer nu_diag
      common/iodiag/nu_diag
      save firstcall
      data firstcall/.true./,nold/-1/,dyearold/-1/,dindex
     */0,31,59,90,120,151,181,212,243,273,304,334/

      flat=.false.
      dyear=int(idate/10000)
      newyear=dyearold.ne.dyear
      
      dmonth=int((idate-10000*dyear)/100)
      dday=int(idate-10000*dyear-100*dmonth)
      
      n=dindex(dmonth)+dday
c      if(4*int(dyear/4).eq.4.and.dmonth.gt.2)n=n+1
      
c      Write(nu_diag,*)'In RHumid'
      
	if(firstcall)then
	
      allocate(rxbox(imt,jmt),rybox(imt,jmt),i0box(imt,jmt),
     *   j0box(imt,jmt),i1box(imt,jmt),j1box(imt,jmt),nodetype(imt,jmt),
     *   PP(imt,jmt))

         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'rhum/rhum.sig995.'//
     *      cyear//'.nc'
         ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'RHum Error: reading NetCDF file'''//
     *         trim(fname)//''''
            call exit_coupler
         endif

      beg(1)=1
      beg(2)=1
      beg(3)=1
      cnt(1)=IDAT
      cnt(2)=1
      cnt(3)=1
      
      ierr=nf_inq_varid(ncid,'lon',vid)
      ierr=nf_get_vara_real(ncid,vid,beg,cnt,Xdat0(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'RHum Error: Error reading lon'
         call exit_coupler
      endif

      cnt(1)=JDAT

      ierr=nf_inq_varid(ncid,'lat',vid)
      ierr=nf_get_vara_real(ncid,vid,beg,cnt,Ydat0(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'RHum Error: Error reading lat'
         call exit_coupler
      endif
      ierr=nf_close(ncid)

	   do i=1,IDAT
	      Xdat(i)=dble(Xdat0(i))
	   enddo
	   do j=1,JDAT
	      Ydat(j)=dble(Ydat0(JDAT-j+1))
           enddo

	   do i=1,imt
	      do j=1,jmt
	         xx=X(i,j)
		 if(xx.lt.0)xx=xx+360.0
		 yy=Y(i,j)
		 mm=mask(i,j)
		    i0=0
		    do ii=1,IDAT-1
		       if(Xdat(ii).le.xx.and.Xdat(ii+1).ge.xx)then
		          i0=ii
			  i1=ii+1
			  rxbox(i,j)=(xx-Xdat(i0))/(Xdat(i1)-Xdat(i0))
			  exit
		       endif
		    enddo
		    if(i0.eq.0)then
		       if(Xdat(IDAT).le.xx.and.(Xdat(1)+360).ge.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0))/
     *  		    (Xdat(i1)+360-Xdat(i0))
                   elseif(Xdat(1).ge.xx.and.(Xdat(IDAT)-360).le.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0)+360)/
     *  		    (Xdat(i1)+360-Xdat(i0))
		       else
                          Write(nu_diag,*) 'Rhumid: Grid from NetCDF is not covering'
			  Write(nu_diag,*) xx,Xdat(1),Xdat(IDAT)
                          call exit_coupler
		       endif
		    endif
		    j0=0
		    if(yy.le.Ydat(2))then
		       j0=2
		       j1=2
		       nodetype(i,j)=1
                       rybox(i,j)=(yy+90.0)/(Ydat(j1)+90.0)
		    else if(yy.ge.Ydat(JDAT-1))then
		       j0=JDAT-1
		       j1=JDAT-1
		       nodetype(i,j)=1
                       rybox(i,j)=(yy-Ydat(j0))/(90.0-Ydat(j0))
		    else
		    do jj=1,JDAT-1
		       if(Ydat(jj).le.yy.and.Ydat(jj+1).ge.yy)then
		          j0=jj
			  j1=jj+1
			  rybox(i,j)=(yy-Ydat(j0))/(Ydat(j1)-Ydat(j0))
			  exit
		       endif
		    enddo
		    nodetype(i,j)=0
		    endif
		    i0box(i,j)=i0
		    i1box(i,j)=i1
		    j0box(i,j)=j0
		    j1box(i,j)=j1
	      enddo
	   enddo
	   
        endif
	
      if(n.ne.nold)then
      
      if(newyear)then

         if(dyearold.ge.0)then      
	   ierr=nf_close(ncid)
	 endif
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'rhum/rhum.sig995.'//
     *      cyear//'.nc'
c	 Write(nu_diag,*)'Reading file: ',fname
	   ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
	   if(ierr.NE.NF_NOERR)then
              Write(nu_diag,*) 'RHum Error: reading NetCDF file'
              call exit_coupler
c           else
c              Write(nu_diag,'(a,i3)')'Open CDF ', ncid
           endif
	   dyearold=dyear
      endif

c      Write(nu_diag,'(a,i3)')'Read record no.',n,' from CDF ', ncid

      beg(1)=1
      beg(2)=1
      beg(3)=n
      cnt(1)=IDAT
      cnt(2)=JDAT
      cnt(3)=1
      
      ierr=nf_inq_varid(ncid,'rhum',vid)
      if(dyear.lt.2014)then
         ierr=nf_get_att_real(ncid,vid,'add_offset',attrib)
         offset=dble(attrib)
         ierr=nf_get_att_real(ncid,vid,'scale_factor',attrib)
         scale=dble(attrib)
c         Write(nu_diag,*)'Offset=',offset
c         Write(nu_diag,*)'Factor=',scale
         ierr=nf_get_vara_int2(ncid,vid,beg,cnt,PB0(1,1))

         Do i=1,IDAT
            Do j=1,JDAT
               PB(i,j)=dble(PB0(i,JDAT-j+1))
            Enddo
         Enddo
      else
         offset=0.0
         scale=1.0
         ierr=nf_get_vara_real(ncid,vid,beg,cnt,PB4(1,1))

         Do i=1,IDAT
            Do j=1,JDAT
               PB(i,j)=dble(PB4(i,JDAT-j+1))
            Enddo
         Enddo
      endif

c      Write(nu_diag,'(a,f7.2,a)')'Read Control: P(115,31)= ',
c     *   (PB(115,31)*scale+offset)/100.0,' hPa'

      PMIN=1e+10
      PMAX=-1e+10
      do j=1,jmt
         do i=1,imt
               i0=i0box(i,j)
	       i1=i1box(i,j)
	       j0=j0box(i,j)
	       j1=j1box(i,j)
	       rx0=rxbox(i,j)
	       rx1=1.0-rx0
	       ry0=rybox(i,j)
	       ry1=1.0-ry0
              if(nodetype(i,j).eq.0)then
               PP(i,j)= scale*(
     *            ry1*(rx1*PB(i0,j0)+rx0*PB(i1,j0))+
     *            ry0*(rx1*PB(i0,j1)+rx0*PB(i1,j1))) +
     *            offset
              else
               PP(i,j)=0.0
	       rr=1.0/IDAT
               if(flat)then
	       do ii=1,IDAT
                  PP(i,j)=PP(i,j)+rr*PB(ii,j0)
	       enddo
	       PP(i,j)=scale*PP(i,j)+offset
               else
c        Components at the pole as mean values
                  PPOLE=0.0
                  do ii=1,IDAT
                     uur=PB(ii,j0)*scale+offset
                     PPOLE=PPOLE+rr*uur
                  enddo

                  uu0=PB(i0,j0)*scale+offset
                  uu1=PB(i1,j0)*scale+offset

                  uu=rx1*uu0+rx0*uu1

                  PP(i,j)=ry1*uu+ry0*PPOLE
               endif
            endif
	    if(PMAX.lt.PP(i,j))PMAX=PP(i,j)
	    if(PMIN.gt.PP(i,j))PMIN=PP(i,j)
         enddo
      enddo

c      Write(nu_diag,*)'Max pressure: ', PMAX/100.0
c      Write(nu_diag,*)'Min pressure: ', PMIN/100.0
      
      endif
      
      nold=n
      
      do j=1,jmt
         do i=1,imt
	    PBOT(i,j)=0.01*PP(i,j)
         enddo
      enddo

      
      firstcall=.false.

      return
      end

      subroutine ciaf_pressure(idate, sec, imt, jmt, X, Y, mask, PBOT, 
     *                       vers)
      
      include 'netcdf.inc'
      
      integer*4, intent(in) :: idate, imt, jmt, mask(imt,jmt), sec, vers
      real*8, intent(in) :: X(imt,jmt), Y(imt,jmt)
      real*8, intent(out) :: PBOT(imt,jmt)
      real*8, allocatable, dimension(:,:) ::
     &  rxbox,rybox,PP
      integer, allocatable, dimension(:,:) ::
     &  i0box,j0box, i1box,j1box, nodetype
      integer dindex(12)
      integer i,j,n,ncid,vid,did,beg(3),cnt(3),ierr,nold
      integer dday, dmonth, dyear, dyearold, dyear0
      character*80 attstr
      character*80 fname
      character, allocatable :: str(:)
      integer,parameter :: IDAT=192
      integer,parameter :: JDAT=94
      real*8,parameter :: XSTEPDAT=360.0/IDAT
      real*8,parameter :: YSTEPDAT=2*(90.0/(JDAT+1))
      real*8,parameter :: pi=3.14159265/180.0
      real*8 ::  Xdat(IDAT), Ydat(JDAT)
      Real*4,dimension(IDAT,JDAT) ::  PB
      Real*4 offset,scale
      character*80 input_path
      common/input_data/input_path
      character*4 cyear

      save

      logical firstcall, newyear, flat
      integer nu_diag
      common/iodiag/nu_diag
      save firstcall
      data firstcall/.true./,nold/-1/,dyearold/-1/,dindex
     */0,31,59,90,120,151,181,212,243,273,304,334/

      flat=.false.
      dyear0=int(idate/10000)
c  Data are started from 1958, so all previous years are considered as 1958
      dyear=max(dyear0,1958)
      newyear=dyearold.ne.dyear
      
      dmonth=int((idate-10000*dyear0)/100)
      dday=int(idate-10000*dyear0-100*dmonth)
      
      n=(dindex(dmonth)+dday-1)*4+int((1.0*sec+3*3600)/(6.0*3600.0))
      n=max(n,1)
      n=min(n,1460)
c      if(4*int(dyear/4).eq.4.and.dmonth.gt.2)n=n+1
      
c      Write(nu_diag,*)'In CIAF AirHum'
      
	if(firstcall)then
	
      allocate(rxbox(imt,jmt),rybox(imt,jmt),i0box(imt,jmt),
     *   j0box(imt,jmt),i1box(imt,jmt),j1box(imt,jmt),nodetype(imt,jmt),
     *   PP(imt,jmt))

      if(vers.eq.1)then
         if(dyear.lt.1958)dyear=1958
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_1.0/slp.'//
     *      cyear//'.2004_08_03.nc'
         if(dyear.gt.2000)
     *      fname=trim(input_path)//'CIAF_1.0/slp.'//
     *      cyear//'.2006_04_06.nc'
      elseif(vers.eq.2)then
         if(dyear.lt.1948)dyear=1948
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_2.0/slp.'//
     *      cyear//'.05APR2010.nc'
      elseif(vers.eq.0)then
         fname=trim(input_path)//'CIAF_2.0/slp.'//
     *      '15JUNE2009.nc'
      endif
         ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
         if(ierr.NE.NF_NOERR)then
            Write(nu_diag,*) 'CIAF Pressure Error: reading NetCDF file'
            call exit_coupler
         endif

      beg(1)=1
      beg(2)=1
      beg(3)=1
      cnt(1)=IDAT
      cnt(2)=1
      cnt(3)=1
      
      ierr=nf_inq_varid(ncid,'LON',vid)
      ierr=nf_get_vara_double(ncid,vid,beg,cnt,Xdat(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'CIAF Pressure Error: Error reading lon'
         call exit_coupler
      endif

      cnt(1)=JDAT

      ierr=nf_inq_varid(ncid,'LAT',vid)
      ierr=nf_get_vara_double(ncid,vid,beg,cnt,Ydat(1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'CIAF Pressure Error: Error reading lat'
         call exit_coupler
      endif
      ierr=nf_close(ncid)

	   do i=1,imt
	      do j=1,jmt
	         xx=X(i,j)
		 if(xx.lt.0)xx=xx+360.0
		 yy=Y(i,j)
		 mm=mask(i,j)
		    i0=0
		    do ii=1,IDAT-1
		       if(Xdat(ii).le.xx.and.Xdat(ii+1).ge.xx)then
		          i0=ii
			  i1=ii+1
			  rxbox(i,j)=(xx-Xdat(i0))/(Xdat(i1)-Xdat(i0))
			  exit
		       endif
		    enddo
		    if(i0.eq.0)then
		       if(Xdat(IDAT).le.xx.and.(Xdat(1)+360).ge.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0))/
     *  		    (Xdat(i1)+360-Xdat(i0))
                   elseif(Xdat(1).ge.xx.and.(Xdat(IDAT)-360).le.xx)then
		          i0=IDAT
			  i1=1
			  rxbox(i,j)=(xx-Xdat(i0)+360)/
     *  		    (Xdat(i1)+360-Xdat(i0))
		       else
                          Write(nu_diag,*) 'Ciaf_slp: Grid from NetCDF is not covering'
			  Write(nu_diag,*) xx,Xdat(1),Xdat(IDAT)
                          call exit_coupler
		       endif
		    endif
		    j0=0
		    if(yy.le.Ydat(2))then
		       j0=2
		       j1=2
		       nodetype(i,j)=1
                       rybox(i,j)=(yy+90.0)/(Ydat(j1)+90.0)
		    else if(yy.ge.Ydat(JDAT-1))then
		       j0=JDAT-1
		       j1=JDAT-1
		       nodetype(i,j)=1
                       rybox(i,j)=(yy-Ydat(j0))/(90.0-Ydat(j0))
		    else
		    do jj=1,JDAT-1
		       if(Ydat(jj).le.yy.and.Ydat(jj+1).ge.yy)then
		          j0=jj
			  j1=jj+1
			  rybox(i,j)=(yy-Ydat(j0))/(Ydat(j1)-Ydat(j0))
			  exit
		       endif
		    enddo
		    nodetype(i,j)=0
		    endif
		    i0box(i,j)=i0
		    i1box(i,j)=i1
		    j0box(i,j)=j0
		    j1box(i,j)=j1
	      enddo
	   enddo
	   
        endif
	
      if(n.ne.nold)then
      
      if(newyear)then

         if(dyearold.ge.0)then      
	   ierr=nf_close(ncid)
	 endif
      if(vers.eq.1)then
         if(dyear.lt.1958)dyear=1958
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_1.0/slp.'//
     *      cyear//'.2004_08_03.nc'
         if(dyear.gt.2000)
     *      fname=trim(input_path)//'CIAF_1.0/slp.'//
     *      cyear//'.2006_04_06.nc'
      elseif(vers.eq.2)then
         if(dyear.lt.1948)dyear=1948
         Write(cyear,'(i4)')dyear
         fname=trim(input_path)//'CIAF_2.0/slp.'//
     *      cyear//'.05APR2010.nc'
      elseif(vers.eq.0)then
         fname=trim(input_path)//'CIAF_2.0/slp.'//
     *      '15JUNE2009.nc'
      endif
c	 Write(nu_diag,*)'Reading file: ',fname
	   ierr=nf_open(trim(fname),NF_NOWRITE,ncid)
	   if(ierr.NE.NF_NOERR)then
              Write(nu_diag,*) 'AirHum Error: reading NetCDF file'
              call exit_coupler
c           else
c              Write(nu_diag,'(a,i3)')'Open CDF ', ncid
           endif
	   dyearold=dyear
      endif

c      Write(nu_diag,'(a,i4)')'Read record no.',n,' from CDF ', ncid

      beg(1)=1
      beg(2)=1
      beg(3)=n
      cnt(1)=IDAT
      cnt(2)=JDAT
      cnt(3)=1
      
      ierr=nf_inq_varid(ncid,'SLP',vid)
      offset=0.0
      scale=1.0
      ierr=nf_get_vara_real(ncid,vid,beg,cnt,PB(1,1))
      if(ierr.NE.NF_NOERR)then
         Write(nu_diag,*) 'CIAF Pressure Error: Error reading slp'
         call exit_coupler
      endif
      
c      Write(nu_diag,'(a,f6.2,a)')'Read Control: QA(115,31)= ',
c     *   PB(115,31)*scale+offset,' K'

      PMIN=1e+10
      PMAX=-1e+10
      do j=1,jmt
         do i=1,imt
             i0=i0box(i,j)
             i1=i1box(i,j)
             j0=j0box(i,j)
             j1=j1box(i,j)
             rx0=rxbox(i,j)
             rx1=1.0-rx0
             ry0=rybox(i,j)
             ry1=1.0-ry0
	       if(nodetype(i,j).eq.0)then
               PP(i,j)= scale*(
     *            ry1*(rx1*dble(PB(i0,j0))+rx0*dble(PB(i1,j0)))+
     *            ry0*(rx1*dble(PB(i0,j1))+rx0*dble(PB(i1,j1)))) +
     *            offset
              else
               PP(i,j)=0.0
	       rr=1.0/IDAT
               if(flat)then
	       do ii=1,IDAT
                  PP(i,j)=PP(i,j)+rr*dble(PB(ii,j0))
	       enddo
	       PP(i,j)=scale*PP(i,j)+offset
               else
c        Components at the pole as mean values
                  PPOLE=0.0
                  do ii=1,IDAT
                     uur=dble(PB(ii,j0))*scale+offset
                     PPOLE=PPOLE+rr*uur
                  enddo

                  uu0=dble(PB(i0,j0))*scale+offset
                  uu1=dble(PB(i1,j0))*scale+offset

                  uu=rx1*uu0+rx0*uu1

                  PP(i,j)=ry1*uu+ry0*PPOLE
               endif
            endif
	    if(PMAX.lt.PP(i,j))PMAX=PP(i,j)
	    if(PMIN.gt.PP(i,j))PMIN=PP(i,j)
         enddo
      enddo

c      Write(nu_diag,*)'Max humidity: ', PMAX
c      Write(nu_diag,*)'Min humidity: ', PMIN
c      Write(nu_diag,*)'Pole humidity: ', PP(88,136)
      
      else

c      Write(nu_diag,*)'Old values not updated.'
      
      endif
      
      nold=n
      
      do j=1,jmt
         do i=1,imt
	    PBOT(i,j)=PP(i,j)
         enddo
      enddo

      
      firstcall=.false.

      return
      end

