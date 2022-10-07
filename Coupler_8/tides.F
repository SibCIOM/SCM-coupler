!     $Id: tides.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $

      subroutine tides(idate, sec, imt, jmt, mask, X,Y, 
     *   PBOT, UBOT, VBOT, fname)
      
      integer*4, intent(in) :: idate, imt, jmt, mask(imt,jmt), sec
      real*8, intent(in) :: X(imt,jmt),Y(imt,jmt)
      real*8, intent(out) :: PBOT(imt,jmt),UBOT(imt,jmt),VBOT(imt,jmt)

      integer dday, dmonth, dyear, dindex(12),constitNum_data(21),
     * ispec_data(21),c(10),ispec,constitNum
      character*(*) fname
      real*8, allocatable, dimension(:,:,:) ::
     &  amp,gph,ampu,gphu,ampv,gphv
      real*8,parameter :: pi=3.14159265358979/180.0
      real*8 alpha_data(21),omega_data(21),phase_data(22),amp_data(21)
     *,  hcr(8),hci(8)

      save

      logical firstcall
      integer nu_diag
      common/iodiag/nu_diag
      save firstcall
      data firstcall/.true./,dindex
     */0,31,59,90,120,151,181,212,243,273,304,334/
      data constitNum_data/1,2,5,6,3,7,4,8,0,0,0,0,0,0,0,0,0,0,0,0,0/,
     *ispec_data/2,2,1,1,2,1,2,1,2,2,2,2,2,1,1,1,1,0,0,0,0/,
     *alpha_data/0.693,0.693,0.736,0.695,0.693,0.706,0.693,0.695,
     *   0.693,0.693,0.693,0.693,0.693,0.693,0.693,0.693,
     *   0.693,0.693,0.693,0.693,0.693/,
     *omega_data/1.405189e-04,1.454441e-04,7.292117e-05,6.759774e-05,
     *   1.378797e-04,7.252295e-05,1.458423e-04,6.495854e-05,
     *   1.352405e-04,1.355937e-04,1.382329e-04,1.431581e-04,
     *   1.452450e-04,7.556036e-05,7.028195e-05,7.824458e-05,
     *   6.531174e-05,0.053234e-04,0.026392e-04,0.003982e-04,
     *   2.81038e-04/,
     *phase_data/1.731557546,0.000000000,0.173003674,1.558553872,
     *   6.050721243,6.110181633,3.487600001,5.877717569,
     *   4.086699633,3.463115091,5.427136701,0.553986502,
     *   0.052841931,2.137025284,2.436575100,1.929046130,
     *   5.254133027,1.756042456,1.964021610,3.487600001,
     *   3.463115091,3.46372e00/,
     *amp_data/0.2441,0.112743,0.141565,0.100661,0.046397,0.046848,
     *   0.030684,0.019273,0.006141,0.007408,0.008811,0.006931,
     *   0.006608,0.007915,0.007915,0.004338,0.003661,0.042041,
     *   0.022191,0.019567,0./,
     *c/1,2,3,4,5,6,7,8,18,19/,ncc/8/

      dyear=int(idate/10000)
      
      dmonth=int((idate-10000*dyear)/100)
      dday=int(idate-10000*dyear-100*dmonth)
      
      rday=693962.0+(dyear-1900)*365.0+dble(int((dyear-1901)*0.25))
      rday0=693962.0+92*365.0+dble(int(91*0.25))
      rday=rday+dble(dindex(dmonth)+dday)+dble(sec)/86400.0
      if(4*int(dyear/4).eq.4.and.dyear.ne.1900.and.dmonth.gt.2)
     *   rday=rday+1.0
      time=rday-rday0
      
      Write(nu_diag,*)'In tides: day=',time
      
	if(firstcall)then
	
         allocate(amp(10,imt,jmt),ampu(10,imt,jmt),ampv(10,imt,jmt),
     *      gph(10,imt,jmt),gphu(10,imt,jmt),gphv(10,imt,jmt))

c         fname='../../MassStore/tides/tides.dat'
         open(792,file=trim(fname),form='unformatted',access='direct',
     *      recl=(imt-2)*(jmt-2)*8,status='old')

         do nc=1,10
            read(792,rec=(nc*2-1))((amp(nc,i,j),i=2,imt-1),
     *         j=2,jmt-1)
            read(792,rec=(nc*2))((gph(nc,i,j),i=2,imt-1),
     *         j=2,jmt-1)
         enddo
         do nc=1,10
            read(792,rec=(20+nc*4-3))((ampu(nc,i,j),i=2,imt-1),
     *         j=2,jmt-1)
            read(792,rec=(20+nc*4-2))((gphu(nc,i,j),i=2,imt-1),
     *         j=2,jmt-1)
            read(792,rec=(20+nc*4-1))((ampv(nc,i,j),i=2,imt-1),
     *         j=2,jmt-1)
            read(792,rec=(20+nc*4))((gphv(nc,i,j),i=2,imt-1),
     *         j=2,jmt-1)
         enddo

         close(792)

         firstcall=.false.

      endif

      do j=1,jmt
         do i=1,imt
            PBOT(i,j)=0.0
            if(mask(i,j).gt.0)then
               hhat=0.0
               do nc=1,ncc
                  kk=c(nc)
                  constitNum=constitNum_data(kk)
                  ispec=ispec_data(kk)
                  am=amp_data(kk)
                  alpha=alpha_data(kk)
                  omega=omega_data(kk)
                  ph=phase_data(kk)
                  call nodal(time+48622.0,nc,pu,pf)
                  hcr(kk)=amp(nc,i,j)*cos(-gph(nc,i,j)*pi)
                  hci(kk)=amp(nc,i,j)*sin(-gph(nc,i,j)*pi)
                  arg=pf*hcr(kk)*cos(omega*time*86400.0+ph+pu)-
     *               pf*hci(kk)*sin(omega*time*86400.0+ph+pu)
                  hhat=hhat+arg
               enddo
               call InferMinor(hcr,hci,time,dh)
               PBOT(i,j)=hhat+dh
            endif
         enddo
      enddo

      do j=1,jmt
         do i=1,imt
            UBOT(i,j)=0.0
            if(mask(i,j).gt.0)then
               hhat=0.0
               do nc=1,ncc
                  kk=c(nc)
                  constitNum=constitNum_data(kk)
                  ispec=ispec_data(kk)
                  am=amp_data(kk)
                  alpha=alpha_data(kk)
                  omega=omega_data(kk)
                  ph=phase_data(kk)
                  call nodal(time+48622.0,nc,pu,pf)
                  hcr(kk)=ampu(nc,i,j)*cos(-gphu(nc,i,j)*pi)
                  hci(kk)=ampu(nc,i,j)*sin(-gphu(nc,i,j)*pi)
                  arg=pf*hcr(kk)*cos(omega*time*86400.0+ph+pu)-
     *               pf*hci(kk)*sin(omega*time*86400.0+ph+pu)
                  hhat=hhat+arg
               enddo
               call InferMinor(hcr,hci,time,dh)
               UBOT(i,j)=(hhat+dh)*0.01
            endif
         enddo
      enddo

      do j=1,jmt
         do i=1,imt
            VBOT(i,j)=0.0
            if(mask(i,j).gt.0)then
               hhat=0.0
               do nc=1,ncc
                  kk=c(nc)
                  constitNum=constitNum_data(kk)
                  ispec=ispec_data(kk)
                  am=amp_data(kk)
                  alpha=alpha_data(kk)
                  omega=omega_data(kk)
                  ph=phase_data(kk)
                  call nodal(time+48622.0,nc,pu,pf)
                  hcr(kk)=ampv(nc,i,j)*cos(-gphv(nc,i,j)*pi)
                  hci(kk)=ampv(nc,i,j)*sin(-gphv(nc,i,j)*pi)
                  arg=pf*hcr(kk)*cos(omega*time*86400.0+ph+pu)-
     *               pf*hci(kk)*sin(omega*time*86400.0+ph+pu)
                  hhat=hhat+arg
               enddo
               call InferMinor(hcr,hci,time,dh)
               VBOT(i,j)=(hhat+dh)*0.01
            endif
         enddo
      enddo
      
      do j=1,jmt
         do i=1,imt
            if(mask(i,j).gt.0)then
               if(Y(i,j).gt.60)then
                  sn=sin(X(i,j)*pi)
                  cs=cos(X(i,j)*pi)
                  if(i.ne.1.and.i.ne.imt)then
                     xx1=(90-Y(i+1,j))*sin(X(i+1,j)*pi)
                     yy1=-(90-Y(i+1,j))*cos(X(i+1,j)*pi)
                     xx2=(90-Y(i-1,j))*sin(X(i-1,j)*pi)
                     yy2=-(90-Y(i-1,j))*cos(X(i-1,j)*pi)
                  else
                     if(i.eq.1)then
                        xx1=(90-Y(i+1,j))*sin(X(i+1,j)*pi)
                        yy1=-(90-Y(i+1,j))*cos(X(i+1,j)*pi)
                        xx2=(90-Y(i,j))*sn
                        yy2=-(90-Y(i,j))*cs
                     else
                        xx1=(90-Y(i,j))*sn
                        yy1=-(90-Y(i,j))*cs
                        xx2=(90-Y(i-1,j))*sin(X(i-1,j)*pi)
                        yy2=-(90-Y(i-1,j))*cos(X(i-1,j)*pi)
                     endif
                  endif
                  dx=xx1-xx2
                  dy=yy1-yy2
                  alph=atan2(dy,dx)
                  uu0=UBOT(i,j)*cs-VBOT(i,j)*sn
                  vv0=UBOT(i,j)*sn+VBOT(i,j)*cs
                  UBOT(i,j)=uu0*cos(alph)+vv0*sin(alph)
                  VBOT(i,j)=-uu0*sin(alph)+vv0*cos(alph)
                  alph=0.0
               endif
            endif
         enddo
      enddo

      return
      end

      subroutine nodal(time,cid,pu,pf);
c
      real*8 arg(53),f(53),u(53)
      integer cid
      integer index(21)
      data index/
     *30,35,19,12,27,17,37,10,25,26,28,33,34,23,14,24,11,5,3,2,45/

c     Determine equilibrium arguments
c     -------------------------------
      pp=282.94 ! solar perigee at epoch 2000
      rad=3.14159265358979/180
      T=time-51544.4993
      s=mod(218.3164+13.17639648*T,360.0)
      h=mod(280.4661+0.98564736*T,360.0)
      p=mod(83.3535+0.11140353*T,360.0)
      omega=mod(125.0445-0.05295377*T,360.0)
      if(s.lt.0.0)s=s+360
      if(h.lt.0.0)h=h+360
      if(p.lt.0.0)p=p+360
      if(omega.lt.0.0)omega=omega+360
      hour = (time - int(time))*24.0
      t1 = 15.0*hour
      t2 = 30.0*hour

      arg( 1) = h - pp                        ! Sa
      arg( 2) = 2*h                           ! Ssa
      arg( 3) = s - p                         ! Mm
      arg( 4) = 2*s - 2*h                     ! MSf
      arg( 5) = 2*s                           ! Mf
      arg( 6) = 3*s - p                       ! Mt
      arg( 7) = t1 - 5*s + 3*h + p - 90       ! alpha1
      arg( 8) = t1 - 4*s + h + 2*p - 90       ! 2Q1
      arg( 9) = t1 - 4*s + 3*h - 90           ! sigma1
      arg(10) = t1 - 3*s + h + p - 90         ! q1
      arg(11) = t1 - 3*s + 3*h - p - 90       ! rho1
      arg(12) = t1 - 2*s + h - 90             ! o1
      arg(13) = t1 - 2*s + 3*h + 90           ! tau1
      arg(14) = t1 - s + h + 90               ! M1
      arg(15) = t1 - s + 3*h - p + 90         ! chi1
      arg(16) = t1 - 2*h + pp - 90            ! pi1
      arg(17) = t1 - h - 90                   ! p1
      arg(18) = t1 + 90                       ! s1
      arg(19) = t1 + h + 90                   ! k1
      arg(20) = t1 + 2*h - pp + 90            ! psi1
      arg(21) = t1 + 3*h + 90                 ! phi1
      arg(22) = t1 + s - h + p + 90           ! theta1
      arg(23) = t1 + s + h - p + 90           ! J1
      arg(24) = t1 + 2*s + h + 90             ! OO1
      arg(25) = t2 - 4*s + 2*h + 2*p          ! 2N2
      arg(26) = t2 - 4*s + 4*h                ! mu2
      arg(27) = t2 - 3*s + 2*h + p            ! n2
      arg(28) = t2 - 3*s + 4*h - p            ! nu2
      arg(29) = t2 - 2*s + h + pp             ! M2a
      arg(30) = t2 - 2*s + 2*h                ! M2
      arg(31) = t2 - 2*s + 3*h - pp           ! M2b
      arg(32) = t2 - s + p + 180.             ! lambda2
      arg(33) = t2 - s + 2*h - p + 180.       ! L2
      arg(34) = t2 - h + pp                   ! T2
      arg(35) = t2                            ! S2
      arg(36) = t2 + h - pp + 180             ! R2
      arg(37) = t2 + 2*h                      ! K2
      arg(38) = t2 + s + 2*h - pp             ! eta2
      arg(39) = t2 - 5*s + 4.0*h + p          ! MNS2
      arg(40) = t2 + 2*s - 2*h                ! 2SM2
      arg(41) = 1.5*arg(30)                   ! M3
      arg(42) = arg(19) + arg(30)             ! MK3
      arg(43) = 3*t1                          ! S3
      arg(44) = arg(27) + arg(30)             ! MN4
      arg(45) = 2*arg(30)                     ! M4
      arg(46) = arg(30) + arg(35)             ! MS4
      arg(47) = arg(30) + arg(37)             ! MK4
      arg(48) = 4*t1                          ! S4
      arg(49) = 5*t1                          ! S5
      arg(50) = 3*arg(30)                     ! M6
      arg(51) = 3*t2                          ! S6
      arg(52) = 7.0*t1                        ! S7
      arg(53) = 4*t2                          ! S8
!
!     determine nodal corrections f and u 
!     -----------------------------------
      sinn = sin(omega*rad);
      cosn = cos(omega*rad);
      sin2n = sin(2*omega*rad);
      cos2n = cos(2*omega*rad);
      sin3n = sin(3*omega*rad);
!!
      f( 1) = 1                                        ! Sa
      f( 2) = 1                                        ! Ssa
      f( 3) = 1 - 0.130*cosn                           ! Mm
      f( 4) = 1                                        ! MSf
      f( 5) = 1.043 + 0.414*cosn                       ! Mf
      f( 6) = sqrt((1+.203*cosn+.040*cos2n)**2 +
     *   (.203*sinn+.040*sin2n)**2)                    ! Mt
      f( 7) = 1                                        ! alpha1
      f( 8) = sqrt((1.+.188*cosn)**2+(.188*sinn)**2)   ! 2Q1
      f( 9) = f(8)                                     ! sigma1
      f(10) = f(8)                                     ! q1
      f(11) = f(8)                                     ! rho1
      f(12) = sqrt((1.0+0.189*cosn-0.0058*cos2n)**2 +
     *   (0.189*sinn-0.0058*sin2n)**2)                 ! O1
      f(13) = 1                                        ! tau1
! tmp1  = 2.*cos(p*rad)+.4*cos((p-omega)*rad);
! tmp2  = sin(p*rad)+.2*sin((p-omega)*rad);            ! Doodson's
      tmp1 = 1.36*cos(p*rad)+.267*cos((p-omega)*rad)   ! Ray's
      tmp2 = 0.64*sin(p*rad)+.135*sin((p-omega)*rad)
      f(14) = sqrt(tmp1**2 + tmp2**2)                  ! M1
      f(15) = sqrt((1.+.221*cosn)**2+(.221*sinn)**2)   ! chi1
      f(16) = 1                                        ! pi1
      f(17) = 1                                        ! P1
      f(18) = 1                                        ! S1
      f(19) = sqrt((1.+.1158*cosn-.0029*cos2n)**2 +
     *   (.1554*sinn-.0029*sin2n)**2)                  ! K1
      f(20) = 1                                        ! psi1
      f(21) = 1                                        ! phi1
      f(22) = 1                                        ! theta1
      f(23) = sqrt((1.+.169*cosn)**2+(.227*sinn)**2)   ! J1
      f(24) = sqrt((1.0+0.640*cosn+0.134*cos2n)**2 +
     *   (0.640*sinn+0.134*sin2n)**2 )                 ! OO1
      f(25) = sqrt((1.-.03731*cosn+.00052*cos2n)**2 +
     *   (.03731*sinn-.00052*sin2n)**2)                ! 2N2
      f(26) = f(25)                                    ! mu2
      f(27) = f(25)                                    ! N2
      f(28) = f(25)                                    ! nu2
      f(29) = 1                                        ! M2a
      f(30) = f(25)                                    ! M2
      f(31) = 1                                        ! M2b
      f(32) = 1                                        ! lambda2
      temp1 = 1.-0.25*cos(2*p*rad)-0.11*cos((2*p-omega)*rad)-0.04*cosn
      temp2 = 0.25*sin(2*p*rad)+0.11*sin((2*p-omega)*rad)+ 0.04*sinn
      f(33) = sqrt(temp1**2 + temp2**2)                ! L2
      f(34) = 1                                        ! T2
      f(35) = 1                                        ! S2
      f(36) = 1                                        ! R2
      f(37) = sqrt((1.+.2852*cosn+.0324*cos2n)**2 +
     *   (.3108*sinn+.0324*sin2n)**2)                  ! K2
      f(38) = sqrt((1.+.436*cosn)**2+(.436*sinn)**2)   ! eta2
      f(39) = f(30)**2                                 ! MNS2
      f(40) = f(30)                                    ! 2SM2
      f(41) = 1      ! wrong                           ! M3
      f(42) = f(19)*f(30)                              ! MK3
      f(43) = 1                                        ! S3
      f(44) = f(30)**2                                 ! MN4
      f(45) = f(44)                                    ! M4
      f(46) = f(44)                                    ! MS4
      f(47) = f(30)*f(37)                              ! MK4
      f(48) = 1                                        ! S4
      f(49) = 1                                        ! S5
      f(50) = f(30)**3                                 ! M6
      f(51) = 1                                        ! S6
      f(52) = 1                                        ! S7
      f(53) = 1                                        ! S8
!
      u( 1) = 0                                       ! Sa
      u( 2) = 0                                       ! Ssa
      u( 3) = 0                                       ! Mm
      u( 4) = 0                                       ! MSf
      u( 5) = -23.7*sinn + 2.7*sin2n - 0.4*sin3n      ! Mf
      u( 6) = atan(-(.203*sinn+.040*sin2n)/
     *   (1+.203*cosn+.040*cos2n))/rad                ! Mt
      u( 7) = 0                                       ! alpha1
      u( 8) = atan(.189*sinn/(1.+.189*cosn))/rad      ! 2Q1
      u( 9) = u(8)                                    ! sigma1
      u(10) = u(8)                                    ! q1
      u(11) = u(8)                                    ! rho1
      u(12) = 10.8*sinn - 1.3*sin2n + 0.2*sin3n       ! O1
      u(13) = 0                                       ! tau1
      u(14) = atan2(tmp2,tmp1)/rad                    ! M1
      u(15) = atan(-.221*sinn/(1.+.221*cosn))/rad     ! chi1
      u(16) = 0                                       ! pi1
      u(17) = 0                                       ! P1
      u(18) = 0                                       ! S1
      u(19) = atan((-.1554*sinn+.0029*sin2n)/
     *   (1.+.1158*cosn-.0029*cos2n))/rad             ! K1
      u(20) = 0                                       ! psi1
      u(21) = 0                                       ! phi1
      u(22) = 0                                       ! theta1
      u(23) = atan(-.227*sinn/(1.+.169*cosn))/rad     ! J1
      u(24) = atan(-(.640*sinn+.134*sin2n)/
     *   (1.+.640*cosn+.134*cos2n))/rad               ! OO1
      u(25) = atan((-.03731*sinn+.00052*sin2n)/ 
     *   (1.-.03731*cosn+.00052*cos2n))/rad           ! 2N2
      u(26) = u(25)                                   ! mu2
      u(27) = u(25)                                   ! N2
      u(28) = u(25)                                   ! nu2
      u(29) = 0                                       ! M2a
      u(30) = u(25)                                   ! M2
      u(31) = 0                                       ! M2b
      u(32) = 0                                       ! lambda2
      u(33) = atan(-temp2/temp1)/rad                  ! L2
      u(34) = 0                                       ! T2
      u(35) = 0                                       ! S2
      u(36) = 0                                       ! R2
      u(37) = atan(-(.3108*sinn+.0324*sin2n)/ 
     *   (1.+.2852*cosn+.0324*cos2n))/rad             ! K2
      u(38) = atan(-.436*sinn/(1.+.436*cosn))/rad     ! eta2
      u(39) = u(30)*2                                 ! MNS2
      u(40) = u(30)                                   ! 2SM2
      u(41) = 1.5d0*u(30)                             ! M3
      u(42) = u(30) + u(19)                           ! MK3
      u(43) = 0                                       ! S3
      u(44) = u(30)*2                                 ! MN4
      u(45) = u(44)                                   ! M4
      u(46) = u(30)                                   ! MS4
      u(47) = u(30)+u(37)                             ! MK4
      u(48) = 0                                       ! S4
      u(49) = 0                                       ! S5
      u(50) = u(30)*3                                 ! M6
      u(51) = 0                                       ! S6
      u(52) = 0                                       ! S7
      u(53) = 0                                       ! S8
! set correspondence between given constituents and supported in OTIS
      pu=u(index(cid))
      pf=f(index(cid))
      return
      end

! Lana Erofeeva, re-make for matlab OCT 2004
! usage:
! [dh]=InferMinor(zmaj,cid,SDtime)
!
! Based on Richard Ray's code perth2
! Return correction for the 16 minor constituents
! zeros, if not enough input constituents for inference
! Input:
! cid(ncon,4)   - GIVEN constituents
! zmaj(ncon, - Complex HC for GIVEN constituents/points
! SDtime(    - time expressed in Serial Days (see help datenum)
! Modes:
!      Time series: zmaj(ncon,1),  SDtime(nt,1)  ->Output: dh(nt,1)
!      Drift Track: zmaj(ncon,nt), SDtime(nt,1)  ->Output: dh(nt,1)
!      Map:         zmaj(ncon,N,M),SDtime(1,1)   ->Output: dh(N,M)
! -------------------------------------------------------------------
      subroutine InferMinor(zmajr,zmaji,time,dh)

      real*8 zmajr(8),zmaji(8),zr8(8),zi8(8)
      real*8 zminr(18),zmini(18),arg(18),f(18),u(18)
      integer index(8)

      data index/6,7,4,2,5,3,8,1/

      rad=3.14159265358979/180
      PP=282.8
      do k=1,8
         zr8(index(k))=zmajr(k)
         zi8(index(k))=zmaji(k)
      enddo
      dh=0.0
      time_mjd=48622+time
      zminr(1)  = 0.263 *zr8(1) - 0.0252*zr8(2)   !2Q1
      zminr(2)  = 0.297 *zr8(1) - 0.0264*zr8(2)   !sigma1
      zminr(3)  = 0.164 *zr8(1) + 0.0048*zr8(2)   !rho1 +
      zminr(4)  = 0.0140*zr8(2) + 0.0101*zr8(4)   !M1
      zminr(5)  = 0.0389*zr8(2) + 0.0282*zr8(4)   !M1
      zminr(6)  = 0.0064*zr8(2) + 0.0060*zr8(4)   !chi1
      zminr(7)  = 0.0030*zr8(2) + 0.0171*zr8(4)   !pi1
      zminr(8)  =-0.0015*zr8(2) + 0.0152*zr8(4)   !phi1
      zminr(9)  =-0.0065*zr8(2) + 0.0155*zr8(4)   !theta1
      zminr(10) =-0.0389*zr8(2) + 0.0836*zr8(4)   !J1 +
      zminr(11) =-0.0431*zr8(2) + 0.0613*zr8(4)   !OO1 +
      zminr(12) = 0.264 *zr8(5) - 0.0253*zr8(6)   !2N2 +
      zminr(13) = 0.298 *zr8(5) - 0.0264*zr8(6)   !mu2 +
      zminr(14) = 0.165 *zr8(5) + 0.00487*zr8(6)  !nu2 +
      zminr(15) = 0.0040*zr8(6) + 0.0074*zr8(7)   !lambda2
      zminr(16) = 0.0131*zr8(6) + 0.0326*zr8(7)   !L2 +
      zminr(17) = 0.0033*zr8(6) + 0.0082*zr8(7)   !L2 +
      zminr(18) = 0.0585*zr8(7)                   !t2 + 

      zmini(1)  = 0.263 *zi8(1) - 0.0252*zi8(2)   !2Q1
      zmini(2)  = 0.297 *zi8(1) - 0.0264*zi8(2)   !sigma1
      zmini(3)  = 0.164 *zi8(1) + 0.0048*zi8(2)   !rho1 +
      zmini(4)  = 0.0140*zi8(2) + 0.0101*zi8(4)   !M1
      zmini(5)  = 0.0389*zi8(2) + 0.0282*zi8(4)   !M1
      zmini(6)  = 0.0064*zi8(2) + 0.0060*zi8(4)   !chi1
      zmini(7)  = 0.0030*zi8(2) + 0.0171*zi8(4)   !pi1
      zmini(8)  =-0.0015*zi8(2) + 0.0152*zi8(4)   !phi1
      zmini(9)  =-0.0065*zi8(2) + 0.0155*zi8(4)   !theta1
      zmini(10) =-0.0389*zi8(2) + 0.0836*zi8(4)   !J1 +
      zmini(11) =-0.0431*zi8(2) + 0.0613*zi8(4)   !OO1 +
      zmini(12) = 0.264 *zi8(5) - 0.0253*zi8(6)   !2N2 +
      zmini(13) = 0.298 *zi8(5) - 0.0264*zi8(6)   !mu2 +
      zmini(14) = 0.165 *zi8(5) + 0.00487*zi8(6)  !nu2 +
      zmini(15) = 0.0040*zi8(6) + 0.0074*zi8(7)   !lambda2
      zmini(16) = 0.0131*zi8(6) + 0.0326*zi8(7)   !L2 +
      zmini(17) = 0.0033*zi8(6) + 0.0082*zi8(7)   !L2 +
      zmini(18) = 0.0585*zi8(7)                   !t2 + 
!
      hour = (time - int(time))*24.
      t1 = 15*hour
      t2 = 30*hour
      T=time-51544.4993
      S=mod(218.3164+13.17639648*T,360.0)
      H=mod(280.4661+0.98564736*T,360.0)
      P=mod(83.3535+0.11140353*T,360.0)
      omega=mod(125.0445-0.05295377*T,360.0)
      if(S.lt.0.0)S=S+360
      if(H.lt.0.0)H=H+360
      if(P.lt.0.0)P=P+360
      if(omega.lt.0.0)omega=omega+360
!
      arg(1) = t1 - 4*S + H + 2*P - 90.     ! 2Q1
      arg(2) = t1 - 4*S + 3*H - 90.         ! sigma1
      arg(3) = t1 - 3*S + 3*H - P - 90.     ! rho1
      arg(4) = t1 - S + H - P + 90.         ! M1
      arg(5) = t1 - S + H + P + 90.         ! M1
      arg(6) = t1 - S + 3*H - P + 90.       ! chi1
      arg(7) = t1 - 2*H + PP - 90.          ! pi1
      arg(8) = t1 + 3*H + 90.               ! phi1
      arg(9) = t1 + S - H + P + 90.         ! theta1
      arg(10) = t1 + S + H - P + 90.        ! J1
      arg(11) = t1 + 2*S + H + 90.          ! OO1
      arg(12) = t2 - 4*S + 2*H + 2*P        ! 2N2
      arg(13) = t2 - 4*S + 4*H              ! mu2
      arg(14) = t2 - 3*S + 4*H - P          ! nu2
      arg(15) = t2 - S + P + 180.D0         ! lambda2
      arg(16) = t2 - S + 2*H - P + 180.D0   ! L2
      arg(17) = t2 - S + 2*H + P            ! L2
      arg(18) = t2 - H + PP                  ! t2
!
!     determine nodal corrections f and u
      sinn = sin(omega*rad)
      cosn = cos(omega*rad)
      sin2n = sin(2*omega*rad)
      cos2n = cos(2*omega*rad)

      f(1) = sqrt((1.0 + 0.189*cosn - 0.0058*cos2n)**2 +
     *   (0.189*sinn - 0.0058*sin2n)**2)
      f(2) = f(1)
      f(3) = f(1)
      f(4) = sqrt((1.0 + 0.185*cosn)**2 + (0.185*sinn)**2)
      f(5) = sqrt((1.0 + 0.201*cosn)**2 + (0.201*sinn)**2)
      f(6) = sqrt((1.0 + 0.221*cosn)**2 + (0.221*sinn)**2)
      f(7) = 1.0
      f(8) = 1.0
      f(9) = 1.0
      f(10) = sqrt((1.0 + 0.198*cosn)**2 + (0.198*sinn)**2)
      f(11) = sqrt((1.0 + 0.640*cosn + 0.134*cos2n)**2 +
     *    (0.640*sinn + 0.134*sin2n)**2 )
      f(12) = sqrt((1.0 - 0.0373*cosn)**2 + (0.0373*sinn)**2)
      f(13) = f(12)
      f(14) = f(12)
      f(15) = 1.0
      f(16) = f(12)
      f(17) = sqrt((1.0 + 0.441*cosn)**2 + (0.441*sinn)**2)
      f(18) = 1.0
!
      u(1) = atan2(0.189*sinn - 0.0058*sin2n,
     *   1.0 + 0.189*cosn - 0.0058*sin2n)/rad
      u(2) = u(1)
      u(3) = u(1)
      u(4) = atan2( 0.185*sinn, 1.0 + 0.185*cosn)/rad
      u(5) = atan2(-0.201*sinn, 1.0 + 0.201*cosn)/rad
      u(6) = atan2(-0.221*sinn, 1.0 + 0.221*cosn)/rad
      u(7) = 0.0
      u(8) = 0.0
      u(9) = 0.0
      u(10) = atan2(-0.198*sinn, 1.0 + 0.198*cosn)/rad
      u(11) = atan2(-0.640*sinn - 0.134*sin2n,
     *   1.0 + 0.640*cosn + 0.134*cos2n)/rad
      u(12) = atan2(-0.0373*sinn, 1.0 - 0.0373*cosn)/rad
      u(13) = u(12)
      u(14) = u(12)
      u(15) = 0.0
      u(16) = u(12)
      u(17) = atan2(-0.441*sinn, 1.0 + 0.441*cosn)/rad
      u(18) = 0.0
!     sum over all tides
!     ------------------
      do k=1,18
         tmp=zminr(k)*f(k)*cos((arg(k) + u(k))*rad)-
     *       zmini(k)*f(k)*sin((arg(k)+u(k))*rad)
         dh = dh + tmp
      enddo
      return
      end
