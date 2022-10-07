!     $Id: conversion.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $

      module conversion
      use struct

      integer*2, allocatable, save, dimension(:,:,:,:) :: 
     *   idx_a2o,     ! node indexes of source grid to be used to derive
     *   idx_a2i,     ! interpolant at target grid node
     *   idx_a2l,
     *   idx_o2a,
     *   idx_o2i,
     *   idx_o2l,
     *   idx_i2a,
     *   idx_i2o,
     *   idx_i2l,
     *   idx_l2a,
     *   idx_l2o,
     *   idx_l2i

      real*8, allocatable, save, dimension(:,:,:) ::
     *   wt_a2o,      ! weights of nodes of source grid to be used to calculate
     *   wt_a2i,      ! interpolant at target grid node
     *   wt_a2l,
     *   wt_o2a,
     *   wt_o2i,
     *   wt_o2l,
     *   wt_i2a,
     *   wt_i2o,
     *   wt_i2l,
     *   wt_l2a,
     *   wt_l2o,
     *   wt_l2i

      real*8, allocatable, save, dimension(:,:) ::
     *   ang_a2o,     ! angles to rotate source grid lines to make them parallel
     *   ang_a2i,     ! to the target grid lines (locally) (positive is counterclockewise)
     *   ang_a2l,
     *   ang_o2a,
     *   ang_o2i,
     *   ang_o2l,
     *   ang_i2a,
     *   ang_i2o,
     *   ang_i2l,
     *   ang_l2a,
     *   ang_l2o,
     *   ang_l2i

      integer ncv,cnv_type
      integer 
     *   ncv_a2o, a2o_type, ncv_a2i, a2i_type, ncv_a2l, a2l_type,
     *   ncv_o2a, o2a_type, ncv_o2i, o2i_type, ncv_o2l, o2l_type,
     *   ncv_i2a, i2a_type, ncv_i2o, i2o_type, ncv_i2l, i2l_type,
     *   ncv_l2a, l2a_type, ncv_l2o, l2o_type, ncv_l2i, l2i_type

      integer, parameter :: 
     *   cnv_exp=1,   ! use exponential weights
     *   cnv_krig=2   ! use inverse-squared weights
      real*8, parameter :: yno=50.0  ! latitude to ... equatorial area and polar
      common/convrt/
     *   ncv,         ! number of neighbouring nodes to use in interpolation
     *   cnv_type,    ! type of weight-function (cnv_exp or cnv_krig)
     *   ncv_a2o, a2o_type,
     *   ncv_a2i, a2i_type,
     *   ncv_a2l, a2l_type,
     *   ncv_o2a, o2a_type,
     *   ncv_o2i, o2i_type,
     *   ncv_o2l, o2l_type,
     *   ncv_i2a, i2a_type,
     *   ncv_i2o, i2o_type,
     *   ncv_i2l, i2l_type,
     *   ncv_l2a, l2a_type,
     *   ncv_l2o, l2o_type,
     *   ncv_l2i, l2i_type

      contains

      Subroutine convert_scalar(typ, A, Agr, Ast, B, Bgr, Bst, Reset)

      Type (state) Ast, Bst
      Type (grid) Agr, Bgr
      real*8 A(Ast%imt,Ast%jmt), B(Bst%imt,Bst%jmt)
      logical, optional :: Reset
      logical first
      integer nu_diag, typ
      common/iodiag/nu_diag

      if (present(Reset)) then
         first=Reset
      else
         first=.false.
      endif

      if(ncv.eq.0)ncv=4
      if(cnv_type.eq.0)cnv_type=cnv_exp

      select case(typ)
         case(a2o)
            ncv = ncv_a2o
            cnv_type = a2o_type
            if(.not.allocated(idx_a2o).or.first)then
               if(.not.allocated(idx_a2o)) 
     *            allocate(idx_a2o(Bst%imt,Bst%jmt,ncv,2),
     *            wt_a2o(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Atmosphere --> Ocean'
               Open(183,file='a2o.con',form='binary',status='old',err=1)
               Read(183)imti,jmti,imtf,jmtf,ncvf
               if(.not.(Ast%imt.eq.imti.and.Ast%jmt.eq.jmti.and.
     *            Bst%imt.eq.imtf.and.Bst%jmt.eq.jmtf.and.ncv.eq.ncvf))then
                  write(*,*)'New grid convertion coeffs are needed'
                  Close(183)
                  goto 1
               endif
               Read(183)idx_a2o,wt_a2o
               Close(183)
               goto 2
    1          call coeffs(idx_a2o,wt_a2o,Agr,Ast,Bgr,Bst)
               Open(183,file='a2o.con',form='binary',status='unknown')
               Write(183)Ast%imt,Ast%jmt,Bst%imt,Bst%jmt,ncv,idx_a2o,wt_a2o
               Close(183)
    2          continue
            endif
            call cnv_scalar(idx_a2o, wt_a2o, A, Ast, B, Bst,Bgr)
         case(a2i)
            ncv = ncv_a2i
            cnv_type = a2i_type
            if(.not.allocated(idx_a2i).or.first)then
               if(.not.allocated(idx_a2i)) 
     *            allocate(idx_a2i(Bst%imt,Bst%jmt,ncv,2),
     *            wt_a2i(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Atmosphere --> Ice'
               Open(183,file='a2i.con',form='binary',status='old',err=3)
               Read(183)imti,jmti,imtf,jmtf,ncvf
               if(.not.(Ast%imt.eq.imti.and.Ast%jmt.eq.jmti.and.
     *            Bst%imt.eq.imtf.and.Bst%jmt.eq.jmtf.and.ncv.eq.ncvf))then
                  write(*,*)'New grid convertion coeffs are needed'
                  Close(183)
                  goto 3
               endif
               Read(183)idx_a2i,wt_a2i
               Close(183)
               goto 4
    3          call coeffs(idx_a2i,wt_a2i,Agr,Ast,Bgr,Bst)
               Open(183,file='a2i.con',form='binary',status='unknown')
               Write(183)Ast%imt,Ast%jmt,Bst%imt,Bst%jmt,ncv,idx_a2i,wt_a2i
               Close(183)
    4          continue
            endif
            call cnv_scalar(idx_a2i, wt_a2i, A, Ast, B, Bst, Bgr)
         case(a2l)
            ncv = ncv_a2l
            cnv_type = a2l_type
            if(.not.allocated(idx_a2l).or.first)then
               if(.not.allocated(idx_a2l)) 
     *            allocate(idx_a2l(Bst%imt,Bst%jmt,ncv,2),
     *            wt_a2l(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Atmosphere --> Land'
               Open(183,file='a2l.con',form='binary',status='old',err=5)
               Read(183)imti,jmti,imtf,jmtf,ncvf
               if(.not.(Ast%imt.eq.imti.and.Ast%jmt.eq.jmti.and.
     *            Bst%imt.eq.imtf.and.Bst%jmt.eq.jmtf.and.ncv.eq.ncvf))then
                  write(*,*)'New grid convertion coeffs are needed'
                  Close(183)
                  goto 5
               endif
               Read(183)idx_a2l,wt_a2l
               Close(183)
               goto 5
    5          call coeffs(idx_a2l,wt_a2l,Agr,Ast,Bgr,Bst)
               Open(183,file='a2l.con',form='binary',status='unknown')
               Write(183)Ast%imt,Ast%jmt,Bst%imt,Bst%jmt,ncv,idx_a2l,wt_a2l
               Close(183)
    6          continue
            endif
            call cnv_scalar(idx_a2l, wt_a2l, A, Ast, B, Bst, Bgr)
         case(o2a)
            ncv = ncv_o2a
            cnv_type = o2a_type
            if(.not.allocated(idx_o2a).or.first)then
               if(.not.allocated(idx_o2a)) 
     *            allocate(idx_o2a(Bst%imt,Bst%jmt,ncv,2),
     *            wt_o2a(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Ocean --> Atmosphere'
               Open(183,file='o2a.con',form='binary',status='old',err=7)
               Read(183)imti,jmti,imtf,jmtf,ncvf
               if(.not.(Ast%imt.eq.imti.and.Ast%jmt.eq.jmti.and.
     *            Bst%imt.eq.imtf.and.Bst%jmt.eq.jmtf.and.ncv.eq.ncvf))then
                  write(*,*)'New grid convertion coeffs are needed'
                  Close(183)
                  goto 7
               endif
               Read(183)idx_o2a,wt_o2a
               Close(183)
               goto 8
    7          call coeffs(idx_o2a,wt_o2a,Agr,Ast,Bgr,Bst)
               Open(183,file='o2a.con',form='binary',status='unknown')
               Write(183)Ast%imt,Ast%jmt,Bst%imt,Bst%jmt,ncv,idx_o2a,wt_o2a
               Close(183)
    8          continue
            endif
            call cnv_scalar(idx_o2a, wt_o2a, A, Ast, B, Bst, Bgr)
         case(o2i)
            ncv = ncv_o2i
            cnv_type = o2i_type
            if(.not.allocated(idx_o2i).or.first)then
               if(.not.allocated(idx_o2i)) 
     *            allocate(idx_o2i(Bst%imt,Bst%jmt,ncv,2),
     *            wt_o2i(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Ocean --> Ice'
               Open(183,file='o2i.con',form='binary',status='old',err=9)
               Read(183)imti,jmti,imtf,jmtf,ncvf
               if(.not.(Ast%imt.eq.imti.and.Ast%jmt.eq.jmti.and.
     *            Bst%imt.eq.imtf.and.Bst%jmt.eq.jmtf.and.ncv.eq.ncvf))then
                  write(*,*)'New grid convertion coeffs are needed'
                  Close(183)
                  goto 9
               endif
               Read(183)idx_o2i,wt_o2i
               Close(183)
               goto 10
    9          call coeffs(idx_o2i,wt_o2i,Agr,Ast,Bgr,Bst)
               Open(183,file='o2i.con',form='binary',status='unknown')
               Write(183)Ast%imt,Ast%jmt,Bst%imt,Bst%jmt,ncv,idx_o2i,wt_o2i
               Close(183)
   10          continue
            endif
            call cnv_scalar(idx_o2i, wt_o2i, A, Ast, B, Bst, Bgr)
         case(o2l)
            ncv = ncv_o2l
            cnv_type = o2l_type
            if(.not.allocated(idx_o2l).or.first)then
               if(.not.allocated(idx_o2l)) 
     *            allocate(idx_o2l(Bst%imt,Bst%jmt,ncv,2),
     *            wt_o2l(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Ocean --> Land'
               Open(183,file='o2l.con',form='binary',status='old',err=11)
               Read(183)imti,jmti,imtf,jmtf,ncvf
               if(.not.(Ast%imt.eq.imti.and.Ast%jmt.eq.jmti.and.
     *            Bst%imt.eq.imtf.and.Bst%jmt.eq.jmtf.and.ncv.eq.ncvf))then
                  write(*,*)'New grid convertion coeffs are needed'
                  Close(183)
                  goto 11
               endif
               Read(183)idx_o2l,wt_o2l
               Close(183)
               goto 12
   11          call coeffs(idx_o2l,wt_o2l,Agr,Ast,Bgr,Bst)
               Open(183,file='o2l.con',form='binary',status='unknown')
               Write(183)Ast%imt,Ast%jmt,Bst%imt,Bst%jmt,ncv,idx_o2l,wt_o2l
               Close(183)
   12          continue
            endif
            call cnv_scalar(idx_o2l, wt_o2l, A, Ast, B, Bst, Bgr)
         case(i2a)
            ncv = ncv_i2a
            cnv_type = i2a_type
            if(.not.allocated(idx_i2a).or.first)then
               if(.not.allocated(idx_i2a)) 
     *            allocate(idx_i2a(Bst%imt,Bst%jmt,ncv,2),
     *            wt_i2a(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Ice --> Atmosphere'
               Open(183,file='i2a.con',form='binary',status='old',err=13)
               Read(183)imti,jmti,imtf,jmtf,ncvf
               if(.not.(Ast%imt.eq.imti.and.Ast%jmt.eq.jmti.and.
     *            Bst%imt.eq.imtf.and.Bst%jmt.eq.jmtf.and.ncv.eq.ncvf))then
                  write(*,*)'New grid convertion coeffs are needed'
                  Close(183)
                  goto 13
               endif
               Read(183)idx_i2a,wt_i2a
               Close(183)
               goto 14
   13          call coeffs(idx_i2a,wt_i2a,Agr,Ast,Bgr,Bst)
               Open(183,file='i2a.con',form='binary',status='unknown')
               Write(183)Ast%imt,Ast%jmt,Bst%imt,Bst%jmt,ncv,idx_i2a,wt_i2a
               Close(183)
   14          continue
            endif
            call cnv_scalar(idx_i2a, wt_i2a, A, Ast, B, Bst, Bgr)
         case(i2o)
            ncv = ncv_i2o
            cnv_type = i2o_type
            if(.not.allocated(idx_i2o).or.first)then
               if(.not.allocated(idx_i2o)) 
     *            allocate(idx_i2o(Bst%imt,Bst%jmt,ncv,2),
     *            wt_i2o(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Ice --> Ocean'
               Open(183,file='i2o.con',form='binary',status='old',err=15)
               Read(183)imti,jmti,imtf,jmtf,ncvf
               if(.not.(Ast%imt.eq.imti.and.Ast%jmt.eq.jmti.and.
     *            Bst%imt.eq.imtf.and.Bst%jmt.eq.jmtf.and.ncv.eq.ncvf))then
                  write(*,*)'New grid convertion coeffs are needed'
                  Close(183)
                  goto 15
               endif
               Read(183)idx_i2o,wt_i2o
               Close(183)
               goto 16
   15          call coeffs(idx_i2o,wt_i2o,Agr,Ast,Bgr,Bst)
               Open(183,file='i2o.con',form='binary',status='unknown')
               Write(183)Ast%imt,Ast%jmt,Bst%imt,Bst%jmt,ncv,idx_i2o,wt_i2o
               Close(183)
   16          continue
            endif
            call cnv_scalar(idx_i2o, wt_i2o, A, Ast, B, Bst, Bgr)
         case(i2l)
            ncv = ncv_i2l
            cnv_type = i2l_type
            if(.not.allocated(idx_i2l).or.first)then
               if(.not.allocated(idx_i2l)) 
     *            allocate(idx_i2l(Bst%imt,Bst%jmt,ncv,2),
     *            wt_i2l(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Ice --> Land'
               Open(183,file='i2l.con',form='binary',status='old',err=17)
               Read(183)imti,jmti,imtf,jmtf,ncvf
               if(.not.(Ast%imt.eq.imti.and.Ast%jmt.eq.jmti.and.
     *            Bst%imt.eq.imtf.and.Bst%jmt.eq.jmtf.and.ncv.eq.ncvf))then
                  write(*,*)'New grid convertion coeffs are needed'
                  Close(183)
                  goto 17
               endif
               Read(183)idx_i2l,wt_i2l
               Close(183)
               goto 18
   17          call coeffs(idx_i2l,wt_i2l,Agr,Ast,Bgr,Bst)
               Open(183,file='i2l.con',form='binary',status='unknown')
               Write(183)Ast%imt,Ast%jmt,Bst%imt,Bst%jmt,ncv,idx_i2l,wt_i2l
               Close(183)
   18          continue
            endif
            call cnv_scalar(idx_i2l, wt_i2l, A, Ast, B, Bst, Bgr)
         case(l2a)
            ncv = ncv_l2a
            cnv_type = l2a_type
            if(.not.allocated(idx_l2a).or.first)then
               if(.not.allocated(idx_l2a)) 
     *            allocate(idx_l2a(Bst%imt,Bst%jmt,ncv,2),
     *            wt_l2a(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Land --> Atmosphere'
               Open(183,file='l2a.con',form='binary',status='old',err=19)
               Read(183)imti,jmti,imtf,jmtf,ncvf
               if(.not.(Ast%imt.eq.imti.and.Ast%jmt.eq.jmti.and.
     *            Bst%imt.eq.imtf.and.Bst%jmt.eq.jmtf.and.ncv.eq.ncvf))then
                  write(*,*)'New grid convertion coeffs are needed'
                  Close(183)
                  goto 19
               endif
               Read(183)idx_l2a,wt_l2a
               Close(183)
               goto 20
   19          call coeffs(idx_l2a,wt_l2a,Agr,Ast,Bgr,Bst)
               Open(183,file='l2a.con',form='binary',status='unknown')
               Write(183)Ast%imt,Ast%jmt,Bst%imt,Bst%jmt,ncv,idx_l2a,wt_l2a
               Close(183)
   20          continue
            endif
            call cnv_scalar(idx_l2a, wt_l2a, A, Ast, B, Bst, Bgr)
         case(l2o)
            ncv = ncv_l2o
            cnv_type = l2o_type
            if(.not.allocated(idx_l2o).or.first)then
               if(.not.allocated(idx_l2o)) 
     *            allocate(idx_l2o(Bst%imt,Bst%jmt,ncv,2),
     *            wt_l2o(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Land --> Ocean'
               Open(183,file='l2o.con',form='binary',status='old',err=21)
               Read(183)imti,jmti,imtf,jmtf,ncvf
               if(.not.(Ast%imt.eq.imti.and.Ast%jmt.eq.jmti.and.
     *            Bst%imt.eq.imtf.and.Bst%jmt.eq.jmtf.and.ncv.eq.ncvf))then
                  write(*,*)'New grid convertion coeffs are needed'
                  Close(183)
                  goto 21
               endif
               Read(183)idx_l2o,wt_l2o
               Close(183)
               goto 22
   21          call coeffs(idx_l2o,wt_l2o,Agr,Ast,Bgr,Bst)
               Open(183,file='l2o.con',form='binary',status='unknown')
               Write(183)Ast%imt,Ast%jmt,Bst%imt,Bst%jmt,ncv,idx_l2o,wt_l2o
               Close(183)
   22          continue
            endif
            call cnv_scalar(idx_l2o, wt_l2o, A, Ast, B, Bst, Bgr)
         case(l2i)
            ncv = ncv_l2i
            cnv_type = l2i_type
            if(.not.allocated(idx_l2i).or.first)then
               if(.not.allocated(idx_l2i)) 
     *            allocate(idx_l2i(Bst%imt,Bst%jmt,ncv,2),
     *            wt_l2i(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Land --> Ice'
               Open(183,file='l2i.con',form='binary',status='old',err=23)
               Read(183)imti,jmti,imtf,jmtf,ncvf
               if(.not.(Ast%imt.eq.imti.and.Ast%jmt.eq.jmti.and.
     *            Bst%imt.eq.imtf.and.Bst%jmt.eq.jmtf.and.ncv.eq.ncvf))then
                  write(*,*)'New grid convertion coeffs are needed'
                  Close(183)
                  goto 23
               endif
               Read(183)idx_l2i,wt_l2i
               Close(183)
               goto 24
   23          call coeffs(idx_l2i,wt_l2i,Agr,Ast,Bgr,Bst)
               Open(183,file='l2i.con',form='binary',status='unknown')
               Write(183)Ast%imt,Ast%jmt,Bst%imt,Bst%jmt,ncv,idx_l2i,wt_l2i
               Close(183)
   24          continue
            endif
            call cnv_scalar(idx_l2i, wt_l2i, A, Ast, B, Bst, Bgr)
      end select

      Return
      End subroutine convert_scalar

      Subroutine convert_vector(typ, Ax, Ay, Agr, Ast, 
     *                               Bx, By, Bgr, Bst, Reset)
! function [Bx,By]=convert_vector(Ax,Ay,Agr,Ast,Bgr,Bst)
! global idx0 wt0

      Type (state) Ast, Bst
      Type (grid) Agr, Bgr
      real*8 Ax(Ast%imt,Ast%jmt), Bx(Bst%imt,Bst%jmt)
      real*8 Ay(Ast%imt,Ast%jmt), By(Bst%imt,Bst%jmt)
      logical, optional :: Reset
      logical first
      integer nu_diag, typ
      common/iodiag/nu_diag

      if (present(Reset)) then
         first=Reset
      else
         first=.false.
      endif

      if(ncv.eq.0)ncv=4
      if(cnv_type.eq.0)cnv_type=cnv_exp

      select case(typ)
         case(a2ov)
            ncv = ncv_a2o
            cnv_type = a2o_type
            if(.not.allocated(idx_a2o).or.first)then
               if(.not.allocated(idx_a2o)) 
     *            allocate(idx_a2o(Bst%imt,Bst%jmt,ncv,2),
     *                   wt_a2o(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Atmosphere --> Ocean'
               call coeffs(idx_a2o,wt_a2o,Agr,Ast,Bgr,Bst)
! if isempty(idx0) | isempty(wt0),
!    [idx0,wt0]=coeffs(Agr,Ast,Bgr,Bst);
! end
! idx=idx0; wt=wt0;
            endif
            if(.not.allocated(ang_a2o))then
               allocate(ang_a2o(Bst%imt,Bst%jmt))
               call angulate(a2o,ang_a2o,Agr,Ast,Bgr,Bst)
! ang=angulate(Agr,Ast,Bgr,Bst);
            endif
            call cnv_vector(idx_a2o, wt_a2o, ang_a2o, 
     *         Ax, Ay, Ast, Agr, Bx, By, Bst, Bgr)
! [Bx,By]=cnv_vector(idx, wt, ang, Ax, Ay, Ast, Bst, Bgr);
         case(a2iv)
            ncv = ncv_a2i
            cnv_type = a2i_type
            if(.not.allocated(idx_a2i).or.first)then
               if(.not.allocated(idx_a2i)) 
     *            allocate(idx_a2i(Bst%imt,Bst%jmt,ncv,2),
     *                   wt_a2i(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Atmosphere --> Ice'
               call coeffs(idx_a2i,wt_a2i,Agr,Ast,Bgr,Bst)
            endif
            if(.not.allocated(ang_a2i))then
               allocate(ang_a2i(Bst%imt,Bst%jmt))
               call angulate(a2i,ang_a2i,Agr,Ast,Bgr,Bst)
            endif
            call cnv_vector(idx_a2i, wt_a2i, ang_a2i, 
     *         Ax, Ay, Ast, Agr, Bx, By, Bst, Bgr)
         case(a2lv)
            ncv = ncv_a2l
            cnv_type = a2l_type
            if(.not.allocated(idx_a2l).or.first)then
               if(.not.allocated(idx_a2l)) 
     *            allocate(idx_a2l(Bst%imt,Bst%jmt,ncv,2),
     *                   wt_a2l(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Atmosphere --> Land'
               call coeffs(idx_a2l,wt_a2l,Agr,Ast,Bgr,Bst)
            endif
            if(.not.allocated(ang_a2l))then
               allocate(ang_a2l(Bst%imt,Bst%jmt))
               call angulate(a2l,ang_a2l,Agr,Ast,Bgr,Bst)
            endif
            call cnv_vector(idx_a2l, wt_a2l, ang_a2l, 
     *         Ax, Ay, Ast, Agr, Bx, By, Bst, Bgr)
         case(o2av)
            ncv = ncv_o2a
            cnv_type = o2a_type
            if(.not.allocated(idx_o2a).or.first)then
               if(.not.allocated(idx_o2a)) 
     *            allocate(idx_o2a(Bst%imt,Bst%jmt,ncv,2),
     *                   wt_o2a(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Ocean --> Atmosphere'
               call coeffs(idx_o2a,wt_o2a,Agr,Ast,Bgr,Bst)
            endif
            if(.not.allocated(ang_o2a))then
               allocate(ang_o2a(Bst%imt,Bst%jmt))
               call angulate(o2a,ang_o2a,Agr,Ast,Bgr,Bst)
            endif
            call cnv_vector(idx_o2a, wt_o2a, ang_o2a,
     *         Ax, Ay, Ast, Agr, Bx, By, Bst, Bgr)
         case(o2iv)
            ncv = ncv_o2i
            cnv_type = o2i_type
            if(.not.allocated(idx_o2i).or.first)then
               if(.not.allocated(idx_o2i)) 
     *            allocate(idx_o2i(Bst%imt,Bst%jmt,ncv,2),
     *                   wt_o2i(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Ocean --> Ice'
               call coeffs(idx_o2i,wt_o2i,Agr,Ast,Bgr,Bst)
            endif
            if(.not.allocated(ang_o2i))then
               allocate(ang_o2i(Bst%imt,Bst%jmt))
               call angulate(o2i,ang_o2i,Agr,Ast,Bgr,Bst)
            endif
            call cnv_vector(idx_o2i, wt_o2i, ang_o2i, 
     *         Ax, Ay, Ast, Agr, Bx, By, Bst, Bgr)
         case(o2lv)
            ncv = ncv_o2l
            cnv_type = o2l_type
            if(.not.allocated(idx_o2l).or.first)then
               if(.not.allocated(idx_o2l)) 
     *            allocate(idx_o2l(Bst%imt,Bst%jmt,ncv,2),
     *                   wt_o2l(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Ocean --> Land'
               call coeffs(idx_o2l,wt_o2l,Agr,Ast,Bgr,Bst)
            endif
            if(.not.allocated(ang_o2l))then
               allocate(ang_o2l(Bst%imt,Bst%jmt))
               call angulate(o2l,ang_o2l,Agr,Ast,Bgr,Bst)
            endif
            call cnv_vector(idx_o2l, wt_o2l, ang_o2l, 
     *         Ax, Ay, Ast, Agr, Bx, By, Bst, Bgr)
         case(i2av)
            ncv = ncv_i2a
            cnv_type = i2a_type
            if(.not.allocated(idx_i2a).or.first)then
               if(.not.allocated(idx_i2a)) 
     *            allocate(idx_i2a(Bst%imt,Bst%jmt,ncv,2),
     *                   wt_i2a(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Ice --> Atmosphere'
               call coeffs(idx_i2a,wt_i2a,Agr,Ast,Bgr,Bst)
            endif
            if(.not.allocated(ang_i2a))then
               allocate(ang_i2a(Bst%imt,Bst%jmt))
               call angulate(i2a,ang_i2a,Agr,Ast,Bgr,Bst)
            endif
            call cnv_vector(idx_i2a, wt_i2a, ang_i2a,
     *         Ax, Ay, Ast, Agr, Bx, By, Bst, Bgr)
         case(i2ov)
            ncv = ncv_i2o
            cnv_type = i2o_type
            if(.not.allocated(idx_i2o).or.first)then
               if(.not.allocated(idx_i2o)) 
     *            allocate(idx_i2o(Bst%imt,Bst%jmt,ncv,2),
     *                   wt_i2o(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Ice --> Ocean'
               call coeffs(idx_i2o,wt_i2o,Agr,Ast,Bgr,Bst)
            endif
            if(.not.allocated(ang_i2o))then
               allocate(ang_i2o(Bst%imt,Bst%jmt))
               call angulate(i2o,ang_i2o,Agr,Ast,Bgr,Bst)
            endif
            call cnv_vector(idx_i2o, wt_i2o, ang_i2o,
     *         Ax, Ay, Ast, Agr, Bx, By, Bst, Bgr)
         case(i2lv)
            ncv = ncv_i2l
            cnv_type = i2l_type
            if(.not.allocated(idx_i2l).or.first)then
               if(.not.allocated(idx_i2l)) 
     *            allocate(idx_i2l(Bst%imt,Bst%jmt,ncv,2),
     *                   wt_i2l(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Ice --> Land'
               call coeffs(idx_i2l,wt_i2l,Agr,Ast,Bgr,Bst)
            endif
            if(.not.allocated(ang_i2l))then
               allocate(ang_i2l(Bst%imt,Bst%jmt))
               call angulate(i2l,ang_i2l,Agr,Ast,Bgr,Bst)
            endif
            call cnv_vector(idx_i2l, wt_i2l, ang_i2l,
     *         Ax, Ay, Ast, Agr, Bx, By, Bst, Bgr)
         case(l2av)
            ncv = ncv_l2a
            cnv_type = l2a_type
            if(.not.allocated(idx_l2a).or.first)then
               if(.not.allocated(idx_l2a)) 
     *            allocate(idx_l2a(Bst%imt,Bst%jmt,ncv,2),
     *                   wt_l2a(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Land --> Atmosphere'
               call coeffs(idx_l2a,wt_l2a,Agr,Ast,Bgr,Bst)
            endif
            if(.not.allocated(ang_l2a))then
               allocate(ang_l2a(Bst%imt,Bst%jmt))
               call angulate(l2a,ang_l2a,Agr,Ast,Bgr,Bst)
            endif
            call cnv_vector(idx_l2a, wt_l2a, ang_l2a,
     *         Ax, Ay, Ast, Agr, Bx, By, Bst, Bgr)
         case(l2ov)
            ncv = ncv_l2o
            cnv_type = l2o_type
            if(.not.allocated(idx_l2o).or.first)then
               if(.not.allocated(idx_l2o)) 
     *            allocate(idx_l2o(Bst%imt,Bst%jmt,ncv,2),
     *                   wt_l2o(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Land --> Ocean'
               call coeffs(idx_l2o,wt_l2o,Agr,Ast,Bgr,Bst)
            endif
            if(.not.allocated(ang_l2o))then
               allocate(ang_l2o(Bst%imt,Bst%jmt))
               call angulate(l2o,ang_l2o,Agr,Ast,Bgr,Bst)
            endif
            call cnv_vector(idx_l2o, wt_l2o, ang_l2o,
     *         Ax, Ay, Ast, Agr, Bx, By, Bst, Bgr)
         case(l2iv)
            ncv = ncv_l2i
            cnv_type = l2i_type
            if(.not.allocated(idx_l2i).or.first)then
               if(.not.allocated(idx_l2i)) 
     *            allocate(idx_l2i(Bst%imt,Bst%jmt,ncv,2),
     *                   wt_l2i(Bst%imt,Bst%jmt,ncv))
               write(*,*)'Land --> Ice'
               call coeffs(idx_l2i,wt_l2i,Agr,Ast,Bgr,Bst)
            endif
            if(.not.allocated(ang_l2i))then
               allocate(ang_l2i(Bst%imt,Bst%jmt))
               call angulate(l2i,ang_l2i,Agr,Ast,Bgr,Bst)
            endif
            call cnv_vector(idx_l2i, wt_l2i, ang_l2i,
     *         Ax, Ay, Ast, Agr, Bx, By, Bst, Bgr)
      end select

      Return
      End subroutine convert_vector

      Subroutine coeffs(idx,wt,Agr,Ast,Bgr,Bst)
! function [idx,wt]=coeffs(Agr,Ast,Bgr,Bst)
! global ncv yno
      Type (state) Ast, Bst
      Type (grid) Agr, Bgr
      integer*2 idx(Bst%imt,Bst%jmt,ncv,2)
      real*8    wt(Bst%imt,Bst%jmt,ncv)

      real*8, allocatable :: w(:,:)
      real*8 pi0,earth,r,rad,hg
      integer ind(2)

      write(*,*)'coeffs is started'
      allocate(w(Ast%imt,Ast%jmt))
! w=zeros(Ast.imt,Ast.jmt);

      pi0=4.0*datan2(1.0,1.0)/180.0
! pi0=pi/180.0;
      earth=6.371d+06
! earth=6.371e+06;
      hg=huge(hg)
! hg=realmax;
      radmin=sqrt(minval(Agr%tarea,MASK=Agr%mask.NE.0))
      radmin=max(radmin,sqrt(minval(Bgr%tarea,MASK=Bgr%mask.NE.0)))/pi0/earth

      Do i=1,Bst%imt
         if(100*int(i/100).eq.i) write(*,*)'    i=',i,'     of ',Bst%imt
         Do j=1,Bst%jmt
! for i=1:Bst.imt,
!     disp(['i=',int2str(i),' of ',int2str(Bst.imt)]);
!     for j=1:Bst.jmt,
            if(Bgr%mask(i,j).ne.0)then
!         if Bgr.mask(i,j)~=0,
               if(dabs(Bgr%lat(i,j)).le.yno)then
!             if abs(Bgr.lat(i,j))<=yno,
                  Do m=1,Ast%imt
                     Do n=1,Ast%jmt
!                 for m=1:Ast.imt,
!                     for n=1:Ast.jmt,
                        w(m,n)=hg
!                         w(m,n)=hg;
                        if(Agr%mask(m,n).ne.0)then
!                         if Agr.mask(m,n)~=0,
                           if(dabs(Agr%lat(m,n)).le.yno+5.0)then
!                             if abs(Agr.lat(m,n))<=yno+5,
                              xrad=sin(Agr%lon(m,n)*pi0)*cos(Agr%lat(m,n)*pi0)-
     *                             sin(Bgr%lon(i,j)*pi0)*cos(Bgr%lat(i,j)*pi0)
!                                 xrad=sin(Agr.lon(m,n)*pi0)*cos(Agr.lat(m,n)*pi0)- ...
!                                     sin(Bgr.lon(i,j)*pi0)*cos(Bgr.lat(i,j)*pi0);
                              yrad=cos(Bgr%lon(i,j)*pi0)*cos(Bgr%lat(i,j)*pi0)-
     *                             cos(Agr%lon(m,n)*pi0)*cos(Agr%lat(m,n)*pi0)
!                                 yrad=cos(Bgr.lon(i,j)*pi0)*cos(Bgr.lat(i,j)*pi0)- ...
!                                     cos(Agr.lon(m,n)*pi0)*cos(Agr.lat(m,n)*pi0);
                              zrad=sin(Agr%lat(m,n)*pi0)-sin(Bgr%lat(i,j)*pi0)
!                                 zrad=sin(Agr.lat(m,n)*pi0)-sin(Bgr.lat(i,j)*pi0);
                              w(m,n)=(xrad**2+yrad**2+zrad**2)/pi0**2
!                                 w(m,n)=(xrad^2+yrad^2+zrad^2)/pi0^2;
                           endif
                        endif
!                             end
!                         end
                     Enddo
                  Enddo
!                     end
!                 end
                  r=0.0
!                 r=0.0;
                  Do k=1,ncv
!                 for k=1:ncv,
                     ind=minloc(w)
                     m=ind(1)
                     n=ind(2)
!                     ind=find(w==min(w(:)));
!                     [m,n]=ind2sub(size(w),ind(1));
                     wt(i,j,k)=w(m,n)
!                     wt(i,j,k)=w(m,n);
                     r=r+sqrt(w(m,n))
!                     if k<=4, r=r+sqrt(w(m,n)); end
                     idx(i,j,k,1)=m
!                     idx(i,j,k,1)=m;
                     idx(i,j,k,2)=n
!                     idx(i,j,k,2)=n;
                     w(m,n)=hg
!                     w(m,n)=hg;
                  Enddo
!                 end
                  r=r/dble(ncv)
!                 r=r/2.0;
!                  rad=min(sqrt(2.0*Bgr%tarea(i,j))/pi0/earth,r)
!                  rad=max(rad,radmin)
                  rad=r
!                 rad=min([sqrt(2.0*Bgr.tarea(i,j))/pi0/earth,r]);
                  r=0.0
!                 r=0.0;
                  do k=1,ncv
!                 for k=1:ncv,
                     wt(i,j,k)=weight_fun(sqrt(wt(i,j,k)),rad)
!                     wt(i,j,k)=weight_fun(sqrt(wt(i,j,k)),rad);
                     if(wt(i,j,k).eq.hg)then
!                     if wt(i,j,k)==hg,
                        r=hg
!                         r=hg;
                     else
!                     else
                        r=r+wt(i,j,k)
!                         r=r+wt(i,j,k);
                     endif
!                     end
                  enddo
!                 end
                  do k=1,ncv
!                 for k=1:ncv,
                     if(r.eq.hg)then
!                     if r==hg,
                        if(wt(i,j,k).eq.hg)then
!                         if wt(i,j,k)==hg,
                           wt(i,j,k)=1.0
!                             wt(i,j,k)=1.0;
                        else
!                         else
                           wt(i,j,k)=0.0
!                             wt(i,j,k)=0.0;
                        endif
!                         end
                     else
!                     else
                        if(r.ne.0.0)then
                           wt(i,j,k)=wt(i,j,k)/r
                        else
                           wt(i,j,k)=1.0/dble(ncv)
                        endif
!                         if r~=0, wt(i,j,k)=wt(i,j,k)/r;
!                         else wt(i,j,k)=1/ncv; end
                     endif
!                     end
                  enddo
!                 end
               else
!             else
                  ref=90.0
!                 ref=90.0;
                  if(Bgr%lat(i,j).lt.0.0)ref=-90.0
!                 if Bgr.lat(i,j)<0.0, ref=-90.0; end
                  rr=dabs(ref-Bgr%lat(i,j))
!                 rr=abs(ref-Bgr.lat(i,j));
                  xi=rr*sin(Bgr%lon(i,j)*pi0)
!                 xi=rr*sin(Bgr.lon(i,j)*pi0);
                  yi=-rr*cos(Bgr%lon(i,j)*pi0)
!                 yi=-rr*cos(Bgr.lon(i,j)*pi0);
                  Do m=1,Ast%imt
                     Do n=1,Ast%jmt
!                 for m=1:Ast.imt,
!                     for n=1:Ast.jmt,
                        w(m,n)=hg
!                         w(m,n)=hg;
                        if(Agr%mask(m,n).ne.0)then
!                         if Agr.mask(m,n)~=0,
                           if(dabs(Agr%lat(m,n)).ge.yno-5.0)then
!                             if abs(Agr.lat(m,n))>=yno-5,
                              rr=dabs(ref-Agr%lat(m,n))
!                                 rr=abs(ref-Agr.lat(m,n));
                              xm=rr*sin(Agr%lon(m,n)*pi0)
!                                 xm=rr*sin(Agr.lon(m,n)*pi0);
                              ym=-rr*cos(Agr%lon(m,n)*pi0)
!                                 ym=-rr*cos(Agr.lon(m,n)*pi0);
                              w(m,n)=(xm-xi)**2+(ym-yi)**2
!                                 w(m,n)=(xm-xi)^2+(ym-yi)^2;
                           endif
!                             end
                        endif
!                         end
                     Enddo
!                     end
                  Enddo
!                 end
                  r=0.0
!                 r=0.0;
                  ilim=ncv/2
                  if(2*ilim.lt.ncv)ilim=ilim+1
                  jlim=ncv/2
                  if(2*jlim.lt.ncv)jlim=jlim+1
                  Do k=1,ncv
!                 for k=1:ncv,
                     ind=minloc(w)
                     m=ind(1)
                     n=ind(2)
                     if(k.eq.1)then
                         m0=m; n0=n
                     endif
!                     ind=find(w==min(w(:)));
!                     [m,n]=ind2sub(size(w),ind(1));
                     wt(i,j,k)=w(m,n)
!                     wt(i,j,k)=w(m,n);
                     r=r+sqrt(w(m,n))
!                     if k<=4, r=r+sqrt(w(m,n)); end
                     idx(i,j,k,1)=m
!                     idx(i,j,k,1)=m;
                     idx(i,j,k,2)=n
!                     idx(i,j,k,2)=n;
                     w(m,n)=hg
!                     w(m,n)=hg;
                     if(m.eq.m0)ilim=ilim-1
                     if(ilim.le.0) w(m0,:)=hg
                     if(n.eq.n0)jlim=jlim-1
                     if(jlim.le.0) w(:,n0)=hg
                  Enddo
!                 end
                  r=r/dble(ncv)
!                 r=r/2.0;
!                  rad=max(sqrt(2.0*Bgr%tarea(i,j))/pi0/earth,r)
!                  rad=max(rad,radmin)
                  rad=r
!                 rad=max([sqrt(2*Bgr.tarea(i,j))/pi0/earth,r]);
                  r=0.0
!                 r=0.0;
                  do k=1,ncv
!                 for k=1:ncv,
                     wt(i,j,k)=weight_fun(sqrt(wt(i,j,k)),rad)
!                     wt(i,j,k)=weight_fun(sqrt(wt(i,j,k)),rad);
                     if(wt(i,j,k).eq.hg)then
!                     if wt(i,j,k)==hg,
                        r=hg
!                         r=hg;
                     else
!                     else
                        r=r+wt(i,j,k)
!                         r=r+wt(i,j,k);
                     endif
!                     end
                  enddo
!                 end
                  do k=1,ncv
!                 for k=1:ncv,
                     if(r.eq.hg)then
!                     if r==hg,
                        if(wt(i,j,k).eq.hg)then
!                         if wt(i,j,k)==hg,
                           wt(i,j,k)=1.0
!                             wt(i,j,k)=1.0;
                        else
!                         else
                           wt(i,j,k)=0.0
!                             wt(i,j,k)=0.0;
                        endif
!                         end
                     else
!                     else
                        if(r.ne.0.0)then
                           wt(i,j,k)=wt(i,j,k)/r
                        else
                           wt(i,j,k)=1.0/dble(ncv)
                        endif
!                         if r~=0, wt(i,j,k)=wt(i,j,k)/r;
!                         else wt(i,j,k)=1/ncv; end
                     endif
!                     end
                  enddo
!                 end
               endif
!             end
            else
!         else
               do k=1,ncv
!             for k=1:ncv,
                  idx(i,j,k,1)=0
!                 idx(i,j,k,1)=0;
                  idx(i,j,k,2)=0
!                 idx(i,j,k,2)=0;
                  wt(i,j,k)=0.0
!                 wt(i,j,k)=0.0;
               enddo
!             end
            endif
!         end
         Enddo
!     end
      Enddo
! end

      deallocate(w)
      write(*,*)'coeffs is done'

      Return
      End subroutine coeffs

      real*8 function weight_fun(r,rad)
! function wf=weight_fun(r,rad)
! global cnv_type cnv_exp cnv_krig
      real*8 r, rad, xx

      xx=r/rad
! xx=r/rad;
      select case(cnv_type)
! switch cnv_type,
         case(cnv_exp)
!     case cnv_exp,
            weight_fun=dexp(-xx**2)
!         wf=exp(-xx^2);
         case(cnv_krig)
!     case cnv_krig,
            if(xx.eq.0.0)then
!         if xx==0.0,
               weight_fun=huge(xx)
!             wf=realmax;
            else
!         else
               weight_fun=1.0/xx**2
!             wf=1.0/xx^2;
            endif
!         end
      end select
! end

      end function weight_fun

      Subroutine angulate(typ,ang,Agr,Ast,Bgr,Bst)
! function ang=angulate(Agr,Ast,Bgr,Bst)
! global yno
      Type (state) Ast, Bst
      Type (grid) Agr, Bgr
      real*8    ang(Bst%imt,Bst%jmt)
      integer typ

      real*8, allocatable, dimension(:,:) ::
     *   dx,
     *   dy,
     *   dxm,
     *   dym,
     *   dxp,
     *   dyp
      integer, allocatable, dimension(:,:) ::
     *   mask,     ! source mask for interpolation
     *   mask0,    ! space to store original source mask
     *   mask1,    ! destination mask for interpolation
     *   mask01    ! space to store original destination mask

      logical postv

      real*8 pi0,earth,r,rad,hg

      pi0=4.0*datan2(1.0,1.0)/180.0
      earth=6.371d+06
      hg=huge(hg)
! pi0=pi/180.0;
! earth=6.371e+06;
! hg=realmax;

      allocate(dx(Bst%imt,Bst%jmt),dy(Bst%imt,Bst%jmt),
     *   dxp(Bst%imt,Bst%jmt),dyp(Bst%imt,Bst%jmt),
     *   dxm(Ast%imt,Ast%jmt),dym(Ast%imt,Ast%jmt),
     *   mask(Ast%imt,Ast%jmt),mask0(Ast%imt,Ast%jmt),
     *   mask1(Bst%imt,Bst%jmt),mask01(Bst%imt,Bst%jmt))

      mask=0
! mask=zeros(Ast.imt,Ast.jmt);
      mask1=0
! mask1=zeros(Bst.imt,Bst.jmt);
      mask0=Agr%mask
! mask0=Agr.mask;
      mask01=Bgr%mask
! mask01=Bgr.mask;
      dxm=0.0
! dxm=zeros(Ast.imt,Ast.jmt);
      dym=0.0
! dym=zeros(Ast.imt,Ast.jmt);

      Do m=1,Ast%imt
! for m=1:Ast.imt,
         mp1=min(m+1,Ast%imt)
         mm1=max(m-1,1)
!     mp1=min([m+1,Ast.imt]);
!     mm1=max([m-1,1]);
         Do n=1,Ast%jmt

!     for n=1:Ast.jmt,
            if(dabs(Agr%lat(m,n)).le. yno+5.0)then
!         if abs(Agr.lat(m,n))<=yno+5,
               xx1=sin(Agr%lon(mm1,n)*pi0)
!             xx1=sin(Agr.lon(mm1,n)*pi0);
               yy1=-cos(Agr%lon(mm1,n)*pi0)
!             yy1=-cos(Agr.lon(mm1,n)*pi0);
               xx2=sin(Agr%lon(mp1,n)*pi0)
!             xx2=sin(Agr.lon(mp1,n)*pi0);
               yy2=-cos(Agr%lon(mp1,n)*pi0)
!             yy2=-cos(Agr.lon(mp1,n)*pi0);
               dxm(m,n)=datan2(xx1*yy2-yy1*xx2,xx1*xx2+yy1*yy2)/pi0/dble(mp1-mm1)
!             dxm(m,n)=atan2(xx1*yy2-yy1*xx2,xx1*xx2+yy1*yy2)/pi0/(mp1-mm1);
               dym(m,n)=(Agr%lat(mp1,n)-Agr%lat(mm1,n))/dble(mp1-mm1)
!             dym(m,n)=(Agr.lat(mp1,n)-Agr.lat(mm1,n))/(mp1-mm1);
               mask(m,n)=1
!             mask(m,n)=1;
            endif
!         end
         enddo
!     end
      enddo
! end

!      Agr%mask=mask
! Agr.mask=mask;

      Do i=1,Bst%imt
! for i=1:Bst.imt,
         Do j=1,Bst%jmt
!     for j=1:Bst.jmt,
            if(dabs(Bgr%lat(i,j)).le. yno)then
!         if abs(Bgr.lat(i,j))<=yno,
               mask1(i,j)=1
!             mask1(i,j)=1;
            endif
!         end
         enddo
!     end
      enddo
! end

      Bgr%mask=mask1
! Bgr.mask=mask1;

      call convert_scalar(typ,dxm,Agr,Ast,dx,Bgr,Bst)
! dx=convert_scalar(dxm,Agr,Ast,Bgr,Bst);
      call convert_scalar(typ,dym,Agr,Ast,dy,Bgr,Bst)
! dy=convert_scalar(dym,Agr,Ast,Bgr,Bst);

      mask=0
! mask=zeros(Ast.imt,Ast.jmt);
      mask1=0
! mask1=zeros(Bst.imt,Bst.jmt);
      dxm=0.0
! dxm=zeros(Ast.imt,Ast.jmt);
      dym=0.0
! dym=zeros(Ast.imt,Ast.jmt);

      Do m=1,Ast%imt
! for m=1:Ast.imt,
         mp1=min(m+1,Ast%imt)
!     mp1=min([m+1,Ast.imt]);
         mm1=max(m-1,1)
!     mm1=max([m-1,1]);
         Do n=1,Ast%jmt
!     for n=1:Ast.jmt,
            if(Agr%lat(m,n).gt.yno-5.0)then
!         if Agr.lat(m,n)>yno-5,
               ref=90.0
!             ref=90;
               rr=ref-Agr%lat(mp1,n)
!             rr=ref-Agr.lat(mp1,n);
               xx2=rr*sin(Agr%lon(mp1,n)*pi0)
!             xx2=rr*sin(Agr.lon(mp1,n)*pi0);
               yy2=-rr*cos(Agr%lon(mp1,n)*pi0)
!             yy2=-rr*cos(Agr.lon(mp1,n)*pi0);

               rr=ref-Agr%lat(mm1,n)
!             rr=ref-Agr.lat(mm1,n);
               xx1=rr*sin(Agr%lon(mm1,n)*pi0)
!             xx1=rr*sin(Agr.lon(mm1,n)*pi0);
               yy1=-rr*cos(Agr%lon(mm1,n)*pi0)
!             yy1=-rr*cos(Agr.lon(mm1,n)*pi0);

               dxm(m,n)=(xx2-xx1)/dble(mp1-mm1)
!             dxm(m,n)=(xx2-xx1)/(mp1-mm1);
               dym(m,n)=(yy2-yy1)/dble(mp1-mm1)
!             dym(m,n)=(yy2-yy1)/(mp1-mm1);
               mask(m,n)=1
!             mask(m,n)=1;
            endif
!         end
         enddo
!     end
      enddo
! end

!      Agr%mask=mask
! Agr.mask=mask;

      Do i=1,Bst%imt
! for i=1:Bst.imt,
         Do j=1,Bst%jmt
!     for j=1:Bst.jmt,
            if(Bgr%lat(i,j).gt. yno)then
!         if Bgr.lat(i,j)>yno,
               mask1(i,j)=1
!             mask1(i,j)=1;
            endif
!         end
         enddo
!     end
      enddo
! end

      Bgr%mask=mask1
! Bgr.mask=mask1;

      call convert_scalar(typ,dxm,Agr,Ast,dxp,Bgr,Bst)
! dxp=convert_scalar(dxm,Agr,Ast,Bgr,Bst);
      call convert_scalar(typ,dym,Agr,Ast,dyp,Bgr,Bst)
! dyp=convert_scalar(dym,Agr,Ast,Bgr,Bst);

      Do i=1,Bst%imt
! for i=1:Bst.imt,
         Do j=1,Bst%jmt
!     for j=1:Bst.jmt,
            if(Bgr%mask(i,j).ne.0)then
c            dx(i,j)=dx(i,j)+dxp(i,j)
            dx(i,j)=dxp(i,j)
!         dx(i,j)=dx(i,j)+dxp(i,j);
c            dy(i,j)=dy(i,j)+dyp(i,j)
            dy(i,j)=dyp(i,j)
!         dy(i,j)=dy(i,j)+dyp(i,j);
            endif
         enddo
!     end
      enddo
! end

      mask=0
! mask=zeros(Ast.imt,Ast.jmt);
      mask1=0
! mask1=zeros(Bst.imt,Bst.jmt);
      dxm=0.0
! dxm=zeros(Ast.imt,Ast.jmt);
      dym=0.0
! dym=zeros(Ast.imt,Ast.jmt);

      Do m=1,Ast%imt
! for m=1:Ast.imt,
         mp1=min(m+1,Ast%imt)
!     mp1=min([m+1,Ast.imt]);
         mm1=max(m-1,1)
!     mm1=max([m-1,1]);
         Do n=1,Ast%jmt
!     for n=1:Ast.jmt,
            if(Agr%lat(m,n).lt.-yno+5.0)then
!         if Agr.lat(m,n)<-yno+5,
               ref=90.0
!             ref=90;
               rr=ref+Agr%lat(mp1,n)
!             rr=ref+Agr.lat(mp1,n);
               xx2=rr*sin(Agr%lon(mp1,n)*pi0)
!             xx2=rr*sin(Agr.lon(mp1,n)*pi0);
               yy2=-rr*cos(Agr%lon(mp1,n)*pi0)
!             yy2=-rr*cos(Agr.lon(mp1,n)*pi0);

               rr=ref+Agr%lat(mm1,n)
!             rr=ref+Agr.lat(mm1,n);
               xx1=rr*sin(Agr%lon(mm1,n)*pi0)
!             xx1=rr*sin(Agr.lon(mm1,n)*pi0);
               yy1=-rr*cos(Agr%lon(mm1,n)*pi0)
!             yy1=-rr*cos(Agr.lon(mm1,n)*pi0);

               dxm(m,n)=(xx2-xx1)/dble(mp1-mm1)
!             dxm(m,n)=(xx2-xx1)/(mp1-mm1);
               dym(m,n)=(yy2-yy1)/dble(mp1-mm1)
!             dym(m,n)=(yy2-yy1)/(mp1-mm1);
               mask(m,n)=1
!             mask(m,n)=1;
            endif
!         end
         enddo
!     end
      enddo
! end

!      Agr%mask=mask
! Agr.mask=mask;

      Do i=1,Bst%imt
! for i=1:Bst.imt,
         Do j=1,Bst%jmt
!     for j=1:Bst.jmt,
            if(Bgr%lat(i,j).lt.-yno)then
!         if Bgr.lat(i,j)<-yno,
               mask1(i,j)=1
!             mask1(i,j)=1;
            endif
!         end
         enddo
!     end
      enddo
! end

      Bgr%mask=mask1
! Bgr.mask=mask1;

      call convert_scalar(typ,dxm,Agr,Ast,dxp,Bgr,Bst)
! dxp=convert_scalar(dxm,Agr,Ast,Bgr,Bst);
      call convert_scalar(typ,dym,Agr,Ast,dyp,Bgr,Bst)
! dyp=convert_scalar(dym,Agr,Ast,Bgr,Bst);

      Do i=1,Bst%imt
! for i=1:Bst.imt,
         Do j=1,Bst%jmt
!     for j=1:Bst.jmt,
            dx(i,j)=dx(i,j)+dxp(i,j)
!         dx(i,j)=dx(i,j)+dxp(i,j);
            dy(i,j)=dy(i,j)+dyp(i,j)
!         dy(i,j)=dy(i,j)+dyp(i,j);
         enddo
!     end
      enddo
! end

      Agr%mask=mask0
! Agr.mask=mask0;
      Bgr%mask=mask01
! Bgr.mask=mask01;
      
      Do i=1,Bst%imt
! for i=1:Bst.imt,
         ip1=min(i+1,Bst%imt)
!     ip1=min([i+1,Bst.imt]);
         im1=max(i-1,1)
!     im1=max([i-1,1]);
         Do j=1,Bst%jmt
!     for j=1:Bst.jmt,
            if(dabs(Bgr%lat(i,j)).le. yno)then
!         if abs(Bgr.lat(i,j))<=yno,
               xx1=sin(Bgr%lon(im1,j)*pi0)
!             xx1=sin(Bgr.lon(im1,j)*pi0);
               yy1=-cos(Bgr%lon(im1,j)*pi0)
!             yy1=-cos(Bgr.lon(im1,j)*pi0);
               xx2=sin(Bgr%lon(ip1,j)*pi0)
!             xx2=sin(Bgr.lon(ip1,j)*pi0);
               yy2=-cos(Bgr%lon(ip1,j)*pi0)
!             yy2=-cos(Bgr.lon(ip1,j)*pi0);
               dxi=datan2(xx1*yy2-yy1*xx2,xx1*xx2+yy1*yy2)/pi0/dble(ip1-im1)
!             dxi=atan2(xx1*yy2-yy1*xx2,xx1*xx2+yy1*yy2)/pi0/(ip1-im1);
               dyi=(Bgr%lat(ip1,j)-Bgr%lat(im1,j))/dble(ip1-im1)
!             dyi=(Bgr.lat(ip1,j)-Bgr.lat(im1,j))/(ip1-im1);
            else
!         else
               ref=90.0
!             ref=90.0;
               if(Bgr%lat(i,j).lt. 0.0) ref=-90.0
!             if Bgr.lat(i,j)<0, ref=-90; end
               rr=dabs(ref-Bgr%lat(ip1,j))
!             rr=abs(ref-Bgr.lat(ip1,j));
               xx2=rr*sin(Bgr%lon(ip1,j)*pi0)
!             xx2=rr*sin(Bgr.lon(ip1,j)*pi0);
               yy2=-rr*cos(Bgr%lon(ip1,j)*pi0)
!             yy2=-rr*cos(Bgr.lon(ip1,j)*pi0);

               rr=dabs(ref-Bgr%lat(im1,j))
!             rr=abs(ref-Bgr.lat(im1,j));
               xx1=rr*sin(Bgr%lon(im1,j)*pi0)
!             xx1=rr*sin(Bgr.lon(im1,j)*pi0);
               yy1=-rr*cos(Bgr%lon(im1,j)*pi0)
!             yy1=-rr*cos(Bgr.lon(im1,j)*pi0);
               dxi=(xx2-xx1)/dble(ip1-im1)
!             dxi=(xx2-xx1)/(ip1-im1);
               dyi=(yy2-yy1)/dble(ip1-im1)
!             dyi=(yy2-yy1)/(ip1-im1);
            endif
!         end
            ang(i,j)=datan2(dxi*dy(i,j)-dyi*dx(i,j),dxi*dx(i,j)+dyi*dy(i,j))
!         ang(i,j)=atan2(dxi*dy(i,j)-dyi*dx(i,j),dxi*dx(i,j)+dyi*dy(i,j));
         enddo
!     end
      enddo
! end
      deallocate(dx,dy,dxp,dyp,dxm,dym,mask,mask0,mask1,mask01)

      return
      end subroutine angulate

      subroutine cnv_scalar(idx, wt, A, Ast, B, Bst,Bgr)
! function B=cnv_scalar(idx,wt,A,Ast,Bst,Bgr)
! global ncv
      Type (state) Ast, Bst
      Type (grid) Bgr
      integer*2   idx(Bst%imt,Bst%jmt,ncv,2)
      real*8    wt(Bst%imt,Bst%jmt,ncv)
      real*8    A(Ast%imt,Ast%jmt)
      real*8    B(Bst%imt,Bst%jmt)

      do i=1,Bst%imt
! for i=1:Bst.imt,
         do j=1,Bst%jmt
!     for j=1:Bst.jmt,
            B(i,j)=0.0
!         B(i,j)=0.0;
            if(Bgr%mask(i,j).ne.0)then
!         if Bgr.mask(i,j)~=0,
               do k=1,ncv
!             for k=1:ncv,
                  m=idx(i,j,k,1)
!                 m=idx(i,j,k,1);
                  n=idx(i,j,k,2)
!                 n=idx(i,j,k,2);
                  if(m*n.ne.0)then
!                 if m*n~=0,
                     B(i,j)=B(i,j)+A(m,n)*wt(i,j,k)
!                     B(i,j)=B(i,j)+A(m,n)*wt(i,j,k);
                  endif
!                 end
               enddo
!             end
            if(isnan(B(i,j)))B(i,j)=0.0
            endif
!         end
         enddo
!     end
      enddo
! end

      return
      end subroutine cnv_scalar

      subroutine cnv_vector(idx, wt, ang, Ax, Ay, Ast, Agr,
     *  Bx, By, Bst, Bgr)
! function [Bx,By]=cnv_vector(idx,wt,ang,Ax,Ay,Ast,Bst,Bgr)
! global ncv
      Type (state) Ast, Bst
      Type (grid) Agr, Bgr
      integer*2   idx(Bst%imt,Bst%jmt,ncv,2)
      real*8    wt(Bst%imt,Bst%jmt,ncv)
      real*8    ang(Bst%imt,Bst%jmt)
      real*8    Ax(Ast%imt,Ast%jmt)
      real*8    Ay(Ast%imt,Ast%jmt)
      real*8    Bx(Bst%imt,Bst%jmt)
      real*8    By(Bst%imt,Bst%jmt)
      real*8    lon0,dlon

      pi0=4.0*datan2(1.0,1.0)/180.0
! pi0=pi/180.0;
      do i=1,Bst%imt
! for i=1:Bst.imt,
         do j=1,Bst%jmt
!     for j=1:Bst.jmt,
            Bx(i,j)=0.0
            By(i,j)=0.0
!         Bx(i,j)=0.0;
!         By(i,j)=0.0;
            if(Bgr%mask(i,j).ne.0)then
            lon0=Bgr%lon(i,j)
            sn=sin(Bgr%lat(i,j)*pi0)
!         if Bgr.mask(i,j)~=0,
               do k=1,ncv
!             for k=1:ncv,
                  m=idx(i,j,k,1)
!                 m=idx(i,j,k,1);
                  n=idx(i,j,k,2)
!                 n=idx(i,j,k,2);
                  if(m*n.ne.0)then
!                 if m*n~=0,
                     dlon=(Agr%lon(m,n)-lon0)*pi0
                     Bx(i,j)=Bx(i,j)+(Ax(m,n)*cos(dlon)-Ay(m,n)*sin(dlon)*sn)*wt(i,j,k)
                     By(i,j)=By(i,j)+(Ax(m,n)*sin(dlon)*sn+Ay(m,n)*cos(dlon))*wt(i,j,k)
!                     Bx(i,j)=Bx(i,j)+Ax(m,n)*wt(i,j,k);
!                     By(i,j)=By(i,j)+Ay(m,n)*wt(i,j,k);
                  endif
!                 end
               enddo
            endif
!             end
!         end

            uu=Bx(i,j)*cos(ang(i,j))-By(i,j)*sin(ang(i,j))
!         uu=Bx(i,j)*cos(ang(i,j))-By(i,j)*sin(ang(i,j));
            vv=Bx(i,j)*sin(ang(i,j))+By(i,j)*cos(ang(i,j))
!         vv=Bx(i,j)*sin(ang(i,j))+By(i,j)*cos(ang(i,j));
            Bx(i,j)=uu
!         Bx(i,j)=uu;
            By(i,j)=vv
!         By(i,j)=vv;
            if(isnan(Bx(i,j)))Bx(i,j)=0.
            if(isnan(By(i,j)))By(i,j)=0.
         enddo
!     end
      enddo
! end
      
      return
      end subroutine cnv_vector

      end module conversion

