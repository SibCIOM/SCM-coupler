!===========================================================================================================
!     $Id: mpi_tools.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $
!===========================================================================================================
      subroutine setup_mpi(model)

      IMPLICIT none

      character*(*) model

      include "mpif.h"                  ! MPI library definitions

      integer ::
     &  cpl_task
     &, my_task           ! master task for ice
     &, MPI_COMM_MY
     &, my_cur_task
     &, nu_diag
      common/tasks/my_task, cpl_task, MPI_COMM_MY,my_cur_task
      common/iodiag/ nu_diag
      integer master_task, ierr

      master_task = 0

      call MPI_INIT(ierr)
      Write(nu_diag,*)'Starting model '//model
      call mpi_coupled (model, cpl_task, my_task, MPI_COMM_MY)
      
      call MPI_COMM_RANK (MPI_COMM_MY, my_cur_task, ierr)
      Write(nu_diag,*)'Model '//model//' was started: ',my_task,my_cur_task

      end subroutine setup_mpi

      subroutine mpi_coupled (in_model_name, 
     &                        cpl_task, model_task, model_comm)

      include "mpif.h"                  ! MPI library definitions

      character (3), intent(in) :: in_model_name   
                        ! 3-letter identifier (atm, lnd, ocn, ice, cpl)
                        ! for the model calling this routine

      integer, intent(out) ::
     &  cpl_task           ! master task of coupler
     &, model_task         ! master task of model (in MPI_COMM_WORLD)
     &, model_comm         ! communicator for internal model comms
      integer cpl_id,ice_id,ocn_id,atm_id,lnd_id,reg_id(20)
      common/allids/cpl_id, ice_id, ocn_id, atm_id, lnd_id, reg_id
      common/iodiag/ nu_diag

      character (3) :: cmodel   ! model name temporary

      integer, dimension(3) :: range  ! array for creating groups for
                                      !  each coupled model

      integer :: 
     &  MPI_GROUP_WORLD  ! group id for MPI_COMM_WORLD
     &, MPI_GROUP_CPL    ! group of processors assigned to cpl
     &, MPI_GROUP_OCN    ! group of processors assigned to ocn
     &, MPI_GROUP_ICE    ! group of processors assigned to ice
     &, MPI_GROUP_ATM    ! group of processors assigned to atm
     &, MPI_GROUP_LND    ! group of processors assigned to lnd
     &, MPI_GROUP_REG(20)! groups of processors assigned to regions
     &, MPI_COMM_REG(20) ! communicators for processors assigned to regions
     &, MPI_COMM_CPL     ! communicator for processors assigned to cpl
     &, MPI_COMM_OCN     ! communicator for processors assigned to ocn
     &, MPI_COMM_ICE     ! communicator for processors assigned to ice
     &, MPI_COMM_ATM     ! communicator for processors assigned to atm
     &, MPI_COMM_LND     ! communicator for processors assigned to lnd
     &, n                ! dummy loop counter
     &, ierr             ! error flag for MPI comms
     &, nprocs_cpl       ! total processor count
     &, my_task_coupled  ! rank of process in coupled domain
     &, cpl_rank_min, cpl_rank_max   ! processor range for each
     &, ocn_rank_min, ocn_rank_max   !  component model
     &, ice_rank_min, ice_rank_max
     &, atm_rank_min, atm_rank_max
     &, lnd_rank_min, lnd_rank_max
     &, reg_rank_min(20), reg_rank_max(20),
     &  regmap(20)
     

      !-----------------------------------------------------------------
      !     determine processor rank in coupled domain
      !-----------------------------------------------------------------

      call MPI_COMM_RANK (MPI_COMM_WORLD, my_task_coupled, ierr)

      !-----------------------------------------------------------------
      !     determine which group of processes assigned to each model
      !-----------------------------------------------------------------

      call MPI_COMM_SIZE (MPI_COMM_WORLD, nprocs_cpl, ierr)

      ocn_rank_min = nprocs_cpl
      ocn_rank_max = 0
      ice_rank_min = nprocs_cpl
      ice_rank_max = 0
      atm_rank_min = nprocs_cpl
      atm_rank_max = 0
      lnd_rank_min = nprocs_cpl
      lnd_rank_max = 0
      cpl_rank_min = nprocs_cpl
      cpl_rank_max = 0
      reg_rank_min = nprocs_cpl
      reg_rank_max = 0
      regmap=0

      !***
      !*** each processor broadcasts its model to all the processors
      !*** in the coupled domain
      !***

      do n=0,nprocs_cpl-1
        if (n == my_task_coupled) then
          cmodel = in_model_name
        else
          cmodel = 'unk'
        endif

        call MPI_BCAST(cmodel, 3, MPI_CHARACTER, n,
     &                            MPI_COMM_WORLD, ierr)

c        print *,'mpi_tool: cmodel= ''',cmodel,''''
        select case(cmodel)
        case ('ocn')
          ocn_rank_min = min(ocn_rank_min, n)
          ocn_rank_max = max(ocn_rank_max, n)
        case ('ice')
          ice_rank_min = min(ice_rank_min, n)
          ice_rank_max = max(ice_rank_max, n)
        case ('atm')
          atm_rank_min = min(atm_rank_min, n)
          atm_rank_max = max(atm_rank_max, n)
        case ('lnd')
          lnd_rank_min = min(lnd_rank_min, n)
          lnd_rank_max = max(lnd_rank_max, n)
        case ('cpl')
          cpl_rank_min = min(cpl_rank_min, n)
          cpl_rank_max = max(cpl_rank_max, n)
        case default
          if(cmodel(1:1).eq.'r')then
             read(cmodel(2:3),'(i2)')kreg
             reg_rank_min(kreg) = min(reg_rank_min(kreg), n)
             reg_rank_max(kreg) = max(reg_rank_max(kreg), n)
             regmap(kreg)=1
c             print *,'mpi_tool: kreg=',kreg,' rank=',n
          else
            write (nu_diag,*)'Unknown model ',cmodel,
     *        ' in coupled model domain'
            write (nu_diag,*)'Model must be cpl, ice, ocn, atm, rNN'
        call sys_abort ('(mpi_coupled) Unknown model in coupled domain')
          endif
        end select

      end do

      !-----------------------------------------------------------------
      !     create subroup and communicators for each models internal 
      !     communciations, note that MPI_COMM_CREATE must be called by
      !     all processes in MPI_COMM_WORLD so this must be done by all
      !     models consistently and in the same order.
      !-----------------------------------------------------------------

      call MPI_COMM_GROUP(MPI_COMM_WORLD, MPI_GROUP_WORLD, ierr)

      range(3) = 1

      range(1) = ocn_rank_min
      range(2) = ocn_rank_max
      call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range,
     &                          MPI_GROUP_OCN, ierr)

      range(1) = ice_rank_min
      range(2) = ice_rank_max
      call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range,
     &                          MPI_GROUP_ICE, ierr)

      range(1) = atm_rank_min
      range(2) = atm_rank_max
      call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range,
     &                          MPI_GROUP_ATM, ierr)

      range(1) = lnd_rank_min
      range(2) = lnd_rank_max
      call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range,
     &                          MPI_GROUP_LND, ierr)

      do kreg=1,20

         if(regmap(kreg).eq.1)then
            range(1) = reg_rank_min(kreg)
            range(2) = reg_rank_max(kreg)
      call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range,
     &                          MPI_GROUP_REG(kreg), ierr)
         endif

      enddo

      range(1) = cpl_rank_min
      range(2) = cpl_rank_max

      call MPI_GROUP_RANGE_INCL(MPI_GROUP_WORLD, 1, range,
     &                          MPI_GROUP_CPL, ierr)
      call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_OCN,
     &                      MPI_COMM_OCN, ierr)
      call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_ICE,
     &                      MPI_COMM_ICE, ierr)
      call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_ATM,
     &                      MPI_COMM_ATM, ierr)
      call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_LND,
     &                      MPI_COMM_LND, ierr)
      call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_CPL,
     &                      MPI_COMM_CPL, ierr)
      do kreg=1,20
         if(regmap(kreg).eq.1)then
      call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_GROUP_REG(kreg),
     &                      MPI_COMM_REG(kreg), ierr)
         reg_id(kreg) = reg_rank_min(kreg)
         endif

      enddo
 
      !-----------------------------------------------------------------
      !     determine coupler process and model processes
      !     assume the first processor in each domain is the task that 
      !     will communicate coupled model messages
      !-----------------------------------------------------------------

      cpl_task = cpl_rank_min
      cpl_id = cpl_rank_min
      ocn_id = ocn_rank_min
      ice_id = ice_rank_min
      atm_id = atm_rank_min
      lnd_id = lnd_rank_min

c        print *,'mpi_tool: in_model_name= ''',cmodel,''''

      select case(in_model_name)
      case ('ocn')
        model_task = ocn_rank_min
        model_comm = MPI_COMM_OCN
      case ('ice')
        model_task = ice_rank_min
        model_comm = MPI_COMM_ICE
      case ('atm')
        model_task = atm_rank_min
        model_comm = MPI_COMM_ATM
      case ('lnd')
        model_task = lnd_rank_min
        model_comm = MPI_COMM_LND
      case ('cpl')
        model_task = cpl_rank_min
        model_comm = MPI_COMM_CPL
      case default
        if(in_model_name(1:1).eq.'r')then
          read(in_model_name(2:3),'(i2)')kreg
          cpl_task = ocn_rank_min
          model_task = reg_rank_min(kreg)
          model_comm = MPI_COMM_REG(kreg)
c          print *,'mpi_tool: kreg=',kreg,' master=',model_task
        else
          write (nu_diag,*)'Unknown model ',in_model_name,
     *   ' in coupled model'
          write (nu_diag,*)'Model must be cpl, ice, ocn'
        call sys_abort ('(mpi_coupled) Unknown model in coupled domain')
        endif
      end select

      end subroutine mpi_coupled

      subroutine recv_i(cbuffi, ncbuffi, msgtype, task_id, ierr)

      include "mpif.h"         ! MPI library definitions

      integer (kind=4) ::
     &   ncbuffi            ! size of integer control buffer
     &,  ierr               ! error flag
      integer (kind=4) ::
     &   cbuffi(ncbuffi)    ! control buffer from cpl

      integer msgtype, task_id
      integer, dimension(MPI_STATUS_SIZE,2) :: 
     &   status             ! status array for communications
      integer nu_diag
      common/iodiag/ nu_diag

      call MPI_RECV(cbuffi, ncbuffi, MPI_INTEGER, task_id,
     &                       msgtype, MPI_COMM_WORLD, status, ierr)

      if (ierr /= MPI_SUCCESS ) then
        write (nu_diag,*) '(recv_i) ERROR after integer recv'
        call exit_coupler
        stop
      endif

c      write(nu_diag,*) '(recv_i) Received integer buffer from ',task_id

      return
      end
c--------------------------------
      subroutine recv_d(cbuffi, ncbuffi, msgtype, task_id, ierr)

      include "mpif.h"         ! MPI library definitions

      integer (kind=4) ::
     &   ncbuffi            ! size of integer control buffer
     &,  ierr               ! error flag
      real (kind=8) ::
     &   cbuffi(ncbuffi)    ! control buffer from cpl

      integer msgtype, task_id
      integer, dimension(MPI_STATUS_SIZE,2) :: 
     &   status             ! status array for communications
      integer nu_diag
      common/iodiag/ nu_diag

      call MPI_RECV(cbuffi, ncbuffi, MPI_DOUBLE_PRECISION, task_id,
     &                       msgtype, MPI_COMM_WORLD, status, ierr)

      if (ierr /= MPI_SUCCESS ) then
        write (nu_diag,*) '(recv_i) ERROR after real8 recv'
        call exit_coupler
        stop
      endif

c      write(nu_diag,*) '(recv_i) Received real8 buffer from ',task_id

      return
      end
c--------------------------------
      subroutine SEND_I(cbuffi, ncbuffi,msgtype, task_id, ierr)

      include "mpif.h"         ! MPI library definitions

      integer (kind=4) ::
     &   ncbuffi            ! size of integer control buffer
     &,  ierr               ! error flag
      integer (kind=4) ::
     &   cbuffi(ncbuffi)    ! control buffer from cpl

      integer msgtype, task_id
      integer, dimension(MPI_STATUS_SIZE,2) :: 
     &   status             ! status array for communications
      integer nu_diag
      common/iodiag/ nu_diag

      call MPI_SEND(cbuffi, ncbuffi, MPI_INTEGER, task_id,
     &                     msgtype, MPI_COMM_WORLD, ierr)

      if (ierr /= MPI_SUCCESS ) then
         write (nu_diag,*)'(send_i) ERROR after integer send'
         call exit_coupler
         stop
      endif

c      write(nu_diag,*) '(send_i) Sent integer buffer to ',task_id

      return
      end
c--------------------------------
      subroutine SEND_D(cbuffi, ncbuffi,msgtype, task_id, ierr)

      include "mpif.h"         ! MPI library definitions

      integer (kind=4) ::
     &   ncbuffi            ! size of integer control buffer
     &,  comm               ! so far unknown parameter
     &,  ierr               ! error flag
      real (kind=8) ::
     &   cbuffi(ncbuffi)    ! control buffer from cpl

      integer msgtype, task_id
      integer, dimension(MPI_STATUS_SIZE,2) :: 
     &   status             ! status array for communications
      integer nu_diag
      common/iodiag/ nu_diag

      call MPI_SEND(cbuffi, ncbuffi, MPI_DOUBLE_PRECISION, task_id,
     &                     msgtype, MPI_COMM_WORLD, ierr)

      if (ierr /= MPI_SUCCESS ) then
         write (nu_diag,*)'(send_d) ERROR after real8 send'
         call exit_coupler
         stop
      endif

c      write(nu_diag,*) '(send_d) Sent real8 buffer to ',task_id

      return
      end

      subroutine wait_mpi_partners

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      return
      end subroutine wait_mpi_partners

      subroutine exit_coupler
      include "mpif.h"         ! MPI library definitions

      integer ierr  ! error flag

      call MPI_FINALIZE(ierr)
      call sys_abort('Finished')
      
      Print *, 'Module is finished'
      stop

      end subroutine exit_coupler

      subroutine sys_abort(string)

c      use iflport

      IMPLICIT none

      include "mpif.h"                  ! MPI library definitions

      character*(*),optional :: string    ! error message string

      integer ::
     &  cpl_task
     &, my_task           ! master task for ice
     &, MPI_COMM_MY
     &, nu_diag
      common/tasks/my_task, cpl_task, MPI_COMM_MY
      common/iodiag/ nu_diag

   !----- local -----
      integer              :: rcode,ierr
      logical              :: flag

   !----- i/o formats -----
      character(len=*),parameter :: F00 = "('(sys_abort) ',a)"

!-------------------------------------------------------------------------------
! PURPOSE: consistent stopping mechanism
!-------------------------------------------------------------------------------

      call flush(nu_diag)
      if (len_trim(string) > 0) write(nu_diag,F00)
     *   'ERROR: '//trim(string)
      write(nu_diag,F00) 'WARNING: calling mpi_abort() and stopping'
      call flush(nu_diag)
      call mpi_initialized(flag,ierr)
      if (flag) call mpi_abort(MPI_COMM_WORLD,rcode,ierr)
      call flush(nu_diag)
      if(nu_diag.ne.6)Close(nu_diag)
      Print *,'Module is finished'
      stop

      end subroutine sys_abort
      
