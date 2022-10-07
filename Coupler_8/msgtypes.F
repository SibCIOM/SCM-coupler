!     $Id: msgtypes.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $

      module msgtypes

      integer, parameter ::
     &   msgtype_c2i  = 30  ! message type for cpl->ice
     &,  msgtype_c2ii = 31  ! message type for cpl->ice init
     &,  msgtype_i2ci = 40  ! message type for ice->cpl init
     &,  msgtype_i2c  = 41  ! message type for ice->cpl
     &,  msgtype_c2o  = 50  ! message type for cpl->ocn
     &,  msgtype_c2oi = 51  ! message type for cpl->ocn init
     &,  msgtype_o2ci = 60  ! message type for ocn->cpl init
     &,  msgtype_o2c  = 61  ! message type for ocn->cpl
     &,  msgtype_c2a  = 70  ! message type for cpl->atm
     &,  msgtype_c2ai = 71  ! message type for cpl->atm init
     &,  msgtype_a2ci = 80  ! message type for atm->cpl init
     &,  msgtype_a2c  = 81  ! message type for atm->cpl
     &,  msgtype_c2l  = 75  ! message type for cpl->lnd
     &,  msgtype_c2li = 76  ! message type for cpl->lnd init
     &,  msgtype_l2ci = 85  ! message type for lnd->cpl init
     &,  msgtype_l2c  = 86  ! message type for lnd->cpl

      end module msgtypes