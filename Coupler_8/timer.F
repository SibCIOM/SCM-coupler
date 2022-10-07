!     $Id: timer.F 7 2015-03-31 07:44:54Z plat $
!     $HeadURL: $

      module timer

      integer*8 timer1,timer0,oldtimer,timescale
      real*8 timers(50),ptimers(50),runtime
      integer*8 reftime,newtime
      data oldtimer/0/,timescale/0/
      
      contains

      subroutine clock_ini
      
         timers=0
         runtime=0.0
         reftime=0
         newtime=0
         call system_clock(count=reftime,count_rate=timescale)
         
      return
      end subroutine clock_ini

      subroutine clock_starttimer(n)
      
         call system_clock(count=timer1)
         if(oldtimer.ne.0)then
            if(timer0.gt.timer1)timer0=timer0-86400
            timers(oldtimer)=timers(oldtimer)+dble(timer1-timer0)
            runtime=runtime+dble(timer1-timer0)
            ptimers(oldtimer)=100.0*timers(oldtimer)/runtime
         endif
         timer0=timer1
         oldtimer=n
         
      return
      end subroutine clock_starttimer

      subroutine clock_start_time_table
      
      call system_clock(count=reftime)

      return
      end subroutine clock_start_time_table

      subroutine clock_print_time_table(preface)
      character*(*) preface
      integer nu_diag
      common/iodiag/nu_diag
      
      call system_clock(count=newtime)

c      if(timescale.ne.0)
c     *   Write(nu_diag,'(a,f7.2)')preface,dble(newtime-reftime)/dble(timescale)

      return
      end subroutine clock_print_time_table

      end module timer
