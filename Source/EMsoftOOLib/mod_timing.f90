! ###################################################################
! Copyright (c) 2013-2020, Marc De Graef Research Group/Carnegie Mellon University
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
!     - Redistributions of source code must retain the above copyright notice, this list 
!        of conditions and the following disclaimer.
!     - Redistributions in binary form must reproduce the above copyright notice, this 
!        list of conditions and the following disclaimer in the documentation and/or 
!        other materials provided with the distribution.
!     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names 
!        of its contributors may be used to endorse or promote products derived from 
!        this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
! USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ###################################################################

!--------------------------------------------------------------------------
! EMsoft:mod_timing.f90
!--------------------------------------------------------------------------
!
! MODULE: mod_timing
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Provides a few timing routines
! 
!> @date   11/19/01 MDG 1.0 original
!> @date   06/04/13 MDG 2.0 rewrite
!> @date   06/05/14 MDG 3.0 replaced globals by timetype argument
!--------------------------------------------------------------------------
module mod_timing

use mod_kinds
use mod_io
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit

IMPLICIT NONE

private
public :: T_TimingClass

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
  character ( len = 3 ), dimension(12) :: month = (/ 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
                                                     'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' /)

  type, public   ::  T_TimingClass
    private 
      integer(kind=irg)             :: numT
      integer(kind=irg),allocatable :: Tstart(:)
      integer(kind=irg),allocatable :: Tstop(:)
      integer(kind=irg),allocatable :: Tinterval(:)
      character(len = 11)           :: datestring
      character(len = 15)           :: timestring
      character(len = 26)           :: timestamp

    contains
      private 

        procedure, pass(self), public :: Time_tick
        procedure, pass(self), public :: Time_tock
        procedure, pass(self), public :: Time_reset
        procedure, pass(self), public :: getInterval
        procedure, pass(self), public :: getDateString
        procedure, pass(self), public :: getTimeString
        procedure, pass(self), public :: printTimeStamp
        procedure, pass(self)         :: makeTimeStamp
 
   end type T_TimingClass

! the constructor routine for this class 
   interface T_TimingClass
     module procedure :: Timing_constructor
   end interface T_TimingClass        


contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! we begin with the functions/subroutines that are public in this class
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! FUNCTION: Timing_constructor
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief initialize the Timing class with the current date and time strings, print if timestamp given
!
!> @date  01/01/20 MDG 1.0 new function
!--------------------------------------------------------------------------
type(T_TimingClass) function Timing_constructor( showDateTime, nCounters ) result(Timing)

use mod_io 

IMPLICIT NONE

  logical, INTENT(IN), OPTIONAL            :: showDateTime
  integer(kind=irg), INTENT(IN), OPTIONAL  :: nCounters

  type(T_IOClass)                          :: Message

! set the maximum number of time intervals to be stored
  if (present(nCounters)) then 
    Timing % numT = nCounters
  else
    Timing % numT = 10   ! 10 is the default number of timer slots
  end if 
  allocate( Timing % Tstart(Timing % numT), Timing % Tstop(Timing % numT), &
            Timing % Tinterval(Timing % numT) )

! initialize date and time strings
  call Timing % makeTimeStamp() 

! and print them if requested 
  if (present(showDateTime)) then
    if (showDateTime) then 
      Message = T_IOClass()
      call Message % printMessage( Timing % timestamp, frm="(A/)")
    end if 
  end if

end function Timing_constructor


!--------------------------------------------------------------------------
!
! SUBROUTINE: Time_tick
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief start time recording using system_clock
!
!> @date   03/17/16 MDG 1.0 original
!--------------------------------------------------------------------------
recursive subroutine Time_tick(self, n)

IMPLICIT NONE

class(T_TimingClass)                    :: self 
integer(kind=irg), intent(IN), OPTIONAL :: n

integer(kind=irg)                       :: i, t

i = 1
if (present(n)) i = n 

call system_clock(t)

self % Tstart(i) = t 

end subroutine Time_tick

!--------------------------------------------------------------------------
!
! SUBROUTINE: Time_tock
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief stop time recording using system_clock
!
!> @date   03/17/16 MDG 1.0 original
!--------------------------------------------------------------------------
recursive subroutine Time_tock(self, n) 

IMPLICIT NONE

class(T_TimingClass)                    :: self 
integer(kind=irg), intent(in), OPTIONAL :: n

integer(kind=irg)                       :: i, now, clock_rate

i = 1
if (present(n)) i = n 

call system_clock(now,clock_rate)
self % Tstop(i) = now 

self % Tinterval(i) = real(now - self % Tstart(i))/real(clock_rate)

end subroutine Time_tock

!--------------------------------------------------------------------------
!
! SUBROUTINE: Time_reset
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief reset one or all counters 
!
!> @date   01/01/20 MDG 1.0 original
!--------------------------------------------------------------------------
recursive subroutine Time_reset(self, n) 

IMPLICIT NONE

class(T_TimingClass)                    :: self 
integer(kind=irg), intent(in), OPTIONAL :: n

integer(kind=irg)                       :: i

if (present(n)) then 
  i = n 
  self % Tstart(i) = 0
  self % Tstop(i) = 0
  self % Tinterval(i) = 0
else
  self % Tstart = 0
  self % Tstop = 0
  self % Tinterval = 0
end if 

end subroutine Time_reset


!--------------------------------------------------------------------------
!
! FUNCTION: getInterval
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief stop time recording using system_clock
!
!> @date   01/01/20 MDG 1.0 original
!--------------------------------------------------------------------------
recursive function getInterval(self, n) result(t)

IMPLICIT NONE

class(T_TimingClass)                    :: self 
integer(kind=irg), intent(in), OPTIONAL :: n

real(kind=sgl)                          :: t

integer(kind=irg)                       :: i

i = 1
if (present(n)) i = n 

t = self % Tinterval(i)

end function getInterval

!--------------------------------------------------------------------------
!
! FUNCTION: getDateString
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief return the current date string 
!
!> @date   01/01/20 MDG 1.0 original
!--------------------------------------------------------------------------
function getDateString(self) result(t)

IMPLICIT NONE

class(T_TimingClass)                    :: self 

character(len=11)                       :: t

t = self % datestring

end function getDateString

!--------------------------------------------------------------------------
!
! FUNCTION: getTimeString
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief return the current time string 
!
!> @date   01/01/20 MDG 1.0 original
!--------------------------------------------------------------------------
function getTimeString(self) result(t)

IMPLICIT NONE

class(T_TimingClass)                    :: self 

character(len=15)                       :: t

t = self % timestring

end function getTimeString

!--------------------------------------------------------------------------
!
! SUBROUTINE: printTimeStamp
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief print the current date-time string to std output or redirected  
!
!> @date   01/01/20 MDG 1.0 original
!--------------------------------------------------------------------------
subroutine printTimeStamp(self, redirect) 

use mod_io

IMPLICIT NONE

class(T_TimingClass)                    :: self 
integer(kind=irg),INTENT(IN),OPTIONAL   :: redirect

type(T_IOClass)                         :: Message
integer(kind=irg)                       :: unit

Message = T_IOClass()

unit = stdout
if (present(redirect)) unit = redirect

call self % makeTimeStamp() 

if (unit.eq.stdout) then 
  call Message % printMessage( self % timestamp, frm="(/A/)")
else
  call Message % printMessage( self % timestamp, frm="(/A/)", redirect = unit)
end if

end subroutine printTimeStamp



! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE: Time_reset
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief reset time recording
! !
! !> @date   06/04/01 MDG 1.0 original
! !> @date   06/04/13 MDG 2.0 rewrite
! !> @date   06/05/14 MDG 3.0 added TT as argument
! !--------------------------------------------------------------------------
! recursive subroutine Time_reset(TT)
! !DEC$ ATTRIBUTES DLLEXPORT :: Time_reset

! IMPLICIT NONE

! type(timetype),INTENT(INOUT)	:: TT
! !f2py intent(in,out) ::  TT

! TT%TIME_t_count = 0.0
! TT%TIME_unit_count = 0.0
! TT%TIME_count = 0
! TT%TIME_newcount = 0
! TT%TIME_count_rate = 0
! TT%TIME_count_max = HUGE(0)
! TT%TIME_old = 0
! TT%TIME_loops = 0

! end subroutine Time_reset

! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE: Time_report
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief report time recording
! !
! !> @param TT time structure
! !> @param interval interval for reporting
! !
! !> @date   06/04/01 MDG 1.0 original
! !> @date   06/04/13 MDG 2.0 rewrite
! !> @date   06/05/14 MDG 3.0 added TT as argument
! !--------------------------------------------------------------------------
! recursive subroutine Time_report(TT, interval)
! !DEC$ ATTRIBUTES DLLEXPORT :: Time_report

! IMPLICIT NONE

! type(timetype),INTENT(INOUT)		:: TT
! !f2py intent(in,out) ::  TT
! real(kind=sgl),intent(IN)   		:: interval

!  TT%TIME_interval = interval
!  TT%TIME_fraction = TT%TIME_interval
!  call Message('Starting computation', frm = "(/A)")

! end subroutine Time_report

! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE: Time_start
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief start time recording
! !
! !> @date   06/04/01 MDG 1.0 original
! !> @date   06/04/13 MDG 2.0 rewrite
! !> @date   06/05/14 MDG 3.0 added TT as argument
! !--------------------------------------------------------------------------
! recursive subroutine Time_start(TT)
! !DEC$ ATTRIBUTES DLLEXPORT :: Time_start

! IMPLICIT NONE

! type(timetype),INTENT(INOUT)		:: TT
! !f2py intent(in,out) ::  TT

! ! start the timing of the computation
!  call Time_reset(TT)
!  call system_clock(TT%TIME_count,TT%TIME_count_rate,TT%TIME_count_max)

! end subroutine Time_start

! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE: Time_estimate
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief estimare remaining time
! !
! !> @param TT time structure
! !> @param numk number of idividual computations
! !
! !> @date   06/04/01 MDG 1.0 original
! !> @date   06/04/13 MDG 2.0 rewrite
! !> @date   06/05/14 MDG 3.0 added TT as argument
! !--------------------------------------------------------------------------
! recursive subroutine Time_estimate(TT, numk)
! !DEC$ ATTRIBUTES DLLEXPORT :: Time_estimate

! IMPLICIT NONE

! type(timetype),INTENT(INOUT)		:: TT
! !f2py intent(in,out) ::  TT
! integer(kind=irg),intent(IN)     	:: numk

! integer(kind=irg)      		:: TIME_nc
! real(kind=sgl)				:: io_real(1)

! ! get the current time
!  call system_clock(TIME_nc, TT%TIME_count_rate, TT%TIME_count_max)
!  TT%TIME_newcount = TIME_nc
!  TT%TIME_t_count = float(TT%TIME_newcount-TT%TIME_count)/float(TT%TIME_count_rate)
!  TT%TIME_unit_count = TT%TIME_t_count
!  io_real(1) = TT%TIME_unit_count
!  call WriteValue(' Time for first computation step [s, typically overestimate] :', io_real, 1, frm = "(F10.5)")
!  call Message('  Anticipated total computation time :', frm = "(A)",advance="no")
!  call PrintTime(TT%TIME_unit_count*float(numk))
 
! end subroutine Time_estimate

! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE: Time_estimate
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief estimate remaining time
! !
! !> @param TT time structure
! !> @param ik current computation
! !> @param numk number of idividual computations
! !
! !> @date   06/04/01 MDG 1.0 original
! !> @date   06/04/13 MDG 2.0 rewrite
! !> @date   06/05/14 MDG 3.0 added TT as argument
! !--------------------------------------------------------------------------
! recursive subroutine Time_remaining(TT, ik, numk)
! !DEC$ ATTRIBUTES DLLEXPORT :: Time_remaining

! IMPLICIT NONE

! type(timetype),INTENT(INOUT)		:: TT
! !f2py intent(in,out) ::  TT
! integer(kind=irg),intent(IN)   	:: ik
! integer(kind=irg),intent(IN)   	:: numk

! integer(kind=irg)    			:: TIME_nc, io_int(1)
! real(kind=sgl)				:: io_real(1)


!  TT%TIME_fraction = TT%TIME_fraction + TT%TIME_interval

! ! get the current time
!  call system_clock(TIME_nc, TT%TIME_count_rate, TT%TIME_count_max)

! ! correct for the resetting of TIME_nc when TIME_count_max is reached
!  if (TIME_nc.lt.TT%TIME_newcount) then     ! we've looped through the entire cycle
!    TT%TIME_loops = TT%TIME_loops+1
!    TT%TIME_count = 0
!  end if 
!  TT%TIME_newcount = TIME_nc

! ! and print it
!  TT%TIME_t_count = (float(TT%TIME_loops)*float(TT%TIME_count_max)+float(TT%TIME_newcount-TT%TIME_count))/float(TT%TIME_count_rate)

! ! reset the time per unit
!  TT%TIME_unit_count = TT%TIME_t_count/float(ik)

! ! print estimated remaining time
!  io_int(1) = nint(100.0*TT%TIME_t_count/(TT%TIME_t_count+TT%TIME_unit_count*(float(numk)-float(ik))))
!  call WriteValue (' ',io_int, 1, frm = "(1x,I3,' % completed; ')",advance="no") 
!  io_real(1) = TT%TIME_t_count
!  call WriteValue(' Total computation time [s] ', io_real, 1, frm = "(F)",advance="no")
!  call Message(';  Estimated remaining time : ', frm = "(A)",advance="no")
!  call PrintTime(TT%TIME_unit_count*(float(numk)-float(ik)))
! !
! end subroutine Time_remaining

! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE: PrintTime
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief print  time
! !
! !> @param tm time variable
! !
! !> @date   06/04/01 MDG 1.0 original
! !> @date   06/04/13 MDG 2.0 rewrite
! !> @date   06/05/14 MDG 3.0 changed IO
! !--------------------------------------------------------------------------
! recursive subroutine PrintTime(tm)
! !DEC$ ATTRIBUTES DLLEXPORT :: PrintTime

! IMPLICIT NONE

! real(kind=sgl),INTENT(IN)		:: tm

! integer(kind=irg)    			:: days, hours, minutes, seconds, io_int(4)
! real(kind=sgl)       			:: secs

!   secs = tm
!   days = 0
!   hours = 0
!   minutes = 0
!   if (secs.gt.86400.0) then
!     days = int(secs)/86400
!     secs = mod(secs,86400.0)
!   end if
!   if (secs.gt.3600.0) then
!     hours = int(secs)/3600
!     secs = mod(secs,3600.0)
!   end if
!   if (secs.gt.60.0) then
!     minutes = int(secs)/60
!     secs = mod(secs,60.0)
!   end if
!   seconds = int(secs)

!   io_int(1:4) = (/ days, hours, minutes, seconds /)
!   call WriteValue(' ',io_int, 4, frm = "(1x,I3,' d,',I3,' h,',I3,' m,',I3,' s')")

! end subroutine PrintTime

! !--------------------------------------------------------------------------
! !
! ! SUBROUTINE: Time_stop
! !
! !> @author Marc De Graef, Carnegie Mellon University
! !
! !> @brief stop time recording
! !
! !> @param TT time structure
! !> @param numk total number of computations
! !
! !> @date   06/04/01 MDG 1.0 original
! !> @date   06/04/13 MDG 2.0 rewrite
! !> @date   06/05/14 MDG 3.0 added TT; changed IO
! !--------------------------------------------------------------------------
! recursive subroutine Time_stop(TT, numk)
! !DEC$ ATTRIBUTES DLLEXPORT :: Time_stop

! IMPLICIT NONE

! type(timetype),INTENT(INOUT)		:: TT
! !f2py intent(in,out) ::  TT
! integer(kind=irg),INTENT(IN)  		:: numk

! real(kind=sgl)				:: io_real(1)


!   call system_clock(TT%TIME_newcount, TT%TIME_count_rate, TT%TIME_count_max)
!   call Message('  Total computation time [s] ', frm = "(A)",advance="no")
!   call PrintTime((float(TT%TIME_loops)*float(TT%TIME_count_max)+float(TT%TIME_newcount-TT%TIME_count))/float(TT%TIME_count_rate))
!   io_real(1)= float(TT%TIME_loops)*float(TT%TIME_count_max)+float(TT%TIME_newcount-TT%TIME_count)
!   io_real(1) = io_real(1)/float(TT%TIME_count_rate)/float(numk)
!   call WriteValue(' Time per step/pixel [s] ', io_real, 1, frm = "(F)")
! end subroutine Time_stop


!--------------------------------------------------------------------------
!
! subroutine: makeTimeStamp
!
!> @brief sets the current YMDHMS date as a time stamp.
!
!> @date 01/01/20 MDG based on timestamp original routine by J. Burkardt, but adapted for timing class
!>                    <https://people.sc.fsu.edu/~jburkardt/f_src/timestamp/timestamp.f90>
!--------------------------------------------------------------------------
subroutine makeTimeStamp (self)

  IMPLICIT NONE

  class(T_TimingClass),intent(inout)    :: self

  integer(kind=irg)                     :: d, h, mo, mm, n, s, values(8), y
  character ( len = 8 )                 :: ampm, date
  character ( len = 10 )                :: time
  character ( len = 5 )                 :: zone

! call the intrinsic routine 
  call date_and_time ( date, time, zone, values )

! extract values  (based on JB code)
  y = values(1)
  mo = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if (h.lt.12) then
    ampm = 'AM'
  else if (h.eq.12) then
    if ((n.eq.0).and.(s.eq.0)) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if (h.lt.12) then
      ampm = 'PM'
    else if (h.eq.12) then
      if ((n.eq.0).and.(s.eq.0)) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if
! end of original JB code

  write (self % datestring, '(a,1x,i2,1x,i4)' ) month(mo), d, y
  write (self % timestring, '(i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) h,':',n,':',s,'.',mm,trim(ampm)
  self % timestamp = self % datestring // self % timestring

end subroutine makeTimeStamp

end module mod_timing
