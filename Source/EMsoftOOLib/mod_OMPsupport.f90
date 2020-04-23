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

module mod_OMPsupport
  !! author: MDG 
  !! version: 1.0 
  !! date: 04/15/20
  !!
  !! a few miscellaneous OpenMP support routines

use mod_kinds
use mod_global
use omp_lib
use mod_io

IMPLICIT NONE 

integer(kind=irg), private :: maxOMPthreads

contains

!--------------------------------------------------------------------------
subroutine OMP_showAvailableThreads() 
!! author: MDG 
!! version: 1.0 
!! date: 04/15/20
!!
!! display the available number of OpenMP threads 
 
IMPLICIT NONE

type(IO_T)        :: Message 
integer(kind=irg) :: io_int(1) 

maxOMPthreads = OMP_GET_MAX_THREADS()
io_int(1) = maxOMPthreads
call Message%WriteValue(' Maximum available number of OpenMP threads : ', io_int, 1)

end subroutine OMP_showAvailableThreads

!--------------------------------------------------------------------------
subroutine OMP_showAllocatedThreads() 
!! author: MDG 
!! version: 1.0 
!! date: 04/15/20
!!
!! display the available number of OpenMP threads 
 
IMPLICIT NONE

type(IO_T)        :: Message 
integer(kind=irg) :: io_int(1) 

io_int(1) = OMP_GET_NUM_THREADS()
call Message%WriteValue(' Allocated number of OpenMP threads : ', io_int, 1)

end subroutine OMP_showAllocatedThreads

!--------------------------------------------------------------------------
subroutine OMP_setNThreads(n) 
!! author: MDG 
!! version: 1.0 
!! date: 04/15/20
!!
!! set the number of OpenMP threads; check that they are available
 
IMPLICIT NONE

integer(kind=irg),INTENT(IN)  :: n

type(IO_T)                    :: Message 
integer(kind=irg)             :: io_int(1) 

if (n.eq.0) then 
  io_int(1) = maxOMPthreads
  call Message%WriteValue(' Number of OpenMP threads set to maximum available : ', io_int, 1)
  call OMP_SET_NUM_THREADS(n) 
else
  if (maxOMPthreads.lt.n) then 
    io_int(1) = n 
    call Message%WriteValue(' Number of OpenMP threads requested : ', io_int, 1)
    io_int(1) =  maxOMPthreads
    call Message%WriteValue(' Number of OpenMP threads available : ', io_int, 1)
    call Message%printMessage(' --> Setting number of threads to maximum available ')
    call OMP_SET_NUM_THREADS(maxOMPthreads) 
  else
    io_int(1) = n 
    call Message%WriteValue(' Number of OpenMP threads set to ', io_int, 1)
    call OMP_SET_NUM_THREADS(n) 
  end if
end if

end subroutine OMP_setNThreads

end module mod_OMPsupport