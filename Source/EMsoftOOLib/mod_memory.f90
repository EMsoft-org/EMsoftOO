! ###################################################################
! Copyright (c) 2014-2021, Marc De Graef Research Group/Carnegie Mellon University
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
module mod_memory
  !! author: MDG
  !! version: 1.0
  !! date: 12/11/20
  !!
  !! Various memory allocation/deallocation/initialization routines.  These 
  !! routines replace the ordinary f90 allocate and deallocate calls to 
  !! make sure that arrays are only allocated if they don't already exist,
  !! and that they are properly initalized to zero or whatever the default value
  !! might be (optional argument).  This will hopefully help alleviate hard
  !! to track runtime bugs, in particular on the Windows platform.
  !!
  !! Note: March 2021, started to add support for memory allocation inside 
  !! OpenMP parallel regions, where each thread needs to allocate its own
  !! arrays. This requires that the constructor be called with the number of 
  !! OMP threads.
  !!
  !! This class is not meant to be used to track all memory use for every program 
  !! and module... it should only be used within a single program and will allow 
  !! the developer to track general memory use from explicit array allocations as 
  !! well as array allocations inside parallel program sections.  Memory allocated 
  !! inside a routine that is called from the main program will not be tracked
  !! (unless the routine is properly instrumented).
  !!
  !! Usage example:  [without threading]
  !! program t 
  !!
  !! use mod_kinds
  !! use mod_global
  !! use mod_memory 
  !!
  !! IMPLICIT NONE 
  !! 
  !! type(memory_T)             :: mem 
  !! real(kind=sgl),allocatable :: ar(:,:)
  !!
  !! ! initiate the memory class 
  !! mem = memory_T() 
  !!
  !! ! allocate the array and initialize with value 10.0
  !! call mem%alloc( ar, (/ 20, 30 /), 'ar', 10.0)
  !! 
  !! ! find out how much memory has been allocated and what type 
  !! call mem%allocated_memory_use()
  !! 
  !! ! deallocate the array 
  !! call mem%dealloc(ar, 'ar')
  !!
  !! ! make sure it was indeed deallocated 
  !! call mem%allocated_memory_use()
  !!
  !! end program 
  !!--------------------------------------------------
  !!
  !! usage example: [with threading]
  !! program t 
  !!
  !! use mod_kinds
  !! use mod_global
  !! use omp_lib
  !! use mod_memory 
  !!
  !! IMPLICIT NONE 
  !! 
  !! type(memory_T)                 :: memth 
  !! real(kind=sgl),allocatable     :: ar(:,:)
  !! complex(kind=dbl),allocatable  :: cc(:)
  !! integer(kind=irg)              :: nthreads, TID 
  !!
  !! nthreads = 3 
  !!
  !! ! initiate the memory class 
  !! memth = memory_T( nt = nthreads ) 
  !!
  !! ! inside the parallel OMP region, allocate the array for each thread 
  !! ! we'll assume that TID contains the thread ID 
  !! call memth%alloc( ar, (/ 10, 10 /), 'ar', 15.0, TID=TID)
  !! call memth%alloc( cc, (/ 10 + 5*TID /), 'cc', TID=TID)
  !!
  !! ! print memory usage information
  !! if (TID.eq.0) call memth%thread_memory_use()
  !! 
  !! ! do stuff with the arrays
  !! ! then just before closing the parallel section 
  !! call memth%dealloc(ar, 'ar', TID=TID)
  !! call memth%dealloc(cc, 'cc', TID=TID)
  !! ! end parallel region 
  !! 
  !! end program 
  !! 
  !!--------------------------------------------------

use mod_kinds
use mod_global
use mod_io

IMPLICIT NONE

type, public :: memory_T
  private 
  ! we'll keep track of all memory that has been allocated by calls to 
  ! this class; so, we'll have memory usage reporting as well.
    integer(kind=irg), allocatable :: totmem_ish(:) 
    integer(kind=irg), allocatable :: totmem_irg(:)
    integer(kind=irg), allocatable :: totmem_ill(:) 
    integer(kind=irg), allocatable :: totmem_sgl(:) 
    integer(kind=irg), allocatable :: totmem_dbl(:) 
    integer(kind=irg), allocatable :: totmem_cmplx(:)
    integer(kind=irg), allocatable :: totmem_cmplxd(:)
    integer(kind=irg), allocatable :: totmem_char(:)
    integer(kind=irg), allocatable :: totmem_logical(:)
    integer(kind=irg), allocatable :: current_memory_allocated(:)
    integer(kind=irg) :: bytes_ish = 2
    integer(kind=irg) :: bytes_irg = 4
    integer(kind=irg) :: bytes_ill = 8
    integer(kind=irg) :: bytes_sgl = 4
    integer(kind=irg) :: bytes_dbl = 8
    integer(kind=irg) :: bytes_cmplx = 8
    integer(kind=irg) :: bytes_cmplxd = 16
    integer(kind=irg) :: bytes_char = 1
    integer(kind=irg) :: bytes_logical = 1
    integer(kind=irg), public :: nthreads = 1
    logical, public   :: verbose

  contains
  private 

    procedure, pass(self) :: alloc_char1_
    procedure, pass(self) :: alloc_char2_
    procedure, pass(self) :: alloc_ish1_
    procedure, pass(self) :: alloc_ish2_
    procedure, pass(self) :: alloc_ish3_
    procedure, pass(self) :: alloc_ish4_
    procedure, pass(self) :: alloc_ish5_
    procedure, pass(self) :: alloc_ish6_
    procedure, pass(self) :: alloc_irg1_
    procedure, pass(self) :: alloc_irg2_
    procedure, pass(self) :: alloc_irg3_
    procedure, pass(self) :: alloc_irg4_
    procedure, pass(self) :: alloc_irg5_
    procedure, pass(self) :: alloc_irg6_
    procedure, pass(self) :: alloc_ill1_
    procedure, pass(self) :: alloc_ill2_
    procedure, pass(self) :: alloc_ill3_
    procedure, pass(self) :: alloc_ill4_
    procedure, pass(self) :: alloc_ill5_
    procedure, pass(self) :: alloc_ill6_
    procedure, pass(self) :: alloc_sgl1_
    procedure, pass(self) :: alloc_sgl2_
    procedure, pass(self) :: alloc_sgl3_
    procedure, pass(self) :: alloc_sgl4_
    procedure, pass(self) :: alloc_sgl5_
    procedure, pass(self) :: alloc_sgl6_
    procedure, pass(self) :: alloc_dbl1_
    procedure, pass(self) :: alloc_dbl2_
    procedure, pass(self) :: alloc_dbl3_
    procedure, pass(self) :: alloc_dbl4_
    procedure, pass(self) :: alloc_dbl5_
    procedure, pass(self) :: alloc_dbl6_
    procedure, pass(self) :: alloc_cmplx1_
    procedure, pass(self) :: alloc_cmplx2_
    procedure, pass(self) :: alloc_cmplx3_
    procedure, pass(self) :: alloc_cmplx4_
    procedure, pass(self) :: alloc_cmplx5_
    procedure, pass(self) :: alloc_cmplx6_
    procedure, pass(self) :: alloc_cmplxd1_
    procedure, pass(self) :: alloc_cmplxd2_
    procedure, pass(self) :: alloc_cmplxd3_
    procedure, pass(self) :: alloc_cmplxd4_
    procedure, pass(self) :: alloc_cmplxd5_
    procedure, pass(self) :: alloc_cmplxd6_
    procedure, pass(self) :: alloc_logical1_
    procedure, pass(self) :: alloc_logical2_
    procedure, pass(self) :: alloc_logical3_

    procedure, pass(self) :: dealloc_char1_
    procedure, pass(self) :: dealloc_char2_
    procedure, pass(self) :: dealloc_ish1_
    procedure, pass(self) :: dealloc_ish2_
    procedure, pass(self) :: dealloc_ish3_
    procedure, pass(self) :: dealloc_ish4_
    procedure, pass(self) :: dealloc_ish5_
    procedure, pass(self) :: dealloc_ish6_
    procedure, pass(self) :: dealloc_irg1_
    procedure, pass(self) :: dealloc_irg2_
    procedure, pass(self) :: dealloc_irg3_
    procedure, pass(self) :: dealloc_irg4_
    procedure, pass(self) :: dealloc_irg5_
    procedure, pass(self) :: dealloc_irg6_
    procedure, pass(self) :: dealloc_ill1_
    procedure, pass(self) :: dealloc_ill2_
    procedure, pass(self) :: dealloc_ill3_
    procedure, pass(self) :: dealloc_ill4_
    procedure, pass(self) :: dealloc_ill5_
    procedure, pass(self) :: dealloc_ill6_
    procedure, pass(self) :: dealloc_sgl1_
    procedure, pass(self) :: dealloc_sgl2_
    procedure, pass(self) :: dealloc_sgl3_
    procedure, pass(self) :: dealloc_sgl4_
    procedure, pass(self) :: dealloc_sgl5_
    procedure, pass(self) :: dealloc_sgl6_
    procedure, pass(self) :: dealloc_dbl1_
    procedure, pass(self) :: dealloc_dbl2_
    procedure, pass(self) :: dealloc_dbl3_
    procedure, pass(self) :: dealloc_dbl4_
    procedure, pass(self) :: dealloc_dbl5_
    procedure, pass(self) :: dealloc_dbl6_
    procedure, pass(self) :: dealloc_cmplx1_
    procedure, pass(self) :: dealloc_cmplx2_
    procedure, pass(self) :: dealloc_cmplx3_
    procedure, pass(self) :: dealloc_cmplx4_
    procedure, pass(self) :: dealloc_cmplx5_
    procedure, pass(self) :: dealloc_cmplx6_
    procedure, pass(self) :: dealloc_cmplxd1_
    procedure, pass(self) :: dealloc_cmplxd2_
    procedure, pass(self) :: dealloc_cmplxd3_
    procedure, pass(self) :: dealloc_cmplxd4_
    procedure, pass(self) :: dealloc_cmplxd5_
    procedure, pass(self) :: dealloc_cmplxd6_
    procedure, pass(self) :: dealloc_logical1_
    procedure, pass(self) :: dealloc_logical2_
    procedure, pass(self) :: dealloc_logical3_

    procedure, pass(self) :: compute_size_
    procedure, pass(self) :: show_allocated_memory_use_ 
    procedure, pass(self) :: show_thread_memory_use_
    procedure, pass(self) :: update_total_memory_use_ 
    procedure, pass(self) :: toggle_verbose_ 
    final :: memory_destructor 

! overload all the allocation and deallocation routines into simple alloc and dealloc methods
    generic, public :: alloc => alloc_char1_, alloc_char2_, alloc_ish1_, alloc_irg1_, alloc_ill1_, &
                                alloc_sgl1_, alloc_dbl1_, alloc_cmplx1_, &
                                alloc_cmplxd1_, alloc_ish2_, alloc_irg2_, alloc_ill2_, &
                                alloc_sgl2_, alloc_dbl2_, alloc_cmplx2_, &
                                alloc_cmplxd2_, alloc_ish3_, alloc_irg3_, alloc_ill3_, &
                                alloc_sgl3_, alloc_dbl3_, alloc_cmplx3_, &
                                alloc_cmplxd3_, alloc_ish4_, alloc_irg4_, alloc_ill4_, &
                                alloc_sgl4_, alloc_dbl4_, alloc_cmplx4_, &
                                alloc_cmplxd4_, alloc_ish5_, alloc_irg5_, alloc_ill5_, &
                                alloc_sgl5_, alloc_dbl5_, alloc_cmplx5_, &
                                alloc_cmplxd5_, alloc_ish6_, alloc_irg6_, alloc_ill6_, &
                                alloc_sgl6_, alloc_dbl6_, alloc_cmplx6_, &
                                alloc_cmplxd6_, alloc_logical1_, alloc_logical2_, alloc_logical3_ 

    generic, public :: dealloc => dealloc_char1_, dealloc_char2_, dealloc_ish1_, dealloc_irg1_, dealloc_ill1_, &
                                  dealloc_sgl1_, dealloc_dbl1_, dealloc_cmplx1_, &
                                  dealloc_cmplxd1_, dealloc_ish2_, dealloc_irg2_, dealloc_ill2_, &
                                  dealloc_sgl2_, dealloc_dbl2_, dealloc_cmplx2_, &
                                  dealloc_cmplxd2_,dealloc_ish3_, dealloc_irg3_, dealloc_ill3_, &
                                  dealloc_sgl3_, dealloc_dbl3_, dealloc_cmplx3_, &
                                  dealloc_cmplxd3_, dealloc_ish4_, dealloc_irg4_, dealloc_ill4_, &
                                  dealloc_sgl4_, dealloc_dbl4_, dealloc_cmplx4_, &
                                  dealloc_cmplxd4_,dealloc_ish5_, dealloc_irg5_, dealloc_ill5_, &
                                  dealloc_sgl5_, dealloc_dbl5_, dealloc_cmplx5_, &
                                  dealloc_cmplxd5_,dealloc_ish6_, dealloc_irg6_, dealloc_ill6_, &
                                  dealloc_sgl6_, dealloc_dbl6_, dealloc_cmplx6_, &
                                  dealloc_cmplxd6_, dealloc_logical1_, dealloc_logical2_, dealloc_logical3_ 

    generic, public :: getsize => compute_size_
    generic, public :: allocated_memory_use => show_allocated_memory_use_
    generic, public :: thread_memory_use => show_thread_memory_use_

    generic, public :: toggle_verbose => toggle_verbose_
end type memory_T

! constructor for this class
interface memory_T 
  module procedure memory_constructor 
end interface memory_T 

contains

!--------------------------------------------------------------------------
type(memory_T) function memory_constructor( nt, silent )  result(mem) 
!DEC$ ATTRIBUTES DLLEXPORT :: memory_constructor
!! author: MDG
!! version: 1.0
!! date: 12/11/20
!!
!! constructor for the memory_T Class

use mod_io
use mod_EMsoft

IMPLICIT NONE

integer(kind=irg), INTENT(IN), OPTIONAL  :: nt 
logical, INTENT(IN), OPTIONAL            :: silent 

type(IO_T)                               :: Message
type(EMsoft_T)                           :: EMsoft
character(fnlen)                         :: dummy 

! deallocate any old memory usage counter arrays 
if (allocated(mem%totmem_ish)) deallocate(mem%totmem_ish) 
if (allocated(mem%totmem_irg)) deallocate(mem%totmem_irg) 
if (allocated(mem%totmem_ill)) deallocate(mem%totmem_ill) 
if (allocated(mem%totmem_sgl)) deallocate(mem%totmem_sgl) 
if (allocated(mem%totmem_dbl)) deallocate(mem%totmem_dbl) 
if (allocated(mem%totmem_cmplx)) deallocate(mem%totmem_cmplx) 
if (allocated(mem%totmem_cmplxd)) deallocate(mem%totmem_cmplxd) 
if (allocated(mem%totmem_char)) deallocate(mem%totmem_char)
if (allocated(mem%totmem_logical)) deallocate(mem%totmem_logical)
if (allocated(mem%current_memory_allocated)) deallocate(mem%current_memory_allocated) 

! allocate for nthreads OpenMP threads or just a single scalar 
if (present(nt)) then 
    mem%nthreads = nt
    allocate(mem%totmem_ish(0:nt-1))
    allocate(mem%totmem_irg(0:nt-1))
    allocate(mem%totmem_ill(0:nt-1))
    allocate(mem%totmem_sgl(0:nt-1))
    allocate(mem%totmem_dbl(0:nt-1))
    allocate(mem%totmem_cmplx(0:nt-1))
    allocate(mem%totmem_cmplxd(0:nt-1))
    allocate(mem%totmem_char(0:nt-1))
    allocate(mem%totmem_logical(0:nt-1))
    allocate(mem%current_memory_allocated(0:nt-1))
else 
    mem%nthreads = 1
    allocate(mem%totmem_ish(1))
    allocate(mem%totmem_irg(1))
    allocate(mem%totmem_ill(1))
    allocate(mem%totmem_sgl(1))
    allocate(mem%totmem_dbl(1))
    allocate(mem%totmem_cmplx(1))
    allocate(mem%totmem_cmplxd(1))
    allocate(mem%totmem_char(1))
    allocate(mem%totmem_logical(1))
    allocate(mem%current_memory_allocated(1))
end if

! set all the memory usage counters to zero
mem%totmem_ish  = 0
mem%totmem_irg = 0
mem%totmem_ill  = 0
mem%totmem_sgl  = 0
mem%totmem_dbl  = 0
mem%totmem_cmplx  = 0
mem%totmem_cmplxd  = 0
mem%totmem_char  = 0
mem%totmem_logical  = 0
mem%current_memory_allocated = 0

mem%verbose = .FALSE.
if (.not.present(silent)) then 
    dummy = ''
    EMsoft = EMsoft_T(dummy,dummy,silent=.TRUE.)
    if (EMsoft%getConfigParameter('EMsoftAllocatetest').eq.'Yes') then 
      call Message%printMessage('>>>>>>>>>>>>>>>  Turning on verbosity in mod_memory class !!!!!!!!!')
      mem%verbose = .TRUE.
    else 
      mem%verbose = .FALSE.
    end if 
end if 

end function memory_constructor

!--------------------------------------------------------------------------
subroutine memory_destructor(self)
!DEC$ ATTRIBUTES DLLEXPORT :: memory_destructor
!! author: MDG
!! version: 1.0
!! date: 03/25/21
!!
!! destructor for the memory_T Class

IMPLICIT NONE

type(memory_T), INTENT(INOUT)  :: self

! deallocate any old memory usage counter arrays 
if (allocated(self%totmem_ish)) deallocate(self%totmem_ish) 
if (allocated(self%totmem_irg)) deallocate(self%totmem_irg) 
if (allocated(self%totmem_ill)) deallocate(self%totmem_ill) 
if (allocated(self%totmem_sgl)) deallocate(self%totmem_sgl) 
if (allocated(self%totmem_dbl)) deallocate(self%totmem_dbl) 
if (allocated(self%totmem_cmplx)) deallocate(self%totmem_cmplx) 
if (allocated(self%totmem_cmplxd)) deallocate(self%totmem_cmplxd) 
if (allocated(self%totmem_char)) deallocate(self%totmem_char) 
if (allocated(self%current_memory_allocated)) deallocate(self%current_memory_allocated) 

call reportDestructor('memory_T')

end subroutine memory_destructor

!--------------------------------------------------------------------------
subroutine show_allocated_memory_use_( self, expl )
!DEC$ ATTRIBUTES DLLEXPORT :: show_allocated_memory_use_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! print out the current memory usage 

IMPLICIT NONE

class(memory_T), INTENT(INOUT)      :: self
character(*),INTENT(IN),OPTIONAL    :: expl 

type(IO_T)                          :: Message 

integer(kind=irg)                   :: io_int(1) 

call Message%printMessage(' Allocated Memory Status (in bytes) : ',"(/A)")
if (present(expl)) then 
  call Message%printMessage(' [ '//trim(expl)//' ]',"(A)")
end if 
call Message%printMessage(' ----------------------------------   ',"(A)")
io_int = sum(self%totmem_ish)
call Message%WriteValue(' integer(kind=ish)     : ', io_int, 1)
io_int = sum(self%totmem_irg)
call Message%WriteValue(' integer(kind=irg)     : ', io_int, 1)
io_int = sum(self%totmem_ill)
call Message%WriteValue(' integer(kind=ill)     : ', io_int, 1)
io_int = sum(self%totmem_sgl)
call Message%WriteValue(' real(kind=sgl)        : ', io_int, 1)
io_int = sum(self%totmem_dbl)
call Message%WriteValue(' real(kind=dbl)        : ', io_int, 1)
io_int = sum(self%totmem_cmplx)
call Message%WriteValue(' complex(kind=sgl)     : ', io_int, 1)
io_int = sum(self%totmem_cmplxd)
call Message%WriteValue(' complex(kind=dbl)     : ', io_int, 1)
io_int = sum(self%totmem_char)
call Message%WriteValue(' character             : ', io_int, 1)


call Message%printMessage(' ')
io_int = sum(self%current_memory_allocated)
call Message%WriteValue(' total allocated bytes (unthreaded) : ', io_int, 1)

end subroutine show_allocated_memory_use_

!--------------------------------------------------------------------------
subroutine show_thread_memory_use_( self, expl )
!DEC$ ATTRIBUTES DLLEXPORT :: show_thread_memory_use_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! print out the current memory usage 

IMPLICIT NONE

class(memory_T), INTENT(INOUT)      :: self
character(*),INTENT(IN),OPTIONAL    :: expl 

type(IO_T)                          :: Message 

integer(kind=irg)                   :: io_int(1), i 

call Message%printMessage(' Allocated Memory Status (in bytes) : ',"(/A)")
if (present(expl)) then 
  call Message%printMessage(' [ '//trim(expl)//' ]',"(A)")
end if 
call Message%printMessage(' ----------------------------------   ',"(A)")
if (self%nthreads.gt.1) then 
    do i=0, self%nthreads-1
        io_int = i
        call Message%WriteValue(' Allocations for OpenMP thread ', io_int, 1)
        io_int = self%totmem_ish(i)
        call Message%WriteValue('    integer(kind=ish) : ', io_int, 1)
        io_int = self%totmem_irg(i)
        call Message%WriteValue('    integer(kind=irg) : ', io_int, 1)
        io_int = self%totmem_ill(i)
        call Message%WriteValue('    integer(kind=ill) : ', io_int, 1)
        io_int = self%totmem_sgl(i)
        call Message%WriteValue('    real(kind=sgl)    : ', io_int, 1)
        io_int = self%totmem_dbl(i)
        call Message%WriteValue('    real(kind=dbl)    : ', io_int, 1)
        io_int = self%totmem_cmplx(i)
        call Message%WriteValue('    complex(kind=sgl) : ', io_int, 1)
        io_int = self%totmem_cmplxd(i)
        call Message%WriteValue('    complex(kind=dbl) : ', io_int, 1)
        io_int = self%totmem_char(i)
        call Message%WriteValue(' character             : ', io_int, 1)
        io_int = self%totmem_logical(i)
        call Message%WriteValue(' logical               : ', io_int, 1)
        call Message%printMessage(' ')
    end do 
else
    io_int = self%totmem_ish(1)
    call Message%WriteValue('    integer(kind=ish) : ', io_int, 1)
    io_int = self%totmem_irg(1)
    call Message%WriteValue('    integer(kind=irg) : ', io_int, 1)
    io_int = self%totmem_ill(1)
    call Message%WriteValue('    integer(kind=ill) : ', io_int, 1)
    io_int = self%totmem_sgl(1)
    call Message%WriteValue('    real(kind=sgl)    : ', io_int, 1)
    io_int = self%totmem_dbl(1)
    call Message%WriteValue('    real(kind=dbl)    : ', io_int, 1)
    io_int = self%totmem_cmplx(1)
    call Message%WriteValue('    complex(kind=sgl) : ', io_int, 1)
    io_int = self%totmem_cmplxd(1)
    call Message%WriteValue('    complex(kind=dbl) : ', io_int, 1)
    io_int = self%totmem_char(1)
    call Message%WriteValue(' character             : ', io_int, 1)
    io_int = self%totmem_logical(1)
    call Message%WriteValue(' logical               : ', io_int, 1)
    call Message%printMessage(' ')
end if 
    
io_int = sum(self%current_memory_allocated)
call Message%WriteValue(' total allocated bytes in threads : ', io_int, 1)

end subroutine show_thread_memory_use_

!--------------------------------------------------------------------------
subroutine update_total_memory_use_( self, value, TID )
!DEC$ ATTRIBUTES DLLEXPORT :: update_total_memory_use_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, updated 03/24/21
!!
!! updates the total memory use

IMPLICIT NONE

class(memory_T), INTENT(INOUT)          :: self
integer(kind=irg),INTENT(IN)            :: value 
integer(kind=irg),INTENT(IN),OPTIONAL   :: TID

integer(kind=irg)                       :: LID 

LID = 1
if (present(TID)) LID = TID

self%current_memory_allocated(LID) = self%current_memory_allocated(LID) + value

end subroutine update_total_memory_use_

!--------------------------------------------------------------------------
function compute_size_(self, dims, tp) result(nbytes)
!DEC$ ATTRIBUTES DLLEXPORT :: compute_size_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! determine the number of bytes in an array of type ish without allocating it 

IMPLICIT NONE

class(memory_T), INTENT(IN)         :: self
integer(kind=irg), INTENT(IN)       :: dims(:)
character(6), INTENT(IN)            :: tp

integer(kind=irg)                   :: nbytes

select case(trim(tp))
    case('ish') 
        nbytes = product(dims) * self%bytes_ish
    case('irg') 
        nbytes = product(dims) * self%bytes_irg
    case('ill') 
        nbytes = product(dims) * self%bytes_ill
    case('sgl') 
        nbytes = product(dims) * self%bytes_sgl
    case('dbl') 
        nbytes = product(dims) * self%bytes_dbl
    case('cmplx') 
        nbytes = product(dims) * self%bytes_cmplx
    case('cmplxd') 
        nbytes = product(dims) * self%bytes_cmplxd
    case('char') 
        nbytes = dims(1) * self%bytes_char
    case('logical') 
        nbytes = dims(1) * self%bytes_logical
end select 

end function compute_size_

!--------------------------------------------------------------------------
subroutine toggle_verbose_(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: toggle_verbose_
!! author: MDG
!! version: 1.0
!! date: 12/14/20
!!
!! toggle the verbose parameter 

IMPLICIT NONE

class(memory_T), INTENT(INOUT)      :: self

if (self%verbose.eqv..TRUE.) then 
    self%verbose = .FALSE.
else
    self%verbose = .TRUE.
end if

end subroutine toggle_verbose_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Here we have the individual routines for all relevant data types 
! Since f90 does not have templates, we need separate routines for all types.
! 
! For the character type, we only do 1D and 2D arrays; it is unlikely that we 
! will ever need a 3D array of strings...
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
recursive subroutine alloc_char1_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_char1_
!! author: MDG
!! version: 1.0
!! date: 03/28/21
!!
!! allocate an array of type char

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
character(fnlen), INTENT(INOUT), allocatable     :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
character(*), INTENT(IN), OPTIONAL               :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(1)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(1) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ish
    self%totmem_ish(LID) = self%totmem_ish(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_char1_:', &
      ' Unable to allocate character(fnlen) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_char(LID) = self%totmem_char(LID) + (dims(1)-startdims(1)+1)*fnlen
  call self%update_total_memory_use_((dims(1)-startdims(1)+1)*fnlen, LID)
else
  allocate(ar(dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_char1_:', &
        ' Unable to allocate character(fnlen) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_char(LID) = self%totmem_char(LID) + dims(1)*fnlen
  call self%update_total_memory_use_(dims(1)*fnlen, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = trim(initval)
else 
    ar = ''
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) (dims(1)-startdims(1)+1)*fnlen
  else 
    write (szstr,*) dims(1)*fnlen
  end if
  if (present(initval)) then 
    write (initstr,*) trim(initval)
  else 
    write (initstr,*) ''
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_char1_

!--------------------------------------------------------------------------
subroutine dealloc_char1_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_char1_
!! author: MDG
!! version: 1.0
!! date: 03/28/21
!!
!! deallocate an array of type char

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
character(fnlen), INTENT(INOUT), allocatable     :: ar(:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar)*fnlen 
    self%totmem_char(LID) = self%totmem_char(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_char1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_char1_

!--------------------------------------------------------------------------
recursive subroutine alloc_char2_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_char2_
!! author: MDG
!! version: 1.0
!! date: 04/14/21
!!
!! allocate an array of type char

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
character(fnlen), INTENT(INOUT), allocatable     :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
character(*), INTENT(IN), OPTIONAL               :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(2)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(2) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ish
    self%totmem_ish(LID) = self%totmem_ish(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_char2_:', &
      ' Unable to allocate character(fnlen) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_char(LID) = self%totmem_char(LID) + (dims(1)-startdims(1)+1)*(dims(2)-startdims(2)+1)*fnlen
  call self%update_total_memory_use_((dims(1)-startdims(1)+1)*(dims(2)-startdims(2)+1)*fnlen, LID)
else
  allocate(ar(dims(1),dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_char2_:', &
        ' Unable to allocate character(fnlen) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_char(LID) = self%totmem_char(LID) + dims(1)*dims(2)*fnlen
  call self%update_total_memory_use_(dims(1)*dims(2)*fnlen, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = trim(initval)
else 
    ar = ''
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) (dims(1)-startdims(1)+1)*(dims(2)-startdims(2)+1)*fnlen
  else 
    write (szstr,*) dims(1)*dims(2)*fnlen
  end if
  if (present(initval)) then 
    write (initstr,*) trim(initval)
  else 
    write (initstr,*) ''
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_char2_

!--------------------------------------------------------------------------
subroutine dealloc_char2_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_char2_
!! author: MDG
!! version: 1.0
!! date: 04/14/21
!!
!! deallocate an array of type char

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
character(fnlen), INTENT(INOUT), allocatable     :: ar(:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar)*fnlen 
    self%totmem_char(LID) = self%totmem_char(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_char2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_char2_

!--------------------------------------------------------------------------
recursive subroutine alloc_ish1_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ish1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type ish1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable       :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
integer(kind=ish), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(1)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(1) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ish
    self%totmem_ish(LID) = self%totmem_ish(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_ish1_:', &
      ' Unable to allocate integer(kind=ish) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_ish(LID) = self%totmem_ish(LID) + product(dims-startdims+1)*self%bytes_ish 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_ish, LID)
else
  allocate(ar(dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_ish1_:', &
        ' Unable to allocate integer(kind=ish) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_ish(LID) = self%totmem_ish(LID) + product(dims)*self%bytes_ish 
  call self%update_total_memory_use_(product(dims)*self%bytes_ish, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ish
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_ish
  else 
    write (szstr,*) product(dims)*self%bytes_ish
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_ish
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_ish1_

!--------------------------------------------------------------------------
subroutine dealloc_ish1_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ish1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type ish1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable    :: ar(:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ish
    self%totmem_ish(LID) = self%totmem_ish(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ish1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ish1_

!--------------------------------------------------------------------------
recursive subroutine alloc_ish2_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ish2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type ish2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable       :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
integer(kind=ish), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(2)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(2) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ish
    self%totmem_ish(LID) = self%totmem_ish(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_ish2_:', &
      ' Unable to allocate integer(kind=ish) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_ish(LID) = self%totmem_ish(LID) + product(dims-startdims+1)*self%bytes_ish 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_ish, LID)
else
  allocate(ar(dims(1),dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_ish2_:', &
        ' Unable to allocate integer(kind=ish) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_ish(LID) = self%totmem_ish(LID) + product(dims)*self%bytes_ish 
  call self%update_total_memory_use_(product(dims)*self%bytes_ish, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ish
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_ish
  else 
    write (szstr,*) product(dims)*self%bytes_ish
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_ish
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_ish2_

!--------------------------------------------------------------------------
subroutine dealloc_ish2_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ish2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type ish2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable    :: ar(:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ish
    self%totmem_ish(LID) = self%totmem_ish(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ish2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ish2_

!--------------------------------------------------------------------------
recursive subroutine alloc_ish3_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ish3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type ish3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable       :: ar(:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(3)
character(*),INTENT(IN)                          :: varname
integer(kind=ish), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(3)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(3) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ish
    self%totmem_ish(LID) = self%totmem_ish(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_ish3_:', &
      ' Unable to allocate integer(kind=ish) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_ish(LID) = self%totmem_ish(LID) + product(dims-startdims+1)*self%bytes_ish 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_ish, LID)
else
  allocate(ar(dims(1),dims(2),dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_ish3_:', &
        ' Unable to allocate integer(kind=ish) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_ish(LID) = self%totmem_ish(LID) + product(dims)*self%bytes_ish 
  call self%update_total_memory_use_(product(dims)*self%bytes_ish, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ish
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_ish
  else 
    write (szstr,*) product(dims)*self%bytes_ish
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_ish
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_ish3_

!--------------------------------------------------------------------------
subroutine dealloc_ish3_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ish3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type ish3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable    :: ar(:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ish
    self%totmem_ish(LID) = self%totmem_ish(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ish3_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ish3_

!--------------------------------------------------------------------------
recursive subroutine alloc_ish4_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ish4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type ish4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable       :: ar(:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(4)
character(*),INTENT(IN)                          :: varname
integer(kind=ish), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(4)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(4) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ish
    self%totmem_ish(LID) = self%totmem_ish(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_ish4_:', &
      ' Unable to allocate integer(kind=ish) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_ish(LID) = self%totmem_ish(LID) + product(dims-startdims+1)*self%bytes_ish 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_ish, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_ish4_:', &
        ' Unable to allocate integer(kind=ish) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_ish(LID) = self%totmem_ish(LID) + product(dims)*self%bytes_ish 
  call self%update_total_memory_use_(product(dims)*self%bytes_ish, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ish
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_ish
  else 
    write (szstr,*) product(dims)*self%bytes_ish
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_ish
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_ish4_

!--------------------------------------------------------------------------
subroutine dealloc_ish4_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ish4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type ish4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable    :: ar(:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ish
    self%totmem_ish(LID) = self%totmem_ish(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ish4_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ish4_

!--------------------------------------------------------------------------
recursive subroutine alloc_ish5_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ish5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type ish5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(5)
character(*),INTENT(IN)                          :: varname
integer(kind=ish), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(5)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(5) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ish
    self%totmem_ish(LID) = self%totmem_ish(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4),startdims(5):dims(5)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_ish5_:', &
      ' Unable to allocate integer(kind=ish) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_ish(LID) = self%totmem_ish(LID) + product(dims-startdims+1)*self%bytes_ish 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_ish, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_ish5_:', &
        ' Unable to allocate integer(kind=ish) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_ish(LID) = self%totmem_ish(LID) + product(dims)*self%bytes_ish 
  call self%update_total_memory_use_(product(dims)*self%bytes_ish, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ish
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_ish
  else 
    write (szstr,*) product(dims)*self%bytes_ish
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_ish
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_ish5_

!--------------------------------------------------------------------------
subroutine dealloc_ish5_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ish5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type ish5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ish
    self%totmem_ish(LID) = self%totmem_ish(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ish5_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ish5_

!--------------------------------------------------------------------------
recursive subroutine alloc_ish6_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ish6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type ish6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(6)
character(*),INTENT(IN)                          :: varname
integer(kind=ish), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(6)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(6) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ish
    self%totmem_ish(LID) = self%totmem_ish(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4),&
              startdims(5):dims(5),startdims(6):dims(6)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_ish6_:', &
      ' Unable to allocate integer(kind=ish) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_ish(LID) = self%totmem_ish(LID) + product(dims-startdims+1)*self%bytes_ish 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_ish, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_ish6_:', &
        ' Unable to allocate integer(kind=ish) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_ish(LID) = self%totmem_ish(LID) + product(dims)*self%bytes_ish 
  call self%update_total_memory_use_(product(dims)*self%bytes_ish, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ish
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_ish
  else 
    write (szstr,*) product(dims)*self%bytes_ish
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_ish
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_ish6_

!--------------------------------------------------------------------------
subroutine dealloc_ish6_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ish6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type ish6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ish
    self%totmem_ish(LID) = self%totmem_ish(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ish6_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ish6_

!--------------------------------------------------------------------------
recursive subroutine alloc_irg1_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_irg1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type irg1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable       :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(1)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(1) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_irg
    self%totmem_irg(LID) = self%totmem_irg(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_irg1_:', &
      ' Unable to allocate integer(kind=irg) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_irg(LID) = self%totmem_irg(LID) + product(dims-startdims+1)*self%bytes_irg 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_irg, LID)
else
  allocate(ar(dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_irg1_:', &
        ' Unable to allocate integer(kind=irg) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_irg(LID) = self%totmem_irg(LID) + product(dims)*self%bytes_irg 
  call self%update_total_memory_use_(product(dims)*self%bytes_irg, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_irg
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_irg
  else 
    write (szstr,*) product(dims)*self%bytes_irg
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_irg
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_irg1_

!--------------------------------------------------------------------------
subroutine dealloc_irg1_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_irg1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type irg1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable    :: ar(:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_irg
    self%totmem_irg(LID) = self%totmem_irg(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_irg1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_irg1_

!--------------------------------------------------------------------------
recursive subroutine alloc_irg2_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_irg2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type irg2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable       :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(2)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(2) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_irg
    self%totmem_irg(LID) = self%totmem_irg(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_irg2_:', &
      ' Unable to allocate integer(kind=irg) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_irg(LID) = self%totmem_irg(LID) + product(dims-startdims+1)*self%bytes_irg 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_irg, LID)
else
  allocate(ar(dims(1),dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_irg2_:', &
        ' Unable to allocate integer(kind=irg) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_irg(LID) = self%totmem_irg(LID) + product(dims)*self%bytes_irg 
  call self%update_total_memory_use_(product(dims)*self%bytes_irg, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_irg
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_irg
  else 
    write (szstr,*) product(dims)*self%bytes_irg
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_irg
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_irg2_

!--------------------------------------------------------------------------
subroutine dealloc_irg2_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_irg2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type irg2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable    :: ar(:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_irg
    self%totmem_irg(LID) = self%totmem_irg(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_irg2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_irg2_

!--------------------------------------------------------------------------
recursive subroutine alloc_irg3_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_irg3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type irg3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable       :: ar(:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(3)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(3)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(3) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_irg
    self%totmem_irg(LID) = self%totmem_irg(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_irg3_:', &
      ' Unable to allocate integer(kind=irg) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_irg(LID) = self%totmem_irg(LID) + product(dims-startdims+1)*self%bytes_irg 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_irg, LID)
else
  allocate(ar(dims(1),dims(2),dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_irg3_:', &
        ' Unable to allocate integer(kind=irg) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_irg(LID) = self%totmem_irg(LID) + product(dims)*self%bytes_irg 
  call self%update_total_memory_use_(product(dims)*self%bytes_irg, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_irg
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_irg
  else 
    write (szstr,*) product(dims)*self%bytes_irg
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_irg
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_irg3_

!--------------------------------------------------------------------------
subroutine dealloc_irg3_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_irg3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type irg3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable    :: ar(:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_irg
    self%totmem_irg(LID) = self%totmem_irg(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_irg3_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_irg3_

!--------------------------------------------------------------------------
recursive subroutine alloc_irg4_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_irg4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type irg4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable       :: ar(:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(4)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(4)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(4) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_irg
    self%totmem_irg(LID) = self%totmem_irg(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_irg4_:', &
      ' Unable to allocate integer(kind=irg) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_irg(LID) = self%totmem_irg(LID) + product(dims-startdims+1)*self%bytes_irg 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_irg, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_irg4_:', &
        ' Unable to allocate integer(kind=irg) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_irg(LID) = self%totmem_irg(LID) + product(dims)*self%bytes_irg 
  call self%update_total_memory_use_(product(dims)*self%bytes_irg, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_irg
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_irg
  else 
    write (szstr,*) product(dims)*self%bytes_irg
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_irg
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_irg4_

!--------------------------------------------------------------------------
subroutine dealloc_irg4_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_irg4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type irg4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable    :: ar(:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_irg
    self%totmem_irg(LID) = self%totmem_irg(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_irg4_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_irg4_

!--------------------------------------------------------------------------
recursive subroutine alloc_irg5_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_irg5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type irg5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(5)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(5)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(5) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_irg
    self%totmem_irg(LID) = self%totmem_irg(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4),startdims(5):dims(5)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_irg5_:', &
      ' Unable to allocate integer(kind=irg) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_irg(LID) = self%totmem_irg(LID) + product(dims-startdims+1)*self%bytes_irg 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_irg, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_irg5_:', &
        ' Unable to allocate integer(kind=irg) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_irg(LID) = self%totmem_irg(LID) + product(dims)*self%bytes_irg 
  call self%update_total_memory_use_(product(dims)*self%bytes_irg, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_irg
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_irg
  else 
    write (szstr,*) product(dims)*self%bytes_irg
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_irg
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_irg5_

!--------------------------------------------------------------------------
subroutine dealloc_irg5_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_irg5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type irg5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_irg
    self%totmem_irg(LID) = self%totmem_irg(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_irg5_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_irg5_

!--------------------------------------------------------------------------
recursive subroutine alloc_irg6_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_irg6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type irg6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(6)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(6)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(6) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_irg
    self%totmem_irg(LID) = self%totmem_irg(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4),& 
              startdims(5):dims(5),startdims(6):dims(6)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_irg6_:', &
      ' Unable to allocate integer(kind=irg) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_irg(LID) = self%totmem_irg(LID) + product(dims-startdims+1)*self%bytes_irg 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_irg, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_irg6_:', &
        ' Unable to allocate integer(kind=irg) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_irg(LID) = self%totmem_irg(LID) + product(dims)*self%bytes_irg 
  call self%update_total_memory_use_(product(dims)*self%bytes_irg, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_irg
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_irg
  else 
    write (szstr,*) product(dims)*self%bytes_irg
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_irg
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_irg6_

!--------------------------------------------------------------------------
subroutine dealloc_irg6_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_irg6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type irg6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_irg
    self%totmem_irg(LID) = self%totmem_irg(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_irg6_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_irg6_

!--------------------------------------------------------------------------
recursive subroutine alloc_ill1_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ill1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type ill1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable       :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
integer(kind=ill), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(1)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(1) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ill
    self%totmem_ill(LID) = self%totmem_ill(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_ill1_:', &
      ' Unable to allocate integer(kind=ill) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_ill(LID) = self%totmem_ill(LID) + product(dims-startdims+1)*self%bytes_ill 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_ill, LID)
else
  allocate(ar(dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_ill1_:', &
        ' Unable to allocate integer(kind=ill) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_ill(LID) = self%totmem_ill(LID) + product(dims)*self%bytes_ill 
  call self%update_total_memory_use_(product(dims)*self%bytes_ill, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ill
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_ill
  else 
    write (szstr,*) product(dims)*self%bytes_ill
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_ill
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_ill1_

!--------------------------------------------------------------------------
subroutine dealloc_ill1_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ill1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type ill1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable    :: ar(:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ill
    self%totmem_ill(LID) = self%totmem_ill(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ill1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ill1_

!--------------------------------------------------------------------------
recursive subroutine alloc_ill2_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ill2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type ill2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable       :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
integer(kind=ill), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(2)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(2) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ill
    self%totmem_ill(LID) = self%totmem_ill(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_ill2_:', &
      ' Unable to allocate integer(kind=ill) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_ill(LID) = self%totmem_ill(LID) + product(dims-startdims+1)*self%bytes_ill 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_ill, LID)
else
  allocate(ar(dims(1),dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_ill2_:', &
        ' Unable to allocate integer(kind=ill) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_ill(LID) = self%totmem_ill(LID) + product(dims)*self%bytes_ill 
  call self%update_total_memory_use_(product(dims)*self%bytes_ill, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ill
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_ill
  else 
    write (szstr,*) product(dims)*self%bytes_ill
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_ill
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_ill2_

!--------------------------------------------------------------------------
subroutine dealloc_ill2_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ill2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type ill2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable    :: ar(:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ill
    self%totmem_ill(LID) = self%totmem_ill(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ill2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ill2_

!--------------------------------------------------------------------------
recursive subroutine alloc_ill3_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ill3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type ill3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable       :: ar(:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(3)
character(*),INTENT(IN)                          :: varname
integer(kind=ill), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(3)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(3) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ill
    self%totmem_ill(LID) = self%totmem_ill(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_ill3_:', &
      ' Unable to allocate integer(kind=ill) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_ill(LID) = self%totmem_ill(LID) + product(dims-startdims+1)*self%bytes_ill 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_ill, LID)
else
  allocate(ar(dims(1),dims(2),dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_ill3_:', &
        ' Unable to allocate integer(kind=ill) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_ill(LID) = self%totmem_ill(LID) + product(dims)*self%bytes_ill 
  call self%update_total_memory_use_(product(dims)*self%bytes_ill, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ill
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_ill
  else 
    write (szstr,*) product(dims)*self%bytes_ill
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_ill
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_ill3_

!--------------------------------------------------------------------------
subroutine dealloc_ill3_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ill3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type ill3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable    :: ar(:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ill
    self%totmem_ill(LID) = self%totmem_ill(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ill3_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ill3_

!--------------------------------------------------------------------------
recursive subroutine alloc_ill4_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ill4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type ill4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable       :: ar(:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(4)
character(*),INTENT(IN)                          :: varname
integer(kind=ill), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(4)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(4) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ill
    self%totmem_ill(LID) = self%totmem_ill(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_ill4_:', &
      ' Unable to allocate integer(kind=ill) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_ill(LID) = self%totmem_ill(LID) + product(dims-startdims+1)*self%bytes_ill 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_ill, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_ill4_:', &
        ' Unable to allocate integer(kind=ill) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_ill(LID) = self%totmem_ill(LID) + product(dims)*self%bytes_ill 
  call self%update_total_memory_use_(product(dims)*self%bytes_ill, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ill
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_ill
  else 
    write (szstr,*) product(dims)*self%bytes_ill
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_ill
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_ill4_

!--------------------------------------------------------------------------
subroutine dealloc_ill4_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ill4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type ill4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable    :: ar(:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ill
    self%totmem_ill(LID) = self%totmem_ill(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ill4_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ill4_

!--------------------------------------------------------------------------
recursive subroutine alloc_ill5_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ill5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type ill5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(5)
character(*),INTENT(IN)                          :: varname
integer(kind=ill), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(5)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(5) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ill
    self%totmem_ill(LID) = self%totmem_ill(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4),startdims(5):dims(5)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_ill5_:', &
      ' Unable to allocate integer(kind=ill) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_ill(LID) = self%totmem_ill(LID) + product(dims-startdims+1)*self%bytes_ill 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_ill, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_ill5_:', &
        ' Unable to allocate integer(kind=ill) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_ill(LID) = self%totmem_ill(LID) + product(dims)*self%bytes_ill 
  call self%update_total_memory_use_(product(dims)*self%bytes_ill, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ill
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_ill
  else 
    write (szstr,*) product(dims)*self%bytes_ill
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_ill
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_ill5_

!--------------------------------------------------------------------------
subroutine dealloc_ill5_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ill5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type ill5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ill
    self%totmem_ill(LID) = self%totmem_ill(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ill5_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ill5_

!--------------------------------------------------------------------------
recursive subroutine alloc_ill6_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ill6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type ill6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(6)
character(*),INTENT(IN)                          :: varname
integer(kind=ill), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(6)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(6) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ill
    self%totmem_ill(LID) = self%totmem_ill(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4), &
              startdims(5):dims(5),startdims(6):dims(6)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_ill6_:', &
      ' Unable to allocate integer(kind=ill) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_ill(LID) = self%totmem_ill(LID) + product(dims-startdims+1)*self%bytes_ill 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_ill, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_ill6_:', &
        ' Unable to allocate integer(kind=ill) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_ill(LID) = self%totmem_ill(LID) + product(dims)*self%bytes_ill 
  call self%update_total_memory_use_(product(dims)*self%bytes_ill, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ill
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_ill
  else 
    write (szstr,*) product(dims)*self%bytes_ill
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0_ill
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_ill6_

!--------------------------------------------------------------------------
subroutine dealloc_ill6_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ill6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type ill6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ill
    self%totmem_ill(LID) = self%totmem_ill(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ill6_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ill6_

!--------------------------------------------------------------------------
recursive subroutine alloc_sgl1_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_sgl1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type sgl1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable       :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
real(kind=sgl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(1)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(1) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_sgl
    self%totmem_sgl(LID) = self%totmem_sgl(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_sgl1_:', &
      ' Unable to allocate real(kind=sgl) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_sgl(LID) = self%totmem_sgl(LID) + product(dims-startdims+1)*self%bytes_sgl 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_sgl, LID)
else
  allocate(ar(dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_sgl1_:', &
        ' Unable to allocate real(kind=sgl) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_sgl(LID) = self%totmem_sgl(LID) + product(dims)*self%bytes_sgl 
  call self%update_total_memory_use_(product(dims)*self%bytes_sgl, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.0_sgl
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_sgl
  else 
    write (szstr,*) product(dims)*self%bytes_sgl
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0.0_sgl
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_sgl1_

!--------------------------------------------------------------------------
subroutine dealloc_sgl1_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_sgl1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type sgl1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable    :: ar(:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_sgl
    self%totmem_sgl(LID) = self%totmem_sgl(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_sgl1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_sgl1_

!--------------------------------------------------------------------------
recursive subroutine alloc_sgl2_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_sgl2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type sgl2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
real(kind=sgl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(2)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(2) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_sgl
    self%totmem_sgl(LID) = self%totmem_sgl(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_sgl2_:', &
      ' Unable to allocate real(kind=sgl) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_sgl(LID) = self%totmem_sgl(LID) + product(dims-startdims+1)*self%bytes_sgl 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_sgl, LID)
else
  allocate(ar(dims(1),dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_sgl2_:', &
        ' Unable to allocate real(kind=sgl) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_sgl(LID) = self%totmem_sgl(LID) + product(dims)*self%bytes_sgl 
  call self%update_total_memory_use_(product(dims)*self%bytes_sgl, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.0_sgl
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_sgl
  else 
    write (szstr,*) product(dims)*self%bytes_sgl
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0.0_sgl
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_sgl2_

!--------------------------------------------------------------------------
subroutine dealloc_sgl2_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_sgl2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type sgl2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_sgl
    self%totmem_sgl(LID) = self%totmem_sgl(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_sgl2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_sgl2_

!--------------------------------------------------------------------------
recursive subroutine alloc_sgl3_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_sgl3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type sgl3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(3)
character(*),INTENT(IN)                          :: varname
real(kind=sgl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(3)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(3) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_sgl
    self%totmem_sgl(LID) = self%totmem_sgl(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_sgl3_:', &
      ' Unable to allocate real(kind=sgl) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_sgl(LID) = self%totmem_sgl(LID) + product(dims-startdims+1)*self%bytes_sgl 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_sgl, LID)
else
  allocate(ar(dims(1),dims(2),dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_sgl3_:', &
        ' Unable to allocate real(kind=sgl) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_sgl(LID) = self%totmem_sgl(LID) + product(dims)*self%bytes_sgl 
  call self%update_total_memory_use_(product(dims)*self%bytes_sgl, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.0_sgl
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_sgl
  else 
    write (szstr,*) product(dims)*self%bytes_sgl
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0.0_sgl
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_sgl3_

!--------------------------------------------------------------------------
subroutine dealloc_sgl3_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_sgl3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type sgl3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_sgl
    self%totmem_sgl(LID) = self%totmem_sgl(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_sgl3_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_sgl3_

!--------------------------------------------------------------------------
recursive subroutine alloc_sgl4_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_sgl4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type sgl4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(4)
character(*),INTENT(IN)                          :: varname
real(kind=sgl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(4)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(4) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_sgl
    self%totmem_sgl(LID) = self%totmem_sgl(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_sgl4_:', &
      ' Unable to allocate real(kind=sgl) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_sgl(LID) = self%totmem_sgl(LID) + product(dims-startdims+1)*self%bytes_sgl 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_sgl, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_sgl4_:', &
        ' Unable to allocate real(kind=sgl) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_sgl(LID) = self%totmem_sgl(LID) + product(dims)*self%bytes_sgl 
  call self%update_total_memory_use_(product(dims)*self%bytes_sgl, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.0_sgl
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_sgl
  else 
    write (szstr,*) product(dims)*self%bytes_sgl
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0.0_sgl
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_sgl4_

!--------------------------------------------------------------------------
subroutine dealloc_sgl4_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_sgl4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type sgl4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_sgl
    self%totmem_sgl(LID) = self%totmem_sgl(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_sgl4_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_sgl4_

!--------------------------------------------------------------------------
recursive subroutine alloc_sgl5_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_sgl5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type sgl5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(5)
character(*),INTENT(IN)                          :: varname
real(kind=sgl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(5)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(5) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_sgl
    self%totmem_sgl(LID) = self%totmem_sgl(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4),startdims(5):dims(5)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_sgl5_:', &
      ' Unable to allocate real(kind=sgl) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_sgl(LID) = self%totmem_sgl(LID) + product(dims-startdims+1)*self%bytes_sgl 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_sgl, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_sgl5_:', &
        ' Unable to allocate real(kind=sgl) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_sgl(LID) = self%totmem_sgl(LID) + product(dims)*self%bytes_sgl 
  call self%update_total_memory_use_(product(dims)*self%bytes_sgl, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.0_sgl
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_sgl
  else 
    write (szstr,*) product(dims)*self%bytes_sgl
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0.0_sgl
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_sgl5_

!--------------------------------------------------------------------------
subroutine dealloc_sgl5_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_sgl5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type sgl5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_sgl
    self%totmem_sgl(LID) = self%totmem_sgl(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_sgl5_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_sgl5_

!--------------------------------------------------------------------------
recursive subroutine alloc_sgl6_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_sgl6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type sgl6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(6)
character(*),INTENT(IN)                          :: varname
real(kind=sgl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(6)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(6) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_sgl
    self%totmem_sgl(LID) = self%totmem_sgl(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4), &
              startdims(5):dims(5),startdims(6):dims(6)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_sgl6_:', &
      ' Unable to allocate real(kind=sgl) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_sgl(LID) = self%totmem_sgl(LID) + product(dims-startdims+1)*self%bytes_sgl 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_sgl, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_sgl6_:', &
        ' Unable to allocate real(kind=sgl) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_sgl(LID) = self%totmem_sgl(LID) + product(dims)*self%bytes_sgl 
  call self%update_total_memory_use_(product(dims)*self%bytes_sgl, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.0_sgl
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_sgl
  else 
    write (szstr,*) product(dims)*self%bytes_sgl
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0.0_sgl
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_sgl6_

!--------------------------------------------------------------------------
subroutine dealloc_sgl6_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_sgl6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type sgl6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_sgl
    self%totmem_sgl(LID) = self%totmem_sgl(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_sgl6_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_sgl6_

!--------------------------------------------------------------------------
recursive subroutine alloc_dbl1_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_dbl1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type dbl1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable       :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
real(kind=dbl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(1)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(1) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_dbl
    self%totmem_dbl(LID) = self%totmem_dbl(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_dbl1_:', &
      ' Unable to allocate real(kind=dbl) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_dbl(LID) = self%totmem_dbl(LID) + product(dims-startdims+1)*self%bytes_dbl 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_dbl, LID)
else
  allocate(ar(dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_dbl1_:', &
        ' Unable to allocate real(kind=dbl) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_dbl(LID) = self%totmem_dbl(LID) + product(dims)*self%bytes_dbl 
  call self%update_total_memory_use_(product(dims)*self%bytes_dbl, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.D0
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_dbl
  else 
    write (szstr,*) product(dims)*self%bytes_dbl
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0.D0
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_dbl1_

!--------------------------------------------------------------------------
subroutine dealloc_dbl1_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_dbl1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type dbl1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable    :: ar(:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_dbl
    self%totmem_dbl(LID) = self%totmem_dbl(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_dbl1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_dbl1_

!--------------------------------------------------------------------------
recursive subroutine alloc_dbl2_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_dbl2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type dbl2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
real(kind=dbl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(2)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(2) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_dbl
    self%totmem_dbl(LID) = self%totmem_dbl(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_dbl2_:', &
      ' Unable to allocate real(kind=dbl) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_dbl(LID) = self%totmem_dbl(LID) + product(dims-startdims+1)*self%bytes_dbl 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_dbl, LID)
else
  allocate(ar(dims(1),dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_dbl2_:', &
        ' Unable to allocate real(kind=dbl) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_dbl(LID) = self%totmem_dbl(LID) + product(dims)*self%bytes_dbl 
  call self%update_total_memory_use_(product(dims)*self%bytes_dbl, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.D0
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_dbl
  else 
    write (szstr,*) product(dims)*self%bytes_dbl
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0.D0
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_dbl2_

!--------------------------------------------------------------------------
subroutine dealloc_dbl2_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_dbl2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type dbl2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_dbl
    self%totmem_dbl(LID) = self%totmem_dbl(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_dbl2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_dbl2_

!--------------------------------------------------------------------------
recursive subroutine alloc_dbl3_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_dbl3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type dbl3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(3)
character(*),INTENT(IN)                          :: varname
real(kind=dbl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(3)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(3) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_dbl
    self%totmem_dbl(LID) = self%totmem_dbl(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_dbl3_:', &
      ' Unable to allocate real(kind=dbl) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_dbl(LID) = self%totmem_dbl(LID) + product(dims-startdims+1)*self%bytes_dbl 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_dbl, LID)
else
  allocate(ar(dims(1),dims(2),dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_dbl3_:', &
        ' Unable to allocate real(kind=dbl) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_dbl(LID) = self%totmem_dbl(LID) + product(dims)*self%bytes_dbl 
  call self%update_total_memory_use_(product(dims)*self%bytes_dbl, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.D0
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_dbl
  else 
    write (szstr,*) product(dims)*self%bytes_dbl
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0.D0
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_dbl3_

!--------------------------------------------------------------------------
subroutine dealloc_dbl3_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_dbl3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type dbl3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_dbl
    self%totmem_dbl(LID) = self%totmem_dbl(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_dbl3_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_dbl3_

!--------------------------------------------------------------------------
recursive subroutine alloc_dbl4_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_dbl4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type dbl4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(4)
character(*),INTENT(IN)                          :: varname
real(kind=dbl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(4)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(4) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_dbl
    self%totmem_dbl(LID) = self%totmem_dbl(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_dbl4_:', &
      ' Unable to allocate real(kind=dbl) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_dbl(LID) = self%totmem_dbl(LID) + product(dims-startdims+1)*self%bytes_dbl 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_dbl, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_dbl4_:', &
        ' Unable to allocate real(kind=dbl) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_dbl(LID) = self%totmem_dbl(LID) + product(dims)*self%bytes_dbl 
  call self%update_total_memory_use_(product(dims)*self%bytes_dbl, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.D0
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_dbl
  else 
    write (szstr,*) product(dims)*self%bytes_dbl
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0.D0
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_dbl4_

!--------------------------------------------------------------------------
subroutine dealloc_dbl4_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_dbl4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type dbl4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_dbl
    self%totmem_dbl(LID) = self%totmem_dbl(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_dbl4_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_dbl4_

!--------------------------------------------------------------------------
recursive subroutine alloc_dbl5_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_dbl5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type dbl5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(5)
character(*),INTENT(IN)                          :: varname
real(kind=dbl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(5)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(5) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_dbl
    self%totmem_dbl(LID) = self%totmem_dbl(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4),startdims(5):dims(5)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_dbl5_:', &
      ' Unable to allocate real(kind=dbl) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_dbl(LID) = self%totmem_dbl(LID) + product(dims-startdims+1)*self%bytes_dbl 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_dbl, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_dbl5_:', &
        ' Unable to allocate real(kind=dbl) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_dbl(LID) = self%totmem_dbl(LID) + product(dims)*self%bytes_dbl 
  call self%update_total_memory_use_(product(dims)*self%bytes_dbl, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.D0
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_dbl
  else 
    write (szstr,*) product(dims)*self%bytes_dbl
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0.D0
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_dbl5_

!--------------------------------------------------------------------------
subroutine dealloc_dbl5_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_dbl5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type dbl5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_dbl
    self%totmem_dbl(LID) = self%totmem_dbl(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_dbl5_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_dbl5_

!--------------------------------------------------------------------------
recursive subroutine alloc_dbl6_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_dbl6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type dbl6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(6)
character(*),INTENT(IN)                          :: varname
real(kind=dbl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(6)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(6) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_dbl
    self%totmem_dbl(LID) = self%totmem_dbl(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4), &
              startdims(5):dims(5),startdims(6):dims(6)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_dbl6_:', &
      ' Unable to allocate real(kind=dbl) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_dbl(LID) = self%totmem_dbl(LID) + product(dims-startdims+1)*self%bytes_dbl 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_dbl, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_dbl6_:', &
        ' Unable to allocate real(kind=dbl) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_dbl(LID) = self%totmem_dbl(LID) + product(dims)*self%bytes_dbl 
  call self%update_total_memory_use_(product(dims)*self%bytes_dbl, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.D0
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_dbl
  else 
    write (szstr,*) product(dims)*self%bytes_dbl
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) 0.D0
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_dbl6_

!--------------------------------------------------------------------------
subroutine dealloc_dbl6_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_dbl6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type dbl6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_dbl
    self%totmem_dbl(LID) = self%totmem_dbl(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_dbl6_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_dbl6_

!--------------------------------------------------------------------------
recursive subroutine alloc_cmplx1_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplx1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type cmplx1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable       :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
complex(kind=sgl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(1)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(1) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplx
    self%totmem_cmplx(LID) = self%totmem_cmplx(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_cmplx1_:', &
      ' Unable to allocate complex(kind=cmplx) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_cmplx(LID) = self%totmem_cmplx(LID) + product(dims-startdims+1)*self%bytes_cmplx 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_cmplx, LID)
else
  allocate(ar(dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_cmplx1_:', &
        ' Unable to allocate complex(kind=cmplx) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_cmplx(LID) = self%totmem_cmplx(LID) + product(dims)*self%bytes_cmplx 
  call self%update_total_memory_use_(product(dims)*self%bytes_cmplx, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.0_sgl,0.0_sgl)
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_cmplx
  else 
    write (szstr,*) product(dims)*self%bytes_cmplx
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) (0.0_sgl,0.0_sgl)
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_cmplx1_

!--------------------------------------------------------------------------
subroutine dealloc_cmplx1_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplx1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type cmplx1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable    :: ar(:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplx
    self%totmem_cmplx(LID) = self%totmem_cmplx(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplx1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplx1_

!--------------------------------------------------------------------------
recursive subroutine alloc_cmplx2_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplx2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type cmplx2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
complex(kind=sgl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(2)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(2) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplx
    self%totmem_cmplx(LID) = self%totmem_cmplx(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_cmplx2_:', &
      ' Unable to allocate complex(kind=cmplx) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_cmplx(LID) = self%totmem_cmplx(LID) + product(dims-startdims+1)*self%bytes_cmplx 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_cmplx, LID)
else
  allocate(ar(dims(1),dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_cmplx2_:', &
        ' Unable to allocate complex(kind=cmplx) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_cmplx(LID) = self%totmem_cmplx(LID) + product(dims)*self%bytes_cmplx 
  call self%update_total_memory_use_(product(dims)*self%bytes_cmplx, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.0_sgl,0.0_sgl)
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_cmplx
  else 
    write (szstr,*) product(dims)*self%bytes_cmplx
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) (0.0_sgl,0.0_sgl)
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_cmplx2_

!--------------------------------------------------------------------------
subroutine dealloc_cmplx2_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplx2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type cmplx2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplx
    self%totmem_cmplx(LID) = self%totmem_cmplx(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplx2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplx2_

!--------------------------------------------------------------------------
recursive subroutine alloc_cmplx3_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplx3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type cmplx3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(3)
character(*),INTENT(IN)                          :: varname
complex(kind=sgl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(3)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(3) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplx
    self%totmem_cmplx(LID) = self%totmem_cmplx(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_cmplx3_:', &
      ' Unable to allocate complex(kind=cmplx) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_cmplx(LID) = self%totmem_cmplx(LID) + product(dims-startdims+1)*self%bytes_cmplx 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_cmplx, LID)
else
  allocate(ar(dims(1),dims(2),dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_cmplx3_:', &
        ' Unable to allocate complex(kind=cmplx) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_cmplx(LID) = self%totmem_cmplx(LID) + product(dims)*self%bytes_cmplx 
  call self%update_total_memory_use_(product(dims)*self%bytes_cmplx, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.0_sgl,0.0_sgl)
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_cmplx
  else 
    write (szstr,*) product(dims)*self%bytes_cmplx
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) (0.0_sgl,0.0_sgl)
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_cmplx3_

!--------------------------------------------------------------------------
subroutine dealloc_cmplx3_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplx3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type cmplx3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplx
    self%totmem_cmplx(LID) = self%totmem_cmplx(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplx3_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplx3_

!--------------------------------------------------------------------------
recursive subroutine alloc_cmplx4_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplx4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type cmplx4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(4)
character(*),INTENT(IN)                          :: varname
complex(kind=sgl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(4)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(4) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplx
    self%totmem_cmplx(LID) = self%totmem_cmplx(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_cmplx4_:', &
      ' Unable to allocate complex(kind=cmplx) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_cmplx(LID) = self%totmem_cmplx(LID) + product(dims-startdims+1)*self%bytes_cmplx 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_cmplx, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_cmplx4_:', &
        ' Unable to allocate complex(kind=cmplx) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_cmplx(LID) = self%totmem_cmplx(LID) + product(dims)*self%bytes_cmplx 
  call self%update_total_memory_use_(product(dims)*self%bytes_cmplx, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.0_sgl,0.0_sgl)
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_cmplx
  else 
    write (szstr,*) product(dims)*self%bytes_cmplx
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) (0.0_sgl,0.0_sgl)
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_cmplx4_

!--------------------------------------------------------------------------
subroutine dealloc_cmplx4_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplx4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type cmplx4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplx
    self%totmem_cmplx(LID) = self%totmem_cmplx(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplx4_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplx4_

!--------------------------------------------------------------------------
recursive subroutine alloc_cmplx5_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplx5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type cmplx5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(5)
character(*),INTENT(IN)                          :: varname
complex(kind=sgl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(5)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(5) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplx
    self%totmem_cmplx(LID) = self%totmem_cmplx(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4),startdims(5):dims(5)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_cmplx5_:', &
      ' Unable to allocate complex(kind=cmplx) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_cmplx(LID) = self%totmem_cmplx(LID) + product(dims-startdims+1)*self%bytes_cmplx 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_cmplx, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_cmplx5_:', &
        ' Unable to allocate complex(kind=cmplx) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_cmplx(LID) = self%totmem_cmplx(LID) + product(dims)*self%bytes_cmplx 
  call self%update_total_memory_use_(product(dims)*self%bytes_cmplx, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.0_sgl,0.0_sgl)
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_cmplx
  else 
    write (szstr,*) product(dims)*self%bytes_cmplx
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) (0.0_sgl,0.0_sgl)
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_cmplx5_

!--------------------------------------------------------------------------
subroutine dealloc_cmplx5_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplx5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type cmplx5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplx
    self%totmem_cmplx(LID) = self%totmem_cmplx(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplx5_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplx5_

!--------------------------------------------------------------------------
recursive subroutine alloc_cmplx6_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplx6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type cmplx6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(6)
character(*),INTENT(IN)                          :: varname
complex(kind=sgl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(6)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(6) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplx
    self%totmem_cmplx(LID) = self%totmem_cmplx(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4), &
              startdims(5):dims(5),startdims(6):dims(6)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_cmplx6_:', &
      ' Unable to allocate complex(kind=cmplx) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_cmplx(LID) = self%totmem_cmplx(LID) + product(dims-startdims+1)*self%bytes_cmplx 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_cmplx, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_cmplx6_:', &
        ' Unable to allocate complex(kind=cmplx) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_cmplx(LID) = self%totmem_cmplx(LID) + product(dims)*self%bytes_cmplx 
  call self%update_total_memory_use_(product(dims)*self%bytes_cmplx, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.0_sgl,0.0_sgl)
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_cmplx
  else 
    write (szstr,*) product(dims)*self%bytes_cmplx
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) (0.0_sgl,0.0_sgl)
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_cmplx6_

!--------------------------------------------------------------------------
subroutine dealloc_cmplx6_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplx6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type cmplx6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplx
    self%totmem_cmplx(LID) = self%totmem_cmplx(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplx6_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplx6_

!--------------------------------------------------------------------------
recursive subroutine alloc_cmplxd1_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplxd1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type cmplxd1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable       :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
complex(kind=dbl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(1)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(1) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplxd
    self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_cmplxd1_:', &
      ' Unable to allocate complex(kind=cmplxd) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) + product(dims-startdims+1)*self%bytes_cmplxd 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_cmplxd, LID)
else
  allocate(ar(dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_cmplxd1_:', &
        ' Unable to allocate complex(kind=cmplxd) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) + product(dims)*self%bytes_cmplxd 
  call self%update_total_memory_use_(product(dims)*self%bytes_cmplxd, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.D0,0.D0)
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_cmplxd
  else 
    write (szstr,*) product(dims)*self%bytes_cmplxd
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) (0.D0,0.D0)
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_cmplxd1_

!--------------------------------------------------------------------------
subroutine dealloc_cmplxd1_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplxd1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type cmplxd1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable    :: ar(:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplxd
    self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplxd1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplxd1_

!--------------------------------------------------------------------------
recursive subroutine alloc_cmplxd2_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplxd2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type cmplxd2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
complex(kind=dbl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(2)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(2) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplxd
    self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_cmplxd2_:', &
      ' Unable to allocate complex(kind=cmplxd) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) + product(dims-startdims+1)*self%bytes_cmplxd 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_cmplxd, LID)
else
  allocate(ar(dims(1),dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_cmplxd2_:', &
        ' Unable to allocate complex(kind=cmplxd) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) + product(dims)*self%bytes_cmplxd 
  call self%update_total_memory_use_(product(dims)*self%bytes_cmplxd, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.D0,0.D0)
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_cmplxd
  else 
    write (szstr,*) product(dims)*self%bytes_cmplxd
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) (0.D0,0.D0)
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_cmplxd2_

!--------------------------------------------------------------------------
subroutine dealloc_cmplxd2_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplxd2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type cmplxd2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplxd
    self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplxd2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplxd2_

!--------------------------------------------------------------------------
recursive subroutine alloc_cmplxd3_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplxd3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type cmplxd3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(3)
character(*),INTENT(IN)                          :: varname
complex(kind=dbl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(3)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(3) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplxd
    self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_cmplxd3_:', &
      ' Unable to allocate complex(kind=cmplxd) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) + product(dims-startdims+1)*self%bytes_cmplxd 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_cmplxd, LID)
else
  allocate(ar(dims(1),dims(2),dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_cmplxd3_:', &
        ' Unable to allocate complex(kind=cmplxd) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) + product(dims)*self%bytes_cmplxd 
  call self%update_total_memory_use_(product(dims)*self%bytes_cmplxd, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.D0,0.D0)
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_cmplxd
  else 
    write (szstr,*) product(dims)*self%bytes_cmplxd
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) (0.D0,0.D0)
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_cmplxd3_

!--------------------------------------------------------------------------
subroutine dealloc_cmplxd3_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplxd3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type cmplxd3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplxd
    self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplxd3_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplxd3_

!--------------------------------------------------------------------------
recursive subroutine alloc_cmplxd4_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplxd4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type cmplxd4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(4)
character(*),INTENT(IN)                          :: varname
complex(kind=dbl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(4)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(4) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplxd
    self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_cmplxd4_:', &
      ' Unable to allocate complex(kind=cmplxd) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) + product(dims-startdims+1)*self%bytes_cmplxd 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_cmplxd, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_cmplxd4_:', &
        ' Unable to allocate complex(kind=cmplxd) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) + product(dims)*self%bytes_cmplxd 
  call self%update_total_memory_use_(product(dims)*self%bytes_cmplxd, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.D0,0.D0)
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_cmplxd
  else 
    write (szstr,*) product(dims)*self%bytes_cmplxd
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) (0.D0,0.D0)
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_cmplxd4_

!--------------------------------------------------------------------------
subroutine dealloc_cmplxd4_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplxd4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type cmplxd4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplxd
    self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplxd4_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplxd4_

!--------------------------------------------------------------------------
recursive subroutine alloc_cmplxd5_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplxd5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type cmplxd5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(5)
character(*),INTENT(IN)                          :: varname
complex(kind=dbl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(5)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(5) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplxd
    self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4),startdims(5):dims(5)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_cmplxd5_:', &
      ' Unable to allocate complex(kind=cmplxd) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) + product(dims-startdims+1)*self%bytes_cmplxd 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_cmplxd, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_cmplxd5_:', &
        ' Unable to allocate complex(kind=cmplxd) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) + product(dims)*self%bytes_cmplxd 
  call self%update_total_memory_use_(product(dims)*self%bytes_cmplxd, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.D0,0.D0)
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_cmplxd
  else 
    write (szstr,*) product(dims)*self%bytes_cmplxd
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) (0.D0,0.D0)
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_cmplxd5_

!--------------------------------------------------------------------------
subroutine dealloc_cmplxd5_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplxd5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type cmplxd5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplxd
    self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplxd5_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplxd5_

!--------------------------------------------------------------------------
recursive subroutine alloc_cmplxd6_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplxd6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! allocate an array of type cmplxd6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(6)
character(*),INTENT(IN)                          :: varname
complex(kind=dbl), INTENT(IN), OPTIONAL             :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(6)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(6) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplxd
    self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3),startdims(4):dims(4),&
              startdims(5):dims(5),startdims(6):dims(6)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_cmplxd6_:', &
      ' Unable to allocate complex(kind=cmplxd) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) + product(dims-startdims+1)*self%bytes_cmplxd 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_cmplxd, LID)
else
  allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_cmplxd6_:', &
        ' Unable to allocate complex(kind=cmplxd) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) + product(dims)*self%bytes_cmplxd 
  call self%update_total_memory_use_(product(dims)*self%bytes_cmplxd, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.D0,0.D0)
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_cmplxd
  else 
    write (szstr,*) product(dims)*self%bytes_cmplxd
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) (0.D0,0.D0)
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_cmplxd6_

!--------------------------------------------------------------------------
subroutine dealloc_cmplxd6_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplxd6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20, modified for OpenMP use on 03/24/21
!!
!! deallocate an array of type cmplxd6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplxd
    self%totmem_cmplxd(LID) = self%totmem_cmplxd(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplxd6_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplxd6_

!--------------------------------------------------------------------------
subroutine alloc_logical1_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_logical1_
!! author: MDG
!! version: 1.0
!! date: 04/14/21
!!
!! allocate an array of type logical1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
logical, INTENT(INOUT), allocatable              :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
logical, INTENT(IN), OPTIONAL                    :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(1)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(1) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_logical
    self%totmem_logical(LID) = self%totmem_logical(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_logical1_:', &
      ' Unable to allocate logical(kind=logical) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_logical(LID) = self%totmem_logical(LID) + product(dims-startdims+1)*self%bytes_logical 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_logical, LID)
else
  allocate(ar(dims(1)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_logical1_:', &
        ' Unable to allocate logical(kind=logical) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_logical(LID) = self%totmem_logical(LID) + product(dims)*self%bytes_logical 
  call self%update_total_memory_use_(product(dims)*self%bytes_logical, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = .FALSE.
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_logical
  else 
    write (szstr,*) product(dims)*self%bytes_logical
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) .FALSE.
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_logical1_

!--------------------------------------------------------------------------
subroutine dealloc_logical1_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_logical1_
!! author: MDG
!! version: 1.0
!! date: 04/14/21
!!
!! deallocate an array of type logical1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
logical, INTENT(INOUT), allocatable              :: ar(:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_logical
    self%totmem_logical(LID) = self%totmem_logical(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_logical1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_logical1_

!--------------------------------------------------------------------------
subroutine alloc_logical2_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_logical2_
!! author: MDG
!! version: 1.0
!! date: 04/14/21
!!
!! allocate an array of type logical2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
logical, INTENT(INOUT), allocatable              :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
logical, INTENT(IN), OPTIONAL                    :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(2)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(2) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_logical
    self%totmem_logical(LID) = self%totmem_logical(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_logical2_:', &
      ' Unable to allocate logical(kind=logical) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_logical(LID) = self%totmem_logical(LID) + product(dims-startdims+1)*self%bytes_logical 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_logical, LID)
else
  allocate(ar(dims(1),dims(2)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_logical2_:', &
        ' Unable to allocate logical(kind=logical) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_logical(LID) = self%totmem_logical(LID) + product(dims)*self%bytes_logical 
  call self%update_total_memory_use_(product(dims)*self%bytes_logical, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = .FALSE.
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_logical
  else 
    write (szstr,*) product(dims)*self%bytes_logical
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) .FALSE.
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_logical2_

!--------------------------------------------------------------------------
subroutine dealloc_logical2_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_logical2_
!! author: MDG
!! version: 1.0
!! date: 04/14/21
!!
!! deallocate an array of type logical2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
logical, INTENT(INOUT), allocatable              :: ar(:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_logical
    self%totmem_logical(LID) = self%totmem_logical(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_logical2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_logical2_

!--------------------------------------------------------------------------
subroutine alloc_logical3_(self, ar, dims, varname, initval, TID, startdims)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_logical3_
!! author: MDG
!! version: 1.0
!! date: 04/14/21
!!
!! allocate an array of type logical3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
logical, INTENT(INOUT), allocatable              :: ar(:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(3)
character(*),INTENT(IN)                          :: varname
logical, INTENT(IN), OPTIONAL                    :: initval
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID
integer(kind=irg), INTENT(IN), OPTIONAL          :: startdims(3)

type(IO_T)                                       :: Message
character(fnlen)                                 :: estr, estr2, outstr, szstr, initstr
integer(kind=irg)                                :: i, sz, err, LID, szar(3) 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_logical
    self%totmem_logical(LID) = self%totmem_logical(LID) - sz
    deallocate(ar)
endif

! allocate the array 
! use the startdims array if present 
if (present(startdims)) then 
  allocate(ar(startdims(1):dims(1),startdims(2):dims(2),startdims(3):dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    estr2 = ' '
    write (estr,*) dims
    write (estr2,*) startdims
    outstr = trim(estr)//':'//trim(estr2)
    call Message%printError('mod_memory:alloc_logical3_:', &
      ' Unable to allocate logical(kind=logical) array '//trim(varname)//' of dimension '//trim(outstr))
  end if
  self%totmem_logical(LID) = self%totmem_logical(LID) + product(dims-startdims+1)*self%bytes_logical 
  call self%update_total_memory_use_(product(dims-startdims+1)*self%bytes_logical, LID)
else
  allocate(ar(dims(1),dims(2),dims(3)), stat = err)
  if (err.ne.0) then 
    estr = ' '
    write (estr,*) dims
    call Message%printError('mod_memory:alloc_logical3_:', &
        ' Unable to allocate logical(kind=logical) array'//trim(varname)//' of dimension '//trim(estr))
  end if
  self%totmem_logical(LID) = self%totmem_logical(LID) + product(dims)*self%bytes_logical 
  call self%update_total_memory_use_(product(dims)*self%bytes_logical, LID)
end if

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = .FALSE.
end if

if (self%verbose) then 
  if (present(startdims)) then 
    write (szstr,*) product(dims-startdims+1)*self%bytes_logical
  else 
    write (szstr,*) product(dims)*self%bytes_logical
  end if
  if (present(initval)) then 
    write (initstr,*) initval
  else 
    write (initstr,*) .FALSE.
  end if
  outstr = ' -> allocated array '//trim(varname)//' of size '//trim(szstr)//'; initialized to '//trim(initstr)
  call Message%printMessage(outstr)
end if

end subroutine alloc_logical3_

!--------------------------------------------------------------------------
subroutine dealloc_logical3_(self, ar, varname, TID)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_logical3_
!! author: MDG
!! version: 1.0
!! date: 04/14/21
!!
!! deallocate an array of type logical3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
logical, INTENT(INOUT), allocatable              :: ar(:,:,:)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL          :: TID

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err, LID, sz 

! set the local thread identifier
LID = 1
if (present(TID)) LID = TID

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_logical
    self%totmem_logical(LID) = self%totmem_logical(LID) - sz
    call self%update_total_memory_use_( -sz, LID )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_logical3_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_logical3_




end module mod_memory
