! ###################################################################
! Copyright (c) 2014-2020, Marc De Graef Research Group/Carnegie Mellon University
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
  !!--------------------------------------------------

use mod_kinds
use mod_global
use mod_io

IMPLICIT NONE

type, public :: memory_T
  private 
  ! we'll keep track of all memory that has been allocated by calls to 
  ! this class; so, we'll have memory usage reporting as well.
    integer(kind=irg) :: totmem_ish 
    integer(kind=irg) :: totmem_irg
    integer(kind=irg) :: totmem_ill 
    integer(kind=irg) :: totmem_sgl 
    integer(kind=irg) :: totmem_dbl 
    integer(kind=irg) :: totmem_cmplx
    integer(kind=irg) :: totmem_cmplxd
    integer(kind=irg) :: current_memory_allocated
    integer(kind=irg) :: bytes_ish = 2
    integer(kind=irg) :: bytes_irg = 4
    integer(kind=irg) :: bytes_ill = 8
    integer(kind=irg) :: bytes_sgl = 4
    integer(kind=irg) :: bytes_dbl = 8
    integer(kind=irg) :: bytes_cmplx = 8
    integer(kind=irg) :: bytes_cmplxd = 16
    logical, public   :: verbose

  contains
  private 

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

    procedure, pass(self) :: compute_size_
    procedure, pass(self) :: show_allocated_memory_use_ 
    procedure, pass(self) :: update_total_memory_use_ 
    procedure, pass(self) :: toggle_verbose_ 
    ! final :: memory_destructor 

    generic, public :: alloc1 => alloc_ish1_, alloc_irg1_, alloc_ill1_, &
                                 alloc_sgl1_, alloc_dbl1_, alloc_cmplx1_, &
                                 alloc_cmplxd1_

    generic, public :: alloc2 => alloc_ish2_, alloc_irg2_, alloc_ill2_, &
                                 alloc_sgl2_, alloc_dbl2_, alloc_cmplx2_, &
                                 alloc_cmplxd2_

    generic, public :: alloc3 => alloc_ish3_, alloc_irg3_, alloc_ill3_, &
                                 alloc_sgl3_, alloc_dbl3_, alloc_cmplx3_, &
                                 alloc_cmplxd3_

    generic, public :: alloc4 => alloc_ish4_, alloc_irg4_, alloc_ill4_, &
                                 alloc_sgl4_, alloc_dbl4_, alloc_cmplx4_, &
                                 alloc_cmplxd4_

    generic, public :: alloc5 => alloc_ish5_, alloc_irg5_, alloc_ill5_, &
                                 alloc_sgl5_, alloc_dbl5_, alloc_cmplx5_, &
                                 alloc_cmplxd5_

    generic, public :: alloc6 => alloc_ish6_, alloc_irg6_, alloc_ill6_, &
                                 alloc_sgl6_, alloc_dbl6_, alloc_cmplx6_, &
                                 alloc_cmplxd6_

    generic, public :: dealloc1 => dealloc_ish1_, dealloc_irg1_, dealloc_ill1_, &
                                   dealloc_sgl1_, dealloc_dbl1_, dealloc_cmplx1_, &
                                   dealloc_cmplxd1_

    generic, public :: dealloc2 => dealloc_ish2_, dealloc_irg2_, dealloc_ill2_, &
                                    dealloc_sgl2_, dealloc_dbl2_, dealloc_cmplx2_, &
                                    dealloc_cmplxd2_

    generic, public :: dealloc3 => dealloc_ish3_, dealloc_irg3_, dealloc_ill3_, &
                                    dealloc_sgl3_, dealloc_dbl3_, dealloc_cmplx3_, &
                                    dealloc_cmplxd3_

    generic, public :: dealloc4 => dealloc_ish4_, dealloc_irg4_, dealloc_ill4_, &
                                    dealloc_sgl4_, dealloc_dbl4_, dealloc_cmplx4_, &
                                    dealloc_cmplxd4_

    generic, public :: dealloc5 => dealloc_ish5_, dealloc_irg5_, dealloc_ill5_, &
                                    dealloc_sgl5_, dealloc_dbl5_, dealloc_cmplx5_, &
                                    dealloc_cmplxd5_

    generic, public :: dealloc6 => dealloc_ish6_, dealloc_irg6_, dealloc_ill6_, &
                                    dealloc_sgl6_, dealloc_dbl6_, dealloc_cmplx6_, &
                                    dealloc_cmplxd6_

    generic, public :: getsize => compute_size_
    generic, public :: allocated_memory_use => show_allocated_memory_use_

    generic, public :: toggle_verbose => toggle_verbose_
end type memory_T

! constructor for this class
interface memory_T 
  module procedure memory_constructor 
end interface memory_T 

contains

!--------------------------------------------------------------------------
type(memory_T) function memory_constructor( )  result(mem) 
!DEC$ ATTRIBUTES DLLEXPORT :: memory_constructor
!! author: MDG
!! version: 1.0
!! date: 12/11/20
!!
!! constructor for the memory_T Class

IMPLICIT NONE

! set all the memory usage counters to zero
mem%totmem_ish  = 0
mem%totmem_irg = 0
mem%totmem_ill  = 0
mem%totmem_sgl  = 0
mem%totmem_dbl  = 0
mem%totmem_cmplx  = 0
mem%totmem_cmplxd  = 0
mem%current_memory_allocated = 0
mem%verbose = .FALSE.

end function memory_constructor

!--------------------------------------------------------------------------
subroutine show_allocated_memory_use_( self )
!DEC$ ATTRIBUTES DLLEXPORT :: show_allocated_memory_use_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! print out the current memory usage 

IMPLICIT NONE

class(memory_T), INTENT(INOUT)      :: self
type(IO_T)                          :: Message 

integer(kind=irg)                   :: io_int(1) 

call Message%printMessage(' Allocated Memory Status (in bytes) : ',"(/A/)")
io_int = self%totmem_ish
call Message%WriteValue(' integer(kind=ish)     : ', io_int, 1)
io_int = self%totmem_irg
call Message%WriteValue(' integer(kind=irg)     : ', io_int, 1)
io_int = self%totmem_ill
call Message%WriteValue(' integer(kind=ill)     : ', io_int, 1)
io_int = self%totmem_sgl
call Message%WriteValue(' real(kind=sgl)        : ', io_int, 1)
io_int = self%totmem_dbl
call Message%WriteValue(' real(kind=dbl)        : ', io_int, 1)
io_int = self%totmem_cmplx
call Message%WriteValue(' complex(kind=sgl)     : ', io_int, 1)
io_int = self%totmem_cmplxd
call Message%WriteValue(' complex(kind=dbl)     : ', io_int, 1)

call Message%printMessage(' ')
io_int = self%current_memory_allocated
call Message%WriteValue(' total allocated bytes : ', io_int, 1)

end subroutine show_allocated_memory_use_

!--------------------------------------------------------------------------
subroutine update_total_memory_use_( self, value )
!DEC$ ATTRIBUTES DLLEXPORT :: update_total_memory_use_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! updates the total memory use

IMPLICIT NONE

class(memory_T), INTENT(INOUT)      :: self
integer(kind=irg),INTENT(IN)        :: value 

self%current_memory_allocated = self%current_memory_allocated + value

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
subroutine alloc_ish1_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ish1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type ish1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable       :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
integer(kind=ish), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(1) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ish
    self%totmem_ish = self%totmem_ish - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_ish1_:', &
        ' Unable to allocate integer(kind=ish) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_ish = self%totmem_ish + product(dims)*self%bytes_ish 
call self%update_total_memory_use_(product(dims)*self%bytes_ish)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ish
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_ish
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
subroutine dealloc_ish1_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ish1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type ish1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable    :: ar(:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ish
    self%totmem_ish = self%totmem_ish - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ish1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ish1_

!--------------------------------------------------------------------------
subroutine alloc_ish2_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ish2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type ish2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable       :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
integer(kind=ish), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(2) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ish
    self%totmem_ish = self%totmem_ish - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_ish2_:', &
        ' Unable to allocate integer(kind=ish) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_ish = self%totmem_ish + product(dims)*self%bytes_ish 
call self%update_total_memory_use_(product(dims)*self%bytes_ish)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ish
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_ish
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
subroutine dealloc_ish2_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ish2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type ish2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable    :: ar(:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ish
    self%totmem_ish = self%totmem_ish - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ish2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ish2_

!--------------------------------------------------------------------------
subroutine alloc_ish3_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ish3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type ish3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable       :: ar(:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(3)
character(*),INTENT(IN)                          :: varname
integer(kind=ish), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(3) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ish
    self%totmem_ish = self%totmem_ish - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_ish3_:', &
        ' Unable to allocate integer(kind=ish) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_ish = self%totmem_ish + product(dims)*self%bytes_ish 
call self%update_total_memory_use_(product(dims)*self%bytes_ish)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ish
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_ish
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
subroutine dealloc_ish3_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ish3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type ish3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable    :: ar(:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ish
    self%totmem_ish = self%totmem_ish - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ish3_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ish3_

!--------------------------------------------------------------------------
subroutine alloc_ish4_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ish4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type ish4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable       :: ar(:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(4)
character(*),INTENT(IN)                          :: varname
integer(kind=ish), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(4) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ish
    self%totmem_ish = self%totmem_ish - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_ish4_:', &
        ' Unable to allocate integer(kind=ish) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_ish = self%totmem_ish + product(dims)*self%bytes_ish 
call self%update_total_memory_use_(product(dims)*self%bytes_ish)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ish
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_ish
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
subroutine dealloc_ish4_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ish4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type ish4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable    :: ar(:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ish
    self%totmem_ish = self%totmem_ish - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ish4_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ish4_

!--------------------------------------------------------------------------
subroutine alloc_ish5_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ish5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type ish5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(5)
character(*),INTENT(IN)                          :: varname
integer(kind=ish), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(5) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ish
    self%totmem_ish = self%totmem_ish - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_ish5_:', &
        ' Unable to allocate integer(kind=ish) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_ish = self%totmem_ish + product(dims)*self%bytes_ish 
call self%update_total_memory_use_(product(dims)*self%bytes_ish)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ish
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_ish
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
subroutine dealloc_ish5_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ish5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type ish5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ish
    self%totmem_ish = self%totmem_ish - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ish5_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ish5_

!--------------------------------------------------------------------------
subroutine alloc_ish6_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ish6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type ish6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(6)
character(*),INTENT(IN)                          :: varname
integer(kind=ish), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(6) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ish
    self%totmem_ish = self%totmem_ish - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_ish6_:', &
        ' Unable to allocate integer(kind=ish) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_ish = self%totmem_ish + product(dims)*self%bytes_ish 
call self%update_total_memory_use_(product(dims)*self%bytes_ish)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ish
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_ish
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
subroutine dealloc_ish6_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ish6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type ish6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ish), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ish
    self%totmem_ish = self%totmem_ish - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ish6_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ish6_

!--------------------------------------------------------------------------
subroutine alloc_irg1_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_irg1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type irg1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable       :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(1) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_irg
    self%totmem_irg = self%totmem_irg - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_irg1_:', &
        ' Unable to allocate integer(kind=irg) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_irg = self%totmem_irg + product(dims)*self%bytes_irg 
call self%update_total_memory_use_(product(dims)*self%bytes_irg)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_irg
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_irg
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
subroutine dealloc_irg1_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_irg1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type irg1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable    :: ar(:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_irg
    self%totmem_irg = self%totmem_irg - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_irg1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_irg1_

!--------------------------------------------------------------------------
subroutine alloc_irg2_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_irg2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type irg2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable       :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(2) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_irg
    self%totmem_irg = self%totmem_irg - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_irg2_:', &
        ' Unable to allocate integer(kind=irg) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_irg = self%totmem_irg + product(dims)*self%bytes_irg 
call self%update_total_memory_use_(product(dims)*self%bytes_irg)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_irg
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_irg
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
subroutine dealloc_irg2_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_irg2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type irg2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable    :: ar(:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_irg
    self%totmem_irg = self%totmem_irg - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_irg2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_irg2_

!--------------------------------------------------------------------------
subroutine alloc_irg3_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_irg3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type irg3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable       :: ar(:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(3)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(3) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_irg
    self%totmem_irg = self%totmem_irg - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_irg3_:', &
        ' Unable to allocate integer(kind=irg) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_irg = self%totmem_irg + product(dims)*self%bytes_irg 
call self%update_total_memory_use_(product(dims)*self%bytes_irg)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_irg
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_irg
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
subroutine dealloc_irg3_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_irg3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type irg3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable    :: ar(:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_irg
    self%totmem_irg = self%totmem_irg - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_irg3_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_irg3_

!--------------------------------------------------------------------------
subroutine alloc_irg4_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_irg4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type irg4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable       :: ar(:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(4)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(4) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_irg
    self%totmem_irg = self%totmem_irg - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_irg4_:', &
        ' Unable to allocate integer(kind=irg) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_irg = self%totmem_irg + product(dims)*self%bytes_irg 
call self%update_total_memory_use_(product(dims)*self%bytes_irg)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_irg
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_irg
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
subroutine dealloc_irg4_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_irg4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type irg4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable    :: ar(:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_irg
    self%totmem_irg = self%totmem_irg - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_irg4_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_irg4_

!--------------------------------------------------------------------------
subroutine alloc_irg5_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_irg5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type irg5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(5)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(5) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_irg
    self%totmem_irg = self%totmem_irg - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_irg5_:', &
        ' Unable to allocate integer(kind=irg) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_irg = self%totmem_irg + product(dims)*self%bytes_irg 
call self%update_total_memory_use_(product(dims)*self%bytes_irg)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_irg
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_irg
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
subroutine dealloc_irg5_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_irg5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type irg5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_irg
    self%totmem_irg = self%totmem_irg - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_irg5_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_irg5_

!--------------------------------------------------------------------------
subroutine alloc_irg6_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_irg6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type irg6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(6)
character(*),INTENT(IN)                          :: varname
integer(kind=irg), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(6) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_irg
    self%totmem_irg = self%totmem_irg - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_irg6_:', &
        ' Unable to allocate integer(kind=irg) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_irg = self%totmem_irg + product(dims)*self%bytes_irg 
call self%update_total_memory_use_(product(dims)*self%bytes_irg)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_irg
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_irg
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
subroutine dealloc_irg6_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_irg6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type irg6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=irg), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_irg
    self%totmem_irg = self%totmem_irg - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_irg6_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_irg6_

!--------------------------------------------------------------------------
subroutine alloc_ill1_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ill1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type ill1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable       :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
integer(kind=ill), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(1) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ill
    self%totmem_ill = self%totmem_ill - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_ill1_:', &
        ' Unable to allocate integer(kind=ill) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_ill = self%totmem_ill + product(dims)*self%bytes_ill 
call self%update_total_memory_use_(product(dims)*self%bytes_ill)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ill
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_ill
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
subroutine dealloc_ill1_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ill1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type ill1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable    :: ar(:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ill
    self%totmem_ill = self%totmem_ill - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ill1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ill1_

!--------------------------------------------------------------------------
subroutine alloc_ill2_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ill2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type ill2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable       :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
integer(kind=ill), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(2) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ill
    self%totmem_ill = self%totmem_ill - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_ill2_:', &
        ' Unable to allocate integer(kind=ill) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_ill = self%totmem_ill + product(dims)*self%bytes_ill 
call self%update_total_memory_use_(product(dims)*self%bytes_ill)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ill
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_ill
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
subroutine dealloc_ill2_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ill2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type ill2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable    :: ar(:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ill
    self%totmem_ill = self%totmem_ill - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ill2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ill2_

!--------------------------------------------------------------------------
subroutine alloc_ill3_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ill3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type ill3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable       :: ar(:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(3)
character(*),INTENT(IN)                          :: varname
integer(kind=ill), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(3) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ill
    self%totmem_ill = self%totmem_ill - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_ill3_:', &
        ' Unable to allocate integer(kind=ill) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_ill = self%totmem_ill + product(dims)*self%bytes_ill 
call self%update_total_memory_use_(product(dims)*self%bytes_ill)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ill
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_ill
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
subroutine dealloc_ill3_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ill3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type ill3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable    :: ar(:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ill
    self%totmem_ill = self%totmem_ill - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ill3_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ill3_

!--------------------------------------------------------------------------
subroutine alloc_ill4_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ill4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type ill4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable       :: ar(:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(4)
character(*),INTENT(IN)                          :: varname
integer(kind=ill), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(4) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ill
    self%totmem_ill = self%totmem_ill - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_ill4_:', &
        ' Unable to allocate integer(kind=ill) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_ill = self%totmem_ill + product(dims)*self%bytes_ill 
call self%update_total_memory_use_(product(dims)*self%bytes_ill)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ill
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_ill
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
subroutine dealloc_ill4_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ill4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type ill4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable    :: ar(:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ill
    self%totmem_ill = self%totmem_ill - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ill4_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ill4_

!--------------------------------------------------------------------------
subroutine alloc_ill5_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ill5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type ill5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(5)
character(*),INTENT(IN)                          :: varname
integer(kind=ill), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(5) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ill
    self%totmem_ill = self%totmem_ill - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_ill5_:', &
        ' Unable to allocate integer(kind=ill) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_ill = self%totmem_ill + product(dims)*self%bytes_ill 
call self%update_total_memory_use_(product(dims)*self%bytes_ill)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ill
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_ill
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
subroutine dealloc_ill5_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ill5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type ill5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ill
    self%totmem_ill = self%totmem_ill - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ill5_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ill5_

!--------------------------------------------------------------------------
subroutine alloc_ill6_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_ill6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type ill6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(6)
character(*),INTENT(IN)                          :: varname
integer(kind=ill), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(6) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_ill
    self%totmem_ill = self%totmem_ill - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_ill6_:', &
        ' Unable to allocate integer(kind=ill) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_ill = self%totmem_ill + product(dims)*self%bytes_ill 
call self%update_total_memory_use_(product(dims)*self%bytes_ill)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0_ill
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_ill
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
subroutine dealloc_ill6_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_ill6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type ill6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
integer(kind=ill), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_ill
    self%totmem_ill = self%totmem_ill - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_ill6_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_ill6_

!--------------------------------------------------------------------------
subroutine alloc_sgl1_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_sgl1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type sgl1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable       :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
real(kind=sgl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(1) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_sgl
    self%totmem_sgl = self%totmem_sgl - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_sgl1_:', &
        ' Unable to allocate real(kind=sgl) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_sgl = self%totmem_sgl + product(dims)*self%bytes_sgl 
call self%update_total_memory_use_(product(dims)*self%bytes_sgl)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.0_sgl
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_sgl
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
subroutine dealloc_sgl1_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_sgl1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type sgl1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable    :: ar(:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_sgl
    self%totmem_sgl = self%totmem_sgl - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_sgl1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_sgl1_

!--------------------------------------------------------------------------
subroutine alloc_sgl2_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_sgl2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type sgl2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
real(kind=sgl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(2) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_sgl
    self%totmem_sgl = self%totmem_sgl - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_sgl2_:', &
        ' Unable to allocate real(kind=sgl) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_sgl = self%totmem_sgl + product(dims)*self%bytes_sgl 
call self%update_total_memory_use_(product(dims)*self%bytes_sgl)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.0_sgl
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_sgl
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
subroutine dealloc_sgl2_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_sgl2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type sgl2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_sgl
    self%totmem_sgl = self%totmem_sgl - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_sgl2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_sgl2_

!--------------------------------------------------------------------------
subroutine alloc_sgl3_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_sgl3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type sgl3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(3)
character(*),INTENT(IN)                          :: varname
real(kind=sgl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(3) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_sgl
    self%totmem_sgl = self%totmem_sgl - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_sgl3_:', &
        ' Unable to allocate real(kind=sgl) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_sgl = self%totmem_sgl + product(dims)*self%bytes_sgl 
call self%update_total_memory_use_(product(dims)*self%bytes_sgl)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.0_sgl
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_sgl
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
subroutine dealloc_sgl3_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_sgl3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type sgl3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_sgl
    self%totmem_sgl = self%totmem_sgl - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_sgl3_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_sgl3_

!--------------------------------------------------------------------------
subroutine alloc_sgl4_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_sgl4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type sgl4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(4)
character(*),INTENT(IN)                          :: varname
real(kind=sgl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(4) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_sgl
    self%totmem_sgl = self%totmem_sgl - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_sgl4_:', &
        ' Unable to allocate real(kind=sgl) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_sgl = self%totmem_sgl + product(dims)*self%bytes_sgl 
call self%update_total_memory_use_(product(dims)*self%bytes_sgl)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.0_sgl
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_sgl
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
subroutine dealloc_sgl4_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_sgl4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type sgl4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_sgl
    self%totmem_sgl = self%totmem_sgl - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_sgl4_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_sgl4_

!--------------------------------------------------------------------------
subroutine alloc_sgl5_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_sgl5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type sgl5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(5)
character(*),INTENT(IN)                          :: varname
real(kind=sgl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(5) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_sgl
    self%totmem_sgl = self%totmem_sgl - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_sgl5_:', &
        ' Unable to allocate real(kind=sgl) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_sgl = self%totmem_sgl + product(dims)*self%bytes_sgl 
call self%update_total_memory_use_(product(dims)*self%bytes_sgl)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.0_sgl
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_sgl
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
subroutine dealloc_sgl5_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_sgl5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type sgl5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_sgl
    self%totmem_sgl = self%totmem_sgl - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_sgl5_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_sgl5_

!--------------------------------------------------------------------------
subroutine alloc_sgl6_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_sgl6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type sgl6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(6)
character(*),INTENT(IN)                          :: varname
real(kind=sgl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(6) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_sgl
    self%totmem_sgl = self%totmem_sgl - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_sgl6_:', &
        ' Unable to allocate real(kind=sgl) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_sgl = self%totmem_sgl + product(dims)*self%bytes_sgl 
call self%update_total_memory_use_(product(dims)*self%bytes_sgl)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.0_sgl
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_sgl
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
subroutine dealloc_sgl6_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_sgl6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type sgl6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_sgl
    self%totmem_sgl = self%totmem_sgl - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_sgl6_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_sgl6_

!--------------------------------------------------------------------------
subroutine alloc_dbl1_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_dbl1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type dbl1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable       :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
real(kind=dbl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(1) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_dbl
    self%totmem_dbl = self%totmem_dbl - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_dbl1_:', &
        ' Unable to allocate real(kind=dbl) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_dbl = self%totmem_dbl + product(dims)*self%bytes_dbl 
call self%update_total_memory_use_(product(dims)*self%bytes_dbl)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.D0
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_dbl
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
subroutine dealloc_dbl1_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_dbl1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type dbl1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable    :: ar(:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_dbl
    self%totmem_dbl = self%totmem_dbl - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_dbl1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_dbl1_

!--------------------------------------------------------------------------
subroutine alloc_dbl2_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_dbl2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type dbl2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
real(kind=dbl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(2) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_dbl
    self%totmem_dbl = self%totmem_dbl - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_dbl2_:', &
        ' Unable to allocate real(kind=dbl) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_dbl = self%totmem_dbl + product(dims)*self%bytes_dbl 
call self%update_total_memory_use_(product(dims)*self%bytes_dbl)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.D0
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_dbl
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
subroutine dealloc_dbl2_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_dbl2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type dbl2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_dbl
    self%totmem_dbl = self%totmem_dbl - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_dbl2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_dbl2_

!--------------------------------------------------------------------------
subroutine alloc_dbl3_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_dbl3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type dbl3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(3)
character(*),INTENT(IN)                          :: varname
real(kind=dbl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(3) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_dbl
    self%totmem_dbl = self%totmem_dbl - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_dbl3_:', &
        ' Unable to allocate real(kind=dbl) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_dbl = self%totmem_dbl + product(dims)*self%bytes_dbl 
call self%update_total_memory_use_(product(dims)*self%bytes_dbl)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.D0
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_dbl
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
subroutine dealloc_dbl3_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_dbl3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type dbl3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_dbl
    self%totmem_dbl = self%totmem_dbl - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_dbl3_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_dbl3_

!--------------------------------------------------------------------------
subroutine alloc_dbl4_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_dbl4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type dbl4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(4)
character(*),INTENT(IN)                          :: varname
real(kind=dbl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(4) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_dbl
    self%totmem_dbl = self%totmem_dbl - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_dbl4_:', &
        ' Unable to allocate real(kind=dbl) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_dbl = self%totmem_dbl + product(dims)*self%bytes_dbl 
call self%update_total_memory_use_(product(dims)*self%bytes_dbl)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.D0
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_dbl
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
subroutine dealloc_dbl4_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_dbl4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type dbl4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_dbl
    self%totmem_dbl = self%totmem_dbl - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_dbl4_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_dbl4_

!--------------------------------------------------------------------------
subroutine alloc_dbl5_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_dbl5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type dbl5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(5)
character(*),INTENT(IN)                          :: varname
real(kind=dbl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(5) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_dbl
    self%totmem_dbl = self%totmem_dbl - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_dbl5_:', &
        ' Unable to allocate real(kind=dbl) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_dbl = self%totmem_dbl + product(dims)*self%bytes_dbl 
call self%update_total_memory_use_(product(dims)*self%bytes_dbl)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.D0
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_dbl
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
subroutine dealloc_dbl5_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_dbl5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type dbl5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_dbl
    self%totmem_dbl = self%totmem_dbl - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_dbl5_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_dbl5_

!--------------------------------------------------------------------------
subroutine alloc_dbl6_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_dbl6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type dbl6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(6)
character(*),INTENT(IN)                          :: varname
real(kind=dbl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(6) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_dbl
    self%totmem_dbl = self%totmem_dbl - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_dbl6_:', &
        ' Unable to allocate real(kind=dbl) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_dbl = self%totmem_dbl + product(dims)*self%bytes_dbl 
call self%update_total_memory_use_(product(dims)*self%bytes_dbl)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = 0.D0
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_dbl
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
subroutine dealloc_dbl6_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_dbl6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type dbl6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
real(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_dbl
    self%totmem_dbl = self%totmem_dbl - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_dbl6_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_dbl6_

!--------------------------------------------------------------------------
subroutine alloc_cmplx1_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplx1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type cmplx1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable       :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
complex(kind=sgl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(1) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplx
    self%totmem_cmplx = self%totmem_cmplx - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_cmplx1_:', &
        ' Unable to allocate complex(kind=cmplx) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_cmplx = self%totmem_cmplx + product(dims)*self%bytes_cmplx 
call self%update_total_memory_use_(product(dims)*self%bytes_cmplx)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.0_sgl,0.0_sgl)
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_cmplx
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
subroutine dealloc_cmplx1_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplx1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type cmplx1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable    :: ar(:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplx
    self%totmem_cmplx = self%totmem_cmplx - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplx1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplx1_

!--------------------------------------------------------------------------
subroutine alloc_cmplx2_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplx2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type cmplx2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
complex(kind=sgl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(2) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplx
    self%totmem_cmplx = self%totmem_cmplx - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_cmplx2_:', &
        ' Unable to allocate complex(kind=cmplx) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_cmplx = self%totmem_cmplx + product(dims)*self%bytes_cmplx 
call self%update_total_memory_use_(product(dims)*self%bytes_cmplx)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.0_sgl,0.0_sgl)
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_cmplx
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
subroutine dealloc_cmplx2_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplx2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type cmplx2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplx
    self%totmem_cmplx = self%totmem_cmplx - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplx2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplx2_

!--------------------------------------------------------------------------
subroutine alloc_cmplx3_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplx3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type cmplx3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(3)
character(*),INTENT(IN)                          :: varname
complex(kind=sgl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(3) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplx
    self%totmem_cmplx = self%totmem_cmplx - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_cmplx3_:', &
        ' Unable to allocate complex(kind=cmplx) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_cmplx = self%totmem_cmplx + product(dims)*self%bytes_cmplx 
call self%update_total_memory_use_(product(dims)*self%bytes_cmplx)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.0_sgl,0.0_sgl)
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_cmplx
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
subroutine dealloc_cmplx3_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplx3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type cmplx3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplx
    self%totmem_cmplx = self%totmem_cmplx - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplx3_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplx3_

!--------------------------------------------------------------------------
subroutine alloc_cmplx4_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplx4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type cmplx4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(4)
character(*),INTENT(IN)                          :: varname
complex(kind=sgl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(4) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplx
    self%totmem_cmplx = self%totmem_cmplx - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_cmplx4_:', &
        ' Unable to allocate complex(kind=cmplx) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_cmplx = self%totmem_cmplx + product(dims)*self%bytes_cmplx 
call self%update_total_memory_use_(product(dims)*self%bytes_cmplx)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.0_sgl,0.0_sgl)
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_cmplx
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
subroutine dealloc_cmplx4_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplx4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type cmplx4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplx
    self%totmem_cmplx = self%totmem_cmplx - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplx4_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplx4_

!--------------------------------------------------------------------------
subroutine alloc_cmplx5_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplx5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type cmplx5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(5)
character(*),INTENT(IN)                          :: varname
complex(kind=sgl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(5) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplx
    self%totmem_cmplx = self%totmem_cmplx - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_cmplx5_:', &
        ' Unable to allocate complex(kind=cmplx) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_cmplx = self%totmem_cmplx + product(dims)*self%bytes_cmplx 
call self%update_total_memory_use_(product(dims)*self%bytes_cmplx)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.0_sgl,0.0_sgl)
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_cmplx
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
subroutine dealloc_cmplx5_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplx5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type cmplx5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplx
    self%totmem_cmplx = self%totmem_cmplx - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplx5_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplx5_

!--------------------------------------------------------------------------
subroutine alloc_cmplx6_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplx6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type cmplx6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(6)
character(*),INTENT(IN)                          :: varname
complex(kind=sgl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(6) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplx
    self%totmem_cmplx = self%totmem_cmplx - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_cmplx6_:', &
        ' Unable to allocate complex(kind=cmplx) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_cmplx = self%totmem_cmplx + product(dims)*self%bytes_cmplx 
call self%update_total_memory_use_(product(dims)*self%bytes_cmplx)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.0_sgl,0.0_sgl)
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_cmplx
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
subroutine dealloc_cmplx6_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplx6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type cmplx6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=sgl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplx
    self%totmem_cmplx = self%totmem_cmplx - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplx6_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplx6_

!--------------------------------------------------------------------------
subroutine alloc_cmplxd1_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplxd1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type cmplxd1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable       :: ar(:)
integer(kind=irg), INTENT(IN)                    :: dims(1)
character(*),INTENT(IN)                          :: varname
complex(kind=dbl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(1) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplxd
    self%totmem_cmplxd = self%totmem_cmplxd - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_cmplxd1_:', &
        ' Unable to allocate complex(kind=cmplxd) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_cmplxd = self%totmem_cmplxd + product(dims)*self%bytes_cmplxd 
call self%update_total_memory_use_(product(dims)*self%bytes_cmplxd)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.D0,0.D0)
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_cmplxd
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
subroutine dealloc_cmplxd1_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplxd1_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type cmplxd1

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable    :: ar(:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplxd
    self%totmem_cmplxd = self%totmem_cmplxd - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplxd1_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplxd1_

!--------------------------------------------------------------------------
subroutine alloc_cmplxd2_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplxd2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type cmplxd2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:)
integer(kind=irg), INTENT(IN)                    :: dims(2)
character(*),INTENT(IN)                          :: varname
complex(kind=dbl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(2) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplxd
    self%totmem_cmplxd = self%totmem_cmplxd - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_cmplxd2_:', &
        ' Unable to allocate complex(kind=cmplxd) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_cmplxd = self%totmem_cmplxd + product(dims)*self%bytes_cmplxd 
call self%update_total_memory_use_(product(dims)*self%bytes_cmplxd)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.D0,0.D0)
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_cmplxd
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
subroutine dealloc_cmplxd2_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplxd2_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type cmplxd2

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplxd
    self%totmem_cmplxd = self%totmem_cmplxd - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplxd2_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplxd2_

!--------------------------------------------------------------------------
subroutine alloc_cmplxd3_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplxd3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type cmplxd3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(3)
character(*),INTENT(IN)                          :: varname
complex(kind=dbl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(3) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplxd
    self%totmem_cmplxd = self%totmem_cmplxd - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_cmplxd3_:', &
        ' Unable to allocate complex(kind=cmplxd) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_cmplxd = self%totmem_cmplxd + product(dims)*self%bytes_cmplxd 
call self%update_total_memory_use_(product(dims)*self%bytes_cmplxd)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.D0,0.D0)
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_cmplxd
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
subroutine dealloc_cmplxd3_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplxd3_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type cmplxd3

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplxd
    self%totmem_cmplxd = self%totmem_cmplxd - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplxd3_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplxd3_

!--------------------------------------------------------------------------
subroutine alloc_cmplxd4_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplxd4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type cmplxd4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(4)
character(*),INTENT(IN)                          :: varname
complex(kind=dbl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(4) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplxd
    self%totmem_cmplxd = self%totmem_cmplxd - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_cmplxd4_:', &
        ' Unable to allocate complex(kind=cmplxd) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_cmplxd = self%totmem_cmplxd + product(dims)*self%bytes_cmplxd 
call self%update_total_memory_use_(product(dims)*self%bytes_cmplxd)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.D0,0.D0)
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_cmplxd
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
subroutine dealloc_cmplxd4_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplxd4_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type cmplxd4

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplxd
    self%totmem_cmplxd = self%totmem_cmplxd - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplxd4_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplxd4_

!--------------------------------------------------------------------------
subroutine alloc_cmplxd5_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplxd5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type cmplxd5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(5)
character(*),INTENT(IN)                          :: varname
complex(kind=dbl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(5) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplxd
    self%totmem_cmplxd = self%totmem_cmplxd - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_cmplxd5_:', &
        ' Unable to allocate complex(kind=cmplxd) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_cmplxd = self%totmem_cmplxd + product(dims)*self%bytes_cmplxd 
call self%update_total_memory_use_(product(dims)*self%bytes_cmplxd)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.D0,0.D0)
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_cmplxd
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
subroutine dealloc_cmplxd5_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplxd5_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type cmplxd5

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplxd
    self%totmem_cmplxd = self%totmem_cmplxd - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplxd5_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplxd5_

!--------------------------------------------------------------------------
subroutine alloc_cmplxd6_(self, ar, dims, varname, initval)
!DEC$ ATTRIBUTES DLLEXPORT :: alloc_cmplxd6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! allocate an array of type cmplxd6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable       :: ar(:,:,:,:,:,:)
integer(kind=irg), INTENT(IN)                    :: dims(6)
character(*),INTENT(IN)                          :: varname
complex(kind=dbl), INTENT(IN), OPTIONAL             :: initval

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr, outstr, szstr, initstr
integer(kind=irg)                                :: sz, err, szar(6) 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    szar = size(ar)
    sz = product(szar) * self%bytes_cmplxd
    self%totmem_cmplxd = self%totmem_cmplxd - sz
    deallocate(ar)
endif

! allocate the array 
allocate(ar(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)), stat = err)
if (err.ne.0) then 
    errorstr = ' '
    write (errorstr,*) dims
    call Message%printError('mod_memory:alloc_cmplxd6_:', &
        ' Unable to allocate complex(kind=cmplxd) array'//trim(varname)//' of dimension '//errorstr)
end if
self%totmem_cmplxd = self%totmem_cmplxd + product(dims)*self%bytes_cmplxd 
call self%update_total_memory_use_(product(dims)*self%bytes_cmplxd)

! initalize the array to zeros or initval if present  
if (present(initval)) then 
    ar = initval
else 
    ar = (0.D0,0.D0)
end if

if (self%verbose) then 
  write (szstr,*) product(dims)*self%bytes_cmplxd
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
subroutine dealloc_cmplxd6_(self, ar, varname)
!DEC$ ATTRIBUTES DLLEXPORT :: dealloc_cmplxd6_
!! author: MDG
!! version: 1.0
!! date: 12/12/20
!!
!! deallocate an array of type cmplxd6

IMPLICIT NONE

class(memory_T), INTENT(INOUT)                   :: self
complex(kind=dbl), INTENT(INOUT), allocatable    :: ar(:,:,:,:,:,:)
character(*),INTENT(IN)                          :: varname

type(IO_T)                                       :: Message
character(fnlen)                                 :: errorstr
integer(kind=irg)                                :: err 
integer(kind=irg)                                :: sz 

! if already allocated, then deallocate 
if (allocated(ar)) then 
    ! get the size of the allocated array and decrement the correct counter 
    sz = size(ar) * self%bytes_cmplxd
    self%totmem_cmplxd = self%totmem_cmplxd - sz
    call self%update_total_memory_use_( -sz )
    if (self%verbose) write (*,*) '   -> deallocated array '//trim(varname)
    deallocate(ar)
else 
    call Message%printMessage(' mod_memory:dealloc_cmplxd6_:Warning: attempting to deallocate array that is not allocated. ')
endif

end subroutine dealloc_cmplxd6_




end module mod_memory