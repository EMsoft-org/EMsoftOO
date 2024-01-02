! ###################################################################
! Copyright (c) 2013-2024, Marc De Graef Research Group/Carnegie Mellon University
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
! EMsoft:JSONsupport.f90
!--------------------------------------------------------------------------
!
! MODULE: JSONsupport
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief routines for conversion between json and nml files and reading of json files
!
!> @date  08/11/15 MDG 1.0 original
!> @date  08/12/15 MDG 1.1 added all routines currently also in NameListHDFwriters.f90
!> @date  08/12/15 MDG 1.2 replaced all the json_failed stuff by short routine JSON_failtest
!> @date  11/20/15 MDG 1.3 started defect file format
!--------------------------------------------------------------------------
module mod_JSONsupport

use, intrinsic :: iso_fortran_env , only: error_unit, wp => real64
use json_module
use mod_kinds
use mod_global

IMPLICIT NONE

contains

!--------------------------------------------------------------------------
!
! SUBROUTINE:JSON_failtest
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief executes the json_fail routine; mostly to shorten the remaining code a little
!
!> @param error_cnt error counter
!
!> @date 08/12/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
recursive subroutine JSON_failtest(error_cnt)
!DEC$ ATTRIBUTES DLLEXPORT :: JSON_failtest

IMPLICIT NONE

integer(kind=irg),INTENT(INOUT)         :: error_cnt
!f2py intent(in,out) ::  error_cnt

if (json_failed().eqv..TRUE.) then
  call json_print_error_message(error_unit)
  error_cnt = error_cnt + 1
end if

end subroutine JSON_failtest

!--------------------------------------------------------------------------
!
! FUNCTION:JSONgetDouble
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief get vector from a json_value
!
!> @param child json_value structure
!> @param str text with variable name
!> @param v verbose if 1
!
!> @date 11/21/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
recursive function JSONgetDouble(child,str,v) result(oval)
!DEC$ ATTRIBUTES DLLEXPORT :: JSONgetDouble

use mod_io
use, intrinsic :: iso_fortran_env, only: wp => real64

IMPLICIT NONE

type(json_value), pointer,INTENT(IN)            :: child
type(IO_T)                                      :: Message
character(fnlen)                                :: str
integer(kind=irg),INTENT(IN)                    :: v
real(kind=dbl)                                  :: oval

real(kind=wp)                                   :: val
real(kind=sgl)                                  :: io_real(1)

call json_get(child, val)
if (v.eq.1) then
  io_real(1) = val
  call Message%WriteValue(str,io_real,1)
end if
oval = val

end function JSONgetDouble

!--------------------------------------------------------------------------
!
! FUNCTION:JSONgetDoubleVector
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief get vector from a json_value
!
!> @param child json_value structure
!> @param nc number of items to read
!> @param str text with variable name
!> @param v verbose if 1
!
!> @date 11/21/15  MDG 1.0 new routine
!--------------------------------------------------------------------------
recursive function JSONgetDoubleVector(child,nc,str,v) result(ovec)
!DEC$ ATTRIBUTES DLLEXPORT :: JSONgetDoubleVector

use mod_io
use, intrinsic :: iso_fortran_env, only: wp => real64

IMPLICIT NONE

type(json_value), pointer,INTENT(IN)            :: child
type(IO_T)                                      :: Message
integer(kind=irg),INTENT(IN)                    :: nc
character(fnlen)                                :: str
integer(kind=irg),INTENT(IN)                    :: v
real(kind=dbl)                                  :: ovec(nc)

real(kind=wp),dimension(:),allocatable          :: vec
real(kind=sgl)                                  :: io_real(nc)

allocate(vec(nc))
call json_get(child, vec)
if (v.eq.1) then
  io_real(1:nc) = vec(1:nc)
  call Message%WriteValue(str,io_real,nc)
end if
ovec(1:nc) = vec(1:nc)
deallocate(vec)

end function JSONgetDoubleVector

end module mod_JSONsupport
