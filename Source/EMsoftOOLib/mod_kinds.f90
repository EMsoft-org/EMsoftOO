! ###################################################################
! Copyright (c) 2014-2020, Marc De Graef/Carnegie Mellon University
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
! EMsoft:mod_kinds.f90
!--------------------------------------------------------------------------
!
! MODULE: mod_kinds
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief definitions of single and double precision
!
!> @date 12/30/19 MDG 1.0 original
!--------------------------------------------------------------------------
module mod_kinds
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! definitions of single and double precision reals as well as integer types

  use iso_fortran_env, only: int16, int32, int64, real32, real64

  IMPLICIT NONE 

  private

  public :: sgl, dbl, ish, irg, ill 

!> @note This module must be "use"d by every program, subroutine, and function.
!> These are the kind variables used by the EMsoft package.

! Define the "kind" parameters for single and double precision reals,
  integer,parameter                     :: sgl = real32    !  SELECTED_REAL_KIND(p=6,r=37)
   !! single precision real kind parameter
  integer,parameter                     :: dbl = real64    !  SELECTED_REAL_KIND(p=13,r=200)
   !! double precision real kind parameter
!DEC$ ATTRIBUTES DLLEXPORT :: sgl
!DEC$ ATTRIBUTES DLLEXPORT :: dbl

! Define the "kind" parameters for short and regular integers,
  integer,parameter                     :: ish = int16     ! SELECTED_INT_KIND(3)
   !! short integer kind parameter
  integer,parameter                     :: irg = int32     !  SELECTED_INT_KIND(9)
   !! long integer kind parameter
  integer,parameter                     :: ill = int64     !  SELECTED_INT_KIND(12)
   !! long long kind parameter
!DEC$ ATTRIBUTES DLLEXPORT :: ish
!DEC$ ATTRIBUTES DLLEXPORT :: irg
!DEC$ ATTRIBUTES DLLEXPORT :: ill

end module mod_kinds
