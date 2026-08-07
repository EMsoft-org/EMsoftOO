! ###################################################################
! Copyright (c) 2026-2026, Marc De Graef Research Group/Carnegie Mellon University
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

! ###################################################################
! original license text <https://github.com/MRedies/NPY-for-Fortran/tree/master>
!
! MIT License
!
! Copyright (c) 2017 Matthias Redies
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
! ###################################################################

module mod_npz
  !! author: Matthias Redies / adapted for EMsoftOO coding style by MDG 
  !! version: 1.0 
  !! date: 08/07/26
  !!
  !! class definition for writing python .npy and .npz files

use mod_kinds
use mod_global

IMPLICIT NONE 

private 

  integer(kind=irg), parameter        :: p_un = 23
  character, parameter                :: magic_num = achar(147) ! x93
  character, parameter                :: major = achar(2)   !major *.npy version
  character, parameter                :: minor = achar(0)   !minor *.npy version
  character(len=*), parameter         :: zip_flag = "-q0"
  character(len=*), parameter         :: magic_str = "NUMPY"

  public :: save_npy
  interface save_npy
    module procedure write_int64_vec, write_int64_mtx, &
         write_int32_vec, write_int32_mtx, write_int32_3d, &
         write_int16_vec, write_int16_mtx, &
         write_int8_vec, write_int8_mtx, write_int8_3d, &
         write_dbl_vec, write_dbl_mtx, &
         write_sng_vec, write_sng_mtx, &
         write_cmplx_sgn_vec, write_cmplx_sgn_mtx, &
         write_cmplx_dbl_vec, write_cmplx_dbl_mtx, &
         write_sng_3dT, write_dbl_3dT, &
         write_sng_4dT, write_dbl_4dT, &
         write_dbl_5dT, &
         write_cmplx_dbl_3dT, &
         write_cmplx_dbl_4dT, &
         write_cmplx_dbl_5dT, &
         write_cmplx_dbl_6dT, &
         write_int64_mtx_3d, &
         write_int16_mtx_3d, &
         write_cmplx_sgn_mtx_3d, &
         write_int32_scalar, write_str
  end interface save_npy

  public :: add_npz
  interface add_npz
    module procedure addrpl_int8_vec, addrpl_int8_mtx, &
         addrpl_int16_vec, addrpl_int16_mtx, &
         addrpl_int32_vec, addrpl_int32_mtx, &
         addrpl_int64_vec, addrpl_int64_mtx, &
         addrpl_sng_vec, addrpl_sng_mtx, &
         addrpl_dbl_vec, addrpl_dbl_mtx, &
         addrpl_cmplx_dbl_vec, addrpl_cmplx_dbl_mtx, &
         addrpl_cmplx_sng_vec, addrpl_cmplx_sng_mtx, &
         addrpl_int8_mtx_3d, &
         addrpl_int16_mtx_3d, &
         addrpl_int32_mtx_3d, &
         addrpl_int64_mtx_3d, &
         addrpl_sng_mtx_3d, &
         addrpl_dbl_mtx_3d, &
         addrpl_cmplx_dbl_mtx_3d, &
         addrpl_cmplx_sng_mtx_3d, &
         addrpl_sng_4d, addrpl_dbl_4d, &
         addrpl_int32_scalar, addrpl_str
  end interface add_npz

contains
  subroutine run_sys(cmd, stat)
      implicit none
      character(len=*), intent(in)     :: cmd
      integer(4), intent(out)          :: stat

      call execute_command_line(cmd, wait=.True., exitstat=stat)
   end subroutine run_sys

   subroutine addrpl_cmplx_sng_vec(zipfile, var_name, vec)
      implicit none
      complex(4), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_cmplx_sng_vec

   subroutine addrpl_cmplx_sng_mtx(zipfile, var_name, mtx)
      implicit none
      complex(4), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_cmplx_sng_mtx

   subroutine addrpl_cmplx_sng_mtx_3d(zipfile, var_name, mtx)
      implicit none
      complex(4), intent(in)           :: mtx(:, :, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_cmplx_sng_mtx_3d

   subroutine addrpl_cmplx_dbl_vec(zipfile, var_name, vec)
      implicit none
      complex(8), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_cmplx_dbl_vec

   subroutine addrpl_cmplx_dbl_mtx(zipfile, var_name, mtx)
      implicit none
      complex(8), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_cmplx_dbl_mtx

   subroutine addrpl_cmplx_dbl_mtx_3d(zipfile, var_name, mtx)
      implicit none
      complex(8), intent(in)           :: mtx(:, :, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_cmplx_dbl_mtx_3d

   subroutine addrpl_dbl_vec(zipfile, var_name, vec)
      implicit none
      real(8), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_dbl_vec

   subroutine addrpl_dbl_mtx(zipfile, var_name, mtx)
      implicit none
      real(8), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_dbl_mtx

   subroutine addrpl_dbl_mtx_3d(zipfile, var_name, mtx)
      implicit none
      real(8), intent(in)           :: mtx(:, :, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_dbl_mtx_3d

   subroutine addrpl_sng_vec(zipfile, var_name, vec)
      implicit none
      real(4), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_sng_vec

   subroutine addrpl_sng_mtx(zipfile, var_name, mtx)
      implicit none
      real(4), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_sng_mtx

   subroutine addrpl_sng_mtx_3d(zipfile, var_name, mtx)
      implicit none
      real(4), intent(in)           :: mtx(:, :, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_sng_mtx_3d

   subroutine addrpl_int8_vec(zipfile, var_name, vec)
      implicit none
      integer(1), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int8_vec

   subroutine addrpl_int8_mtx(zipfile, var_name, mtx)
      implicit none
      integer(1), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int8_mtx

   subroutine addrpl_int8_mtx_3d(zipfile, var_name, mtx)
      implicit none
      integer(1), intent(in)           :: mtx(:, :, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int8_mtx_3d

   subroutine addrpl_int16_vec(zipfile, var_name, vec)
      implicit none
      integer(2), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int16_vec

   subroutine addrpl_int16_mtx(zipfile, var_name, mtx)
      implicit none
      integer(2), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int16_mtx

   subroutine addrpl_int16_mtx_3d(zipfile, var_name, mtx)
      implicit none
      integer(2), intent(in)           :: mtx(:, :, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int16_mtx_3d

   subroutine addrpl_int32_vec(zipfile, var_name, vec)
      implicit none
      integer(4), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int32_vec

   subroutine addrpl_int32_mtx(zipfile, var_name, mtx)
      implicit none
      integer(4), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int32_mtx

   subroutine addrpl_int32_mtx_3d(zipfile, var_name, mtx)
      implicit none
      integer(4), intent(in)           :: mtx(:, :, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int32_mtx_3d

   subroutine addrpl_int64_vec(zipfile, var_name, vec)
      implicit none
      integer(8), intent(in)           :: vec(:)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, vec)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int64_vec

   subroutine addrpl_int64_mtx(zipfile, var_name, mtx)
      implicit none
      integer(8), intent(in)           :: mtx(:, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int64_mtx

   subroutine addrpl_int64_mtx_3d(zipfile, var_name, mtx)
      implicit none
      integer(8), intent(in)           :: mtx(:, :, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, mtx)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int64_mtx_3d

   subroutine addrpl_sng_4d(zipfile, var_name, tensor)
      implicit none
      real(4), intent(in)              :: tensor(:, :, :, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, tensor)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_sng_4d

   subroutine addrpl_dbl_4d(zipfile, var_name, tensor)
      implicit none
      real(8), intent(in)              :: tensor(:, :, :, :)
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, tensor)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_dbl_4d

   subroutine addrpl_int32_scalar(zipfile, var_name, val)
      implicit none
      integer(4), intent(in)           :: val
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, val)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_int32_scalar

   subroutine addrpl_str(zipfile, var_name, str)
      implicit none
      character(len=*), intent(in)     :: str
      character(len=*), intent(in)     :: zipfile, var_name
      character(len=:), allocatable    :: npy_name
      integer(4)                       :: succ

      npy_name = var_name//".npy"

      call save_npy(npy_name, str)
      ! just store and be quite while zipping
      call run_sys("zip "//zip_flag//" "//zipfile &
                   //" "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute zip command"
      endif

      call run_sys("rm "//npy_name, succ)
      if (succ /= 0) then
         write (*, *) "Can't execute rm command"
      endif
   end subroutine addrpl_str

   Subroutine write_cmplx_sgn_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(4), intent(in)           :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<c8"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor
      write (p_un) header_len
      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_cmplx_sgn_mtx

   Subroutine write_cmplx_sgn_mtx_3d(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(4), intent(in)           :: mtx(:, :, :)
      character(len=*), parameter      :: var_type = "<c8"
      integer(4)                       :: header_len, s_mtx(3), i, j, k

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor
      write (p_un) header_len
      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_cmplx_sgn_mtx_3d

   Subroutine write_cmplx_sgn_vec(filename, vec)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(4), intent(in)           :: vec(:)
      character(len=*), parameter      :: var_type = "<c8"
      integer(4)                       :: header_len, s_vec(1), i

      s_vec = shape(vec)
      header_len = len(dict_str(var_type, s_vec))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor
      write (p_un) header_len

      write (p_un) dict_str(var_type, s_vec)

      write (p_un) vec

      close (unit=p_un)
   End Subroutine write_cmplx_sgn_vec

   Subroutine write_cmplx_dbl_6dT(filename, tensor)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(8), intent(in)           :: tensor(:, :, :, :, :, :)
      character(len=*), parameter      :: var_type = "<c16"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor))
      write (p_un) tensor
      close (unit=p_un)
   End Subroutine write_cmplx_dbl_6dT

   Subroutine write_cmplx_dbl_5dT(filename, tensor)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(8), intent(in)           :: tensor(:, :, :, :, :)
      character(len=*), parameter      :: var_type = "<c16"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor))
      write (p_un) tensor
      close (unit=p_un)
   End Subroutine write_cmplx_dbl_5dT

   Subroutine write_cmplx_dbl_4dT(filename, tensor)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(8), intent(in)           :: tensor(:, :, :, :)
      character(len=*), parameter      :: var_type = "<c16"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor))
      write (p_un) tensor
      close (unit=p_un)
   End Subroutine write_cmplx_dbl_4dT

   Subroutine write_cmplx_dbl_3dT(filename, tensor)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(8), intent(in)           :: tensor(:, :, :)
      character(len=*), parameter      :: var_type = "<c16"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor))
      write (p_un) tensor
      close (unit=p_un)
   End Subroutine write_cmplx_dbl_3dT

   Subroutine write_cmplx_dbl_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(8), intent(in)           :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<c16"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_cmplx_dbl_mtx

   Subroutine write_cmplx_dbl_vec(filename, vec)
      Implicit None
      character(len=*), intent(in)     :: filename
      complex(8), intent(in)           :: vec(:)
      character(len=*), parameter      :: var_type = "<c16"
      integer(4)                       :: header_len, s_vec(1), i

      s_vec = shape(vec)
      header_len = len(dict_str(var_type, s_vec))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_vec)

      write (p_un) vec

      close (unit=p_un)
   End Subroutine write_cmplx_dbl_vec

   Subroutine write_sng_3dT(filename, tensor)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(4), intent(in)              :: tensor(:, :, :)
      character(len=*), parameter      :: var_type = "<f4"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor))
      write (p_un) tensor
      close (unit=p_un)
   End Subroutine write_sng_3dT

   Subroutine write_sng_4dT(filename, tensor)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(4), intent(in)              :: tensor(:, :, :, :)
      character(len=*), parameter      :: var_type = "<f4"
      integer(4)                       :: header_len

      header_len = len(dict_str(var_type, shape(tensor)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor))
      write (p_un) tensor
      close (unit=p_un)
   End Subroutine write_sng_4dT

   Subroutine write_sng_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(4), intent(in)              :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<f4"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_sng_mtx

   Subroutine write_sng_vec(filename, vec)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(4), intent(in)              :: vec(:)
      character(len=*), parameter      :: var_type = "<f4"
      integer(4)                       :: header_len, s_vec(1), i

      s_vec = shape(vec)
      header_len = len(dict_str(var_type, s_vec))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_vec)

      write (p_un) vec

      close (unit=p_un)
   End Subroutine write_sng_vec

   Subroutine write_dbl_3dT(filename, tensor)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(8), intent(in)              :: tensor(:, :, :)
      character(len=*), parameter      :: var_type = "<f8"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor))
      write (p_un) tensor
      close (unit=p_un)
   End Subroutine write_dbl_3dT

   Subroutine write_dbl_4dT(filename, tensor4)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(8), intent(in)              :: tensor4(:, :, :, :)
      character(len=*), parameter      :: var_type = "<f8"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor4)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor4))
      write (p_un) tensor4
      close (unit=p_un)
   End Subroutine write_dbl_4dT

   Subroutine write_dbl_5dT(filename, tensor5)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(8), intent(in)              :: tensor5(:, :, :, :, :)
      character(len=*), parameter      :: var_type = "<f8"
      integer(4)                       :: header_len, i, j, k

      header_len = len(dict_str(var_type, shape(tensor5)))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, shape(tensor5))
      write (p_un) tensor5
      close (unit=p_un)
   End Subroutine write_dbl_5dT

   Subroutine write_dbl_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(8), intent(in)              :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<f8"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_dbl_mtx

   Subroutine write_dbl_vec(filename, vec)
      Implicit None
      character(len=*), intent(in)     :: filename
      real(8), intent(in)              :: vec(:)
      character(len=*), parameter      :: var_type = "<f8"
      integer(4)                       :: header_len, s_vec(1), i

      s_vec = shape(vec)
      header_len = len(dict_str(var_type, s_vec))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_vec)

      write (p_un) vec

      close (unit=p_un)
   End Subroutine write_dbl_vec

   Subroutine write_int64_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(8), intent(in)           :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<i8"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_int64_mtx

   Subroutine write_int64_mtx_3d(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(8), intent(in)           :: mtx(:, :, :)
      character(len=*), parameter      :: var_type = "<i8"
      integer(4)                       :: header_len, s_mtx(3), i, j, k

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_int64_mtx_3d

   Subroutine write_int64_vec(filename, vec)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(8), intent(in)           :: vec(:)
      character(len=*), parameter      :: var_type = "<i8"
      integer(4)                       :: header_len, s_vec(1), i

      s_vec = shape(vec)
      header_len = len(dict_str(var_type, s_vec))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_vec)

      write (p_un) vec

      close (unit=p_un)
   End Subroutine write_int64_vec

   Subroutine write_int32_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(4), intent(in)           :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<i4"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_int32_mtx

   Subroutine write_int32_3d(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(4), intent(in)           :: mtx(:,:,:)
      character(len=*), parameter      :: var_type = "<i4"
      integer(4)                       :: header_len, s_mtx(3), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_int32_3d

   Subroutine write_int32_vec(filename, vec)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(4), intent(in)           :: vec(:)
      character(len=*), parameter      :: var_type = "<i4"
      integer(4)                       :: header_len, s_vec(1), i

      s_vec = shape(vec)
      header_len = len(dict_str(var_type, s_vec))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_vec)

      write (p_un) vec

      close (unit=p_un)
   End Subroutine write_int32_vec

   Subroutine write_int16_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(2), intent(in)           :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<i2"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_int16_mtx

   Subroutine write_int16_mtx_3d(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(2), intent(in)           :: mtx(:, :, :)
      character(len=*), parameter      :: var_type = "<i2"
      integer(4)                       :: header_len, s_mtx(3), i, j, k

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_int16_mtx_3d

   Subroutine write_int16_vec(filename, vec)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(2), intent(in)           :: vec(:)
      character(len=*), parameter      :: var_type = "<i2"
      integer(4)                       :: header_len, s_vec(1), i

      s_vec = shape(vec)
      header_len = len(dict_str(var_type, s_vec))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_vec)

      write (p_un) vec

      close (unit=p_un)
   End Subroutine write_int16_vec

   Subroutine write_int8_mtx(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(1), intent(in)           :: mtx(:, :)
      character(len=*), parameter      :: var_type = "<i1"
      integer(4)                       :: header_len, s_mtx(2), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_int8_mtx

   Subroutine write_int8_3d(filename, mtx)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(1), intent(in)           :: mtx(:,:,:)
      character(len=*), parameter      :: var_type = "<i1"
      integer(4)                       :: header_len, s_mtx(3), i, j

      s_mtx = shape(mtx)
      header_len = len(dict_str(var_type, s_mtx))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_mtx)

      write (p_un) mtx

      close (unit=p_un)
   End Subroutine write_int8_3d

   Subroutine write_int8_vec(filename, vec)
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(1), intent(in)           :: vec(:)
      character(len=*), parameter      :: var_type = "<i1"
      integer(4)                       :: header_len, s_vec(1), i

      s_vec = shape(vec)
      header_len = len(dict_str(var_type, s_vec))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor

      write (p_un) header_len

      write (p_un) dict_str(var_type, s_vec)

      write (p_un) vec

      close (unit=p_un)
   End Subroutine write_int8_vec

   Subroutine write_int32_scalar(filename, val)
      ! writes a rank-0 (0-d) array, which is what numpy stores for a python scalar
      Implicit None
      character(len=*), intent(in)     :: filename
      integer(4), intent(in)           :: val
      character(len=*), parameter      :: var_type = "<i4"
      integer(4)                       :: header_len

      header_len = len(scalar_dict_str(var_type))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor
      write (p_un) header_len
      write (p_un) scalar_dict_str(var_type)

      write (p_un) val

      close (unit=p_un)
   End Subroutine write_int32_scalar

   Subroutine write_str(filename, str)
      ! writes a rank-0 numpy unicode array ('<U<n>'), i.e. what numpy stores for
      ! a python string; the payload is n UCS-4 code points.  The string is used
      ! as given, so the caller should trim() it first if that is what is wanted.
      Implicit None
      character(len=*), intent(in)     :: filename
      character(len=*), intent(in)     :: str
      character(len=:), allocatable    :: var_type, cnt_str
      integer(4)                       :: header_len, i, n
      integer(4), allocatable          :: ucs4(:)

      n = len(str)
      allocate (character(24) :: cnt_str)
      write (cnt_str, "(I0)") n
      var_type = "<U"//trim(cnt_str)

      allocate (ucs4(n))
      do i = 1, n
         ucs4(i) = iachar(str(i:i))
      enddo

      header_len = len(scalar_dict_str(var_type))

      open (unit=p_un, file=filename, form="unformatted", &
            access="stream")
      write (p_un) magic_num, magic_str, major, minor
      write (p_un) header_len
      write (p_un) scalar_dict_str(var_type)

      write (p_un) ucs4

      close (unit=p_un)

      deallocate (ucs4)
   End Subroutine write_str

   function scalar_dict_str(var_type) result(str)
      ! header dictionary for a rank-0 array; 'shape': () has no entries at all,
      ! so this cannot go through shape_str
      implicit none
      character(len=*), intent(in)   :: var_type
      character(len=:), allocatable  :: str

      str = "{'descr': '"//var_type// &
            "', 'fortran_order': True, 'shape': (), }"

      do while (mod(len(str) + 11, 16) /= 0)
         str = str//" "
      enddo

      str = str//achar(10)

   end function scalar_dict_str

   function dict_str(var_type, var_shape) result(str)
      implicit none
      character(len=*), intent(in)   :: var_type
      integer(4), intent(in)         :: var_shape(:)
      character(len=:), allocatable  :: str
      integer(4)                     :: cnt

      cnt = len("{'descr': '")
      cnt = cnt + len(var_type)
      cnt = cnt + len("', 'fortran_order': True, 'shape': (")
      cnt = cnt + len(shape_str(var_shape))
      cnt = cnt + len(",), }")
      do while (mod(cnt + 10, 16) /= 0)
         cnt = cnt + 1
      enddo

      allocate (character(cnt) :: str)

      str = "{'descr': '"//var_type// &
            "', 'fortran_order': True, 'shape': ("// &
            shape_str(var_shape)//"), }"

      do while (mod(len(str) + 11, 16) /= 0)
         str = str//" "
      enddo

      str = str//achar(10)

   end function dict_str

   function shape_str(var_shape) result(fin_str)
      implicit none
      integer(4), intent(in)        :: var_shape(:)
      character(len=:), allocatable :: str, small_str, fin_str
      integer(4)                    :: i, length, start, halt

      length = 14*size(var_shape)
      allocate (character(length) :: str)
      allocate (character(14)     :: small_str)
      str = " "

      do i = 1, size(var_shape)
         start = (i - 1)*length + 1
         halt = i*length + 1
         write (small_str, "(I13,A)") var_shape(i), ","
         str = trim(str)//adjustl(small_str)
      enddo

      fin_str = trim(str)
   end function shape_str

end module mod_npz