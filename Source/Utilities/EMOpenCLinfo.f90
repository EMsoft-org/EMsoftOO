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
! original Copyright information (clfortran's query_platforms_devices.f90)
! -----------------------------------------------------------------------------
!
! Copyright (C) 2013-2014 Company for Advanced Supercomputing Solutions LTD
! Bosmat 2a St.
! Shoham
! Israel 60850
! http://www.cass-hpc.com
!
! Author: Mordechai Butrashvily <support@cass-hpc.com>
!
! -----------------------------------------------------------------------------
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! -----------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! EMsoft:EMOpenCLinfo.f90
!--------------------------------------------------------------------------
!
! PROGRAM: EMOpenCLinfo 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief show information for OpenCL devices (based on clfortran's query_platforms_devices.f90)
!
!> @date   02/18/16 MDG 1.0 original
!--------------------------------------------------------------------------
program EMOpenCLinfo

use mod_EMsoft
use mod_global
use clfortran
use mod_CLsupport
use ISO_C_BINDING

IMPLICIT NONE

character(fnlen)                         :: progname, progdesc

type(OpenCL_T)                           :: CL
type(EMsoft_T)                           :: EMsoft 

progname = 'EMOpenCLinfo.f90'
progdesc = 'List OpenCL platform and device information'
EMsoft = EMsoft_T(progname,progdesc)

! deal with the command line arguments, if any (mostly to generate pdf version of wiki pages)
!call Interpret_Program_Arguments(2,(/ 904, 930 /), progname)

! initialize the OpenCL class and print all OpenCL information
CL = OpenCL_T()
call CL%print_platform_info()

end program EMOpenCLinfo

