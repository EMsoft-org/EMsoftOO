! ###################################################################
! Copyright (c) 2013-2026, Marc De Graef Research Group/Carnegie Mellon University
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

program EMGPUinfo
  !! author: MDG
  !! version: 1.0
  !! date: 01/12/20
  !!
  !! show information for the active GPU backend (OpenCL platforms/devices, or
  !! Apple Metal devices when the Metal backend is enabled).  Both backends
  !! expose the same GPU_T interface, so this single program reports either one.
  !! (formerly EMOpenCLinfo, which is still built as an alias.)

use mod_EMsoft
use mod_global
use mod_GPUsupport

IMPLICIT NONE

character(fnlen)  :: progname, progdesc

type(GPU_T)    :: CL
type(EMsoft_T)    :: EMsoft

progname = 'EMGPUinfo.f90'
progdesc = 'List GPU platform and device information (OpenCL or Apple Metal)'
EMsoft = EMsoft_T(progname, progdesc, tpl = (/ 904, 930 /) )

! initialize the GPU class and print all platform/device information
CL = GPU_T()
call CL%print_platform_info()

end program EMGPUinfo
