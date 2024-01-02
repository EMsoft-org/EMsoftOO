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

module mod_platformsupport
  !! author: MDG 
  !! version: 1.0 
  !! date: 08/29/23
  !!
  !! class definition for the platformsupport module
  !!
  !! This file is automatically edited by cmake at build time to select either the 
  !! routines from the IFPORT module (if using the ifort compiler on Windows) or 
  !! the gfortran extensions (on all other platforms).  CMAKE will replace the 
  !! @CMAKE_USE_IFORT@ and @CMAKE_USE_GFORTRAN@ strings by the appropriate characters
  !! (spaces or comment characters).
 
use mod_kinds
use mod_global

! deal with the IFPORT module that is part of ifort on Windows
@CMAKE_USE_IFORT@use IFPORT 

IMPLICIT NONE 

contains

! there are two sets of routines for the chdir, rename, system, and hostnm commands;
! in the code, they should have the system_ prefix; at configuration time, the command 
! strings in this file are replaced as needed, depending on the operating system


! here are the calls for the Windows ifort compiler:
@CMAKE_USE_IFORT@!--------------------------------------------------------------------------
@CMAKE_USE_IFORT@function system_chdir(c) result(status)
@CMAKE_USE_IFORT@!DEC$ ATTRIBUTES DLLEXPORT :: system_chdir
@CMAKE_USE_IFORT@!! author: MDG
@CMAKE_USE_IFORT@!! version: 1.0
@CMAKE_USE_IFORT@!! date: 08/29/23 
@CMAKE_USE_IFORT@!!
@CMAKE_USE_IFORT@!! change directory wrapper for ifort
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@implicit none
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@character(*), INTENT(IN)  :: c 
@CMAKE_USE_IFORT@integer(kind=irg)         :: status
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@status = chdir(c)
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@end function system_chdir
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@!--------------------------------------------------------------------------
@CMAKE_USE_IFORT@function system_hostnm(c) result(status)
@CMAKE_USE_IFORT@!DEC$ ATTRIBUTES DLLEXPORT :: system_hostnm
@CMAKE_USE_IFORT@!! author: MDG
@CMAKE_USE_IFORT@!! version: 1.0
@CMAKE_USE_IFORT@!! date: 08/29/23 
@CMAKE_USE_IFORT@!!
@CMAKE_USE_IFORT@!! get hostname wrapper for ifort
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@implicit none
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@character(*), INTENT(INOUT)   :: c 
@CMAKE_USE_IFORT@integer(kind=irg)             :: status
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@status = hostnm(c)
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@end function system_hostnm
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@!--------------------------------------------------------------------------
@CMAKE_USE_IFORT@function system_rename(c1, c2) result(status)
@CMAKE_USE_IFORT@!DEC$ ATTRIBUTES DLLEXPORT :: system_rename
@CMAKE_USE_IFORT@!! author: MDG
@CMAKE_USE_IFORT@!! version: 1.0
@CMAKE_USE_IFORT@!! date: 08/29/23 
@CMAKE_USE_IFORT@!!
@CMAKE_USE_IFORT@!! file rename wrapper for ifort
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@implicit none
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@character(*), INTENT(IN)      :: c1
@CMAKE_USE_IFORT@character(*), INTENT(INOUT)   :: c2 
@CMAKE_USE_IFORT@integer(kind=irg)             :: status
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@status = rename(c1, c2)
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@end function system_rename
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@!--------------------------------------------------------------------------
@CMAKE_USE_IFORT@function system_system(c) result(status)
@CMAKE_USE_IFORT@!DEC$ ATTRIBUTES DLLEXPORT :: system_system
@CMAKE_USE_IFORT@!! author: MDG
@CMAKE_USE_IFORT@!! version: 1.0
@CMAKE_USE_IFORT@!! date: 08/29/23 
@CMAKE_USE_IFORT@!!
@CMAKE_USE_IFORT@!! system call wrapper for ifort
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@implicit none
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@character(*), INTENT(IN)      :: c
@CMAKE_USE_IFORT@integer(kind=irg)             :: status
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@status = system(c)
@CMAKE_USE_IFORT@
@CMAKE_USE_IFORT@end function system_system

! and here are the equivalent calls for the gfortran compiler:
@CMAKE_USE_GFORTRAN@!--------------------------------------------------------------------------
@CMAKE_USE_GFORTRAN@function system_chdir(c) result(status)
@CMAKE_USE_GFORTRAN@!DEC$ ATTRIBUTES DLLEXPORT :: system_chdir
@CMAKE_USE_GFORTRAN@!! author: MDG
@CMAKE_USE_GFORTRAN@!! version: 1.0
@CMAKE_USE_GFORTRAN@!! date: 08/29/23 
@CMAKE_USE_GFORTRAN@!!
@CMAKE_USE_GFORTRAN@!! change directory wrapper for ifort
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@implicit none
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@character(*), INTENT(IN)  :: c 
@CMAKE_USE_GFORTRAN@integer(kind=irg)         :: status
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@call chdir(c, status)
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@end function system_chdir
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@!--------------------------------------------------------------------------
@CMAKE_USE_GFORTRAN@function system_hostnm(c) result(status)
@CMAKE_USE_GFORTRAN@!DEC$ ATTRIBUTES DLLEXPORT :: system_hostnm
@CMAKE_USE_GFORTRAN@!! author: MDG
@CMAKE_USE_GFORTRAN@!! version: 1.0
@CMAKE_USE_GFORTRAN@!! date: 08/29/23 
@CMAKE_USE_GFORTRAN@!!
@CMAKE_USE_GFORTRAN@!! get hostname wrapper for ifort
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@implicit none
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@character(*), INTENT(INOUT)   :: c 
@CMAKE_USE_GFORTRAN@integer(kind=irg)             :: status
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@call hostnm(c, status)
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@end function system_hostnm
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@!--------------------------------------------------------------------------
@CMAKE_USE_GFORTRAN@function system_rename(c1, c2) result(status)
@CMAKE_USE_GFORTRAN@!DEC$ ATTRIBUTES DLLEXPORT :: system_rename
@CMAKE_USE_GFORTRAN@!! author: MDG
@CMAKE_USE_GFORTRAN@!! version: 1.0
@CMAKE_USE_GFORTRAN@!! date: 08/29/23 
@CMAKE_USE_GFORTRAN@!!
@CMAKE_USE_GFORTRAN@!! file rename wrapper for ifort
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@implicit none
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@character(*), INTENT(IN)      :: c1
@CMAKE_USE_GFORTRAN@character(*), INTENT(INOUT)   :: c2 
@CMAKE_USE_GFORTRAN@integer(kind=irg)             :: status
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@call rename(c1, c2, status)
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@end function system_rename
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@!--------------------------------------------------------------------------
@CMAKE_USE_GFORTRAN@function system_system(c) result(status)
@CMAKE_USE_GFORTRAN@!DEC$ ATTRIBUTES DLLEXPORT :: system_system
@CMAKE_USE_GFORTRAN@!! author: MDG
@CMAKE_USE_GFORTRAN@!! version: 1.0
@CMAKE_USE_GFORTRAN@!! date: 08/29/23 
@CMAKE_USE_GFORTRAN@!!
@CMAKE_USE_GFORTRAN@!! system call wrapper for ifort
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@implicit none
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@character(*), INTENT(IN)      :: c
@CMAKE_USE_GFORTRAN@integer(kind=irg)             :: status
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@call system(c, status)
@CMAKE_USE_GFORTRAN@
@CMAKE_USE_GFORTRAN@end function system_system
@CMAKE_USE_GFORTRAN@

end module mod_platformsupport