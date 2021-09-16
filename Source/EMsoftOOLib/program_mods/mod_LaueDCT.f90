! ###################################################################
! Copyright (c) 2013-2021, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_LaueDCT
  !! author: MDG 
  !! version: 1.0 
  !! date: 05/27/21
  !!
  !! class definition for the EMLaueDCT program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMLaueDCT program [this was the EMLaueSlit program for single crystals before EMsoft v. 6.0]
type, public :: LaueDCTNameListType
        real(kind=dbl)          :: Lw               ! slit width (mm)
        real(kind=dbl)          :: Lh               ! slit height (mm)
        real(kind=dbl)          :: Lx               ! distance front face of slit to divergent x-ray source (mm)
        real(kind=dbl)          :: Ly               ! slit center x position (mm)
        real(kind=dbl)          :: Lz               ! slit center y position (mm)
        real(kind=dbl)          :: VoltageH         ! highest tube voltage     
        real(kind=dbl)          :: VoltageL         ! lowest tube voltage     
        real(kind=dbl)          :: Sx               ! distance from source to samplefront (mm)
        real(kind=dbl)          :: sampletodetector ! distance sample front to detector face (mm)
        real(kind=dbl)          :: ps               ! detector pixel size (mm)
        integer(kind=irg)       :: Ny               ! number of detector pixels horizontally
        integer(kind=irg)       :: Nz               ! number of detector pixels vertically
        integer(kind=irg)       :: sampleNrotations ! number of rotation steps around sample axis; step size = 2pi/sampleNrotations
        real(kind=dbl)          :: Dx               ! detector pattern center x coordinate  [mm]
        real(kind=dbl)          :: Dy               ! detector pattern center y coordinate  [mm]
        real(kind=dbl)          :: Dz               ! detector pattern center z coordinate  [mm]
        real(kind=dbl)          :: vs               ! size of the voxels that make up the sample (mm)
        real(kind=dbl)          :: absl             ! sample absorption length [mm]
        real(kind=dbl)          :: beamstopatf      ! beam stop attenuation factor
        real(kind=dbl)          :: beamstopwidth    ! beam stop width [mm]
        real(kind=dbl)          :: beamstopheight   ! beam stop height [mm]
        real(kind=sgl)          :: spotw
        real(kind=sgl)          :: gammavalue
        real(kind=dbl)          :: intcutoffratio
        real(kind=dbl)          :: samplescalefactor! scale factor from DREAM.3D units to mm
        real(kind=dbl)          :: samplingstepsize ! step size along the x-ray traces [mm]
        integer(kind=irg)       :: BPx
        integer(kind=irg)       :: nthreads
        logical                 :: binarize
        character(fnlen)        :: backprojection
        character(fnlen)        :: tiffprefix
        character(fnlen)        :: hdfname
        character(fnlen)        :: xtalname
        character(fnlen)        :: DREAM3Dfilename
        character(fnlen)        :: EulerAnglesHDFpath(10)
        character(fnlen)        :: FeatureIDsHDFpath(10)
end type LaueDCTNameListType


type, private :: samplinglisttype 
  real(kind=dbl)                        :: pos(3)
  real(kind=dbl)                        :: ray(3)
  real(kind=dbl)                        :: front(3)
  real(kind=dbl)                        :: back(3)
  real(kind=dbl)                        :: voxelvolume
  integer(kind=irg)                     :: pixely 
  integer(kind=irg)                     :: pixelz 
  type(samplinglisttype), pointer       :: next
end type samplinglisttype

! class definition
type, public :: LaueDCT_T
private 
  character(fnlen)            :: nmldeffile = 'EMLaueDCT.nml'
  type(LaueDCTNameListType)   :: lnl 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: writeHDFNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: LaueDCT_
  procedure, pass(self) :: get_Lw_
  procedure, pass(self) :: get_Lh_
  procedure, pass(self) :: get_Lx_
  procedure, pass(self) :: get_Ly_
  procedure, pass(self) :: get_Lz_
  procedure, pass(self) :: get_VoltageH_
  procedure, pass(self) :: get_VoltageL_
  procedure, pass(self) :: get_Sx_
  procedure, pass(self) :: get_sampletodetector_
  procedure, pass(self) :: get_ps_
  procedure, pass(self) :: get_Ny_
  procedure, pass(self) :: get_Nz_
  procedure, pass(self) :: get_sampleNrotations_
  procedure, pass(self) :: get_Dx_
  procedure, pass(self) :: get_Dy_
  procedure, pass(self) :: get_Dz_
  procedure, pass(self) :: get_vs_
  procedure, pass(self) :: get_absl_
  procedure, pass(self) :: get_beamstopatf_
  procedure, pass(self) :: get_beamstopwidth_
  procedure, pass(self) :: get_beamstopheight_
  procedure, pass(self) :: get_spotw_
  procedure, pass(self) :: get_gammavalue_
  procedure, pass(self) :: get_intcutoffratio_
  procedure, pass(self) :: get_samplescalefactor_
  procedure, pass(self) :: get_samplingstepsize_
  procedure, pass(self) :: get_BPx_
  procedure, pass(self) :: get_nthreads_
  procedure, pass(self) :: get_binarize_
  procedure, pass(self) :: get_backprojection_
  procedure, pass(self) :: get_tiffprefix_
  procedure, pass(self) :: get_hdfname_
  procedure, pass(self) :: get_xtalname_
  procedure, pass(self) :: get_DREAM3Dfilename_
  procedure, pass(self) :: set_Lw_
  procedure, pass(self) :: set_Lh_
  procedure, pass(self) :: set_Lx_
  procedure, pass(self) :: set_Ly_
  procedure, pass(self) :: set_Lz_
  procedure, pass(self) :: set_VoltageH_
  procedure, pass(self) :: set_VoltageL_
  procedure, pass(self) :: set_Sx_
  procedure, pass(self) :: set_sampletodetector_
  procedure, pass(self) :: set_ps_
  procedure, pass(self) :: set_Ny_
  procedure, pass(self) :: set_Nz_
  procedure, pass(self) :: set_sampleNrotations_
  procedure, pass(self) :: set_Dx_
  procedure, pass(self) :: set_Dy_
  procedure, pass(self) :: set_Dz_
  procedure, pass(self) :: set_vs_
  procedure, pass(self) :: set_absl_
  procedure, pass(self) :: set_beamstopatf_
  procedure, pass(self) :: set_beamstopwidth_
  procedure, pass(self) :: set_beamstopheight_
  procedure, pass(self) :: set_spotw_
  procedure, pass(self) :: set_gammavalue_
  procedure, pass(self) :: set_intcutoffratio_
  procedure, pass(self) :: set_samplescalefactor_
  procedure, pass(self) :: set_samplingstepsize_
  procedure, pass(self) :: set_BPx_
  procedure, pass(self) :: set_nthreads_
  procedure, pass(self) :: set_binarize_
  procedure, pass(self) :: set_backprojection_
  procedure, pass(self) :: set_tiffprefix_
  procedure, pass(self) :: set_hdfname_
  procedure, pass(self) :: set_xtalname_
  procedure, pass(self) :: set_DREAM3Dfilename_

  generic, public :: getNameList => getNameList_
  generic, public :: writeHDFNameList => writeHDFNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: LaueDCT => LaueDCT_
  generic, public :: get_Lw => get_Lw_
  generic, public :: get_Lh => get_Lh_
  generic, public :: get_Lx => get_Lx_
  generic, public :: get_Ly => get_Ly_
  generic, public :: get_Lz => get_Lz_
  generic, public :: get_VoltageH => get_VoltageH_
  generic, public :: get_VoltageL => get_VoltageL_
  generic, public :: get_Sx => get_Sx_
  generic, public :: get_sampletodetector => get_sampletodetector_
  generic, public :: get_ps => get_ps_
  generic, public :: get_Ny => get_Ny_
  generic, public :: get_Nz => get_Nz_
  generic, public :: get_sampleNrotations => get_sampleNrotations_
  generic, public :: get_Dx => get_Dx_
  generic, public :: get_Dy => get_Dy_
  generic, public :: get_Dz => get_Dz_
  generic, public :: get_vs => get_vs_
  generic, public :: get_absl => get_absl_
  generic, public :: get_beamstopatf => get_beamstopatf_
  generic, public :: get_beamstopwidth => get_beamstopwidth_
  generic, public :: get_beamstopheight => get_beamstopheight_
  generic, public :: get_spotw => get_spotw_
  generic, public :: get_gammavalue => get_gammavalue_
  generic, public :: get_intcutoffratio => get_intcutoffratio_
  generic, public :: get_samplescalefactor => get_samplescalefactor_
  generic, public :: get_samplingstepsize => get_samplingstepsize_
  generic, public :: get_BPx => get_BPx_
  generic, public :: get_nthreads => get_nthreads_
  generic, public :: get_binarize => get_binarize_
  generic, public :: get_backprojection => get_backprojection_
  generic, public :: get_tiffprefix => get_tiffprefix_
  generic, public :: get_hdfname => get_hdfname_
  generic, public :: get_xtalname => get_xtalname_
  generic, public :: get_DREAM3Dfilename => get_DREAM3Dfilename_
  generic, public :: set_Lw => set_Lw_
  generic, public :: set_Lh => set_Lh_
  generic, public :: set_Lx => set_Lx_
  generic, public :: set_Ly => set_Ly_
  generic, public :: set_Lz => set_Lz_
  generic, public :: set_VoltageH => set_VoltageH_
  generic, public :: set_VoltageL => set_VoltageL_
  generic, public :: set_Sx => set_Sx_
  generic, public :: set_sampletodetector => set_sampletodetector_
  generic, public :: set_ps => set_ps_
  generic, public :: set_Ny => set_Ny_
  generic, public :: set_Nz => set_Nz_
  generic, public :: set_sampleNrotations => set_sampleNrotations_
  generic, public :: set_Dx => set_Dx_
  generic, public :: set_Dy => set_Dy_
  generic, public :: set_Dz => set_Dz_
  generic, public :: set_vs => set_vs_
  generic, public :: set_absl => set_absl_
  generic, public :: set_beamstopatf => set_beamstopatf_
  generic, public :: set_beamstopwidth => set_beamstopwidth_
  generic, public :: set_beamstopheight => set_beamstopheight_
  generic, public :: set_spotw => set_spotw_
  generic, public :: set_gammavalue => set_gammavalue_
  generic, public :: set_intcutoffratio => set_intcutoffratio_
  generic, public :: set_samplescalefactor => set_samplescalefactor_
  generic, public :: set_samplingstepsize => set_samplingstepsize_
  generic, public :: set_BPx => set_BPx_
  generic, public :: set_nthreads => set_nthreads_
  generic, public :: set_binarize => set_binarize_
  generic, public :: set_backprojection => set_backprojection_
  generic, public :: set_tiffprefix => set_tiffprefix_
  generic, public :: set_hdfname => set_hdfname_
  generic, public :: set_xtalname => set_xtalname_
  generic, public :: set_DREAM3Dfilename => set_DREAM3Dfilename_
end type LaueDCT_T

! the constructor routine for this class 
interface LaueDCT_T
  module procedure LaueDCT_constructor
end interface LaueDCT_T

contains

!--------------------------------------------------------------------------
type(LaueDCT_T) function LaueDCT_constructor( nmlfile ) result(LaueDCT)
!DEC$ ATTRIBUTES DLLEXPORT :: LaueDCT_constructor 
!! author: MDG 
!! version: 1.0 
!! date: 05/27/21
!!
!! constructor for the LaueDCT_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call LaueDCT%readNameList(nmlfile)

end function LaueDCT_constructor

!--------------------------------------------------------------------------
subroutine LaueDCT_destructor(self) 
!DEC$ ATTRIBUTES DLLEXPORT :: LaueDCT_destructor
!! author: MDG 
!! version: 1.0 
!! date: 05/27/21
!!
!! destructor for the LaueDCT_T Class
 
IMPLICIT NONE

type(LaueDCT_T), INTENT(INOUT)  :: self 

call reportDestructor('LaueDCT_T')

end subroutine LaueDCT_destructor

!--------------------------------------------------------------------------
function get_Lw_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Lw_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get Lw from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%Lw

end function get_Lw_

!--------------------------------------------------------------------------
subroutine set_Lw_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Lw_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set Lw in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%Lw = inp

end subroutine set_Lw_

!--------------------------------------------------------------------------
function get_Lh_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Lh_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get Lh from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%Lh

end function get_Lh_

!--------------------------------------------------------------------------
subroutine set_Lh_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Lh_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set Lh in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%Lh = inp

end subroutine set_Lh_

!--------------------------------------------------------------------------
function get_Lx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Lx_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get Lx from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%Lx

end function get_Lx_

!--------------------------------------------------------------------------
subroutine set_Lx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Lx_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set Lx in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%Lx = inp

end subroutine set_Lx_

!--------------------------------------------------------------------------
function get_Ly_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Ly_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get Ly from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%Ly

end function get_Ly_

!--------------------------------------------------------------------------
subroutine set_Ly_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Ly_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set Ly in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%Ly = inp

end subroutine set_Ly_

!--------------------------------------------------------------------------
function get_Lz_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Lz_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get Lz from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%Lz

end function get_Lz_

!--------------------------------------------------------------------------
subroutine set_Lz_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Lz_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set Lz in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%Lz = inp

end subroutine set_Lz_

!--------------------------------------------------------------------------
function get_VoltageH_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_VoltageH_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get VoltageH from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%VoltageH

end function get_VoltageH_

!--------------------------------------------------------------------------
subroutine set_VoltageH_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_VoltageH_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set VoltageH in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%VoltageH = inp

end subroutine set_VoltageH_

!--------------------------------------------------------------------------
function get_VoltageL_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_VoltageL_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get VoltageL from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%VoltageL

end function get_VoltageL_

!--------------------------------------------------------------------------
subroutine set_VoltageL_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_VoltageL_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set VoltageL in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%VoltageL = inp

end subroutine set_VoltageL_

!--------------------------------------------------------------------------
function get_Sx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Sx_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get Sx from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%Sx

end function get_Sx_

!--------------------------------------------------------------------------
subroutine set_Sx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Sx_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set Sx in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%Sx = inp

end subroutine set_Sx_

!--------------------------------------------------------------------------
function get_sampletodetector_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_sampletodetector_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get sampletodetector from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%sampletodetector

end function get_sampletodetector_

!--------------------------------------------------------------------------
subroutine set_sampletodetector_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_sampletodetector_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set sampletodetector in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%sampletodetector = inp

end subroutine set_sampletodetector_

!--------------------------------------------------------------------------
function get_ps_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ps_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get ps from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%ps

end function get_ps_

!--------------------------------------------------------------------------
subroutine set_ps_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_ps_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set ps in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%ps = inp

end subroutine set_ps_

!--------------------------------------------------------------------------
function get_Ny_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Ny_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get Ny from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%lnl%Ny

end function get_Ny_

!--------------------------------------------------------------------------
subroutine set_Ny_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Ny_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set Ny in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%lnl%Ny = inp

end subroutine set_Ny_

!--------------------------------------------------------------------------
function get_Nz_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Nz_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get Nz from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%lnl%Nz

end function get_Nz_

!--------------------------------------------------------------------------
subroutine set_Nz_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Nz_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set Nz in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%lnl%Nz = inp

end subroutine set_Nz_

!--------------------------------------------------------------------------
function get_sampleNrotations_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_sampleNrotations_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get sampleNrotations from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%lnl%sampleNrotations

end function get_sampleNrotations_

!--------------------------------------------------------------------------
subroutine set_sampleNrotations_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_sampleNrotations_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set sampleNrotations in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%lnl%sampleNrotations = inp

end subroutine set_sampleNrotations_

!--------------------------------------------------------------------------
function get_Dx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Dx_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get Dx from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%Dx

end function get_Dx_

!--------------------------------------------------------------------------
subroutine set_Dx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Dx_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set Dx in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%Dx = inp

end subroutine set_Dx_

!--------------------------------------------------------------------------
function get_Dy_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Dy_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get Dy from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%Dy

end function get_Dy_

!--------------------------------------------------------------------------
subroutine set_Dy_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Dy_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set Dy in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%Dy = inp

end subroutine set_Dy_

!--------------------------------------------------------------------------
function get_Dz_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_Dz_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get Dz from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%Dz

end function get_Dz_

!--------------------------------------------------------------------------
subroutine set_Dz_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_Dz_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set Dz in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%Dz = inp

end subroutine set_Dz_

!--------------------------------------------------------------------------
function get_vs_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_vs_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get vs from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%vs

end function get_vs_

!--------------------------------------------------------------------------
subroutine set_vs_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_vs_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set vs in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%vs = inp

end subroutine set_vs_

!--------------------------------------------------------------------------
function get_absl_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_absl_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get absl from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%absl

end function get_absl_

!--------------------------------------------------------------------------
subroutine set_absl_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_absl_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set absl in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%absl = inp

end subroutine set_absl_

!--------------------------------------------------------------------------
function get_beamstopatf_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_beamstopatf_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get beamstopatf from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%beamstopatf

end function get_beamstopatf_

!--------------------------------------------------------------------------
subroutine set_beamstopatf_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_beamstopatf_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set beamstopatf in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%beamstopatf = inp

end subroutine set_beamstopatf_

!--------------------------------------------------------------------------
function get_beamstopwidth_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_beamstopwidth_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get beamstopwidth from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%beamstopwidth

end function get_beamstopwidth_

!--------------------------------------------------------------------------
subroutine set_beamstopwidth_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_beamstopwidth_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set beamstopwidth in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%beamstopwidth = inp

end subroutine set_beamstopwidth_

!--------------------------------------------------------------------------
function get_beamstopheight_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_beamstopheight_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get beamstopheight from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%beamstopheight

end function get_beamstopheight_

!--------------------------------------------------------------------------
subroutine set_beamstopheight_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_beamstopheight_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set beamstopheight in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%beamstopheight = inp

end subroutine set_beamstopheight_

!--------------------------------------------------------------------------
function get_spotw_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_spotw_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get spotw from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=sgl)                      :: out

out = self%lnl%spotw

end function get_spotw_

!--------------------------------------------------------------------------
subroutine set_spotw_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_spotw_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set spotw in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)          :: inp

self%lnl%spotw = inp

end subroutine set_spotw_

!--------------------------------------------------------------------------
function get_gammavalue_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_gammavalue_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get gammavalue from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=sgl)                      :: out

out = self%lnl%gammavalue

end function get_gammavalue_

!--------------------------------------------------------------------------
subroutine set_gammavalue_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_gammavalue_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set gammavalue in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=sgl), INTENT(IN)          :: inp

self%lnl%gammavalue = inp

end subroutine set_gammavalue_

!--------------------------------------------------------------------------
function get_intcutoffratio_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_intcutoffratio_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get intcutoffratio from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%intcutoffratio

end function get_intcutoffratio_

!--------------------------------------------------------------------------
subroutine set_intcutoffratio_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_intcutoffratio_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set intcutoffratio in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%intcutoffratio = inp

end subroutine set_intcutoffratio_

!--------------------------------------------------------------------------
function get_samplescalefactor_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_samplescalefactor_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get samplescalefactor from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%samplescalefactor

end function get_samplescalefactor_

!--------------------------------------------------------------------------
subroutine set_samplescalefactor_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_samplescalefactor_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set samplescalefactor in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%samplescalefactor = inp

end subroutine set_samplescalefactor_

!--------------------------------------------------------------------------
function get_samplingstepsize_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_samplingstepsize_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get samplingstepsize from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl)                      :: out

out = self%lnl%samplingstepsize

end function get_samplingstepsize_

!--------------------------------------------------------------------------
subroutine set_samplingstepsize_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_samplingstepsize_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set samplingstepsize in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
real(kind=dbl), INTENT(IN)          :: inp

self%lnl%samplingstepsize = inp

end subroutine set_samplingstepsize_

!--------------------------------------------------------------------------
function get_BPx_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_BPx_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get BPx from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%lnl%BPx

end function get_BPx_

!--------------------------------------------------------------------------
subroutine set_BPx_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_BPx_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set BPx in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%lnl%BPx = inp

end subroutine set_BPx_

!--------------------------------------------------------------------------
function get_nthreads_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nthreads_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get nthreads from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
integer(kind=irg)                   :: out

out = self%lnl%nthreads

end function get_nthreads_

!--------------------------------------------------------------------------
subroutine set_nthreads_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_nthreads_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set nthreads in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
integer(kind=irg), INTENT(IN)       :: inp

self%lnl%nthreads = inp

end subroutine set_nthreads_

!--------------------------------------------------------------------------
function get_binarize_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_binarize_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get binarize from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
logical                             :: out

out = self%lnl%binarize

end function get_binarize_

!--------------------------------------------------------------------------
subroutine set_binarize_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_binarize_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set binarize in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
logical, INTENT(IN)                 :: inp

self%lnl%binarize = inp

end subroutine set_binarize_

!--------------------------------------------------------------------------
function get_backprojection_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_backprojection_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get backprojection from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
character(fnlen)                    :: out

out = self%lnl%backprojection

end function get_backprojection_

!--------------------------------------------------------------------------
subroutine set_backprojection_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_backprojection_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set backprojection in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)        :: inp

self%lnl%backprojection = inp

end subroutine set_backprojection_

!--------------------------------------------------------------------------
function get_tiffprefix_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_tiffprefix_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get tiffprefix from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
character(fnlen)                    :: out

out = self%lnl%tiffprefix

end function get_tiffprefix_

!--------------------------------------------------------------------------
subroutine set_tiffprefix_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_tiffprefix_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set tiffprefix in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)        :: inp

self%lnl%tiffprefix = inp

end subroutine set_tiffprefix_

!--------------------------------------------------------------------------
function get_hdfname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_hdfname_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get hdfname from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
character(fnlen)                    :: out

out = self%lnl%hdfname

end function get_hdfname_

!--------------------------------------------------------------------------
subroutine set_hdfname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_hdfname_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set hdfname in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)        :: inp

self%lnl%hdfname = inp

end subroutine set_hdfname_

!--------------------------------------------------------------------------
function get_xtalname_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_xtalname_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get xtalname from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
character(fnlen)                    :: out

out = self%lnl%xtalname

end function get_xtalname_

!--------------------------------------------------------------------------
subroutine set_xtalname_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_xtalname_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set xtalname in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)        :: inp

self%lnl%xtalname = inp

end subroutine set_xtalname_

!--------------------------------------------------------------------------
function get_DREAM3Dfilename_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: get_DREAM3Dfilename_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! get DREAM3Dfilename from the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
character(fnlen)                    :: out

out = self%lnl%DREAM3Dfilename

end function get_DREAM3Dfilename_

!--------------------------------------------------------------------------
subroutine set_DREAM3Dfilename_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: set_DREAM3Dfilename_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/2021
!!
!! set DREAM3Dfilename in the LaueDCT_T class

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)     :: self
character(fnlen), INTENT(IN)        :: inp

self%lnl%DREAM3Dfilename = inp

end subroutine set_DREAM3Dfilename_

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/21
!!
!! read the namelist from an nml file for the LaueDCT_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)      :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

real(kind=dbl)                       :: Lw               ! slit width (mm)
real(kind=dbl)                       :: Lh               ! slit height (mm)
real(kind=dbl)                       :: Lx               ! distance front face of slit to divergent x-ray source (mm)
real(kind=dbl)                       :: Ly               ! slit center x position (mm)
real(kind=dbl)                       :: Lz               ! slit center y position (mm)
real(kind=dbl)                       :: VoltageH         ! highest tube voltage     
real(kind=dbl)                       :: VoltageL         ! lowest tube voltage     
real(kind=dbl)                       :: Sx               ! distance from source to samplefront (mm)
real(kind=dbl)                       :: sampletodetector ! distance sample front to detector face (mm)
real(kind=dbl)                       :: ps               ! detector pixel size (mm)
integer(kind=irg)                    :: Ny               ! number of detector pixels horizontally
integer(kind=irg)                    :: Nz               ! number of detector pixels vertically
integer(kind=irg)                    :: sampleNrotations ! number of rotation steps around sample axis
real(kind=dbl)                       :: Dx               ! detector pattern center x coordinate  [mm]
real(kind=dbl)                       :: Dy               ! detector pattern center y coordinate  [mm]
real(kind=dbl)                       :: Dz               ! detector pattern center z coordinate  [mm]
real(kind=dbl)                       :: vs               ! size of the voxels that make up the sample (mm)
real(kind=dbl)                       :: absl             ! sample absorption length [mm]
real(kind=dbl)                       :: beamstopatf      ! beam stop attenuation factor
real(kind=dbl)                       :: beamstopwidth    ! beam stop width [mm]
real(kind=dbl)                       :: beamstopheight   ! beam stop height [mm]
real(kind=sgl)                       :: spotw
real(kind=sgl)                       :: gammavalue
real(kind=dbl)                       :: intcutoffratio
real(kind=dbl)                       :: samplescalefactor
real(kind=dbl)                       :: samplingstepsize
integer(kind=irg)                    :: BPx
integer(kind=irg)                    :: nthreads
logical                              :: binarize
character(fnlen)                     :: backprojection
character(fnlen)                     :: tiffprefix
character(fnlen)                     :: hdfname
character(fnlen)                     :: xtalname
character(fnlen)                     :: DREAM3Dfilename
character(fnlen)                     :: EulerAnglesHDFpath(10)
character(fnlen)                     :: FeatureIDsHDFpath(10)

! define the IO namelist to facilitate passing variables to the program.
namelist  / LaueDCTData / Lw,Lh,Lx,Ly,Lz,VoltageH,VoltageL,Sx,sampletodetector, beamstopwidth, beamstopheight,&
                          ps,Ny,Nz,sampleNrotations,Dx,Dy,Dz,vs,absl, binarize, &
                          beamstopatf,spotw,BPx,nthreads,backprojection, intcutoffratio, samplescalefactor, &
                          tiffprefix,hdfname,xtalname, gammavalue, &
                          samplingstepsize, DREAM3Dfilename, EulerAnglesHDFpath, FeatureIDsHDFpath

Lw               = 2.D0    ! slit width (mm)
Lh               = 2.D0    ! slit height (mm)
Lx               = 100.D0  ! distance front face of slit to divergent x-ray source (mm)
Ly               = 0.D0    ! slit center x position (mm)
Lz               = 0.D0    ! slit center y position (mm)
VoltageH         = 60.D0   ! highest tube voltage     
VoltageL         = 40.D0   ! lowest tube voltage     
Sx               = 120.D0  ! distance from source to samplefront (mm)
sampletodetector = 120.D0  ! distance sample front to detector face (mm)
ps               = 0.254D0 ! pixel width (mm)
Ny               = 960     ! number of pixels horizontally
Nz               = 780     ! number of pixels vertically
sampleNrotations = 180     ! number of rotation steps in a 2pi rotation
Dx               = 0.D0    ! pattern center x coordinate 
Dy               = 0.D0    ! pattern center y coordinate 
Dz               = 0.D0    ! pattern center z coordinate 
vs               = 0.10D0  ! size of the voxels that make up the sample (mm)
absl             = 0.5D0   ! absorption length (mm)
beamstopatf      = 0.1D0   ! beam stop attenuation factor
beamstopwidth    = 2.0D0   ! beam stop width [mm]
beamstopheight   = 2.0D0   ! beam stop height [mm]
nthreads         = 1       ! number of parallel threads for pattern computation
BPx              = 300     ! semi-edge length for back projection square Lambert maps
spotw            = 0.1     ! spot size weight factor (1/(2*sigma^2))
gammavalue       = 1.0     ! scaling factor for gamma intensity scaling
intcutoffratio   = 0.0001D0! intensity ratio cut off
samplescalefactor= 10.D0   ! 10 microns per DREAM.3D unit length
samplingstepsize = 0.005D0 ! in mm, step size along ray
binarize         = .FALSE.
backprojection   = 'No'    ! 'Yes' or 'No'; adds backprojections to output file
tiffprefix       = 'undefined'  ! prefix for tiff output files with individual patterns
xtalname         = 'undefined'  ! structure file name
hdfname          = 'undefined'  ! HDF output file name
DREAM3Dfilename  = 'undefined'  ! DREAM.3D input file name (only for polycrystal)
EulerAnglesHDFpath = (/ '', '', '', '', '', '', '', '', '', '' /) ! path to EulerAngles 
FeatureIDsHDFpath  = (/ '', '', '', '', '', '', '', '', '', '' /)  ! path to FeatureIDs 

EulerAnglesHDFpath(1) = 'DataContainers'
EulerAnglesHDFpath(3) = 'CellData'
EulerAnglesHDFpath(4) = 'EulerAngles'

FeatureIDsHDFpath(1) = 'DataContainers'
FeatureIDsHDFpath(3) = 'CellData'
FeatureIDsHDFpath(4) = 'FeatureIDs'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=LaueDCTData)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call Message%printError('readNameList:',' crystal structure file name is undefined in '//nmlfile)
 end if
 if (trim(hdfname).eq.'undefined') then
  call Message%printError('readNameList:',' output file name is undefined in '//nmlfile)
 end if
end if

self%lnl%Lw = Lw               
self%lnl%Lh = Lh               
self%lnl%Lx = Lx               
self%lnl%Ly = Ly               
self%lnl%Lz = Lz               
self%lnl%VoltageH = VoltageH         
self%lnl%VoltageL = VoltageL         
self%lnl%Sx = Sx               
self%lnl%sampletodetector = sampletodetector 
self%lnl%ps = ps               
self%lnl%Ny = Ny               
self%lnl%Nz = Nz               
self%lnl%sampleNrotations = sampleNrotations               
self%lnl%Dx = Dx               
self%lnl%Dy = Dy               
self%lnl%Dz = Dz               
self%lnl%vs = vs               
self%lnl%absl = absl             
self%lnl%beamstopatf = beamstopatf
self%lnl%beamstopwidth = beamstopwidth
self%lnl%beamstopheight = beamstopheight
self%lnl%spotw = spotw
self%lnl%BPx = BPx
self%lnl%nthreads = nthreads
self%lnl%intcutoffratio = intcutoffratio
self%lnl%samplescalefactor = samplescalefactor
self%lnl%samplingstepsize = samplingstepsize
self%lnl%tiffprefix = tiffprefix
self%lnl%hdfname = hdfname
self%lnl%xtalname = xtalname
self%lnl%binarize = binarize
self%lnl%DREAM3Dfilename = DREAM3Dfilename
self%lnl%EulerAnglesHDFpath = EulerAnglesHDFpath
self%lnl%FeatureIDsHDFpath = FeatureIDsHDFpath

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/21
!!
!! pass the namelist for the LaueDCT_T Class to the calling program

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)          :: self
type(LaueDCTNameListType)                :: nml

nml = self%lnl

end function getNameList_

!--------------------------------------------------------------------------
recursive subroutine writeHDFNameList_(self, HDF, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: writeHDFNameList_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/21
!!
!! write namelist to HDF file

use mod_HDFsupport
use mod_HDFnames
use stringconstants 

use ISO_C_BINDING

IMPLICIT NONE

class(LaueDCT_T), INTENT(INOUT)         :: self 
type(HDF_T), INTENT(INOUT)              :: HDF
type(HDFnames_T), INTENT(INOUT)         :: HDFnames

integer(kind=irg),parameter             :: n_int = 5, n_real = 2, n_dbl = 20
integer(kind=irg)                       :: hdferr,  io_int(n_int), i
real(kind=sgl)                          :: io_real(n_real)
real(kind=dbl)                          :: io_dbl(n_dbl)   
character(20)                           :: intlist(n_int), reallist(n_real), dbllist(n_dbl)
character(fnlen)                        :: dataset, sval(1),groupname
character(fnlen,kind=c_char)            :: line2(1)
character(fnlen,kind=c_char)            :: line10(10)

associate( lnl => self%lnl )

! create the group for this namelist
hdferr = HDF%createGroup(HDFnames%get_NMLlist())

! write all the single integers
io_int = (/ lnl%Ny, lnl%Nz, lnl%nthreads, lnl%BPx, lnl%sampleNrotations /)
intlist(1) = 'Ny'
intlist(2) = 'Nz'
intlist(3) = 'nthreads'
intlist(4) = 'BPx'
intlist(5) = 'sampleNrotations'
call HDF%writeNMLintegers(io_int, intlist, n_int)

! write all the single reals
io_real = (/ lnl%spotw, lnl%gammavalue /)
reallist(1) = 'spotw'
reallist(2) = 'gammavalue'
call HDF%writeNMLreals(io_real, reallist, n_real)

! write all the single reals
io_dbl = (/ lnl%Lw, lnl%Lh, lnl%Lx, lnl%Ly, lnl%Lz, lnl%VoltageH, lnl%VoltageL, lnl%Sx, &
            lnl%sampletodetector, lnl%ps, lnl%Dy, &
            lnl%Dz, lnl%vs, lnl%absl, lnl%beamstopatf, lnl%beamstopwidth, lnl%beamstopheight, lnl%intcutoffratio, &
            lnl%samplescalefactor, lnl%samplingstepsize /)
dbllist(1) = 'Lw'
dbllist(2) = 'Lh'
dbllist(3) = 'Lx'
dbllist(4) = 'Ly'
dbllist(5) = 'Lz'
dbllist(6) = 'VoltageH'
dbllist(7) = 'VoltageL'
dbllist(8) = 'Sx'
dbllist(9) = 'sampletodetector'
dbllist(10) = 'ps'
dbllist(11) = 'Dy'
dbllist(12) = 'Dz'
dbllist(13) = 'vs'
dbllist(14) = 'absl'
dbllist(15) = 'beamstopatf'
dbllist(16) = 'beamstopwidth'
dbllist(17) = 'beamstopheight'
dbllist(18) = 'intcutoffratio'
dbllist(19) = 'samplescalefactor'
dbllist(20) = 'samplingstepsize'
call HDF%writeNMLdbles(io_dbl, dbllist, n_dbl)

! write all the strings
dataset = SC_xtalname
line2(1) = trim(lnl%xtalname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create xtalname dataset', hdferr)

dataset = 'hdfname'
line2(1) = trim(lnl%hdfname)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create hdfname dataset', hdferr)

dataset = 'tiffprefix'
line2(1) = trim(lnl%tiffprefix)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create tiffprefix dataset', hdferr)

dataset = 'DREAM3Dfilename'
line2(1) = trim(lnl%DREAM3Dfilename)
hdferr = HDF%writeDatasetStringArray(dataset, line2, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create DREAM3Dfilename dataset', hdferr)

dataset = 'EulerAnglesHDFpath'
do i=1,10 
  line10(i) = trim(lnl%EulerAnglesHDFpath(i))
end do 
hdferr = HDF%writeDatasetStringArray(dataset, line10, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create EulerAnglesHDFpath dataset', hdferr)

dataset = 'FeatureIDsHDFpath'
do i=1,10 
  line10(i) = trim(lnl%FeatureIDsHDFpath(i))
end do 
hdferr = HDF%writeDatasetStringArray(dataset, line10, 1)
if (hdferr.ne.0) call HDF%error_check('writeHDFNameList: unable to create FeatureIDsHDFpath dataset', hdferr)

! and pop this group off the stack
call HDF%pop()

end associate

end subroutine writeHDFNameList_

!--------------------------------------------------------------------------
subroutine LaueDCT_(self, EMsoft, progname, HDFnames)
!DEC$ ATTRIBUTES DLLEXPORT :: LaueDCT_
!! author: MDG 
!! version: 1.0 
!! date: 05/27/21
!!
!! perform the computations
!!
!! This routine uses some of the approaches in the following paper:
!! "A flexible and standalone forward simulation model for laboratory X-ray diffraction contrast tomography"
!! H. Fang, D. Juul Jensen and Y. Zhang, Acta Cryst (2020) A76, 652–663
!!
!! but the microstructure sampling algorithm is different in the sense that a dense grid of points is used
!! with all points lying on the rays connecting the detector pixels with the source.
!!
!! As we develop this routine further we will add more of the contributions to the intensity to make 
!! the simulation more realistic  
!!
!!

use mod_EMsoft
use mod_initializers
use mod_xrd
use mod_symmetry
use mod_crystallography
use mod_gvectors
use mod_kvectors
use mod_diffraction
use mod_rotations
use mod_quaternions
use mod_so3
use mod_LaueSupport
use mod_io
use mod_math
use mod_timing
use mod_Lambert
use HDF5
use mod_HDFsupport
use mod_HDFnames
use mod_memory
use ISO_C_BINDING
use omp_lib
use mod_notifications
use stringconstants
use mod_image
use mod_Laue
use mod_DREAM3D
use, intrinsic :: iso_fortran_env

IMPLICIT NONE 

class(LaueDCT_T), INTENT(INOUT)            :: self
type(EMsoft_T), INTENT(INOUT)              :: EMsoft
character(fnlen), INTENT(INOUT)            :: progname 
type(HDFnames_T), INTENT(INOUT)            :: HDFnames 

type(Cell_T)                               :: cell
type(Timing_T)                             :: timer
type(IO_T)                                 :: Message
type(HDF_T)                                :: HDF
type(SpaceGroup_T)                         :: SG
type(Diffraction_T)                        :: Diff
type(DynType)                              :: Dyn
type(memory_T)                             :: mem, memth
type(so3_T)                                :: SO
type(QuaternionArray_T)                    :: qAR, orlist
type(Quaternion_T)                         :: quat
type(LaueReflist_T)                        :: LaueReflist
type(q_T)                                  :: qu
type(microstructure)                       :: microstr
type(samplinglisttype), pointer            :: samplinglist, stmp, stmpb

integer(kind=irg)                          :: numangles, numbatches, remainder, ii, jj, pid, tickstart, shadow(2,4), np

integer(kind=irg)                          :: i, j, icnt, numvox, hdferr, npx, npy, refcnt, io_int(4), Lstart, bsw, bsh, gcnt, &
                                              g(3), gr(3), rf, NUMTHREADS, TID, BPnpx, BPnpy, m, betamin, betamax, NNy, NNz
real(kind=sgl)                             :: l, kouter, kinner, tstart, tstop, mi, ma, lambdamin, lambdamax, kv(3), &
                                              scl, kv2(3), shortg, info, gg, mps, dmin, intfactor, io_real(3), fpar(20) 
real(kind=sgl),allocatable                 :: pp(:,:), pattern(:,:), patternsum(:,:,:), bppatterns(:,:,:), bp(:,:)
real(kind=dbl)                             :: qq(4), pre, dtp
logical,allocatable                        :: line(:)

logical                                    :: verbose, f_exists, g_exists, insert=.TRUE., overwrite=.TRUE.

type(Laue_grow_list),pointer               :: reflist, rltmp          

real(kind=dbl)                             :: dt, pixel(3), k, disc, tt, tm, tp, newt, th, dpar(20), previous 
real(kind=dbl),allocatable                 :: slist(:,:), rotated_slist(:,:), rays(:,:), frontpts(:,:,:), backpts(:,:,:), voxvol(:)
integer(kind=irg)                          :: it, nt, npoints, ipar(20)
integer(kind=irg),allocatable              :: spppix(:,:), nst(:,:,:)


character(fnlen)                           :: hdfname, groupname, datagroupname, attributename, dataset, fname, &
                                              TIFF_filename, BPmode
character(11)                              :: dstr
character(15)                              :: tstrb
character(15)                              :: tstre
character(4)                               :: pnum
character(fnlen)                           :: HDF_FileVersion
integer(HSIZE_T)                           :: dims3(3), cnt3(3), offset3(3)
character(fnlen,kind=c_char)               :: line2(1)

! declare variables for use in object oriented image module
integer                                    :: iostat
character(len=128)                         :: iomsg
logical                                    :: isInteger
type(image_t)                              :: im
integer(int8)                              :: i8 (3,4)
integer(int8), allocatable                 :: TIFF_image(:,:)

! Legendre lattitude arrays
real(kind=dbl),allocatable                 :: LegendreArray(:), upd(:), diagonal(:)

! parameters for the slit model 
real(kind=dbl)                             :: ds, dsvec(3), d0, d0vec(3), d, dvec(3), t, slitc(3), sw, sh, slitcorners(3,4), &
                                              scuvec(3,4), scdet(3,4), scsbp(3,4), samplecenter(3), dx, dy, dz, kuvec(3)
integer(kind=irg)                          :: minvy, minvz, maxvy, maxvz, numvx, ix, iy, iz  
real(kind=sgl),allocatable                 :: kvecs(:,:), kvox(:,:), kinpre(:) 

! initialize the HDF class 
call openFortranHDFInterface()
HDF = HDF_T()

! simplify the notation a little
associate( lnl => self%lnl )

! initialize the timing routines
timer = Timing_T()
tstrb = timer%getTimeString()
call timer%Time_tick(1)

! initialize the memory class 
mem = memory_T()

! rotations in double precision
call setRotationPrecision('d')

! read a DREAM.3D microstructure file or a list of orientations for a single crystal ?
call ReadDREAM3Dfile(EMsoft, lnl%DREAM3Dfilename, microstr, lnl%EulerAnglesHDFpath, lnl%FeatureIDsHDFpath)
microstr%samplescalefactor = lnl%samplescalefactor
call Message%printMessage(' ------------------- ')
io_int(1:3) = microstr%dimensions
call Message%writeValue('  # of voxels           :', io_int, 3)
io_real(1:3) = microstr%origin
call Message%writeValue('  Origin                :', io_real, 3)
io_real(1:3) = microstr%gridspacing
call Message%writeValue('  Spacing               :', io_real, 3)
gr = microstr%Quaternions%getQnumber()
io_int(1:4) = (/ 4, gr(1), gr(2), gr(3) /)
call Message%writeValue('  Quaternion array shape:', io_int, 4)
io_int(1:3) = shape(microstr%FeatureIDs)
call Message%writeValue('  FeatureIDs array shape:', io_int, 3)
io_int(1) = maxval(microstr%FeatureIDs)
call Message%writeValue('  # of grains           :', io_int, 1)
numangles = lnl%sampleNrotations 

! the sample thickness is the number of non-zero voxels along the x-direction at the center of the y-z plane
call mem%alloc(line, (/ int(microstr%dimensions(1)) /), 'line')
line = ( microstr%FeatureIDs(:, microstr%dimensions(2)/2, microstr%dimensions(3)/2) .ne. 0)
t = 0.5 * dble(count(line)) * microstr%gridspacing(1) * lnl%samplescalefactor
call mem%dealloc(line, 'line')

io_real(1) = t 
call Message%writeValue('  Cylinder radius       :', io_real, 1)

! compute the limiting wave numbers for the outer and inner Ewald spheres
kouter = getXRDwavenumber(sngl(lnl%VoltageH))
kinner = getXRDwavenumber(sngl(lnl%VoltageL))
lambdamin = 1.0/kouter
lambdamax = 1.0/kinner

!=============================================
!=============================================
! crystallography section 
!allocate(cell)        
verbose = .TRUE.

call cell%setFileName(lnl%xtalname)
call Diff%setrlpmethod('XR')

dmin = 0.05
call Diff%setV(dble(1.0)) ! any value will work except 0.0
call Initialize_Cell(cell, Diff, SG, Dyn, EMsoft, dmin, verbose, useHDF=HDF, noLUT=.TRUE.)

if ((SG%getSpaceGroupXtalSystem().eq.5).and.(cell%getLatParm('b').eq.cell%getLatParm('c'))) then
    call Message%printMessage( (/ &
    '                                                                         ', &
    ' ========Program Aborted========                                         ', &
    ' The Laue pattern simulation for rhombohedral/trigonal structures        ', &
    ' requires that the structure be described using the hexagonal reference  ', &
    ' frame.  Please re-enter the crystal structure in this setting and re-run', &
    ' the program.                                                            '/) )
    stop
end if

!=============================================
!=============================================
! compute possible reflection list with kinematical intensities, and intensities
intfactor = 0.0001D0   ! default intensity cutoff factor (from EMLauemaster program)
LaueReflist = LaueReflist_T( grow = .TRUE. )
call LaueReflist%Init_Unit_Reflist(cell, SG, Diff, lambdamin, intfactor, gcnt, verbose, shortg)

! we need to compute the correct value for Lstart... setting to 8 for now
! (this is used for the backprojections, and should probably become a namelist parameter)
Lstart = 8

!=============================================
!=============================================
! start creation of the output file

! Open a new file
  hdfname = trim(EMsoft%generateFilePath('EMdatapathname',lnl%hdfname))

! but delete it first if it already exists
  inquire(file=trim(hdfname), exist=f_exists)
  if (f_exists) then
    open(unit=dataunit, file=trim(hdfname), status='old',form='unformatted')
    close(unit=dataunit, status='delete')
  end if

  hdferr =  HDF%createFile(hdfname)

! write the EMheader to the file
  datagroupname = trim(HDFnames%get_ProgramData())
  call HDF%writeEMheader(EMsoft, dstr, tstrb, tstre, progname, datagroupname)

! add the Duration field to the EMheader group
  hdferr = HDF%openGroup(HDFnames%get_EMheader())
  hdferr = HDF%openGroup(HDFnames%get_ProgramData())

dataset = SC_Duration
  call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
  tstop = 0
  if (g_exists) then
    hdferr = HDF%writeDatasetFloat(dataset, tstop, overwrite)
  else
    hdferr = HDF%writeDatasetFloat(dataset, tstop)
  end if
  call HDF%pop()
  call HDF%pop()

groupname = SC_NMLfiles
  hdferr = HDF%createGroup(HDFnames%get_NMLfiles())

! read the text file and write the array to the file
dataset = trim(HDFnames%get_NMLfilename())
  hdferr = HDF%writeDatasetTextFile(dataset, EMsoft%nmldeffile)

! leave this group
  call HDF%pop()

! create a namelist group to write all the namelist files into
  hdferr = HDF%createGroup(HDFnames%get_NMLparameters())

  call self%writeHDFNameList(HDF, HDFnames)

! leave this group
  call HDF%pop()

! then the remainder of the data in a EMData group
  hdferr = HDF%createGroup(HDFnames%get_EMData())
  hdferr = HDF%createGroup(HDFnames%get_ProgramData())

! create the Lauemaster group and add a HDF_FileVersion attribbute to it 
  HDF_FileVersion = '4.0'
  HDF_FileVersion = cstringify(HDF_FileVersion)
  attributename = SC_HDFFileVersion
  hdferr = HDF%addStringAttributeToGroup(attributename, HDF_FileVersion)

! finally, write all the necessary data:  orientations and simulated patterns along with geometrical parameters
dataset = 'kouter'
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then 
      hdferr = HDF%writeDatasetFloat(dataset, kouter, overwrite)
    else
      hdferr = HDF%writeDatasetFloat(dataset, kouter)
    end if

dataset = 'kinner'
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then 
      hdferr = HDF%writeDatasetFloat(dataset, kinner, overwrite)
    else
      hdferr = HDF%writeDatasetFloat(dataset, kinner)
    end if

dataset = 'numangles'
    call H5Lexists_f(HDF%getobjectID(),trim(dataset),g_exists, hdferr)
    if (g_exists) then 
      hdferr = HDF%writeDatasetInteger(dataset, numangles, overwrite)
    else
      hdferr = HDF%writeDatasetInteger(dataset, numangles)
    end if

! leave the output HDF5 file open so that we can write the results


! create the hyperslabs and write zeroes to them for now
call mem%alloc(patternsum, (/ lnl%Ny, lnl%Nz, numangles /), 'patternsum', initval = 0.0)

! dataset = 'LauePatterns'
!     dims3 = (/ lnl%Ny, lnl%Nz, numangles /)
!     cnt3 = (/ lnl%Ny, lnl%Nz, numangles /)
!     offset3 = (/ 0, 0, 0 /)
!     hdferr = HDF_writeHyperslabFloatArray3D(dataset, patternsum, dims3, offset3, cnt3, HDF_head)
! deallocate(patternsum)

! should we add the backprojections to the same file ?
if (trim(lnl%backprojection).eq.'Yes') then 
  BPnpx = 2*lnl%BPx+1
  BPnpy = 2*lnl%BPx+1
  call mem%alloc(bppatterns, (/ BPnpx, BPnpy, numangles /), 'bppaterns', initval = 0.0)

! dataset = 'backprojections'
!     dims3 = (/ BPnpx, BPnpy, numangles /)
!     cnt3 = (/ BPnpx, BPnpy, numangles /)
!     offset3 = (/ 0, 0, 0 /)
!     hdferr = HDF_writeHyperslabFloatArray3D(dataset, bppatterns, dims3, offset3, cnt3, HDF_head)
end if


! projectionmode = 'T'
  call Message%printMessage(' Initializing sampling geometry for Laue Transmission mode')
  
! all computations are done in mm
  ds = lnl%Lx 
  dsvec = (/ ds, 0.D0, 0.D0 /)   ! vector to the slit plane in x-direction 
  d0 = lnl%Sx 
  d0vec = (/ d0, 0.D0, 0.D0 /)   ! vector to the sample front surface in x-direction
  d = lnl%sampletodetector 
  dvec = (/ d0+d, 0.D0, 0.D0 /)  ! vector to the detector in x-direction
  slitc = (/ 0.D0, lnl%Ly, lnl%Lz /)  ! location of slit center in slit plane
! half width and height of the slit 
  sw = lnl%Lw*0.5D0
  sh = lnl%Lh*0.5D0

! vectors to the slit corners
  slitcorners(1:3,1) = dsvec + slitc + (/ 0.D0, sw, sh /)
  slitcorners(1:3,2) = dsvec + slitc + (/ 0.D0,-sw, sh /)
  slitcorners(1:3,3) = dsvec + slitc + (/ 0.D0,-sw,-sh /)
  slitcorners(1:3,4) = dsvec + slitc + (/ 0.D0, sw,-sh /)

! unit vectors to the slit corners
  do i=1,4 
    scuvec(1:3,i) = slitcorners(1:3,i) / vecnorm(slitcorners(1:3,i))
  end do

! slit corner vectors extended to the detector plane; delineates the projected image of the slit
  do i=1, 4
    l = dvec(1) / scuvec(1,i)
    scdet(1:3,i) = l * scuvec(1:3,i)
  end do

! the sample is cylindrical, so we set up the sampling points inside the cylinder on rays from the 
! source to the detector; one ray per detector pixel.  

! the following are the corners of the slit projection on the detector in units of number of pixels
  do i=1,4 
    shadow(1:2,i) = nint( scdet(2:3,i) / lnl%ps )
  end do 

! for each detector pixel, we have a ray and along this ray we have multiple sampling points that lie
! inside the sample cylinder.  Since we do not know the number of points ahead of time we'll use a 
! linked list and convert it to a regular array when we're done.
  nullify(samplinglist)
  allocate(samplinglist)
  stmp => samplinglist 
  nullify(stmp%next)
  npoints = 0
  dt = lnl%samplingstepsize 
  call mem%alloc(nst, (/ 2, shadow(1,1)-1, shadow(2,2)-1 /), 'nst', initval = 0, startdims = (/ 1, shadow(1,2), shadow(2,3) /) )

! we need to keep the sampling positions as well as the ray vectors for those rays that 
! intersect the cylindrical sample; we offset the detector pixels by 1/2 so that the origin falls
! between pixels; we do the same thing inside the cylinder along each ray so that the first 
! sampling point falls 0.5*dt from the cylinder surface.
  pre = (lnl%ps/(lnl%Sx + lnl%sampletodetector))**2
  do iz=shadow(2,3),shadow(2,2)-1
    do iy=shadow(1,2),shadow(1,1)-1
      pixel = (/ scdet(1,1), (dble(iy)+0.5D0)*lnl%ps, (dble(iz)+0.5D0)*lnl%ps /) / scdet(1,1) 
      k = pixel(2)**2
      disc = (1.D0+k)*0.25D0*t*t - k*lnl%Sx 
      if (disc.ge.0.D0) then 
        if (disc.eq.0.D0) then ! this will rarely happen, so we don't need to include these points 
          ! npoints = npoints + 1  
          ! nst(1:2,iy,iz) = (/ npoints, npoints /)
          ! tt = lnl%Sx / (1.D0+k)
          ! stmp%pos = (/ tt, pixel(2)*tt, pixel(3)*tt /)
          ! stmp%ray = pixel
          ! stmp%front = stmp%pos
          ! stmp%back = stmp%pos
          ! stmp%pixely = iy 
          ! stmp%pixelz = iz
          ! stmp%voxelvolume = 0.D0    ! surface voxels are not counted
          ! allocate( stmp%next )
          ! stmp => stmp%next
          ! nullify(stmp%next)
        else 
          tm = (lnl%Sx-sqrt(disc))/(1.D0+k)
          tp = (lnl%Sx+sqrt(disc))/(1.D0+k)
          stmp%front = (/ tm, pixel(2)*tm, pixel(3)*tm /)
          stmp%back = (/ tp, pixel(2)*tp, pixel(3)*tp /)
          nt = int( (tp-tm)/dt )+1 
          dtp = (tp-tm)/nt
          previous = tm
          do it=0,nt 
            npoints = npoints + 1  
            if (it.eq.0) nst(1,iy,iz) = npoints
            if (it.eq.nt) nst(2,iy,iz) = npoints
            newt = tm+(dble(it)+0.5D0)*dtp
            stmp%pos = (/ newt, pixel(2)*newt, pixel(3)*newt /)
            stmp%ray = pixel
            stmp%pixely = iy 
            stmp%pixelz = iz
            stmp%voxelvolume = pre * dtp * (0.75D0*dtp*dtp + newt*newt)
            allocate( stmp%next )
            stmp => stmp%next
            nullify(stmp%next)
          end do 
        end if 
      end if 
    end do 
  end do 

! convert the linked list to regular array 
  call mem%alloc(slist, (/ 3, npoints /), 'slist', initval = 0.D0)
  call mem%alloc(rays, (/ 3, npoints /), 'rays', initval = 0.D0)
  call mem%alloc(voxvol, (/ npoints /), 'voxvol', initval = 0.D0)
  call mem%alloc(rotated_slist, (/ 3, npoints /), 'rotated_slist', initval = 0.D0)
  call mem%alloc(spppix, (/ shadow(1,1)-1, shadow(2,2)-1 /), 'spppix', initval = 0, startdims = (/ shadow(1,2), shadow(2,3) /) )
  call mem%alloc(frontpts, (/ 3, shadow(1,1)-1, shadow(2,2)-1 /), 'frontpts', initval = 0.D0, &
                 startdims = (/ 1, shadow(1,2), shadow(2,3) /) )
  call mem%alloc(backpts, (/ 3, shadow(1,1)-1, shadow(2,2)-1 /), 'backpts', initval = 0.D0, &
                 startdims = (/ 1, shadow(1,2), shadow(2,3) /) )

  stmp => samplinglist
  do it=1,npoints 
    spppix( stmp%pixely, stmp%pixelz ) = spppix( stmp%pixely, stmp%pixelz ) + 1 
    if (maxval(frontpts(1:3,stmp%pixely, stmp%pixelz)).gt.0.D0) then 
      frontpts(1:3,stmp%pixely, stmp%pixelz) = stmp%front(1:3)
      backpts(1:3,stmp%pixely, stmp%pixelz) = stmp%back(1:3)
    end if
    slist(1:3,it) = stmp%pos
    rays(1:3,it) = stmp%ray
    voxvol(it) = stmp%voxelvolume
    stmp => stmp%next
  end do 

! we'll transfer the sampling points to the sample reference frame by a translation along x
! this way we can perform the sample rotations and do the interpolations needed for proper sampling
! The other data points remain in the source reference frame so that we can compute the incident wave vector
  slist(1,:) = slist(1,:) - lnl%sampletodetector

! delete the linked list completely 
  stmp => samplinglist 
  stmpb => stmp%next 
  do  
    if (associated(stmp)) deallocate(stmp)
    if (.not.associated(stmpb)) EXIT 
    stmp => stmpb 
    stmpb => stmp%next 
  end do 
  nullify(samplinglist, stmp, stmpb)

! next we need the Legendre lattitudes for the back projector 
if (trim(lnl%backprojection).eq.'Yes') then 
  call Message%printMessage(' Computing Legendre lattitudinal grid values')
  call mem%alloc(diagonal, (/ BPnpx /), 'diagonal', initval=0.D0)
  call mem%alloc(upd, (/ BPnpx /), 'upd')
  upd = (/ (dble(i) / dsqrt(4.D0 * dble(i)**2 - 1.D0), i=1,BPnpx) /)
  call dsterf(BPnpx-2, diagonal, upd, info) 
  ! the eigenvalues are stored from smallest to largest and we need them in the opposite direction
  call mem%alloc(LegendreArray, (/ BPnpx-1 /), 'LegendreArray', startdims = (/ 0 /) )
  LegendreArray(0:BPnpx-1) = diagonal(BPnpx:1:-1)
  ! set the center eigenvalue to 0
  LegendreArray((BPnpx-1)/2) = 0.D0
  call mem%dealloc(diagonal, 'diagonal') 
  call mem%dealloc(upd, 'upd')
end if 

! generate the sample rotation QuaternionArray_T class 
qAR = QuaternionArray_T( n = numangles, s = 'd' )
do it=1,numangles
  th = ( dble(it-1)*2.D0*cPi/dble(numangles) ) / 2.D0 ! need the semi-angle
  quat = Quaternion_T( qd = (/ cos(th), 0.D0, 0.D0, sin(th) /) )
  call qAR%insertQuatinArray(it, quat)
end do 

! define the ipar, fpar, and dpar arrays to pass parameters to the actual computation
! we'll give them 20 entries each to make sure that we allow for future additions
ipar(1) = lnl%Ny 
ipar(2) = lnl%Nz 
ipar(3) = gcnt 

dpar(1) = lnl%Dy 
dpar(2) = lnl%Dz 
dpar(3) = lnl%ps 
dpar(4) = lnl%samplingstepsize 
dpar(5) = lnl%absl  
dpar(6) = lnl%sampletodetector 
dpar(7) = lnl%spotw  
dpar(8) = t 
dpar(9) = lambdamin 
dpar(10)= lambdamax
! dpar(9) = lnl% 



!=============================================
!=============================================
! Here we perform the actual simulations; the threads cover all the sample voxels for a 
! given sample orientation and the resulting patterns are summed together.
! Then we go to the next orientation.
! Write them to the HDF5 output file, along with (optionally) the pattern tiff files

memth = Memory_T( nt=lnl%nthreads )

! set the number of OpenMP threads 
  call OMP_SET_NUM_THREADS(lnl%nthreads)
  io_int(1) = lnl%nthreads
  call Message%WriteValue(' Attempting to set number of threads to ',io_int, 1, frm = "(I4)")

! outer loop ... 
  do ii = 1, numangles
! rotate all the sampling points except for the first angle which is zero
    if (ii.eq.1) then 
      rotated_slist = slist 
    else 
      quat = qAR%getQuatfromArray( ii ) 
      rotated_slist = quat%quat_Lp_vecarray( npoints, slist )
    end if 
! transform the list to the microstructure reference frame
    rotated_slist(1,:) =  rotated_slist(1,:) / microstr%gridspacing(1) + microstr%dimensions(1)/2
    rotated_slist(2,:) =  rotated_slist(2,:) / microstr%gridspacing(2) + microstr%dimensions(2)/2
    rotated_slist(3,:) =  rotated_slist(3,:) / microstr%gridspacing(3) + microstr%dimensions(3)/2

! use OpenMP to run on multiple cores ... 
! We loop over all the detector pixels that are illuminated by the slit. For each such pixel,
! we pass on the list of rotated sampling points along with their number; the routine then returns
! the pattern for that pixel, including all the diffracted beams along the ray. Each thread works on
! a horizontal row of pixels before returning the accumulated pattern and adding it to the overall pattern.

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(TID, iy, pattern, pp, np, orlist)

      NUMTHREADS = OMP_GET_NUM_THREADS()
      TID = OMP_GET_THREAD_NUM()
      call memth%alloc(pattern, (/ lnl%Ny, lnl%Nz /), 'pattern', initval=0.0, TID = TID)
      call memth%alloc(pp, (/ lnl%Ny, lnl%Nz /), 'pp', initval=0.0, TID = TID)

!$OMP DO SCHEDULE(DYNAMIC)
      do iz = shadow(2,3),shadow(2,2)-1
        do iy= shadow(1,2),shadow(1,1)-1
          np = nst(2, iy, iz) - nst(1, iy, iz) + 1 
! for this ray, get the sequence of orientations from the microstructure 
          call getRayOrientations( microstr, rotated_slist(1:3,nst(1,iy,iz):nst(2,iy,iz)), np, orlist )
! pass only the relevant portions for this pixel of the slist and rotated_slist arrays etc ... 
          pp = LaueReflist%getnewLaueDCTPattern(ipar, fpar, dpar, np, slist(1:3,nst(1,iy,iz):nst(2,iy,iz)), &
                                                rays(1:3,nst(1,iy,iz)), voxvol(nst(1,iy,iz):nst(2,iy,iz)), & 
                                                quat, orlist )
          call orlist%deleteArray() ! explicitly clean up the orientations array
          pattern = pattern + pp 
        end do 

!$OMP CRITICAL
          patternsum(:, :, ii) = patternsum(:, :, ii) + pattern(:, :) 
!$OMP END CRITICAL

     end do 
!$OMP END DO

      if (TID.eq.0) then 
        io_int(1) = ii 
        call Message%WriteValue(' patterns completed for batch ', io_int, 1) 
      end if 
      call memth%dealloc(pattern, 'pattern', TID = TID)
      call memth%dealloc(pp, 'pp', TID = TID)

! end of OpenMP portion
!$OMP END PARALLEL

    if (lnl%binarize.eqv..TRUE.) then 
      mps = maxval(patternsum(:,:,ii))
      where (patternsum(:,:,ii) .gt. mps*0.01)
        patternsum(:,:,ii) = 1.0
      end where
      where (patternsum(:,:,ii) .le. mps*0.01)
        patternsum(:,:,ii) = 0.0
      end where
    end if 
end do ! outer loop
! 

  call timer%Time_tock(1)
  tstop = timer%getInterval(1)
  io_int(1) = tstop
  call Message%WriteValue('Execution time [s]: ',io_int,1)

! apply the rectangular beam stop to all patterns; we assume that the beam stop is 
! centered.  First convert the beam stop dimensions to pixels, then multiply all pixels 
! inside the beam stop by the attenuation factor
  bsw = int( 0.5 * lnl%beamstopwidth / lnl%ps )
  bsh = int( 0.5 * lnl%beamstopheight / lnl%ps )
  NNy = lnl%Ny/2
  NNz = lnl%Nz/2
  io_int(1:2) = (/ 2*bsw, 2*bsh /)
  call Message%WriteValue(' ---> Applying beam stop attenuation for stop size ', io_int, 2)

  patternsum(NNy-bsw:NNy+bsw,NNz-bsh:NNz+bsh,1:numangles) = & 
      patternsum(NNy-bsw:NNy+bsw,NNz-bsh:NNz+bsh,1:numangles) * lnl%beamstopatf


! write the hyperslab to the HDF5 file 
dataset = 'LauePatterns'
  hdferr = HDF%writeDatasetFloatArray(dataset, patternsum, lnl%Ny, lnl%Nz, numangles)

  if (trim(lnl%backprojection).eq.'Yes') then 
    call Message%printMessage('Starting pattern back projection computation')
    call mem%alloc(bp, (/ BPnpx, BPnpy /), 'bp')
    BPmode = 'forward'
    do ii=1,numangles
      bp = backprojectLauePattern( (/kouter, kinner/), sngl(lnl%ps), sngl(lnl%sampletodetector), Lstart, (/lnl%Ny, lnl%Nz/), &
                                   (/lnl%BPx, lnl%BPx/), patternsum(1:lnl%Ny,1:lnl%Nz,ii), BPmode, LegendreArray)
      ! write (*,*) ' max intensity in backprojection : ', maxval(bp), maxval(patternsum(1:lnl%Ny,1:lnl%Nz,ii))
      bppatterns(:,:,ii) = bp
    end do
    call mem%dealloc(bp, 'bp')

dataset = 'backprojections'
    hdferr = HDF%writeDatasetFloatArray(dataset, bppatterns, BPnpx, BPnpy, numangles)
  end if

 call HDF%pop()
 call HDF%pop()
 tstre = timer%getTimeString()

! update the time string
groupname = SC_EMheader
  hdferr = HDF%openGroup(groupname)
  hdferr = HDF%openGroup(datagroupname)

dataset = SC_StopTime
  line2(1) = dstr//', '//tstre
  hdferr = HDF%writeDatasetStringArray(dataset, line2, 1, overwrite)

dataset = SC_Duration
  hdferr = HDF%writeDatasetFloat(dataset, tstop, overwrite)

 call HDF%pop(.TRUE.)

! and close the fortran hdf interface
 call closeFortranHDFInterface()

 call Message%printMessage("LaueDCT: patterns stored in "//trim(HDFname))

if (trim(lnl%tiffprefix).ne.'undefined') then 
  npx = lnl%Ny 
  npy = lnl%Nz 

! use the same intensity scaling factors for all patterns 
  ma = maxval(patternsum)
  mi = minval(patternsum)

! optionally, write the individual tiff image files 
  do ii=1,numangles
    write (pnum,"(I4.4)") ii
    ! fname = trim(EMsoft_getEMdatapathname())//trim(lnl%tiffprefix)//'_'//pnum//'.tiff'
    fname = EMsoft%generateFilePath('EMdatapathname',trim(lnl%tiffprefix)//'_'//pnum//'.tiff')
    TIFF_filename = trim(fname)

! allocate memory for image
    allocate(TIFF_image(npx,npy))

! fill the image with whatever data you have (between 0 and 255)
    TIFF_image = 255 - int(255 * (patternsum(:,:,ii)-mi)/(ma-mi))

! set up the image_t structure
    im = image_t(TIFF_image)
    if(im%empty()) call Message%printMessage("LaueDCT","failed to convert array to image")

! create the file
    call im%write(trim(TIFF_filename), iostat, iomsg) ! format automatically detected from extension
    if(0.ne.iostat) then
        call Message%printMessage("failed to write image to file : "//iomsg)
    end if 
    deallocate(TIFF_image)
  end do 
end if 

end associate 

end subroutine LaueDCT_



end module mod_LaueDCT
