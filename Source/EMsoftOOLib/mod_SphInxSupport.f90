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

module mod_SphInxSupport
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/22/20
  !!
  !! class definition for the EMSphInxSupport program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMSphInxSupport program
! namelist for the EMSphInx program
type, public :: SphInxNameListType
! Indexing Parameters
    integer(kind=irg)       :: bw             ! spherical harmonic bandwidth
    logical                 :: normed         ! normalization flag
    logical                 :: refine         ! refinement flag
    logical                 :: flipy          ! pattern vertical flip flag
    integer(kind=irg)       :: ROImask(4)     ! to be supported in C++ version
    character(fnlen)        :: ROIfile        ! mask file [for f90 version only]
    integer(kind=irg)       :: nregions       ! for AHE
    integer(kind=irg)       :: nthread        ! parallel threads to be used [if 0, then omp_get_max_threads()]
    integer(kind=irg)       :: batchsize      ! number of patterns to dispatch to a thread at once [indexer.hpp, BatchEstimate]

! Scan Parameters and Camera Calibration
    real(kind=sgl)          :: scandims(4)    ! [ipf_wd, ipf_ht, stepx, stepy]   
    integer(kind=sgl)       :: patdims(2)     ! [numsx, numsy] pattern dimensions
    real(kind=sgl)          :: delta          ! detector pixel size [microns]
    real(kind=sgl)          :: pctr(3)        ! [xpc, ypc, L] pattern center coordinates in any vendor convention
    character(fnlen)        :: vendor         ! ['EMsoft', 'EDAX/TSL', 'Oxford', 'Bruker']
    real(kind=sgl)          :: thetac         ! detector tilt angle from vertical
    integer(kind=irg)       :: binning        ! pattern binning factor
    logical                 :: circmask       ! replaces 'maskpattern' and 'maskradius' (radius == smallest side/2)
!   logical                 :: flipy          ! should patterns be flipped upside down ? [C++ only ?]

! Input Data
    character(fnlen)        :: masterfile     ! master pattern file
    character(fnlen)        :: patfile        ! experimental file
    character(fnlen)        :: HDFstrings(10) ! this is 'patdset' in C++ version
    character(fnlen)        :: inputtype      ! this is determined by extension in the C++ version

! Output Data
    character(fnlen)        :: datafile       ! hdf5, must be defined
    character(fnlen)        :: ctffile        ! not currently supported in C++ version
    character(fnlen)        :: angfile        ! not currently supported in C++ version
end type SphInxNameListType

end module mod_SphInxSupport