! ###################################################################
! Copyright (c) 2016-2023, Marc De Graef Research Group/Carnegie Mellon University
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

program EMTKDmaster
  !! author: MDG
  !! version: 1.0 
  !! date: 03/18/20
  !!
  !! TKD Energy-dependent Master Pattern Simulation 

use mod_kinds
use mod_global
use mod_EMsoft
use mod_TKDmaster
use mod_HDFnames
use stringconstants
use mod_EBSDmaster

IMPLICIT NONE

character(fnlen)      :: progname = 'EMTKDmaster.f90'
character(fnlen)      :: progdesc = 'TKD Energy-dependent Master Pattern Simulation'

type(EMsoft_T)        :: EMsoft
type(TKDmaster_T)     :: TKD 
type(EBSDmaster_T)    :: MP 
type(HDFnames_T)      :: HDFnames
real(kind=sgl)        :: thickness

! print the EMsoft header and handle any command line arguments  
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 29, 0 /) )

! deal with the namelist stuff
TKD = TKDmaster_T(EMsoft%nmldeffile)

! We will use the EBSDmaster code to compute the TKD master pattern to avoid
! duplicating large amounts of code... So we need to copy the TKD name list
! parameters into the EBSD master name list (they are pretty much identical)
call MP%set_npx(TKD%get_npx())
call MP%set_Esel(TKD%get_Esel())
call MP%set_nthreads(TKD%get_nthreads())
call MP%set_dmin(TKD%get_dmin())
call MP%set_copyfromenergyfile(TKD%get_copyfromenergyfile())
call MP%set_h5copypath(TKD%get_h5copypath())
call MP%set_energyfile(TKD%get_energyfile())
call MP%set_BetheParametersFile(TKD%get_BetheParametersFile())
call MP%set_combinesites(TKD%get_combinesites())
call MP%set_restart(TKD%get_restart())
call MP%set_uniform(TKD%get_uniform())
call MP%set_Notify(TKD%get_Notify())
call MP%set_kinematical(TKD%get_kinematical())
thickness = TKD%get_thickness()

! set the HDFnames class to TKD mode  
HDFnames = HDFnames_T() 
call HDFnames%set_ProgramData(SC_TKDmaster) 
call HDFnames%set_NMLlist(SC_TKDmasterNameList) 
call HDFnames%set_NMLfilename(SC_TKDmasterNML) 
call HDFnames%set_Variable(SC_MCOpenCL) 

! perform the computations (we use the standard EBSD routine to do this, but with 
! altered HDFnames to make sure we generate the correct HDF5 output file)
call MP%EBSDmaster(EMsoft, progname, HDFnames, thickness)

end program EMTKDmaster
