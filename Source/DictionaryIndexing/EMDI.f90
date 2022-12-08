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

program EMDI
  !! author: MDG
  !! version: 1.0 
  !! date: 04/02/20
  !!
  !! Indexing of EBSD/ECP/TKD patterns using a dynamically generated pattern dictionary
  !!
  !! This program has a different structure from most other programs because the indexing
  !! routine must be callable from C++ (EMsoftWorkBench).  So the regular class structure
  !! is handled a bit differently, and the name list is not initialized in the main program
  !! at all. The DIdriver routine is also not part of a class, but is called as a regular 
  !! routine. 

use mod_kinds
use mod_global
use mod_EMsoft
use mod_DI
use mod_HDFsupport
use ISO_C_BINDING
use mod_FitOrientation

IMPLICIT NONE

character(fnlen)        :: progname
character(fnlen)        :: progdesc = 'Indexing of EBSD/ECP/TKD patterns using a dynamically calculated pattern dictionary'

type(EMsoft_T)          :: EMsoft
character(kind=c_char)  :: Cprogname(fnlen)
character(kind=c_char)  :: Cnmldeffile(fnlen)
type(DIfile_T)          :: DIFT
type(FitOrientation_T)  :: FitOr

progname = 'EMDI.f90'

! print the EMsoft header and handle any command line arguments  
EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 80 /) )

! call the DIdriver routine to take care of the entire indexing process 
Cnmldeffile = carstringify(EMsoft%nmldeffile)
Cprogname = carstringify(progname)
call DIdriver(Cnmldeffile, Cprogname, C_NULL_FUNPTR, C_NULL_FUNPTR, C_NULL_FUNPTR, 0_c_size_t)

! next, get the nml file to see if we also need to start up the refinement
DIFT = DIfile_T(EMsoft%nmldeffile)

if (DIFT%getrefinementfilename().ne.'undefined') then 
    progname = 'EMFitOrientation.f90'
    progdesc = 'Refine orientations by searching orientation space about a point including the pseudosymmetric variant(s)'
    call EMsoft%printEMsoftHeader(progname, progdesc)

! deal with the namelist stuff
    call EMsoft%setnmldeffile(trim(DIFT%getrefinementfilename()))
    FitOr = FitOrientation_T(EMsoft%nmldeffile)

! perform the computations
    call FitOr%FitOrientation(EMsoft, progname)
end if 

end program EMDI
