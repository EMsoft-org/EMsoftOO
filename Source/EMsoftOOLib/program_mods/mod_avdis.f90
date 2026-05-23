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

module mod_avdis
  !! author: MDG 
  !! version: 1.0 
  !! date: 11/15/25
  !!
  !! class definition for the EMavdis program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMavdis program
type, public :: avdisNameListType
  integer(kind=irg) :: cnumstart
  integer(kind=irg) :: cstepsize
  integer(kind=irg) :: cnsteps
  character(fnlen)  :: avcsvfile
  character(2)      :: sampletype
  character(4)      :: sampleoffset
end type avdisNameListType

! class definition
type, public :: avdis_T
private 
  character(fnlen)          :: nmldeffile = 'EMavdis.nml'
  type(avdisNameListType)   :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: setcnumstart_
  procedure, pass(self) :: getcnumstart_
  procedure, pass(self) :: setcstepsize_
  procedure, pass(self) :: getcstepsize_
  procedure, pass(self) :: setcnsteps_
  procedure, pass(self) :: getcnsteps_
  procedure, pass(self) :: setavcsvfile_
  procedure, pass(self) :: getavcsvfile_
  procedure, pass(self) :: setsampletype_
  procedure, pass(self) :: getsampletype_
  procedure, pass(self) :: setsampleoffset_
  procedure, pass(self) :: getsampleoffset_
  procedure, pass(self) :: avdis_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: setcnumstart => setcnumstart_
  generic, public :: getcnumstart => getcnumstart_
  generic, public :: setcstepsize => setcstepsize_
  generic, public :: getcstepsize => getcstepsize_
  generic, public :: setcnsteps => setcnsteps_
  generic, public :: getcnsteps => getcnsteps_
  generic, public :: setavcsvfile => setavcsvfile_
  generic, public :: getavcsvfile => getavcsvfile_
  generic, public :: setsampletype => setsampletype_
  generic, public :: getsampletype => getsampletype_
  generic, public :: setsampleoffset => setsampleoffset_
  generic, public :: getsampleoffset => getsampleoffset_
  generic, public :: avdis => avdis_

end type avdis_T

! the constructor routine for this class 
interface avdis_T
  module procedure avdis_constructor
end interface avdis_T

contains

!--------------------------------------------------------------------------
type(avdis_T) function avdis_constructor( nmlfile ) result(avdis)
!! author: MDG 
!! version: 1.0 
!! date: 11/15/25
!!
!! constructor for the avdis_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call avdis%readNameList(nmlfile)

end function avdis_constructor

!--------------------------------------------------------------------------
subroutine avdis_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 11/15/25
!!
!! destructor for the avdis_T Class
 
IMPLICIT NONE

type(avdis_T), INTENT(INOUT)  :: self 

call reportDestructor('avdis_T')

end subroutine avdis_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 11/15/25
!!
!! read the namelist from an nml file for the avdis_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(avdis_T), INTENT(INOUT)        :: self
character(fnlen),INTENT(IN)          :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)          :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                       :: EMsoft 
type(IO_T)                           :: Message       
logical                              :: skipread = .FALSE.

integer(kind=irg)                    :: cnumstart
integer(kind=irg)                    :: cstepsize
integer(kind=irg)                    :: cnsteps
character(fnlen)                     :: avcsvfile
character(2)                         :: sampletype
character(4)                         :: sampleoffset

namelist  / avdis / avcsvfile, cnumstart, cstepsize, cnsteps, sampletype, sampleoffset

avcsvfile = 'undefined'
cnumstart = 20
cstepsize = 20
cnsteps = 16
sampletype = 'cP'
sampleoffset =  'none'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
    open(UNIT=dataunit,FILE=trim(nmlfile),DELIM='apostrophe',STATUS='old')
    read(UNIT=dataunit,NML=avdis)
    close(UNIT=dataunit,STATUS='keep')

! check for required entries
    if (trim(avcsvfile).eq.'undefined') then
        call Message%printError('readNameList:',' avcsvfile file name is undefined in '//nmlfile)
    end if
end if

self%nml%avcsvfile = avcsvfile
self%nml%cnumstart = cnumstart
self%nml%cstepsize = cstepsize
self%nml%cnsteps = cnsteps
self%nml%sampletype = sampletype
self%nml%sampleoffset =sampleoffset

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 11/15/25
!!
!! pass the namelist for the avdis_T Class to the calling program

IMPLICIT NONE 

class(avdis_T), INTENT(INOUT)          :: self
type(avdisNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine setcnumstart_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setcnumstart_
!! author: MDG
!! version: 1.0
!! date: 11/15/25
!!
!! set cnumstart in the avdis_T class

IMPLICIT NONE

class(avdis_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN) :: inp

self%nml%cnumstart = inp

end subroutine setcnumstart_

!--------------------------------------------------------------------------
function getcnumstart_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getcnumstart_
!! author: MDG
!! version: 1.0
!! date: 11/15/25
!!
!! get cnumstart from the avdis_T class

IMPLICIT NONE

class(avdis_T), INTENT(INOUT) :: self
integer(kind=irg)             :: out

out = self%nml%cnumstart

end function getcnumstart_

!--------------------------------------------------------------------------
subroutine setcstepsize_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setcstepsize_
!! author: MDG
!! version: 1.0
!! date: 11/15/25
!!
!! set cstepsize in the avdis_T class

IMPLICIT NONE

class(avdis_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN) :: inp

self%nml%cstepsize = inp

end subroutine setcstepsize_

!--------------------------------------------------------------------------
function getcstepsize_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getcstepsize_
!! author: MDG
!! version: 1.0
!! date: 11/15/25
!!
!! get cstepsize from the avdis_T class

IMPLICIT NONE

class(avdis_T), INTENT(INOUT) :: self
integer(kind=irg)             :: out

out = self%nml%cstepsize

end function getcstepsize_

!--------------------------------------------------------------------------
subroutine setcnsteps_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setcnsteps_
!! author: MDG
!! version: 1.0
!! date: 11/15/25
!!
!! set cnsteps in the avdis_T class

IMPLICIT NONE

class(avdis_T), INTENT(INOUT) :: self
integer(kind=irg), INTENT(IN) :: inp

self%nml%cnsteps = inp

end subroutine setcnsteps_

!--------------------------------------------------------------------------
function getcnsteps_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getcnsteps_
!! author: MDG
!! version: 1.0
!! date: 11/15/25
!!
!! get cnsteps from the avdis_T class

IMPLICIT NONE

class(avdis_T), INTENT(INOUT) :: self
integer(kind=irg)             :: out

out = self%nml%cnsteps

end function getcnsteps_

!--------------------------------------------------------------------------
subroutine setavcsvfile_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setavcsvfile_
!! author: MDG
!! version: 1.0
!! date: 11/15/25
!!
!! set avcsvfile in the avdis_T class

IMPLICIT NONE

class(avdis_T), INTENT(INOUT) :: self
character(fnlen), INTENT(IN)  :: inp

self%nml%avcsvfile = trim(inp)

end subroutine setavcsvfile_

!--------------------------------------------------------------------------
function getavcsvfile_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getavcsvfile_
!! author: MDG
!! version: 1.0
!! date: 11/15/25
!!
!! get avcsvfile from the avdis_T class

IMPLICIT NONE

class(avdis_T), INTENT(INOUT) :: self
character(fnlen)              :: out

out = trim(self%nml%avcsvfile)

end function getavcsvfile_

!--------------------------------------------------------------------------
subroutine setsampletype_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setsampletype_
!! author: MDG
!! version: 1.0
!! date: 11/15/25
!!
!! set sampletype in the avdis_T class

IMPLICIT NONE

class(avdis_T), INTENT(INOUT)     :: self
character(2), INTENT(IN)          :: inp

self%nml%sampletype = trim(inp)

end subroutine setsampletype_

!--------------------------------------------------------------------------
function getsampletype_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getsampletype_
!! author: MDG
!! version: 1.0
!! date: 11/15/25
!!
!! get sampletype from the avdis_T class

IMPLICIT NONE

class(avdis_T), INTENT(INOUT)     :: self
character(2)                      :: out

out = trim(self%nml%sampletype)

end function getsampletype_

!--------------------------------------------------------------------------
subroutine setsampleoffset_(self,inp)
!DEC$ ATTRIBUTES DLLEXPORT :: setsampleoffset_
!! author: MDG
!! version: 1.0
!! date: 11/15/25
!!
!! set sampleoffset in the avdis_T class

IMPLICIT NONE

class(avdis_T), INTENT(INOUT)     :: self
character(4), INTENT(IN)          :: inp

self%nml%sampleoffset = trim(inp)

end subroutine setsampleoffset_

!--------------------------------------------------------------------------
function getsampleoffset_(self) result(out)
!DEC$ ATTRIBUTES DLLEXPORT :: getsampleoffset_
!! author: MDG
!! version: 1.0
!! date: 11/15/25
!!
!! get sampleoffset from the avdis_T class

IMPLICIT NONE

class(avdis_T), INTENT(INOUT)     :: self
character(4)                      :: out

out = trim(self%nml%sampleoffset)

end function getsampleoffset_

!--------------------------------------------------------------------------
subroutine avdis_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: avdis_
!! author: MDG 
!! version: 1.0 
!! date: 11/15/25
!!
!! perform the computations

use mod_EMsoft
use mod_so3
use mod_rotations
use mod_quaternions
use mod_io

IMPLICIT NONE 

class(avdis_T), INTENT(INOUT)           :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(so3_T)                             :: SO, SO2 
type(c_T)                               :: cu
type(c_T), allocatable                  :: NNlist(:)
type(q_T)                               :: qu, qu2
type(r_T)                               :: ro
type(Quaternion_T)                      :: q1, q2, Mu
type(IO_T)                              :: Message 

integer(kind=irg)                       :: idelta, ic, jc, kc, igroup, n, NN, nump, numpk, pgnum(11), nvalid, ival
integer(kind=irg)                       :: histogram(1000,11), io_int(1), ntotal
real(kind=dbl)                          :: delta, coor(3), coor2(3), Mu1(4), dhist, mval
real(kind=dbl),allocatable              :: disor(:,:,:)
character(fnlen)                        :: fname

associate( nml=>self%nml )
call setRotationPrecision('d')

pgnum = (/ 1, 3, 6, 9, 12, 16, 18, 21, 24, 28, 30 /)

n = 6 
if (nml%sampletype.eq.'cI') n = 14
if (nml%sampletype.eq.'cF') n = 12
allocate( NNlist(n) )
dhist = 0.01D0
histogram = 0 

fname = EMsoft%generateFilePath('EMdatapathname',trim(nml%avcsvfile))
open(unit=dataunit, file=trim(fname),status='unknown',form='unformatted')
call Message%printMessage(' created output file : '//trim(fname))

! outer loop is over all the delta values
do idelta=1, nml%cnsteps 
  ntotal = 0
  nump = nml%cnumstart + (idelta-1) * nml%cstepsize
  delta = 0.5D0 * LPs%ap / dble( nump )
  if (nml%sampletype.ne.'cP') then  ! allow for centered lattices
    numpk = nump*2
  else
    numpk = nump 
  end if 
  allocate ( disor( -nump:nump, -nump:nump, -numpk:numpk ) )
  disor = 0.D0 
  SO2 = so3_T( pgnum(1), initshifts=nml%sampletype )
   
  write (*,*) ' idelta: ', idelta, shape(disor), nump, numpk

  ! next loop over all the cubochoric points and determine which RFZ they are part of
  do ic = -nump,nump
    do jc = -nump,nump
      do kc = -numpk,numpk ! this allows for handling of cI and cF lattices
        if (nml%sampletype.eq.'cP') then 
          coor = (/ dble(ic), dble(jc), dble(kc) /) * delta
        end if 
        if ((nml%sampletype.eq.'cF').and.(mod(kc,2).eq.0)) then 
          coor = (/ dble(ic), dble(jc), dble(kc)/2.D0 /) * delta
        end if 
        if (nml%sampletype.eq.'cI') then 
          if (mod(kc,2).eq.0) then 
            coor = (/ dble(ic), dble(jc), dble(kc)/2.D0 /) * delta
          else
            coor = (/ dble(ic)+0.5D0, dble(jc)+0.5D0, dble(kc)/2.D0 /) * delta
          end if 
        end if 
        cu = c_T( cdinp=coor )
        qu = cu%cq()
        Mu1 = qu%q_copyd()
        Mu1(2:4) = -Mu1(2:4)
        q1 = Quaternion_T( qd = Mu1 )
        call SO2%getcuboNN( cu, delta, nml%sampletype, NNlist, n, nvalid )
        do NN=1,nvalid
          qu2 = NNlist(NN)%cq()
          Mu1 = qu2%q_copyd()
          q2 = Quaternion_T( qd = Mu1 )
          Mu = q2 * q1
          Mu1 = Mu%get_quatd()
          disor(ic, jc, kc) = disor(ic, jc, kc) + 2.D0*acos(Mu1(1))
        end do
        disor(ic, jc, kc) = disor(ic, jc, kc) / dble(nvalid)
        if (nml%sampletype.eq.'cF') then 
          if (mod(kc,2).eq.0) then ! add a second site to this layer
            coor = (/ dble(ic)+0.5D0, dble(jc)+0.5D0, dble(kc)/2.D0 /) * delta
            cu = c_T( cdinp=coor )
            qu = cu%cq()
            Mu1 = qu%q_copyd()
            Mu1(2:4) = -Mu1(2:4)
            q1 = Quaternion_T( qd = Mu1 )
            call SO2%getcuboNN( cu, delta, nml%sampletype, NNlist, n, nvalid )
            do NN=1,nvalid
              qu2 = NNlist(NN)%cq()
              Mu1 = qu2%q_copyd()
              q2 = Quaternion_T( qd = Mu1 )
              Mu = q2 * q1
              Mu1 = Mu%get_quatd()
              disor(ic, jc, kc) = disor(ic, jc, kc) + 2.D0*acos(Mu1(1))
            end do
            disor(ic, jc, kc) = disor(ic, jc, kc) / dble(nvalid)
          else
            coor = (/ dble(ic)+0.5D0, dble(jc), dble(kc-1)/2.D0+0.5D0 /) * delta  ! first site
            cu = c_T( cdinp=coor )
            qu = cu%cq()
            Mu1 = qu%q_copyd()
            Mu1(2:4) = -Mu1(2:4)
            q1 = Quaternion_T( qd = Mu1 )
            call SO2%getcuboNN( cu, delta, nml%sampletype, NNlist, n, nvalid )
            do NN=1,nvalid
              qu2 = NNlist(NN)%cq()
              Mu1 = qu2%q_copyd()
              q2 = Quaternion_T( qd = Mu1 )
              Mu = q2 * q1
              Mu1 = Mu%get_quatd()
              disor(ic, jc, kc) = disor(ic, jc, kc) + 2.D0*acos(Mu1(1))
            end do
            disor(ic, jc, kc) = disor(ic, jc, kc) / dble(nvalid)
!
            coor = (/ dble(ic), dble(jc)+0.5D0, dble(kc-1)/2.D0+0.5D0 /) * delta  ! second site
            cu = c_T( cdinp=coor )
            qu = cu%cq()
            Mu1 = qu%q_copyd()
            Mu1(2:4) = -Mu1(2:4)
            q1 = Quaternion_T( qd = Mu1 )
            call SO2%getcuboNN( cu, delta, nml%sampletype, NNlist, n, nvalid )
            do NN=1,nvalid
              qu2 = NNlist(NN)%cq()
              Mu1 = qu2%q_copyd()
              q2 = Quaternion_T( qd = Mu1 )
              Mu = q2 * q1
              Mu1 = Mu%get_quatd()
              disor(ic, jc, kc) = disor(ic, jc, kc) + 2.D0*acos(Mu1(1))
            end do
            disor(ic, jc, kc) = disor(ic, jc, kc) / dble(nvalid)
          end if
        
        end if 
      end do 
    end do 
  end do
  disor = 100.D0 * disor / dtor  
  io_int(1) = (2*nump+1)**2*(2*numpk+1)
  call Message%WriteValue(' total number of points generated : ', io_int, 1)
  
! convert this to a histogram but only if the point is inside the RFZ 
! loop uses the current set of disorientations and goes over the rotational point groups
  do igroup=1,11
    SO = so3_T( pgnum(igroup) )

    write (*,*) ' checking point group ', pgnum(igroup)

    do ic = -nump,nump
      do jc = -nump,nump
        do kc = -nump,nump
          coor = (/ dble(ic), dble(jc), dble(kc) /) * delta
          cu = c_T( cdinp=coor )
          ro = cu%cr()
          if (SO%IsinsideFZ( ro ).eqv..TRUE.) then
            mval = disor(ic, jc, kc)
            if ((mval.lt.1000.D0).and.(mval.gt.0.D0)) then 
              ival = int(mval)
              histogram(ival,igroup) = histogram(ival,igroup) + 1
            end if  
          end if 
        end do
      end do 
    end do 
  end do
  deallocate( disor )

! write this histogram array to a binary datafile; to be replaced with HDF5 file

  write (dataunit) histogram
  histogram = 0

  io_int(1) = idelta 
  call Message%WriteValue(' completed computations for ncubochoric = ',io_int,1) 

end do 

close(unit=dataunit, status='keep')
call Message%printMessage(' -> computation complete.')

end associate

end subroutine avdis_



end module mod_avdis