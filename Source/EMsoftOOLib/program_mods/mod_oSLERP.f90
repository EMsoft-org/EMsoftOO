! ###################################################################
! Copyright (c) 2013-2025, Marc De Graef Research Group/Carnegie Mellon University
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

module mod_oSLERP
  !! author: MDG 
  !! version: 1.0 
  !! date: 07/18/25
  !!
  !! class definition for the EMoSLERP program

use mod_kinds
use mod_global

IMPLICIT NONE 

! namelist for the EMoSLERP program
type, public :: oSLERPNameListType
  integer(kind=irg)       :: framesize
  integer(kind=irg)       :: nthreads
  real(kind=dbl)          :: qm(4)
  real(kind=dbl)          :: mA(3)
  real(kind=dbl)          :: mC(3)
  real(kind=dbl)          :: o1(8)
  real(kind=dbl)          :: o2(8)
  real(kind=dbl)          :: dOmega
  logical                 :: silentrender
  character(fnlen)        :: rendermode
  character(fnlen)        :: GBmode
  character(fnlen)        :: xtalname
  character(fnlen)        :: povrayfile
  character(fnlen)        :: framefolder
  character(fnlen)        :: moviename
  character(fnlen)        :: metric
  character(fnlen)        :: PVincludepath
  character(fnlen)        :: PVexec 
  character(fnlen)        :: PVffmpeg
end type oSLERPNameListType

! class definition
type, public :: oSLERP_T
private 
  character(fnlen)          :: nmldeffile = 'EMoSLERP.nml'
  type(oSLERPNameListType)  :: nml 

contains
private 
  procedure, pass(self) :: readNameList_
  procedure, pass(self) :: getNameList_
  procedure, pass(self) :: oSLERP_

  generic, public :: getNameList => getNameList_
  generic, public :: readNameList => readNameList_
  generic, public :: oSLERP => oSLERP_

end type oSLERP_T

! the constructor routine for this class 
interface oSLERP_T
  module procedure oSLERP_constructor
end interface oSLERP_T

contains

!--------------------------------------------------------------------------
type(oSLERP_T) function oSLERP_constructor( nmlfile ) result(oSLERP)
!! author: MDG 
!! version: 1.0 
!! date: 07/18/25
!!
!! constructor for the oSLERP_T Class; reads the name list 
 
IMPLICIT NONE

character(fnlen), OPTIONAL   :: nmlfile 

call oSLERP%readNameList(nmlfile)

end function oSLERP_constructor

!--------------------------------------------------------------------------
subroutine oSLERP_destructor(self) 
!! author: MDG 
!! version: 1.0 
!! date: 07/18/25
!!
!! destructor for the oSLERP_T Class
 
IMPLICIT NONE

type(oSLERP_T), INTENT(INOUT)  :: self 

call reportDestructor('oSLERP_T')

end subroutine oSLERP_destructor

!--------------------------------------------------------------------------
subroutine readNameList_(self, nmlfile, initonly)
!DEC$ ATTRIBUTES DLLEXPORT :: readNameList_
!! author: MDG 
!! version: 1.0 
!! date: 07/18/25
!!
!! read the namelist from an nml file for the oSLERP_T Class 

use mod_io 
use mod_EMsoft

IMPLICIT NONE 

class(oSLERP_T), INTENT(INOUT)        :: self
character(fnlen),INTENT(IN)           :: nmlfile
 !! full path to namelist file 
logical,OPTIONAL,INTENT(IN)           :: initonly
 !! fill in the default values only; do not read the file

type(EMsoft_T)                        :: EMsoft 
type(IO_T)                            :: Message       
logical                               :: skipread = .FALSE.

integer(kind=irg)                     :: framesize
integer(kind=irg)                     :: nthreads
real(kind=dbl)                        :: qm(4)
real(kind=dbl)                        :: mA(3)
real(kind=dbl)                        :: mC(3)
real(kind=dbl)                        :: o1(8)
real(kind=dbl)                        :: o2(8)
logical                               :: silentrender
real(kind=dbl)                        :: dOmega
character(fnlen)                      :: GBmode
character(fnlen)                      :: rendermode
character(fnlen)                      :: xtalname 
character(fnlen)                      :: povrayfile
character(fnlen)                      :: framefolder
character(fnlen)                      :: moviename
character(fnlen)                      :: PVincludepath
character(fnlen)                      :: PVexec 
character(fnlen)                      :: PVffmpeg
character(fnlen)                      :: metric

namelist /oSLERPlist/ framesize, nthreads, qm, mA, mC, o1, o2, dOmega, GBmode, xtalname, povrayfile, framefolder, &
                      rendermode, moviename, PVincludepath, PVexec, PVffmpeg, silentrender, metric

framesize = 1024
nthreads = 1
PVincludepath = 'undefined'
PVexec = 'undefined'
PVffmpeg = 'undefined'
! if GBmode = 'normal'
qm = (/ 1.D0, 0.D0, 0.D0, 0.D0 /)
mA = (/ 1.D0, 0.D0, 0.D0 /)
mC = (/ 1.D0, 0.D0, 0.D0 /)
! if GBmode = 'octonion'
o1 = (/ 1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0 /)   ! normalization will be done by program
o2 = (/ 1.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0 /)   ! normalization will be done by program

dOmega = 0.25D0
GBmode = 'normal'      ! 'normal' for (mA, qm) description; 'octonion' for (qA, qB) description
rendermode = 'cubes'
silentrender = .TRUE.
xtalname = 'undefined'
povrayfile = 'underfined'
framefolder = 'frames'
moviename = 'render.mp4'
metric = 'octonion'

if (present(initonly)) then
  if (initonly) skipread = .TRUE.
end if

if (.not.skipread) then
! read the namelist file
 open(UNIT=dataunit,FILE=trim(EMsoft%toNativePath(nmlfile)),DELIM='apostrophe',STATUS='old')
 read(UNIT=dataunit,NML=oSLERPlist)
 close(UNIT=dataunit,STATUS='keep')

! check for required entries
 if (trim(xtalname).eq.'undefined') then
  call Message%printError('EMoSLERP:',' xtal input file name is undefined in '//nmlfile)
 end if
 if (trim(povrayfile).eq.'undefined') then
  call Message%printError('EMoSLERP:',' POVray output file name is undefined in '//nmlfile)
 end if
end if

self%nml%framesize = framesize
self%nml%qm = qm
self%nml%mA = mA
self%nml%mC = mC
self%nml%o1 = o1
self%nml%o2 = o2
self%nml%qm = qm
self%nml%dOmega = dOmega
self%nml%silentrender = silentrender
self%nml%rendermode = rendermode
self%nml%xtalname = trim(xtalname)
self%nml%povrayfile = trim(povrayfile) 
self%nml%framefolder = trim(framefolder)
self%nml%moviename = trim(moviename)
self%nml%PVincludepath = trim(PVincludepath)
self%nml%PVexec = trim(PVexec)
self%nml%PVffmpeg = trim(PVffmpeg)
self%nml%metric = metric

end subroutine readNameList_

!--------------------------------------------------------------------------
function getNameList_(self) result(nml)
!DEC$ ATTRIBUTES DLLEXPORT :: getNameList_
!! author: MDG 
!! version: 1.0 
!! date: 07/18/25
!!
!! pass the namelist for the oSLERP_T Class to the calling program

IMPLICIT NONE 

class(oSLERP_T), INTENT(INOUT)          :: self
type(oSLERPNameListType)                :: nml

nml = self%nml

end function getNameList_

!--------------------------------------------------------------------------
subroutine oSLERP_(self, EMsoft, progname)
!DEC$ ATTRIBUTES DLLEXPORT :: oSLERP_
!! author: MDG 
!! version: 1.0 
!! date: 07/18/25
!!
!! performs the interpolation and creates all the required POVray files for rendering

use mod_EMsoft
use mod_HDFsupport
use mod_io 
use mod_crystallography
use mod_symmetry
use mod_rotations
use mod_quaternions
use mod_octonions
use mod_GBoctonions
use mod_povray
use omp_lib
use mod_dirstats
use mod_math
use mod_platformsupport

IMPLICIT NONE 

class(oSLERP_T), INTENT(INOUT)          :: self
type(EMsoft_T), INTENT(INOUT)           :: EMsoft
character(fnlen), INTENT(INOUT)         :: progname 

type(Cell_T)                            :: Cell 
type(SpaceGroup_T)                      :: SG
type(IO_T)                              :: Message 
type(PoVRay_T)                          :: PoVRay
type(Quaternion_T)                      :: qmc, rhoA, rhoA2, rhoC, rhoC2, qint, qint1, qint2, quat 
type(GBoctonion_T)                      :: oct1, oct2
type(Octonion_T)                        :: o
type(DirStat_T)                         :: DS
type(QuaternionArray_T)                 :: qAR
type(q_T)                               :: qu
type(o_T)                               :: om

integer(kind=irg)                       :: hdferr, pgnum, i, sgnum, numf, status, io_int(1) 
character(fnlen)                        :: fname, pvcmd, subfolder, povname, str, dirstring, xname
real(kind=dbl)                          :: Omega, dOmega, qn1(4), qn2(4), qmat(3,3), OB(3,3), ON(3,3), p, pA(4), pB(4), pC(4), &
                                           qinter(4), ointer(8), pD(4), qq(3), phiA, phiC, msA(3), msC(3), pp(3), io_real(1)
logical                                 :: dexists, fexists, frames_generated
real(kind=dbl),allocatable              :: tval(:)
character(4)                            :: number4
character(3)                            :: number3
character(2)                            :: number2
character(1)                            :: number1

call openFortranHDFInterface()

associate( nml=>self%nml )

! initialize the PoVRay class with a dummy file name
fname = 'dummy.txt'
PoVRay = PoVRay_T( EMsoft, fname, nofile=.TRUE. )

! In order for the PoVRay program to properly function, two include files are needed;
! here we copy those files from the EMsoftOO resources folder into the current folder. 
! the first file contains a set of PoVRay macros by F. Lohmueller,
! the other file contains specific macros for the oSLERP rendering
call PoVRay%get_incfile(EMsoft, 'analytical_g.inc')
call PoVRay%get_incfile(EMsoft, 'octonionSLERP.inc')

! get the point group number for this crystal structure
cell = cell_T()
xname = EMsoft%generateFilePath('EMXtalFolderpathname',trim(nml%xtalname))
call cell%readDataHDF(SG, EMsoft, useXtalName=xname)

pgnum = 0
sgnum = SG%getSpaceGroupNumber()
do i=1,32
  if (SGPG(i).le.sgnum) pgnum = i
end do

! if GBmode = 'normal', we must compute the relevant octonions
if (trim(nml%GBmode).eq.'normal') then 
    ! grains B and D are unrotated
    ! so A and C are rotated by the misorientation quaternion qm
    qmc = Quaternion_T( qd = nml%qm )
    qmc = conjg(qmc)
    ! determine the boundary normals in the sample reference frame
    msA = qmc%quat_Lp(nml%mA)
    msC = qmc%quat_Lp(nml%mC)
    ! get the rotations to bring the normals to the z-axis
    phiA = acos(abs(msA(3))) * 0.5D0
    pp = (/ msA(2), -msA(1), 0.D0 /)
    pp = sin(phiA) * pp / vecnorm(pp)
    rhoA = Quaternion_T( qd = (/ cos(phiA), pp(1), pp(2), pp(3) /) )
    phiC = acos(abs(msC(3))) * 0.5D0
    pp = (/ msC(2), -msC(1), 0.D0 /)
    pp = sin(phiC) * pp / vecnorm(pp)
    rhoC = Quaternion_T( qd = (/ cos(phiC), pp(1), pp(2), pp(3) /) )
    ! transform the quaternions and form the octonions
    rhoA2 = qmc * rhoA
    rhoC2 = qmc * rhoC
    oct1 = GBoctonion_T( rhoA2, rhoA )
    oct2 = GBoctonion_T( rhoC2, rhoC )
else
    o = Octonion_T( od = nml%o1 )
    oct1 = GBoctonion_T( oct = o )
    call oct1%oct_print(' input oct1 (normalized) : ')
    o = Octonion_T( od = nml%o2 )
    oct2 = GBoctonion_T( oct = o )
    call oct2%oct_print(' input oct2 (normalized) : ')
end if

! define the subfolder name for the frames
subfolder = EMsoft%generateFilePath('EMdatapathname',trim(nml%framefolder))

! make sure that the frame folder exists
pvcmd = trim(subfolder)
inquire(file=trim(subfolder),exist=dexists)
if (dexists.eqv..FALSE.) then 
    call Message%printMessage(' --> creating frames folder')
    status = system_mkdir(trim(subfolder))
    if (status.ne.0) call Message%printMessage(' --> error while creating frames folder')
else
    call Message%printMessage(' --> frames folder exists; removing existing files ')
    status = system_deletefile(trim(subfolder)//'/*')
    if (status.ne.0) call Message%printMessage(' --> error while deleting files in frames folder')
end if

! initialize the symmetry operators
DS = DirStat_T( PGnum = pgnum )

! compute the misorientation angle between the grain boundary octonions
Omega = oct1%GBO_Omega_Symmetric(oct2, DS, metric = trim(nml%metric))
io_real(1) = Omega * 180.D0 / cPi
call Message%WriteValue('--> octonion misorientation angle (degrees) ', io_real, 1)

if (Omega.eq.0.D0) then 
    call Message%printMessage(' The misorientation angles between the octonions is zero.')
    call Message%printMessage(' ---> No movie will be generated.')
    stop 'Program run aborted'
end if

! determine the number of movie frames
dOmega = nml%dOmega * cPi/180.D0
numf = nint(Omega/dOmega/2.0)
io_int(1) = numf
call Message%WriteValue('Number of movie frames requested = ',io_int,1)

! construct the t values array for the interpolation parameter
allocate(tval(0:numf+1))
do i=0,numf+1
    tval(i) = float(i)/float(numf+1)
end do

! convert the octonions to the correct parameters for POVray visualization
quat = conjg(oct1%GBO_get_q(1))
qn1 = quat%get_quatd()
qn1 = qn1/vecnorm(qn1)
quat = conjg(oct2%GBO_get_q(1))
qn2 = quat%get_quatd()
qn2 = qn2/vecnorm(qn2)

! create the PoVRay frame parameter files
do i=0,numf+1
    ! interpolate the octonions
    ointer = oct1%GBO_SLERP(oct1%get_octd(), oct2%get_octd(), Omega, tval(i), 8)
    qint1 = conjg( Quaternion_T( qd = ointer(1:4) ) )
    qint2 = Quaternion_T( qd = ointer(5:8) )
    qint = qint2 * qint1
    call qint%quat_normalize()
    qu = q_T( qdinp = qint%get_quatd() )
    om = qu%qo()
    qmat = om%o_copyd()
    OB = transpose(PoVRay%flipRotationMatrix(qmat))

    ! interpolate quaternions to get the grain boundary normal
    qinter = oct1%GBO_SLERP(qn1, qn2, Omega, tval(i), 4) 
    qu = q_T( qdinp = qinter(1:4) )
    om = qu%qo()
    qmat = om%o_copyd()
    ON = transpose(PoVRay%flipRotationMatrix(qmat))

    ! create the filename prefix
    if (numf.gt.1000) then 
        write (number4,"(I4.4)") i
        povname = trim(subfolder)//'/parameter'//number4//'.pov'
    else 
        write (number3,"(I3.3)") i
        povname = trim(subfolder)//'/parameter'//number3//'.pov'
    end if
        
    open(unit=dataunit,file=trim(povname),status='unknown',form='formatted')
    write(dataunit,"('// this file contains the orientations for grain B and for the interface normal')")
    write(dataunit,"('// this file can be generated by any script or code and must contain two orientation')")
    write(dataunit,"('// in the form of 4x4 matrices')")
    write(dataunit,"(' ')")
    write(dataunit,"('#macro GrainMatrix()  ')")
    write(dataunit,"('matrix < ',3(F10.6,','))") OB(1,1), OB(1,2), OB(1,3)
    write(dataunit,"(3(F10.6,','))") OB(2,1), OB(2,2), OB(2,3)
    write(dataunit,"(3(F10.6,','))") OB(3,1), OB(3,2), OB(3,3)
    write(dataunit,"('0.0, 0.0, 0.0 > ')")
    write(dataunit,"('#end')")
    write(dataunit,"(' ')")
    write(dataunit,"('#macro NormalMatrix()')")
    write(dataunit,"('matrix < ',3(F10.6,','))") ON(1,1), ON(1,2), ON(1,3)
    write(dataunit,"(3(F10.6,','))") ON(2,1), ON(2,2), ON(2,3)
    write(dataunit,"(3(F10.6,','))") ON(3,1), ON(3,2), ON(3,3)
    write(dataunit,"('0, 0, 0  > ')")
    write(dataunit,"('#end')")
    close(unit=dataunit,status='keep')
end do
call Message%printMessage(' --> All PoVRay frame parameter files created ')

! generate the main POVray scene file with the correct SubFolder variable
if (trim(nml%rendermode).eq.'spheres') then
    open(unit=dataunit,file='sphere-scene.pov',status='unknown',form='formatted')
else
    open(unit=dataunit,file='cube-scene.pov',status='unknown',form='formatted')
end if
write(dataunit,"('// POV-Ray 3.7 Scene File for visualization of grain boundary configurations')")
write(dataunit,"('// produced by EMsoftOO 6.0')")
write(dataunit,"('// email: degraef@cmu.edu')")
write(dataunit,"('// homepage: https://www.mse.engineering.cmu.edu/directory/bios/degraef-marc.html')")
write(dataunit,"('//')")
write(dataunit,"('')")
write(dataunit,"('// load all the definitions')")
write(dataunit,"('#version 3.7;')")
write(dataunit,"('#include ""octonionSLERP.inc""')")
write(dataunit,"('')")
write(dataunit,"('#declare SubFolder = ""',A,'/""')") trim(subfolder)
write(dataunit,"('')")
write(dataunit,"('#declare scl = 2.0;')")
write(dataunit,"('#declare he =  2.8;')")
write(dataunit,"('#declare FileName = concat(SubFolder,""/parameter"",str(clock,-3,0),"".pov"")')")
write(dataunit,"('#include FileName')")
write(dataunit,"('')")
write(dataunit,"('// the fixed grain A reference frame')")
write(dataunit,"('#declare scl = 3.0;')")
write(dataunit,"('#declare sclA = 3.0;')")
write(dataunit,"('#declare sclB = 2.0;')")

if (trim(nml%rendermode).eq.'spheres') then
    write(dataunit,"('// external coordinate axes (sample reference frame)')")
    write(dataunit,"('object{ AxisXYZ( scl, scl, scl, Texture_A_Red, Texture_A_Blue, Texture_A_Green) ')")
    write(dataunit,"('scale <0.5,0.5,0.5>  translate <-he, -he, he> }')")
    write(dataunit,"('')")
    write(dataunit,"('object{ AxisXYZgrainA( scl, scl, scl, Texture_A_Red, Texture_A_Blue,Texture_A_Green)  }')")
    write(dataunit,"('')")
    write(dataunit,"('// the rotatable grain B reference frame')")
    write(dataunit,"('#declare scl = 2.0;')")
    write(dataunit,"('object{ AxisXYZgrainB( 1.5*scl, 1.5*scl, 1.5*scl, Texture_A_Red, Texture_A_Blue, Texture_A_Green, scl)')") 
    write(dataunit,"('GrainMatrix() }')")
    write(dataunit,"('')")
    write(dataunit,"('// and the grain boundary plane normal')")
    write(dataunit,"('object { normals translate <0, R1, 0> NormalMatrix() }')")
else
    write(dataunit,"('#declare grainB = ')")
    write(dataunit,"('union {')")
    write(dataunit,"('difference {')")
    write(dataunit,"('object { AxisXYZcubeB( 1.5*sclB, 1.5*sclB, 1.5*sclB,Texture_A_Red,Texture_A_Blue,Texture_A_Green, sclB)')") 
    write(dataunit,"('GrainMatrix() }')")
    write(dataunit,"('object { TruncationB() NormalMatrix() }')")
    write(dataunit,"('}')")
    write(dataunit,"('object { normals NormalMatrix() }')")
    write(dataunit,"('} ')")
    write(dataunit,"('')")
    write(dataunit,"('#declare grainA = ')")
    write(dataunit,"('difference {')")
    write(dataunit,"('object { AxisXYZcubeA( sclA, sclA, sclA, Texture_A_Red, Texture_A_Blue,Texture_A_Green) }')")
    write(dataunit,"('object { TruncationA() NormalMatrix() }')")
    write(dataunit,"('}')")
    write(dataunit,"('')")
    write(dataunit,"('object { grainA }')")
    write(dataunit,"('object { grainB }')")
end if
close(unit=dataunit,status='keep')
call Message%printMessage(' --> PoVRay scene file sphere-scene.pov created ')

! next we generate the PoVRay.ini file with the rendering instructions
open(unit=dataunit,file='povray.ini',status='unknown',form='formatted')
if (trim(nml%rendermode).eq.'spheres') then
    write(dataunit,"('Input_File_Name=sphere-scene.pov')")
else
    write(dataunit,"('Input_File_Name=cube-scene.pov')")
end if
write(dataunit,"('Output_File_Name=',A,'/frame')") trim(subfolder)

str = '+L'//trim(nml%PVincludepath)
write(dataunit,"(A)") trim(str)

if (nml%framesize.ge.1000) then
   write(dataunit,"('+W',I4,' +H',I4)") nml%framesize, nml%framesize
else if (nml%framesize.ge.100) then 
   write(dataunit,"('+W',I3,' +H',I3)") nml%framesize, nml%framesize
else 
   write(dataunit,"('+W',I2,' +H',I2)") nml%framesize, nml%framesize
end if
write(dataunit,"('Initial_Clock=1')")
write(dataunit,"('Initial_Frame=1')")
if (numf.gt.1000) then 
    write(dataunit,"('Final_Clock=',I4)") numf+1
    write(dataunit,"('Final_Frame=',I4)") numf+1
else if (numf.gt.100) then
    write(dataunit,"('Final_Clock=',I3)") numf+1
    write(dataunit,"('Final_Frame=',I3)") numf+1
else if (numf.gt.10) then
    write(dataunit,"('Final_Clock=',I2)") numf+1
    write(dataunit,"('Final_Frame=',I2)") numf+1
else 
    write(dataunit,"('Final_Clock=',I1)") numf+1
    write(dataunit,"('Final_Frame=',I1)") numf+1
end if
write(dataunit,"('Work_Threads=6')")
close(unit=dataunit,status='keep')
call Message%printMessage(' --> povray.ini file created ')


! finally, execute the rendering programs if they can be found; 
! otherwise, print out some suggestions on how to proceed with the scene files
pvcmd = trim(nml%PVexec)
inquire(file=trim(pvcmd),exist=fexists)
frames_generated = .FALSE.

if (fexists.eqv..TRUE.) then
    call Message%printMessage('')
    if (nml%silentrender.eqv..TRUE.) then
      call Message%printMessage('Found PovRay command line executable; rendering frames in silent mode (may take a while)')
      call Message%printMessage('Executing '//trim(pvcmd)//' povray.ini >/dev/null 2>/dev/null')
      call system(trim(pvcmd)//' povray.ini >/dev/null 2>/dev/null')
    else
      call Message%printMessage('Found PovRay command line executable; rendering frames')
      call Message%printMessage('Executing '//trim(pvcmd)//' povray.ini')
      call system(trim(pvcmd)//' povray.ini')
    end if 
    frames_generated = .TRUE.
else
    call Message%printMessage(' =============================== ')
    call Message%printMessage('')
    call Message%printMessage('To render the individual frames, please run the command line PoVRay executable')
    call Message%printMessage('on the povray.ini file.  On the Mac platform, this can be accomplished as follows:')
    call Message%printMessage('')
    call Message%printMessage('shell> ...path.../PovrayCommandLineMacV2/Povray37UnofficialMacCmd povray >/dev/null 2>/dev/null')
    call Message%printMessage('')
    call Message%printMessage('Make sure that the second line of the povray.ini file has the correct path to the include folder.')
    call Message%printMessage('')
    if (trim(nml%rendermode).eq.'spheres') then
      call Message%printMessage('You can also simply start the interactive PoVRay program and render the sphere-scene.pov file.')
    else
      call Message%printMessage('You can also simply start the interactive PoVRay program and render the cube-scene.pov file.')
    end if
call Message%printMessage('')
end if

! if the frames have been generated, then look for the ffmpeg program; if found, try to execute it, otherwise print
! an explanation on what to do next 
pvcmd = trim(nml%PVffmpeg)
inquire(file=trim(pvcmd),exist=fexists)

if ((frames_generated.eqv..FALSE.).OR.(fexists.eqv..FALSE.)) then
    call Message%printMessage(' =============================== ')
    call Message%printMessage('')
    call Message%printMessage('The individual frames can be turned into a movie in a number of different ways;  if you have the ')
    call Message%printMessage('ffmpeg program installed, you can run it as follows to generate an mp4 file:')
    call Message%printMessage('')
    if (numf.gt.1000) then 
      str = 'shell> ...path.../ffmpeg -i frames/frame%04d.png -y -c:v libx264 -pix_fmt yuv420p '//trim(nml%moviename)
    else if (numf.gt.100) then  
      str = 'shell> ...path.../ffmpeg -i frames/frame%03d.png -y -c:v libx264 -pix_fmt yuv420p '//trim(nml%moviename)
    else if (numf.gt.10) then  
      str = 'shell> ...path.../ffmpeg -i frames/frame%02d.png -y -c:v libx264 -pix_fmt yuv420p '//trim(nml%moviename)
    else 
      str = 'shell> ...path.../ffmpeg -i frames/frame%01d.png -y -c:v libx264 -pix_fmt yuv420p '//trim(nml%moviename)
    end if
    call Message%printMessage(trim(str))
    call Message%printMessage('')
    call Message%printMessage('This should generate an mp4 file in your current folder.')
    call Message%printMessage('')
else if (fexists.eqv..TRUE.) then 
    if (numf.gt.1000) then 
     pvcmd = trim(pvcmd)//' -i frames/frame%04d.png -y -c:v libx264 -pix_fmt yuv420p '//trim(nml%moviename)
    else if (numf.gt.100) then  
     pvcmd = trim(pvcmd)//' -i frames/frame%03d.png -y -c:v libx264 -pix_fmt yuv420p '//trim(nml%moviename)
    else if (numf.gt.10) then  
     pvcmd = trim(pvcmd)//' -i frames/frame%02d.png -y -c:v libx264 -pix_fmt yuv420p '//trim(nml%moviename)
    else 
     pvcmd = trim(pvcmd)//' -i frames/frame%01d.png -y -c:v libx264 -pix_fmt yuv420p '//trim(nml%moviename)
    end if
    call Message%printMessage('')
    call Message%printMessage('Found ffmpeg executable; silently converting frames into mp4 file...')
    call Message%printMessage('Executing '//trim(pvcmd)//' >/dev/null 2>/dev/null')
    call system(trim(pvcmd)//' >/dev/null 2>/dev/null')
end if
end associate

call closeFortranHDFInterface()

end subroutine oSLERP_





end module mod_oSLERP