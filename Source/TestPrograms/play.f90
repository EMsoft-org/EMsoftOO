! ###################################################################
! Copyright (c) 2016-2024, Marc De Graef Research Group/Carnegie Mellon University
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

program EMplay
  !! author: MDG
  !! version: 1.0 
  !! date: 02/26/24
  !!
  !! test program to look at triple junctions and fundamental zones

use mod_kinds
use mod_global
use mod_EMsoft
use mod_rotations
use mod_quaternions
use mod_povray
use mod_so3
use mod_io
use mod_symmetry
use mod_crystallography
use HDF5
use mod_HDFsupport


IMPLICIT NONE

character(fnlen)        :: progname = 'EMplay.f90'
character(fnlen)        :: progdesc = 'test program'

type(EMsoft_T)          :: EMsoft
type(QuaternionArray_T) :: qAR, qsym, qdummy
type(Cell_T)            :: cell
type(SpaceGroup_T)      :: SG
type(so3_T)             :: SO
type(IO_T)              :: Message
type(q_T)               :: qu1, qu2, qu3, disor, disor12, disor13, disor23
type(a_T)               :: disax12, disax13, disax23
type(Quaternion_T)      :: q1, q2, q3, disor1
type(r_T)               :: roFZ

integer(kind=irg)       :: FZtype, FZorder
type(FZpointd),pointer  :: FZhead, FZtmp

character(fnlen)        :: xtalname, inputfile, fname
integer(kind=irg)       :: io_int(3), pgnum, norient, grain1, grain2, i1, i2, i3, nbins, bin 
real(kind=dbl)          :: quat(4)
integer(kind=irg),allocatable :: hist(:)

! the main idea of this program is to read in a list of orientations in the RFZ, pick one at random
! and generate all possible pairs of grain boundaries that have this one as grain 1. For each 
! misorientation we determine the disorientation, and then we combine these two into a new 
! misorientation for the third boundary.  The output of this program will then be a 3D plot
! of the corresponding disorientations; presumably/possibly there will be regions of the RFZ 
! that remain unpopulated ?


! print the EMsoft header and handle any command line arguments  
EMsoft = EMsoft_T( progname, progdesc )

call setRotationPrecision('d')

call openFortranHDFInterface()

! get the point group number from the crystal file
! xtalname = 'Ti-alpha.xtal'
! call cell%getCrystalData(xtalname, SG, EMsoft)

! define the symmetry operators
! pgnum = SG%getPGnumber()

pgnum = 6
call qdummy%QSym_init( pgnum, qsym )

call qsym%quat_print()

! define the fundamental zone class
SO = so3_T( pgnum )

! print some output
call SO%getFZtypeandorder(FZtype, FZorder)
io_int = (/ pgnum, FZtype, FZorder /)
call Message%WriteValue('  Point group, type, and order: ',io_int,3)

! read the input orientations
inputfile = 'play/orientations.txt'
fname = EMsoft%generateFilePath('EMdatapathname',inputfile)
call Message%printMessage(' - read orientation data from file '//trim(fname), frm="(/A)")
call SO%getOrientationsfromFile(fname)

norient = SO%getListCount('FZ')
io_int(1) = norient
call Message%WriteValue('   -> # orientations read from file : ', io_int, 1)

! apply the reduction to the RFZ
nbins = qsym%getQnumber()
allocate(hist(nbins))
call SO%ReducelisttoRFZ( qsym )

! transform them into an array of quaternions
call SO%listtoQuaternionArray( qAR, 'FZ' )

! ask for a number
! call Message%ReadValue(' enter a number for the grain 1 quaternion : ', io_int, 1)
! grain1 = io_int(1)
! q1 = qAR%getQuatfromArray(grain1)
! qu1 = q_T( qdinp = q1%get_quatd() )

! call Message%ReadValue(' enter a number for the grain 2 quaternion : ', io_int, 1)
! grain2 = io_int(1)
! q2 = qAR%getQuatfromArray(grain2)
! qu2 = q_T( qdinp = q2%get_quatd() )

! loop over all other orientations, avoiding grain1, and compute the 
! disorientation quaternions for grain boundaries 1-2 and 1-3; then 
! compute the disorientation for the third boundary and collect those 
! together in a new list that can then be visualized

! open(unit=dataunit,file='output.txt',status='unknown',form='formatted')
! write (dataunit,"(A2)") 'qu'
! write (dataunit,"(A3)") 'xxx'
do i1=1,norient 
! get the grain 2 orientation from the list
    q1 = qAR%getQuatfromArray(i1)
    ! qu1 = q_T( qdinp = q1%get_quatd() )
    do i2=1,norient 
        if (i1.ne.i2) then 
            q2 = qAR%getQuatfromArray(i2)
            ! qu2 = q_T( qdinp = q2%get_quatd() )
! and compute the disorientation between grain 1 and grain 2
            disor1 = q1 * conjg(q2)
            call disor1%quat_pos() 
            disor = q_T( qdinp = disor1%get_quatd() )
            call SO%ReduceOrientationtoRFZ(disor, qsym, roFZ, bin=bin)
            hist(bin) = hist(bin)+1
! and write this to the output file 
            ! quat = disor23%q_copyd()
            ! if (quat(1).lt.0.D0) quat = -quat
            ! write (dataunit, "(3(F10.7,','),F10.7)") quat(1:4)
        end if 
    end do
    ! if (mod(i2,250).eq.0) write (*,*) ' i2 = ', i2
  ! end if 
end do
! close(dataunit,status='keep')

call Message%printMessage(' computation completed.')

write (*,*) ' histogram count : ', hist

write (*,*) ' total # tested  : ', sum(hist)


end program EMplay
