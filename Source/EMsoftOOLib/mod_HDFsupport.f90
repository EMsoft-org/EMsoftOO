!###################################################################
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

module mod_HDFsupport
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! EMsoft HDF5 helper routines
  !!
  !! original history (roughly through version 5.x)
  !! 03/17/15 MDG 1.0 original
  !! 03/27/15 MDG 1.1 added integer and real write routines
  !! 03/29/15 MDG 1.2 removed all h5lt routines
  !! 03/31/15 MDG 1.3 added support for arrays of c_chars
  !! 04/07/15 MDG 1.4 added hyperslab routines for char, integer, float and double in 2, 3, and 4D
  !! 04/08/15 MDG 1.5 removed HDF_tail pointer as it was no longer needed
  !! 04/08/15 MDG 1.6 added optional overwrite keyword to all writeDataset... routines
  !! 01/06/15 MDG 1.7 changed order of array declarations so that dimensions are declared first (needed for Intel Fortran compiler)
  !! 10/31/16 MDG 1.8 added C_NULL_CHAR to every string that is an argument of a core HDF5 routine (via cstringify function)
  !! 11/01/16 MDG 1.9 fixed dimensional error in hyperslab read routines
  !! 12/14/16 MDG 2.0 added logical switch to flag DREAM.3D-generated files which require different string handling (fixed length vs variable length)
  !! 12/16/16 MDG 3.0 completely reworked HDF error handling; introduced h5open_EMsoft to initialize the fortran HDF interface
  !! 08/30/19 MDG 4.0 modified HDF_head definition for python f90wrap compatibility
  !! 09/30/19 MAJ 4.1 initial mods of allocations that caused Mac OSX/ifort issues in write routines 
  !! 10/01/19 MDG 4.2 additional mods to make ifort work on Mac OS X 
  !! 11/08/19 MDG 4.3 replaced individual dims parameters by single dims array in multiple routines

use mod_kinds
use mod_global
use HDF5
use stringconstants

IMPLICIT NONE 
  private

! type definition for HDF push-pop stack to keep track of open objects
  type HDFobjectStackType   
    character(LEN=1)                      :: objectType
    character(fnlen)                      :: objectName
    integer(HID_T)                        :: objectID
    type(HDFobjectStackType),pointer      :: next
  end type HDFobjectStackType

  logical, save                           :: FixedLengthflag, HDFverbose, dumpHDFstack
  logical, public, save                           :: interfaceOpen = .FALSE.

  public :: cstringify

  type, public :: HDF_T 
   !! HDF Class definition (this class takes over from the original self structure)
   private
     type(HDFobjectStackType)             :: head 

   contains 
   private 
! to simplify calling these routines as class methods, and simultaneously minimize the 
! amount of work involved in making these changes everywhere, we add an underscore to the name
! of the routine if there is no overload of the generic name; that then allows us to continue
! using the normal call scheme, but with HDF_foo replaced by HDF%foo everywhere. In the absence
! of procedure overload, the internal name of that routine will then be foo_. In cases where 
! there are several overloaded options, the 1D, 2D, etc portions of the routine names will 
! serve the role of the underscore.

! we start with the push-pop-stack management stuff and error handling routines 
     procedure, pass(self) :: push_
     procedure, pass(self) :: pop_
     procedure, pass(self) :: stackdump_
     procedure, pass(self) :: getObjectID_
     ! procedure, pass(self) :: toggleStackDump_
     procedure, pass(self) :: error_check_
     procedure, pass(self) :: associatedHead_
! open and closing routines for files, groups, and datasets  
     procedure, pass(self) :: createFile_
     procedure, pass(self) :: openFile_
     procedure, pass(self) :: createGroup_
     procedure, pass(self) :: openGroup_
     procedure, pass(self) :: openDataset_
! general purpose dataset writing routines 
     procedure, pass(self) :: writeNMLintegers_
     procedure, pass(self) :: writeNMLreals_
     procedure, pass(self) :: writeNMLdbles_
     procedure, pass(self) :: writeEMheader_
     procedure, pass(self) :: writeDatasetTextFile_
     procedure, pass(self) :: extractDatasetTextfile_
     procedure, pass(self) :: writeDatasetStringArray_
     procedure, pass(self) :: writeDatasetCharArray1D
     procedure, pass(self) :: writeDatasetCharArray2D
     procedure, pass(self) :: writeDatasetCharArray3D
     procedure, pass(self) :: writeDatasetCharArray4D
     procedure, pass(self) :: writeDatasetInteger_
     procedure, pass(self) :: writeDatasetInteger1byteArray1D_
     procedure, pass(self) :: writeDatasetIntegerArray1D
     procedure, pass(self) :: writeDatasetIntegerArray2D
     procedure, pass(self) :: writeDatasetIntegerArray3D
     procedure, pass(self) :: writeDatasetIntegerArray4D
     procedure, pass(self) :: writeDatasetFloat_
     procedure, pass(self) :: writeDatasetDouble_
     procedure, pass(self) :: writeDatasetFloatArray1D
     procedure, pass(self) :: writeDatasetFloatArray2D
     procedure, pass(self) :: writeDatasetFloatArray3D
     procedure, pass(self) :: writeDatasetFloatArray4D
     procedure, pass(self) :: writeDatasetFloatArray6D
     procedure, pass(self) :: writeDatasetDoubleArray1D
     procedure, pass(self) :: writeDatasetDoubleArray2D
     procedure, pass(self) :: writeDatasetDoubleArray3D
     procedure, pass(self) :: writeDatasetDoubleArray4D
     procedure, pass(self) :: writeHyperslabCharArray2D
     procedure, pass(self) :: writeHyperslabCharArray3D
     procedure, pass(self) :: writeHyperslabCharArray4D
     procedure, pass(self) :: writeHyperslabIntegerArray2D
     procedure, pass(self) :: writeHyperslabIntegerArray3D
     procedure, pass(self) :: writeHyperslabIntegerArray4D
     procedure, pass(self) :: writeHyperslabFloatArray2D
     procedure, pass(self) :: writeHyperslabFloatArray3D
     procedure, pass(self) :: writeHyperslabFloatArray4D
     procedure, pass(self) :: writeHyperslabDoubleArray2D
     procedure, pass(self) :: writeHyperslabDoubleArray3D
     procedure, pass(self) :: writeHyperslabDoubleArray4D
! general purpose dataset reading routines
     procedure, pass(self) :: readfromTextfile_
     procedure, pass(self) :: readDatasetStringArray_
     procedure, pass(self) :: readDatasetCharArray1D
     procedure, pass(self) :: readDatasetCharArray2D
     procedure, pass(self) :: readDatasetCharArray3D
     procedure, pass(self) :: readDatasetCharArray4D
     procedure, pass(self) :: readDatasetInteger_
     procedure, pass(self) :: readDatasetIntegerArray1D
     procedure, pass(self) :: readDatasetIntegerArray2D
     procedure, pass(self) :: readDatasetIntegerArray3D
     procedure, pass(self) :: readDatasetIntegerArray4D
     procedure, pass(self) :: readDatasetFloat_
     procedure, pass(self) :: readDatasetFloatArray1D
     procedure, pass(self) :: readDatasetFloatArray2D
     procedure, pass(self) :: readDatasetFloatArray3D
     procedure, pass(self) :: readDatasetFloatArray4D
     procedure, pass(self) :: readDatasetDouble_
     procedure, pass(self) :: readDatasetDoubleArray1D
     procedure, pass(self) :: readDatasetDoubleArray2D
     procedure, pass(self) :: readDatasetDoubleArray3D
     procedure, pass(self) :: readDatasetDoubleArray4D
     procedure, pass(self) :: readHyperslabCharArray2D_
     procedure, pass(self) :: readHyperslabCharArray3D_
     procedure, pass(self) :: readHyperslabCharArray4D_
     procedure, pass(self) :: readHyperslabIntegerArray2D_
     procedure, pass(self) :: readHyperslabIntegerArray3D_
     procedure, pass(self) :: readHyperslabIntegerArray4D_
     procedure, pass(self) :: readHyperslabFloatArray2D_
     procedure, pass(self) :: readHyperslabFloatArray3D_
     procedure, pass(self) :: readHyperslabFloatArray4D_
     procedure, pass(self) :: readHyperslabDoubleArray2D_
     procedure, pass(self) :: readHyperslabDoubleArray3D_
     procedure, pass(self) :: readHyperslabDoubleArray4D_
! general purpose attribute reading and writing routines 
     procedure, pass(self) :: addStringAttributeToGroup_
     procedure, pass(self) :: getStringAttributeFromGroup_
! EMsoft-specific predefined routines (e.g., writing the EMheader)
     ! procedure, pass(self) :: CrystalData
     ! procedure, pass(self) :: SaveDataHDF
     ! procedure, pass(self) :: ReadDataHDF
     procedure, pass(self) :: h5_write_pseudo_bse_image_
     procedure, pass(self) :: h5_tsl_read_ebsd_pattern_
     procedure, pass(self) :: h5_read_integer_dataset_
! miscellaneous auxiliary routines 
     procedure, pass(self) :: read2DImage_
     procedure, pass(self) :: CheckFixedLengthflag_
     procedure, pass(self) :: resetFixedLengthflag_
! finally, define the destructor routine 
     final :: HDF_destructor

! generic (public) function definitions and overloads

     generic, public :: push => push_
     generic, public :: pop => pop_
     generic, public :: stackdump => stackdump_
     generic, public :: getObjectID => getObjectID_
     generic, public :: associatedHead => associatedHead_
     ! generic, public :: togglestackdump => togglestackdump_
 
     generic, public :: createFile => createFile_
     generic, public :: openFile => openFile_
     generic, public :: createGroup => createGroup_
     generic, public :: openGroup => openGroup_
     generic, public :: openDataset => openDataset_
     generic, public :: error_check => error_check_
 
     generic, public :: writeNMLintegers => writeNMLintegers_
     generic, public :: writeNMLreals => writeNMLreals_
     generic, public :: writeNMLdbles => writeNMLdbles_
     generic, public :: writeEMheader => writeEMheader_
     generic, public :: writeDatasetTextFile => writeDatasetTextFile_
     generic, public :: extractDatasetTextfile => extractDatasetTextfile_
     generic, public :: writeDatasetStringArray => writeDatasetStringArray_
     generic, public :: writeDatasetCharArray => writeDatasetCharArray1D, writeDatasetCharArray2D, &
                                                 writeDatasetCharArray3D, writeDatasetCharArray4D
     generic, public :: writeDatasetInteger => writeDatasetInteger_
     generic, public :: writeDatasetInteger1byteArray1D => writeDatasetInteger1byteArray1D_
     generic, public :: writeDatasetIntegerArray => writeDatasetIntegerArray1D, writeDatasetIntegerArray2D, &
                                                    writeDatasetIntegerArray3D, writeDatasetIntegerArray4D
     generic, public :: writeDatasetFloat => writeDatasetFloat_
     generic, public :: writeDatasetDouble => writeDatasetDouble_
     generic, public :: writeDatasetFloatArray => writeDatasetFloatArray1D, writeDatasetFloatArray2D, &
                                                  writeDatasetFloatArray3D, writeDatasetFloatArray4D, &
                                                  writeDatasetFloatArray6D
     generic, public :: writeDatasetDoubleArray => writeDatasetDoubleArray1D, writeDatasetDoubleArray2D, &
                                                   writeDatasetDoubleArray3D, writeDatasetDoubleArray4D
     generic, public :: writeHyperslabCharArray => writeHyperslabCharArray2D, writeHyperslabCharArray3D, &
                                                   writeHyperslabCharArray4D
     generic, public :: writeHyperslabIntegerArray => writeHyperslabIntegerArray2D, writeHyperslabIntegerArray3D, &
                                                      writeHyperslabIntegerArray4D
     generic, public :: writeHyperslabFloatArray => writeHyperslabFloatArray2D,  writeHyperslabFloatArray3D, &
                                                    writeHyperslabFloatArray4D
     generic, public :: writeHyperslabDoubleArray => writeHyperslabDoubleArray2D, writeHyperslabDoubleArray3D, &
                                                     writeHyperslabDoubleArray4D
     generic, public :: readfromTextfile => readfromTextfile_
     generic, public :: readDatasetStringArray => readDatasetStringArray_
     generic, public :: readDatasetCharArray => readDatasetCharArray1D, readDatasetCharArray2D, &
                                                readDatasetCharArray3D, readDatasetCharArray4D
     generic, public :: readDatasetInteger => readDatasetInteger_
     generic, public :: readDatasetIntegerArray => readDatasetIntegerArray1D, readDatasetIntegerArray2D, &
                                                   readDatasetIntegerArray3D, readDatasetIntegerArray4D
     generic, public :: readDatasetFloat => readDatasetFloat_
     generic, public :: readDatasetFloatArray => readDatasetFloatArray1D, readDatasetFloatArray2D, &
                                                 readDatasetFloatArray3D, readDatasetFloatArray4D
     generic, public :: readDatasetDouble => readDatasetDouble_
     generic, public :: readDatasetDoubleArray => readDatasetDoubleArray1D, readDatasetDoubleArray2D, &
                                                  readDatasetDoubleArray3D, readDatasetDoubleArray4D
     generic, public :: readHyperslabCharArray2D => readHyperslabCharArray2D_ 
     generic, public :: readHyperslabCharArray3D => readHyperslabCharArray3D_
     generic, public :: readHyperslabCharArray4D => readHyperslabCharArray4D_
     generic, public :: readHyperslabIntegerArray2D => readHyperslabIntegerArray2D_
     generic, public :: readHyperslabIntegerArray3D =>readHyperslabIntegerArray3D_
     generic, public :: readHyperslabIntegerArray4D =>readHyperslabIntegerArray4D_
     generic, public :: readHyperslabFloatArray2D => readHyperslabFloatArray2D_
     generic, public :: readHyperslabFloatArray3D => readHyperslabFloatArray3D_
     generic, public :: readHyperslabFloatArray4D =>  readHyperslabFloatArray4D_
     generic, public :: readHyperslabDoubleArray2D => readHyperslabDoubleArray2D_
     generic, public :: readHyperslabDoubleArray3D => readHyperslabDoubleArray3D_
     generic, public :: readHyperslabDoubleArray4D => readHyperslabDoubleArray4D_
     generic, public :: addStringAttributeToGroup => addStringAttributeToGroup_
     generic, public :: getStringAttributeFromGroup => getStringAttributeFromGroup_
     generic, public :: read2DImage => read2DImage_ 
     ! generic, public :: CrystalData => CrystalData_
     ! generic, public :: SaveDataHDF => SaveDataHDF_
     ! generic, public :: ReadDataHDF => ReadDataHDF_
     generic, public :: h5_write_pseudo_bse_image => h5_write_pseudo_bse_image_
     generic, public :: h5_tsl_read_ebsd_pattern => h5_tsl_read_ebsd_pattern_
     generic, public :: h5_read_integer_dataset => h5_read_integer_dataset_
     generic, public :: CheckFixedLengthflag => CheckFixedLengthflag_
     generic, public :: resetFixedLengthflag => resetFixedLengthflag_

  end type HDF_T

!DEC$ ATTRIBUTES DLLEXPORT :: push
!DEC$ ATTRIBUTES DLLEXPORT :: pop
!DEC$ ATTRIBUTES DLLEXPORT :: stackdump
!DEC$ ATTRIBUTES DLLEXPORT :: getObjectID
! !DEC$ ATTRIBUTES DLLEXPORT :: togglestackdump
!DEC$ ATTRIBUTES DLLEXPORT :: error_check
!DEC$ ATTRIBUTES DLLEXPORT :: createFile
!DEC$ ATTRIBUTES DLLEXPORT :: openFile
!DEC$ ATTRIBUTES DLLEXPORT :: openFortranHDFInterface
!DEC$ ATTRIBUTES DLLEXPORT :: closeFortranHDFInterface
!DEC$ ATTRIBUTES DLLEXPORT :: createGroup
!DEC$ ATTRIBUTES DLLEXPORT :: openGroup
!DEC$ ATTRIBUTES DLLEXPORT :: openDataset
!DEC$ ATTRIBUTES DLLEXPORT :: writeDatasetTextFile
!DEC$ ATTRIBUTES DLLEXPORT :: extractDatasetTextfile
!DEC$ ATTRIBUTES DLLEXPORT :: writeDatasetStringArray
!DEC$ ATTRIBUTES DLLEXPORT :: writeDatasetCharArray
!DEC$ ATTRIBUTES DLLEXPORT :: writeDatasetInteger
!DEC$ ATTRIBUTES DLLEXPORT :: writeDatasetInteger1byteArray1D
!DEC$ ATTRIBUTES DLLEXPORT :: writeDatasetIntegerArray
!DEC$ ATTRIBUTES DLLEXPORT :: writeDatasetFloat
!DEC$ ATTRIBUTES DLLEXPORT :: writeDatasetDouble
!DEC$ ATTRIBUTES DLLEXPORT :: writeDatasetFloatArray
!DEC$ ATTRIBUTES DLLEXPORT :: writeDatasetDoubleArray
!DEC$ ATTRIBUTES DLLEXPORT :: writeHyperslabCharArray
!DEC$ ATTRIBUTES DLLEXPORT :: writeHyperslabIntegerArray
!DEC$ ATTRIBUTES DLLEXPORT :: writeHyperslabFloatArray
!DEC$ ATTRIBUTES DLLEXPORT :: writeHyperslabDoubleArray
!DEC$ ATTRIBUTES DLLEXPORT :: readfromTextfile
!DEC$ ATTRIBUTES DLLEXPORT :: readDatasetStringArray
!DEC$ ATTRIBUTES DLLEXPORT :: readDatasetCharArray
!DEC$ ATTRIBUTES DLLEXPORT :: readDatasetInteger
!DEC$ ATTRIBUTES DLLEXPORT :: readDatasetIntegerArray
!DEC$ ATTRIBUTES DLLEXPORT :: readDatasetFloat
!DEC$ ATTRIBUTES DLLEXPORT :: readDatasetFloatArray
!DEC$ ATTRIBUTES DLLEXPORT :: readDatasetDouble
!DEC$ ATTRIBUTES DLLEXPORT :: readDatasetDoubleArray
!DEC$ ATTRIBUTES DLLEXPORT :: readHyperslabCharArray2D
!DEC$ ATTRIBUTES DLLEXPORT :: readHyperslabCharArray3D
!DEC$ ATTRIBUTES DLLEXPORT :: readHyperslabCharArray4D
!DEC$ ATTRIBUTES DLLEXPORT :: readHyperslabIntegerArray2D
!DEC$ ATTRIBUTES DLLEXPORT :: readHyperslabIntegerArray3D
!DEC$ ATTRIBUTES DLLEXPORT :: readHyperslabIntegerArray4D
!DEC$ ATTRIBUTES DLLEXPORT :: readHyperslabFloatArray2D
!DEC$ ATTRIBUTES DLLEXPORT :: readHyperslabFloatArray3D
!DEC$ ATTRIBUTES DLLEXPORT :: readHyperslabFloatArray4D
!DEC$ ATTRIBUTES DLLEXPORT :: readHyperslabDoubleArray2D
!DEC$ ATTRIBUTES DLLEXPORT :: readHyperslabDoubleArray3D
!DEC$ ATTRIBUTES DLLEXPORT :: readHyperslabDoubleArray4D
!DEC$ ATTRIBUTES DLLEXPORT :: addStringAttributeToGroup
!DEC$ ATTRIBUTES DLLEXPORT :: getStringAttributeFromGroup
! !DEC$ ATTRIBUTES DLLEXPORT :: CrystalData
! !DEC$ ATTRIBUTES DLLEXPORT :: SaveDataHDF
! !DEC$ ATTRIBUTES DLLEXPORT :: ReadDataHDF
!DEC$ ATTRIBUTES DLLEXPORT :: h5_write_pseudo_bse_image
!DEC$ ATTRIBUTES DLLEXPORT :: h5_tsl_read_ebsd_pattern
!DEC$ ATTRIBUTES DLLEXPORT :: h5_read_integer_dataset
!DEC$ ATTRIBUTES DLLEXPORT :: CheckFixedLengthflag
!DEC$ ATTRIBUTES DLLEXPORT :: resetFixedLengthflag

! the constructor routine for this class 
  interface HDF_T
    module procedure HDF_constructor
  end interface HDF_T

public :: openFortranHDFInterface, closeFortranHDFInterface 

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! We begin with the functions/subroutines that are public in this class
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
type(HDF_T) function HDF_constructor( ) result(HDF)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! constructor for the HDF Class 

use mod_io 

IMPLICIT NONE

  type(IO_T)                    :: Message
  integer(kind=irg)             :: hdferr
  integer(kind=irg)             :: printonoff

! silent HDF operation or somewhat verbose ?
  HDFverbose = .FALSE. 
  dumpHDFstack = .FALSE.     ! set with the switchStackDump method

! make sure that the fortran HDF interface is open; if not call an error 
  if (HDFinterfaceOpen.eqv..FALSE.) then 
    call Message%printError('HDF_constructor', 'The HDF fortran interface is closed; must be opened in calling program')
  end if 

! then we nullify the HDF%head private variable (again, the user does not get direct access
! to this linked list... it is used internally and for debugging purposes).
  nullify(HDF%head%next)

end function HDF_constructor
  
!--------------------------------------------------------------------------
subroutine HDF_destructor(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! destructor (final) for the HDF Class 

  type(HDF_T)                         :: self 

  integer(kind=irg)                   :: hdferr
  integer(kind=irg)                   :: printonoff
  type(HDFobjectStackType), pointer   :: tmp

  call reportDestructor('HDF_T')

! for now, there is nothing to destruct since the push-pop stack has already been 
! nullified in the pop(.TRUE.) call.

end subroutine HDF_destructor

!--------------------------------------------------------------------------
subroutine openFortranHDFInterface(verbose)
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/13/20
  !!
  !! open the HDF interface 

use mod_io 

IMPLICIT NONE 

  logical, INTENT(IN), OPTIONAL       :: verbose

  type(IO_T)                          :: Message 
  integer(kind=irg)                   :: hdferr, io_int(1)
  integer(kind=irg)                   :: printonoff

! open the interface; note that we can not use the error_check routine here
! because the HDF_T class has not been instantiated yet !!!
    call h5open_f(hdferr)
    if (hdferr.lt.0) then 
      io_int(1) = hdferr 
      call Message%WriteValue('Error code : ',io_int,1)
      call Message%printMessage('   returned by routine h5open_f',frm="(A)")
    end if 

! and turn standard error reporting off; we'll handle errors our way...
    printonoff = 0
    call h5eset_auto_f(printonoff,hdferr)
    if (hdferr.lt.0) then 
      io_int(1) = hdferr 
      call Message%WriteValue('Error code : ',io_int,1)
      call Message%printMessage('   returned by routine h5eset_auto_f',frm="(A)")
    end if 

    HDFinterfaceOpen = .TRUE.

    if (present(verbose)) then 
        if (verbose.eqv..TRUE.) call Message%printMessage(' *** Fortran HDF interface OPEN *** ')
    end if

end subroutine openFortranHDFInterface

!--------------------------------------------------------------------------
subroutine closeFortranHDFInterface(verbose)
  !! author: MDG 
  !! version: 1.0 
  !! date: 02/13/20
  !!
  !! close the HDF interface 

use mod_io 

IMPLICIT NONE 

  logical, INTENT(IN), OPTIONAL       :: verbose

  type(IO_T)                          :: Message 
  integer(kind=irg)                   :: hdferr, io_int(1)
  integer(kind=irg)                   :: printonoff

! turn standard error reporting back on; note that we can not use the error_check routine here
! because the HDF_T class has not been instantiated yet !!! 
    printonoff = 1
    call h5eset_auto_f(printonoff, hdferr)
    if (hdferr.lt.0) then 
      io_int(1) = hdferr 
      call Message%WriteValue('Error code : ',io_int,1)
      call Message%printMessage('   returned by routine h5eset_auto_f',frm="(A)")
    end if 

! close the HDF fortran interface
    call h5close_f(hdferr)
     if (hdferr.lt.0) then 
      io_int(1) = hdferr 
      call Message%WriteValue('Error code : ',io_int,1)
      call Message%printMessage('   returned by routine h5close_f',frm="(A)")
    end if 

    HDFinterfaceOpen = .FALSE.

    if (present(verbose)) then 
        if (verbose.eqv..TRUE.) call Message%printMessage(' *** Fortran HDF interface CLOSED *** ')
    end if

end subroutine closeFortranHDFInterface

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! to avoid user usage of file, group, dataset, etc, IDs, we use a push_-pop_
! stack to keep track of the open items and close them again.  
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine push_(self, oT, oID, oName)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! push_ an HDF object to the stack

use mod_io 

IMPLICIT NONE

class(HDF_T), INTENT(INOUT)              :: self
character(LEN=1),INTENT(IN)             :: oT
 !! object type character
integer(HID_T),INTENT(IN)               :: oID 
 !! oID object identifier
character(fnlen),INTENT(IN)             :: oName
 !! oName name

type(HDFobjectStackType),pointer        :: node
integer(kind=irg)                       :: istat
type(IO_T)                              :: Message

! the stack always exists but we never use the top level
if (.not.associated(self%head%next)) then 
   allocate(self%head%next,stat=istat)     ! allocate new value
   call error_check_(self, 'push_: unable to allocate self pointer', istat, .TRUE.)

   nullify(self%head%next%next)            ! nullify next in list
   if (HDFverbose.eqv..TRUE.) call Message%printMessage('  -> creating self linked list', frm = "(A)")
else
   allocate(node,stat=istat)          ! allocate new value
   call error_check_(self, 'push_: unable to allocate node pointer', istat, .TRUE.)

   node%next => self%head%next
   self%head%next => node
end if

! set the values
self%head%next%objectType = oT
self%head%next%objectID   = oID
self%head%next%objectname = trim(oName)

if (dumpHDFstack.eqv..TRUE.) call stackdump_(self)

end subroutine push_

!--------------------------------------------------------------------------
recursive subroutine pop_(self, closeall)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! pop_ an HDF object from the stack and close it

use mod_io

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
!f2py intent(in,out) ::  self
logical,INTENT(IN),optional                :: closeall

integer                                    :: error, istat
type(HDFobjectStackType),pointer           :: tmp

nullify(tmp)

if (PRESENT(closeall)) then  
! this would be called if an error arises that forces a complete shutdown of the program, or at the end of a regular program
  do while (associated(self%head%next)) 
! close the current object 
    error = HDF_close_level(self%head%next%objectType, self%head%next%objectID)
    call error_check_(self, 'pop_:unable to close requested level for object type '//self%head%next%objectType, error,.TRUE.)

! and re-point the stack head
    tmp => self%head%next
    self%head%next => self%head%next%next
! delete the old entry
    deallocate(tmp)
  end do
  nullify(self%head%next)
else
! close the current object 
  error = HDF_close_level(self%head%next%objectType, self%head%next%objectID)
  call error_check_(self, 'pop_:unable to close requested level for object type '//self%head%next%objectType, error,.TRUE.)

! and re-point the stack head
  tmp => self%head%next
  self%head%next => self%head%next%next
! delete the old entry
  deallocate(tmp)

  if (dumpHDFstack.eqv..TRUE.) call stackdump_(self)
end if

contains

  recursive function HDF_close_level(oT, oID) result(error)
  
  IMPLICIT NONE

  character(LEN=1),INTENT(IN)   :: oT
  integer(HID_T),INTENT(IN)     :: oID 
  integer(kind=irg)             :: error

  select case(oT)
  case ('f') 
    call h5fclose_f(oID, error)  ! close the file
    call error_check_(self, 'pop_:h5fclose_f', error)


  case ('g') 
    call h5gclose_f(oID, error)  ! close the group
    call error_check_(self, 'pop_:h5gclose_f', error)


  case ('d') 
    call h5dclose_f(oID, error)  ! close the data set
    call error_check_(self, 'pop_:h5dclose_f', error)


  case ('a') 
    call h5aclose_f(oID, error)  ! close the attribute
    call error_check_(self, 'pop_:h5aclose_f', error)


  case ('t') 
    call h5tclose_f(oID, error)  ! close the data type
    call error_check_(self, 'pop_:h5tclose_f', error)


  case ('s') 
    call h5sclose_f(oID, error)  ! close the data space
    call error_check_(self, 'pop_:h5sclose_f', error)


  case DEFAULT
  end select

end function HDF_close_level

end subroutine pop_

!--------------------------------------------------------------------------
recursive subroutine stackdump_(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! print out the entire stack for debugging purposes

use mod_io

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)              :: self

type(HDFobjectStackType),pointer        :: tmp
integer(kind=irg)                       :: io_int(1)
type(IO_T)                              :: Message

tmp => self%head%next
if (.not.associated(tmp)) then
  call Message%WriteValue('stackdump_:',' stack is empty')

else
  call Message%WriteValue('stackdump_:','HDF stack entries')

  do
    if (.not.associated(tmp)) EXIT
    call Message%WriteValue('','>'//tmp%objectType//'<  >'//trim(tmp%objectName)//'<', frm = "(A)",advance="no") 

    io_int(1) = tmp%objectID
    call Message%WriteValue('',io_int,1,frm="(I12)")

    tmp => tmp%next
  end do
end if

end subroutine stackdump_

!--------------------------------------------------------------------------
recursive function associatedHead_(self) result(assoc)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/14/20
  !!
  !! check whether or not the head pointer is associated

class(HDF_T), INTENT(INOUT)     :: self 
logical                         :: assoc

assoc = .FALSE.
if (associated(self%head%next)) then 
  assoc = .TRUE.
end if 

end function associatedHead_

!--------------------------------------------------------------------------
recursive function getObjectID_(self) result(ObjID)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/15/20
  !!
  !! returns an object identifier

class(HDF_T), INTENT(INOUT)     :: self 
integer(HID_T)                  :: ObjID

ObjID = 0
if (associated(self%head%next)) then 
  ObjID = self%head%next%ObjectID
end if 

end function getObjectID_


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! from here on, we have basic HDF support routines
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive function createFile_(self, HDFname) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! Create a new HDF file (this also opens the file)

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: HDFname

integer(kind=irg)                          :: success
integer(HID_T)                             :: file_id ! file identifier
integer                                    :: hdferr  ! hdferr flag

success = 0

! Create a new file using default properties.
call h5fcreate_f(cstringify(HDFname), H5F_ACC_TRUNC_F, file_id, hdferr)
call error_check_(self, 'createFile_:h5fcreate_f:'//trim(HDFname), hdferr)

if (hdferr.lt.0) then
  success = -1
else
  call self%push_('f', file_id, HDFname)
end if

FixedLengthflag = .FALSE.

end function createFile_

!--------------------------------------------------------------------------
recursive function openFile_(self, HDFname, readonly) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! Open an HDF file

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: HDFname

logical,INTENT(IN),OPTIONAL                :: readonly
integer(kind=irg)                          :: success

integer(HID_T)                             :: file_id ! file identifier
integer                                    :: hdferr  ! hdferr flag

success = 0
if (present(readonly)) then 
  call H5Fopen_f(cstringify(HDFname), H5F_ACC_RDONLY_F, file_id, hdferr)
  call error_check_(self, 'openFile_:H5Fopen_f:'//trim(HDFname)//':readonly', hdferr)
else
  call H5Fopen_f(cstringify(HDFname), H5F_ACC_RDWR_F, file_id, hdferr)
  call error_check_(self, 'openFile_:H5Fopen_f:'//trim(HDFname), hdferr)
end if

if (hdferr.lt.0) then 
  success = -1
else
  call self%push_('f', file_id, HDFname)
end if

FixedLengthflag = .FALSE.

end function openFile_

!--------------------------------------------------------------------------
recursive function createGroup_(self, groupname) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! create and open a new group; if it already exists, just open it

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: groupname

integer(kind=irg)                          :: success
integer(HID_T)                             :: group_id!  identifier
integer                                    :: hdferr  ! hdferr flag
logical                                    :: g_exists

success = 0

call H5Lexists_f(self%head%next%objectID,cstringify(groupname),g_exists, hdferr)
call error_check_(self, 'createGroup_:H5Lexists_f:'//trim(groupname), hdferr)

if (g_exists) then 
  call H5Gopen_f(self%head%next%objectID, cstringify(groupname), group_id, hdferr)
  call error_check_(self, 'createGroup_:H5Gopen_f:'//trim(groupname), hdferr)

  if (hdferr.lt.0) then
    success = -1
  else
  ! put the group_id onto the HDF_stack
    call self%push_('g', group_id, groupname)
  end if
else 
  call H5Gcreate_f(self%head%next%objectID, cstringify(groupname), group_id, hdferr)
  call error_check_(self, 'createGroup_:H5Gcreate_f:'//trim(groupname), hdferr)

  if (hdferr.lt.0) then
    success = -1
  else
  ! and put the group_id onto the HDF_stack
    call self%push_('g', group_id, groupname)
  end if
end if

end function createGroup_

!--------------------------------------------------------------------------
recursive function openGroup_(self, groupname) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! open an existing group

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: groupname

integer(kind=irg)                          :: success
integer(HID_T)                             :: group_id !  identifier
integer                                    :: hdferr  ! hdferr flag

success = 0

call H5Gopen_f(self%head%next%objectID, cstringify(groupname), group_id, hdferr)
call error_check_(self, 'openGroup_:H5Gopen_f:'//trim(groupname), hdferr)

if (hdferr.lt.0) then
  success = -1
else
! put the group_id onto the HDF_stack
  call self%push_('g', group_id, groupname)
end if

end function openGroup_

!--------------------------------------------------------------------------
recursive function openDataset_(self, dataname) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! open an existing dataset

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname

integer(kind=irg)                          :: success
integer(HID_T)                             :: data_id !  identifier
integer                                    :: hdferr  ! hdferr flag

success = 0

call H5dopen_f(self%head%next%objectID, cstringify(dataname), data_id, hdferr)
call error_check_(self, 'openDataset_:H5dopen_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
else
! put the data_id onto the HDF_stack
  call self%push_('d', data_id, dataname)
end if

end function openDataset_

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! various text writing routines 
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive function writeDatasetTextFile_(self, dataname, filename) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! write a text file data set to the current file or group ID

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                              :: self
character(fnlen),INTENT(IN)                             :: dataname
character(fnlen),INTENT(IN)                             :: filename

integer(kind=irg)                                       :: success

character(len=fnlen, KIND=c_char),allocatable, TARGET   :: stringarray(:) 
integer(kind=irg)                                       :: nlines

integer(HSIZE_T)                                        :: dim0 
integer(SIZE_T)                                         :: sdim 
integer(HID_T)                                          :: filetype, space, dset ! Handles
integer                                                 :: hdferr, i, rnk
integer(HSIZE_T), DIMENSION(1:1)                        :: dims
integer(HSIZE_T), DIMENSION(1:2)                        :: maxdims
logical                                                 :: g_exists

TYPE(C_PTR), ALLOCATABLE, TARGET                        :: wdata(:)
TYPE(C_PTR)                                             :: f_ptr

success = 0

stringarray = readfromTextfile_(self, filename, nlines) 

! first, convert the stringarray to an array of C-pointers, with each string
! terminated by a C_NULL_CHAR.
dims(1) = nlines
allocate(wdata(1:dims(1)))
do i=1,dims(1)
  wdata(i) = C_LOC(stringarray(i))
end do

! then we write this C_ptr to the HDF file in the proper data set

! first create the memory data type (filetype)
call H5Tcopy_f(H5T_STRING, filetype, hdferr)
call error_check_(self, 'writeDatasetTextFile_:H5Tcopy_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetTextFile_:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
call H5Lexists_f(self%head%next%objectID,cstringify(dataname),g_exists, hdferr)
call error_check_(self, 'writeDatasetTextFile_:H5Lexists_f:'//trim(dataname), hdferr)

if (g_exists) then 
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetTextFile_:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), filetype, space, dset, hdferr)
  call error_check_(self, 'writeDatasetTextFile_:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

f_ptr = C_LOC(wdata(1))
call h5dwrite_f(dset, filetype, f_ptr, hdferr)
call error_check_(self, 'writeDatasetTextFile_:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetTextFile_:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetTextFile_:h5sclose_f:'//trim(dataname), hdferr)

call h5tclose_f(filetype, hdferr)
call error_check_(self, 'writeDatasetTextFile_:h5tclose_f:'//trim(dataname), hdferr)

deallocate(wdata)

! that's it

end function writeDatasetTextFile_

!--------------------------------------------------------------------------
recursive subroutine readDatasetStringArray_(self, dataname, nlines, hdferr, stringarray)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! read a string array from a data set

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                              :: self
character(fnlen),INTENT(IN)                             :: dataname
integer(kind=irg),INTENT(OUT)                           :: nlines

integer(kind=irg),INTENT(OUT)                           :: hdferr
character(len=fnlen, KIND=c_char),allocatable, TARGET, INTENT(OUT)   :: stringarray(:) 

integer(HID_T)                                          :: filetype, space, memtype ! Handles
integer                                                 :: i, length
integer(HSIZE_T), DIMENSION(1:1)                        :: dims
integer(HSIZE_T), DIMENSION(1:1)                        :: maxdims
integer(SIZE_T)                                         :: size

character(len = fnlen, kind=c_char),  POINTER           :: pfstr ! A pointer to a Fortran string
TYPE(C_PTR), DIMENSION(:), ALLOCATABLE, TARGET          :: rdata ! Read buffer
character(len=fnlen), TARGET                            :: fl_rdata 
TYPE(C_PTR)                                             :: f_ptr

! Open dataset.
!
hdferr = self%openDataset_(dataname)
!
! Get the datatype.
!
call H5Dget_type_f(self%head%next%objectID, filetype, hdferr)
call error_check_(self, 'readDatasetStringArray_:H5Dget_type_f:'//trim(dataname), hdferr)

if (FixedLengthflag.eqv..TRUE.) then ! this option is only set up to read one single string into stringarray...
  call H5Tget_size_f(filetype, size, hdferr)
  call error_check_(self, 'readDatasetStringArray_:H5Tget_size_f:'//trim(dataname), hdferr)

! Get dataspace and allocate memory for read buffer.
!
  call H5Dget_space_f(self%head%next%objectID, space, hdferr)
  call error_check_(self, 'readDatasetStringArray_:H5Dget_space_f:'//trim(dataname), hdferr)
 
  call H5Tcopy_f(H5T_FORTRAN_S1, memtype, hdferr)
  call error_check_(self, 'readDatasetStringArray_:H5Tcopy_f:'//trim(dataname), hdferr)
  
  call H5Tset_size_f(memtype, size-1, hdferr)
  call error_check_(self, 'readDatasetStringArray_:H5Tset_size_f:'//trim(dataname), hdferr)

!
! Read the data.
!
  f_ptr = C_LOC(fl_rdata(1:1))
  call h5dread_f(self%head%next%objectID, memtype, f_ptr, hdferr) !, space)
  call error_check_(self, 'readDatasetStringArray_:h5dread_f:'//trim(dataname), hdferr)

  allocate(stringarray(1))
  do i=1,size-1
    stringarray(1)(i:i) = fl_rdata(i:i)
  end do
  nlines = 1
else ! there could be multiple variable length strings to be read...
! Get dataspace and allocate memory for read buffer.
  call H5Dget_space_f(self%head%next%objectID, space, hdferr)
  call error_check_(self, 'readDatasetStringArray_:H5Dget_space_f:'//trim(dataname), hdferr)

  ! this routine returns the rank of the data set in the hdferr variable when successful, otherwise -1
  call H5Sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
  if (hdferr.lt.0) then 
    call error_check_(self, 'readDatasetStringArray_:H5Sget_simple_extent_dims_f:'//trim(dataname), hdferr)
  end if

  ALLOCATE(rdata(1:dims(1)), stringarray(1:dims(1)))
!
! Read the data.
!
  f_ptr = C_LOC(rdata(1))
  call h5dread_f(self%head%next%objectID, H5T_STRING, f_ptr, hdferr)
  ! call h5dread_f(self%head%next%objectID, H5T_NATIVE_CHARACTER, f_ptr, hdferr)
  call error_check_(self, 'readDatasetStringArray_:h5dread_f:'//trim(dataname), hdferr)

!
! convert the data to a string array
!
  DO i = 1, dims(1)
    call C_F_POINTER(rdata(i), pfstr)

    length = 0
    DO
     IF(pfstr(length+1:length+1).EQ.C_NULL_CHAR.OR.length.GE.fnlen) EXIT
     length = length + 1
    ENDDO
    stringarray(i) = pfstr(1:length)
  END DO

  nlines = dims(1)
  DEALLOCATE(rdata)
end if

call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetStringArray_:h5sclose_f:'//trim(dataname), hdferr)

call H5Tclose_f(filetype, hdferr)
call error_check_(self, 'readDatasetStringArray_:H5Tclose_f:'//trim(dataname), hdferr)

! close the dataset
call self%pop_()

end subroutine readDatasetStringArray_

!--------------------------------------------------------------------------
recursive function extractDatasetTextfile_(self, dataname, textfile) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! read a string array from a data set and stores it as a text file

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                              :: self
character(fnlen),INTENT(IN)                             :: dataname
character(fnlen),INTENT(IN)                             :: textfile

integer(kind=irg)                                       :: success

integer(HID_T)                                          :: filetype, space ! Handles
integer                                                 :: hdferr, i, length, dt
integer(HSIZE_T), DIMENSION(1:1)                        :: dims
integer(HSIZE_T), DIMENSION(1:1)                        :: maxdims

character(len = fnlen, kind=c_char),  POINTER           :: pfstr ! A pointer to a Fortran string
TYPE(C_PTR), DIMENSION(:), ALLOCATABLE, TARGET          :: rdata ! Read buffer
TYPE(C_PTR)                                             :: f_ptr

! Open dataset.
!
hdferr = self%openDataset_(dataname)
!
! Get the datatype.
!
call H5Dget_type_f(self%head%next%objectID, filetype, hdferr)
call error_check_(self, 'extractDatasetTextfile_:H5Dget_type_f:'//trim(dataname), hdferr)

! Get dataspace and allocate memory for read buffer.
!
call H5Dget_space_f(self%head%next%objectID, space, hdferr)
call error_check_(self, 'extractDatasetTextfile_:H5Dget_space_f:'//trim(dataname), hdferr)

call H5Sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
if (hdferr.lt.0) then 
  call error_check_(self, 'extractDatasetTextfile_:H5Sget_simple_extent_dims_f:'//trim(dataname), hdferr)
end if

ALLOCATE(rdata(1:dims(1)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1))
call h5dread_f(self%head%next%objectID, H5T_STRING, f_ptr, hdferr)
call error_check_(self, 'extractDatasetTextfile_:h5dread_f:'//trim(dataname), hdferr)

!
! store the datain a textfile
!
dt = 55
open(unit=dt, file=trim(textfile), status='unknown', form='formatted')
DO i = 1, dims(1)
  call C_F_POINTER(rdata(i), pfstr)

  length = 0
  DO
     IF(pfstr(length+1:length+1).EQ.C_NULL_CHAR.OR.length.GE.fnlen) EXIT
     length = length + 1
  ENDDO
  write(dt,"(A)") pfstr(1:length)
END DO
close(unit=dt,status='keep')

DEALLOCATE(rdata)

call h5sclose_f(space, hdferr)
call error_check_(self, 'extractDatasetTextfile_:h5sclose_f:'//trim(dataname), hdferr)

call H5Tclose_f(filetype, hdferr)
call error_check_(self, 'extractDatasetTextfile_:H5Tclose_f:'//trim(dataname), hdferr)

! close the dataset
call pop_(self)

success = 0

end function extractDatasetTextfile_

!--------------------------------------------------------------------------
recursive function writeDatasetStringArray_(self, dataname, inputarray, nlines, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/08/20
  !!
  !! write a string array

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: nlines
character(len=fnlen),INTENT(IN)            :: inputarray(nlines) 

logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success

integer(kind=irg),parameter                :: fnlenp = fnlen+1
character(len=fnlenp, KIND=c_char),TARGET  :: stringarray(nlines) 
integer(HSIZE_T)                           :: dim0 
integer(SIZE_T)                            :: sdim 
integer(HID_T)                             :: filetype, space, dset ! Handles
integer                                    :: hdferr, i, rnk, l
integer(HSIZE_T), DIMENSION(1:1)           :: dims
! integer(HSIZE_T), DIMENSION(1:2)           :: maxdims

TYPE(C_PTR), ALLOCATABLE, TARGET           :: wdata(:)
TYPE(C_PTR)                                :: f_ptr

success = 0

stringarray = ''
! first, convert the stringarray to an array of C-pointers, with each string
! terminated by a C_NULL_CHAR.
dims(1) = nlines
allocate(wdata(1:dims(1)))
do i=1,dims(1)
  l = len(trim(inputarray(i)))+1
  if (l.gt.fnlenp) l = fnlenp           ! if a string is too long for some reason, we just truncate it
  stringarray(i) = trim(inputarray(i))
  stringarray(i)(l:l) = C_NULL_CHAR
  wdata(i) = C_LOC(stringarray(i)(1:1))
end do

! then we write this C_ptr to the HDF file in the proper data set

! first create the memory data type (filetype)
call H5Tcopy_f(H5T_STRING, filetype, hdferr)
call error_check_(self, 'writeDatasetStringArray_:H5Tcopy_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
if (hdferr.lt.0) then
  call error_check_(self, 'writeDatasetStringArray_:h5screate_simple_f:'//trim(dataname), hdferr)
end if

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then 
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetStringArray_:h5dopen_f:'//trim(dataname)//':overwrite', hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), filetype, space, dset, hdferr)
  call error_check_(self, 'writeDatasetStringArray_:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

f_ptr = C_LOC(wdata(1))
call h5dwrite_f(dset, filetype, f_ptr, hdferr )
call error_check_(self, 'writeDatasetStringArray_:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetStringArray_:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetStringArray_:h5sclose_f:'//trim(dataname), hdferr)

call H5Tclose_f(filetype, hdferr)
call error_check_(self, 'writeDatasetStringArray_:H5Tclose_f:'//trim(dataname), hdferr)

deallocate(wdata)

! that's it

end function writeDatasetStringArray_

!--------------------------------------------------------------------------
recursive function writeDatasetCharArray1D(self, dataname, chararray, dims, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 1D array of c_chars    [variable type H5T_STD_U8LE] 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T), INTENT(IN)               :: dims(1)
character(len=1),TARGET                    :: chararray(dims(1)) 
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, i, rnk, l

TYPE(C_PTR), dimension(1:1), TARGET        :: wdata
TYPE(C_PTR)                                :: f_ptr

success = 0

wdata(1) = C_LOC(chararray(1))

! then we write this C_ptr to the HDF file in the proper data set

! Create dataspace.
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetCharArray1D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the c_char data to it.
!
if (present(overwrite)) then 
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetCharArray1D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_U8LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetCharArray1D:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_STD_U8LE, wdata(1), hdferr )
call error_check_(self, 'writeDatasetCharArray1D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetCharArray1D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetCharArray1D:h5sclose_f:'//trim(dataname), hdferr)

! that's it

end function writeDatasetCharArray1D

!--------------------------------------------------------------------------
recursive function writeDatasetCharArray2D(self, dataname, chararray, dims, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 2D array of c_chars    [variable type H5T_STD_U8LE] 


use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T), INTENT(IN)               :: dims(2)
character(len=1),TARGET                    :: chararray(dims(1), dims(2)) 
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, i, rnk, l

TYPE(C_PTR), dimension(1:1), TARGET        :: wdata

success = 0

wdata(1) = C_LOC(chararray(1,1))
! then we write this C_ptr to the HDF file in the proper data set

! Create dataspace.
rnk = 2
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetCharArray2D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the c_char data to it.
!
if (present(overwrite)) then 
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetCharArray2D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_U8LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetCharArray2D:h5dcreate_f:'//trim(dataname), hdferr)
end if 

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_STD_U8LE, wdata(1), hdferr )
call error_check_(self, 'writeDatasetCharArray2D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetCharArray2D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetCharArray2D:h5sclose_f:'//trim(dataname), hdferr)

! that's it

end function writeDatasetCharArray2D

!--------------------------------------------------------------------------
recursive function writeDatasetCharArray3D(self, dataname, chararray, dims, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 3D array of c_chars    [variable type H5T_STD_U8LE] 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T), INTENT(IN)               :: dims(3)
character(len=1),TARGET                    :: chararray(dims(1), dims(2), dims(3)) 
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, i, rnk, l

TYPE(C_PTR), dimension(1:3), TARGET        :: wdata

success = 0

wdata(1) = C_LOC(chararray(1,1,1))
! then we write this C_ptr to the HDF file in the proper data set

! Create dataspace.
rnk = 3
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetCharArray3D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the c_char data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetCharArray3D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_U8LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetCharArray3D:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_STD_U8LE, wdata(1), hdferr )
call error_check_(self, 'writeDatasetCharArray3D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetCharArray3D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetCharArray3D:h5sclose_f:'//trim(dataname), hdferr)

! that's it

end function writeDatasetCharArray3D

!--------------------------------------------------------------------------
recursive function writeDatasetCharArray4D(self, dataname, chararray, dims, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 4D array of c_chars    [variable type H5T_STD_U8LE] 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T), INTENT(IN)               :: dims(4)
character(len=1),TARGET                    :: chararray(dims(1), dims(2), dims(3), dims(4)) 
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, i, rnk, l

TYPE(C_PTR), dimension(1:4), TARGET        :: wdata

success = 0

wdata(1) = C_LOC(chararray(1,1,1,1))

! then we write this C_ptr to the HDF file in the proper data set

! Create dataspace.
rnk = 4
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetCharArray4D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the c_char data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetCharArray4D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_U8LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetCharArray4D:h5dcreate_f:'//trim(dataname), hdferr)
end if 

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_STD_U8LE, wdata(1), hdferr )
call error_check_(self, 'writeDatasetCharArray4D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetCharArray4D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetCharArray4D:h5sclose_f:'//trim(dataname), hdferr)

! that's it

end function writeDatasetCharArray4D

!--------------------------------------------------------------------------
recursive function writeDatasetInteger_(self, dataname, intval, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write an integer data set to the current file or group ID

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: intval

logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:1)           :: dims

integer, dimension(1:1), TARGET            :: wdata
TYPE(C_PTR)                                :: f_ptr

success = 0

dims(1) = 1
wdata(1) = intval

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1))

! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetInteger_:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetInteger_:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_I32LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetInteger_:h5dcreate_f:'//trim(dataname), hdferr)
end if 

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdferr )
call error_check_(self, 'writeDatasetInteger_:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetInteger_:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetInteger_:h5sclose_f:'//trim(dataname), hdferr)

! that's it

end function writeDatasetInteger_

!--------------------------------------------------------------------------
recursive function writeDatasetInteger1byteArray1D_(self, dataname, intarr, dim0, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 1D integer array data set to the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: dim0
integer(kind=1),INTENT(IN)                 :: intarr(dim0)
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:1)           :: dims

integer, dimension(1:dim0), TARGET         :: wdata
TYPE(C_PTR)                                :: f_ptr

success = 0

dims(1) = dim0
wdata = intarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1))

! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetInteger1byteArray1D_:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetInteger1byteArray1D_:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_I8LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetInteger1byteArray1D_:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdferr )
call error_check_(self, 'writeDatasetInteger1byteArray1D_:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetInteger1byteArray1D_:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetInteger1byteArray1D_:h5sclose_f:'//trim(dataname), hdferr)

! that's it

end function writeDatasetInteger1byteArray1D_

!--------------------------------------------------------------------------
recursive function writeDatasetIntegerArray1D(self, dataname, intarr, dim0, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 1D integer array data set to the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: dim0
integer(kind=irg),INTENT(IN)               :: intarr(dim0)
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success, istat

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:1)           :: dims

integer(kind=4),allocatable,TARGET         :: wdata(:)

TYPE(C_PTR)                                :: f_ptr

success = 0

allocate(wdata(dim0), stat=istat)

wdata = intarr
dims(1) = dim0

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1))

! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetIntegerArray1D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetIntegerArray1D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_I32LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetIntegerArray1D:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdferr )
call error_check_(self, 'writeDatasetIntegerArray1D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetIntegerArray1D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetIntegerArray1D:h5sclose_f:'//trim(dataname), hdferr)
DEALLOCATE(wdata, stat=istat)

! that's it

end function writeDatasetIntegerArray1D

!--------------------------------------------------------------------------
recursive function writeDatasetIntegerArray2D(self, dataname, intarr, dim0, dim1, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 2D integer array data set to the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: dim0
integer(kind=irg),INTENT(IN)               :: dim1
integer(kind=irg),INTENT(IN)               :: intarr(dim0, dim1)
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success, istat

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)           :: dims

integer(kind=4),allocatable,TARGET         :: wdata(:,:)

TYPE(C_PTR)                                :: f_ptr

success = 0

allocate(wdata(dim0, dim1), stat=istat)

wdata = intarr

dims(1:2) = (/ dim0, dim1 /)


! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1))

! Create dataspace.
!
rnk = 2
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetIntegerArray2D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetIntegerArray2D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_I32LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetIntegerArray2D:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdferr )
call error_check_(self, 'writeDatasetIntegerArray2D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetIntegerArray2D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetIntegerArray2D:h5sclose_f:'//trim(dataname), hdferr)
DEALLOCATE(wdata, stat=istat)

! that's it

end function writeDatasetIntegerArray2D

!--------------------------------------------------------------------------
recursive function writeDatasetIntegerArray3D(self, dataname, intarr, dim0, dim1, dim2, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 3D integer array data set to the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: dim0
integer(kind=irg),INTENT(IN)               :: dim1
integer(kind=irg),INTENT(IN)               :: dim2
integer(kind=irg),INTENT(IN),TARGET        :: intarr(dim0, dim1, dim2)
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success, istat

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:3)           :: dims

integer(kind=4),allocatable,TARGET         :: wdata(:,:,:)
TYPE(C_PTR)                                :: f_ptr

success = 0

allocate(wdata(dim0, dim1, dim2), stat=istat)

dims(1:3) = (/ dim0, dim1, dim2 /)
wdata = intarr

! get a C pointer to the integer array
f_ptr = C_LOC(intarr(1, 1, 1))

! Create dataspace.
!
rnk = 3
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetIntegerArray3D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetIntegerArray3D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_I32LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetIntegerArray3D:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdferr )
call error_check_(self, 'writeDatasetIntegerArray3D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetIntegerArray3D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetIntegerArray3D:h5sclose_f:'//trim(dataname), hdferr)

DEALLOCATE(wdata, stat=istat)

! that's it

end function writeDatasetIntegerArray3D

!--------------------------------------------------------------------------
recursive function writeDatasetIntegerArray4D(self, dataname, intarr, dim0, dim1, dim2, dim3, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 4D integer array data set to the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: dim0
integer(kind=irg),INTENT(IN)               :: dim1
integer(kind=irg),INTENT(IN)               :: dim2
integer(kind=irg),INTENT(IN)               :: dim3
integer(kind=irg),INTENT(IN)               :: intarr(dim0, dim1, dim2, dim3)
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success, istat

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:4)           :: dims

integer(kind=4),allocatable,TARGET         :: wdata(:,:,:,:)
TYPE(C_PTR)                                :: f_ptr

success = 0

allocate(wdata(dim0, dim1, dim2, dim3), stat=istat)

wdata = intarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1,1,1))

! Create dataspace.
!
rnk = 4
dims(1:4) = (/ dim0, dim1, dim2, dim3 /)

call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetIntegerArray4D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetIntegerArray4D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_I32LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetIntegerArray4D:h5dcreate_f:'//trim(dataname), hdferr)
end if 

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdferr )
call error_check_(self, 'writeDatasetIntegerArray4D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetIntegerArray4D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetIntegerArray4D:h5sclose_f:'//trim(dataname), hdferr)

DEALLOCATE(wdata, stat=istat)

! that's it

end function writeDatasetIntegerArray4D

!--------------------------------------------------------------------------
recursive function writeDatasetFloat_(self, dataname, fltval, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a single precision float data set to the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
real(kind=sgl),INTENT(IN)                  :: fltval

logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)
integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:1)           :: dims

real(real_kind4), dimension(1:1), TARGET   :: wdata
TYPE(C_PTR)                                :: f_ptr

success = 0

dims(1) = 1
wdata(1) = fltval

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1))

! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetFloat_:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetFloat_:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetFloat_:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_REAL, f_ptr, hdferr )
call error_check_(self, 'writeDatasetFloat_:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetFloat_:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetFloat_:h5sclose_f:'//trim(dataname), hdferr)

! that's it

end function writeDatasetFloat_

!--------------------------------------------------------------------------
recursive function writeDatasetDouble_(self, dataname, dblval, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a double precision float data set to the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
real(kind=dbl),INTENT(IN)                  :: dblval

logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)
integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:1)           :: dims

real(real_kind8), dimension(1:1), TARGET   :: wdata
TYPE(C_PTR)                                :: f_ptr

success = 0

dims(1) = 1
wdata(1) = dblval

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1))

! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetDouble_:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetDouble_:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetDouble_:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
call error_check_(self, 'writeDatasetDouble_:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetDouble_:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetDouble_:h5sclose_f:'//trim(dataname), hdferr)

! that's it

end function writeDatasetDouble_

!--------------------------------------------------------------------------
recursive function writeDatasetFloatArray1D(self, dataname, fltarr, dim0, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 1D single precision float array data set to the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: dim0
real(kind=sgl),INTENT(IN)                  :: fltarr(dim0)
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success, istat

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)
integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:1)           :: dims

real(real_kind4),allocatable,TARGET        :: wdata(:)

TYPE(C_PTR)                                :: f_ptr

success = 0

dims(1) = dim0
allocate(wdata(dim0), stat=istat)

wdata = fltarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1))

! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetFloatArray1D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetFloatArray1D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetFloatArray1D:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_REAL, f_ptr, hdferr )
call error_check_(self, 'writeDatasetFloatArray1D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetFloatArray1D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetFloatArray1D:h5sclose_f:'//trim(dataname), hdferr)
DEALLOCATE(wdata, stat=istat)

! that's it

end function writeDatasetFloatArray1D

!--------------------------------------------------------------------------
recursive function writeDatasetFloatArray2D(self, dataname, fltarr, dim0, dim1, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 2D single precision float array data set to the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: dim0
integer(kind=irg),INTENT(IN)               :: dim1
real(kind=sgl),INTENT(IN)                  :: fltarr(dim0, dim1)
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success, istat

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)
integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)           :: dims

real(real_kind4),allocatable,TARGET        :: wdata(:,:)

TYPE(C_PTR)                                :: f_ptr

success = 0
allocate(wdata(dim0, dim1), stat=istat)

wdata = fltarr


dims(1:2) = (/ dim0, dim1 /)

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1))

! Create dataspace.
!
rnk = 2
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetFloatArray2D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetFloatArray2D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetFloatArray2D:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_REAL, f_ptr, hdferr )
call error_check_(self, 'writeDatasetFloatArray2D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetFloatArray2D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetFloatArray2D:h5sclose_f:'//trim(dataname), hdferr)
DEALLOCATE(wdata, stat=istat)

! that's it

end function writeDatasetFloatArray2D

!--------------------------------------------------------------------------
recursive function writeDatasetFloatArray3D(self, dataname, fltarr, dim0, dim1, dim2, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 3D single precision float array data set to the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: dim0
integer(kind=irg),INTENT(IN)               :: dim1
integer(kind=irg),INTENT(IN)               :: dim2
real(kind=sgl),INTENT(IN)                  :: fltarr(dim0, dim1, dim2)
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success, istat

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)
integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:3)           :: dims

real(real_kind4),allocatable,TARGET        :: wdata(:,:,:)

TYPE(C_PTR)                                :: f_ptr

success = 0

allocate(wdata(dim0, dim1, dim2), stat=istat)

wdata = fltarr

dims(1:3) = (/ dim0, dim1, dim2 /)

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1,1))

! Create dataspace.
!
rnk = 3
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetFloatArray3D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetFloatArray3D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetFloatArray3D:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_REAL, f_ptr, hdferr )
call error_check_(self, 'writeDatasetFloatArray3D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetFloatArray3D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetFloatArray3D:h5sclose_f:'//trim(dataname), hdferr)
DEALLOCATE(wdata, stat=istat)

! that's it

end function writeDatasetFloatArray3D

!--------------------------------------------------------------------------
recursive function writeDatasetFloatArray4D(self, dataname, fltarr, dim0, dim1, dim2, dim3, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 4D single precision float array data set to the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: dim0
integer(kind=irg),INTENT(IN)               :: dim1
integer(kind=irg),INTENT(IN)               :: dim2
integer(kind=irg),INTENT(IN)               :: dim3
real(kind=sgl),INTENT(IN)                  :: fltarr(dim0, dim1, dim2, dim3)
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success, istat

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)
integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:4)           :: dims

real(real_kind4),allocatable,TARGET        :: wdata(:,:,:,:)

TYPE(C_PTR)                                :: f_ptr

success = 0

allocate(wdata(dim0, dim1, dim2, dim3), stat=istat)

wdata = fltarr

dims(1:4) = (/ dim0, dim1, dim2, dim3 /)

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1,1,1))

! Create dataspace.
!
rnk = 4
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetFloatArray4D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetFloatArray4D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetFloatArray4D:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_REAL, f_ptr, hdferr )
call error_check_(self, 'writeDatasetFloatArray4D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetFloatArray4D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetFloatArray4D:h5sclose_f:'//trim(dataname), hdferr)
DEALLOCATE(wdata, stat=istat)

! that's it

end function writeDatasetFloatArray4D

!--------------------------------------------------------------------------
recursive function writeDatasetFloatArray6D(self, dataname, fltarr, dim0, dim1, dim2, dim3, dim4, dim5, &
                                            overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 6D single precision float array data set to the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: dim0
integer(kind=irg),INTENT(IN)               :: dim1
integer(kind=irg),INTENT(IN)               :: dim2
integer(kind=irg),INTENT(IN)               :: dim3
integer(kind=irg),INTENT(IN)               :: dim4
integer(kind=irg),INTENT(IN)               :: dim5
real(kind=sgl),INTENT(IN)                  :: fltarr(dim0, dim1, dim2, dim3, dim4, dim5)
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success, istat

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)
integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:6)           :: dims

TYPE(C_PTR)                                :: f_ptr
real(real_kind4),allocatable,TARGET        :: wdata(:,:,:,:,:,:)

success = 0

allocate(wdata(dim0,dim1,dim2,dim3,dim4,dim5), stat=istat)
dims(1:6) = (/ dim0, dim1, dim2, dim3, dim4, dim5 /)
wdata = fltarr

 ! get a C pointer to the float array
f_ptr = C_LOC(wdata(1,1,1,1,1,1))

! Create dataspace.
!
rnk = 6
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetFloatArray6D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the 6D float data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetFloatArray6D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetFloatArray6D:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_REAL, f_ptr, hdferr )
call error_check_(self, 'writeDatasetFloatArray6D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetFloatArray6D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetFloatArray6D:h5sclose_f:'//trim(dataname), hdferr)

! that's it

end function writeDatasetFloatArray6D

!--------------------------------------------------------------------------
recursive function writeDatasetDoubleArray1D(self, dataname, dblarr, dim0, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 1D double precision float array data set to the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: dim0
real(kind=dbl),INTENT(IN)                  :: dblarr(dim0)
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success, istat

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)
integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:1)           :: dims

real(real_kind8),allocatable,TARGET        :: wdata(:)

TYPE(C_PTR)                                :: f_ptr

success = 0

allocate(wdata(dim0), stat=istat)
dims(1) = dim0
wdata = dblarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1))

! Create dataspace.
!
rnk = 1
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetDoubleArray1D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetDoubleArray1D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetDoubleArray1D:h5dcreate_f:'//trim(dataname), hdferr)
end if 

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
call error_check_(self, 'writeDatasetDoubleArray1D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetDoubleArray1D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetDoubleArray1D:h5sclose_f:'//trim(dataname), hdferr)

! that's it

end function writeDatasetDoubleArray1D

!--------------------------------------------------------------------------
recursive function writeDatasetDoubleArray2D(self, dataname, dblarr, dim0, dim1, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 2D double precision float array data set to the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: dim0
integer(kind=irg),INTENT(IN)               :: dim1
real(kind=dbl),INTENT(IN)                  :: dblarr(dim0, dim1)
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success, istat

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)
integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:2)           :: dims

real(real_kind8),allocatable,TARGET        :: wdata(:,:)
TYPE(C_PTR)                                :: f_ptr

success = 0

allocate(wdata(dim0,dim1), stat=istat)
dims(1:2) = (/ dim0, dim1 /)
wdata = dblarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1))

! Create dataspace.
!
rnk = 2
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetDoubleArray2D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetDoubleArray2D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetDoubleArray2D:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
call error_check_(self, 'writeDatasetDoubleArray2D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetDoubleArray2D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetDoubleArray2D:h5sclose_f:'//trim(dataname), hdferr)

! that's it

end function writeDatasetDoubleArray2D

!--------------------------------------------------------------------------
recursive function writeDatasetDoubleArray3D(self, dataname, dblarr, dim0, dim1, dim2, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 3D double precision float array data set to the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: dim0
integer(kind=irg),INTENT(IN)               :: dim1
integer(kind=irg),INTENT(IN)               :: dim2
real(kind=dbl),INTENT(IN)                  :: dblarr(dim0, dim1, dim2)
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success, istat

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)
integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:3)           :: dims

real(real_kind8),allocatable,TARGET        :: wdata(:,:,:)

TYPE(C_PTR)                                :: f_ptr

success = 0

allocate(wdata(dim0,dim1,dim2), stat=istat)
dims(1:3) = (/ dim0, dim1, dim2 /)
wdata = dblarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1,1))

! Create dataspace.
!
rnk = 3
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetDoubleArray3D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetDoubleArray3D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetDoubleArray3D:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
call error_check_(self, 'writeDatasetDoubleArray3D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetDoubleArray3D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetDoubleArray3D:h5sclose_f:'//trim(dataname), hdferr)

! that's it

end function writeDatasetDoubleArray3D

!--------------------------------------------------------------------------
recursive function writeDatasetDoubleArray4D(self, dataname, dblarr, dim0, dim1, dim2, dim3, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! write a 4D double precision float array data set to the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg),INTENT(IN)               :: dim0
integer(kind=irg),INTENT(IN)               :: dim1
integer(kind=irg),INTENT(IN)               :: dim2
integer(kind=irg),INTENT(IN)               :: dim3
real(kind=dbl),INTENT(IN)                  :: dblarr(dim0, dim1, dim2, dim3)
logical,INTENT(IN),OPTIONAL                :: overwrite
integer(kind=irg)                          :: success, istat

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)
integer(HID_T)                             :: space, dset ! Handles
integer                                    :: hdferr, rnk
integer(HSIZE_T), DIMENSION(1:4)           :: dims

real(real_kind8),allocatable,TARGET        :: wdata(:,:,:,:)
TYPE(C_PTR)                                :: f_ptr

success = 0

allocate(wdata(dim0,dim1,dim2,dim3), stat=istat)
dims(1:4) = (/ dim0, dim1, dim2, dim3 /)
wdata = dblarr

! get a C pointer to the integer array
f_ptr = C_LOC(wdata(1,1,1,1))

! Create dataspace.
!
rnk = 4
call h5screate_simple_f(rnk, dims, space, hdferr)
call error_check_(self, 'writeDatasetDoubleArray4D:h5screate_simple_f:'//trim(dataname), hdferr)

!
! Create the dataset and write the variable-length string data to it.
!
if (present(overwrite)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeDatasetDoubleArray4D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
  call error_check_(self, 'writeDatasetDoubleArray4D:h5dcreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr )
call error_check_(self, 'writeDatasetDoubleArray4D:h5dwrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5dclose_f(dset , hdferr)
call error_check_(self, 'writeDatasetDoubleArray4D:h5dclose_f:'//trim(dataname), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'writeDatasetDoubleArray4D:h5sclose_f:'//trim(dataname), hdferr)

! that's it

end function writeDatasetDoubleArray4D


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------READ subroutines below--------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine readDatasetCharArray1D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 1D char array data set from the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(1)
integer(kind=irg), INTENT(OUT)             :: hdferr
character(len=1), dimension(:), allocatable, TARGET, INTENT(OUT)     :: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:1)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetCharArray1D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetCharArray1D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetCharArray1D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1))
call h5dread_f( dset, H5T_STD_U8LE, f_ptr, hdferr)
call error_check_(self, 'readDatasetCharArray1D:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetCharArray1D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetCharArray1D:h5dclose_f:'//trim(dataname), hdferr)


! that's it

end subroutine readDatasetCharArray1D

!--------------------------------------------------------------------------
recursive subroutine readDatasetCharArray2D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 2D char array data set from the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(2)
integer(kind=irg), INTENT(OUT)             :: hdferr
character(len=1), dimension(:,:), allocatable, TARGET, INTENT(OUT)   :: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:2)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetCharArray2D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetCharArray2D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetCharArray2D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1),1:dims(2)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1))
call h5dread_f( dset, H5T_STD_U8LE, f_ptr, hdferr)
call error_check_(self, 'readDatasetCharArray2D:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetCharArray2D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetCharArray2D:h5dclose_f:'//trim(dataname), hdferr)


! that's it

end subroutine readDatasetCharArray2D

!--------------------------------------------------------------------------
recursive subroutine readDatasetCharArray3D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 3D char array data set from the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(3)
integer(kind=irg), INTENT(OUT)             :: hdferr
character(len=1), dimension(:,:,:), allocatable, TARGET, INTENT(OUT) :: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:3)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetCharArray3D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetCharArray3D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetCharArray3D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1),1:dims(2),1:dims(3)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1))
call h5dread_f( dset, H5T_STD_U8LE, f_ptr, hdferr)
call error_check_(self, 'readDatasetCharArray3D:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetCharArray3D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetCharArray3D:h5dclose_f:'//trim(dataname), hdferr)


! that's it

end subroutine readDatasetCharArray3D

!--------------------------------------------------------------------------
recursive subroutine readDatasetCharArray4D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 4D char array data set from the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(4)
integer(kind=irg), INTENT(OUT)             :: hdferr
character(len=1), dimension(:,:,:,:), allocatable, TARGET, INTENT(OUT) :: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:4)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetCharArray4D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetCharArray4D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetCharArray4D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1),1:dims(2),1:dims(3),1:dims(4)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1,1))
call h5dread_f( dset, H5T_STD_U8LE, f_ptr, hdferr)
call error_check_(self, 'readDatasetCharArray4D:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetCharArray4D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetCharArray4D:h5dclose_f:'//trim(dataname), hdferr)


! that's it

end subroutine readDatasetCharArray4D

!--------------------------------------------------------------------------
recursive subroutine readDatasetInteger_(self, dataname, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns an integer data set from the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg), INTENT(OUT)             :: hdferr
integer,  TARGET, INTENT(OUT)              :: rdata

integer(HID_T)                             :: space, dset ! Handles

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetInteger_:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetInteger_:h5dget_space_f:'//trim(dataname), hdferr)

!
! Read the data.
!
f_ptr = C_LOC(rdata)
call h5dread_f( dset, H5T_NATIVE_INTEGER, f_ptr, hdferr)
call error_check_(self, 'readDatasetInteger_:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetInteger_:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetInteger_:h5dclose_f:'//trim(dataname), hdferr)

! that's it

end subroutine readDatasetInteger_

!--------------------------------------------------------------------------
recursive subroutine readDatasetIntegerArray1D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 1D integer array data set from the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(1)
integer(kind=irg), INTENT(OUT)             :: hdferr
integer, dimension(:), allocatable, TARGET, INTENT(OUT)              :: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:1)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetIntegerArray1D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetIntegerArray1D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetIntegerArray1D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1))
call h5dread_f( dset, H5T_NATIVE_INTEGER, f_ptr, hdferr)
call error_check_(self, 'readDatasetIntegerArray1D:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetIntegerArray1D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetIntegerArray1D:h5dclose_f:'//trim(dataname), hdferr)


! that's it

end subroutine readDatasetIntegerArray1D

!--------------------------------------------------------------------------
recursive subroutine readDatasetIntegerArray2D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 2D integer array data set from the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(2)
integer(kind=irg), INTENT(OUT)             :: hdferr
integer, dimension(:,:), allocatable, TARGET, INTENT(OUT)            :: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:2)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetIntegerArray2D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetIntegerArray2D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetIntegerArray2D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1),1:dims(2)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1))
call h5dread_f( dset, H5T_NATIVE_INTEGER, f_ptr, hdferr)
call error_check_(self, 'readDatasetIntegerArray2D:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetIntegerArray2D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetIntegerArray2D:h5dclose_f:'//trim(dataname), hdferr)


! that's it

end subroutine readDatasetIntegerArray2D

!--------------------------------------------------------------------------
recursive subroutine readDatasetIntegerArray3D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 3D integer array data set from the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(3)
integer(kind=irg), INTENT(OUT)             :: hdferr
integer, dimension(:,:,:), allocatable, TARGET, INTENT(OUT)          :: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:3)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetIntegerArray3D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetIntegerArray3D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetIntegerArray3D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1),1:dims(2),1:dims(3)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1))
call h5dread_f( dset, H5T_NATIVE_INTEGER, f_ptr, hdferr)
call error_check_(self, 'readDatasetIntegerArray3D:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetIntegerArray3D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetIntegerArray3D:h5dclose_f:'//trim(dataname), hdferr)


! that's it

end subroutine readDatasetIntegerArray3D

!--------------------------------------------------------------------------
recursive subroutine readDatasetIntegerArray4D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 4D integer array data set from the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(4)
integer(kind=irg), INTENT(OUT)             :: hdferr
!integer, dimension(:,:,:,:), allocatable, TARGET, INTENT(OUT)        :: rdata
integer, dimension(:,:,:,:), allocatable, TARGET        :: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:4)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetIntegerArray4D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetIntegerArray4D:h5dget_space_f:'//trim(dataname), hdferr)


call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetIntegerArray4D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1),1:dims(2),1:dims(3),1:dims(4)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1,1))
call h5dread_f( dset, H5T_NATIVE_INTEGER, f_ptr, hdferr)
call error_check_(self, 'readDatasetIntegerArray4D:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetIntegerArray4D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetIntegerArray4D:h5dclose_f:'//trim(dataname), hdferr)


! that's it

end subroutine readDatasetIntegerArray4D

!--------------------------------------------------------------------------
recursive subroutine readDatasetFloat_(self, dataname, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a float data set from the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg), INTENT(OUT)             :: hdferr
real(real_kind4), TARGET, INTENT(OUT)      :: rdata

integer(HID_T)                             :: space, dset ! Handles

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetFloat_:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetFloat_:h5dget_space_f:'//trim(dataname), hdferr)


! Read the data.
!
f_ptr = C_LOC(rdata)
call h5dread_f( dset, H5T_NATIVE_REAL, f_ptr, hdferr)
call error_check_(self, 'readDatasetFloat_:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetFloat_:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetFloat_:h5dclose_f:'//trim(dataname), hdferr)

! that's it

end subroutine readDatasetFloat_

!--------------------------------------------------------------------------
recursive subroutine readDatasetFloatArray1D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 1D float array data set from the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(1)
integer(kind=irg), INTENT(OUT)             :: hdferr
real(real_kind4), dimension(:), allocatable, TARGET, INTENT(OUT)      :: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:1)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetFloatArray1D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetFloatArray1D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetFloatArray1D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1))
call h5dread_f( dset, H5T_NATIVE_REAL, f_ptr, hdferr)
call error_check_(self, 'readDatasetFloatArray1D:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetFloatArray1D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetFloatArray1D:h5dclose_f:'//trim(dataname), hdferr)


! that's it

end subroutine readDatasetFloatArray1D

!--------------------------------------------------------------------------
recursive subroutine readDatasetFloatArray2D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 2D float array data set from the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(2)
integer(kind=irg), INTENT(OUT)             :: hdferr
real(real_kind4), dimension(:,:), allocatable, TARGET, INTENT(OUT)    :: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:2)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetFloatArray2D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetFloatArray2D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetFloatArray2D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1),1:dims(2)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1))
call h5dread_f( dset, H5T_NATIVE_REAL, f_ptr, hdferr)
call error_check_(self, 'readDatasetFloatArray2D:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetFloatArray2D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetFloatArray2D:h5dclose_f:'//trim(dataname), hdferr)


! that's it

end subroutine readDatasetFloatArray2D

!--------------------------------------------------------------------------
recursive subroutine readDatasetFloatArray3D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 3D float array data set from the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(3)
integer(kind=irg), INTENT(OUT)             :: hdferr
real(real_kind4), dimension(:,:,:), allocatable, TARGET, INTENT(OUT)  :: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:3)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetFloatArray3D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetFloatArray3D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetFloatArray3D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1),1:dims(2),1:dims(3)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1))
call h5dread_f( dset, H5T_NATIVE_REAL, f_ptr, hdferr)
call error_check_(self, 'readDatasetFloatArray3D:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetFloatArray3D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetFloatArray3D:h5dclose_f:'//trim(dataname), hdferr)


! that's it

end subroutine readDatasetFloatArray3D

!--------------------------------------------------------------------------
recursive subroutine readDatasetFloatArray4D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 4D float array data set from the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(4)
integer(kind=irg), INTENT(OUT)             :: hdferr
real(real_kind4), dimension(:,:,:,:), allocatable, TARGET, INTENT(OUT):: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:4)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetFloatArray4D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetFloatArray4D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetFloatArray4D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1),1:dims(2),1:dims(3),1:dims(4)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1,1))
call h5dread_f( dset, H5T_NATIVE_REAL, f_ptr, hdferr)
call error_check_(self, 'readDatasetFloatArray4D:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetFloatArray4D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetFloatArray4D:h5dclose_f:'//trim(dataname), hdferr)

! that's it

end subroutine readDatasetFloatArray4D

!--------------------------------------------------------------------------
recursive subroutine readDatasetDouble_(self, dataname, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a double data set from the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(kind=irg), INTENT(OUT)             :: hdferr
real(real_kind8), TARGET, INTENT(OUT)      :: rdata

integer(HID_T)                             :: space, dset ! Handles

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetDouble_:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetDouble_:h5dget_space_f:'//trim(dataname), hdferr)


! Read the data.
!
f_ptr = C_LOC(rdata)
call h5dread_f( dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr)
call error_check_(self, 'readDatasetDouble_:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetDouble_:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetDouble_:h5dclose_f:'//trim(dataname), hdferr)

! that's it

end subroutine readDatasetDouble_

!--------------------------------------------------------------------------
recursive subroutine readDatasetDoubleArray1D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 1D double array data set from the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(1)
integer(kind=irg), INTENT(OUT)             :: hdferr
real(real_kind8), dimension(:), allocatable, TARGET, INTENT(OUT)      :: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:1)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetDoubleArray1D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetDoubleArray1D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetDoubleArray1D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1))
call h5dread_f( dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr)
call error_check_(self, 'readDatasetDoubleArray1D:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetDoubleArray1D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetDoubleArray1D:h5dclose_f:'//trim(dataname), hdferr)

! that's it

end subroutine readDatasetDoubleArray1D

!--------------------------------------------------------------------------
recursive subroutine readDatasetDoubleArray2D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 2D double array data set from the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(2)
integer(kind=irg), INTENT(OUT)             :: hdferr
real(real_kind8), dimension(:,:), allocatable, TARGET, INTENT(OUT)    :: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:2)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetDoubleArray2D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetDoubleArray2D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetDoubleArray2D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1),1:dims(2)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1))
call h5dread_f( dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr)
call error_check_(self, 'readDatasetDoubleArray2D:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetDoubleArray2D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetDoubleArray2D:h5dclose_f:'//trim(dataname), hdferr)


! that's it

end subroutine readDatasetDoubleArray2D

!--------------------------------------------------------------------------
recursive subroutine readDatasetDoubleArray3D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 3D double array data set from the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(3)
integer(kind=irg), INTENT(OUT)             :: hdferr
real(real_kind8), dimension(:,:,:), allocatable, TARGET, INTENT(OUT)  :: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:3)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetDoubleArray3D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetDoubleArray3D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetDoubleArray3D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1),1:dims(2),1:dims(3)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1))
call h5dread_f( dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr)
call error_check_(self, 'readDatasetDoubleArray3D:h5dread_f:'//trim(dataname), hdferr)


!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetDoubleArray3D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetDoubleArray3D:h5dclose_f:'//trim(dataname), hdferr)


! that's it

end subroutine readDatasetDoubleArray3D

!--------------------------------------------------------------------------
recursive subroutine readDatasetDoubleArray4D(self, dataname, dims, hdferr, rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 4D double array data set from the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(OUT)               :: dims(4)
integer(kind=irg), INTENT(OUT)             :: hdferr
real(real_kind8), dimension(:,:,:,:), allocatable, TARGET, INTENT(OUT):: rdata

integer(HID_T)                             :: space, dset ! Handles
integer                                    :: rnk
integer(HSIZE_T), DIMENSION(1:4)           :: maxdims

TYPE(C_PTR)                                :: f_ptr

! open the data set
call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readDatasetDoubleArray4D:h5dopen_f:'//trim(dataname), hdferr)

! get dataspace and allocate memory for read buffer 
call h5dget_space_f(dset,space, hdferr)
call error_check_(self, 'readDatasetDoubleArray4D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, dims, maxdims, hdferr)
call error_check_(self, 'readDatasetDoubleArray4D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


allocate(rdata(1:dims(1),1:dims(2),1:dims(3),1:dims(4)))
!
! Read the data.
!
f_ptr = C_LOC(rdata(1,1,1,1))
call h5dread_f( dset, H5T_NATIVE_DOUBLE, f_ptr, hdferr)
call error_check_(self, 'readDatasetDoubleArray4D:h5dread_f:'//trim(dataname), hdferr)

!
! Close and release resources.
!
call h5sclose_f(space, hdferr)
call error_check_(self, 'readDatasetDoubleArray4D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset , hdferr)
call error_check_(self, 'readDatasetDoubleArray4D:h5dclose_f:'//trim(dataname), hdferr)


! that's it

end subroutine readDatasetDoubleArray4D

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! hyperslab read and write routines
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive function writeHyperslabCharArray2D(self, dataname, wdata, hdims, offset, &
                                       dims, insert) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! writes a 2D hyperslab to the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: hdims(2)
integer(HSIZE_T),INTENT(IN)                :: offset(2)
integer(HSIZE_T),INTENT(IN)                :: dims(2)
character(kind=c_char),INTENT(IN),TARGET   :: wdata(dims(1),dims(2))
logical, OPTIONAL, INTENT(IN)              :: insert
integer(kind=irg)                          :: success

integer(HID_T)                             :: memspace, space, dset ! Handles
integer                                    :: hdferr, rnk

TYPE(C_PTR)                                :: f_ptr

success = 0

rnk = 2
call h5screate_simple_f(rnk, hdims, space, hdferr)
call error_check_(self, 'writeHyperslabCharArray2D:h5screate_simple_f:'//trim(dataname), hdferr)

if (present(insert)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeHyperslabCharArray2D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_U8LE, space, dset, hdferr)
  call error_check_(self, 'writeHyperslabCharArray2D:h5dcreate_f:'//trim(dataname), hdferr)
end if

f_ptr = C_LOC(wdata(1,1))

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call error_check_(self, 'writeHyperslabCharArray2D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'writeHyperslabCharArray2D:h5screate_simple_f:'//trim(dataname), hdferr)

call h5dwrite_f(dset, H5T_STD_U8LE, f_ptr, hdferr, memspace, space)
call error_check_(self, 'writeHyperslabCharArray2D:h5dwrite_f:'//trim(dataname), hdferr)


call h5dclose_f(dset, hdferr)
call error_check_(self, 'writeHyperslabCharArray2D:h5dclose_f:'//trim(dataname), hdferr)


end function writeHyperslabCharArray2D

!--------------------------------------------------------------------------
recursive function writeHyperslabCharArray3D(self, dataname, wdata, hdims, offset, &
                                                 dims, insert) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! writes a 3D hyperslab to the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: hdims(3)
integer(HSIZE_T),INTENT(IN)                :: offset(3)
integer(HSIZE_T),INTENT(IN)                :: dims(3)
character(kind=c_char),INTENT(IN)          :: wdata(dims(1),dims(2),dims(3))
logical, OPTIONAL, INTENT(IN)              :: insert
integer(kind=irg)                          :: success

integer(HID_T)                             :: memspace, space, dset ! Handles
integer                                    :: hdferr, rnk

success = 0

rnk = 3
call h5screate_simple_f(rnk, hdims, space, hdferr)
call error_check_(self, 'writeHyperslabCharArray3D:h5screate_simple_f:'//trim(dataname), hdferr)

if (present(insert)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeHyperslabCharArray3D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_U8LE, space, dset, hdferr)
  call error_check_(self, 'writeHyperslabCharArray3D:h5dcreate_f:'//trim(dataname), hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call error_check_(self, 'writeHyperslabCharArray3D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'writeHyperslabCharArray3D:h5screate_simple_f:'//trim(dataname), hdferr)

call h5dwrite_f(dset, H5T_STD_U8LE, wdata, dims, hdferr, memspace, space)
call error_check_(self, 'writeHyperslabCharArray3D:h5dwrite_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'writeHyperslabCharArray3D:h5dclose_f:'//trim(dataname), hdferr)


end function writeHyperslabCharArray3D

!--------------------------------------------------------------------------
recursive function writeHyperslabCharArray4D(self, dataname, wdata, hdims, offset, &
                                                 dims, insert) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! writes a 4D hyperslab to the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: hdims(4)
integer(HSIZE_T),INTENT(IN)                :: offset(4)
integer(HSIZE_T),INTENT(IN)                :: dims(4)
character(kind=c_char),INTENT(IN)          :: wdata(dims(1),dims(2),dims(3),dims(4))
logical, OPTIONAL, INTENT(IN)              :: insert
integer(kind=irg)                          :: success

integer(HID_T)                             :: memspace, space, dset ! Handles
integer                                    :: hdferr, rnk

success = 0

rnk = 4
call h5screate_simple_f(rnk, hdims, space, hdferr)
call error_check_(self, 'writeHyperslabCharArray4D:h5screate_simple_f:'//trim(dataname), hdferr)

if (present(insert)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeHyperslabCharArray4D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_U8LE, space, dset, hdferr)
  call error_check_(self, 'writeHyperslabCharArray4D:h5dcreate_f:'//trim(dataname), hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call error_check_(self, 'writeHyperslabCharArray4D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'writeHyperslabCharArray4D:h5screate_simple_f:'//trim(dataname), hdferr)

call h5dwrite_f(dset, H5T_STD_U8LE, wdata, dims, hdferr, memspace, space)
call error_check_(self, 'writeHyperslabCharArray4D:h5dwrite_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'writeHyperslabCharArray4D:h5dclose_f:'//trim(dataname), hdferr)


end function writeHyperslabCharArray4D

!--------------------------------------------------------------------------
recursive function writeHyperslabIntegerArray2D(self, dataname, wdata, hdims, offset, &
                                                    dims, insert) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! writes a 2D hyperslab to the current file or group ID 

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: hdims(2)
integer(HSIZE_T),INTENT(IN)                :: offset(2)
integer(HSIZE_T),INTENT(IN)                :: dims(2)
integer(kind=irg),INTENT(IN)               :: wdata(dims(1),dims(2))
logical, OPTIONAL, INTENT(IN)              :: insert
integer(kind=irg)                          :: success

integer(HID_T)                             :: memspace, space, dset ! Handles
integer                                    :: hdferr, rnk

success = 0

rnk = 2
call h5screate_simple_f(rnk, hdims, space, hdferr)
call error_check_(self, 'writeHyperslabIntegerArray2D:h5screate_simple_f:'//trim(dataname), hdferr)

if (present(insert)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeHyperslabIntegerArray2D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_I32LE, space, dset, hdferr)
  call error_check_(self, 'writeHyperslabIntegerArray2D:h5dcreate_f:'//trim(dataname), hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call error_check_(self, 'writeHyperslabIntegerArray2D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'writeHyperslabIntegerArray2D:h5screate_simple_f:'//trim(dataname), hdferr)

call h5dwrite_f(dset, H5T_NATIVE_INTEGER, wdata, dims, hdferr, memspace, space)
call error_check_(self, 'writeHyperslabIntegerArray2D:h5dwrite_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'writeHyperslabIntegerArray2D:h5dclose_f:'//trim(dataname), hdferr)


end function writeHyperslabIntegerArray2D

!--------------------------------------------------------------------------
recursive function writeHyperslabIntegerArray3D(self, dataname, wdata, hdims, offset, &
                                                    dims, insert) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! writes a 3D hyperslab to the current file or group ID  

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: hdims(3)
integer(HSIZE_T),INTENT(IN)                :: offset(3)
integer(HSIZE_T),INTENT(IN)                :: dims(3)
integer(kind=irg),INTENT(IN)               :: wdata(dims(1),dims(2),dims(3))
logical, OPTIONAL, INTENT(IN)              :: insert
integer(kind=irg)                          :: success

integer(HID_T)                             :: memspace, space, dset ! Handles
integer                                    :: hdferr, rnk

success = 0

rnk = 3
call h5screate_simple_f(rnk, hdims, space, hdferr)
call error_check_(self, 'writeHyperslabIntegerArray3D:h5screate_simple_f:'//trim(dataname), hdferr)

if (present(insert)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeHyperslabIntegerArray3D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_I32LE, space, dset, hdferr)
  call error_check_(self, 'writeHyperslabIntegerArray3D:h5dcreate_f:'//trim(dataname), hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call error_check_(self, 'writeHyperslabIntegerArray3D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'writeHyperslabIntegerArray3D:h5screate_simple_f:'//trim(dataname), hdferr)

call h5dwrite_f(dset, H5T_NATIVE_INTEGER, wdata, dims, hdferr, memspace, space)
call error_check_(self, 'writeHyperslabIntegerArray3D:h5dwrite_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'writeHyperslabIntegerArray3D:h5dclose_f:'//trim(dataname), hdferr)


end function writeHyperslabIntegerArray3D

!--------------------------------------------------------------------------
recursive function writeHyperslabIntegerArray4D(self, dataname, wdata, hdims, offset, &
                                                    dims, insert) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! writes a 4D hyperslab to the current file or group ID  

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: hdims(4)
integer(HSIZE_T),INTENT(IN)                :: offset(4)
integer(HSIZE_T),INTENT(IN)                :: dims(4)
integer(kind=irg),INTENT(IN)               :: wdata(dims(1),dims(2),dims(3),dims(4))
logical, OPTIONAL, INTENT(IN)              :: insert
integer(kind=irg)                          :: success

integer(HID_T)                             :: memspace, space, dset ! Handles
integer                                    :: hdferr, rnk

success = 0

rnk = 4
call h5screate_simple_f(rnk, hdims, space, hdferr)
call error_check_(self, 'writeHyperslabintegerArray4D:h5screate_simple_f:'//trim(dataname), hdferr)

if (present(insert)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeHyperslabintegerArray4D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_STD_I32LE, space, dset, hdferr)
  call error_check_(self, 'writeHyperslabintegerArray4D:h5dcreate_f:'//trim(dataname), hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call error_check_(self, 'writeHyperslabintegerArray4D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'writeHyperslabintegerArray4D:h5screate_simple_f:'//trim(dataname), hdferr)

call h5dwrite_f(dset, H5T_NATIVE_INTEGER, wdata, dims, hdferr, memspace, space)
call error_check_(self, 'writeHyperslabintegerArray4D:h5dwrite_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'writeHyperslabintegerArray4D:h5dclose_f:'//trim(dataname), hdferr)


end function writeHyperslabIntegerArray4D

!--------------------------------------------------------------------------
recursive function writeHyperslabFloatArray2D(self, dataname, wdata, hdims, offset, &
                                                  dims, insert) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !!  writes a 2D hyperslab to the current file or group ID  

IMPLICIT NONE

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: hdims(2)
integer(HSIZE_T),INTENT(IN)                :: offset(2)
integer(HSIZE_T),INTENT(IN)                :: dims(2)
real(real_kind4),INTENT(IN)                :: wdata(dims(1),dims(2))
logical, OPTIONAL, INTENT(IN)              :: insert
integer(kind=irg)                          :: success

integer(HID_T)                             :: memspace, space, dset ! Handles
integer                                    :: hdferr, rnk

success = 0

rnk = 2
call h5screate_simple_f(rnk, hdims, space, hdferr)
call error_check_(self, 'writeHyperslabFloatArray2D:h5screate_simple_f:'//trim(dataname), hdferr)

if (present(insert)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeHyperslabFloatArray2D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
  call error_check_(self, 'writeHyperslabFloatArray2D:h5dcreate_f:'//trim(dataname), hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call error_check_(self, 'writeHyperslabFloatArray2D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'writeHyperslabFloatArray2D:h5screate_simple_f:'//trim(dataname), hdferr)

call h5dwrite_f(dset, H5T_NATIVE_REAL, wdata, dims, hdferr, memspace, space)
call error_check_(self, 'writeHyperslabFloatArray2D:h5dwrite_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'writeHyperslabFloatArray2D:h5dclose_f:'//trim(dataname), hdferr)


end function writeHyperslabFloatArray2D

!--------------------------------------------------------------------------
recursive function writeHyperslabFloatArray3D(self, dataname, wdata, hdims, offset, &
                                                  dims, insert) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! writes a 3D hyperslab to the current file or group ID  

IMPLICIT NONE

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: hdims(3)
integer(HSIZE_T),INTENT(IN)                :: offset(3)
integer(HSIZE_T),INTENT(IN)                :: dims(3)
real(real_kind4),INTENT(IN)                :: wdata(dims(1),dims(2),dims(3))
logical, OPTIONAL, INTENT(IN)              :: insert
integer(kind=irg)                          :: success

integer(HID_T)                             :: memspace, space, dset ! Handles
integer                                    :: hdferr, rnk

success = 0

rnk = 3
call h5screate_simple_f(rnk, hdims, space, hdferr)
call error_check_(self, 'writeHyperslabFloatArray3D:h5screate_simple_f:'//trim(dataname), hdferr)

if (present(insert)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeHyperslabFloatArray3D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
  call error_check_(self, 'writeHyperslabFloatArray3D:h5dcreate_f:'//trim(dataname), hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call error_check_(self, 'writeHyperslabFloatArray3D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'writeHyperslabFloatArray3D:h5screate_simple_f:'//trim(dataname), hdferr)

call h5dwrite_f(dset, H5T_NATIVE_REAL, wdata, dims, hdferr, memspace, space)
call error_check_(self, 'writeHyperslabFloatArray3D:h5dwrite_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'writeHyperslabFloatArray3D:h5dclose_f:'//trim(dataname), hdferr)


end function writeHyperslabFloatArray3D

!--------------------------------------------------------------------------
recursive function writeHyperslabFloatArray4D(self, dataname, wdata, hdims, offset, &
                                                  dims, insert) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! writes a 4D hyperslab to the current file or group ID  

IMPLICIT NONE

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: hdims(4)
integer(HSIZE_T),INTENT(IN)                :: offset(4)
integer(HSIZE_T),INTENT(IN)                :: dims(4)
real(real_kind4),INTENT(IN)                :: wdata(dims(1),dims(2),dims(3),dims(4))
logical, OPTIONAL, INTENT(IN)              :: insert
integer(kind=irg)                          :: success

integer(HID_T)                             :: memspace, space, dset ! Handles
integer                                    :: hdferr, rnk

success = 0

rnk = 4
call h5screate_simple_f(rnk, hdims, space, hdferr)
call error_check_(self, 'writeHyperslabFloatArray4D:h5screate_simple_f:'//trim(dataname), hdferr)

if (present(insert)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeHyperslabFloatArray4D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F32LE, space, dset, hdferr)
  call error_check_(self, 'writeHyperslabFloatArray4D:h5dcreate_f:'//trim(dataname), hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call error_check_(self, 'writeHyperslabFloatArray4D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'writeHyperslabFloatArray4D:h5screate_simple_f:'//trim(dataname), hdferr)

call h5dwrite_f(dset, H5T_NATIVE_REAL, wdata, dims, hdferr, memspace, space)
call error_check_(self, 'writeHyperslabFloatArray4D:h5dwrite_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'writeHyperslabFloatArray4D:h5dclose_f:'//trim(dataname), hdferr)


end function writeHyperslabFloatArray4D

!--------------------------------------------------------------------------
recursive function writeHyperslabDoubleArray2D(self, dataname, wdata, hdims, offset, &
                                                   dims, insert) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !!  writes a 2D hyperslab to the current file or group ID  

IMPLICIT NONE

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: hdims(2)
integer(HSIZE_T),INTENT(IN)                :: offset(2)
integer(HSIZE_T),INTENT(IN)                :: dims(2)
real(real_kind8),INTENT(IN)                :: wdata(dims(1),dims(2))
logical, OPTIONAL, INTENT(IN)              :: insert
integer(kind=irg)                          :: success

integer(HID_T)                             :: memspace, space, dset ! Handles
integer                                    :: hdferr, rnk

success = 0

rnk = 2
call h5screate_simple_f(rnk, hdims, space, hdferr)
call error_check_(self, 'writeHyperslabDoubleArray2D:h5screate_simple_f:'//trim(dataname), hdferr)

if (present(insert)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeHyperslabDoubleArray2D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
  call error_check_(self, 'writeHyperslabDoubleArray2D:h5dcreate_f:'//trim(dataname), hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call error_check_(self, 'writeHyperslabDoubleArray2D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'writeHyperslabDoubleArray2D:h5screate_simple_f:'//trim(dataname), hdferr)

call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, wdata, dims, hdferr, memspace, space)
call error_check_(self, 'writeHyperslabDoubleArray2D:h5dwrite_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'writeHyperslabDoubleArray2D:h5dclose_f:'//trim(dataname), hdferr)


end function writeHyperslabDoubleArray2D

!--------------------------------------------------------------------------
recursive function writeHyperslabDoubleArray3D(self, dataname, wdata, hdims, offset, &
                                                   dims, insert) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! writes a 3D hyperslab to the current file or group ID  

IMPLICIT NONE

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: hdims(3)
integer(HSIZE_T),INTENT(IN)                :: offset(3)
integer(HSIZE_T),INTENT(IN)                :: dims(3)
real(real_kind8),INTENT(IN)                :: wdata(dims(1),dims(2),dims(3))
logical, OPTIONAL, INTENT(IN)              :: insert
integer(kind=irg)                          :: success

integer(HID_T)                             :: memspace, space, dset ! Handles
integer                                    :: hdferr, rnk

success = 0

rnk = 3
call h5screate_simple_f(rnk, hdims, space, hdferr)
call error_check_(self, 'writeHyperslabDoubleArray3D:h5screate_simple_f:'//trim(dataname), hdferr)

if (present(insert)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeHyperslabDoubleArray3D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
  call error_check_(self, 'writeHyperslabDoubleArray3D:h5dcreate_f:'//trim(dataname), hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call error_check_(self, 'writeHyperslabDoubleArray3D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'writeHyperslabDoubleArray3D:h5screate_simple_f:'//trim(dataname), hdferr)

call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, wdata, dims, hdferr, memspace, space)
call error_check_(self, 'writeHyperslabDoubleArray3D:h5dwrite_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'writeHyperslabDoubleArray3D:h5dclose_f:'//trim(dataname), hdferr)


end function writeHyperslabDoubleArray3D

!--------------------------------------------------------------------------
recursive function writeHyperslabDoubleArray4D(self, dataname, wdata, hdims, offset, &
                                                   dims, insert) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! writes a 4D hyperslab to the current file or group ID  

IMPLICIT NONE

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: hdims(4)
integer(HSIZE_T),INTENT(IN)                :: offset(4)
integer(HSIZE_T),INTENT(IN)                :: dims(4)
real(real_kind8),INTENT(IN)                :: wdata(dims(1),dims(2),dims(3),dims(4))
logical, OPTIONAL, INTENT(IN)              :: insert
integer(kind=irg)                          :: success

integer(HID_T)                             :: memspace, space, dset ! Handles
integer                                    :: hdferr, rnk

success = 0

rnk = 4
call h5screate_simple_f(rnk, hdims, space, hdferr)
call error_check_(self, 'writeHyperslabDoubleArray4D:h5screate_simple_f:'//trim(dataname), hdferr)

if (present(insert)) then
  call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
  call error_check_(self, 'writeHyperslabDoubleArray4D:h5dopen_f:'//trim(dataname), hdferr)
else
  call h5dcreate_f(self%head%next%objectID, cstringify(dataname), H5T_IEEE_F64LE, space, dset, hdferr)
  call error_check_(self, 'writeHyperslabDoubleArray4D:h5dcreate_f:'//trim(dataname), hdferr)
end if

call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr)
call error_check_(self, 'writeHyperslabDoubleArray4D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'writeHyperslabDoubleArray4D:h5screate_simple_f:'//trim(dataname), hdferr)

call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, wdata, dims, hdferr, memspace, space)
call error_check_(self, 'writeHyperslabDoubleArray4D:h5dwrite_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'writeHyperslabDoubleArray4D:h5dclose_f:'//trim(dataname), hdferr)


end function writeHyperslabDoubleArray4D

!--------------------------------------------------------------------------
recursive function readHyperslabCharArray2D_(self, dataname, offset, dims) result(rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 2D hyperslab from the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: offset(2)
integer(HSIZE_T),INTENT(IN)                :: dims(2)
character(len=1,kind=c_char), dimension(:,:), allocatable, TARGET   :: rdata

integer(HID_T)                             :: memspace, space, dset ! Handles
integer(HSIZE_T)                           :: hdims(2), max_dims(2)
integer                                    :: hdferr, rnk

allocate(rdata(1:dims(1),1:dims(2)))

call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readHyperslabCharArray2D:h5dopen_f:'//trim(dataname), hdferr)

call h5dget_space_f(dset, space, hdferr)
call error_check_(self, 'readHyperslabCharArray2D:h5dget_space_f:'//trim(dataname), hdferr)


call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)
call error_check_(self, 'readHyperslabCharArray2D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


rnk = 2
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call error_check_(self, 'readHyperslabCharArray2D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'readHyperslabCharArray2D:h5screate_simple_f:'//trim(dataname), hdferr)


call h5dread_f(dset, H5T_STD_U8LE, rdata, hdims, hdferr, memspace, space)
call error_check_(self, 'readHyperslabCharArray2D:h5dread_f:'//trim(dataname), hdferr)


call h5sclose_f(space, hdferr)
call error_check_(self, 'readHyperslabCharArray2D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'readHyperslabCharArray2D:h5dclose_f:'//trim(dataname), hdferr)


end function readHyperslabCharArray2D_

!--------------------------------------------------------------------------
recursive function readHyperslabCharArray3D_(self, dataname3, offset3, dims3) result(rdata3)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 3D hyperslab from the current file or group ID 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname3
integer(HSIZE_T),INTENT(IN)                :: offset3(3)
integer(HSIZE_T),INTENT(IN)                :: dims3(3)
character(kind=c_char), dimension(:,:,:), allocatable, TARGET   :: rdata3

integer(HID_T)                             :: memspace, space, dset ! Handles
integer(HSIZE_T)                           :: hdims(3), max_dims(3)
integer                                    :: hdferr, rnk

allocate(rdata3(1:dims3(1),1:dims3(2),1:dims3(3)))

call h5dopen_f(self%head%next%objectID, cstringify(dataname3), dset, hdferr)
call error_check_(self, 'readHyperslabCharArray3D:h5dopen_f:'//trim(dataname3), hdferr)

call h5dget_space_f(dset, space, hdferr)
call error_check_(self, 'readHyperslabCharArray3D:h5dget_space_f:'//trim(dataname3), hdferr)

call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)
call error_check_(self, 'readHyperslabCharArray3D:h5sget_simple_extent_dims_f:'//trim(dataname3), hdferr)

rnk = 3
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset3, dims3, hdferr) 
call error_check_(self, 'readHyperslabCharArray3D:h5sselect_hyperslab_f:'//trim(dataname3), hdferr)

call h5screate_simple_f(rnk, dims3, memspace, hdferr)
call error_check_(self, 'readHyperslabCharArray3D:h5screate_simple_f:'//trim(dataname3), hdferr)

call h5dread_f(dset, H5T_STD_U8LE, rdata3, hdims, hdferr, memspace, space)
call error_check_(self, 'readHyperslabCharArray3D:h5dread_f:'//trim(dataname3), hdferr)

call h5sclose_f(space, hdferr)
call error_check_(self, 'readHyperslabCharArray3D:h5sclose_f:'//trim(dataname3), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'readHyperslabCharArray3D:h5dclose_f:'//trim(dataname3), hdferr)

end function readHyperslabCharArray3D_

!--------------------------------------------------------------------------
recursive function readHyperslabCharArray4D_(self, dataname4, offset4, dims4) result(rdata4)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 4D hyperslab from the current file or group ID  

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname4
integer(HSIZE_T),INTENT(IN)                :: offset4(4)
integer(HSIZE_T),INTENT(IN)                :: dims4(4)
character(len=1,kind=c_char), dimension(:,:,:,:), allocatable, TARGET   :: rdata4

integer(HID_T)                             :: memspace, space, dset ! Handles
integer(HSIZE_T)                           :: hdims(4), max_dims(4)
integer                                    :: hdferr, rnk

allocate(rdata4(1:dims4(1),1:dims4(2),1:dims4(3),1:dims4(4)))

call h5dopen_f(self%head%next%objectID, cstringify(dataname4), dset, hdferr)
call error_check_(self, 'readHyperslabCharArray4D:h5dopen_f:'//trim(dataname4), hdferr)

call h5dget_space_f(dset, space, hdferr)
call error_check_(self, 'readHyperslabCharArray4D:h5dget_space_f:'//trim(dataname4), hdferr)


call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)
call error_check_(self, 'readHyperslabCharArray4D:h5sget_simple_extent_dims_f:'//trim(dataname4), hdferr)


rnk = 4
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset4, dims4, hdferr) 
call error_check_(self, 'readHyperslabCharArray4D:h5sselect_hyperslab_f:'//trim(dataname4), hdferr)

call h5screate_simple_f(rnk, dims4, memspace, hdferr)
call error_check_(self, 'readHyperslabCharArray4D:h5screate_simple_f:'//trim(dataname4), hdferr)


call h5dread_f(dset, H5T_STD_U8LE, rdata4, hdims, hdferr, memspace, space)
call error_check_(self, 'readHyperslabCharArray4D:h5dread_f:'//trim(dataname4), hdferr)


call h5sclose_f(space, hdferr)
call error_check_(self, 'readHyperslabCharArray4D:h5sclose_f:'//trim(dataname4), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'readHyperslabCharArray4D:h5dclose_f:'//trim(dataname4), hdferr)


end function readHyperslabCharArray4D_

!--------------------------------------------------------------------------
recursive function readHyperslabIntegerArray2D_(self, dataname, offset, dims) result(rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 2D hyperslab from the current file or group ID  

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: offset(2)
integer(HSIZE_T),INTENT(IN)                :: dims(2)
integer, dimension(:,:), allocatable, TARGET :: rdata

integer(HID_T)                             :: memspace, space, dset ! Handles
integer(HSIZE_T)                           :: hdims(2), max_dims(2)
integer                                    :: hdferr, rnk

call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readHyperslabIntegerArray2D:h5dopen_f:'//trim(dataname), hdferr)

call h5dget_space_f(dset, space, hdferr)
call error_check_(self, 'readHyperslabIntegerArray2D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)
call error_check_(self, 'readHyperslabIntegerArray2D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


rnk = 2
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call error_check_(self, 'readHyperslabIntegerArray2D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'readHyperslabIntegerArray2D:h5screate_simple_f:'//trim(dataname), hdferr)


allocate(rdata(dims(1),dims(2)))
call h5dread_f(dset, H5T_NATIVE_INTEGER, rdata, hdims, hdferr, memspace, space)
call error_check_(self, 'readHyperslabIntegerArray2D:h5dread_f:'//trim(dataname), hdferr)


call h5sclose_f(space, hdferr)
call error_check_(self, 'readHyperslabIntegerArray2D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'readHyperslabIntegerArray2D:h5dclose_f:'//trim(dataname), hdferr)


end function readHyperslabIntegerArray2D_

!--------------------------------------------------------------------------
recursive function readHyperslabIntegerArray3D_(self, dataname, offset, dims) result(rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 3D hyperslab from the current file or group ID  

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: offset(3)
integer(HSIZE_T),INTENT(IN)                :: dims(3)
integer, dimension(:,:,:), allocatable, TARGET :: rdata

integer(HID_T)                             :: memspace, space, dset ! Handles
integer(HSIZE_T)                           :: hdims(3), max_dims(3)
integer                                    :: hdferr, rnk

call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readHyperslabIntegerArray3D:h5dopen_f:'//trim(dataname), hdferr)

call h5dget_space_f(dset, space, hdferr)
call error_check_(self, 'readHyperslabIntegerArray3D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)
call error_check_(self, 'readHyperslabIntegerArray3D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


rnk = 3
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call error_check_(self, 'readHyperslabIntegerArray3D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'readHyperslabIntegerArray3D:h5screate_simple_f:'//trim(dataname), hdferr)


allocate(rdata(dims(1),dims(2),dims(3)))
call h5dread_f(dset, H5T_NATIVE_INTEGER, rdata, hdims, hdferr, memspace, space)
call error_check_(self, 'readHyperslabIntegerArray3D:h5dread_f:'//trim(dataname), hdferr)


call h5sclose_f(space, hdferr)
call error_check_(self, 'readHyperslabIntegerArray3D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'readHyperslabIntegerArray3D:h5dclose_f:'//trim(dataname), hdferr)


end function readHyperslabIntegerArray3D_

!--------------------------------------------------------------------------
recursive function readHyperslabIntegerArray4D_(self, dataname, offset, dims) result(rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 4D hyperslab from the current file or group ID  

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: offset(4)
integer(HSIZE_T),INTENT(IN)                :: dims(4)
integer, dimension(:,:,:,:), allocatable, TARGET :: rdata

integer(HID_T)                             :: memspace, space, dset ! Handles
integer(HSIZE_T)                           :: hdims(4), max_dims(4)
integer                                    :: hdferr, rnk

call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readHyperslabIntegerArray4D:h5dopen_f:'//trim(dataname), hdferr)

call h5dget_space_f(dset, space, hdferr)
call error_check_(self, 'readHyperslabIntegerArray4D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)
call error_check_(self, 'readHyperslabIntegerArray4D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


rnk = 4
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call error_check_(self, 'readHyperslabIntegerArray4D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'readHyperslabIntegerArray4D:h5screate_simple_f:'//trim(dataname), hdferr)


allocate(rdata(dims(1),dims(2),dims(3),dims(4)))
call h5dread_f(dset, H5T_NATIVE_INTEGER, rdata, hdims, hdferr, memspace, space)
call error_check_(self, 'readHyperslabIntegerArray4D:h5dread_f:'//trim(dataname), hdferr)


call h5sclose_f(space, hdferr)
call error_check_(self, 'readHyperslabIntegerArray4D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'readHyperslabIntegerArray4D:h5dclose_f:'//trim(dataname), hdferr)


end function readHyperslabIntegerArray4D_

!--------------------------------------------------------------------------
recursive function readHyperslabFloatArray2D_(self, dataname, offset, dims) result(rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 2D hyperslab from the current file or group ID  

IMPLICIT NONE

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: offset(2)
integer(HSIZE_T),INTENT(IN)                :: dims(2)
real(real_kind4), dimension(:,:), allocatable, TARGET    :: rdata

integer(HID_T)                             :: memspace, space, dset ! Handles
integer(HSIZE_T)                           :: hdims(2), max_dims(2)
integer                                    :: hdferr, rnk

call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readHyperslabFloatArray2D:h5dopen_f:'//trim(dataname), hdferr)

call h5dget_space_f(dset, space, hdferr)
call error_check_(self, 'readHyperslabFloatArray2D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)
call error_check_(self, 'readHyperslabFloatArray2D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


rnk = 2
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call error_check_(self, 'readHyperslabFloatArray2D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'readHyperslabFloatArray2D:h5screate_simple_f:'//trim(dataname), hdferr)


allocate(rdata(dims(1),dims(2)))
call h5dread_f(dset, H5T_NATIVE_REAL, rdata, hdims, hdferr, memspace, space)
call error_check_(self, 'readHyperslabFloatArray2D:h5dread_f:'//trim(dataname), hdferr)


call h5sclose_f(space, hdferr)
call error_check_(self, 'readHyperslabFloatArray2D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'readHyperslabFloatArray2D:h5dclose_f:'//trim(dataname), hdferr)


end function readHyperslabFloatArray2D_

!--------------------------------------------------------------------------
recursive function readHyperslabFloatArray3D_(self, dataname, offset, dims) result(rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 3D hyperslab from the current file or group ID 

IMPLICIT NONE

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: offset(3)
integer(HSIZE_T),INTENT(IN)                :: dims(3)
real(real_kind4), dimension(:,:,:), allocatable, TARGET  :: rdata

integer(HID_T)                             :: memspace, space, dset ! Handles
integer(HSIZE_T)                           :: hdims(3), max_dims(3)
integer                                    :: hdferr, rnk

call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readHyperslabFloatArray3D:h5dopen_f:'//trim(dataname), hdferr)

call h5dget_space_f(dset, space, hdferr)
call error_check_(self, 'readHyperslabFloatArray3D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)
call error_check_(self, 'readHyperslabFloatArray3D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


rnk = 3
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call error_check_(self, 'readHyperslabFloatArray3D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'readHyperslabFloatArray3D:h5screate_simple_f:'//trim(dataname), hdferr)


allocate(rdata(dims(1),dims(2),dims(3)))
call h5dread_f(dset, H5T_NATIVE_REAL, rdata, hdims, hdferr, memspace, space)
call error_check_(self, 'readHyperslabFloatArray3D:h5dread_f:'//trim(dataname), hdferr)


call h5sclose_f(space, hdferr)
call error_check_(self, 'readHyperslabFloatArray3D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'readHyperslabFloatArray3D:h5dclose_f:'//trim(dataname), hdferr)


end function readHyperslabFloatArray3D_

!--------------------------------------------------------------------------
recursive function readHyperslabFloatArray4D_(self, dataname, offset, dims) result(rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !!  reads and returns a 4D hyperslab from the current file or group ID  

IMPLICIT NONE

integer,parameter                          :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: offset(4)
integer(HSIZE_T),INTENT(IN)                :: dims(4)
real(real_kind4), dimension(:,:,:,:), allocatable, TARGET:: rdata

integer(HID_T)                             :: memspace, space, dset ! Handles
integer(HSIZE_T)                           :: hdims(4), max_dims(4)
integer                                    :: hdferr, rnk

call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readHyperslabFloatArray4D:h5dopen_f:'//trim(dataname), hdferr)

call h5dget_space_f(dset, space, hdferr)
call error_check_(self, 'readHyperslabFloatArray4D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)
call error_check_(self, 'readHyperslabFloatArray4D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


rnk = 4
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call error_check_(self, 'readHyperslabFloatArray4D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'readHyperslabFloatArray4D:h5screate_simple_f:'//trim(dataname), hdferr)


allocate(rdata(dims(1),dims(2),dims(3),dims(4)))
call h5dread_f(dset, H5T_NATIVE_REAL, rdata, hdims, hdferr, memspace, space)
call error_check_(self, 'readHyperslabFloatArray4D:h5dread_f:'//trim(dataname), hdferr)


call h5sclose_f(space, hdferr)
call error_check_(self, 'readHyperslabFloatArray4D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'readHyperslabFloatArray4D:h5dclose_f:'//trim(dataname), hdferr)


end function readHyperslabFloatArray4D_

!--------------------------------------------------------------------------
recursive function readHyperslabDoubleArray2D_(self, dataname, offset, dims) result(rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 2D hyperslab from the current file or group ID  

IMPLICIT NONE

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: offset(2)
integer(HSIZE_T),INTENT(IN)                :: dims(2)
real(real_kind8), dimension(:,:), allocatable, TARGET    :: rdata

integer(HID_T)                             :: memspace, space, dset ! Handles
integer(HSIZE_T)                           :: hdims(2), max_dims(2)
integer                                    :: hdferr, rnk

call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readHyperslabDoubleArray2D:h5dopen_f:'//trim(dataname), hdferr)

call h5dget_space_f(dset, space, hdferr)
call error_check_(self, 'readHyperslabDoubleArray2D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)
call error_check_(self, 'readHyperslabDoubleArray2D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


rnk = 2
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call error_check_(self, 'readHyperslabDoubleArray2D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'readHyperslabDoubleArray2D:h5screate_simple_f:'//trim(dataname), hdferr)


allocate(rdata(dims(1),dims(2)))
call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, hdims, hdferr, memspace, space)
call error_check_(self, 'readHyperslabDoubleArray2D:h5dread_f:'//trim(dataname), hdferr)


call h5sclose_f(space, hdferr)
call error_check_(self, 'readHyperslabDoubleArray2D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'readHyperslabDoubleArray2D:h5dclose_f:'//trim(dataname), hdferr)


end function readHyperslabDoubleArray2D_

!--------------------------------------------------------------------------
recursive function readHyperslabDoubleArray3D_(self, dataname, offset, dims) result(rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 3D hyperslab from the current file or group ID  

IMPLICIT NONE

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: offset(3)
integer(HSIZE_T),INTENT(IN)                :: dims(3)
real(real_kind8), dimension(:,:,:), allocatable, TARGET  :: rdata

integer(HID_T)                             :: memspace, space, dset ! Handles
integer(HSIZE_T)                           :: hdims(3), max_dims(3)
integer                                    :: hdferr, rnk

call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'readHyperslabDoubleArray3D:h5dopen_f:'//trim(dataname), hdferr)

call h5dget_space_f(dset, space, hdferr)
call error_check_(self, 'readHyperslabDoubleArray3D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)
call error_check_(self, 'readHyperslabDoubleArray3D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


rnk = 3
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call error_check_(self, 'readHyperslabDoubleArray3D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'readHyperslabDoubleArray3D:h5screate_simple_f:'//trim(dataname), hdferr)


allocate(rdata(dims(1),dims(2),dims(3)))
call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, hdims, hdferr, memspace, space)
call error_check_(self, 'readHyperslabDoubleArray3D:h5dread_f:'//trim(dataname), hdferr)


call h5sclose_f(space, hdferr)
call error_check_(self, 'readHyperslabDoubleArray3D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'readHyperslabDoubleArray3D:h5dclose_f:'//trim(dataname), hdferr)


end function readHyperslabDoubleArray3D_

!--------------------------------------------------------------------------
recursive function readHyperslabDoubleArray4D_(self, dataname, offset, dims) result(rdata)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! reads and returns a 4D hyperslab from the current file or group ID  

IMPLICIT NONE

integer,parameter                          :: real_kind8 = SELECTED_REAL_KIND(Fortran_REAL_8)

class(HDF_T),INTENT(INOUT)                 :: self
character(fnlen),INTENT(IN)                :: dataname
integer(HSIZE_T),INTENT(IN)                :: offset(4)
integer(HSIZE_T),INTENT(IN)                :: dims(4)
real(real_kind8), dimension(:,:,:,:), allocatable, TARGET:: rdata

integer(HID_T)                             :: memspace, space, dset ! Handles
integer(HSIZE_T)                           :: hdims(4), max_dims(4)
integer                                    :: hdferr, rnk

call h5dopen_f(self%head%next%objectID, cstringify(dataname), dset, hdferr)
call error_check_(self, 'hdf_readHyperslabDoubleArray4D:h5dopen_f:'//trim(dataname), hdferr)

call h5dget_space_f(dset, space, hdferr)
call error_check_(self, 'hdf_readHyperslabDoubleArray4D:h5dget_space_f:'//trim(dataname), hdferr)

call h5sget_simple_extent_dims_f(space, hdims, max_dims, hdferr)
call error_check_(self, 'hdf_readHyperslabDoubleArray4D:h5sget_simple_extent_dims_f:'//trim(dataname), hdferr)


rnk = 4
call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offset, dims, hdferr) 
call error_check_(self, 'hdf_readHyperslabDoubleArray4D:h5sselect_hyperslab_f:'//trim(dataname), hdferr)

call h5screate_simple_f(rnk, dims, memspace, hdferr)
call error_check_(self, 'hdf_readHyperslabDoubleArray4D:h5screate_simple_f:'//trim(dataname), hdferr)


allocate(rdata(dims(1),dims(2),dims(3),dims(4)))
call h5dread_f(dset, H5T_NATIVE_DOUBLE, rdata, hdims, hdferr, memspace, space)
call error_check_(self, 'hdf_readHyperslabDoubleArray4D:h5dread_f:'//trim(dataname), hdferr)


call h5sclose_f(space, hdferr)
call error_check_(self, 'hdf_readHyperslabDoubleArray4D:h5sclose_f:'//trim(dataname), hdferr)

call h5dclose_f(dset, hdferr)
call error_check_(self, 'hdf_readHyperslabDoubleArray4D:h5dclose_f:'//trim(dataname), hdferr)


end function readHyperslabDoubleArray4D_


!--------------------------------------------------------------------------
recursive function CheckFixedLengthflag_(self, dataset) result(itis)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! returns TRUE if there is a FixedLength=1 data set in the current EMheader 

use HDF5
use h5lt
 
IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                      :: self
character(fnlen),INTENT(IN)                     :: dataset

logical                                         :: itis

integer(kind=irg)                               :: FL, hdferr, i

! we assume this file has variable length strings (the default for EMsoft)
itis = .FALSE.

! look for the data set
i = h5ltfind_dataset_f(self%head%next%ObjectID, trim(dataset))

if (i.eq.1) then 
  call readDatasetInteger_(self, dataset, hdferr, FL)
  call error_check_(self, 'CheckFixedLengthflag:readDatasetInteger:'//trim(dataset), hdferr)

  if (FL.eq.1) then 
    itis = .TRUE.
  end if
end if

! and set the saved FixedLength flag (private variable)
FixedLengthflag = itis

end function CheckFixedLengthflag_

!--------------------------------------------------------------------------
recursive subroutine resetFixedLengthflag_(self)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! sets the FixedLengthflag to .FALSE. 

IMPLICIT NONE

class(HDF_T),INTENT(INOUT):: self

FixedLengthflag = .FALSE.

end subroutine resetFixedLengthflag_

!--------------------------------------------------------------------------
recursive SUBROUTINE h5_write_pseudo_bse_image_(self, fname, dsetnm, hdferr, wdata)
  !! author: Patrick G. Callahan, UCSB 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! Write an hdf5 file containing psuedo/virtual bse images 

  USE ISO_C_BINDING
  IMPLICIT NONE

  INTEGER,PARAMETER                                       :: real_kind = SELECTED_REAL_KIND(4)

  class(HDF_T),INTENT(INOUT)                              :: self
  CHARACTER(LEN=fnlen), INTENT(IN)                        :: fname  !NAME OF THE HDF5 FILE YOU WANT TO CREATE
  CHARACTER(LEN=fnlen), INTENT(IN)                        :: dsetnm !NAME OF THE DATASET YOU WANT TO WRITE IN THE HDF5 FILE
  REAL(KIND=4), DIMENSION(:,:,:),  TARGET, INTENT(INOUT)  :: wdata

  INTEGER                                                 :: rnk
                          
  INTEGER,INTENT(OUT)                                     :: hdferr   !hdferr flag    
  INTEGER(HID_T)                                          :: fid, dsetid, spaceid, memspace !FILE ID, DATASET ID, DATASPACE ID, MEMORY SPACE
  INTEGER(HSIZE_T), ALLOCATABLE                           :: dims(:)!, max_dims(2)    
  TYPE(C_PTR)                                             :: buff

  buff = C_LOC(wdata(1,1,1))
  dims = shape(wdata)
  rnk = rank(wdata)
  write(*,*) 'dims =', dims!, 'rank = ',rnk

  !CREATE THE FILE
  CALL h5fcreate_f(fname, H5F_ACC_TRUNC_F, fid, hdferr)
  !CREATE THE DATASPACE
  CALL h5screate_simple_f(rnk, dims, spaceid, hdferr)
  !CREATE THE DATASET
  CALL h5dcreate_f(fid, dsetnm, H5T_IEEE_F32LE, spaceid, dsetid, hdferr)
  !WRITE THE DATASET
  CALL h5dwrite_f(dsetid, H5T_NATIVE_REAL, buff, hdferr)
  

  !CLEANUP
  !CLOSE THE DATASPACE
  CALL h5sclose_f(spaceid, hdferr)
  !CLOSE THE DATASET
  CALL h5dclose_f(dsetid, hdferr)
!  !CLOSE THE DATASPACE
!  CALL h5sclose_f(spaceid, hdferr)
  !CLOSE THE FILE
  CALL h5fclose_f(fid, hdferr)

END SUBROUTINE h5_write_pseudo_bse_image_

!--------------------------------------------------------------------------
recursive SUBROUTINE h5_tsl_read_ebsd_pattern_(self, fname, dsetnm, hdferr, rdata, offset, szx, szy)
  !! author: Patrick G. Callahan, UCSB 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! Read a single EBSD pattern from a TSL hdf5 file
  !!
  !! Helper routine for h5tslpbse. Reads a single pattern because
  !! many of our datasets are 100s of GB so can't load them all in memory.
  !! This will also be useful for dictionary indexing of our scans that are
  !! large. In the future we should read a number of patterns at once depending 
  !! on memory limits and work with them, then move to the next set of patterns.

  USE ISO_C_BINDING
  IMPLICIT NONE

  class(HDF_T),INTENT(INOUT)                        :: self
  CHARACTER(LEN=fnlen), INTENT(IN)                  :: fname  !NAME OF THE HDF5 FILE YOU ARE READING FROM
  CHARACTER(LEN=fnlen), INTENT(IN)                  :: dsetnm !NAME OF THE DATASET YOU WANT TO READ
  INTEGER,INTENT(OUT)                               :: hdferr   !hdferr flag
  INTEGER(HID_T)                                    :: fid, dsetid, spaceid, memspace !FILE ID, DATASET ID, DATASPACE ID, MEMORY SPACE
  INTEGER, DIMENSION(:,:), TARGET, INTENT(INOUT)    :: rdata  ! variable for data to be read into
  INTEGER(KIND=8),INTENT(IN)                        :: offset ! 1d index of pattern in nRow x nColumns sized vector
  INTEGER,INTENT(IN)                                :: szx, szy ! nColumns, nRows in the scan

  TYPE(C_PTR)                                       :: buff
  INTEGER(KIND=8)                                   :: start(3), cnt(3) 
  INTEGER(KIND=4)                                   :: rnk 

  !! start is the offset of the starting element of the specificied hyperslab
  start(1) = 0
  start(2) = 0
  start(3) = offset

  !! how many blocks to select from the dataspace in each dimension
  cnt(1) = szx
  cnt(2) = szy
  cnt(3) = 1

  !Open the hdf5 file
  CALL h5fopen_f(fname, H5F_ACC_RDONLY_F, fid, hdferr)
  !Open an existing dataset
  CALL h5dopen_f(fid, dsetnm, dsetid, hdferr)  
  !GET DATASPACE AND ALLOCATE MEMORY FOR READ BUFFER
  CALL h5dget_space_f(dsetid, spaceid, hdferr)
  !READ THE DATA
  buff = C_LOC(rdata)
  ! Select the region in the hyperslab 
  CALL h5sselect_hyperslab_f( spaceid, H5S_SELECT_SET_F, start, cnt, hdferr)
  rnk=3
  !CREATE A SIMPLE DATASPACE AND OPEN IT FOR ACCESS
  CALL h5screate_simple_f( rnk, cnt, memspace, hdferr)
  CALL h5dread_f(dsetid, H5T_NATIVE_INTEGER, rdata, cnt, hdferr, memspace, spaceid)
  
  !CLOSE MEMSPACE
  CALL h5sclose_f(memspace,hdferr)
  !CLOSE THE DATASPACE
  CALL h5sclose_f(spaceid,hdferr) 
  !CLOSE THE DATASET
  CALL h5dclose_f(dsetid,hdferr)

END SUBROUTINE h5_tsl_read_ebsd_pattern_

!--------------------------------------------------------------------------
recursive SUBROUTINE h5_read_integer_dataset_(self, fname, dsetnm, hdferr, rdata)
  !! author: Patrick G. Callahan, UCSB 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! Read a single EBSD pattern from a TSL hdf5 file 

  USE ISO_C_BINDING
  IMPLICIT NONE

  class(HDF_T),INTENT(INOUT)        :: self
  CHARACTER(LEN=fnlen), INTENT(IN)  :: fname  !NAME OF THE HDF5 FILE YOU ARE READING FROM
  CHARACTER(LEN=fnlen), INTENT(IN)  :: dsetnm !NAME OF THE DATASET YOU WANT TO READ
  INTEGER,INTENT(OUT)               :: hdferr   !hdferr flag
  INTEGER(HID_T)                    :: fid      !FILE ID
  INTEGER(HID_T)                    :: dsetid   !DATASET ID
  INTEGER(HID_T)                    :: spaceid  !DATASPACE ID
  INTEGER, TARGET, INTENT(OUT)      :: rdata
  TYPE(C_PTR)                       :: buff

  !Open the hdf5 file
  CALL h5fopen_f(fname, H5F_ACC_RDONLY_F, fid, hdferr)

  !Open an existing dataset
  CALL h5dopen_f(fid, dsetnm, dsetid, hdferr)  

  !GET DATASPACE AND ALLOCATE MEMORY FOR READ BUFFER
  CALL h5dget_space_f(dsetid, spaceid, hdferr)

  !READ THE DATA
  buff = C_LOC(rdata)
  CALL h5dread_f(dsetid,H5T_NATIVE_INTEGER, buff, hdferr)
  !CLOSE THE DATASPACE
  CALL h5sclose_f(spaceid,hdferr)     
  !CLOSE THE DATASET
  CALL h5dclose_f(dsetid,hdferr)

END SUBROUTINE h5_read_integer_dataset_

!--------------------------------------------------------------------------
recursive function addStringAttributeToGroup_(self, dataname, stratt, overwrite) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! add a string attribute to the current level in the HDF file 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                              :: self
character(fnlen),INTENT(IN)                             :: dataname
character(len=fnlen, KIND=c_char),INTENT(INOUT)         :: stratt 
logical,INTENT(IN),OPTIONAL                             :: overwrite
integer(kind=irg)                                       :: success

integer,parameter                                       :: real_kind4 = SELECTED_REAL_KIND(Fortran_REAL_4)
integer(HID_T)                                          :: aspace_id, dset, atype_id, attr_id ! Handles
integer                                                 :: hdferr, rnk
integer(SIZE_T)                                         :: attrlen
integer(HSIZE_T), DIMENSION(1:1)                        :: dims
integer(HSIZE_T), DIMENSION(1)                          :: data_dims

success = 0

dims(1) = 1

attrlen = len_trim(stratt)
attrlen = attrlen+1
stratt(attrlen:attrlen) = C_NULL_CHAR

! Create dataspace.
rnk = 1
call h5screate_simple_f(rnk, dims, aspace_id, hdferr)
call error_check_(self, 'addStringAttributeToGroup_:h5screate_simple_f:'//trim(dataname), hdferr)
!
! Create datatype for the attribute.
!
call h5tcopy_f(H5T_NATIVE_CHARACTER, atype_id, hdferr)
call error_check_(self, 'addStringAttributeToGroup_:h5tcopy_f:'//trim(dataname), hdferr)

call h5tset_size_f(atype_id, attrlen, hdferr)
call error_check_(self, 'addStringAttributeToGroup_:h5tset_size_f:'//trim(dataname), hdferr)
!
! Create the attribute and write the string data to it.
!
if (present(overwrite)) then
  call h5aopen_f(self%head%next%objectID, cstringify(dataname), attr_id, hdferr)
  call error_check_(self, 'addStringAttributeToGroup_:h5aopen_f:'//trim(dataname), hdferr)
else
  call h5acreate_f(self%head%next%objectID, cstringify(dataname), atype_id, aspace_id, attr_id, hdferr)
  call error_check_(self, 'addStringAttributeToGroup_:h5acreate_f:'//trim(dataname), hdferr)
end if

if (hdferr.lt.0) then
  success = -1
end if

data_dims(1) = 1
call h5awrite_f(attr_id, atype_id, stratt, data_dims, hdferr )
call error_check_(self, 'addStringAttributeToGroup_:h5awrite_f:'//trim(dataname), hdferr)

if (hdferr.lt.0) then
  success = -1
end if
!
! Close and release resources.
!
call h5aclose_f(attr_id , hdferr)
call error_check_(self, 'addStringAttributeToGroup_:h5aclose_f:'//trim(dataname), hdferr)

call h5sclose_f(aspace_id, hdferr)
call error_check_(self, 'addStringAttributeToGroup_:h5sclose_f:'//trim(dataname), hdferr)

! that's it

end function addStringAttributeToGroup_

!--------------------------------------------------------------------------
recursive function getStringAttributeFromGroup_(self, dataname, stratt, slen) result(success)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! read a string attribute from the current level in the HDF file 

use ISO_C_BINDING

class(HDF_T),INTENT(INOUT)                              :: self
character(fnlen),INTENT(IN)                             :: dataname
integer(SIZE_T),INTENT(IN)                              :: slen
character(len=slen, KIND=c_char),INTENT(INOUT)          :: stratt 
integer(kind=irg)                                       :: success

integer(HID_T)                                          :: aspace_id, filetype, atype_id, attr_id, memtype ! Handles
integer                                                 :: hdferr, rnk
integer(SIZE_T)                                         :: attrlen
INTEGER(HSIZE_T), DIMENSION(1:1)                        :: maxdims
INTEGER(hsize_t), DIMENSION(1:1)                        :: dims
CHARACTER(LEN=slen), DIMENSION(:), ALLOCATABLE, TARGET  :: rdata
INTEGER(SIZE_T)                                         :: size

INTEGER, DIMENSION(:), POINTER                          :: ptr_r 
TYPE(C_PTR)                                             :: f_ptr
  

dims(1) = slen

success = 0

! open the attribute for this group
  call h5aopen_f(self%head%next%objectID, cstringify(dataname), attr_id, hdferr)
  call error_check_(self, 'HDF_getStringAttributeFromGroup:h5aopen_f:'//trim(dataname), hdferr)

  ! Get the datatype and its size.
  !
  CALL H5Aget_type_f(attr_id, filetype, hdferr)
  CALL H5Tget_size_f(filetype, size, hdferr)

  ! Make sure the declared length is large enough
  IF(size.GT.slen+1)THEN
     PRINT*,'ERROR:mod_HDFsupport:getStringAttributeFromGroup_: Character LEN is too small'
     STOP
  ENDIF
  !
  ! Get dataspace and allocate memory for read buffer.
  ! 
  CALL H5Aget_space_f(attr_id, aspace_id, hdferr)
  CALL H5Sget_simple_extent_dims_f(aspace_id, dims, maxdims, hdferr)

  ALLOCATE(rdata(1:dims(1)))
  !
  ! Create the memory datatype.
  !
  CALL H5Tcopy_f(H5T_FORTRAN_S1, memtype, hdferr)
  CALL H5Tset_size_f(memtype, slen, hdferr)
  !
  ! Read the data.
  !
  f_ptr = C_LOC(rdata(1)(1:1))
  CALL H5Aread_f(attr_id, memtype, f_ptr, hdferr)
  stratt = trim(rdata(1))
  !
  ! Close and release resources.
  !
  CALL H5Aclose_f(attr_id, hdferr)
  CALL H5Sclose_f(aspace_id, hdferr)
  CALL H5Tclose_f(memtype, hdferr)

end function getStringAttributeFromGroup_

!--------------------------------------------------------------------------
subroutine read2DImage_(self, dataset, image, numx, numy)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! read a gray scale image from the HDF5 file  

use h5im
use h5lt

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)             :: self
character(fnlen),INTENT(IN)            :: dataset
integer(kind=irg),INTENT(IN)           :: numx
integer(kind=irg),INTENT(IN)           :: numy
integer(kind=irg),INTENT(INOUT)        :: image(numx,numy)

integer(kind=irg),allocatable          :: vec(:)
integer(kind=irg)                      :: hdferr

! read the image from the file
allocate(vec(numx*numy))
call h5imread_image_f(self%head%next%objectID,dataset,vec,hdferr)

! reorganize it into a regular image
image = reshape( vec, (/ numx, numy/) )

end subroutine read2DImage_


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! EMsoft-specific procedures 
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive subroutine error_check_(self, OffendingRoutine, error, Fatal)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! deal with an HDF error 

use mod_io

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)   :: self
character(LEN=*),INTENT(IN)  :: OffendingRoutine   ! name of offending routine + message
integer(kind=irg),INTENT(IN) :: error              ! returned error code
logical,OPTIONAL,INTENT(IN)  :: Fatal              ! if true, then report the error and stop the program

type(IO_T)                   :: Message 
integer(kind=irg)            :: io_int(1)

if (error.lt.0) then 
  io_int(1) = error
  call Message%WriteValue('Error code : ',io_int,1)
  
  call Message%printMessage('   returned by routine '//OffendingRoutine,frm="(A)")
  
  if (present(Fatal)) then
    if (Fatal.eqv..TRUE.) STOP 'Unrecoverable Error' ! this is not very graceful, but it'll do the job for now ...
  end if
end if

end subroutine error_check_

!--------------------------------------------------------------------------
recursive subroutine writeEMheader_(self, dstr, tstrb, tstre, prn, dataname)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! The EMheader is an HDF group that contains the following basic dataset strings
  !!
  !! EMsoft version       : scversion from local.f90
  !! execution date       : dstr
  !! start time           : tstr1
  !! end time             : tstr2
  !! program name         : prn
  !! user name            : EMsoft%getConfigParameter('Username') 
  !! user location        : EMsoft%getConfigParameter('Userlocation') 
  !! user email           : EMsoft%getConfigParameter('Useremail') 
  !! computer name        : read via system call hostnm() 

use mod_EMsoft
use mod_timing
use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)               :: self
character(11),INTENT(INOUT)              :: dstr
character(15),INTENT(IN)                 :: tstrb
character(15),INTENT(IN)                 :: tstre
character(fnlen),INTENT(IN)              :: prn
character(fnlen),INTENT(IN),OPTIONAL     :: dataname

type(EMsoft_T)                           :: EMsoft 
type(Timing_T)                           :: Timing
integer                                  :: hdferr ! error flag
integer                                  :: i,ic,nlen 
character(100)                           :: c
character(fnlen)                         :: line, groupname
character(fnlen,kind=c_char)             :: line2(1)
logical                                  :: g_exists, overwrite=.TRUE.

! create and open the EMheader group
groupname = SC_EMheader
hdferr = self%createGroup_(groupname)
call error_check_(self, 'writeEMheader_:HDF_createGroup_:'//trim(groupname), hdferr)

! create and open the dataname group to allow for different data sets in the same file
if (PRESENT(dataname).eqv..TRUE.) then 
  groupname = trim(dataname)
  hdferr = self%createGroup_(groupname)
  call error_check_(self, 'writeEMheader_:HDF_createGroup_:'//trim(groupname), hdferr)
end if

! version number /EMheader/Version 'character'
line = 'Version'
line2(1) = trim(EMsoft%getConfigParameter('EMsoftversion'))
call H5Lexists_f(self%head%next%objectID,trim(line),g_exists, hdferr)
call error_check_(self, 'writeEMheader_:H5Lexists_f:'//trim(line), hdferr)

if (g_exists) then 
  hdferr = self%writeDatasetStringArray(line, line2, 1, overwrite)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line)//':overwrite', hdferr)
else
  hdferr = self%writeDatasetStringArray(line, line2, 1)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line), hdferr)
end if

! execution data /EMheader/Date 'character'
line = 'Date'
line2(1) = dstr
call H5Lexists_f(self%head%next%objectID,trim(line),g_exists, hdferr)
call error_check_(self, 'writeEMheader_:H5Lexists_f:'//trim(line), hdferr)

if (g_exists) then 
  hdferr = self%writeDatasetStringArray(line, line2, 1, overwrite)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line)//':overwrite', hdferr)
else
  hdferr = self%writeDatasetStringArray(line, line2, 1)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line), hdferr)
end if

! start time /EMheader/StartTime 'character'
line = 'StartTime'
line2(1) = dstr//', '//tstrb
call H5Lexists_f(self%head%next%objectID,trim(line),g_exists, hdferr)
call error_check_(self, 'writeEMheader_:H5Lexists_f:'//trim(line), hdferr)

if (g_exists) then 
  hdferr = self%writeDatasetStringArray(line, line2, 1, overwrite)

else
  hdferr = self%writeDatasetStringArray(line, line2, 1)

end if

! stop time /EMheader/StopTime 'character'; this is often updated at the end of a run
! since the end date can be different from the start date, especially for long runs, we get a new dstr string
dstr = Timing%getDateString()
line = 'StopTime'
line2(1) = dstr//', '//tstre
call H5Lexists_f(self%head%next%objectID,trim(line),g_exists, hdferr)
call error_check_(self, 'writeEMheader_:H5Lexists_f:'//trim(line), hdferr)

if (g_exists) then 
  hdferr = self%writeDatasetStringArray(line, line2, 1, overwrite)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line)//':overwrite', hdferr)
else
  hdferr = self%writeDatasetStringArray(line, line2, 1)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line), hdferr)
end if

! program name /EMheader/ProgramName 'character'
line = 'ProgramName'
line2(1) = prn 
call H5Lexists_f(self%head%next%objectID,trim(line),g_exists, hdferr)
call error_check_(self, 'writeEMheader_:H5Lexists_f:'//trim(line), hdferr)

if (g_exists) then 
  hdferr = self%writeDatasetStringArray(line, line2, 1, overwrite)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line)//':overwrite', hdferr)
else
  hdferr = self%writeDatasetStringArray(line, line2, 1)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line), hdferr)
end if

! user name /EMheader/UserName 'character'
line = 'UserName'
line2(1) = trim(EMsoft%getConfigParameter('UserName'))
call H5Lexists_f(self%head%next%objectID,trim(line),g_exists, hdferr)
call error_check_(self, 'writeEMheader_:H5Lexists_f:'//trim(line), hdferr)

if (g_exists) then 
  hdferr = self%writeDatasetStringArray(line, line2, 1, overwrite)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line)//':overwrite', hdferr)
else
  hdferr = self%writeDatasetStringArray(line, line2, 1)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line), hdferr)
end if

! user location /EMheader/UserLocation 'character'
line = 'UserLocation'
line2(1) = trim(EMsoft%getConfigParameter('UserLocation'))
call H5Lexists_f(self%head%next%objectID,trim(line),g_exists, hdferr)
call error_check_(self, 'writeEMheader_:H5Lexists_f:'//trim(line), hdferr)

if (g_exists) then 
  hdferr = self%writeDatasetStringArray(line, line2, 1, overwrite)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line)//':overwrite', hdferr)
else
  hdferr = self%writeDatasetStringArray(line, line2, 1)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line), hdferr)
end if

! user email /EMheader/UserEmail 'character'
line = 'UserEmail'
line2(1) = trim(EMsoft%getConfigParameter('UserEmail'))
call H5Lexists_f(self%head%next%objectID,trim(line),g_exists, hdferr)
call error_check_(self, 'writeEMheader_:H5Lexists_f:'//trim(line), hdferr)

if (g_exists) then 
  hdferr = self%writeDatasetStringArray(line, line2, 1, overwrite)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line)//':overwrite', hdferr)
else
  hdferr = self%writeDatasetStringArray(line, line2, 1)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line), hdferr)
end if

! hostname /EMheader/HostName 'character'
call hostnm(c)

! lowercase it
nlen = len(c) 
do i=1,nlen 
   ic = ichar(c(i:i)) 
   if (ic >= 65 .and. ic < 90) c(i:i) = char(ic+32) 
end do 
line = 'HostName'
line2(1) = c
call H5Lexists_f(self%head%next%objectID,trim(line),g_exists, hdferr)
call error_check_(self, 'writeEMheader_:H5Lexists_f:'//trim(line), hdferr)

if (g_exists) then 
  hdferr = self%writeDatasetStringArray(line, line2, 1, overwrite)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line)//':overwrite', hdferr)
else
  hdferr = self%writeDatasetStringArray(line, line2, 1)
  call error_check_(self, 'writeEMheader_:writeDatasetStringArray_:'//trim(line), hdferr)
end if

! close the dataname group, if present
if (PRESENT(dataname).eqv..TRUE.) call pop_(self)

! and close this group
call self%pop()

end subroutine writeEMheader_


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! utility functions and subroutines 
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
recursive function readfromTextfile_(self, filename, nlines) result(stringarray)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! read a text file and return it in a C_NULL_CHAR terminated string array 

use ISO_C_BINDING

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                                   :: self
character(fnlen),INTENT(IN)                                  :: filename
integer(kind=irg),INTENT(OUT)                                :: nlines
character(len=fnlen, KIND=c_char), allocatable               :: stringarray(:) 

integer(kind=irg)                                            :: i, j, dt
character(len=fnlen, KIND=c_char), DIMENSION(1)              :: line 

dt = 55
! read the file first to determine the number of lines
open(unit=dt,file=trim(filename),action='read',form='formatted',status='old')
nlines = 0
do
  read (dt,"(A)",end=10) line(1)
  nlines = nlines + 1
end do
10 close(unit=dt,status='keep')

! then re-read the file and store all the lines in the wdata array
allocate(stringarray(1:nlines))
open(unit=dt,file=trim(filename),action='read',form='formatted',status='old')
do i=1,nlines
! initialize the line to null characters before each read
  do j=1,fnlen
    line(1)(j:j) = char(0)
  end do
! read the line
  read (dt,"(A)") line(1)
! find the string length and put the next character equal to C_NULL_CHAR
  j = len(trim(line(1)))+1
! truncate a line if it has more than fnlen characters
  if (j.gt.fnlen) j = fnlen
  line(1)(j:j) = C_NULL_CHAR
! store the line in the array
  stringarray(i) = line(1)
end do
close(unit=dt,status='keep')

end function readfromTextfile_

!--------------------------------------------------------------------------
pure recursive function cstringify(strin) result(cstrout)
  !! author: MDG 
  !! version: 1.0 
  !! date: 01/09/20
  !!
  !! turn a fortran string into a null-terminated c-string 

use ISO_C_BINDING

IMPLICIT NONE

character(fnlen),INTENT(IN)                   :: strin
character(len=len_trim(strin)+1,kind=c_char)  :: cstrout

integer(kind=irg)                             :: slen, i

slen = len_trim(strin)
do i=1,slen
  cstrout(i:i) = strin(i:i)
end do
slen = slen+1
cstrout(slen:slen) = C_NULL_CHAR 

end function cstringify

!--------------------------------------------------------------------------
recursive subroutine writeNMLintegers_(self, io_int, intlist, n_int)
!DEC$ ATTRIBUTES DLLEXPORT :: HDF_writeNMLintegers
 !! author: MDG
 !! version: 1.0 
 !! date: 02/04/20
 !!
 !! write a series of integer namelist entries to an HDF file

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                :: self
integer(kind=irg),INTENT(IN)              :: n_int
integer(kind=irg),INTENT(IN)              :: io_int(n_int)
character(20),INTENT(IN)                  :: intlist(n_int)

integer(kind=irg)                         :: hdferr, i
character(fnlen)                          :: dataset
logical                                   :: g_exists, overwrite=.TRUE.

do i=1,n_int
  dataset = intlist(i)
  call H5Lexists_f(self%head%next%objectID,trim(dataset),g_exists, hdferr)
  if (g_exists) then 
    hdferr = self%writeDatasetInteger(dataset, io_int(i), overwrite)
   else
    hdferr = self%writeDatasetInteger(dataset, io_int(i))
  end if
  if (hdferr.ne.0) call error_check_(self,'writeNMLintegers: unable to create '//trim(intlist(i))//' dataset', hdferr)
end do

end subroutine writeNMLintegers_

!--------------------------------------------------------------------------
recursive subroutine writeNMLreals_(self, io_real, reallist, n_real)
 !! author: MDG
 !! version: 1.0 
 !! date: 02/04/20
 !!
 !! write a series of real namelist entries to an HDF file

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                :: self
integer(kind=irg),INTENT(IN)              :: n_real
real(kind=sgl),INTENT(IN)                 :: io_real(n_real)
character(20),INTENT(IN)                  :: reallist(n_real)

integer(kind=irg)                         :: hdferr, i
character(fnlen)                          :: dataset
logical                                   :: g_exists, overwrite=.TRUE.

do i=1,n_real
  dataset = reallist(i)
  call H5Lexists_f(self%head%next%objectID,trim(dataset),g_exists, hdferr)
  if (g_exists) then 
    hdferr = self%writeDatasetFloat(dataset, io_real(i), overwrite)
  else
    hdferr = self%writeDatasetFloat(dataset, io_real(i))
  end if
  if (hdferr.ne.0) call error_check_(self,'writeNMLreals: unable to create '//trim(reallist(i))//' dataset',hdferr)
end do

end subroutine writeNMLreals_

!--------------------------------------------------------------------------
recursive subroutine writeNMLdbles_(self, io_real, reallist, n_real)
 !! author: MDG
 !! version: 1.0 
 !! date: 02/04/20
 !!
 !! write a series of double precision namelist entries to an HDF file

IMPLICIT NONE

class(HDF_T),INTENT(INOUT)                :: self
integer(kind=irg),INTENT(IN)              :: n_real
real(kind=dbl),INTENT(IN)                 :: io_real(n_real)
character(20),INTENT(IN)                  :: reallist(n_real)

integer(kind=irg)                         :: hdferr, i
character(fnlen)                          :: dataset
logical                                   :: g_exists, overwrite=.TRUE.

do i=1,n_real
  dataset = reallist(i)
  call H5Lexists_f(self%head%next%objectID,trim(dataset),g_exists, hdferr)
  if (g_exists) then 
    hdferr = self%writeDatasetDouble(dataset, io_real(i), overwrite)
  else
    hdferr = self%writeDatasetDouble(dataset, io_real(i))
  end if
  if (hdferr.ne.0) call error_check_(self,'writeNMLdbles: unable to create '//trim(reallist(i))//' dataset', hdferr)
end do

end subroutine writeNMLdbles_

end module mod_HDFsupport
