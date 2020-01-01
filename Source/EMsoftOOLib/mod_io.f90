! ###################################################################
! Copyright (c) 2013-2019, Marc De Graef Research Group/Carnegie Mellon University
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

!--------------------------------------------------------------------------
! EMsoft:mod_io.f90
!--------------------------------------------------------------------------
!
! MODULE: mod_io
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief message and error handling routines
!
!> @details  Ideally, this should be the only module that has explicit write statements
!> in it.  
!
!> @date 12/31/19 MDG 1.0 original
!--------------------------------------------------------------------------

module mod_io

use mod_global 
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit

IMPLICIT NONE 

private
public :: T_IOClass 


  type, public  :: T_IOClass
    private
      character(fnlen)  :: origin 
      character(fnlen)  :: message 

    contains
    private

      procedure, pass(self) :: ReadValueIntShort
      procedure, pass(self) :: ReadValueIntLong
      procedure, pass(self) :: ReadValueRealSingle
      procedure, pass(self) :: ReadValueRealDouble
      procedure, pass(self) :: ReadValueString
      procedure, pass(self) :: ReadValueStringArray

      procedure, pass(self) :: WriteValueIntShort
      procedure, pass(self) :: WriteValueIntLong
      procedure, pass(self) :: WriteValueIntLongLong
      procedure, pass(self) :: WriteValueRealSingle
      procedure, pass(self) :: WriteValueRealDouble
      procedure, pass(self) :: WriteValueRealComplex
      procedure, pass(self) :: WriteValueString

      procedure, pass(self) :: printShortError
      procedure, pass(self) :: printErrorStatus

      procedure, pass(self), public :: printWarning
      procedure, pass(self), public :: printMessage

      generic, public :: ReadValue => ReadValueIntShort, ReadValueIntLong, ReadValueRealSingle, &
                                      ReadValueRealDouble, ReadValueString, ReadValueStringArray
      generic, public :: WriteValue => WriteValueIntShort, WriteValueIntLong, WriteValueIntLongLong, &
                                       WriteValueRealSingle, WriteValueRealDouble, WriteValueRealComplex, &
                                       WriteValueString
      generic, public :: printError => printShortError, printErrorStatus

  end type T_IOClass

  ! the constructor routine for this class 
  interface T_IOClass
    module procedure :: Message_constructor
  end interface T_IOClass

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! we begin with the functions/subroutines that are public in this class
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!
! FUNCTION: Message_constructor
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief initialize the IO class; this doesn't really do anything yet 
!
!> @date  12/31/19 MDG 1.0 new function
!--------------------------------------------------------------------------
type(T_IOClass) function Message_constructor() result(Message)

IMPLICIT NONE

end function Message_constructor

!--------------------------------------------------------------------------
!
! SUBROUTINE: printMessage
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief dump a message to standard output
!
!> @details Simple routine to print a string on the standard output, with optional formatting
!> instructions, for instance if one wants an empty line before (frm='(/A)') or after (frm='(A/)') 
!> the string.  Note that one can include the name of the optional variable in the subroutine
!> call, as in:
!> call Message('this is a string', frm='(//A//)' , stdout = 22)
!> this makes it clear that frm and stdout are optional variables.
! 
!> @param mess message string
!> @param frm optional string formatting command
!> @param advance  print a new line character ?
!
!> @date 12/31/19 MDG 1.0 original
!--------------------------------------------------------------------------
subroutine printMessage(self, mess, frm, advance, redirect)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*),INTENT(IN)                 :: mess         !< message string
  character(*),OPTIONAL,INTENT(IN)        :: frm          !< optional formatting string
  character(*),OPTIONAL,INTENT(IN)        :: advance      !< optional formatting string
  integer(kind=irg),OPTIONAL,INTENT(IN)   :: redirect     !< redirect to this unit

  integer(kind=irg)                       :: unit 

  unit = stdout
  if (present(redirect)) unit = redirect 

! default format or not ?
  if (PRESENT(frm)) then
   if (present(advance)) then
     write (unit,fmt=frm,advance="no") trim(mess)
   else 
     write (unit,fmt=frm) trim(mess)
   end if
  else    ! default output format: a simple string
   write (unit,fmt="(A)") trim(mess)
  end if 

end subroutine printMessage

!--------------------------------------------------------------------------
!
! SUBROUTINE: printShortError
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Write error message and abort program
!
!> @param s1 routine name string
!> @param s2 explanation string
!
!> @date   12/31/19 MDG 1.0 original
! ###################################################################
subroutine printShortError(self, s1, s2)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*), INTENT(IN)  :: s1  !< first part of error message (routine name)
  character(*), INTENT(IN)  :: s2  !< second part of error message (brief explanation)

  call self % printMessage(' ----> Fatal error in routine '//s1//': '//s2, frm='(//A/)', redirect=stderr) 
  stop '  Progam ended abnormally'

end subroutine printShortError

!--------------------------------------------------------------------------
!
! SUBROUTINE: printErrorStatus
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Write error message with status number and abort program
!
!> @param s1 string
!> @param s2 optional string 
!
!> @date   12/31/19 MDG 1.0 original
! ###################################################################
subroutine printErrorStatus(self, s1, status, s2)

IMPLICIT NONE

  class(T_IOClass),intent(inout)      :: self

  character(*), INTENT(IN)            :: s1      !< first part of error message (routine name)
  integer(kind=irg),INTENT(IN)        :: status  !< error identifier
  character(*), INTENT(IN),OPTIONAL   :: s2(:)   !< second part of error message (brief explanation)

  integer(kind=irg)                   :: io_int(1), ss2(1), i

  ss2 = shape(s2)

  call self % printMessage(' EMsoft error encountered:', frm='(//A)', redirect=stderr) 
  call self % printMessage('  '//s1, frm='(A)', redirect=stderr) 
  if (present(s2)) then 
    do i=1,ss2(1) 
      call self % printMessage('  '//s2(i), frm='(A)', redirect=stderr)
    end do 
  end if 
  io_int(1) = status
  call self % WriteValue(' Error Status :',io_int,1,"(I7/)", redirect=stderr)
  stop '  Progam ended abnormally'

end subroutine printErrorStatus


!--------------------------------------------------------------------------
!
! SUBROUTINE: printWarning
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Write warning message
!
!> @param s1 routine name string
!> @param s2 explanation string
!> @param stdout optional output unit identifier
!
!> @date   12/31/19 MDG 1.0 original
! ###################################################################
subroutine printWarning(self, s1, s2)

IMPLICIT NONE

  class(T_IOClass),intent(inout)      :: self

  character(*), INTENT(IN)            :: s1     !< first part of error message (routine name)
  character(*), INTENT(IN),OPTIONAL   :: s2(:)  !< second part of error message (brief explanation)

  integer(kind=irg)                   :: ss2(1), i

  ss2 = shape(s2)

  call self % printMessage(' EMsoft warning encountered:', frm='(//A)', redirect=stderr) 
  call self % printMessage('  '//s1, frm='(A)', redirect=stderr) 
  if (present(s2)) then 
    do i=1,ss2(1) 
      call self % printMessage('  '//s2(i), frm='(A)', redirect=stderr)
    end do 
  end if 

end subroutine printWarning


! ###################################################################
! reading routines
! ###################################################################


! ###################################################################
! 
!  subroutine ReadValueString   
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read a string from standard input (unit = 5)
!
!> @param Qstring question string
!> @param rd_string string to be read
!> @param frm optional format string

!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
subroutine ReadValueString(self, Qstring, rd_string, frm)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*),INTENT(IN)                         :: Qstring
  character(*),INTENT(OUT)                        :: rd_string
  character(*),INTENT(IN),OPTIONAL                :: frm

  call self % printMessage(Qstring, frm = "(' ',A,' ')", advance="no")

  if (PRESENT(frm)) then
    read (stdin, fmt=frm) rd_string
  else
    read (stdin,*) rd_string
  end if

end subroutine ReadValueString

! ###################################################################
! 
!  subroutine ReadValueStringArray   
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read an array of strings from standard input (unit = 5)
!
!> @param Qstring question string
!> @param rd_string string to be read
!> @param num number of strings in array
!> @param frm optional format string

!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
subroutine ReadValueStringArray(self, Qstring, rd_string, num, frm)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*),INTENT(IN)                         :: Qstring
  character(1),INTENT(OUT)                        :: rd_string(num)
  integer(kind=irg),INTENT(IN)                    :: num
  character(*),INTENT(IN),OPTIONAL                :: frm

  integer(kind=irg)                               :: i

  call self % printMessage(Qstring, frm = "(' ',A,' ')",advance="no")

  if (PRESENT(frm)) then 
    do i=1,num
      read (stdin, fmt=frm) rd_string(i)
    end do
  else  
    do i=1,num
      read (stdin,*) rd_string(i)
    end do
  end if 

end subroutine ReadValueStringArray

! ###################################################################
! 
!  subroutine ReadValueIntShort   
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read one or more short integers
!
!> @param Qstring question string
!> @param rd_int integer to be read
!> @param num optional number of integers to be read

!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
subroutine ReadValueIntShort(self, Qstring, rd_int, num)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*), INTENT(IN)                        :: Qstring
  integer(kind=ish),INTENT(OUT)                   :: rd_int(*)
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num

  integer(kind=irg)                               :: i

  call self % printMessage(Qstring, frm = "(' ',A,' ')",advance="no")

  ! one or more than one values expected ?
  if (PRESENT(num)) then
    read (stdin,*) (rd_int(i),i=1,num)
  else
    read (stdin,*) rd_int(1)
  end if
  
end subroutine ReadValueIntShort

! ###################################################################
! 
!  subroutine ReadValueIntLong 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read one or more regular (4-byte) integers
!
!> @param Qstring question string
!> @param rd_int integer to be read
!> @param num optional number of integers to be read

!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
subroutine ReadValueIntLong(self, Qstring, rd_int, num)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*), INTENT(IN)                        :: Qstring
  integer(kind=irg),INTENT(OUT)                   :: rd_int(*)
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num

  integer(kind=irg)                               :: i

  call self % printMessage(Qstring, frm = "(' ',A,' ')",advance="no")

  ! one or more than one values expected ?
  if (PRESENT(num)) then
    read (stdin,*) (rd_int(i),i=1,num)
  else
    read (stdin,*) rd_int(1)
  end if
  
end subroutine ReadValueIntLong

! ###################################################################
! 
!  subroutine ReadValueRealSingle 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read one or more regular (4-byte) reals
!
!> @param Qstring question string
!> @param rd_real integer to be read
!> @param num optional number of integers to be read

!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
subroutine ReadValueRealSingle(self, Qstring, rd_real, num)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*), INTENT(IN)                        :: Qstring
  real(kind=sgl),INTENT(OUT)                      :: rd_real(*)
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num

  integer(kind=irg)                               :: i

  call self % printMessage(Qstring, frm = "(' ',A,' ')",advance="no")

  ! one or more than one values expected ?
  if (PRESENT(num)) then
    read (stdin,*) (rd_real(i),i=1,num)
  else
    read (stdin,*) rd_real(1)
  end if

end subroutine ReadValueRealSingle

! ###################################################################
! 
!  subroutine ReadValueRealDouble 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief read one or more regular (4-byte) reals
!
!> @param Qstring question string
!> @param rd_real integer to be read
!> @param num optional number of integers to be read
!> @param stdout optional output unit identifier

!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
subroutine ReadValueRealDouble(self, Qstring, rd_real, num)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*), INTENT(IN)                        :: Qstring
  real(kind=dbl),INTENT(OUT)                      :: rd_real(*)
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num

  integer(kind=irg)                               :: i

  call self % printMessage(Qstring, frm = "(' ',A,' ')",advance="no")

  ! one or more than one values expected ?
  if (PRESENT(num)) then
    read (stdin,*) (rd_real(i),i=1,num)
  else
    read (stdin,*) rd_real(1)
  end if

end subroutine ReadValueRealDouble


! ###################################################################
! writing routines
! ###################################################################

! ###################################################################
! 
!  subroutine WriteValueString 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write a string
!
!> @param Qstring question string
!> @param out_string output string
!> @param frm optional formatting argument
!
!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
subroutine WriteValueString(self, Qstring, out_string, frm, advance, redirect)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*),INTENT(IN)                         :: Qstring 
  character(*),INTENT(IN)                         :: out_string
  character(*),INTENT(IN),OPTIONAL                :: frm
  character(*),INTENT(IN),OPTIONAL                :: advance
  integer(kind=irg),INTENT(IN),OPTIONAL           :: redirect

  ! send Qstring to the output only if it is non-zero length
  if (len(Qstring).ne.0) then
    if (present(redirect)) then 
      call self % printMessage(Qstring, frm = "(A)",advance="no",redirect=redirect)
    else
      call self % printMessage(Qstring, frm = "(A)",advance="no")
    end if 
  end if


  if (PRESENT(frm)) then 
    if (present(advance)) then 
      if (present(redirect)) then 
        call self % printMessage(out_string, frm = frm, advance="no", redirect=redirect)
      else
        call self % printMessage(out_string, frm = frm, advance="no")
      end if
    else 
      if (present(redirect)) then 
        call self % printMessage(out_string, frm = frm, redirect=redirect)
      else
        call self % printMessage(out_string, frm = frm)
      end if 
    end if
  else
    if (present(redirect)) then 
      call self % printMessage(out_string, frm = "(A)", redirect=redirect)
    else
      call self % printMessage(out_string, frm = "(A)")
    end if
  end if

end subroutine WriteValueString

! ###################################################################
! 
!  subroutine WriteValueIntShort 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write one or more short integers
!
!> @param Qstring question string
!> @param out_int output string
!> @param num optional number of integers
!> @param frm optional formatting argument
!
!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
subroutine WriteValueIntShort(self, Qstring, out_int, num, frm, advance, redirect)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self
 
  character(*), INTENT(IN)                        :: Qstring
  integer(kind=ish),INTENT(IN)                    :: out_int(*)
  character(*),INTENT(IN),OPTIONAL                :: frm
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
  character(*),INTENT(IN),OPTIONAL                :: advance
  integer(kind=irg),INTENT(IN),OPTIONAL           :: redirect

  integer(kind=irg)                               :: i, unit  

  unit = stdout
  if (present(redirect)) unit = redirect

  if (len(Qstring).ne.0) then 
    if (present(redirect)) then 
      call self % printMessage(Qstring, frm = "(A)",advance="no",redirect=redirect)
    else
      call self % printMessage(Qstring, frm = "(A)",advance="no")
    end if 
  end if 

  ! one or more than one values expected ?
  if (PRESENT(num)) then
   if (PRESENT(frm)) then
     if (present(advance)) then 
      write (unit, fmt=frm, advance="no") (out_int(i),i=1,num)
    else
      write (unit, fmt=frm) (out_int(i),i=1,num)
    end if 
   else
    write (unit,*) (out_int(i),i=1,num)
   end if
  else
   if (PRESENT(frm)) then
    if (present(advance)) then 
      write (unit, fmt=frm, advance="no") out_int(1)
    else
      write (unit, fmt=frm) out_int(1)
    end if
   else
    write (unit,*) out_int(1)
   end if
  end if

end subroutine WriteValueIntShort

! ###################################################################
! 
!  subroutine WriteValueIntLong 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write one or more 4-byte integers
!
!> @param Qstring question string
!> @param out_int output string
!> @param num optional number of integers
!> @param frm optional formatting argument
!
!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
subroutine WriteValueIntLong(self, Qstring, out_int, num, frm, advance, redirect)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self
 
  character(*), INTENT(IN)                        :: Qstring
  integer(kind=irg),INTENT(IN)                    :: out_int(*)
  character(*),INTENT(IN),OPTIONAL                :: frm
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
  character(*),INTENT(IN),OPTIONAL                :: advance
  integer(kind=irg),INTENT(IN),OPTIONAL           :: redirect

  integer(kind=irg)                               :: i, unit  

  unit = stdout
  if (present(redirect)) unit = redirect

  if (len(Qstring).ne.0) then 
    if (present(redirect)) then 
      call self % printMessage(Qstring, frm = "(A)",advance="no",redirect=redirect)
    else
      call self % printMessage(Qstring, frm = "(A)",advance="no")
    end if 
  end if 

  ! one or more than one values expected ?
  if (PRESENT(num)) then
   if (PRESENT(frm)) then
    if (present(advance)) then 
      write (unit, fmt=frm, advance="no") (out_int(i),i=1,num)
    else
      write (unit, fmt=frm) (out_int(i),i=1,num)
    end if 
   else
    write (unit,*) (out_int(i),i=1,num)
   end if
  else
   if (PRESENT(frm)) then
    if (present(advance)) then 
      write (unit, fmt=frm, advance="no") out_int(1)
    else
      write (unit, fmt=frm) out_int(1)
    end if
   else
    write (unit,*) out_int(1)
   end if
  end if

end subroutine WriteValueIntLong

! ###################################################################
! 
!  subroutine WriteValueIntLongLong
!
!> @author Saransh, Carnegie Mellon University
!
!> @brief write one or more 8-byte integers
!
!> @param Qstring question string
!> @param out_int output string
!> @param num optional number of integers
!> @param frm optional formatting argument
!
!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
subroutine WriteValueIntLongLong(self, Qstring, out_int, num, frm, advance, redirect)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self
 
  character(*), INTENT(IN)                        :: Qstring
  integer(kind=ill),INTENT(IN)                    :: out_int(*)
  character(*),INTENT(IN),OPTIONAL                :: frm
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
  character(*),INTENT(IN),OPTIONAL                :: advance
  integer(kind=irg),INTENT(IN),OPTIONAL           :: redirect

  integer(kind=irg)                               :: i, unit  

  unit = stdout
  if (present(redirect)) unit = redirect

  if (len(Qstring).ne.0) then 
    if (present(redirect)) then 
      call self % printMessage(Qstring, frm = "(A)",advance="no",redirect=redirect)
    else
      call self % printMessage(Qstring, frm = "(A)",advance="no")
    end if 
  end if 

  ! one or more than one values expected ?
  if (PRESENT(num)) then
   if (PRESENT(frm)) then
    if (present(advance)) then 
      write (unit, fmt=frm, advance="no") (out_int(i),i=1,num)
    else
      write (unit, fmt=frm) (out_int(i),i=1,num)
    end if 
   else
    write (unit,*) (out_int(i),i=1,num)
   end if
  else
   if (PRESENT(frm)) then
      if (present(advance)) then 
      write (unit, fmt=frm, advance="no") out_int(1)
    else
      write (unit, fmt=frm) out_int(1)
    end if
   else
    write (unit,*) out_int(1)
   end if
  end if

end subroutine WriteValueIntLongLong


! ###################################################################
! 
!  subroutine WriteValueRealSingle 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write one or more single precision reals
!
!> @param Qstring question string
!> @param out_real output string
!> @param num optional number of integers
!> @param frm optional formatting argument
!> @param stdout optional output unit identifier
!
!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
subroutine WriteValueRealSingle(self, Qstring, out_real, num, frm, advance, redirect)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self
 
  character(*), INTENT(IN)                        :: Qstring
  real(kind=sgl),INTENT(IN)                       :: out_real(*)
  character(*),INTENT(IN),OPTIONAL                :: frm
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
  character(*),INTENT(IN),OPTIONAL                :: advance
  integer(kind=irg),INTENT(IN),OPTIONAL           :: redirect

  integer(kind=irg)                               :: i, unit  

  unit = stdout
  if (present(redirect)) unit = redirect

  if (len(Qstring).ne.0) then 
    if (present(redirect)) then 
      call self % printMessage(Qstring, frm = "(A)",advance="no",redirect=redirect)
    else
      call self % printMessage(Qstring, frm = "(A)",advance="no")
    end if 
  end if 

  ! one or more than one values expected ?
  if (PRESENT(num)) then
   if (PRESENT(frm)) then
    if (present(advance)) then 
      write (unit, fmt=frm, advance="no") (out_real(i),i=1,num)
    else
      write (unit, fmt=frm) (out_real(i),i=1,num)
    end if  
   else
    write (unit,*) (out_real(i),i=1,num)
   end if
  else
   if (PRESENT(frm)) then
    if (present(advance)) then 
      write (unit, fmt=frm, advance="no") out_real(1)
    else
      write (unit, fmt=frm) out_real(1)
    end if
   else
    write (unit,*) out_real(1)
   end if
  end if

end subroutine WriteValueRealSingle



! ###################################################################
! 
!  subroutine WriteValueRealDouble 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write one or more double precision reals
!
!> @param Qstring question string
!> @param out_real output string
!> @param num optional number of integers
!> @param frm optional formatting argument
!
!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
subroutine WriteValueRealDouble(self, Qstring, out_real, num, frm, advance, redirect)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self
 
  character(*), INTENT(IN)                        :: Qstring
  real(kind=dbl),INTENT(IN)                       :: out_real(*)
  character(*),INTENT(IN),OPTIONAL                :: frm
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
  character(*),INTENT(IN),OPTIONAL                :: advance
  integer(kind=irg),INTENT(IN),OPTIONAL           :: redirect

  integer(kind=irg)                               :: i, unit  

  unit = stdout
  if (present(redirect)) unit = redirect

  if (len(Qstring).ne.0) then 
    if (present(redirect)) then 
      call self % printMessage(Qstring, frm = "(A)",advance="no",redirect=redirect)
    else
      call self % printMessage(Qstring, frm = "(A)",advance="no")
    end if 
  end if 

  ! one or more than one values expected ?
  if (PRESENT(num)) then
   if (PRESENT(frm)) then
      if (present(advance)) then 
      write (unit, fmt=frm, advance="no") (out_real(i),i=1,num)
    else
      write (unit, fmt=frm) (out_real(i),i=1,num)
    end if  
   else
    write (unit,*) (out_real(i),i=1,num)
   end if
  else
   if (PRESENT(frm)) then
    if (present(advance)) then 
      write (unit, fmt=frm, advance="no") out_real(1)
    else
      write (unit, fmt=frm) out_real(1)
    end if
   else
    write (unit,*) out_real(1)
   end if
  end if

end subroutine WriteValueRealDouble


! ###################################################################
! 
!  subroutine WriteValueRealComplex 
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief write one or more single precision complex numbers
!
!> @param Qstring question string
!> @param out_real output string
!> @param num optional number of integers
!> @param frm optional formatting argument
!
!> @date 03/19/13 MDG 1.0 new routine
!> @date 06/05/14 MDG 2.0 changed io handling
!> @date 03/29/18 MDG 4.1 removed stdout argument
! ###################################################################
subroutine WriteValueRealComplex(self, Qstring, out_cmplx, num, frm, advance, redirect)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self
 
  character(*), INTENT(IN)                        :: Qstring
  complex(kind=sgl),INTENT(IN)                    :: out_cmplx(*)
  character(*),INTENT(IN),OPTIONAL                :: frm
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
  character(*),INTENT(IN),OPTIONAL                :: advance
  integer(kind=irg),INTENT(IN),OPTIONAL           :: redirect

  integer(kind=irg)                               :: i, unit

   unit = stdout
  if (present(redirect)) unit = redirect

  if (len(Qstring).ne.0) then 
    if (present(redirect)) then 
      call self % printMessage(Qstring, frm = "(A)",advance="no",redirect=redirect)
    else
      call self % printMessage(Qstring, frm = "(A)",advance="no")
    end if 
  end if 

  ! one or more than one values expected ?
  if (PRESENT(num)) then
   if (PRESENT(frm)) then
      if (present(advance)) then 
      write (unit, fmt=frm, advance="no") (out_cmplx(i),i=1,num)
    else
      write (unit, fmt=frm) (out_cmplx(i),i=1,num)
    end if  
   else
    write (unit,*) (out_cmplx(i),i=1,num)
   end if
  else
   if (PRESENT(frm)) then
    if (present(advance)) then 
      write (unit, fmt=frm, advance="no") out_cmplx(1)
    else
      write (unit, fmt=frm) out_cmplx(1)
    end if
   else
    write (unit,*) out_cmplx(1)
   end if
  end if

end subroutine WriteValueRealComplex



end module mod_io

