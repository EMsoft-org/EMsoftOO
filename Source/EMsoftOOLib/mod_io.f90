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
!> @date 12/31/19 MDG 1.0 original
!--------------------------------------------------------------------------

module mod_io
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! Message and error handling routines.
  !!
  !! We try to eliminate calls to *write* and *read* as much as possible from
  !! all programs and replace them by the routines in this module (for any IO
  !! that involves the command line, not for data files)

use mod_global 
use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                          stdout=>output_unit, &
                                          stderr=>error_unit

IMPLICIT NONE 

private
public :: T_IOClass 


  type, public  :: T_IOClass
    private
      character(fnlen)  :: message 
       !! a simple string of length fnlen

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
      procedure, pass(self) :: printMessageSingle
      procedure, pass(self) :: printMessageMultiple

      procedure, pass(self), public :: printWarning

      generic, public :: ReadValue => ReadValueIntShort, ReadValueIntLong, ReadValueRealSingle, &
                                      ReadValueRealDouble, ReadValueString, ReadValueStringArray
      generic, public :: WriteValue => WriteValueIntShort, WriteValueIntLong, WriteValueIntLongLong, &
                                       WriteValueRealSingle, WriteValueRealDouble, WriteValueRealComplex, &
                                       WriteValueString
      generic, public :: printError => printShortError, printErrorStatus
      generic, public :: printMessage => printMessageSingle, printMessageMultiple

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
!> @brief initialize the IO class; 
!
!> @date  12/31/19 MDG 1.0 new function
!--------------------------------------------------------------------------
type(T_IOClass) function Message_constructor( m ) result(Message)
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! constructor for the IO Class 
  
IMPLICIT NONE

character(fnlen), INTENT(IN), OPTIONAL      :: m 
 !! input character string

if (present(m)) then 
  Message % message = trim(m)
else 
  Message % message = ''
end if 

end function Message_constructor

!--------------------------------------------------------------------------
!
! SUBROUTINE: printMessageSingle
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief dump a message to standard output
! 
!> @date 12/31/19 MDG 1.0 original
!--------------------------------------------------------------------------
subroutine printMessageSingle(self, mess, frm, advance, redirect)
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! Simple routine to print a string to standard output (default), with optional formatting
  !! instructions; for instance, if one wants an empty line before (frm='(/A)') or after (frm='(A/)') 
  !! the string.  Note that one can include the name of the optional variable in the subroutine
  !! call, as in:
  !!
  !! call self % printMessage('this is a string', frm='(//A//)' , redirect = 22)


IMPLICIT NONE

  class(T_IOClass),intent(inout)          :: self

  character(*),INTENT(IN)                 :: mess         
   !! message string
  character(*),OPTIONAL,INTENT(IN)        :: frm          
   !! optional formatting string
  character(*),OPTIONAL,INTENT(IN)        :: advance      
   !! optional keyword to omit linefeed character
  integer(kind=irg),OPTIONAL,INTENT(IN)   :: redirect     
   !! optional redirect to this unit (stdout by default)

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

end subroutine printMessageSingle

!--------------------------------------------------------------------------
!
! SUBROUTINE: printMessageMultiple
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief dump a message with multiple lines to standard output
!
!> @date 12/31/19 MDG 1.0 original
!--------------------------------------------------------------------------
subroutine printMessageMultiple(self, mess, frm, redirect)
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! Simple routine to print one or more strings to standard output (default), with optional formatting
  !! instructions; for instance, if one wants an empty line before (frm='(/A)') or after (frm='(A/)') 
  !! the string.  Note that one can include the name of the optional variable in the subroutine
  !! call, as in:
  !!
  !! call self % printMessage( (/'this is a string       ', &
  !!                             'and this is another one'/), redirect = 10)
  !! Note that *mess* is a string array, so all component strings MUST have the same length!
   
IMPLICIT NONE

  class(T_IOClass),intent(inout)          :: self

  character(*),INTENT(IN)                 :: mess(:)      
   !! message array of strings
  character(*),OPTIONAL,INTENT(IN)        :: frm          
   !! optional formatting string
  integer(kind=irg),OPTIONAL,INTENT(IN)   :: redirect     
   !! redirect to this unit

  integer(kind=irg)                       :: unit, ss(1), i 

  ss = shape(mess) 

  unit = stdout
  if (present(redirect)) unit = redirect 

! default format or not ?
  if (PRESENT(frm)) then
    do i=1,ss(1)
      write (unit,fmt=frm) trim(mess(i))
    end do
  else    ! default output format: a simple string
    do i=1,ss(1)
      write (unit,fmt="(A)") trim(mess(i))
    end do
  end if 

end subroutine printMessageMultiple

!--------------------------------------------------------------------------
!
! SUBROUTINE: printShortError
!
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief Write error message and abort program
!
!> @date   12/31/19 MDG 1.0 original
! ###################################################################
subroutine printShortError(self, s1, s2)
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! write an error message (routine_name: message) and abort the program

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*), INTENT(IN)  :: s1  
   !! first part of error message (routine name)
  character(*), INTENT(IN)  :: s2  
   !! second part of error message (brief explanation)

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
!> @date   12/31/19 MDG 1.0 original
! ###################################################################
subroutine printErrorStatus(self, s1, status, s2)
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! write an error message with status number and abort program

IMPLICIT NONE

  class(T_IOClass),intent(inout)      :: self

  character(*), INTENT(IN)            :: s1      
   !! first part of error message (routine name)
  integer(kind=irg),INTENT(IN)        :: status  
   !! error identifier
  character(*), INTENT(IN),OPTIONAL   :: s2(:)   
   !! optional second part of error message (brief explanation); can have multiple lines
 
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
!> @date   12/31/19 MDG 1.0 original
! ###################################################################
subroutine printWarning(self, s1, s2)
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! write an warning message, potentially multiple lines 

IMPLICIT NONE

  class(T_IOClass),intent(inout)      :: self

  character(*), INTENT(IN)            :: s1     
   !! first part of error message (routine name)
  character(*), INTENT(IN),OPTIONAL   :: s2(:)  
   !! second part of error message (brief explanation)

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
!> @brief read a string from standard input (stdin)
!
!> @date 12/31/19 MDG 1.0 new routine
! ###################################################################
subroutine ReadValueString(self, Qstring, rd_string, frm)
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! read a string from standard input (stdin)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*),INTENT(IN)                         :: Qstring
   !! user prompt string 
  character(*),INTENT(OUT)                        :: rd_string
   !! string to be read 
  character(*),INTENT(IN),OPTIONAL                :: frm
   !! optional formatting string

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
!> @brief read an array of strings from standard input (stdin)
!
!> @date 12/31/19 MDG 1.0 new routine
! ###################################################################
subroutine ReadValueStringArray(self, Qstring, rd_string, num, frm)
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! read an array of strings from standard input (stdin)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*),INTENT(IN)                         :: Qstring
   !! user prompt string 
  character(1),INTENT(OUT)                        :: rd_string(num)
   !! array of strings to be read
  integer(kind=irg),INTENT(IN)                    :: num
   !! number of strings to read 
  character(*),INTENT(IN),OPTIONAL                :: frm
   !! optional formatting string

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
!> @date 12/31/19 MDG 1.0 new routine
! ###################################################################
subroutine ReadValueIntShort(self, Qstring, rd_int, num)
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! read one or more short integers from standard input (stdin)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*), INTENT(IN)                        :: Qstring
   !! user prompt string 
  integer(kind=ish),INTENT(OUT)                   :: rd_int(*)
   !! output array of short integers
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
   !! number of integers to read

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
!> @date 12/31/19 MDG 1.0 new routine
! ###################################################################
subroutine ReadValueIntLong(self, Qstring, rd_int, num)
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! read one or more long integers from standard input (stdin)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*), INTENT(IN)                        :: Qstring
   !! user prompt string
  integer(kind=irg),INTENT(OUT)                   :: rd_int(*)
   !! array to hold integers 
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
   !! number of integers to read

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
!> @date 12/31/19 MDG 1.0 new routine
! ###################################################################
subroutine ReadValueRealSingle(self, Qstring, rd_real, num)
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! read one or more 4-byte reals from standard input (stdin)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*), INTENT(IN)                        :: Qstring
   !! user prompt string 
  real(kind=sgl),INTENT(OUT)                      :: rd_real(*)
   !! array to hold single precision reals
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
   !! number of reals to read

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
!> @date 12/31/19 MDG 1.0 new routine
! ###################################################################
subroutine ReadValueRealDouble(self, Qstring, rd_real, num)
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! read one or more 8-byte reals from standard input (stdin)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*), INTENT(IN)                        :: Qstring
   !! user prompt string
  real(kind=dbl),INTENT(OUT)                      :: rd_real(*)
   !! array to hold double precision reals
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
   !! number of doubles to read 

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
!> @date 12/31/19 MDG 1.0 new routine
! ###################################################################
subroutine WriteValueString(self, Qstring, out_string, frm, advance, redirect)
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! write a string to standard input (stdin)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*),INTENT(IN)                         :: Qstring 
   !! comment string
  character(*),INTENT(IN)                         :: out_string
   !! output string
  character(*),INTENT(IN),OPTIONAL                :: frm
   !! optional formatting string
  character(*),INTENT(IN),OPTIONAL                :: advance
   !! optional hold on linefeed
  integer(kind=irg),INTENT(IN),OPTIONAL           :: redirect
   !! optional redirect to a different output unit

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
!> @date 12/31/19 MDG 1.0 new routine
! ###################################################################
subroutine WriteValueIntShort(self, Qstring, out_int, num, frm, advance, redirect)
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! write one or more short integers to output (stdout or redirect)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self
 
  character(*), INTENT(IN)                        :: Qstring
   !! comment string
  integer(kind=ish),INTENT(IN)                    :: out_int(*)
   !! one or more output short integers
  character(*),INTENT(IN),OPTIONAL                :: frm
   !! optional formatting string
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
   !! optional number of integers to write 
  character(*),INTENT(IN),OPTIONAL                :: advance
   !! optional hold on linefeed
  integer(kind=irg),INTENT(IN),OPTIONAL           :: redirect
   !! optional redirect to other output unit 

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
!> @date 12/31/19 MDG 1.0 new routine
! ###################################################################
subroutine WriteValueIntLong(self, Qstring, out_int, num, frm, advance, redirect)
  !! author: MDG 
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! write one or more regular integers to output (stdout or redirect)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self

  character(*), INTENT(IN)                        :: Qstring
   !! comment string
  integer(kind=irg),INTENT(IN)                    :: out_int(*)
   !! one or more output short integers
  character(*),INTENT(IN),OPTIONAL                :: frm
   !! optional formatting string
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
   !! optional number of integers to write 
  character(*),INTENT(IN),OPTIONAL                :: advance
   !! optional hold on linefeed
  integer(kind=irg),INTENT(IN),OPTIONAL           :: redirect
   !! optional redirect to other output unit 

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
!> @author Saransh Singh, Carnegie Mellon University
!
!> @brief write one or more 8-byte integers
!
!> @date 12/31/19 MDG 1.0 new routine
! ###################################################################
subroutine WriteValueIntLongLong(self, Qstring, out_int, num, frm, advance, redirect)
  !! author: Saransh Singh, revised by MDG
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! write one or more 8-byte integers to output (stdout or redirect)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self
 
  character(*), INTENT(IN)                        :: Qstring
   !! comment string
  integer(kind=ill),INTENT(IN)                    :: out_int(*)
   !! one or more output short integers
  character(*),INTENT(IN),OPTIONAL                :: frm
   !! optional formatting string
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
   !! optional number of integers to write 
  character(*),INTENT(IN),OPTIONAL                :: advance
   !! optional hold on linefeed
  integer(kind=irg),INTENT(IN),OPTIONAL           :: redirect
   !! optional redirect to other output unit 

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
!> @date 12/31/19 MDG 1.0 new routine
! ###################################################################
subroutine WriteValueRealSingle(self, Qstring, out_real, num, frm, advance, redirect)
  !! author: MDG
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! write one or more single precision reals to output (stdout or redirect)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self
 
  character(*), INTENT(IN)                        :: Qstring
   !! comment string
  real(kind=sgl),INTENT(IN)                       :: out_real(*)
   !! one or more output single precision reals
  character(*),INTENT(IN),OPTIONAL                :: frm
   !! optional formatting string
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
   !! optional number of reals to write 
  character(*),INTENT(IN),OPTIONAL                :: advance
   !! optional hold on linefeed
  integer(kind=irg),INTENT(IN),OPTIONAL           :: redirect
   !! optional redirect to other output unit 

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
!> @date 12/31/19 MDG 1.0 new routine
! ###################################################################
subroutine WriteValueRealDouble(self, Qstring, out_real, num, frm, advance, redirect)
  !! author: MDG
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! write one or more double precision reals to output (stdout or redirect)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self
 
  character(*), INTENT(IN)                        :: Qstring
   !! comment string
  real(kind=dbl),INTENT(IN)                       :: out_real(*)
   !! one or more output double precision reals
  character(*),INTENT(IN),OPTIONAL                :: frm
   !! optional formatting string
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
   !! optional number of reals to write 
  character(*),INTENT(IN),OPTIONAL                :: advance
   !! optional hold on linefeed
  integer(kind=irg),INTENT(IN),OPTIONAL           :: redirect
   !! optional redirect to other output unit 

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
!> @date 12/31/19 MDG 1.0 new routine
! ###################################################################
subroutine WriteValueRealComplex(self, Qstring, out_cmplx, num, frm, advance, redirect)
  !! author: MDG
  !! version: 1.0 
  !! date: 12/31/19
  !!
  !! write one or more single precision complex numbers to output (stdout or redirect)

IMPLICIT NONE

  class(T_IOClass),intent(inout) :: self
 
  character(*), INTENT(IN)                        :: Qstring
   !! comment string
  complex(kind=sgl),INTENT(IN)                    :: out_cmplx(*)
   !! one or more output single precision complex numbers
  character(*),INTENT(IN),OPTIONAL                :: frm
   !! optional formatting string
  integer(kind=irg),INTENT(IN),OPTIONAL           :: num
   !! optional number of complex numbers to write 
  character(*),INTENT(IN),OPTIONAL                :: advance
   !! optional hold on linefeed
  integer(kind=irg),INTENT(IN),OPTIONAL           :: redirect
   !! optional redirect to other output unit 

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

