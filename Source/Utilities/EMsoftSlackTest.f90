! ###################################################################
! Copyright (c) 2013-2023, Marc De Graef Research Group/Carnegie Mellon University
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

program EMsoftSlackTest
  !! author: MDG
  !! version: 1.0 
  !! date: 01/26/20
  !!
  !! EMsoftSlackTest sends a simple test message to a Slack channel

use mod_kinds
use mod_global
use mod_EMsoft 
use mod_io
use mod_notifications

IMPLICIT NONE

character(fnlen)              :: progname = 'EMsoftSlackTest'
character(fnlen)              :: progdesc = 'Sends a simple test message to the currently configured Slack channel'
type(EMsoft_T)                :: EMsoft 

type(IO_T)                    :: Message
character(fnlen),ALLOCATABLE  :: MessageLines(:)
integer(kind=irg)             :: NumLines, i
character(fnlen)              :: MessageTitle, line
character(100)                :: c

EMsoft = EMsoft_T( progname, progdesc, tpl = (/ 922 /) )

if (trim(EMsoft%getConfigParameter('Notify')).ne.'Off') then
    NumLines = 1
    allocate(MessageLines(NumLines))
    call Message%printMessage(' Enter a test sentence (between single quotes):')
    call Message%ReadValue(' ---> ', line)

    call hostnm(c)
 
    MessageLines(1) = trim(line)
    MessageTitle = 'EMsoft on '//trim(c)
    i = PostMessage(EMsoft, MessageLines, NumLines, MessageTitle)
else
   call Message%printMessage('Notifications are turned off in your EMsoftConfig.json configuration file')
end if

end program EMsoftSlackTest
