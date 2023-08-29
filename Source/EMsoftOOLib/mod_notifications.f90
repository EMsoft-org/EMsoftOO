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

module mod_notifications
  !! author: MDG
  !! version: 1.0
  !! date: 01/26/20
  !!
  !! notifications module contains routines to update the user on program status
  !!
  !! The user can set the "notifications" parameter in the EMsoftConfig.json file
  !! to "Off", "Email", or "Slack"; each program will also have a "EMNotify" parameter which
  !! can be turned "On" of "Off".  The EMsoftConfig.json file also has the "EMSlackWebHookURL"
  !! parameter which the user needs to set to the correct URL for the "EMSlackChannel" to which
  !! the message should be POSTed.  If the option selected is "Email", then the UserEmail
  !! parameter will be used to send the message via regular email, assuming that email sending
  !! has been set up properly.

use mod_kinds
use mod_global
use mod_platformsupport

IMPLICIT NONE

contains

!--------------------------------------------------------------------------
recursive function PostSlackMessage(EMsoft, MessageLines, NumLines, MessageTitle) result(status)
 !DEC$ ATTRIBUTES DLLEXPORT :: PostSlackMessage
 !! author: MDG
  !! version: 1.0
  !! date: 01/26/20
  !!
  !! function to POST a message to a slack channel; this routine requires that the 'curl' program
  !! be installed.  This still needs to be implemented for Windows 10...

use mod_EMsoft
use mod_io

IMPLICIT NONE

type(EMsoft_T),INTENT(INOUT) :: EMsoft
integer(kind=irg),INTENT(IN) :: NumLines
character(fnlen),INTENT(IN)  :: MessageLines(NumLines)
character(fnlen),INTENT(IN)  :: MessageTitle
integer(kind=irg)            :: status

type(IO_T)                   :: Message
character(fnlen)             :: EMSlackWebHookURL, EMSlackChannel
character(4096)              :: JSONMessage, cmd
integer(kind=irg)            :: j

if (EMsoft%getConfigParameter('EMsoftplatform').ne.'Windows') then
  status = 0

  ! first get the EMSlackWebHookURL parameter
  EMSlackWebHookURL = EMsoft%getConfigParameter('SlackWebHookURL')
  if (len(trim(EMSlackWebHookURL)).eq.0) status = -20

  ! get the Slack Channel name
  EMSlackChannel = EMsoft%getConfigParameter('SlackChannel')
  if (len(trim(EMSlackChannel)).eq.0) status = -21

  if (status.eq.0) then
  ! put the json message together
     JSONMessage = '{"channel": "'
     JSONMessage = trim(JSONMessage)//trim(EMSlackChannel)//'",'
     JSONMessage = trim(JSONMessage)//' "username": "'//trim(MessageTitle)//'",'
     JSONMessage = trim(JSONMessage)//' "icon_url": "http://muri.materials.cmu.edu/wp-content/uploads/2017/08/128x128.png",'
     JSONMessage = trim(JSONMessage)//' "text": "'

     do j=1,NumLines-1
       JSONMessage = trim(JSONMessage)//trim(MessageLines(j))//'\n'
     end do
     JSONMessage = trim(JSONMessage)//trim(MessageLines(NumLines))//'"}'

   ! and generate the curl command string
     cmd = 'curl -s -d ''payload='//trim(JSONMessage)//''' '//trim(EMSlackWebHookURL)//' >/dev/null'
     status = system_system(cmd)
  end if
else
  call Message%printMessage('This option is not yet implemented on Windows 10...')
end if

end function PostSlackMessage

!--------------------------------------------------------------------------
recursive function PostEmailMessage(EMsoft, MessageLines, NumLines, MessageTitle) result(status)
!DEC$ ATTRIBUTES DLLEXPORT :: PostEmailMessage
  !! author: MDG
  !! version: 1.0
  !! date: 01/26/20
  !!
  !! function to POST a message via email; requires a functional mailing system ...

use mod_EMsoft

IMPLICIT NONE

type(EMsoft_T),INTENT(INOUT) :: EMsoft
integer(kind=irg),INTENT(IN) :: NumLines
character(fnlen),INTENT(IN)  :: MessageLines(NumLines)
character(fnlen),INTENT(IN)  :: MessageTitle
integer(kind=irg)            :: status

character(fnlen)             :: UserEmail
character(4096)              :: EmailMessage, cmd
integer(kind=irg)            :: j

status = 0

! first get the user's email address
UserEmail = EMsoft%getConfigParameter('Useremail')

if (len(trim(UserEmail)).eq.0) then
   status = -10
else
! put the email message together
   EmailMessage = '\nEMsoft automated email message service\n\n'

   do j=1,NumLines-1
      EmailMessage = trim(EmailMessage)//trim(MessageLines(j))//'\n'
   end do
   EmailMessage = trim(EmailMessage)//trim(MessageLines(NumLines))

! and generate the email send command string
   cmd = 'echo "'//trim(EmailMessage)//'" | mail -s "Message from '//trim(MessageTitle)//'" '//trim(UserEmail)
   status = system_system(cmd)
end if

end function PostEmailMessage

!--------------------------------------------------------------------------
recursive function PostMessage(EMsoft, MessageLines, NumLines, MessageTitle) result(status)
!DEC$ ATTRIBUTES DLLEXPORT :: PostMessage
  !! author: MDG
  !! version: 1.0
  !! date: 01/26/20
  !!
  !! function to POST a message; function decides whether to use email or slack

use mod_EMsoft
use mod_io

IMPLICIT NONE

type(EMsoft_T),INTENT(INOUT) :: EMsoft
integer(kind=irg),INTENT(IN) :: NumLines
character(fnlen),INTENT(IN)  :: MessageLines(NumLines)
character(fnlen),INTENT(IN)  :: MessageTitle
integer(kind=irg)            :: status

type(IO_T)                   :: Message
character(fnlen)             :: notifymode, platform
integer(kind=irg)            :: ierr

platform = EMsoft%getConfigParameter('EMsoftplatform')

if (trim(platform).ne.'Windows') then

! get the notification mode (should be 'Email' or 'Slack')
   notifymode = EMsoft%getConfigParameter('Notify')

   if (notifymode.eq.'Email') then
     ierr = PostEmailMessage(EMsoft, MessageLines, NumLines, MessageTitle)
   else
     ierr = PostSlackMessage(EMsoft, MessageLines, NumLines, MessageTitle)
   end if

! do some error handling based on the value of ierr
   if (ierr.ne.0) then
     select case (ierr)
       case(-10)
        call Message%printMessage('PostMessage Warning: User email address not set in EMsoftConfig.json configuration file')
       case(-20)
        call Message%printMessage('PostMessage Warning: SlackWebHookURL not set in EMsoftConfig.json configuration file')
       case(-21)
        call Message%printMessage('PostMessage Warning: SlackChannel not set in EMsoftConfig.json configuration file')
       case default
     end select
     call Message%printMessage('Message not posted via '//trim(notifymode))
   end if

else
   call Message%printMessage('PostMessage Warning: notifications are disabled on Windows at this time')
end if

end function PostMessage

end module mod_notifications
