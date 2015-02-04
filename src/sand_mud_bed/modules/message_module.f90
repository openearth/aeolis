module message_module
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011.                                     
!                                                                               
!  This program is free software: you can redistribute it and/or modify         
!  it under the terms of the GNU General Public License as published by         
!  the Free Software Foundation version 3.                                      
!                                                                               
!  This program is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU General Public License for more details.                                 
!                                                                               
!  You should have received a copy of the GNU General Public License            
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.        
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------
!  $Id: message_module.f90 7697 2012-11-16 14:10:17Z boer_aj $
!  $HeadURL: https://svn.oss.deltares.nl/repos/openearthtools/trunk/programs/SandMudBedModule/03_Fortran/example/example/modules/message_module.f90 $
!!--description-----------------------------------------------------------------
!
!    Function: - Keep track of a stack of messages
!
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------

private

!
! constants
!
integer, parameter, public :: message_len = 1024

!
! data types
!
public message_stack

!
! functions and subroutines
!
public message_module_info
public initstack
public isempty
public addmessage
public adderror
public addwarning
public getmessage
public writemessages
public clearstack

!
! A message stack is a chained list of messages; each message contains a single
! string and a pointer to the next message.
!
type message_type
    private
    character(message_len)       :: message
    type(message_type) , pointer :: other_messages
end type message_type

type message_stack
    private
    type(message_type), pointer  :: message_list
end type message_stack

contains
!
!
!
!==============================================================================
subroutine message_module_info(messages)
    !    Function: - Add info about this message module to the message stack
    !
    implicit none
    !
    ! arguments
    !
    type(message_stack) :: messages
    !
    ! body
    !
    call addmessage(messages,'$Id: message_module.f90 7697 2012-11-16 14:10:17Z boer_aj $')
    call addmessage(messages,'$URL: https://svn.oss.deltares.nl/')
end subroutine message_module_info
!
!
!
!==============================================================================
subroutine initstack(stack)
    !    Function: - Create a message stack
    !
    implicit none
    !
    ! Call variables
    !
    type(message_stack)          :: stack
    !
    ! Local variables
    !
    !
    ! body
    !
    nullify(stack%message_list)
end subroutine initstack
!
!
!
!==============================================================================
subroutine addmessage(stack,newmessage)
    !    Function: - Add a new message to the bottom of the message stack
    !
    implicit none
    !
    ! Call variables
    !
    character(*)        :: newmessage
    type(message_stack) :: stack
    !
    ! Local variables
    !
    type(message_type), pointer :: localstack
    !
    ! body
    !
    if (associated(stack%message_list)) then
       localstack => stack%message_list
       do while (associated(localstack%other_messages))
          localstack => localstack%other_messages
       end do
       allocate(localstack%other_messages)
       localstack => localstack%other_messages
       localstack%message = newmessage
       nullify(localstack%other_messages)
    else
       allocate(stack%message_list)
       localstack => stack%message_list
       localstack%message = newmessage
       nullify(localstack%other_messages)
    endif
end subroutine addmessage
!
!
!
!==============================================================================
subroutine adderror(stack,newmessage)
    !    Function: - Add a new error message to the bottom of the message stack
    !
    implicit none
    !
    ! Call variables
    !
    type(message_stack) :: stack
    character(*)        :: newmessage
    !
    ! Local variables
    !
    !
    ! body
    !
    newmessage = '*** ERROR ' // newmessage
    call addmessage(stack,newmessage)
end subroutine adderror
!
!
!
!==============================================================================
subroutine addwarning(stack,newmessage)
    !    Function: - Add a new warning message to the bottom of the message stack
    !
    implicit none
    !
    ! Call variables
    !
    type(message_stack) :: stack
    character(*)        :: newmessage
    !
    ! Local variables
    !
    !
    ! body
    !
    newmessage = '*** WARNING ' // newmessage
    call addmessage(stack,newmessage)
end subroutine addwarning
!
!
!
!==============================================================================
function isempty(stack)
    !    Function: - Checks whether there is any message on the stack
    !
    implicit none
    !
    ! Call variables
    !
    type(message_stack) :: stack
    logical             :: isempty
    !
    ! Local variables
    !
    !
    ! body
    !
    isempty = .not.associated(stack%message_list)
end function isempty
!
!
!
!==============================================================================
subroutine getmessage(stack,message)
    !    Function: - Get a message from the top of the message stack
    !
    implicit none
    !
    ! Call variables
    !
    type(message_stack)    :: stack
    character(message_len) :: message
    !
    ! Local variables
    !
    type(message_type), pointer :: localstack
    !
    ! body
    !
    if (associated(stack%message_list)) then
       localstack => stack%message_list
       message = localstack%message
       stack%message_list => localstack%other_messages
       deallocate(localstack)
    else
       message = 'Trying to get message from empty stack.'
    endif
end subroutine getmessage
!
!
!
!==============================================================================
subroutine writemessages(stack,unit)
    !    Function: - Write all messages to file
    !
    implicit none
    !
    ! Call variables
    !
    type(message_stack)          :: stack
    integer                      :: unit
    !
    ! Local variables
    !
    character(message_len)       :: message
    !
    ! body
    !
    do while (.not.isempty(stack))
       call getmessage(stack,message)
       write(unit,'(A)') trim(message)
    end do
end subroutine writemessages
!
!
!
!==============================================================================
subroutine clearstack(stack)
    !    Function: - Get a message from the top of the message stack
    !
    implicit none
    !
    ! Call variables
    !
    type(message_stack)    :: stack
    character(message_len) :: message
    !
    ! Local variables
    !
    type(message_type), pointer :: localstack
    !
    ! body
    !
    do while (associated(stack%message_list))
       localstack => stack%message_list
       message = localstack%message
       stack%message_list => localstack%other_messages
       deallocate(localstack)
    end do
end subroutine clearstack

end module message_module
