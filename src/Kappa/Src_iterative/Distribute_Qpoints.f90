
! ============================================================================= !
!    Copyright (C) 2022  Soham Mandal                                           !
!                                                                               !
!    This program is free software: you can redistribute it and/or modify       !
!    it under the terms of the GNU General Public License as published by       !
!    the Free Software Foundation, either version 3 of the License, or          !
!    (at your option) any later version.                                        !
!                                                                               !
!    This program is distributed in the hope that it will be useful,            !
!    but WITHOUT ANY WARRANTY; without even the implied warranty of             !
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              !
!    GNU General Public License for more details.                               !
!                                                                               !
!    You should have received a copy of the GNU General Public License          !
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.     !
!                                                                               !
!    e-mail: phy.soham@gmail.com                                                !
! ============================================================================= !

module DistributeQs

    implicit none
    private

    public          :: DistributeQpoints, ShowProcessorDistribution

contains

    subroutine DistributeQpoints(Nqs, my_Qsize, my_offset, my_edge)

        implicit none

        integer, intent(in)                                 :: Nqs
        integer, intent(out)                                :: my_Qsize, my_offset
        integer, dimension(2), intent(out)                  :: my_edge

        ! ================================= Local variable =================================== !

        character(len=128)                                  :: msg
        integer, dimension(:), allocatable                  :: QinEachProcs
        integer                                             :: Nprocs, div_int, rem, istat

        ! ================================= Local variable =================================== !

        Nprocs = num_images()

        allocate( QinEachProcs( Nprocs ) )
        QinEachProcs = 0

        OnlyImage1: if ( this_image() == 1 ) then

            div_int = Nqs / Nprocs
            rem = mod( Nqs, Nprocs )

            Warn: if ( div_int == 0 ) then

                write(*, 25) Nqs, Nprocs
                25 FORMAT('Number of q-points ( ', I5, ' ) is smaller than number of processors ( ', I5, ' ) ')

            end if Warn

            QinEachProcs( : ) = div_int
            QinEachProcs( 1:rem ) = QinEachProcs( 1:rem ) + 1

            Check: if ( sum(QinEachProcs) /= Nqs ) then

                write(*, 35)
                35 FORMAT( 'Something wrong in DistributeQpoints ')
                ERROR STOP

            end if Check

        end if OnlyImage1

        SYNC ALL
        call co_broadcast( QinEachProcs, source_image=1, stat=istat, ERRMSG=msg )
        if ( istat /= 0 ) write(*, "( 'ERROR in co_broadcast : ', A128 )") msg
        SYNC ALL

        my_Qsize = QinEachProcs(this_image())

        my_edge(1) = sum( QinEachProcs( 1 : (this_image() - 1) ) ) + 1
        my_edge(2) = my_edge(1) + my_Qsize - 1

        my_offset = my_edge(1) - 1

        deallocate( QinEachProcs )

    end subroutine DistributeQpoints


    subroutine ShowProcessorDistribution( character_to_print, &
                                        & my_edge, my_Qsize, my_offset )

        implicit none

        character(len=*)                                :: character_to_print
        integer, dimension(2), intent(in)               :: my_edge
        integer, intent(in)                             :: my_Qsize, my_offset

        ! =============================== Local variables =============================== !

        integer                             :: char_len, left, right

        character(len=8)                    :: frmt
        character(len=8)                    :: int1_to_char, int2_to_char
        character(len=256)                  :: write_frmt

        ! =============================== Local variables =============================== !

        char_len = LEN_TRIM( character_to_print )

        left = (111 - char_len - 2) / 2
        right = 111 - (left + char_len + 2)

        frmt = '(I5)'

        write(int1_to_char, frmt) left
        write(int2_to_char, frmt) right

        write_frmt = "(10X, '|', "//trim(adjustl(adjustr(int1_to_char)))//"('='), ' ', '" &
                   & //trim(adjustl(adjustr(character_to_print)))//"', ' ', "&
                   & //trim(adjustl(adjustr(int2_to_char)))//"('='), '|')"

        ! ////////////////////  Show the qpoints distributions among processors \\\\\\\\\\\\\\\\\\\\\ !
        img1_chk1: if ( this_image() == 1 ) then
            write(*, *)
            call execute_command_line(' ')

            write(*, write_frmt)
            call execute_command_line(' ')

            write(*, 25)
            call execute_command_line(' ')

        end if img1_chk1

        SYNC ALL

        write(*, 30) this_image(), my_edge, my_Qsize, my_offset
        call execute_command_line(' ')

        SYNC ALL

        img1_chk2: if ( this_image() == 1 ) then

            write(*, 25)
            call execute_command_line(' ')

            write(*, write_frmt)
            call execute_command_line(' ')
            write(*, *)
            call execute_command_line(' ')

        end if img1_chk2

        25 FORMAT(10X,'|',111X,'|')
        30 FORMAT(10X,'| **', 1X,' Image No. = ', I4, ', Q-points covered: (', I9, '=> ', I9, &
                & '). my_Qsize = ', I9, ', my_offset = ', I7, 1X,'** |')
        ! \\\\\\\\\\\\\\\\\\\\  Show the qpoints distributions among processors ///////////////////// !

    end subroutine ShowProcessorDistribution

end module DistributeQs

