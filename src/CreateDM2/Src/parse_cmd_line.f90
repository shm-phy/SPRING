
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

module parse_cmd_line

    use kinds,              only : dp
    implicit none

    private

    public  ::  get_command_line

contains

subroutine get_command_line(infile, T, mu, alpha)

    implicit none

    character(len=512), parameter       :: expected_cmd_line="./a.out -inp [INPUT FILE] -T [TEMPERATURE] &
                                         & -mu [INTEGER] -alpha [INTEGER]"
    integer, parameter                  :: MAX_CMD_ARG=8

    character(len=*), intent(out)       :: infile
    real(dp), intent(out)               :: T
    integer, intent(out)                :: mu, alpha
    
    character(len=64)                   :: Temp_char, mu_char, &
                                         & alpha_char
    logical, dimension(MAX_CMD_ARG/2)   :: dflt_val

    character(len=1024)                 :: cmd_line
    integer                             :: cmd_line_len
    integer                             :: Narg
    character(len=64)                   :: cmd_name, cmd_arg
    integer                             :: cmd_arg_len
    integer                             :: nn, i, err


    call get_command(cmd_line, cmd_line_len, err)

    err_chk0: if ( err /= 0 ) then

        write(*, '(A)') "Command line arguments are incorrectly given"

        large_len: if ( err == -1 ) then
                write(*, '(A)') "Command line is to large to fit in 1024 len character"
        end if large_len

        write(*, 100) expected_cmd_line
        write(*, 200) cmd_line

        100 FORMAT('Expected command line: ', A)
        200 FORMAT('Given command line: ', A)
        
    end if err_chk0

    Narg = command_argument_count()
    dflt_val(:) = .true.

    chk_no_cmd: if ( Narg > 0 .and. Narg <= MAX_CMD_ARG ) then
        
        nn = 1
        loop_args: do 
                
            call get_command_argument(nn, cmd_name, cmd_arg_len, err)
            nn = nn + 1
            
            err_chk: if ( err == 0) then
            
                chk_cmd_name: if ( cmd_name == "-inp" ) then
            
                    call get_command_argument(nn, cmd_arg, cmd_arg_len, err)

                    if ( err == 0 ) then
                        infile = cmd_arg
                        dflt_val(1) = .false.

                        img1_chk1: if ( this_image() == 1) then 
                            write(*, 500) infile
                            500 FORMAT('Input file to be read: ', A64)
                        end if img1_chk1
                    else
                        call print_err_msg(err, nn, cmd_line, expected_cmd_line)

                    end if
                    
                    nn = nn + 1

                else if ( cmd_name == "-T" ) then chk_cmd_name

                    call get_command_argument(nn, cmd_arg, cmd_arg_len, err)
                    if ( err == 0 ) then
                        Temp_char = cmd_arg
                        dflt_val(2) = .false.
                        read(Temp_char, *) T

                        img1_chk2: if ( this_image() == 1) then 
                            write(*, 600) T
                            600 FORMAT('Temperature: ', F7.2, 'K')
                        end if img1_chk2

                    else
                        call print_err_msg(err, nn, cmd_line, expected_cmd_line)

                    end if

                    nn = nn + 1

                else if ( cmd_name == "-mu" ) then chk_cmd_name

                    call get_command_argument(nn, cmd_arg, cmd_arg_len, err)
                    if ( err == 0 ) then
                        mu_char = cmd_arg
                        dflt_val(3) = .false.
                        read(mu_char, *) mu

                    else
                        call print_err_msg(err, nn, cmd_line, expected_cmd_line)

                    end if

                    nn = nn + 1

                else if ( cmd_name == "-alpha" ) then chk_cmd_name

                    call get_command_argument(nn, cmd_arg, cmd_arg_len, err)
                    if ( err == 0 ) then
                        alpha_char = cmd_arg
                        dflt_val(4) = .false.
                        read(alpha_char, *) alpha

                    else
                        call print_err_msg(err, nn, cmd_line, expected_cmd_line)

                    end if

                    nn = nn + 1

                else chk_cmd_name
                    
                    if ( this_image() == 1 ) write(*, 320) cmd_name 
                    320 FORMAT('WARNING: Unrecognised command name -> ', A64)

                    nn = nn + 1
            
                end if chk_cmd_name

            else err_chk
                call print_err_msg(err, nn, cmd_line, expected_cmd_line)
            
            end if err_chk
            
            if ( nn > Narg ) exit loop_args
        
        end do loop_args

    else if ( Narg > MAX_CMD_ARG ) then chk_no_cmd
        call print_err_msg(-5, 0, cmd_line, expected_cmd_line)

    end if chk_no_cmd

    dflt_chk_loop: do i = 1, size(dflt_val)

        outer_if: if ( dflt_val(i) ) then

            inner_if: if ( i == 1 ) then
                infile = 'input.nml'
                img1_chk3: if ( this_image() == 1) then 
                    write(*, 550) infile
                    550 FORMAT('Input file is set to default: ', A64)
                end if img1_chk3

            else if ( i == 2 ) then inner_if
                T = 300.0_dp
                img1_chk4: if ( this_image() == 1) then 
                    write(*, 650) T
                    650 FORMAT('Temperature is set to default: ', F7.2, 'K')
                end if img1_chk4

            else if ( i == 3 ) then inner_if
                mu = 1
                img1_chk5: if ( this_image() == 1) then 
                    write(*, 350) mu
                    350 FORMAT('mu set to default: ', I3)
                end if img1_chk5

            else if ( i == 4 ) then inner_if
                alpha = 1
                img1_chk6: if ( this_image() == 1) then 
                    write(*, 950) alpha
                    950 FORMAT('alpha is set to default: ', I3)
                end if img1_chk6

            end if inner_if

        end if outer_if

    end do dflt_chk_loop

end subroutine get_command_line


subroutine print_err_msg(err, nn, cmd_line, expected_cmd_line)

    implicit none
    integer, intent(in)                 :: err, nn
    character(len=*), intent(in)        :: expected_cmd_line
    character(len=*), intent(in)        :: cmd_line

    err_chk: if ( err /= 0 ) then

            write(*, 300) nn
            write(*, 100) expected_cmd_line
            write(*, 200) cmd_line

            300 FORMAT('Error in reading argument number ', I3)
            100 FORMAT('Expected command line: ', A64)
            200 FORMAT('Given command line: ', A64)

    end if err_chk

end subroutine print_err_msg

end module parse_cmd_line

