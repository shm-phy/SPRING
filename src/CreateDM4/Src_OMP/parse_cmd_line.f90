
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

    public  ::  get_command_line, ShowWelcomeBanner

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

                         write(*, 500) infile
                         500 FORMAT('Input file to be read: ', A64)

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

                        write(*, 600) T
                        600 FORMAT('Temperature: ', F7.2, 'K')
                        

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

                    write(*, 320) cmd_name
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
                
                write(*, 550) infile
                550 FORMAT('Input file is set to default: ', A64)

            else if ( i == 2 ) then inner_if
                T = 300.0_dp
                
                write(*, 650) T
                650 FORMAT('Temperature is set to default: ', F7.2, 'K')

            else if ( i == 3 ) then inner_if
                mu = 1
                
                write(*, 350) mu
                350 FORMAT('mu set to default: ', I3)

            else if ( i == 4 ) then inner_if
                alpha = 1
                
                write(*, 950) alpha
                950 FORMAT('alpha is set to default: ', I3)

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


    subroutine ShowWelcomeBanner()

        implicit none

        write(*, 10)
        write(*, 20)
        write(*, 30)
#ifdef __SPR_VERSION__
        write(*, 22) __SPR_VERSION__
        22 FORMAT(10X, "Version:", A11)
#endif
        write(*, 40)
#ifdef __COMPILE_TIME__
        write(*, 12) __COMPILE_TIME__
        12 FORMAT( 10X, "Compiled on: ", A20)
#endif
        write(*, 50)
        write(*, 10)

        10 FORMAT(6X, 5('='), 70('-'), 5('='))
        20 FORMAT( &
                   10X, "     _______..______   .______       __  .__   __.   _______ ", /, &
                   10X, "    /       ||   _  \  |   _  \     |  | |  \ |  |  /  _____|", /, &
                   10X, "   |   (----`|  |_)  | |  |_)  |    |  | |   \|  | |  |  __  ", /, &
                   10X, "    \   \    |   ___/  |      /     |  | |  . `  | |  | |_ | ", /, &
                   10X, ".----)   |   |  |      |  |\  \----.|  | |  |\   | |  |__| | ", /, &
                   10X, "|_______/    | _|      | _| `._____||__| |__| \__|  \______| " &
                 )
        30 FORMAT(10X, "Self-consistent Phonon RenormalizING-package")
        40 FORMAT(10X, "Executable: create_disp_mat4_omp.x (Creation of 4th-order displacement matrix)")
        50 FORMAT(10X, "If you use this program for any publication, please cite: ", /, &
                  10X, "https://arxiv.org/abs/2402.02787", / )

    end subroutine ShowWelcomeBanner


end module parse_cmd_line

