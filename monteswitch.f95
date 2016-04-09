!
! By extracting lines beginning with the regular expression '\!\! ?'
! (ignoring leading whitespace), and then removing matches to the
! regular expression, html documentation corresponding to this 
! source code will be created.
!
!! <html>
!! <head>
!!  <title> monteswitch Documentation </title>
!! </head>
!! <body>
!! 
!! <h1> <code>monteswitch</code> Documentation </h1>
!!
!! <h2> Author </h2>
!! <p> Tom Underwood </p>
!!
!! <h2> General Notes </h2>
!! <p>
!! The program performs a simulation by calling the <code>run</code> subroutine in the <code>monteswitch_mod</code>
!! module. The program always creates two output files: 'data' and 'state', which correspond to the <code>datafile</code>
!! and <code>statefile</code> arguments of the <code>run</code> subroutine. 
!! </p>
!! <p>
!! This program terminates with a non-zero exit status of 1 for most errors. If the system has 
!! melted then the exit status is 2.  If the energy has diverged due to precision-related
!! errors then the exit status is 3. Note that exit statuses are not part of the Fortran standard, and may not 
!! work for all operating systems or compilers.
!! </p>
!! <p>
!! Usage of this program is: <code>monteswitch [-seed <i>seed</i>] -new [-wf]</code> or
!! <code>monteswitch [-seed <i>seed</i>] (-resume|-reset)</code>, where the command-line argument options are as follows.
!! </p>
!! <table border="1">
!!  <tr>
!!   <td> <b> Argument </b> </td>
!!   <td> <b> Description </b> </td>
!!  </tr>
!!  <tr>
!!   <td> -seed <i>seed</i> </td>
!!   <td>
!!   Specify the integer seed for the random number generator to be <i>seed</i>. <i>This command line option, if present,
!!   should precede the other possible arguments below</i>. Note that the seed must be non-zero.
!!   </td>
!!  </tr>
!!  <tr>
!!   <td> -new </td>
!!   <td>
!!   Start a new simulation from scratch; the simulation is initialised with the <code>initialise_from_files</code>
!!   subroutine in the <code>monteswitch_mod</code> module. In this case there are three input files: 'params_in', 
!!   'lattices_in' and 'interactions_in', which correspond to the arguments to the <code>initialise_from_files</code>
!!   in an obvious way. For this option the file 'data' if it exists is overwritten.
!!   </td>
!!  </tr>
!!  <tr>
!!   <td> -wf </td>
!!   <td>
!!   When used in conjunction with <code>-new</code>, the initial weight function to be used by the simulation is read 
!!   from a file 'wf_in'. The file must contain <code>M_grid_size</code> (specified in 'params_in') lines, each 
!!   containing two tokens (extra lines and tokens are ignored). Both tokens should be of type <code>real</code>. 
!!   The first token on line <i>i</i> is ignored, and the second is the value of the weight function for macrostate <i>i</i>.
!!   </td>
!!  </tr>
!!  <tr>
!!   <td> -resume </td>
!!   <td>
!!   Resume an old simulation from the file 'state;. In this case the <code>import</code> subroutine in the 
!!   <code>monteswitch_mod</code> module is used to initialise the simulation from the file 'state'. The simulation 
!!   is then 'resumed', and the file 'data' if it exists is ammended.
!!   </td>
!!  </tr>
!!  <tr>
!!   <td> -reset </td>
!!   <td>
!!   Start a new simulation from the file 'state': the program imports the state from the file 'state' using the
!!   <code>import</code> subroutine in the <code>monteswitch_mod</code> module, but resets all counters before
!!   'resuming' the simulation. In this case the file 'data' if it exists is overwritten. Note that <code>trans</code> 
!!   is not regarded as a counter, and hence will not be reset. Similarly for the weight function.
!!   </td>
!!  </tr>
!! </table>
!! </p>
!!
program monteswitch

    !! <h2> Dependencies </h2>
    !! <p> 
    !! <ul>
    !!  <li> <code> kinds_mod </code> </li>
    !!  <li> <code> monteswitch_mod </code> </li>
    !! </ul>
    !! </p>
    use kinds_mod
    use monteswitch_mod

    implicit none
    
    ! This is .true. if the argument "-seed" is present
    logical :: seed_specified
    ! This is .true. if the argument "-wf" is present
    logical :: wf_specified
    ! The RNG seed if it is specified
    integer :: seed
    ! The command line argument we are on for reading: we have read n-1 command line arguments so far
    integer :: n

    ! Unimportant variables
    character(len=20) :: char
    integer :: error, i
    real(rk) :: x


    seed_specified = .false.
    wf_specified = .false.
    n=1

    ! Get the 1st command line argument
    call getarg(n,char)
    n = n+1

    ! Check that there actually is an argument
    if(trim(char)=="") then

        write(0,*) "monteswitch: Error. No command line argument detected."
        stop 1

    end if

    ! If the 1st argument is "-seed" then read the next argument and use it as the RNG seed
    if(trim(char)=="-seed") then

        seed_specified = .true.

        call getarg(n,char)
        n = n+1

        read(char,*,iostat=error) seed

        ! Check that the seed is sensible
        if(error/=0) then
            write(0,*) "monteswitch: Error. Problem reading integer after the command line argument '-seed'"
            stop 1
        end if
        if(seed==0) then
            write(0,*) "monteswitch: Error. The RNG seed must be non-zero."
            stop 1
        end if

        ! Read the next command line argument: we are still expecting one of "-new", "-resume" or "-reset"
        call getarg(n,char)
        n = n+1

    end if


    ! If the 1st argument (or 3rd argument if the 1st argument is "-seed") is "-new", "-resume", 
    ! "-reset" or something else
    select case(trim(char))

    case("-new")


        ! CODE FOR A NEW SIMULATION

        ! Check for the presence of "-wf"
        call getarg(n,char)
        n = n+1
        select case(trim(char))
        case("-wf")
            wf_specified = .true.
            ! Check that there aren't any more command line arguments
            call getarg(n,char)
            if(trim(char)/="") then
                write(0,*) "monteswitch: Error. Command line argument detected after '-wf' where there should be none"
                stop 1
            end if
        case("")            
            ! Do nothing
        case default 
            write(0,*) "monteswitch: Error. Only argument '-wf' is allowed after '-new'"
            stop 1
        end select

        call initialise_from_files("params_in","lattices_in")

        ! Set the initial weight function from the file "wf_in" if "-wf" is present
        if(wf_specified) then

            ! Open the file
            open(unit=10,file="wf_in",iostat=error,status="old")
            if(error/=0) then
                write(0,*) "monteswitch: Error. Problem opening file 'wf_in'"
                stop 1
            end if
            
            ! Read the weight function from the file and set eta_grid accordingly        
            do i=1,M_grid_size
                read(10,*,iostat=error) x, eta_grid(i)
                if(error/=0) then
                    write(0,*) "monteswitch: Error. Problem reading weight function from line ",i," in file 'wf_in'"
                    stop 1
                end if
            end do
            
            close(unit=10)

        end if

        ! Run the simulation
        if(seed_specified) then
            call run("data","state",.false.,seed)
        else
            call run("data","state",.false.)
        end if


    case("-resume")


        ! CODE FOR A RESUMED SIMULATION

        ! Check that there aren't any more command line arguments
        call getarg(n,char)
        if(trim(char)/="") then
            write(0,*) "monteswitch: Error. Command line argument detected after '-resume' where there should be none"
            stop 1
        end if

        call import("state")

        if(seed_specified) then
            call run("data","state",.true.,seed)
        else
            call run("data","state",.true.)
        end if


    case("-reset")


        ! CODE FOR SIMULATION INITIALISED FROM THE STATE FILE BUT WITH COUNTERS RESET

        ! Check that there aren't any more command line arguments
        call getarg(n,char)
        if(trim(char)/="") then
            write(0,*) "monteswitch: Error. Command line argument detected after '-reset' where there should be none"
            stop 1
        end if

        call import("state")
        call initialise_counters()

        if(seed_specified) then
            call run("data","state",.false.,seed)
        else
            call run("data","state",.false.)
        end if

    case("")

        ! Code for missing 3rd argument if "-seed" is present

        write(0,*) "monteswitch: Error. Failed to detect argument '-new', '-resume' or '-reset'"
        stop 1

    case("-wf")

        ! Code for out of place "-wf" argument

        write(0,*) "monteswitch: Error. Incorrect usage: argument '-wf' must follow '-new'"
        stop 1

    case default

        ! Code for unknown 1st (or 3rd if "-seed" is present) argument

        write(0,*) "monteswitch: Error. Unrecognised command line argument '",trim(char),"'"
        stop 1

    end select


end program
!!
!! </body>
!! </html>
