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
!! The command line options are as follows:
!! </p>
!! <table border="1">
!!  <tr>
!!   <td> <b> Argument </b> </td>
!!   <td> <b> Description </b> </td>
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
  
  ! Unimportant variables
  character(len=20) :: char

  ! THE PROGRAM
  

  ! Read the first command line argument 
  call getarg(1,char)
  if(char=="") then
     write(0,*) "monteswitch: Error. No command line argument detected."
     stop 1
  else if(trim(char)=="-new") then

     ! CODE FOR A NEW SIMULATION

     call initialise_from_files("params_in","lattices_in")
     call run("data","state",.false.)

  else if(trim(char)=="-resume") then

     ! CODE FOR A RESUMED SIMULATION
     
     call import("state")
     call run("data","state",.true.)

  else if(trim(char)=="-reset") then
     
     ! CODE FOR SIMULATION INITIALISED FROM THE STATE FILE BUT WITH COUNTERS RESET
     
     call import("state")
     call initialise_counters()
     call run("data","state",.false.)

  else

     ! CODE FOR ANY OTHER COMMAND LINE ARGUMENTS

     write(0,*) "monteswitch: Error. Unrecognised first command line argument."
     stop 1

  end if


end program monteswitch
!!
!! </body>
!! </html>
