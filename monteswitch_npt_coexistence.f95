!
! By extracting lines beginning with the regular expression '\!\! ?'
! (ignoring leading whitespace), and then removing matches to the
! regular expression, html documentation corresponding to this 
! source code will be created.
!
!! <html>
!! <head>
!!  <title> monteswitch_npt_coexistence Documentation </title>
!! </head>
!! <body>
!! 
!! <h1> <code>monteswitch_npt_coexistence</code> Documentation </h1>
!!
!! <h2> Author </h2>
!! <p> Tom Underwood </p>
!!
!! <h2> General Notes </h2>
!! <p>
!! This program determines the coexistence curve for an NPT ensemble for a specified range of pressures.
!! </p>
!! <p>
!! The variables of the simulation are specified in a file 'init_coexistence', and are listed below.
!! </p>
!! <p>
!! This program creates the following output files:
!! 'coexistence_points': contains the coordinates (P,beta) of coexistence points found by the program, with the third
!! token in each line being the uncertainty in beta (specifically, the upper bound on the uncertainty in the Newton-Rhapson
!! change in beta).
!! 'examined_points': contains the coordinates (P,beta) of all points in phase space explored by the program
!! 'wf_P_BETA': contains the weight function vs. order parameter for the P and beta
!! 'vs_P_BETA': contains the number of times each order parameter was visited during a batch simulation. The first
!! column pertains to lattice 1, while the second pertains to lattice 2.
!! 'state_P_BETA': the 'state' file after the batch simulation for the P and beta.
!! </p>
!!
program monteswitch_npt_coexistence

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

  !! <h2> Variables for new simulations which are specified in the file 'init_coexistence'</h2>
  !! <p>
  !! The following variables <i>must</i> be contained within the input file 'init_coexistence', in the order specified below, 
  !! in the following format: Each line corresponds to a variable. Each line contains two tokens, with whitespace
  !! acting as delimitters. The first token is a descriptive string, e.g., <code>rho=</code> (note the lack of space
  !! before the equals sign) for the line pertaining to the variable <code>rho</code>. This string's sole purpose is
  !! to allow the 'init' file to be more easily readable; it is not used by the program. The second token is the
  !! value of the variable corresponding to the line. For example, the first few lines of 'init' could be
  !! <code>
  !!  rho= 1.0901668831
  !!  nx= 2
  !!  ny= 2
  !! </code> 
  !! </p>
  !! <p>
  !! Most of the variables are inheritted from the module <code>monteswitch_mod</code>, but, in addition to playing
  !! the role specified in that module, some will also determine how a simulation will be initialised by this program.
  !! Descriptions are given below only for those variables not inherited from <code>monteswitch_mod</code>, or if the
  !! value of the variable determines how the simulation will be initialised. Note that even variables which are not 
  !! relevant to the particular simulation you wish to run must be inlcuded in the 'init_coexistence' file, and given a value.
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Variable </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> rho_init </code> </td>
  !!   <td> <code> real(rk) </code> </td>
  !!   <td>  The starting density for the lattices. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> nx </code> </td>
  !!   <td> <code> integer(ik) </code> </td>
  !!   <td> 
  !!   The number of unit cells to create in the x direction to form the supercell.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> ny </code> </td>
  !!   <td> <code> integer(ik) </code> </td>
  !!   <td> 
  !!   The number of unit cells to create in the y direction to form the supercell.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> nz </code> </td>
  !!   <td> <code> integer(ik) </code> </td>
  !!   <td> 
  !!   The number of unit cells to create in the z direction to form the supercell.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> init_lattice </code> </td>
  !!   <td> <code> integer(ik) </code> </td>
  !!   <td> 
  !!   The lattice number (1 or 2) which the simulation will be initialised in.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> lj_epsilon </code> </td>
  !!   <td> <code> real(rk) </code> </td>
  !!   <td> 
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> lj_sigma </code> </td>
  !!   <td> <code> real(rk) </code> </td>
  !!   <td> 
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> lj_cutoff </code> </td>
  !!   <td> <code> real(rk) </code> </td>
  !!   <td> 
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> list_cutoff </code> </td>
  !!   <td> <code> real(rk) </code> </td>
  !!   <td> 
  !!   The radius which defines which sites interact: sites within 
  !!   <code>list_cutoff</code> of each other when in the initial perfect bcc and fcc
  !!   positions will interact throughout the simulation.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> M_grid_size </code> </td>
  !!   <td> <code> real(rk) </code> </td>
  !!   <td> 
  !!   The number of 'bins' which the order parameter grid will have.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> M_grid_min_kT </code> </td>
  !!   <td> <code> real(rk) </code> </td>
  !!   <td> 
  !!   The minimum value which the order parameter grid will support, in units of <i>kT</i>, where <i>k</i> denotes the
  !!   Boltzmann constant. In other words, the lower limit of order parameter which the simulation will be
  !!   able to cope with.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> M_grid_max_kT </code> </td>
  !!   <td> <code> real(rk) </code> </td>
  !!   <td> 
  !!   The maximum value which the order parameter grid will support, in units of <i>kT</i>, where <i>k</i> denotes the
  !!   Boltzmann constant. In other words, the upper limit of order parameter which the simulation will be
  !!   able to cope with.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!    <td><code> enable_multicanonical_generation </code></td> 
  !!    <td><code>logical</code></td>
  !!    <td>The value of <code>enable_multicanonical</code> for weight function generation simulations. Its value is
  !!     <code>.true.</code> for batch simulations</td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> P_min </code></td> 
  !!    <td><code>real(rk)</code></td>
  !!    <td>The minimum value of <code>P</code> to consider. </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> P_max </code></td> 
  !!    <td><code>real(rk)</code></td>
  !!    <td>The maximum value of <code>P</code> to consider. </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> Delta_P </code></td> 
  !!    <td><code>real(rk)</code></td>
  !!    <td>The increment between subsequent values of <code>P</code>, i.e. the resolution of the coexistence curve 
  !!    in <code>P</code>-space.</td> 
  !!  </tr>
  !!  <tr>
  !!   <td><code> beta_init </code></td> 
  !!   <td><code>real(rk)</code></td> 
  !!   <td>The initial value of <code>beta</code> to consider. </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> part_select </code></td> 
  !!     <td><code>character(len=30)</code></td> 
  !!     <td> </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> part_step </code></td> 
  !!     <td><code>real(rk)</code></td> 
  !!     <td> </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> enable_COM_frame </code></td> 
  !!     <td><code>logical</code></td> 
  !!     <td> </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> vol_dynamics </code></td> 
  !!     <td><code>character(len=30)</code></td> 
  !!     <td> </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> vol_freq </code></td> 
  !!     <td><code>integer(ik)</code></td> 
  !!     <td> </td> 
  !!  </tr>

  !!  <tr>
  !!    <td><code> vol_step </code></td> 
  !!     <td><code>real(rk)</code></td> 
  !!     <td> </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> stop_sweeps_generation </code></td> 
  !!     <td><code>integer(ik)</code></td> 
  !!     <td>The number of sweeps in a weight function generation simulation. </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> stop_sweeps_batch </code></td> 
  !!     <td><code>integer(ik)</code></td> 
  !!     <td>The number of sweeps in a batch simulation. </td> 
  !!  </tr>
  !! <tr>
  !!  <td> <code>enable_melt_checks</code> </td>
  !!  <td> <code>logical</code> </td>
  !!  <td> </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>melt_sweeps</code> </td>
  !!  <td> <code>integer(ik)</code> </td>
  !!  <td> </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>melt_threshold</code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>melt_option</code> </td>
  !!  <td> <code>character(len=30)</code> </td>
  !!  <td> </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>enable_divergence_checks</code> </td>
  !!  <td> <code>logical</code> </td>
  !!  <td> 
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>divergence_sweeps</code> </td>
  !!  <td> <code>integer(ik)</code> </td>
  !!  <td> 
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>divergence_tol</code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> 
  !!  </td>
  !! </tr>
  !!  <tr>
  !!    <td><code> update_eta_sweeps </code></td> 
  !!     <td><code>integer(ik)</code></td> 
  !!    <td> </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> update_trans_generation </code></td> 
  !!     <td><code>logical</code></td> 
  !!     <td>
  !!     The value of <code>update_trans</code> in weight function generation simulations. 
  !!     If set to <code>.true.</code> then <code>trans</code> is initialised; if not then <code>trans</code>
  !!     is not initialised.
  !!     </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> update_eta_method </code></td> 
  !!     <td><code>logical</code></td> 
  !!     <td></td>
  !!  </tr>
  !!  <tr>
  !!    <td><code> enable_barriers_generation </code></td> 
  !!    <td><code>logical</code></td> 
  !!    <td>
  !!     The value of <code>enable_barriers</code> in weight function generation simulations. 
  !!    If set to <code>.true.</code> then the variables pertaining to barriers are initialised.
  !!    </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> barrier_dynamics </code></td> 
  !!    <td><code>character(len=30)</code></td> 
  !!    <td>
  !!    </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> lock_moves </code></td> 
  !!     <td><code>integer(ik)</code></td> 
  !!     <td> </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> block_sweeps </code></td> 
  !!     <td><code>integer(ik)</code></td> 
  !!     <td> </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> equil_sweeps </code></td> 
  !!     <td><code>integer(ik)</code></td> 
  !!     <td> </td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> nr_max </code></td> 
  !!     <td><code>integer(ik)</code></td> 
  !!     <td>The maximum number of iterations which the Newton-Rhapson loop will perform in an attempt to find a
  !!       coexistence point.</td> 
  !!  </tr>
  !!  <tr>
  !!    <td><code> shell_command </code></td> 
  !!     <td><code>character(len=100)</code></td> 
  !!     <td> 
  !!     This is the command passed to the shell to execute the <code>monteswitch</code> or <code>monteswitch_mpi</code> binary
  !!     with the argument <code>-reset</code.
  !!     For example, it could be <code>"../monteswitch -reset"</code> if the binary is in the directory above the current one; or
  !!     to run the OpenMPI version using 4 processes it could be <code>"mpiexec -n 4 ../monteswitch_mpi -reset"</code>. Note that
  !!     the length of <code>shell_command</code> is limited to 500 characters. If one needs more than 500 characters,
  !!     then defining shell variables and using these in <code>shell_command</code> is one way to get around the limit.
  !!     </td> 
  !!  </tr>
  !! </table>

  ! Variables not inherited from monteswitch_mod
  real(rk) :: rho_init
  integer(ik) :: nx
  integer(ik) :: ny
  integer(ik) :: nz
  integer(ik) :: init_lattice
  real(rk) :: list_cutoff
  real(rk) :: M_grid_min_kT
  real(rk) :: M_grid_max_kT
  logical :: enable_multicanonical_generation
  real(rk) :: P_min
  real(rk) :: P_max
  real(rk) :: Delta_P
  real(rk) :: beta_init
  integer(ik) :: stop_sweeps_generation
  integer(ik) :: stop_sweeps_batch
  logical :: update_trans_generation
  logical :: enable_barriers_generation
  integer(ik) :: nr_max
  character(len=500) :: shell_command
 
  ! The change in beta which gives the coexistence point, as predicted by Newton-Rhapson
  real(rk) :: Delta_beta
  ! Quantities related to the uncertainty in the estimated coexistence beta using Newton-Rhapson
  real(rk) :: sigma_H1H2_upper
  real(rk) :: sigma_H1H2_lower
  real(rk) :: sigma_Delta_beta_upper
  real(rk) :: sigma_Delta_beta_lower

  ! Unimportant variables
  integer(ik) :: error, i, j
  character(len=100) :: char_P, char_beta




  ! INITIALISATION


  ! Import the 'init_coexistence' variables from that file
  call import_init_coexistence()

  ! Set some variables which will be used for all simulations
  enable_lattice_moves= .true.
  enable_part_moves= .true.
  enable_vol_moves= .true.
  output_file_period=-1
  output_file_V= .false.
  output_file_R_1= .false.
  output_file_R_2= .false.
  output_file_u= .false.
  output_file_lattice= .false.
  output_file_E= .false.
  output_file_M= .false.
  output_file_moves_lattice= .false.
  output_file_accepted_moves_lattice= .false.
  output_file_moves_part= .false.
  output_file_accepted_moves_part= .false.
  output_file_moves_vol= .false.
  output_file_accepted_moves_vol= .false.
  output_file_rejected_moves_M_OOB= .false.
  output_file_M_OOB_low= .false.
  output_file_M_OOB_high= .false.
  output_file_barrier_macro_low= .false.
  output_file_barrier_macro_high= .false.
  output_file_rejected_moves_M_barrier= .false.
  output_file_moves_since_lock= .false.
  output_file_melts= .false.
  output_file_equil_DeltaF= .false.
  output_file_sigma_equil_DeltaF= .false.
  output_file_equil_H_1= .false.
  output_file_sigma_equil_H_1= .false.
  output_file_equil_H_2= .false.
  output_file_sigma_equil_H_2= .false.
  output_file_equil_V_1= .false.
  output_file_sigma_equil_V_1= .false.
  output_file_equil_V_2= .false.
  output_file_sigma_equil_V_2= .false.
  output_stdout_period= -1
  output_stdout_V= .false.
  output_stdout_R_1= .false.
  output_stdout_R_2= .false.
  output_stdout_u= .false.
  output_stdout_lattice= .false.
  output_stdout_E= .false.
  output_stdout_M= .false.
  output_stdout_moves_lattice= .false.
  output_stdout_accepted_moves_lattice= .false.
  output_stdout_moves_part= .false.
  output_stdout_accepted_moves_part= .false.
  output_stdout_moves_vol= .false.
  output_stdout_accepted_moves_vol= .false.
  output_stdout_rejected_moves_M_OOB= .false.
  output_stdout_M_OOB_low= .false.
  output_stdout_M_OOB_high= .false.
  output_stdout_barrier_macro_low= .false.
  output_stdout_barrier_macro_high= .false.
  output_stdout_rejected_moves_M_barrier= .false.
  output_stdout_moves_since_lock= .false.
  output_stdout_melts= .false.
  output_stdout_equil_DeltaF= .false.
  output_stdout_sigma_equil_DeltaF= .false.
  output_stdout_equil_H_1= .false.
  output_stdout_sigma_equil_H_1= .false.
  output_stdout_equil_H_2= .false.
  output_stdout_sigma_equil_H_2= .false.
  output_stdout_equil_V_1= .false.
  output_stdout_sigma_equil_V_1= .false.
  output_stdout_equil_V_2= .false.
  output_stdout_sigma_equil_V_2= .false.
  ! The below value of 'checkpoint_period' gives no checkpointing during any simulation, but checkpointing at the
  ! end of every simulation, i.e., we get a 'state' file only at the end of every simulation.  
  checkpoint_period=stop_sweeps_generation+stop_sweeps_batch+1 



  ! DETERMINE THE COEXISTENCE CURVE


  
  ! Wipe the 'coexistence_points' and 'examined_points' files beforehand
  open(unit=11,file="coexistence_points",status="replace")
  close(unit=11)
  open(unit=11,file="examined_points",status="replace")
  close(unit=11)


  ! Set initial beta and P
  beta=beta_init
  P= P_min

  do
     if(P>P_max) then
        write(*,*)
        write(*,*) "Completed all simulations."
        write(*,*)
        stop 0
     end if

     write(*,*)
     write(*,*)
     write(*,*) "Starting search for coexistence beta for P = ",P
     write(*,*)


     ! The below loop is to find a pressure on the coexistence curve at the current P
     i=1
     do
        if(i>nr_max) then
           write(0,*) "monteswitch_npt: Error. We have performed 'nr_max' iterations in the Newton-Rhapson loop."
           stop 1
        end if

        write(*,*) 
        write(*,*) "On iteration ",i," for P = ",P
        write(*,*) 


        ! WEIGHT FUNCTION GENERATION SIMULATION

        ! Create a new simulation 'from scratch' for this beta and P, for weight function generation. Note that all simulations
        ! will start at a density of 'rho_init' and lattice 'init_lattice', with lists determined by 'list_cutoff'
        call initialise_lattices_bcc_fcc(rho_init,nx,ny,nz)
        call initialise_lists(list_cutoff)
        call initialise_cold_microstate(init_lattice)
        ! Initialise variables for weight function generation. Note that we create 'M_grid' and 'eta_grid'; and initialise
        ! the M_counts_1 and M_counts_2 arrays to be of the appropriate size via initialise_counters().
        enable_multicanonical=enable_multicanonical_generation
        stop_sweeps= stop_sweeps_generation
        update_eta=.true.
        update_trans= update_trans_generation
        enable_barriers= enable_barriers_generation
        calc_equil_properties= .false.
        call initialise_M_grid(M_grid_size,M_grid_min_kT/beta,M_grid_max_kT/beta)
        ! Error trap: Check if the current microstate is within the supported order parameter range
        if(M<M_grid(1) .or. M>=M_grid(M_grid_size)+(M_grid(2)-M_grid(1))) then
           write(0,*) "monteswitch_npt_coexistence: Error. Cold microstate outwith supported order parameter range."
           stop 1
        end if
        call initialise_counters()
        if(allocated(eta_grid)) then
           deallocate(eta_grid)
        end if
        allocate(eta_grid(M_grid_size))
        eta_grid=1.0_rk
        if(update_trans) then
           if(allocated(trans)) then
              deallocate(trans)
           end if
           allocate(trans(M_grid_size,M_grid_size))
           trans=0.0_rk
        end if
        if(enable_barriers) then
           call initialise_barriers()
        end if
        ! Export to 'state' for use by the 'shell_command' binary
        call export("state")
        ! Run the simulation to generate the weight function
        write(*,*) "Running simulation to generate weight function at (P,beta) = (",P,",",beta,")"
        call system(shell_command,error)
        if(error/=0) then
           write(*,*) "monteswitch_npt_coexistence: Error in executing shell command ",shell_command
           stop 1
        end if

        
        
        ! BATCH SIMULATION
        
        ! Modify the state file to correspond to a batch simulation
        call import("state")
        stop_sweeps=stop_sweeps_batch
        enable_multicanonical=.true.
        update_trans=.false.
        enable_barriers=.false.
        update_eta=.false.
        calc_equil_properties= .true.
        ! While we have the imported state output the weight function
        write(char_P,*) P
        write(char_beta,*) beta
        open(unit=11,file="wf_"//trim(adjustl(char_P))//"_"//trim(adjustl(char_beta)),status="replace")
        do j=1,M_grid_size
           write(11,*) M_grid(j),eta_grid(j)
        end do
        close(unit=11)
        ! Export the state file for use by the simulation
        call export("state")
        ! Run the batch simulation
        write(*,*) "Running batch simulation at (P,beta) = (",P,",",beta,")"
        call system(shell_command,error)
        if(error/=0) then
           write(*,*) "monteswitch_npt_coexistence: Error in executing shell command ",shell_command
           stop 1
        end if
        ! Append the 'examined points file'
        open(unit=11,file="examined_points",position="append")
        write(11,*) P,beta
        close(unit=11)
        ! Output the visited states
        open(unit=11,file="vs_"//trim(adjustl(char_P))//"_"//trim(adjustl(char_beta)),status="replace")
        do j=1,M_grid_size
           write(11,*) M_grid(j),M_counts_1(j),M_counts_2(j)
        end do
        close(unit=11)


        
        ! CHECK FOR CONVERGENCE (ON THE COEXISTENCE BETA) AND, IF NOT, AMMEND BETA AND TRY AGAIN

        ! Import the state variables for this beta and P
        call import("state")
        ! Export the state variables to the file 'state_P_BETA' for error trapping
        call export("state_"//trim(adjustl(char_P))//"_"//trim(adjustl(char_beta)))
        ! Print some useful information
        write(*,*) "State variables from this simulation:"
        write(*,*) "   (P,beta) = (",P,",",beta,")"
        write(*,*) "   equil_DeltaF = ",equil_DeltaF
        write(*,*) "   sigma_equil_DeltaF = ",sigma_equil_DeltaF
        write(*,*) "   block_counts_DeltaF = ",block_counts_DeltaF
        write(*,*) "   equil_H_1 = ",equil_H_1
        write(*,*) "   sigma_equil_H_1 = ",sigma_equil_H_1
        write(*,*) "   block_counts_H_1 = ",block_counts_H_1
        write(*,*) "   equil_H_2 = ",equil_H_2
        write(*,*) "   sigma_equil_H_2 = ",sigma_equil_H_2
        write(*,*) "   block_counts_H_2 = ",block_counts_H_2
        write(*,*) "   equil_V_1 = ",equil_V_1
        write(*,*) "   sigma_equil_V_1 = ",sigma_equil_V_1
        write(*,*) "   block_counts_V_1 = ",block_counts_V_1
        write(*,*) "   equil_V_2 = ",equil_V_2
        write(*,*) "   sigma_equil_V_2 = ",sigma_equil_V_2
        write(*,*) "   block_counts_V_2 = ",block_counts_V_2
        ! If any of the 'block_counts' variables are 0 then we cannot proceed, since the corresponding 'equil' variables
        ! are meaningless.
        if(block_counts_DeltaF<=0) then
           write(0,*) "monteswitch_npt_coexistence: Error. 'block_counts_DeltaF'<=0. Hence we have no 'equil_DeltaF' ", &
                "for this (P,beta) and cannot proceed."
           stop 1
        else if(block_counts_V_1<=0) then
           write(0,*) "monteswitch_npt_coexistence: Error. 'block_counts_V_1'<=0. Hence we have no 'equil_V_1' ", &
                "for this (P,beta) and cannot proceed."
           stop 1
        else if(block_counts_V_2<=0) then
           write(0,*) "monteswitch_npt_coexistence: Error. 'block_counts_V_2'<=0. Hence we have no 'equil_V_2' ", &
                "for this (P,beta) and cannot proceed."
           stop 1
        else if(block_counts_H_1<=0) then
           write(0,*) "monteswitch_npt_coexistence: Error. 'block_counts_H_1'<=0. Hence we have no 'equil_H_1' ", &
                "for this (P,beta) and cannot proceed."
           stop 1
        else if(block_counts_H_2<=0) then
           write(0,*) "monteswitch_npt_coexistence: Error. 'block_counts_H_2'<=0. Hence we have no 'equil_H_2' ", &
                "for this (P,beta) and cannot proceed."
           stop 1
        end if
        ! Calculate the change in beta predicted by Newton-Rhapson, and its uncertainty. Note that
        ! beta_new=beta-Delta_F(beta)/Delta_F'(beta), where recall that here F is the Gibbs free energy
        Delta_beta=-beta*equil_DeltaF/(equil_H_1-equil_H_2)
        sigma_H1H2_upper=abs(sigma_equil_H_1)+abs(sigma_equil_H_2)
        sigma_H1H2_lower=sqrt(sigma_equil_H_1**2+sigma_equil_H_2**2)
        sigma_Delta_beta_upper=abs(beta*equil_DeltaF/(equil_H_1-equil_H_2)) &
             * ( abs(sigma_equil_DeltaF/equil_DeltaF) + abs(sigma_H1H2_upper/(equil_H_1-equil_H_2)) )
        sigma_Delta_beta_lower=abs(beta*equil_DeltaF/(equil_H_1-equil_H_2)) &
             * sqrt( (sigma_equil_DeltaF/equil_DeltaF)**2 + (sigma_H1H2_lower/(equil_H_1-equil_H_2))**2 )
        ! Output this information
        write(*,*) "Newton-Rhapson predictions:"
        write(*,*) "   Delta_beta = ",Delta_beta
        write(*,*) "   sigma_Delta_beta_upper = ",sigma_Delta_beta_upper
        write(*,*) "   sigma_Delta_beta_lower = ",sigma_Delta_beta_lower
        ! Check for convergence: if Delta_beta<sigma_Delta_beta then the next beta is statistically indistinguishable
        ! from the current one, and we break the loop - we cannot get closer to the coexistence beta, which is
        ! what the Newton-Rhapson algorithm in this situation will converge upon (given the standard assumptions)
        if(abs(Delta_beta)<sigma_Delta_beta_upper) then
           write(*,*) "We have convergence on the coexistence beta. Exiting Newton-Rhapson loop."
           exit
        end if
        write(*,*) "No convergence on coexistence beta yet. Proceeding with loop."
        ! Otherwise update beta to correspond to the Newton-Rhapson prediction, and perform another simulation at this P.
        beta=beta+Delta_beta
        if(beta<0) then
           write(0,*) "monteswitch_npt_coexistence: Error. Estimated beta is < 0."
           stop 1
        end if
        write(*,*)
        ! Export the state variables to the file 'state' in preparation for the next simulation
        call export("state")


        i=i+1
     end do



     ! We have a beta at the current P which corresponds to a point on the coexistence curve. Output this to
     ! a file
     write(*,*) "Coexistence point found at (P,beta) = (",P,",",beta,")"
     write(*,*) "Uncertainty in beta is sigma_Delta_beta_upper = ",sigma_Delta_beta_upper
     open(unit=11,file="coexistence_points",position="append")
     write(11,*) P,beta,sigma_Delta_beta_upper
     close(unit=11)

     

     ! PREPARE FOR SIMULATIONS AT THE NEXT P

     ! Update P for the next loop. Guess beta according to a tangental extrapolation: beta_new=beta+Delta_P*dbeta/dP
     P=P+Delta_P
     beta=beta-Delta_P*beta*(equil_V_1-equil_V_2)/(equil_H_1-equil_H_2)
     write(*,*)
     write(*,*) "Proceeding to next P: P =",P
     write(*,*) "Estimated coexistence beta for this beta using tangental extrapolation is ",beta
     write(*,*) "(Using beta->beta-Delta_P*beta*(equil_V_1-equil_V_2)/(equil_H_1-equil_H_2)"
     if(beta<0) then
        write(0,*) "monteswitch_npt_coexistence: Error. Estimated beta is < 0."
        stop 1
     end if

  end do



contains



  ! This subroutine imports the variables contained within the file 'init_coexistence'
  subroutine import_init_coexistence()
    character(len=20) string

    ! Open the file
    open(unit=10,file="init_coexistence",iostat=error,status="old")
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem opening file 'init_coexistence'."
       stop 1
    end if

    ! Read the variables
    read(10,*,iostat=error) string, rho_init
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'rho_init' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 1 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, nx
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'nx' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 2 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, ny
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'ny' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 3 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, nz
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'nz' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 4 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, init_lattice
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'init_lattice' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 5 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, lj_epsilon
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'lj_epsilon' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 6 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, lj_sigma
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'lj_sigma' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 7 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, lj_cutoff
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'lj_cutoff' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 8 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, list_cutoff
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'list_cutoff' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 9 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, M_grid_size
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'M_grid_size' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 10 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, M_grid_min_kT
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'M_grid_min_kT' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 11 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, M_grid_max_kT
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'M_grid_max_kT' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 12 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, enable_multicanonical_generation
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'enable_multicanonical_generation' from ", &
            "file 'init_coexistence'. Note that this variable should be specified on line 13 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, P_min
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'P_min' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 14 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, P_max
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'P_max' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 15 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, Delta_P
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'Delta_P' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 16 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, beta_init
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'beta_init' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 17 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, part_select
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'part_select' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 18 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, part_step
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'part_step' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 19 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, enable_COM_frame
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'enable_COM_frame' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 20 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, vol_dynamics
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'vol_dynamics' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 21 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, vol_freq
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'vol_freq' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 22 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, vol_step
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'vol_step' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 23 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, stop_sweeps_generation
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'stop_sweeps_generation' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 24 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, stop_sweeps_batch
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'stop_sweeps_generation' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 25 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, enable_melt_checks
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'enable_melt_checks' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 26 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, melt_sweeps
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'melt_sweeps' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 27 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, melt_threshold
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'melt_threshold' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 28 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, melt_option
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'melt_option' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 29 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, enable_divergence_checks
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'enable_divergence_checks' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 30 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, divergence_sweeps
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'divergence_sweeps' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 31 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, divergence_tol
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'divergence_tol' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 32 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, update_eta_sweeps
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'update_eta_sweeps' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 33 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, update_trans_generation
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'update_trans_generation' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 34 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, update_eta_method
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'update_eta_method' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 35 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, enable_barriers_generation
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'enable_barriers_generation' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 36 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, barrier_dynamics
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'barrier_dynamics' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 37 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, lock_moves
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'lock_moves' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 38 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, block_sweeps
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'block_sweeps' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 39 of this file."
       stop 1
    end if 
    read(10,*,iostat=error) string, equil_sweeps
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'equil_sweeps' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 40 of this file."
       stop 1
    end if 
    read(10,*,iostat=error) string, nr_max
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'nr_max' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 41 of this file."
       stop 1
    end if
    read(10,*,iostat=error) string, shell_command
    if(error/=0) then
       write(0,*) "monteswitch_npt_coexistence: Error. Problem reading 'shell_command' from file 'init_coexistence'.", &
            " Note that this variable should be specified on line 42 of this file."
       stop 1
    end if

    close(unit=10)
  end subroutine import_init_coexistence




end program monteswitch_npt_coexistence
