!-----------------------------------------------------------------------
!
! Copyright (c) 2016 Tom L. Underwood
!
! Permission is hereby granted, free of charge, to any person obtaining 
! a copy of this software and associated documentation files (the 
! "Software"), to deal in the Software without restriction, including 
! without limitation the rights to use, copy, modify, merge, publish, 
! distribute, sublicense, and/or sell copies of the Software, and to 
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
! 
! The above copyright notice and this permission notice shall be 
! included in all copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
! BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
! ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
! SOFTWARE.
!
!-----------------------------------------------------------------------
!
! By extracting lines beginning with the regular expression '\!\! ?'
! (ignoring leading whitespace), and then removing matches to the
! regular expression, html documentation corresponding to this 
! source code will be created.
!
! By extracting lines beginning with the regular expression '\!init_params ?'
! (ignoring leading whitespace), and then removing matches to the
! regular expression, a template of an 'init_params' file will be created,
! i.e. a file corresponding to the 'filename_params' argument in the
! subroutine 'initialise_from_files(filename_params,filename_lattice)'.
! However, this template does not include the contributions from the
! 'energy_mod' module, which must be added to the end of the template.
!
!! <html>
!! <head>
!!  <title> monteswitch_mod Documentation </title>
!! </head>
!! <body>
!! 
!! <h1> <code>monteswitch_mod</code> Documentation </h1>
!!
!! <h2> Author </h2>
!! <p> Tom Underwood </p>
!!
!! <h2> General Notes </h2>
!!
!! <h3> Description </h3>
!! <p>
!! This module contains the nuts and bolts for lattice switch Monte Carlo simulations. The variables in this module define the 
!! 'state' of the simulation, which includes - among other things - variables describing the microstate of the
!! system at the current time and variables determining in what manner information will be output during a simulation. 
!! The variables associated with particle interactions are inherited from the module <code>interactions_mod</code>.
!! </p>
!! <p>
!! Particles are moved according to a random walk, with a maximum move of <code>part_step</code> in any Cartesian
!! direction. The volume is moved, if volume moves are used at all, in a manner according the flag
!! <code>vol_dynamics</code>. Details of the various options are provided below.
!! </p>
!! <p>
!! This module terminates the current program with a non-zero exit status of 1 for most errors. If the system has 
!! melted, and <code>melt_option="stop"</code> then the exit status is 2. The exit status is also 2 if the system has 'exploded',
!! i.e., it has negative volume. If the energy has diverged due to precision-related
!! errors then the exit status is 3. Note that exit statuses are not part of the Fortran standard, and may not work for
!! all operating systems or compilers.
!! </p>
!!
!! <h3> Definitions and conventions </h3>
!! <p>
!! <ul>
!!  <li>
!!  The free energy difference is the free energy of phase 1 minus that of phase 2.
!!  </li>
!!  <li>
!!  A "move" is one of either a: particle move, in which one particle is moved; a lattice move, in which
!!  the lattice is switched; and a volume move, in which the volume of the unit cell is altered. Note that
!!  if the move is rejected then it is still deemed to have taken place. For an NVT ensemble with lattice
!!  switch dynamics, the following move "cycle" is performed: particle move, lattice move. For an NPT 
!!  ensemble with lattice switch dynamics the move cycle is: particle move, lattice move, (volume move, 
!!  lattice move), where the set of moves in the brackets occur on average <code>vol_freq</code> times per 
!!  sweep. A different particle is used for each cycle; the particles are either cycled through one by one,
!!  or selected randomly for each cycle, according to the flag <code>part_select</code>. 
!!  If the flag <code>enable_lattice_moves</code> is set to <code>.false</code>
!!  then lattice switch dynamics are not used, and lattice moves are never attempted.
!!  </li>
!!  <li>
!!  A "sweep" comprises a cycle containing <code>n_part</code> particle moves, where <code>n_part</code> is the number
!!  of particles in the system. When the number of sweeps is printed to stdout
!!  or a file, specifically it is the number of <i>completed</i> sweeps. Tasks performed periodically on the scale of
!!  sweeps, e.g. updating the weight function, outputing data to files, checkpointing, are performed at the <i>end</i> of the
!!  sweep. The following table gives a breakdown of the number of moves of each type per sweep for different
!!  combinations of flags.
!!  <table border="1">
!!   <tr>
!!    <td> <b> Dynamics </b> </td>
!!    <td> <b> Particle moves </b> </td>
!!    <td> <b> Volume moves </b> </td>
!!    <td> <b> Lattice moves </b> </td>
!!    <td> <b> Total moves </b> </td>
!!   </tr>
!!   <tr>
!!    <td> 
!!    <code>enable_part_moves=.true.</code>
!!    <code>enable_vol_moves=.true.</code>
!!    <code>enable_lattice_moves=.true.</code>
!!    </td>
!!    <td> <code>n_part</code> </td>
!!    <td> <code>v_freq</code> (on average) </td>
!!    <td> <code>(n_part+vol_freq)</code> (on average)</td>
!!    <td> <code>2*(n_part+vol_freq)</code> (on average)</td>
!!   </tr> 
!!   <tr>
!!    <td> 
!!    <code>enable_part_moves=.true.</code>
!!    <code>enable_vol_moves=.false.</code>
!!    <code>enable_lattice_moves=.true.</code>
!!    </td>
!!    <td> <code>n_part</code> </td>
!!    <td> 0 </td>
!!    <td> <code>n_part</code> </td>
!!    <td> <code>2*n_part</code> </td>
!!   </tr>
!!   </tr>
!!   <tr>
!!    <td> 
!!    <code>enable_part_moves=.true.</code>
!!    <code>enable_vol_moves=.true.</code>
!!    <code>enable_lattice_moves=.false.</code>
!!    </td>
!!    <td> <code>n_part</code> </td>
!!    <td> <code>v_freq</code> (on average) </td>
!!    <td> 0 (on average)</td>
!!    <td> <code>(n_part+vol_freq)</code> (on average)</td>
!!   </tr> 
!!   </tr>
!!   <tr>
!!    <td> 
!!    <code>enable_part_moves=.true.</code>
!!    <code>enable_vol_moves=.false.</code>
!!    <code>enable_lattice_moves=.false.</code>
!!    </td>
!!    <td> <code>n_part</code> </td>
!!    <td> 0 </td>
!!    <td> 0 </td>
!!    <td> <code>n_part</code></td>
!!   </tr> 
!!   <tr>
!!    <td> 
!!    <code>enable_part_moves=.false.</code>
!!    <code>enable_vol_moves=.false.</code>
!!    <code>enable_lattice_moves=.true.</code>
!!    </td>
!!    <td> 0 </td>
!!    <td> 0 </td>
!!    <td> <code>n_part</code></td>
!!    <td> <code>n_part</code></td>
!!   </tr> 
!!   <tr>
!!    <td> 
!!    <code>enable_part_moves=.false.</code>
!!    <code>enable_vol_moves=.true.</code>
!!    <code>enable_lattice_moves=.true.</code>
!!    </td>
!!    <td> 0 </td>
!!    <td> <code>v_freq</code> (on average) </td>
!!    <td> <code>(n_part+vol_freq)</code> (on average)</td>
!!    <td> <code>(n_part+vol_freq)</code> (on average)</td>
!!   </tr> 
!!  </table>
!!  </li>
!!  <li>
!!  The order parameter for a microstate is defined as <code>E_1-E_2+P*(V_1-V_2)-log(V_1/V_2)/beta</code>, 
!!  where <code>E_1</code> is the energy of the lattice-type-1 microstate corresponding to the current displacements, 
!!  <code>E_2</code> is the energy of the lattice-type-2 microstate, <code>V_1</code> is the volume of the lattice-type-1
!!  microstate, <code>V_2</code> is the volume of the lattice-type-2 microstate, and <code>P</code> and <code>beta</code>
!!  have their usual significance.
!!  </li>
!!  <li>
!!  The weight function is such that the probability of microstate 'sigma' being sampled using multicanonical
!!  sampling is exp(-beta*E(sigma))exp(eta(sigma)), where eta(sigma) is the value of the weight fuunction for
!!  sigma.
!!  </li>
!! </ul>
!! </p>
!!
!! <h3> Warnings </h3>
!! <p>
!! <ul>
!!  <li>
!!  The program may crash without any error messages if the Fortran intrinsic function <code>mod(a,p)</code> is
!!  called with <code>p=0</code>. This could occur if any of the variables which describe a 'period' in
!!  sweeps is 0. Therefore do not set any of these to 0 if they are in use. Depending on how I've coded, it
!!  may even be the case that the program crashes if any of these are 0 and they aren't 'in use'. I will have
!!  to test this at some point; to be safe just don't use 0 for such variables.
!!  </li>
!!  <li>
!!  There are no default values for any of the variables. It is therefore recommended that their values are all
!!  explicitly stated at initialisation: one would definitely not want a strange value of, say, 
!!  <code>stop_sweeps</code> creeping into the simulation because its value wasn't initialised.
!!  </li>
!!  <li>
!!  In a similar vein to the above point, if one, say, tells the simulation to not evaluate equilibrium
!!  quantities, then the values of the variables associated with equilibrium quantities will be nonsense, and
!!  do not imply a bug in the program.
!!  </li>
!! </ul>
!! </p>
!!
module monteswitch_mod

  !! <h2> Dependencies </h2>
  !! <p> 
  !! <ul>
  !!  <li> <code> kinds_mod </code> </li>
  !!  <li> <code> rng_mod </code> </li>
  !!  <li> <code> metropolis_mod </code> </li>
  !!  <li> <code> interactions_mod </code> </li>
  !! </ul>
  !! </p>
  use kinds_mod
  use rng_mod
  use metropolis_mod
  use interactions_mod

  implicit none

  !! <h2> Variables </h2>
  
  !! <p>
  !! What follows is a list of the module variables. Note also that variables are inherited from <code>interactions_mod</code>
  !! which contribute to the specification of the state of the system, such as neighbour lists. 
  !! Variables can be loosely categorised as either user-defined or 'internal'. In theory a user need never tamper with 
  !! internal variables - and I would not recommend it unless they <i>really</i> know what they are doing. Internal 
  !! variables' descriptions are in red.
  !! </p>
  !!
  !! <h3> Variables determining the ensemble and dynamics</h3>
  !! <table border="1">
  !! <tr>
  !!  <td> <b> Variable </b> </td>
  !!  <td> <b> Type </b> </td>
  !!  <td> <b> Description </b> </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>enable_multicanonical<code> </td>
  !!  <td> <code>logical</code> </td>
  !!  <td> 
  !!  Flag which determines whether multicanonical sampling is used. If set to <code>.true.</code> then
  !!  multicanonical sampling is used via the weight function <code>eta</code>; if set to <code>.false</code>
  !!  then Boltzmann sampling is used.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>beta<code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> 
  !!  The thermodynamic beta. This should be greater than zero - setting this to zero will yield division by
  !!  zero errors in some quantities, e.g. the free energy difference between phases
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>P<code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> Pressure of the system. </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>enable_lattice_moves<code> </td>
  !!  <td> <code>logical</code> </td>
  !!  <td> 
  !!  Flag which which determines whether lattice switch dynamics is used. If set to <code>.true.</code>,
  !!  then after every particle or volume move a lattice switch is attempted. If set to
  !!  <code>.false.</code> then conventional dynamics are used: the underlying lattice is never changed
  !!  during the course of a simulation.
  !!  </td>
  !! <tr>
  !!  <td> <code>enable_part_moves<code> </td>
  !!  <td> <code>logical</code> </td>
  !!  <td> 
  !!  Flag which which determines whether particles will be moved: if set to <code>.true.</code>,
  !!  then particles will be moved. This exists mainly for debugging purposes.
  !!  </td>
  !! <tr>
  !!  <td> <code>enable_vol_moves<code> </td>
  !!  <td> <code>logical</code> </td>
  !!  <td> 
  !!  Flag which which determines whether volume moves are used, and determines the ensemble to be used:
  !!  if set to <code>.true.</code> then we have an NPT ensemble; if set to <code>.false.</code> then we
  !!  have an NVT ensemble. Note that the Metropolis algorithm differs slightly between NVT and NPT
  !!  ensembles.
  !!  </td> 
  !! <tr>
  !!  <td> <code>part_select<code> </td>
  !!  <td> <code>character(len=30)</code> </td>
  !!  <td> 
  !!  Flag which which determines how a particle is to be selected at the
  !!  beginning of each cycle of moves. Possible values are as follows:
  !!  <code>"cycle"</code> results in selecting particles in order, i.e. particle 
  !!  1 is used first, then particle 2,..., then particle <code>n_part</code>, then 
  !!  particle 1, etc.; <code>"rand"</code> selects a particle at the beginning of each
  !!  cycle at random.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>part_step<code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> 
  !!  Maximum displacement change in any Cartesian direction which any particle can be
  !!  moved by during a particle move. Recall that random walk particle dynamics is used.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>enable_COM_frame<code> </td>
  !!  <td> <code>logical</code> </td>
  !!  <td> 
  !!  Flag which determines whether we use the centre-of-mass frame. If set to <code>.true.</code>
  !!  the particle displacements are ammended after every accepted particle move such that the
  !!  average displacement is 0.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>vol_dynamics<code> </td>
  !!  <td> <code>character(len=30)</code> </td>
  !!  <td> 
  !!  Flag which determines how the supercell is altered for a volume move. Note
  !!  that this is only important for NPT ensembles. Possible values are as follows:
  !!  <code>"FVM"</code> uses a random walk in ln(V) with a maximum step size of <code>vol_step<code>,
  !!  and a fixed aspect ratio, i.e. the proportions of the supercell dimensions are unchanged by a 
  !!  volume move; <code>"UVM"</code> uses a random walk to expand/contract each dimension of the supercell. 
  !!  </td>
  !! <tr>
  !!  <td> <code>vol_freq<code> </td>
  !!  <td> <code>integer(ik)</code> </td>
  !!  <td> 
  !!  The number of volume moves which will occur per sweep on average.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>vol_step<code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td>
  !!  Maximum step size for a volume move random walk. The volume can change by a maximum of <code>vol_step<code>
  !!  each volume move.
  !!  </td>
  !! </tr>
  !! </table>
  logical :: enable_multicanonical
  real(rk) :: beta
  real(rk) :: P
  logical :: enable_lattice_moves
  logical :: enable_part_moves
  logical :: enable_vol_moves
  character(len=30) :: part_select
  character(len=*), parameter :: part_select_cycle="cycle"
  character(len=*), parameter :: part_select_rand="rand"
  real(rk) :: part_step
  logical :: enable_COM_frame
  character(len=30) :: vol_dynamics
  character(len=*), parameter :: vol_dynamics_FVM="FVM"
  character(len=*), parameter :: vol_dynamics_UVM="UVM"
  integer(ik) :: vol_freq
  real(rk) :: vol_step

  !! <h3> Variables determining the length of the simulation and the equilibration time </h3>
  !! <p>
  !! <table border="1">
  !! <tr>
  !!  <td> <b> Variable </b> </td>
  !!  <td> <b> Type </b> </td>
  !!  <td> <b> Description </b> </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>stop_sweeps</code> </td>
  !!  <td> <code>integer(ik)</code> </td>
  !!  <td>
  !!  The number of sweeps to perform during this simulation. If this is set to 0 then the Monte Carlo loop
  !!  is by-passed, but tasks performed periodically during the Monte Carlo loop, e.g. updating the weight
  !!  function and checking whether or not the system has melted/exploded, are done once. Furthermore, if it is set to
  !!  0, the equilibrium properties - the free energy difference between phases, the volume and enthalpy of 
  !!  each phase, and the associated uncertainties of these quantities - are recalculated if they are to be 
  !!  calculated.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>equil_sweeps</code> </td>
  !!  <td> <code>integer(ik)</code> </td>
  !!  <td> The number of sweeps after which equilibration is assumed to have taken place. </td>
  !! </tr>
  !! </table>
  integer(ik) :: stop_sweeps
  integer(ik) :: equil_sweeps

  !! <h3> Variables pertaining to safety checks on the simulation </h3>
  !! <p>
  !! <table border="1">
  !! <tr>
  !!  <td> <b> Variable </b> </td>
  !!  <td> <b> Type </b> </td>
  !!  <td> <b> Description </b> </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>enable_melt_checks</code> </td>
  !!  <td> <code>logical</code> </td>
  !!  <td> 
  !!  Flag determining whether the simulation will periodically check if the particle displacements have
  !!  exceeded 'safe' values, i.e. if the system has melted, and if the system has a negative volume, i.e.,
  !!  the system has 'exploded'.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>melt_sweeps</code> </td>
  !!  <td> <code>integer(ik)</code> </td>
  !!  <td> 
  !!  If <code>enable_melt_checks=.true.</code>, the particle displacements are checked for melting/exploding every
  !!  <code>melt_sweeps</code> sweeps.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>melt_threshold</code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> 
  !!  If any particle displacement is of magnitude <code>>melt_threshold</code> in any Cartesian direction then the system
  !!  is deemed to have melted.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>melt_option</code> </td>
  !!  <td> <code>character(len=30)</code> </td>
  !!  <td> 
  !!  Flag which determines what happens if the system melts. Options are as follows. <code>"zero_1"</code> and 
  !!  <code>"zero_2"</code> move the system to the zero displacement microstate at the current volume in lattice types 1
  !!  and 2 respectively, and then proceeds. <code>"zero_random"</code> does the same but picks a random lattice type, and
  !!  <code>"zero_current"</code> picks the current lattice.
  !!  In all of these cases, the system waits for <code>equil_sweeps</code> after the displacements are 'reset' for the
  !!  system to equilibrate before variables which should not be updated before equilibration are updated. Also, the current
  !!  block with regards to the calculation of equilibrium quantities is disregarded.
  !!  <code>"stop"</code> prints an error message to stderr if the system melts and returns an exit status of 2.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>enable_divergence_checks</code> </td>
  !!  <td> <code>logical</code> </td>
  !!  <td> 
  !!  Flag determining whether the simulation will periodically check if the energy is correct. For particle moves, the energy
  !!  is not calculated exactly for the trial particle displacements because it is very computationally expensive and unnecessary.
  !!  Instead the energy <i>change</i> for the trial displacements is calculated, which involves consideration of the energy 
  !!  stored in the interactions between the particle which has been moved, and its neighbours. This is far less demanding to
  !!  calculate. If the move is accepted this change is ammended to the total energy. However, over time it is possible that this
  !!  'running total' approach will yield incorrect energies due to the finite precision of the computer. Hence one should periodically
  !!  recalculate the energy exactly. This is done during volume moves. If <code>enable_divergence_checks=.true.</code> this is also
  !!  done every <code>divergence_sweeps</code> sweeps, after which the recalculated energy is compared to the 'current' energy, and
  !!  an error is flagged if they are different - outwith a tolerence of <code>divergence_tol</code>. Note that the order parameter
  !!  is also ammended after the energy (for lattice 1 and 2) is recalculated.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>divergence_sweeps</code> </td>
  !!  <td> <code>integer(ik)</code> </td>
  !!  <td> 
  !!  If <code>enable_divergence_checks=.true.</code>, the energy is checked every <code>divergence_sweeps</code> sweeps.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>divergence_tol</code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> 
  !!  The difference between the current energy and the 'true' energy of the system which would count as an error. This should
  !!  be positive.
  !!  </td>
  !! </tr>
  !! </table>
  logical :: enable_melt_checks
  integer(ik) :: melt_sweeps
  real(rk) :: melt_threshold
  character(len=30) :: melt_option
  character(len=30), parameter :: melt_option_zero_1="zero_1"
  character(len=30), parameter :: melt_option_zero_2="zero_2"
  character(len=30), parameter :: melt_option_zero_random="zero_random"
  character(len=30), parameter :: melt_option_zero_current="zero_current"
  character(len=30), parameter :: melt_option_stop="stop"
  logical :: enable_divergence_checks
  integer(ik) :: divergence_sweeps
  real(rk) :: divergence_tol

  !! <h3> Variables determining the nature of the output </h3>
  !! <p>
  !! The flags <code>output_file_X</code>, where <code>X</code> is the name of a simulation variable, when set to <code>.true</code>,
  !! will result in a line consisting of "<code>X: </code>" followed by the number of completed sweeps, followed by the current value of
  !! <code>X</code> being printed to a specified file (<code>datafile</code> in the <code>run</code> subroutine) every 
  !! <code>output_file_period</code> sweeps (at the beginning of the sweep) - where the file is an argument of the <code>run</code> 
  !! subroutine. However, if <code>output_file_period=0</code>, then instead the output is <i>after every move</i>; and if 
  !! <code>output_file_period<0</code> then there is no output to the file. Similar applies to <code>output_stdout_X</code>, but for 
  !! stdout during a simulation.
  !! </p>
  !! <p>
  !! The variable <code>checkpoint_period</code> determines how often (in sweeps) the simulation is checkpointed, i.e., how often the 
  !! state of the simulation is exported to a specified file (<code>statefile</code> in the <code>run</code> subroutine). Note that 
  !! the simulation is automatically checkpointed at the end of every simulation. If <code>checkpoint_period<=0</code>, then there is
  !! no output to the file <code>statefile</code>, i.e., it is empty. If one wants the simulation to be checkpointed only at the end
  !! of a simulation, and not during it, then set <code>checkpoint_period>stop_sweeps</code>.
  !! </p>
  !! <table border="1">
  !! <tr>
  !!  <td> <b> Variable </b> </td>
  !!  <td> <b> Type </b> </td>
  !! </tr>
  !! <tr> <td> <code>output_file_period</code> </td> <td> <code>integer(ik)</code> </td> </tr>  
  !! <tr> <td> <code>output_file_Lx </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_Ly </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_Lz </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_V </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_R_1 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_R_2 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_u </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_lattice </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_E </code> </td> <td> <code>logical</code> </td>  </tr>
  !! <tr> <td> <code>output_file_M </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_macro </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_eta </code> </td> <td> <code>logical</code> </td> </tr>  
  !! <tr> <td> <code>output_file_moves_lattice </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_accepted_moves_lattice </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_moves_part </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_accepted_moves_part </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_moves_vol </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_accepted_moves_vol </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_rejected_moves_M_OOB </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_M_OOB_high </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_barrier_macro_low </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_barrier_macro_high </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_rejected_moves_M_barrier </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_moves_since_lock </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_melts </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_equil_DeltaFs </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_sigma_equil_DeltaFs </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_equil_H_1 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_sigma_equil_H_1 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_equil_H_2 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_sigma_equil_H_2 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_equil_V_1 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_sigma_equil_V_1 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_equil_V_2 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_sigma_equil_V_2 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_equil_umsd_1 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_sigma_equil_umsd_1 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_equil_umsd_2 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_sigma_equil_umsd_2 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_equil_L_1 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_sigma_equil_L_1 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_equil_L_2 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_file_sigma_equil_L_2 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_period</code> </td> <td> <code>integer(ik)</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_Lx </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_Ly </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_Lz </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_V </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_R_1 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_R_2 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_u </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_lattice </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_E </code> </td> <td> <code>logical</code> </td>  </tr>
  !! <tr> <td> <code>output_stdout_M </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_macro </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_eta </code> </td> <td> <code>logical</code> </td> </tr>  
  !! <tr> <td> <code>output_stdout_accepted_moves_lattice </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_accepted_moves_part </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_accepted_moves_vol </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_rejected_moves_M_OOB </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_M_OOB_high </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_barrier_macro_low </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_barrier_macro_high </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_rejected_moves_M_barrier </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_moves_since_lock </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_equil_DeltaFs </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_sigma_equil_DeltaFs </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_equil_H_1s </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_sigma_equil_H_1s </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_equil_H_2s </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_sigma_equil_H_2s </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_equil_V_1s </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_sigma_equil_V_1s </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_equil_V_2s </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_sigma_equil_V_2s </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_melts </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_equil_umsd_1 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_sigma_equil_umsd_1 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_equil_umsd_2 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_sigma_equil_umsd_2 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_equil_L_1 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_sigma_equil_L_1 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_equil_L_2 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>output_stdout_sigma_equil_L_2 </code> </td> <td> <code>logical</code> </td> </tr>
  !! <tr> <td> <code>> <code>checkpoint_period </code> </td> <td> <code>integer(ik)</code> </td> </tr> </code> </td> <td> <code>logical</code> </td> </tr>
  !! </table>
  integer(ik) :: output_file_period
  logical :: output_file_Lx
  logical :: output_file_Ly
  logical :: output_file_Lz
  logical :: output_file_V
  logical :: output_file_R_1
  logical :: output_file_R_2
  logical :: output_file_u
  logical :: output_file_lattice
  logical :: output_file_E
  logical :: output_file_M
  logical :: output_file_macro
  logical :: output_file_eta
  logical :: output_file_moves_lattice
  logical :: output_file_accepted_moves_lattice
  logical :: output_file_moves_part
  logical :: output_file_accepted_moves_part
  logical :: output_file_moves_vol
  logical :: output_file_accepted_moves_vol
  logical :: output_file_rejected_moves_M_OOB
  logical :: output_file_M_OOB_high
  logical :: output_file_M_OOB_low
  logical :: output_file_barrier_macro_low
  logical :: output_file_barrier_macro_high
  logical :: output_file_rejected_moves_M_barrier
  logical :: output_file_moves_since_lock
  logical :: output_file_melts
  logical :: output_file_equil_DeltaF
  logical :: output_file_sigma_equil_DeltaF
  logical :: output_file_equil_H_1
  logical :: output_file_sigma_equil_H_1
  logical :: output_file_equil_H_2
  logical :: output_file_sigma_equil_H_2
  logical :: output_file_equil_V_1
  logical :: output_file_sigma_equil_V_1
  logical :: output_file_equil_V_2
  logical :: output_file_sigma_equil_V_2
  logical :: output_file_equil_umsd_1
  logical :: output_file_sigma_equil_umsd_1
  logical :: output_file_equil_umsd_2
  logical :: output_file_sigma_equil_umsd_2
  logical :: output_file_equil_L_1
  logical :: output_file_sigma_equil_L_1
  logical :: output_file_equil_L_2
  logical :: output_file_sigma_equil_L_2
  integer(ik) :: output_stdout_period
  logical :: output_stdout_Lx
  logical :: output_stdout_Ly
  logical :: output_stdout_Lz
  logical :: output_stdout_V
  logical :: output_stdout_R_1
  logical :: output_stdout_R_2
  logical :: output_stdout_u
  logical :: output_stdout_lattice
  logical :: output_stdout_E
  logical :: output_stdout_M
  logical :: output_stdout_macro
  logical :: output_stdout_eta
  logical :: output_stdout_moves_lattice
  logical :: output_stdout_accepted_moves_lattice
  logical :: output_stdout_moves_part
  logical :: output_stdout_accepted_moves_part
  logical :: output_stdout_moves_vol
  logical :: output_stdout_accepted_moves_vol
  logical :: output_stdout_rejected_moves_M_OOB
  logical :: output_stdout_M_OOB_high
  logical :: output_stdout_M_OOB_low
  logical :: output_stdout_barrier_macro_low
  logical :: output_stdout_barrier_macro_high
  logical :: output_stdout_rejected_moves_M_barrier
  logical :: output_stdout_moves_since_lock
  logical :: output_stdout_melts
  logical :: output_stdout_equil_DeltaF
  logical :: output_stdout_sigma_equil_DeltaF
  logical :: output_stdout_equil_H_1
  logical :: output_stdout_sigma_equil_H_1
  logical :: output_stdout_equil_H_2
  logical :: output_stdout_sigma_equil_H_2
  logical :: output_stdout_equil_V_1
  logical :: output_stdout_sigma_equil_V_1
  logical :: output_stdout_equil_V_2
  logical :: output_stdout_sigma_equil_V_2
  logical :: output_stdout_equil_umsd_1
  logical :: output_stdout_sigma_equil_umsd_1
  logical :: output_stdout_equil_umsd_2
  logical :: output_stdout_sigma_equil_umsd_2
  logical :: output_stdout_equil_L_1
  logical :: output_stdout_sigma_equil_L_1
  logical :: output_stdout_equil_L_2
  logical :: output_stdout_sigma_equil_L_2
  integer(ik) :: checkpoint_period


  !! <h3> Variables defining macrostates </h3>
  !! <p>
  !! In a simulation we need to define macrostates, each of which consists of a range of order parameters. 
  !! <code>M_grid</code> is an array containing <code>M_grid_size</code> equidistant
  !! values. We emphasise that <code>M_grid</code> must have this form. 
  !! Furthermore, <code>M_grid_size</code> must be greater than 1.
  !! <code>M_grid</code> defines macrostates: the 'i'th macrostate is defined as consisting of microstates
  !! whose order parameter is between <code>M_grid(i)</code> (inclusive) and <code>M_grid(i)+(M_grid(2)-M_grid(1))</code> 
  !! (exclusive) - with <code>(M_grid(2)-M_grid(1))</code> being the spacing between successive elements of <code>M_grid</code>. Moves 
  !! which take us to states with <code>M>=M_grid(M_grid_size)+(M_grid(2)-M_grid(1))</code> or <code>M&ltM_grid(1)</code>
  !! are always rejected.
  !! </p>
  !! <table border="1">
  !! <tr>
  !!  <td> <b> Variable </b> </td>
  !!  <td> <b> Type </b> </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>M_grid_size</code> </td>
  !!  <td> <code>integer(ik)</code> </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>M_grid</code> </td>
  !!  <td> <code> real(rk), dimension(:), allocatable</code> </td>
  !! </tr>
  !! </table>
  integer(ik) :: M_grid_size
  real(rk), dimension(:), allocatable :: M_grid

  !! <h3> Variables describing the state of the system </h3>
  !! <table border="1">
  !! <tr>
  !!  <td> <b> Variable </b> </td>
  !!  <td> <b> Type </b> </td>
  !!  <td> <b> Description </b> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>n_part</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  Total number of particles in the system. </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>Lx</code> </font> </td>
  !!  <td> <font color="red">  <code> real(rk), dimension(2)</code> </font> </td>
  !!  <td> <font color="red">  Dimension of supercell in x direction for lattice types 1 and 2. </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>Ly</code> </font> </td>
  !!  <td> <font color="red">  <code> real(rk), dimension(2)</code> </font> </td>
  !!  <td> <font color="red">  Dimension of supercell in y direction for lattice types 1 and 2. </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>Lz</code> </font> </td>
  !!  <td> <font color="red">  <code> real(rk), dimension(2)</code> </font> </td>
  !!  <td> <font color="red">  Dimension of supercell in z direction for lattice types 1 and 2. </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>V</code> </font> </td>
  !!  <td> <font color="red">  <code> real(rk)</code> </font> </td>
  !!  <td> <font color="red">  Current volume of the system (for NPT simulations) </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>spec_1</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik), dimension(:), allocatable</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  <code>spec_1(n)</code> is the species of the nth particle for lattice type 1.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>spec_2</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik), dimension(:), allocatable</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  <code>spec_2(n)</code> is the species of the nth particle for lattice type 2.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>pos_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(:,:), allocatable</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  (<code>pos_1(1,n)</code>,<code>pos_1(2,n)</code>,<code>pos_1(3,n)</code>) is the position 
  !!  of the nth particle for lattice type 1.
  !!  This array should have <code>n_part</code> elements in its first dimension,
  !!  and 3 elements in its second.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>pos_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(:,:), allocatable</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  (<code>pos_2(1,n)</code>,<code>pos_2(2,n)</code>,<code>pos_2(3,n)</code>) is the position 
  !!  of the nth particle for lattice type 2.
  !!  This array should have <code>n_part</code> elements in its first dimension,
  !!  and 3 elements in its second.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>R_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(:,:), allocatable</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  (<code>R_1(1,n)</code>,<code>R_1(2,n)</code>,<code>R_1(3,n)</code>) is the position 
  !!  of the nth particle with zero displacement for lattice type 1.
  !!  This array should have <code>n_part</code> elements in its first dimension,
  !!  and 3 elements in its second.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>R_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(:,:), allocatable</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  (<code>R_2(1,n)</code>,<code>R_2(2,n)</code>,<code>R_2(3,n)</code>) is the position 
  !!  of the nth particle with zero displacement for lattice type 2.
  !!  This array should have <code>n_part</code> elements in its first dimension,
  !!  and 3 elements in its second.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>u</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(:,:), allocatable</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  (<code>u(1,n)</code>,<code>u(2,n)</code>,<code>u(3,n)</code>) is the current displacement vector
  !!  for the nth particle in the system, i.e., if the current lattice
  !!  is lattice type 1, then <code>u(:,n)+R_1(:,n)</code> are the current positions
  !!  of the particles in the system. This array should have <code>n_part</code>
  !!  elements in its first dimension, and 3 elements in its second.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>lattice</code> </font> </td>
  !!  <td> <font color="red">  <code>integer</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  The current lattice type (1 or 2).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>E_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  The energy associated with lattice type 1 for the current set of
  !!  displacements and volume.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>E_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  The energy associated with lattice type 2 for the current set of
  !!  displacements and volume.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>E</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  The current energy. (<code>E=E_1</code> if <code>lattice=1</code> and <code>E=E_2</code> if <code>lattice=2</code>).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>M</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  The current order parameter.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>macro</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  The current macrostate.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>eta</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  The value of the weight function for the current state.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>M_const</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  The value of <code>log(V_1/V_2)/beta</code> determined at the beginning of the simulation, where
  !!  <code>V_1</code> is the volume of box 1. Note that this is a constant throughout a simulation since <code>V_1/V_2</code> 
  !!  is preserved by the lattice switch. This quantity is stored in memory for optimisation purposes - 
  !!  to avoid calling the log function every move.
  !!  </font> </td>
  !! </tr>
  !! </table>
  integer(ik) :: n_part
  real(rk), dimension(2) :: Lx
  real(rk), dimension(2) :: Ly
  real(rk), dimension(2) :: Lz
  real(rk) :: V
  integer(ik), dimension(:), allocatable :: spec_1
  integer(ik), dimension(:), allocatable :: spec_2
  real(rk), dimension(:,:), allocatable :: pos_1
  real(rk), dimension(:,:), allocatable :: pos_2
  real(rk), dimension(:,:), allocatable :: R_1
  real(rk), dimension(:,:), allocatable :: R_2
  real(rk), dimension(:,:), allocatable :: u
  integer :: lattice
  real(rk) :: E_1
  real(rk) :: E_2
  real(rk) :: E
  real(rk) :: M
  integer(ik) :: macro
  real(rk) :: eta
  real(rk) :: M_const

  !! <h3> Variables describing the lattice switch </h3>
  !! <table border="1">
  !! <tr>
  !!  <td> <b> Variable </b> </td>
  !!  <td> <b> Type </b> </td>
  !!  <td> <b> Description </b> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>E</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  The current energy. (<code>E=E_1</code> if <code>lattice=1</code> and <code>E=E_2</code> if <code>lattice=2</code>).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>switchscalex</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  The scalefactor of the supercell in the x-dimension to take the system from lattice 1 to lattice 2
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>switchscaley</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  The scalefactor of the supercell in the y-dimension to take the system from lattice 1 to lattice 2
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>switchscalez</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red"> 
  !!  The scalefactor of the supercell in the z-dimension to take the system from lattice 1 to lattice 2
  !!  </font> </td>
  !! </tr>
  real(rk) :: switchscalex
  real(rk) :: switchscaley
  real(rk) :: switchscalez

  !! <h3> Variables keeping track of move numbers, acceptance rates, etc. </h3>
  !! <table border="1">
  !! <tr>
  !!  <td> <b> Variable </b> </td>
  !!  <td> <b> Type </b> </td>
  !!  <td> <b> Description </b> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>M_counts_1</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik), dimension(:), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  <code>M_counts_1(i)</code> is the number of microstates visited which belong to
  !!   macrostate 'i' and lattice 1. This array must have <code>M_grid_size</code> elements.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>M_counts_2</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik), dimension(:), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  <code>M_counts_2(i)</code> is the number of microstates visited which belong to
  !!   macrostate 'i' and lattice 2. This array must have <code>M_grid_size</code> elements.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>sweeps</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Total number of sweeps which have been made.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>moves</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Total number of moves which have been made.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>moves_lattice</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Total number of lattice moves which have been made.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>accepted_moves_lattice</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Total number of accepted lattice switch moves.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>moves_part</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Total number of particle moves which have been made.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>accepted_moves_part</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Total number of accepted particle moves.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>moves_vol</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Total number of volume moves which have been made.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>accepted_moves_vol</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Total number of accepted volume moves.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>rejected_moves_M_OOB</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  The number of moves rejected because the order parameter was outwith the
  !!  order parameter range supported by <code>M_grid</code>, i.e. it is the number moves 
  !!  rejected because <code>M>=M_grid(M_grid_size)+(M_grid(2)-M_grid(1))</code> or
  !!  <code>M<M_grid(1)</code>. (OOB = out of bounds).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>M_OOB_high</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  The highest order parameter encountered which is outwith the supported range.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>M_OOB_low</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  The lowest order parameter encountered which is outwith the supported range.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>melts</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  The number of times the system has melted.
  !!  </font> </td>
  !! </tr>
  !! </table>
  integer(ik), dimension(:), allocatable :: M_counts_1
  integer(ik), dimension(:), allocatable :: M_counts_2
  integer(ik) :: sweeps
  integer(ik) :: moves
  integer(ik) :: moves_lattice
  integer(ik) :: accepted_moves_lattice
  integer(ik) :: moves_part
  integer(ik) :: accepted_moves_part
  integer(ik) :: moves_vol
  integer(ik) :: accepted_moves_vol
  integer(ik) :: rejected_moves_M_OOB
  real(rk) :: M_OOB_high
  real(rk) :: M_OOB_low
  integer(ik) :: melts

  !! <h3> Variables pertaining to the weight function and weight function generation </h3>
  !! <table border="1">
  !! <tr>
  !!  <td> <b> Variable </b> </td>
  !!  <td> <b> Type </b> </td>
  !!  <td> <b> Description </b> </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>eta_grid</code> </td>
  !!  <td> <code>real(rk), dimension(:), allocatable</code> </td>
  !!  <td> 
  !!  <code>eta_grid(i)</code> is the value of the weight function for the 'i'th macrostate.
  !!  <code>eta_grid</code> must have <code>M_grid_size</code> elements.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>trans</code> </td>
  !!  <td> <code>real(rk), dimension(:,:), allocatable</code> </td>
  !!  <td> 
  !!  <code>trans(1,i)</code>, <code>trans(2,i)</code>, <code>trans(3,i)</code> and <code>trans(4,i)</code>
  !!  are, respectively, the number of infered transitions from macrostate 'i' to macrostate 'i-1', 'i',
  !!  'i+1', and any other macrostate.  <code>trans</code> must have 4 elements in its first and
  !!  <code>M_grid_size</code> elements its second dimension.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>update_eta</code> </td>
  !!  <td> <code>logical</td>
  !!  <td> 
  !!  Flag determining whether the weight function is periodically updated or not. If set to <code>.true.</code>
  !!  then the weight function is updated every <code>update_eta_sweeps</code> sweeps. 
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>update_eta_sweeps</code> </td>
  !!  <td> <code>integer(ik)</td>
  !!  <td> 
  !!  The period of sweeps at which the weight function is updated if it is updated.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>update_eta_method<code> </td>
  !!  <td> <code>character(len=30)</code> </td>
  !!  <td> 
  !!  Flag which determines how the weight function will be updated if it is updated.
  !!  Possible values are as follows:
  !!  <ul>
  !!     <li>
  !!     <code>"VS"</code> uses the visited states method. Here the weight function is updated every 
  !!     <code>update_eta_sweeps</code> sweeps: it uses the order parameter histogram (comprised of the histograms
  !!     <code>M_counts_1</code> and <code>M_counts_2</code> pertaining to both phases) to estimate the
  !!     'ideal weight function' which would result in all order parameter macrostates being
  !!     visited with equal probability with multicanonical sampling. See <code>update_eta_VS()</code> for further
  !!     details. After the weight function is updated,
  !!     <code>M_counts_1</code> and <code>M_counts_2</code> are reset to zero for the next 'block' of 
  !!     <code>update_eta_sweeps</code> sweeps. During this block, the weight function will be clcoser to
  !!     the ideal than it was in the previous block, and hence the order parameter histogram will be 'flatter'
  !!     - since the macrostates will have been visited more uniformly. After this block the weight function
  !!     is recalculated using the 'new' order parameter histogram, the order parameter histograms are reset,
  !!     etc. Eventually the weight function will converge upon the ideal one. Note that this method is not
  !!     suitable for parallelisation, and that if <code>stop_sweeps=0</code> then the weight function is
  !!     immediately recalculated, and the order parameter histograms are reset. Also, be sure not to set 
  !!     <code>update_eta_sweeps</code> sweeps to be too low, in which case the order parameter histogram 
  !!     for each block will be too sparse to be meaningful.
  !!     </li>
  !!     <li>
  !!     <code>"shooting"</code> uses the shooting method. Here <code>trans</code> is used to construct the
  !!     macrostate transition probability matrix (MTPM). This is then used to estimate the probability of the
  !!     system being in each macrostate, i.e. the order parameter histogram, using the 'shooting method':
  !!     Starting with the macrostate corresponding to the lowest order parameters, i.e., macrostate 1, the 
  !!     MTPM is used to infer the relative probability of the system being in macrostate 2 - which corresponds
  !!     to the macrostate 'above' 1 with regards to its supported range of order parameters. This estimate of 
  !!     the probability of the system being in macrostate 2 is in turn used to estimate, via the MTPM, the
  !!     probability of being in macrostate 3, etc. See <code>update_eta_shooting()</code> for further
  !!     details. This method, while very fast, will never converge upon the
  !!     'ideal' weight function; however, it will yield an accurate enough weight function for most practical
  !!     purposes. Note that this method is parallelisable: the <code>trans</code> matrix from different simulations
  !!     can be pooled together to yield a more accurate one. In our implementation, the weight function is
  !!     updated every <code>update_eta_sweeps</code> sweeps.
  !!     </li>
  !!  </ul>
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>update_trans</code> </td>
  !!  <td> <code>logical</td>
  !!  <td> 
  !!  This should be <code>.true.</code> if <code>trans</code> is to be updated each move and 
  !!  <code>.false.</code> if not.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>enable_barriers</code> </td>
  !!  <td> <code>logical</code> </td>
  !!  <td> 
  !!  Flag which determines whether order parameter barriers are to be used or not.
  !!  </td>
  !! </tr>
  !! </tr>
  !! <tr>
  !!  <td> <code>barrier_dynamics</code> </td>
  !!  <td> <code>character(len=30)</code> </td>
  !!  <td> 
  !!  Flag which determines how the order parameter barriers are to change throughout the simulation. 
  !!  <code>"random"</code> evolves the macrostate the system is locked into via a random walk: the
  !!  next macrostate is decided with equal probability to be that above or that below the current
  !!  macrostate. <code>"pong_up"</code> moves to increasingly higher macrostates until the upper limit of the
  !!  supported order parameter range is encountered, at which point it reverses direction and proceeds to 
  !!  increasingly lower macrostates until it reaches the lower limit of the order parameter range, at which
  !!  point it reverses direction... <code>"pong_down"</code> instead moves initially to increasingly lower
  !!  macrostates. 
  !!  (Used if <code>enable_barriers=.true.</code>).
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>lock_moves</code> </td>
  !!  <td> <code>integer(ik)</code> </td>
  !!  <td> 
  !!  The number of moves to lock the system into one macrostate for. This should be greater than 0.
  !!  (Used if <code>enable_barriers=.true.</code>).
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>barrier_macro_low</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  An integer defining the lower position of the order parameter barrier: only microstates
  !!  belonging to macrostates >=<code>barrier_low</code> are accepted. (Used if <code>enable_barriers=.true.</code>).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>barrier_macro_high</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  An integer defining the higher position of the order parameter barrier: only microstates
  !!  belonging to macrostates <=</code>barrier_high</code> are accepted.  (Used if <code>enable_barriers=.true.</code>).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>rejected_moves_M_barrier</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Number of moves rejected because they are outwith the order parameter barriers. Note that trial moves are checked whether they are
  !!  within the range supported by <code>M_grid</code> <i>before</i> they are checked whether they are within the barriers; this counter
  !!  only accounts for moves rejected at the latter stage. (Used if <code>enable_barriers=.true.</code>).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>moves_since_lock</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  The number of moves since we were last locked into a particular macrostate.  (Used if <code>enable_barriers=.true.</code>).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>barrier_to_lock</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  This is: 1 if the low barrier is to be shifted in order to lock the system into the next
  !!  macrostate; 2 if the high barrier is to be shifted. Its value is unimportant before any barriers have
  !!  been shifted.  (Used if <code>enable_barriers=.true.</code>).
  !!  </font> </td>
  !! </tr>
  !! </table>
  real(rk), dimension(:), allocatable :: eta_grid
  real(rk), dimension(:,:), allocatable :: trans
  logical :: update_eta
  integer(ik) :: update_eta_sweeps
  logical :: update_trans
  character(len=30) update_eta_method
  character(len=*), parameter :: update_eta_method_VS="VS"
  character(len=*), parameter :: update_eta_method_shooting="shooting"
  logical :: enable_barriers
  character(len=30) :: barrier_dynamics
  character(len=*), parameter :: barrier_dynamics_random="random"
  character(len=*), parameter :: barrier_dynamics_pong_up="pong_up"
  character(len=*), parameter :: barrier_dynamics_pong_down="pong_down"
  integer(ik) :: lock_moves
  integer(ik) :: barrier_macro_low
  integer(ik) :: barrier_macro_high
  integer(ik) :: rejected_moves_M_barrier
  integer(ik) :: moves_since_lock
  integer(ik) :: barrier_to_lock


  !! <h3> Variables pertaining to equilibrium properties and their calculation </h3>
  !! <table border="1">
  !! <tr>
  !!  <td> <b> Variable </b> </td>
  !!  <td> <b> Type </b> </td>
  !!  <td> <b> Description </b> </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>calc_equil_properties</code> </td>
  !!  <td> <code>logical</code> </td>
  !!  <td>
  !!  Flag determining if equilibrium properties will be calculated using block analysis. The properties in question
  !!  are those described in the rest of this table. Where 'energy/enthalpy' is written, note that it is the energy
  !!  which is considered in NVT ensembles, and the enthalpy which is considered in NPT ensembles.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>block_sweeps</code> </td>
  !!  <td> <code>integer(ik)</code> </td>
  !!  <td>
  !!  The number of sweeps which comprise a 'block' which will be used to evaluate equilibrium properties. After
  !!  each sweep, the property under consideration is evaluated, e.g. free energy difference between phases, volume
  !!  of phase 1. There are <code>block_sweep</code> such quantities in a block, which are then used to determine the
  !!  average of the property under consideration for that block. The equilibrium value of the property is the 
  !!  'average of these averages' over many blocks, and the associated uncertainty is the standard error of the mean 
  !!  over the averages for each block. <code>block_sweeps</code> must be greater than 0.
  !!  </td>
  !! </tr>
  !! <tr> 
  !!  <td> <code>block_counts</code> </td>
  !!  <td> <code>integer(ik)</code> </td>
  !!  <td>
  !!  The number of blocks considered so far. Note that this is, say, 2, after 2 blocks have been <i>completed</i>:
  !!  it does not mean we are on our 2nd block.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>intrablock_counts_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of equilibration properties. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>intrablock_counts_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of equilibration properties. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>equil_DeltaF</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the free energy difference evaluated using block analysis. Note that this
  !!  variable is only meaningful if <code>block_counts_DeltaF>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>sigma_equil_DeltaF</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the uncertainty in the free energy difference evaluated using block analysis. Note that this
  !!  variable is only meaningful if <code>block_counts_DeltaF>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>intrablock_sum_DeltaF</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_DeltaF</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_DeltaF</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_DeltaF</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_DeltaF_sqrd</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>sigma_equil_DeltaF</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>block_counts_DeltaF</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  The number of blocks which have been considered so far for the purposes of evaluating <code>equil_DeltaF</code>
  !!  and  <code>sigma_DeltaF</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>equil_H_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the energy/enthalpy of phase 1 evaluated using block analysis. Note that this
  !!  variable is only meaningful if <code>block_counts_H_1>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>equil_H_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the energy/enthalpy of phase 2 evaluated using block analysis. Note that this
  !!  variable is only meaningful if <code>block_counts_H_2>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>sigma_equil_H_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the uncertainty in the energy/enthalpy of phase 1 evaluated using block analysis. Note that this
  !!  variable is only meaningful if <code>block_counts_H_1>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>sigma_equil_H_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the uncertainty in the energy/enthalpy of phase 2 evaluated using block analysis. Note that this
  !!  variable is only meaningful if <code>block_counts_H_2>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>intrablock_sum_H_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_H_1</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>intrablock_sum_H_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_H_2</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_H_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">
  !!  A running sum used in the evaluation of <code>equil_H_1</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_H_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_H_2</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_H_1_sqrd</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_H_1</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_H_2_sqrd</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>sigma_equil_H_2</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>block_counts_H_1</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  The number of blocks which have been considered so far for the purposes of evaluating <code>equil_H_1</code>
  !!  and  <code>sigma_H_1</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>block_counts_H_2</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  The number of blocks which have been considered so far for the purposes of evaluating <code>equil_H_2</code>
  !!  and  <code>sigma_H_2</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>equil_V_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the volume of phase 1 evaluated using block analysis. Note that this
  !!  variable is only meaningful if <code>block_counts_V_1>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>equil_V_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the volume of phase 2 evaluated using block analysis. Note that this
  !!  variable is only meaningful if <code>block_counts_V_2>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>sigma_equil_V_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the uncertainty in the volume of phase 1 evaluated using block analysis. Note that this
  !!  variable is only meaningful if <code>block_counts_V_1>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>sigma_equil_V_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the uncertainty in the volume of phase 2 evaluated using block analysis. Note that this
  !!  variable is only meaningful if <code>block_counts_V_2>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>intrablock_sum_V_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_V_1</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>intrablock_sum_V_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_V_2</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_V_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_V_1</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_V_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_V_2</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_V_1_sqrd</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>sigma_equil_V_1</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_V_2_sqrd</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>sigma_equil_V_2</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>block_counts_V_1</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  The number of blocks which have been considered so far for the purposes of evaluating <code>equil_V_1</code>
  !!  and  <code>sigma_V_1</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>block_counts_V_2</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  The number of blocks which have been considered so far for the purposes of evaluating <code>equil_V_2</code>
  !!  and  <code>sigma_V_2</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>equil_umsd_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(:), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the mean-squared displacement (MSD) of the particles from their lattice sites for phase
  !!  1 evaluated using block analysis: element <code>n</code> is the MSD for particle <code>n</code>.
  !!  Note that this variable is only meaningful if <code>block_counts_umsd_1>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>equil_umsd_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(:), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the MSD of the particles from their lattice sites for phase
  !!  2 evaluated using block analysis: element <code>n</code> is the MSD for particle <code>n</code>.
  !!  Note that this variable is only meaningful if <code>block_counts_umsd_2>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>sigma_equil_umsd_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(:), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the uncertainty in the MSD of the particles from their lattice sites for phase
  !!  1 evaluated using block analysis: element <code>n</code> pertains to particle <code>n</code>. Note that this variable 
  !!  is only meaningful if <code>block_counts_umsd_1>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>sigma_equil_umsd_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(:), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the uncertainty in the MSD of the particles from their lattice sites for phase
  !!  2 evaluated using block analysis: element <code>n</code> pertains to particle <code>n</code>. Note that this variable 
  !!  is only meaningful if <code>block_counts_umsd_2>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>intrablock_sum_umsd_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(:), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_umsd_1</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>intrablock_sum_umsd_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(:), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_umsd_2</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_umsd_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(:), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_umsd_1</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_umsd_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(:), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_umsd_2</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_umsd_1_sqrd</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(:), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_umsd_1</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_umsd_2_sqrd</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(:), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_umsd_2</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>block_counts_umsd_1</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  The number of blocks which have been considered so far for the purposes of evaluating <code>equil_umsd_1</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>block_counts_umsd_2</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  The number of blocks which have been considered so far for the purposes of evaluating <code>equil_umsd_2</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>equil_L_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(3), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the mean dimensions of the supercell for phase 1: equil_L_1(1) is the mean x-dimension of the supercell for
  !!  phase 1, equil_L_1(2) is the y-dimension, equil_L_1(3) is the z-dimension.
  !!  Note that this variable is only meaningful if <code>block_counts_L_1>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>simga_equil_L_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(3), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the uncertainties in the mean dimensions of the supercell for phase 1: sigma_equil_L_1(1) is the uncertainty
  !!  in equil_L_1(1), etc.
  !!  Note that this variable is only meaningful if <code>block_counts_L_1>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>equil_L_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(3), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the mean dimensions of the supercell for phase 2: equil_L_2(1) is the mean x-dimension of the supercell for
  !!  phase 1, equil_L_2(2) is the y-dimension, equil_L_2(3) is the z-dimension.
  !!  Note that this variable is only meaningful if <code>block_counts_L_2>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>simga_equil_L_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(3), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  Variable containing the uncertainties in the mean dimensions of the supercell for phase 2: sigma_equil_L_2(1) is the uncertainty
  !!  in equil_L_2(1), etc.
  !!  Note that this variable is only meaningful if <code>block_counts_L_2>0</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>intrablock_sum_L_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(3), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_L_1</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>intrablock_sum_L_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(3), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_L_2</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_L_1</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(3), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_L_1</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_L_2</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(3), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_L_2</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_L_1_sqrd</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(3), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_L_1</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>interblock_sum_L_2_sqrd</code> </font> </td>
  !!  <td> <font color="red">  <code>real(rk), dimension(3), allocatable</code> </font> </td>
  !!  <td> <font color="red">  
  !!  A running sum used in the evaluation of <code>equil_L_2</code>. (Details in the source code).
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>block_counts_L_1</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  The number of blocks which have been considered so far for the purposes of evaluating <code>equil_L_1</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>block_counts_L_2</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(ik)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  The number of blocks which have been considered so far for the purposes of evaluating <code>equil_L_2</code>.
  !!  </font> </td>
  !! </tr>
  !! <tr>
  !!  <td> <font color="red">  <code>sweep_equil_reference</code> </font> </td>
  !!  <td> <font color="red">  <code>integer(rk)</code> </font> </td>
  !!  <td> <font color="red">  
  !!  The sweep number which will act as a reference with regards to determining whether equilibration has occured: if
  !!  <code>sweeps>=sweep_equil_reference+equil_sweeps</code>, then equilibration is considered to occur. For 'normal'
  !!  simulations this will be 0. If the simulation has melted, and the <code>melt_option</code> variable is such that
  !!  the simulation does not 'stop', but continues, then this will be set to the sweep at which the system is 'reset',
  !!  e.g. to an ideal lattice 1, such that equilibration occurs before further variables which contribute to 
  !!  equilibrium quantities are ammended.
  !!  </font> </td>
  !! </tr>
  !! </table>
  logical :: calc_equil_properties
  integer(ik) :: block_sweeps
  integer(ik) :: block_counts
  real(rk) :: intrablock_counts_1
  real(rk) :: intrablock_counts_2
  real(rk) :: equil_DeltaF
  real(rk) :: sigma_equil_DeltaF
  real(rk) :: interblock_sum_DeltaF
  real(rk) :: interblock_sum_DeltaF_sqrd
  integer(ik) :: block_counts_DeltaF
  real(rk) :: equil_H_1
  real(rk) :: equil_H_2
  real(rk) :: sigma_equil_H_1
  real(rk) :: sigma_equil_H_2
  real(rk) :: intrablock_sum_H_1
  real(rk) :: intrablock_sum_H_2
  real(rk) :: interblock_sum_H_1
  real(rk) :: interblock_sum_H_2
  real(rk) :: interblock_sum_H_1_sqrd
  real(rk) :: interblock_sum_H_2_sqrd
  integer(ik) :: block_counts_H_1
  integer(ik) :: block_counts_H_2
  real(rk) :: equil_V_1
  real(rk) :: equil_V_2
  real(rk) :: sigma_equil_V_1
  real(rk) :: sigma_equil_V_2
  real(rk) :: intrablock_sum_V_1
  real(rk) :: intrablock_sum_V_2
  real(rk) :: interblock_sum_V_1
  real(rk) :: interblock_sum_V_2
  real(rk) :: interblock_sum_V_1_sqrd
  real(rk) :: interblock_sum_V_2_sqrd
  integer(ik) :: block_counts_V_1
  integer(ik) :: block_counts_V_2
  real(rk), dimension(:), allocatable :: equil_umsd_1
  real(rk), dimension(:), allocatable :: equil_umsd_2
  real(rk), dimension(:), allocatable :: sigma_equil_umsd_1
  real(rk), dimension(:), allocatable :: sigma_equil_umsd_2
  real(rk), dimension(:), allocatable :: intrablock_sum_umsd_1
  real(rk), dimension(:), allocatable :: intrablock_sum_umsd_2
  real(rk), dimension(:), allocatable :: interblock_sum_umsd_1
  real(rk), dimension(:), allocatable :: interblock_sum_umsd_2
  real(rk), dimension(:), allocatable :: interblock_sum_umsd_1_sqrd
  real(rk), dimension(:), allocatable :: interblock_sum_umsd_2_sqrd
  integer(ik) :: block_counts_umsd_1
  integer(ik) :: block_counts_umsd_2
  real(rk), dimension(3) :: equil_L_1
  real(rk), dimension(3) :: sigma_equil_L_1
  real(rk), dimension(3) :: equil_L_2
  real(rk), dimension(3) :: sigma_equil_L_2
  real(rk), dimension(3) :: intrablock_sum_L_1
  real(rk), dimension(3) :: intrablock_sum_L_2
  real(rk), dimension(3) :: interblock_sum_L_1
  real(rk), dimension(3) :: interblock_sum_L_2
  real(rk), dimension(3) :: interblock_sum_L_1_sqrd
  real(rk), dimension(3) :: interblock_sum_L_2_sqrd
  integer(ik) :: block_counts_L_1
  integer(ik) :: block_counts_L_2
  integer(ik) :: sweep_equil_reference


contains




  !! <h2> Procedures </h2>


  !! <h3> Key procedures: <code>initialise_from_files</code>, <code>export</code>, <code>import</code>, <code>run</code> </h3>



  !! <h4> <code> subroutine initialise_from_files(filename_params,filename_lattice) </code> </h4>
  !! <p>
  !! This procedure initialises a simulation according to simulation parameters contained in the file <code>filename_params</code> and
  !! the lattices contained in the file <code>filename_lattice</code>. Be aware that this procedure opens and
  !! closes unit 10.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> filename_params </code> </td>
  !!   <td> <code>  character(*), intent(in) </code> </td>
  !!   <td> The file name containing the simulation parameters required for initialisation. The format of this file can be
  !!   obtained by extracting lines beginning with the regular expression '\!init_params ?' (ignoring leading 
  !!   whitespace) for the source code of this module, and then removing matches to the regular expression.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> filename_lattice </code> </td>
  !!   <td> <code>  character(*), intent(in) </code> </td>
  !!   <td> The file name containing the initial lattices. The format of this file is described in the documentation
  !!   for the <code>initialise_lattices(filename)</code> subroutine.
  !!   </td>
  !!  </tr>
  !! </table>
  subroutine initialise_from_files(filename_params,filename_lattice)
    character(*), intent(in) :: filename_params
    character(*), intent(in) :: filename_lattice
    ! The lattice the system will be initialised in (cold)
    integer(ik) :: init_lattice
    ! The minimum and maximum ranges of the order parameter grid
    real(rk) :: M_grid_min
    real(rk) :: M_grid_max

    ! Be aware that the order of calling the initialisation procedures below is important.
      
    ! Import 'params' and 'lattice' variables
    call import_params()
    call initialise_lattices(filename_lattice)

    call initialise_M_variables(M_grid_size,M_grid_min,M_grid_max)
    call initialise_counters()

    ! Import the 'interactions' variables
    call initialise_interactions(Lx(1),Ly(1),Lz(1),spec_1,R_1,Lx(2),Ly(2),Lz(2),spec_2,R_2)

    call initialise_cold_microstate(init_lattice)
    if(enable_barriers) then
       call initialise_barriers()
    end if


  contains


    ! This nested subroutine imports parameters from the filename 'filename_params'
    subroutine import_params()
      character(len=20) string
      integer(ik) :: error

      ! Open the file
      open(unit=10,file=filename_params,iostat=error,status="old")
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem opening file '",filename_params,"'"
         stop 1
      end if

      ! Read the variables
      !init_params init_lattice= integer
      read(10,*,iostat=error) string, init_lattice
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'init_lattice' from file '",filename_params,"'"
         stop 1
      end if
      !init_params M_grid_size= integer
      read(10,*,iostat=error) string, M_grid_size
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'M_grid_size' from file '",filename_params,"'"
         stop 1
      end if
      !init_params M_grid_min= real
      read(10,*,iostat=error) string, M_grid_min
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'M_grid_min' from file '",filename_params,"'"
         stop 1
      end if
      !init_params M_grid_max= real
      read(10,*,iostat=error) string, M_grid_max
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'M_grid_max' from file '",filename_params,"'"
         stop 1
      end if
      !init_params enable_multicanonical= logical
      read(10,*,iostat=error) string, enable_multicanonical
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'enable_multicanonical' from file '",filename_params,"'"
         stop 1
      end if
      !init_params beta= real
      read(10,*,iostat=error) string, beta
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'beta' from file '",filename_params,"'"
         stop 1
      end if
      !init_params P= real
      read(10,*,iostat=error) string, P
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'P' from file '",filename_params,"'"
         stop 1
      end if
      !init_params enable_lattice_moves= logical
      read(10,*,iostat=error) string, enable_lattice_moves
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'enable_lattice_moves' from file '",filename_params,"'"
         stop 1
      end if
      !init_params enable_part_moves= logical
      read(10,*,iostat=error) string, enable_part_moves
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'enable_part_moves' from file '",filename_params,"'"
         stop 1
      end if
      !init_params enable_vol_moves= logical
      read(10,*,iostat=error) string, enable_vol_moves
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'enable_vol_moves' from file '",filename_params,"'"
         stop 1
      end if
      !init_params part_select= character
      read(10,*,iostat=error) string, part_select
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'part_select' from file '",filename_params,"'"
         stop 1
      end if
      !init_params part_step=  real
      read(10,*,iostat=error) string, part_step
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'part_step' from file '",filename_params,"'"
         stop 1
      end if
      !init_params enable_COM_frame= logical
      read(10,*,iostat=error) string, enable_COM_frame
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'enable_COM_frame' from file '",filename_params,"'"
         stop 1
      end if
      !init_params vol_dynamics= character
      read(10,*,iostat=error) string, vol_dynamics
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'vol_dynamics' from file '",filename_params,"'"
         stop 1
      end if
      !init_params vol_freq= integer
      read(10,*,iostat=error) string, vol_freq
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'vol_freq' from file '",filename_params,"'"
         stop 1
      end if
      !init_params vol_step= real
      read(10,*,iostat=error) string, vol_step
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'vol_step' from file '",filename_params,"'"
         stop 1
      end if
      !init_params stop_sweeps= integer
      read(10,*,iostat=error) string, stop_sweeps
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'stop_sweeps' from file '",filename_params,"'"
         stop 1
      end if
      !init_params equil_sweeps= integer
      read(10,*,iostat=error) string, equil_sweeps
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'equil_sweeps' from file '",filename_params,"'"
         stop 1
      end if
      !init_params enable_melt_checks= logical
      read(10,*,iostat=error) string, enable_melt_checks
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'enable_melt_checks' from file '",filename_params,"'"
         stop 1
      end if
      !init_params melt_sweeps= integer
      read(10,*,iostat=error) string, melt_sweeps
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'melt_sweeps' from file '",filename_params,"'"
         stop 1
      end if
      !init_params melt_threshold= real
      read(10,*,iostat=error) string, melt_threshold
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'melt_threshold' from file '",filename_params,"'"
         stop 1
      end if
      !init_params melt_option= character
      read(10,*,iostat=error) string, melt_option
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'melt_option' from file '",filename_params,"'"
         stop 1
      end if
      !init_params enable_divergence_checks= logical
      read(10,*,iostat=error) string, enable_divergence_checks
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'enable_divergence_checks' from file '",filename_params,"'"
         stop 1
      end if
      !init_params divergence_sweeps= integer
      read(10,*,iostat=error) string, divergence_sweeps
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'divergence_sweeps' from file '",filename_params,"'"
         stop 1
      end if
      !init_params divergence_tol= real
      read(10,*,iostat=error) string, divergence_tol
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'divergence_tol' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_period= integer
      read(10,*,iostat=error) string, output_file_period
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_period' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_Lx= logical
      read(10,*,iostat=error) string, output_file_Lx
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_Lx' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_Ly= logical
      read(10,*,iostat=error) string, output_file_Ly
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_Ly' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_Lz= logical
      read(10,*,iostat=error) string, output_file_Lz
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_Lz' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_V= logical
      read(10,*,iostat=error) string, output_file_V
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_V' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_R_1= logical
      read(10,*,iostat=error) string, output_file_R_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_R_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_R_2= logical
      read(10,*,iostat=error) string, output_file_R_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_R_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_u= logical
      read(10,*,iostat=error) string, output_file_u
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_u' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_lattice= logical
      read(10,*,iostat=error) string, output_file_lattice
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_lattice' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_E= logical
      read(10,*,iostat=error) string, output_file_E
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_E' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_M= logical
      read(10,*,iostat=error) string, output_file_M
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_M' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_macro= logical
      read(10,*,iostat=error) string, output_file_macro
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_macro' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_eta= logical
      read(10,*,iostat=error) string, output_file_eta
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_eta' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_moves_lattice= logical
      read(10,*,iostat=error) string, output_file_moves_lattice
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_moves_lattice' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_accepted_moves_lattice= logical
      read(10,*,iostat=error) string, output_file_accepted_moves_lattice
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_accepted_moves_lattice' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_moves_part= logical
      read(10,*,iostat=error) string, output_file_moves_part
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_moves_part' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_accepted_moves_part= logical
      read(10,*,iostat=error) string, output_file_accepted_moves_part
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_accepted_moves_part' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_moves_vol= logical
      read(10,*,iostat=error) string, output_file_moves_vol
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_moves_vol' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_accepted_moves_vol= logical
      read(10,*,iostat=error) string, output_file_accepted_moves_vol
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_accepted_moves_vol' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_rejected_moves_M_OOB= logical
      read(10,*,iostat=error) string, output_file_rejected_moves_M_OOB
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_rejected_moves_M_OOB' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_M_OOB_high= logical
      read(10,*,iostat=error) string, output_file_M_OOB_high
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_M_OOB_high' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_M_OOB_low= logical
      read(10,*,iostat=error) string, output_file_M_OOB_low
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_M_OOB_low' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_barrier_macro_low= logical
      read(10,*,iostat=error) string, output_file_barrier_macro_low
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_barrier_macro_low' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_barrier_macro_high= logical
      read(10,*,iostat=error) string, output_file_barrier_macro_high
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_barrier_macro_high' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_rejected_moves_M_barrier= logical
      read(10,*,iostat=error) string, output_file_rejected_moves_M_barrier
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_rejected_moves_M_barrier' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_moves_since_lock= logical
      read(10,*,iostat=error) string, output_file_moves_since_lock
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_moves_since_lock' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_melts= logical
      read(10,*,iostat=error) string, output_file_melts
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_melts' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_equil_DeltaF= logical
      read(10,*,iostat=error) string, output_file_equil_DeltaF
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_DeltaF' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_sigma_equil_DeltaF= logical
      read(10,*,iostat=error) string, output_file_sigma_equil_DeltaF
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_DeltaF' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_equil_H_1= logical
      read(10,*,iostat=error) string, output_file_equil_H_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_H_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_sigma_equil_H_1= logical
      read(10,*,iostat=error) string, output_file_sigma_equil_H_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_H_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_equil_H_2= logical
      read(10,*,iostat=error) string, output_file_equil_H_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_H_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_sigma_equil_H_2= logical
      read(10,*,iostat=error) string, output_file_sigma_equil_H_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_H_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_equil_V_1= logical
      read(10,*,iostat=error) string, output_file_equil_V_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_V_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_sigma_equil_V_1= logical
      read(10,*,iostat=error) string, output_file_sigma_equil_V_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_V_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_equil_V_2= logical
      read(10,*,iostat=error) string, output_file_equil_V_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_V_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_sigma_equil_V_2= logical
      read(10,*,iostat=error) string, output_file_sigma_equil_V_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_V_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_equil_umsd_1= logical
      read(10,*,iostat=error) string, output_file_equil_umsd_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_umsd_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_sigma_equil_umsd_1= logical
      read(10,*,iostat=error) string, output_file_sigma_equil_umsd_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_umsd_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_equil_umsd_2= logical
      read(10,*,iostat=error) string, output_file_equil_umsd_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_umsd_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_sigma_equil_umsd_2= logical
      read(10,*,iostat=error) string, output_file_sigma_equil_umsd_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_umsd_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_equil_L_1= logical
      read(10,*,iostat=error) string, output_file_equil_L_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_L_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_sigma_equil_L_1= logical
      read(10,*,iostat=error) string, output_file_sigma_equil_L_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_L_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_equil_L_2= logical
      read(10,*,iostat=error) string, output_file_equil_L_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_L_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_file_sigma_equil_L_2= logical
      read(10,*,iostat=error) string, output_file_sigma_equil_L_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_L_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_period= integer
      read(10,*,iostat=error) string, output_stdout_period
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_period' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_Lx= logical
      read(10,*,iostat=error) string, output_stdout_Lx
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_Lx' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_Ly= logical
      read(10,*,iostat=error) string, output_stdout_Ly
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_Ly' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_Lz= logical
      read(10,*,iostat=error) string, output_stdout_Lz
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_Lz' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_V= logical
      read(10,*,iostat=error) string, output_stdout_V
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_V' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_R_1= logical
      read(10,*,iostat=error) string, output_stdout_R_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_R_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_R_2= logical
      read(10,*,iostat=error) string, output_stdout_R_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_R_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_u= logical
      read(10,*,iostat=error) string, output_stdout_u
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_u' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_lattice= logical
      read(10,*,iostat=error) string, output_stdout_lattice
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_lattice' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_E= logical
      read(10,*,iostat=error) string, output_stdout_E
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_E' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_M= logical
      read(10,*,iostat=error) string, output_stdout_M
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_M' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_macro= logical
      read(10,*,iostat=error) string, output_stdout_macro
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_macro' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_eta= logical
      read(10,*,iostat=error) string, output_stdout_eta
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_eta' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_moves_lattice= logical
      read(10,*,iostat=error) string, output_stdout_moves_lattice
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_moves_lattice' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_accepted_moves_lattice= logical
      read(10,*,iostat=error) string, output_stdout_accepted_moves_lattice
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_accepted_moves_lattice' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_moves_part= logical
      read(10,*,iostat=error) string, output_stdout_moves_part
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_moves_part' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_accepted_moves_part= logical
      read(10,*,iostat=error) string, output_stdout_accepted_moves_part
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_accepted_moves_part' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_moves_vol= logical
      read(10,*,iostat=error) string, output_stdout_moves_vol
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_moves_vol' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_accepted_moves_vol= logical
      read(10,*,iostat=error) string, output_stdout_accepted_moves_vol
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_accepted_moves_vol' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_rejected_moves_M_OOB= logical
      read(10,*,iostat=error) string, output_stdout_rejected_moves_M_OOB
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_rejected_moves_M_OOB' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_M_OOB_high= logical
      read(10,*,iostat=error) string, output_stdout_M_OOB_high
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_M_OOB_high' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_M_OOB_low= logical
      read(10,*,iostat=error) string, output_stdout_M_OOB_low
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_M_OOB_low' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_barrier_macro_low= logical
      read(10,*,iostat=error) string, output_stdout_barrier_macro_low
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_barrier_macro_low' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_barrier_macro_high= logical
      read(10,*,iostat=error) string, output_stdout_barrier_macro_high
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_barrier_macro_high' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_rejected_moves_M_barrier= logical
      read(10,*,iostat=error) string, output_stdout_rejected_moves_M_barrier
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_rejected_moves_M_barrier' from file '", &
              filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_moves_since_lock= logical
      read(10,*,iostat=error) string, output_stdout_moves_since_lock
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_moves_since_lock' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_melts= logical
      read(10,*,iostat=error) string, output_stdout_melts
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_melts' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_equil_DeltaF= logical
      read(10,*,iostat=error) string, output_stdout_equil_DeltaF
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_DeltaF' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_sigma_equil_DeltaF= logical
      read(10,*,iostat=error) string, output_stdout_sigma_equil_DeltaF
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_DeltaF' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_equil_H_1= logical
      read(10,*,iostat=error) string, output_stdout_equil_H_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_H_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_sigma_equil_H_1= logical
      read(10,*,iostat=error) string, output_stdout_sigma_equil_H_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_H_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_equil_H_2= logical
      read(10,*,iostat=error) string, output_stdout_equil_H_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_H_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_sigma_equil_H_2= logical
      read(10,*,iostat=error) string, output_stdout_sigma_equil_H_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_H_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_equil_V_1= logical
      read(10,*,iostat=error) string, output_stdout_equil_V_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_V_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_sigma_equil_V_1= logical
      read(10,*,iostat=error) string, output_stdout_sigma_equil_V_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_V_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_equil_V_2= logical
      read(10,*,iostat=error) string, output_stdout_equil_V_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_V_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_sigma_equil_V_2= logical
      read(10,*,iostat=error) string, output_stdout_sigma_equil_V_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_V_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_equil_umsd_1= logical
      read(10,*,iostat=error) string, output_stdout_equil_umsd_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_umsd_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_sigma_equil_umsd_1= logical
      read(10,*,iostat=error) string, output_stdout_sigma_equil_umsd_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_umsd_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_equil_umsd_2= logical
      read(10,*,iostat=error) string, output_stdout_equil_umsd_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_umsd_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_sigma_equil_umsd_2= logical
      read(10,*,iostat=error) string, output_stdout_sigma_equil_umsd_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_umsd_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_equil_L_1= logical
      read(10,*,iostat=error) string, output_stdout_equil_L_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_L_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_sigma_equil_L_1= logical
      read(10,*,iostat=error) string, output_stdout_sigma_equil_L_1
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_L_1' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_equil_L_2= logical
      read(10,*,iostat=error) string, output_stdout_equil_L_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_L_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params output_stdout_sigma_equil_L_2= logical
      read(10,*,iostat=error) string, output_stdout_sigma_equil_L_2
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_L_2' from file '",filename_params,"'"
         stop 1
      end if
      !init_params checkpoint_period= integer
      read(10,*,iostat=error) string, checkpoint_period
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'checkpoint_period' from file '",filename_params,"'"
         stop 1
      end if
      !init_params update_eta= logical
      read(10,*,iostat=error) string, update_eta
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'update_eta' from file '",filename_params,"'"
         stop 1
      end if
      !init_params update_eta_sweeps= integer
      read(10,*,iostat=error) string, update_eta_sweeps
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'update_eta_sweeps' from file '",filename_params,"'"
         stop 1
      end if
      !init_params update_trans= logical
      read(10,*,iostat=error) string, update_trans
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'update_trans' from file '",filename_params,"'"
         stop 1
      end if
      !init_params update_eta_method= character
      read(10,*,iostat=error) string, update_eta_method
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'update_eta_method' from file '",filename_params,"'"
         stop 1
      end if
      !init_params enable_barriers= logical
      read(10,*,iostat=error) string, enable_barriers
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'enable_barriers' from file '",filename_params,"'"
         stop 1
      end if
      !init_params barrier_dynamics= character
      read(10,*,iostat=error) string, barrier_dynamics
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'barrier_dynamics' from file '",filename_params,"'"
         stop 1
      end if
      !init_params lock_moves= integer
      read(10,*,iostat=error) string, lock_moves
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'lock_moves' from file '",filename_params,"'"
         stop 1
      end if
      !init_params calc_equil_properties= logical
      read(10,*,iostat=error) string, calc_equil_properties
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'calc_equil_properties' from file '",filename_params,"'"
         stop 1
      end if
      !init_params block_sweeps= integer
      read(10,*,iostat=error) string, block_sweeps
      if(error/=0) then
         write(0,*) "monteswitch_mod: Error. Problem reading 'block_sweeps' from file '",filename_params,"'"
         stop 1
      end if

      close(unit=10)

    end subroutine import_params
    
  end subroutine initialise_from_files




  !! <h4> <code> subroutine export(filename) </code> </h4>
  !! <p>
  !! This subroutine exports the state of the simulation to the specified file. Be aware that
  !! this subroutine opens and closes unit 10. Be aware also that if <code>checkpoint_period<=0</code>
  !! then this subroutine creates only an empty file.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> filename </code> </td>
  !!   <td> <code> character(*), intent(in) </code> </td>
  !!   <td> The file to which the state of the simulation will be exported. </td>
  !!  </tr>
  !! </table>
  subroutine export(filename)
    character(*), intent(in) :: filename

    open(unit=10,file=filename,status="replace")

    ! Output something only if checkpoint_period>0; otherwise we have an empty file
    if(checkpoint_period>0) then
       ! Output scalars first: these will be imported first and then used to allocate arrays to the required size
       ! Also, it is easier to read the file if the arrays are at the end
       write(10,*) "enable_multicanonical= ",enable_multicanonical
       write(10,*) "beta= ",beta
       write(10,*) "P= ",P
       write(10,*) "enable_lattice_moves= ",enable_lattice_moves
       write(10,*) "enable_part_moves= ",enable_part_moves
       write(10,*) "enable_vol_moves= ",enable_vol_moves
       write(10,*) "part_select= ",part_select
       write(10,*) "part_step= ",part_step
       write(10,*) "enable_COM_frame= ",enable_COM_frame
       write(10,*) "vol_dynamics= ",vol_dynamics
       write(10,*) "vol_freq= ",vol_freq
       write(10,*) "vol_step= ",vol_step
       write(10,*) "stop_sweeps= ",stop_sweeps
       write(10,*) "equil_sweeps= ",equil_sweeps
       write(10,*) "enable_melt_checks= ",enable_melt_checks
       write(10,*) "melt_sweeps= ",melt_sweeps
       write(10,*) "melt_threshold= ",melt_threshold
       write(10,*) "melt_option= ",melt_option
       write(10,*) "enable_divergence_checks= ",enable_divergence_checks
       write(10,*) "divergence_sweeps= ",divergence_sweeps
       write(10,*) "divergence_tol= ",divergence_tol
       write(10,*) "output_file_period= ",output_file_period
       write(10,*) "output_file_Lx= ",output_file_Lx
       write(10,*) "output_file_Ly= ",output_file_Ly
       write(10,*) "output_file_Lz= ",output_file_Lz
       write(10,*) "output_file_V= ",output_file_V
       write(10,*) "output_file_R_1= ",output_file_R_1
       write(10,*) "output_file_R_2= ",output_file_R_2
       write(10,*) "output_file_u= ",output_file_u
       write(10,*) "output_file_lattice= ",output_file_lattice
       write(10,*) "output_file_E= ",output_file_E
       write(10,*) "output_file_M= ",output_file_M
       write(10,*) "output_file_macro= ",output_file_macro
       write(10,*) "output_file_eta= ",output_file_eta
       write(10,*) "output_file_moves_lattice= ",output_file_moves_lattice
       write(10,*) "output_file_accepted_moves_lattice= ",output_file_accepted_moves_lattice
       write(10,*) "output_file_moves_part= ",output_file_moves_part
       write(10,*) "output_file_accepted_moves_part= ",output_file_accepted_moves_part
       write(10,*) "output_file_moves_vol= ",output_file_moves_vol
       write(10,*) "output_file_accepted_moves_vol= ",output_file_accepted_moves_vol
       write(10,*) "output_file_rejected_moves_M_OOB= ",output_file_rejected_moves_M_OOB
       write(10,*) "output_file_M_OOB_high= ",output_file_M_OOB_high
       write(10,*) "output_file_M_OOB_low= ",output_file_M_OOB_low
       write(10,*) "output_file_barrier_macro_low= ",output_file_barrier_macro_low
       write(10,*) "output_file_barrier_macro_high= ",output_file_barrier_macro_high
       write(10,*) "output_file_rejected_moves_M_barrier= ",output_file_rejected_moves_M_barrier
       write(10,*) "output_file_moves_since_lock= ",output_file_moves_since_lock
       write(10,*) "output_file_melts= ",output_file_melts
       write(10,*) "output_file_equil_DeltaF= ",output_file_equil_DeltaF
       write(10,*) "output_file_sigma_equil_DeltaF= ",output_file_sigma_equil_DeltaF
       write(10,*) "output_file_equil_H_1= ",output_file_equil_H_1
       write(10,*) "output_file_sigma_equil_H_1= ",output_file_sigma_equil_H_1
       write(10,*) "output_file_equil_H_2= ",output_file_equil_H_2
       write(10,*) "output_file_sigma_equil_H_2= ",output_file_sigma_equil_H_2
       write(10,*) "output_file_equil_V_1= ",output_file_equil_V_1
       write(10,*) "output_file_sigma_equil_V_1= ",output_file_sigma_equil_V_1
       write(10,*) "output_file_equil_V_2= ",output_file_equil_V_2
       write(10,*) "output_file_sigma_equil_V_2= ",output_file_sigma_equil_V_2
       write(10,*) "output_file_equil_umsd_1= ",output_file_equil_umsd_1
       write(10,*) "output_file_sigma_equil_umsd_1= ",output_file_sigma_equil_umsd_1
       write(10,*) "output_file_equil_umsd_2= ",output_file_equil_umsd_2
       write(10,*) "output_file_sigma_equil_umsd_2= ",output_file_sigma_equil_umsd_2
       write(10,*) "output_file_equil_L_1= ",output_file_equil_L_1
       write(10,*) "output_file_sigma_equil_L_1= ",output_file_sigma_equil_L_1
       write(10,*) "output_file_equil_L_2= ",output_file_equil_L_2
       write(10,*) "output_file_sigma_equil_L_2= ",output_file_sigma_equil_L_2
       write(10,*) "output_stdout_period= ",output_stdout_period
       write(10,*) "output_stdout_Lx= ",output_stdout_Lx
       write(10,*) "output_stdout_Ly= ",output_stdout_Ly
       write(10,*) "output_stdout_Lz= ",output_stdout_Lz
       write(10,*) "output_stdout_V= ",output_stdout_V
       write(10,*) "output_stdout_R_1= ",output_stdout_R_1
       write(10,*) "output_stdout_R_2= ",output_stdout_R_2
       write(10,*) "output_stdout_u= ",output_stdout_u
       write(10,*) "output_stdout_lattice= ",output_stdout_lattice
       write(10,*) "output_stdout_E= ",output_stdout_E
       write(10,*) "output_stdout_M= ",output_stdout_M
       write(10,*) "output_stdout_macro= ",output_stdout_macro
       write(10,*) "output_stdout_eta= ",output_stdout_eta
       write(10,*) "output_stdout_moves_lattice= ",output_stdout_moves_lattice
       write(10,*) "output_stdout_accepted_moves_lattice= ",output_stdout_accepted_moves_lattice
       write(10,*) "output_stdout_moves_part= ",output_stdout_moves_part
       write(10,*) "output_stdout_accepted_moves_part= ",output_stdout_accepted_moves_part
       write(10,*) "output_stdout_moves_vol= ",output_stdout_moves_vol
       write(10,*) "output_stdout_accepted_moves_vol= ",output_stdout_accepted_moves_vol
       write(10,*) "output_stdout_rejected_moves_M_OOB= ",output_stdout_rejected_moves_M_OOB
       write(10,*) "output_stdout_M_OOB_high= ",output_stdout_M_OOB_high
       write(10,*) "output_stdout_M_OOB_low= ",output_stdout_M_OOB_low
       write(10,*) "output_stdout_barrier_macro_low= ",output_stdout_barrier_macro_low
       write(10,*) "output_stdout_barrier_macro_high= ",output_stdout_barrier_macro_high
       write(10,*) "output_stdout_rejected_moves_M_barrier= ",output_stdout_rejected_moves_M_barrier
       write(10,*) "output_stdout_moves_since_lock= ",output_stdout_moves_since_lock
       write(10,*) "output_stdout_melts= ",output_stdout_melts
       write(10,*) "output_stdout_equil_DeltaF= ",output_stdout_equil_DeltaF
       write(10,*) "output_stdout_sigma_equil_DeltaF= ",output_stdout_sigma_equil_DeltaF
       write(10,*) "output_stdout_equil_H_1= ",output_stdout_equil_H_1
       write(10,*) "output_stdout_sigma_equil_H_1= ",output_stdout_sigma_equil_H_1
       write(10,*) "output_stdout_equil_H_2= ",output_stdout_equil_H_2
       write(10,*) "output_stdout_sigma_equil_H_2= ",output_stdout_sigma_equil_H_2
       write(10,*) "output_stdout_equil_V_1= ",output_stdout_equil_V_1
       write(10,*) "output_stdout_sigma_equil_V_1= ",output_stdout_sigma_equil_V_1
       write(10,*) "output_stdout_equil_V_2= ",output_stdout_equil_V_2
       write(10,*) "output_stdout_sigma_equil_V_2= ",output_stdout_sigma_equil_V_2
       write(10,*) "output_stdout_equil_umsd_1= ",output_stdout_equil_umsd_1
       write(10,*) "output_stdout_sigma_equil_umsd_1= ",output_stdout_sigma_equil_umsd_1
       write(10,*) "output_stdout_equil_umsd_2= ",output_stdout_equil_umsd_2
       write(10,*) "output_stdout_sigma_equil_umsd_2= ",output_stdout_sigma_equil_umsd_2
       write(10,*) "output_stdout_equil_L_1= ",output_stdout_equil_L_1
       write(10,*) "output_stdout_sigma_equil_L_1= ",output_stdout_sigma_equil_L_1
       write(10,*) "output_stdout_equil_L_2= ",output_stdout_equil_L_2
       write(10,*) "output_stdout_sigma_equil_L_2= ",output_stdout_sigma_equil_L_2
       write(10,*) "checkpoint_period= ",checkpoint_period
       write(10,*) "M_grid_size= ",M_grid_size
       write(10,*) "n_part= ",n_part
       write(10,*) "Lx= ",Lx
       write(10,*) "Ly= ",Ly
       write(10,*) "Lz= ",Lz
       write(10,*) "V= ",V
       write(10,*) "lattice= ",lattice
       write(10,*) "E_1= ",E_1
       write(10,*) "E_2= ",E_2
       write(10,*) "E= ",E
       write(10,*) "M= ",M
       write(10,*) "macro= ",macro
       write(10,*) "eta= ",eta
       write(10,*) "M_const= ",M_const
       write(10,*) "switchscalex= ",switchscalex
       write(10,*) "switchscaley= ",switchscaley
       write(10,*) "switchscalez= ",switchscalez
       write(10,*) "sweeps= ",sweeps
       write(10,*) "moves= ",moves
       write(10,*) "moves_lattice= ",moves_lattice
       write(10,*) "accepted_moves_lattice= ",accepted_moves_lattice
       write(10,*) "moves_part= ",moves_part
       write(10,*) "accepted_moves_part= ",accepted_moves_part
       write(10,*) "moves_vol= ",moves_vol
       write(10,*) "accepted_moves_vol= ",accepted_moves_vol
       write(10,*) "rejected_moves_M_OOB= ",rejected_moves_M_OOB
       write(10,*) "M_OOB_high= ",M_OOB_high
       write(10,*) "M_OOB_low= ",M_OOB_low
       write(10,*) "melts= ",melts
       write(10,*) "update_eta= ",update_eta
       write(10,*) "update_eta_sweeps= ",update_eta_sweeps
       write(10,*) "update_eta_method= ",update_eta_method
       write(10,*) "update_trans= ",update_trans
       write(10,*) "enable_barriers= ",enable_barriers
       write(10,*) "barrier_dynamics= ",barrier_dynamics
       write(10,*) "lock_moves= ",lock_moves
       write(10,*) "barrier_macro_low= ",barrier_macro_low
       write(10,*) "barrier_macro_high= ",barrier_macro_high
       write(10,*) "rejected_moves_M_barrier= ",rejected_moves_M_barrier
       write(10,*) "moves_since_lock= ",moves_since_lock
       write(10,*) "barrier_to_lock= ",barrier_to_lock
       write(10,*) "calc_equil_properties= ",calc_equil_properties
       write(10,*) "block_sweeps= ",block_sweeps
       write(10,*) "block_counts= ",block_counts
       write(10,*) "intrablock_counts_1= ",intrablock_counts_1
       write(10,*) "intrablock_counts_2= ",intrablock_counts_2
       write(10,*) "equil_DeltaF= ",equil_DeltaF
       write(10,*) "sigma_equil_DeltaF= ",sigma_equil_DeltaF
       write(10,*) "interblock_sum_DeltaF= ",interblock_sum_DeltaF
       write(10,*) "interblock_sum_DeltaF_sqrd= ",interblock_sum_DeltaF_sqrd
       write(10,*) "block_counts_DeltaF= ",block_counts_DeltaF
       write(10,*) "equil_H_1= ",equil_H_1
       write(10,*) "equil_H_2= ",equil_H_2
       write(10,*) "sigma_equil_H_1= ",sigma_equil_H_1
       write(10,*) "sigma_equil_H_2= ",sigma_equil_H_2
       write(10,*) "intrablock_sum_H_1= ",intrablock_sum_H_1
       write(10,*) "intrablock_sum_H_2= ",intrablock_sum_H_2
       write(10,*) "interblock_sum_H_1= ",interblock_sum_H_1
       write(10,*) "interblock_sum_H_2= ",interblock_sum_H_2
       write(10,*) "interblock_sum_H_1_sqrd= ",interblock_sum_H_1_sqrd
       write(10,*) "interblock_sum_H_2_sqrd= ",interblock_sum_H_2_sqrd
       write(10,*) "block_counts_H_1= ",block_counts_H_1
       write(10,*) "block_counts_H_2= ",block_counts_H_2
       write(10,*) "equil_V_1= ",equil_V_1
       write(10,*) "equil_V_2= ",equil_V_2
       write(10,*) "sigma_equil_V_1= ",sigma_equil_V_1
       write(10,*) "sigma_equil_V_2= ",sigma_equil_V_2
       write(10,*) "intrablock_sum_V_1= ",intrablock_sum_V_1
       write(10,*) "intrablock_sum_V_2= ",intrablock_sum_V_2
       write(10,*) "interblock_sum_V_1= ",interblock_sum_V_1
       write(10,*) "interblock_sum_V_2= ",interblock_sum_V_2
       write(10,*) "interblock_sum_V_1_sqrd= ",interblock_sum_V_1_sqrd
       write(10,*) "interblock_sum_V_2_sqrd= ",interblock_sum_V_2_sqrd
       write(10,*) "block_counts_V_1= ",block_counts_V_1
       write(10,*) "block_counts_V_2= ",block_counts_V_2
       write(10,*) "block_counts_umsd_1= ",block_counts_umsd_1
       write(10,*) "block_counts_umsd_2= ",block_counts_umsd_2
       write(10,*) "equil_L_1= ",equil_L_1
       write(10,*) "equil_L_2= ",equil_L_2
       write(10,*) "sigma_equil_L_1= ",sigma_equil_L_1
       write(10,*) "sigma_equil_L_2= ",sigma_equil_L_2
       write(10,*) "intrablock_sum_L_1= ",intrablock_sum_L_1
       write(10,*) "intrablock_sum_L_2= ",intrablock_sum_L_2
       write(10,*) "interblock_sum_L_1= ",interblock_sum_L_1
       write(10,*) "interblock_sum_L_2= ",interblock_sum_L_2
       write(10,*) "interblock_sum_L_1_sqrd= ",interblock_sum_L_1_sqrd
       write(10,*) "interblock_sum_L_2_sqrd= ",interblock_sum_L_2_sqrd
       write(10,*) "block_counts_L_1= ",block_counts_L_1
       write(10,*) "block_counts_L_2= ",block_counts_L_2
       write(10,*) "sweep_equil_reference= ",sweep_equil_reference
       ! Output allocatable arrays
       write(10,*) "spec_1= ",spec_1
       write(10,*) "spec_2= ",spec_2
       write(10,*) "pos_1= ",pos_1
       write(10,*) "pos_2= ",pos_2
       write(10,*) "R_1= ",R_1
       write(10,*) "R_2= ",R_2
       write(10,*) "u= ",u
       write(10,*) "M_grid= ",M_grid
       write(10,*) "M_counts_1= ",M_counts_1
       write(10,*) "M_counts_2= ",M_counts_2
       write(10,*) "eta_grid= ",eta_grid
       write(10,*) "trans= ",trans
       write(10,*) "equil_umsd_1= ",equil_umsd_1
       write(10,*) "equil_umsd_2= ",equil_umsd_2
       write(10,*) "sigma_equil_umsd_1= ",sigma_equil_umsd_1
       write(10,*) "sigma_equil_umsd_2= ",sigma_equil_umsd_2
       write(10,*) "intrablock_sum_umsd_1= ",intrablock_sum_umsd_1
       write(10,*) "intrablock_sum_umsd_2= ",intrablock_sum_umsd_2
       write(10,*) "interblock_sum_umsd_1= ",interblock_sum_umsd_1
       write(10,*) "interblock_sum_umsd_2= ",interblock_sum_umsd_2
       write(10,*) "interblock_sum_umsd_1_sqrd= ",interblock_sum_umsd_1_sqrd
       write(10,*) "interblock_sum_umsd_2_sqrd= ",interblock_sum_umsd_2_sqrd

       ! Export 'interactions' variables
       call export_interactions()

    end if
    close(unit=10)
  end subroutine export




  !! <h4> <code> subroutine import(filename) </code> </h4>
  !! <p>
  !! This subroutine imports the state of the simulation from the specified file. Be
  !! aware that this subroutine opens and closes unit 10.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> filename </code> </td>
  !!   <td> <code> character(*), intent(in) </code> </td>
  !!   <td> The file from which the state of the simulation will be imported. </td>
  !!  </tr>
  !! </table>
  subroutine import(filename)
    character(*), intent(in) :: filename
    integer(ik) :: error
    character(20) :: string

    open(unit=10,file=filename,iostat=error,status="old")
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem opening file '",trim(filename),"'"
       stop 1
    end if
    
    ! Read the scalar variables. Then use these to allocate arrays before reading the array variables.
    read(10,*,iostat=error) string, enable_multicanonical
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'enable_multicanonical' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, beta
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'beta' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, P
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'P' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, enable_lattice_moves
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'enable_lattice_moves' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, enable_part_moves
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'enable_part_moves' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, enable_vol_moves
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'enable_vol_moves' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, part_select
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'part_select' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, part_step
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'part_step' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, enable_COM_frame
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'enable_COM_frame' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, vol_dynamics
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'vol_dynamics' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, vol_freq
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'vol_freq' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, vol_step
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'vol_step' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, stop_sweeps
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'stop_sweeps' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, equil_sweeps
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'equil_sweeps' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, enable_melt_checks
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'enable_melt_checks' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, melt_sweeps
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'melt_sweeps' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, melt_threshold
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'melt_threshold' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, melt_option
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'melt_option' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, enable_divergence_checks
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'enable_divergence_checks' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, divergence_sweeps
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'divergence_sweeps' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, divergence_tol
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'divergence_tol' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_period
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_period' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_Lx
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_Lx' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_Ly
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_Ly' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_Lz
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_Lz' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_V
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_V' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_R_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_R_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_R_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_R_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_u
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_u' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_lattice
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_lattice' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_E
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_E' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_M
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_M' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_macro
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_macro' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_eta
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_eta' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_moves_lattice
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_moves_lattice' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_accepted_moves_lattice
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_accepted_moves_lattice' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_moves_part
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_moves_part' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_accepted_moves_part
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_accepted_moves_part' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_moves_vol
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_moves_vol' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_accepted_moves_vol
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_accepted_moves_vol' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_rejected_moves_M_OOB
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_rejected_moves_M_OOB' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_M_OOB_high
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_M_OOB_high' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_M_OOB_low
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_M_OOB_low' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_barrier_macro_low
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_barrier_macro_low' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_barrier_macro_high
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_barrier_macro_high' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_rejected_moves_M_barrier
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_rejected_moves_M_barrier' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_moves_since_lock
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_moves_since_lock' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_melts
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_melts' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_equil_DeltaF
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_DeltaF' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_sigma_equil_DeltaF
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_DeltaF' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_equil_H_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_H_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_sigma_equil_H_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_H_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_equil_H_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_H_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_sigma_equil_H_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_H_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_equil_V_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_V_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_sigma_equil_V_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_V_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_equil_V_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_V_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_sigma_equil_V_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_V_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_equil_umsd_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_umsd_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_sigma_equil_umsd_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_umsd_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_equil_umsd_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_umsd_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_sigma_equil_umsd_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_umsd_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_equil_L_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_L_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_sigma_equil_L_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_L_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_equil_L_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_equil_L_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_file_sigma_equil_L_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_file_sigma_equil_L_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_period
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_period' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_Lx
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_Lx' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_Ly
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_Ly' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_Lz
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_Lz' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_V
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_V' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_R_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_R_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_R_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_R_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_u
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_u' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_lattice
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_lattice' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_E
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_E' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_M
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_M' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_macro
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_macro' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_eta
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_eta' from file '",trim(filename)
       stop 1
    end if 
    read(10,*,iostat=error) string, output_stdout_moves_lattice
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_moves_lattice' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_accepted_moves_lattice
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_accepted_moves_lattice' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_moves_part
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_moves_part' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_accepted_moves_part
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_accepted_moves_part' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_moves_vol
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_moves_vol' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_accepted_moves_vol
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_accepted_moves_vol' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_rejected_moves_M_OOB
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_rejected_moves_M_OOB' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_M_OOB_high
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_M_OOB_high' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_M_OOB_low
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_M_OOB_low' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_barrier_macro_low
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_barrier_macro_low' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_barrier_macro_high
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_barrier_macro_high' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_rejected_moves_M_barrier
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_rejected_moves_M_barrier' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_moves_since_lock
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_moves_since_lock' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_melts
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_melts' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_equil_DeltaF
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_DeltaF' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_sigma_equil_DeltaF
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_DeltaF' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_equil_H_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_H_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_sigma_equil_H_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_H_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_equil_H_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_H_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_sigma_equil_H_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_H_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_equil_V_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_V_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_sigma_equil_V_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_V_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_equil_V_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_V_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_sigma_equil_V_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_V_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_equil_umsd_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_umsd_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_sigma_equil_umsd_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_umsd_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_equil_umsd_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_umsd_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_sigma_equil_umsd_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_umsd_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_equil_L_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_L_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_sigma_equil_L_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_L_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_equil_L_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_equil_L_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, output_stdout_sigma_equil_L_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'output_stdout_sigma_equil_L_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, checkpoint_period
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'checkpoint_period' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, M_grid_size
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'M_grid_size' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, n_part
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'n_part' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, Lx
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'Lx' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, Ly
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'Ly' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, Lz
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'Lz' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, V
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'V' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, lattice
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'lattice' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, E_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'E_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, E_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'E_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, E
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'E' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, M
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'M' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, macro
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'macro' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, eta
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'eta' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, M_const
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'M_const' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, switchscalex
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'switchscalex' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, switchscaley
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'switchscaley' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, switchscalez
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'switchscalez' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, sweeps
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'sweeps' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, moves
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'moves' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, moves_lattice
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'moves_lattice' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, accepted_moves_lattice
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'accepted_moves_lattice' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, moves_part
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'moves_part' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, accepted_moves_part
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'accepted_moves_part' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, moves_vol
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'moves_vol' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, accepted_moves_vol
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'accepted_moves_vol' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, rejected_moves_M_OOB
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'rejected_moves_M_OOB' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, M_OOB_high
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'M_OOB_high' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, M_OOB_low
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'M_OOB_low' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, melts
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'melts' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, update_eta
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'update_eta' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, update_eta_sweeps
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'update_eta_sweeps' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, update_eta_method
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'update_eta_method' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, update_trans
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'update_trans' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, enable_barriers
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'enable_barriers' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, barrier_dynamics
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'barrier_dynamics' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, lock_moves
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'lock_moves' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, barrier_macro_low
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'barrier_macro_low' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, barrier_macro_high
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'barrier_macro_high' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, rejected_moves_M_barrier
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'rejected_moves_M_barrier' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, moves_since_lock
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'moves_since_lock' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, barrier_to_lock
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'barrier_to_lock' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, calc_equil_properties
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'calc_equil_properties' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, block_sweeps
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'block_sweeps' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, block_counts
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'block_counts' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, intrablock_counts_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'intrablock_counts_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, intrablock_counts_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'intrablock_counts_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, equil_DeltaF
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'equil_DeltaF' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, sigma_equil_DeltaF
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'sigma_equil_DeltaF' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, interblock_sum_DeltaF
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_DeltaF' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, interblock_sum_DeltaF_sqrd
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_DeltaF_sqrd' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, block_counts_DeltaF
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'block_counts_DeltaF' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, equil_H_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'equil_H_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, equil_H_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'equil_H_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, sigma_equil_H_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'sigma_equil_H_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, sigma_equil_H_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'sigma_equil_H_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, intrablock_sum_H_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'intrablock_sum_H_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, intrablock_sum_H_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'intrablock_sum_H_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, interblock_sum_H_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_H_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, interblock_sum_H_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_H_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, interblock_sum_H_1_sqrd
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_H_1_sqrd' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, interblock_sum_H_2_sqrd
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_H_2_sqrd' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, block_counts_H_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'block_counts_H_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, block_counts_H_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'block_counts_H_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, equil_V_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'equil_V_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, equil_V_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'equil_V_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, sigma_equil_V_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'sigma_equil_V_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, sigma_equil_V_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'sigma_equil_V_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, intrablock_sum_V_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'intrablock_sum_V_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, intrablock_sum_V_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'intrablock_sum_V_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, interblock_sum_V_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_V_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, interblock_sum_V_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_V_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, interblock_sum_V_1_sqrd
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_V_1_sqrd' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, interblock_sum_V_2_sqrd
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_V_2_sqrd' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, block_counts_V_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'block_counts_V_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, block_counts_V_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'block_counts_V_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, block_counts_umsd_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'block_counts_umsd_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, block_counts_umsd_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'block_counts_umsd_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, equil_L_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'equil_L_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, equil_L_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'equil_L_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, sigma_equil_L_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'sigma_equil_L_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, sigma_equil_L_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'sigma_equil_L_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, intrablock_sum_L_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'intrablock_sum_L_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, intrablock_sum_L_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'intrablock_sum_L_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, interblock_sum_L_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_L_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, interblock_sum_L_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_L_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, interblock_sum_L_1_sqrd
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_L_1_sqrd' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, interblock_sum_L_2_sqrd
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_L_2_sqrd' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, block_counts_L_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'block_counts_L_1' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, block_counts_L_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'block_counts_L_2' from file '",trim(filename)
       stop 1
    end if
    read(10,*,iostat=error) string, sweep_equil_reference
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'sweep_equil_reference' from file '",trim(filename)
       stop 1
    end if

    ! Read the array values from the file. Note that if the arrays are 'not in use' then
    ! we read only the first token (as a string) of each line - since there is no second, third etc.
    ! tokens.
    !
    ! spec_1
    if(allocated(spec_1)) then
       deallocate(spec_1)
    end if
    allocate(spec_1(n_part))
    read(10,*,iostat=error) string, spec_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'spec_1' from file '",trim(filename)
       stop 1
    end if
    ! spec_2
    if(allocated(spec_2)) then
       deallocate(spec_2)
    end if
    allocate(spec_2(n_part))
    read(10,*,iostat=error) string, spec_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'spec_2' from file '",trim(filename)
       stop 1
    end if
    ! pos_1
    if(allocated(pos_1)) then
       deallocate(pos_1)
    end if
    allocate(pos_1(3,n_part))
    read(10,*,iostat=error) string, pos_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'pos_1' from file '",trim(filename)
       stop 1
    end if
    ! pos_2
    if(allocated(pos_2)) then
       deallocate(pos_2)
    end if
    allocate(pos_2(3,n_part))
    read(10,*,iostat=error) string, pos_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'pos_2' from file '",trim(filename)
       stop 1
    end if
    ! R_1
    if(allocated(R_1)) then
       deallocate(R_1)
    end if
    allocate(R_1(3,n_part))
    read(10,*,iostat=error) string, R_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'R_1' from file '",trim(filename)
       stop 1
    end if
    ! R_2
    if(allocated(R_2)) then
       deallocate(R_2)
    end if
    allocate(R_2(3,n_part))
    read(10,*,iostat=error) string, R_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'R_2' from file '",trim(filename)
       stop 1
    end if
    ! u
    if(allocated(u)) then
       deallocate(u)
    end if
    allocate(u(3,n_part))
    read(10,*,iostat=error) string, u
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'u' from file '",trim(filename)
       stop 1
    end if
    ! M_grid
    if(allocated(M_grid)) then
       deallocate(M_grid)
    end if
    allocate(M_grid(M_grid_size))
    read(10,*,iostat=error) string, M_grid
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'M_grid' from file '",trim(filename)
       stop 1
    end if
    ! M_counts_1
    if(allocated(M_counts_1)) then
       deallocate(M_counts_1)
    end if
    allocate(M_counts_1(M_grid_size))
    read(10,*,iostat=error) string, M_counts_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'M_counts_1' from file '",trim(filename)
       stop 1
    end if
    ! M_counts_2
    if(allocated(M_counts_2)) then
       deallocate(M_counts_2)
    end if
    allocate(M_counts_2(M_grid_size))
    read(10,*,iostat=error) string, M_counts_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'M_counts_2' from file '",trim(filename)
       stop 1
    end if
    ! eta_grid
    if(allocated(eta_grid)) then
       deallocate(eta_grid)
    end if
    allocate(eta_grid(M_grid_size))
    read(10,*,iostat=error) string, eta_grid
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'eta_grid' from file '",trim(filename)
       stop 1
    end if
    ! trans
    if(allocated(trans)) then
       deallocate(trans)
    end if
    allocate(trans(4,M_grid_size))
    read(10,*,iostat=error) string, trans
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'trans' from file '",trim(filename)
       stop 1
    end if
    ! equil_umsd_1
    if(allocated(equil_umsd_1)) then
       deallocate(equil_umsd_1)
    end if
    allocate(equil_umsd_1(n_part))
    read(10,*,iostat=error) string, equil_umsd_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'equil_umsd_1' from file '",trim(filename)
       stop 1
    end if
    ! equil_umsd_2
    if(allocated(equil_umsd_2)) then
       deallocate(equil_umsd_2)
    end if
    allocate(equil_umsd_2(n_part))
    read(10,*,iostat=error) string, equil_umsd_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'equil_umsd_2' from file '",trim(filename)
       stop 1
    end if
    ! sigma_equil_umsd_1
    if(allocated(sigma_equil_umsd_1)) then
       deallocate(sigma_equil_umsd_1)
    end if
    allocate(sigma_equil_umsd_1(n_part))
    read(10,*,iostat=error) string, sigma_equil_umsd_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'sigma_equil_umsd_1' from file '",trim(filename)
       stop 1
    end if
    ! sigma_equil_umsd_2
    if(allocated(sigma_equil_umsd_2)) then
       deallocate(sigma_equil_umsd_2)
    end if
    allocate(sigma_equil_umsd_2(n_part))
    read(10,*,iostat=error) string, sigma_equil_umsd_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'sigma_equil_umsd_2' from file '",trim(filename)
       stop 1
    end if
    ! intrablock_sum_umsd_1
    if(allocated(intrablock_sum_umsd_1)) then
       deallocate(intrablock_sum_umsd_1)
    end if
    allocate(intrablock_sum_umsd_1(n_part))
    read(10,*,iostat=error) string, intrablock_sum_umsd_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'intrablock_sum_umsd_1' from file '",trim(filename)
       stop 1
    end if
    ! intrablock_sum_umsd_2
    if(allocated(intrablock_sum_umsd_2)) then
       deallocate(intrablock_sum_umsd_2)
    end if
    allocate(intrablock_sum_umsd_2(n_part))
    read(10,*,iostat=error) string, intrablock_sum_umsd_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'intrablock_sum_umsd_2' from file '",trim(filename)
       stop 1
    end if
    ! interblock_sum_umsd_1
    if(allocated(interblock_sum_umsd_1)) then
       deallocate(interblock_sum_umsd_1)
    end if
    allocate(interblock_sum_umsd_1(n_part))
    read(10,*,iostat=error) string, interblock_sum_umsd_1
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_umsd_1' from file '",trim(filename)
       stop 1
    end if
    ! interblock_sum_umsd_2
    if(allocated(interblock_sum_umsd_2)) then
       deallocate(interblock_sum_umsd_2)
    end if
    allocate(interblock_sum_umsd_2(n_part))
    read(10,*,iostat=error) string, interblock_sum_umsd_2
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_umsd_2' from file '",trim(filename)
       stop 1
    end if
    ! interblock_sum_umsd_1_sqrd
    if(allocated(interblock_sum_umsd_1_sqrd)) then
       deallocate(interblock_sum_umsd_1_sqrd)
    end if
    allocate(interblock_sum_umsd_1_sqrd(n_part))
    read(10,*,iostat=error) string, interblock_sum_umsd_1_sqrd
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_umsd_1_sqrd' from file '",trim(filename)
       stop 1
    end if
    ! interblock_sum_umsd_2_sqrd
    if(allocated(interblock_sum_umsd_2_sqrd)) then
       deallocate(interblock_sum_umsd_2_sqrd)
    end if
    allocate(interblock_sum_umsd_2_sqrd(n_part))
    read(10,*,iostat=error) string, interblock_sum_umsd_2_sqrd
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem reading 'interblock_sum_umsd_2_sqrd' from file '",trim(filename)
       stop 1
    end if

    ! Import 'interactions' variables
    call import_interactions(Lx(1), Ly(1), Lz(1), spec_1, pos_1, R_1, Lx(2), Ly(2), Lz(2), spec_2, pos_2, R_2, u)

    close(unit=10)

  end subroutine import




  !! <h4> <code>  subroutine run(datafile,statefile,appenddata,seed) </code> </h4>
  !! <p>
  !! This subroutine runs the lattice switch Monte Carlo simulation using the current simulation variables.
  !! Data is output to the file 'datafile'; if the file already exists then is replaced or appended according
  !! to the <code>appenddata</code> flag. Data is also output to stdout. The state of the simulation is 
  !! periodically output to 'statefile'. Note that the heavy lifting is deferred to the 'move' procedures. 
  !! Be aware that this subroutine opens and closes units 10 and 11 (for the 'state' and 'data' files respectively).
  !! Note also that the variables described in the above section 'Variables determining the nature of the output'
  !! also determine the nature of the output: and certain values for these variables result in empty 'state'
  !! or 'data' files being created, or suppression of output to stdout.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> datafile </code> </td>
  !!   <td> <code> character(*), intent(in) </code> </td>
  !!   <td> The file to which data during the simulation will be output to. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> statefile </code> </td>
  !!   <td> <code> character(*), intent(in) </code> </td>
  !!   <td> 
  !!   The file to which the state of the simulation will be writen to periodically, and at the end of the
  !!   simulation
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> appenddata </code> </td>
  !!   <td> <code> logical, intent(in) </code> </td>
  !!   <td>
  !!   This should be <code>.true.</code> if <code>datafile</code> is to be
  !!   ammended, and <code>.false.</code> if it is to be overwritten. 
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> seed </code> </td>
  !!   <td> <code> integer, intent(in), optional </code> </td>
  !!   <td> A specific seed for the random number generator. If this argument is not present then a 
  !!   random seed is used - which is taken from the system clock.</td>
  !!  </tr>
  !! </table>
  subroutine run(datafile,statefile,appenddata,seed)
    character(*), intent(in) :: datafile
    character(*), intent(in) :: statefile
    logical, intent(in) :: appenddata
    integer, intent(in), optional :: seed
    integer(ik) :: i
    integer(ik) :: n
    integer(ik) :: p
    real(rk) :: vol_move_prob
    ! The time the simulation was started and the current time
    real :: time_start, time

    vol_move_prob=(1.0_rk*vol_freq)/n_part
    ! Initialise a random seed for the random number generator
    if(present(seed)) then
       call init_random_seed(seed)
    else
       call init_random_seed()
    end if

    call cpu_time(time_start)
    
    ! Open the unit corresponding to 'datafile'. Note that unit 11 is used instead
    ! of unit 10 since the subroutine 'export' used in the Monte Carlo loop uses unit
    ! 10. 
    if(appenddata) then
       open(unit=11,file=datafile,position="append")
    else
       open(unit=11,file=datafile,status="replace")
    end if

    ! If stop_sweeps=0 then do not run the Monte Carlo loop, but instead perform tasks which
    ! are normally performed periodically after sweeps, as well as calculate the equilibrium properties. 
    ! Otherwise proceed to the Monte Carlo loop
    if(stop_sweeps==0) then

       ! Check for melting
       if(enable_melt_checks) then
          call check_for_melt()
       end if
       ! Check energy divergence
       if(enable_divergence_checks) then
          call check_for_divergence()
       end if
       ! Update the weight function
       if(update_eta) then
          select case(update_eta_method)
          case(update_eta_method_VS)
             call update_eta_VS()
             M_counts_1=0
             M_counts_2=0
          case(update_eta_method_shooting)
             call update_eta_shooting()
          case default
             write(0,*) "monteswitch_mod: Error. 'update_eta_method' value is not recognised."
             stop 1
          end select
       end if
       ! Calculate equilibrium quantities
       if(sweeps>=sweep_equil_reference+equil_sweeps .and. calc_equil_properties) then
          call set_equil_properties()
       end if
       
    else
       
       ! THE LATTICE SWITCH MONTE CARLO LOOP

       ! Print preamble to stdout if appropriate: if output_stdout_period<0 then we don't.
       if(output_stdout_period>=0) then
          write(6,*)
          write(6,*) "****************************************************************"
          write(6,*)  
          write(6,*) "Starting lattice switch Monte Carlo loop"
          write(6,*)
          write(6,*)
       end if
       
       ! i is the sweep number for this simulation
       do i=1,stop_sweeps
          
          do n=1,n_part
             
             if(enable_part_moves) then
                ! Select the particle number 'p' to use for this cycle
                select case(part_select)
                case(part_select_cycle)
                   p=n
                case(part_select_rand)
                   p=int(get_random_number()*n_part)+1
                case default
                   write(0,*) "monteswitch_mod: Error. 'part_select' value is not recognised."
                   stop 1
                end select
                ! Move the selected particle
                call move_particle(p)
                call output_if_needed()
                
                ! Move the lattice
                if(enable_lattice_moves) then
                   call move_lattice()
                   call output_if_needed()
                end if
             else
                ! Move the lattice
                if(enable_lattice_moves) then
                   call move_lattice()
                   call output_if_needed()
                end if
             end if
             if(enable_vol_moves) then
                ! If we have an NPT ensemble move the volume on average 'vol_freq' times per sweep
                if(get_random_number()<vol_move_prob) then
                   ! Move the volume
                   call move_volume()
                   call output_if_needed()
                   ! Move the lattice
                   if(enable_lattice_moves) then
                      call move_lattice()
                      call output_if_needed()
                   end if
                end if
             end if
          end do
          
          sweeps=sweeps+1
          
          ! Perform some stuff after each sweep

          ! Check for melting when appropriate. Note that this should preceed checkpointing: otherwise
          ! a melted simualation would be checkpointed - which is not what we want.
          ! Note that if the system has melted, and the 'melt_option' variable is such that the simulation
          ! is not stopped, then after check_for_melt() is called the simulation has a different set
          ! of displacements, and 'sweep_equil_reference' has been set to the current 'sweeps'.
          if(enable_melt_checks) then
             if(mod(sweeps,melt_sweeps)==0) then
                call check_for_melt()
             end if
          end if
          ! Check for energy divergence every so often. Again this should precede checkpointing.
          if(enable_divergence_checks) then
             if(mod(sweeps,divergence_sweeps)==0) then
                call check_for_divergence()
             end if
          end if
          ! Print info to unit 11 and stdout when appropriate: if output_file_period=0 then we defer
          ! output to the output_if_needed() subroutine - which outputs every move; if output_file_period<0
          ! then we never output.
          if(output_file_period>0) then
             if(mod(sweeps,output_file_period)==0) then
                call output_file()
             end if
          end if
          if(output_stdout_period>0) then
             if(mod(sweeps,output_stdout_period)==0) then
                call output_stdout()
             end if
          end if
          ! Checkpoint when appropriate
          if(mod(sweeps,checkpoint_period)==0) then
             call export(statefile)
          end if
          ! Calculate equilibrium properties when appropriate (Note > instead of >= is not a typo
          ! even though >= is used elsewhere - which is for DURING a sweep, while > is for AFTER 
          ! a sweep: 'sweeps' was updated already above.
          if(sweeps>sweep_equil_reference+equil_sweeps .and. calc_equil_properties) then
             ! If we are at the end of a block, perform tasks which must be performed at the
             ! end of a block pertaining to the evaluation of equilibrium properties of the
             ! system - if appropriate.
             if(mod(sweeps-sweep_equil_reference-equil_sweeps,block_sweeps)==0) then
                call block_stuff()
             end if
          end if
          ! Update the weight function when appropriate
          if(update_eta) then
             if(mod(sweeps,update_eta_sweeps)==0) then
                select case(update_eta_method)
                case(update_eta_method_VS)
                   call update_eta_VS()
                   M_counts_1=0
                   M_counts_2=0
                case(update_eta_method_shooting)
                   call update_eta_shooting()
                case default
                   write(0,*) "monteswitch_mod: Error. 'update_eta_method' value is not recognised."
                   stop 1
                end select
             end if
          end if
          
       end do

       ! Print final info to stdout if appropriate: if output_stdout_period<0 then we don't.
       if(output_stdout_period>=0) then
          write(6,*) "****************************************************************"
          write(6,*)
          write(6,*) "Finished lattice switch Monte Carlo loop"
          write(6,*)
          call cpu_time(time)
          write(6,*) "CPU time elapsed this simulation (s) = ",(time-time_start)
          write(6,*) "Average time per sweep = ",(time-time_start)/stop_sweeps
          write(6,*)
       end if

    end if
    
    close(unit=11)

    call export(statefile)

  contains

    ! DETAILS OF HOW EQUILIBRIUM PROPERTIES ARE EVALUATED
    !
    ! This nested subroutine performs tasks related to evaluating the equilibrium variables - equil_DeltaF, equil_H_1,
    ! equil_H_2, equil_V_1, equil_V_2, sigma_equil_DeltaF, sigma_equil_H_1, sigma_equil_H_2, sigma_equil_V_1, 
    ! sigma_equil_V_2 - which must be performed at the end of every block.
    !
    ! In order to determine the equilibrium value and associated uncertainty of some property P which pertains to, say,
    ! phase 1, the variables intrablock_sum_P and intrablock_counts_1 are used to determine the value of P for the 
    ! current block. If the microstates were sampled according to the Boltzmann distribution, then if intrablock_sum_P were a running
    ! total of P after every move, and if intrablock_counts_1 were a running total of the number of times phase 1 was
    ! visited in the current block, then the value of P for the block would be
    ! intrablock_sum_P/intrablock_counts_1. However,
    ! if we have multicanonical sampling, then we must take into account the fact that our sampling is biased:
    ! microstate sigma is a factor of exp(eta(sigma))*Z/Z' more likely to be sampled in multicanonical sampling than
    ! Boltzmann sampling, where Z is the Boltzmann partition function and Z' is the multicanonical one. Therefore
    ! we oversample each visited state sigma by a factor exp(eta(sigma))*Z/Z'. This can be accounted for by giving
    ! microstate sigma a weight of exp(-eta(sigma)) in the evaluation, which amounts to ammending intrablock_sum_P
    ! and intrablock_counts_1 as follows (each time we visit state 1): intrablock_sum_P=intrablock_sum_P+exp(-eta(sigma))*P(sigma);
    ! intrablock_counts_1=intrablock_counts_1+exp(-eta(sigma)). Note that it is unnecessary to include the factor 
    ! Z/Z' since it is the same for all microstates, and hence vanishes when the block value of P, i.e. 
    ! intrablock_sum_P/intrablock_counts_P is evaluated.
    !
    ! Consider how the quantity DeltaF is evaluated. This is the difference in the free energies between
    ! phases 1 and 2 (defined as that of 1 minus that of 2), and is given by the formula:
    ! DeltaF=-(1/beta)*ln(p(1)/p(2)), where p(1) is the probability of the system being in phase 1 and p(2) is the
    ! probability of the system being in phase 2. Hence to evaluate DeltaF for a block we must evaluate the ratio
    ! p(1)/p(2) for the block. In Boltzmann sampling intrablock_counts_1 and intrablock_counts_2 can be used to 
    ! evaluate the ratio given that intrablock_counts_1 is the number of microstates visited for phase 1, and
    ! similarly for intrablock_counts_2: p(1)/p(2) is simply intrablock_counts_1/intrablock_counts_2. Again, for
    ! multicanonical sampling, we oversample each microstate sigma by a factor of exp(eta(sigma))*Z/Z'. As
    ! above, this can be accounted for by ammending intrablock_counts_1(2) by exp(-eta(sigma)) each time we visit
    ! a microstate belonging to phase 1(2) - instead of by 1 as is the case for Boltzmann sampling.
    !
    ! At the end of each block the properties pertaining to the block are evaluated as just described. Then they
    ! are added to running totals which are used to evaluate the equilibrium properties and their associated
    ! uncertainties. For property P: interblock_sum_P is a running total of P of the values obtained from all blocks
    ! so far; interblock_sum_P_sqrd is a running total of P*P of the values obtained from all blocks
    ! so far. The former is used to evaluate the equilibrium value of P as follows: equil_P=interblock_sum_P/block_counts_P,
    ! where block_counts_P is the number of blocks considered so far for property P. The latter is used to evaluate the 
    ! uncertainty in equil_P using the standard error of the mean: 
    ! sigma_equil_P=sqrt(interblock_sum_P_sqrt/block_counts_P-equil_P*equil_P)/sqrt(block_counts*block_counts_P).
    ! We only add a block to the running total if it yields a sensible value for P. E.g. if we are considering the enthalpy
    ! of phase 1 (H_1), and during the block we are never in phase 1, then H_1 will be nonsense, which will yield nonsence
    ! in the interblock running totals if used. It is for this reason that we have separate 'block_counts' counters for each 
    ! property.
    !
    ! Of course, all intrablock counters and sums must be reset at the end of each block.
    !
    subroutine block_stuff()
      ! The values of DeltaF, H_1, H_2, V_1, V_2, L_1 and L_2 for the current block
      real(rk) :: block_DeltaF
      real(rk) :: block_H_1
      real(rk) :: block_H_2
      real(rk) :: block_V_1
      real(rk) :: block_V_2
      real(rk), dimension(3) :: block_L_1
      real(rk), dimension(3) :: block_L_2
      ! The values of the arrays umsd_1 and umsd_2 for the current block
      real(rk), dimension(n_part) :: block_umsd_1
      real(rk), dimension(n_part) :: block_umsd_2

      ! Update block_counts
      block_counts=block_counts+1
      ! Code for DeltaF: consider the block only if intrablock_counts_1 and intrablock_counts_2 are >0
      if(intrablock_counts_1>0 .and. intrablock_counts_2>0) then
         block_DeltaF=-log(intrablock_counts_1/intrablock_counts_2)/beta
         block_counts_DeltaF=block_counts_DeltaF+1
         interblock_sum_DeltaF=interblock_sum_DeltaF+block_DeltaF
         interblock_sum_DeltaF_sqrd=interblock_sum_DeltaF_sqrd + block_DeltaF*block_DeltaF
      end if
      ! Code for H_1 and H_2: consider the block only if intrablock_counts_1>0 and intrablock_counts_2>0 respectively
      if(intrablock_counts_1>0) then
         block_H_1=intrablock_sum_H_1/intrablock_counts_1
         block_counts_H_1=block_counts_H_1+1
         interblock_sum_H_1=interblock_sum_H_1+block_H_1
         interblock_sum_H_1_sqrd=interblock_sum_H_1_sqrd + block_H_1*block_H_1
      end if
      if(intrablock_counts_2>0) then
         block_H_2=intrablock_sum_H_2/intrablock_counts_2
         block_counts_H_2=block_counts_H_2+1
         interblock_sum_H_2=interblock_sum_H_2+block_H_2
         interblock_sum_H_2_sqrd=interblock_sum_H_2_sqrd + block_H_2*block_H_2
      end if
      ! Code for V_1 and V_2: consider the block only if intrablock_counts_1>0 and intrablock_counts_2>0 respectively
      if(intrablock_counts_1>0) then
         block_V_1=intrablock_sum_V_1/intrablock_counts_1
         block_counts_V_1=block_counts_V_1+1
         interblock_sum_V_1=interblock_sum_V_1+block_V_1
         interblock_sum_V_1_sqrd=interblock_sum_V_1_sqrd + block_V_1*block_V_1
      end if
      if(intrablock_counts_2>0) then
         block_V_2=intrablock_sum_V_2/intrablock_counts_2
         block_counts_V_2=block_counts_V_2+1
         interblock_sum_V_2=interblock_sum_V_2+block_V_2
         interblock_sum_V_2_sqrd=interblock_sum_V_2_sqrd + block_V_2*block_V_2
      end if
      ! Code for umsd_1 and umsd_2: consider the block only if intrablock_counts_1>0 and intrablock_counts_2>0 respectively
      if(intrablock_counts_1>0) then
         block_umsd_1=intrablock_sum_umsd_1/intrablock_counts_1
         block_counts_umsd_1=block_counts_umsd_1+1
         interblock_sum_umsd_1=interblock_sum_umsd_1+block_umsd_1
         interblock_sum_umsd_1_sqrd=interblock_sum_umsd_1_sqrd + block_umsd_1*block_umsd_1
      end if
      if(intrablock_counts_2>0) then
         block_umsd_2=intrablock_sum_umsd_2/intrablock_counts_2
         block_counts_umsd_2=block_counts_umsd_2+1
         interblock_sum_umsd_2=interblock_sum_umsd_2+block_umsd_2
         interblock_sum_umsd_2_sqrd=interblock_sum_umsd_2_sqrd + block_umsd_2*block_umsd_2
      end if
      ! Code for L_1 and L_2: consider the block only if intrablock_counts_1>0 and intrablock_counts_2>0 respectively
      if(intrablock_counts_1>0) then
         block_L_1=intrablock_sum_L_1/intrablock_counts_1
         block_counts_L_1=block_counts_L_1+1
         interblock_sum_L_1=interblock_sum_L_1+block_L_1
         interblock_sum_L_1_sqrd=interblock_sum_L_1_sqrd + block_L_1*block_L_1
      end if
      if(intrablock_counts_2>0) then
         block_L_2=intrablock_sum_L_2/intrablock_counts_2
         block_counts_L_2=block_counts_L_2+1
         interblock_sum_L_2=interblock_sum_L_2+block_L_2
         interblock_sum_L_2_sqrd=interblock_sum_L_2_sqrd + block_L_2*block_L_2
      end if

      ! Reset sums
      intrablock_counts_1=0.0_rk
      intrablock_counts_2=0.0_rk
      intrablock_sum_H_1=0.0_rk
      intrablock_sum_H_2=0.0_rk
      intrablock_sum_V_1=0.0_rk
      intrablock_sum_V_2=0.0_rk
      intrablock_sum_umsd_1=0.0_rk
      intrablock_sum_umsd_2=0.0_rk
      intrablock_sum_L_1=0.0_rk
      intrablock_sum_L_2=0.0_rk
      ! Calculate equilibrium properties
      call set_equil_properties()

    end subroutine block_stuff
    

    ! This nested subroutine checks whether or not the system has melted and acts accordingly.
    subroutine check_for_melt()

      ! Check for melting first
      if(max(abs(maxval(u)),abs(minval(u)))>melt_threshold) then
         melts=melts+1

         select case(melt_option)
         case(melt_option_stop)
            write(0,*) "monteswitch_mod: Error. System has melted."
            call export("state_ERROR")
            stop 2
         case(melt_option_zero_1,melt_option_zero_2,melt_option_zero_current,melt_option_zero_random)
            u=0.0_rk
            pos_1=R_1
            pos_2=R_2
            E_1=calc_energy_scratch(1,Lx(1),Ly(1),Lz(1),spec_1,pos_1,R_1,u)
            E_2=calc_energy_scratch(2,Lx(2),Ly(2),Lz(2),spec_2,pos_2,R_2,u)
            select case(melt_option)
            case(melt_option_zero_1)
               lattice=1
               E=E_1
            case(melt_option_zero_2)
               lattice=2
               E=E_2
            case(melt_option_zero_current)
               select case(lattice)
               case(1)
                  E=E_1
               case(2)
                  E=E_2
               end select
            case(melt_option_zero_random)
               if(get_random_number()<0.5) then
                  lattice=1
               else
                  lattice=2
               end if
               select case(lattice)
               case(1)
                  E=E_1
               case(2)
                  E=E_2
               end select
            end select
            M=calc_M(E_1, E_2, Lx(1)*Ly(1)*Lz(1), Lx(2)*Ly(2)*Lz(2))
            macro=get_macro(M)
            ! If M is outwith the supported range then get_macro will return -2 or -1.
            if(macro<1) then
                write(0,*) "monteswitch_mod: Error. Order parameter after melt reset is outwith supported range: M = ",M
                stop 1
            end if
            eta=eta_grid(macro)
            if(enable_barriers) then
               call initialise_barriers()
            end if
            ! Set 'sweep_equil_reference' such that the system re-equilibrates, 
            ! and reset intrablock variables
            sweep_equil_reference=sweeps
            intrablock_counts_1=0.0_rk
            intrablock_counts_2=0.0_rk
            intrablock_sum_H_1=0.0_rk
            intrablock_sum_H_2=0.0_rk
            intrablock_sum_V_1=0.0_rk
            intrablock_sum_V_2=0.0_rk
         case default
            write(0,*) "monteswitch_mod: Error. 'melt_option' not recognised."
            stop 1
         end select
      
      end if

    end subroutine check_for_melt

    
    ! This nested subroutine checks whether or not the energy has diverged and resets the energies and order parameter
    ! to their 'exact' values. 
    subroutine check_for_divergence()
      real(rk) :: E_exact_1, E_exact_2
      ! Calculate the exact energy for both lattice types
      E_exact_1=calc_energy_scratch(1,Lx(1),Ly(1),Lz(1),spec_1,pos_1,R_1,u)
      E_exact_2=calc_energy_scratch(2,Lx(2),Ly(2),Lz(2),spec_2,pos_2,R_2,u)
      ! Compare to the current energies and flag any errors
      if(abs(E_1-E_exact_1)>divergence_tol) then
         write(0,*) "monteswitch_mod: Error. Energy of lattice 1 has diverged from exact value by ",abs(E_1-E_exact_1),&
              "; E_1= ",E_1,"; E_exact_1= ",E_exact_1
         call export("state_ERROR")
         stop 3
      end if
      if(abs(E_2-E_exact_2)>divergence_tol) then
         write(0,*) "monteswitch_mod: Error. Energy of lattice 2 has diverged from exact value by ",abs(E_2-E_exact_2),&
              "; E_2= ",E_2,"; E_exact_2= ",E_exact_2
         call export("state_ERROR")
         stop 3
      end if
      ! Set the energies and order parameters to their exact values
      E_1=E_exact_1
      E_2=E_exact_2
      M=calc_M(E_1, E_2, Lx(1)*Ly(1)*Lz(1), Lx(2)*Ly(2)*Lz(2))
      macro=get_macro(M)
      ! If M is outwith the supported range then get_macro will return -2 or -1.
      if(macro<1) then
          write(0,*) "monteswitch_mod: Error. Order parameter after divergence check is outwith supported range: M = ",M
          stop 1
      end if
      eta=eta_grid(macro)
    end subroutine check_for_divergence


    ! This nested subroutine performs any necessary output to be performed after every move, which occurs if
    ! output_file_period==0 or output_stdout_period==0
    subroutine output_if_needed()
      if(output_file_period==0) then
         call output_file()
      end if
      if(output_stdout_period==0) then
         call output_stdout()
      end if
    end subroutine output_if_needed


    ! This nested subroutine outputs the required data to the file 'datafile'
    subroutine output_file()
      if(output_file_Lx) then
         write(11,*) "Lx: ",sweeps,Lx
      end if 
      if(output_file_Ly) then
         write(11,*) "Ly: ",sweeps,Ly
      end if 
      if(output_file_Lz) then
         write(11,*) "Lz: ",sweeps,Lz
      end if
      if(output_file_V) then
         write(11,*) "V: ",sweeps,V
      end if
      if(output_file_R_1) then
         write(11,*) "R_1: ",sweeps,R_1
      end if
      if(output_file_R_2) then
         write(11,*) "R_2: ",sweeps,R_2
      end if
      if(output_file_u) then
         write(11,*) "u: ",sweeps,u
      end if
      if(output_file_lattice) then
         write(11,*) "lattice: ",sweeps,lattice
      end if
      if(output_file_E) then
         write(11,*) "E: ",sweeps,E
      end if
      if(output_file_M) then
         write(11,*) "M: ",sweeps,M
      end if  
      if(output_file_macro) then
         write(11,*) "macro: ",sweeps,macro
      end if  
      if(output_file_eta) then
         write(11,*) "eta: ",sweeps,eta
      end if  
      if(output_file_moves_lattice) then
         write(11,*) "moves_lattice: ",sweeps,moves_lattice
      end if
      if(output_file_accepted_moves_lattice) then
         write(11,*) "accepted_moves_lattice: ",sweeps,accepted_moves_lattice
      end if      
      if(output_file_moves_part) then
         write(11,*) "moves_part: ",sweeps,moves_part
      end if
      if(output_file_accepted_moves_part) then
         write(11,*) "accepted_moves_part: ",sweeps,accepted_moves_part
      end if      
      if(output_file_moves_vol) then
         write(11,*) "moves_vol: ",sweeps,moves_vol
      end if
      if(output_file_accepted_moves_vol) then
         write(11,*) "accepted_moves_vol: ",sweeps,accepted_moves_vol
      end if
      if(output_file_rejected_moves_M_OOB) then
         write(11,*) "rejected_moves_M_OOB: ",sweeps,rejected_moves_M_OOB
      end if
      if(output_file_M_OOB_high) then
         write(11,*) "M_OOB_high: ",sweeps,M_OOB_high
      end if
      if(output_file_M_OOB_low) then
         write(11,*) "M_OOB_low: ",sweeps,M_OOB_low
      end if
      if(output_file_barrier_macro_low) then
         write(11,*) "barrier_macro_low: ",sweeps,barrier_macro_low
      end if
      if(output_file_barrier_macro_high) then
         write(11,*) "barrier_macro_high: ",sweeps,barrier_macro_high
      end if
      if(output_file_rejected_moves_M_barrier) then
         write(11,*) "rejected_moves_M_barrier: ",sweeps,rejected_moves_M_barrier
      end if
      if(output_file_moves_since_lock) then
         write(11,*) "moves_since_lock: ",sweeps,moves_since_lock
      end if
      if(output_file_melts) then
         write(11,*) "melts: ",sweeps,melts
      end if
      if(output_file_equil_DeltaF) then
         write(11,*) "equil_DeltaF: ",sweeps,equil_DeltaF
      end if
      if(output_file_sigma_equil_DeltaF) then
         write(11,*) "sigma_equil_DeltaF: ",sweeps,sigma_equil_DeltaF
      end if
      if(output_file_equil_H_1) then
         write(11,*) "equil_H_1: ",sweeps,equil_H_1
      end if
      if(output_file_sigma_equil_H_1) then
         write(11,*) "sigma_equil_H_1: ",sweeps,sigma_equil_H_1
      end if
      if(output_file_equil_H_2) then
         write(11,*) "equil_H_2: ",sweeps,equil_H_2
      end if
      if(output_file_sigma_equil_H_2) then
         write(11,*) "sigma_equil_H_2: ",sweeps,sigma_equil_H_2
      end if
      if(output_file_equil_V_1) then
         write(11,*) "equil_V_1: ",sweeps,equil_V_1
      end if
      if(output_file_sigma_equil_V_1) then
         write(11,*) "sigma_equil_V_1: ",sweeps,sigma_equil_V_1
      end if
      if(output_file_equil_V_2) then
         write(11,*) "equil_V_2: ",sweeps,equil_V_2
      end if
      if(output_file_sigma_equil_V_2) then
         write(11,*) "sigma_equil_V_2: ",sweeps,sigma_equil_V_2
      end if
      if(output_file_equil_umsd_1) then
         write(11,*) "equil_umsd_1: ",sweeps,equil_umsd_1
      end if
      if(output_file_sigma_equil_umsd_1) then
         write(11,*) "sigma_equil_umsd_1: ",sweeps,sigma_equil_umsd_1
      end if
      if(output_file_equil_umsd_2) then
         write(11,*) "equil_umsd_2: ",sweeps,equil_umsd_2
      end if
      if(output_file_sigma_equil_umsd_2) then
         write(11,*) "sigma_equil_umsd_2: ",sweeps,sigma_equil_umsd_2
      end if
      if(output_file_equil_L_1) then
         write(11,*) "equil_L_1: ",sweeps,equil_L_1
      end if
      if(output_file_sigma_equil_L_1) then
         write(11,*) "sigma_equil_L_1: ",sweeps,sigma_equil_L_1
      end if
      if(output_file_equil_L_2) then
         write(11,*) "equil_L_2: ",sweeps,equil_L_2
      end if
      if(output_file_sigma_equil_L_2) then
         write(11,*) "sigma_equil_L_2: ",sweeps,sigma_equil_L_2
      end if
    end subroutine output_file

    
    ! This nested subroutine outputs the required info to stdout
    subroutine output_stdout()

      call cpu_time(time)
      write(6,*) "****************************************************************"
      write(6,*)
      write(6,*) "Completed ",i," sweeps of ",stop_sweeps," for this simulation"
      write(6,*) "(Completed ",sweeps," sweeps over all simulations)"
      write(6,*)
      write(6,*) "CPU time elapsed this simulation (s) = ",(time-time_start)
      write(6,*) "Average time per sweep so far (s) this simulation = ",(time-time_start)/i
      write(6,*)
      write(6,*) "Current information:"

      if(output_stdout_Lx) then
         write(6,*) "Lx: ",Lx
      end if
      if(output_stdout_Ly) then
         write(6,*) "Ly: ",Ly
      end if
      if(output_stdout_Lz) then
         write(6,*) "Lz: ",Lz
      end if
      if(output_stdout_V) then
         write(6,*) "V: ",V
      end if
      if(output_stdout_R_1) then
         write(6,*) "R_1: ",R_1
      end if
      if(output_stdout_R_2) then
         write(6,*) "R_2: ",R_2
      end if
      if(output_stdout_u) then
         write(6,*) "u: ",u
      end if
      if(output_stdout_lattice) then
         write(6,*) "lattice: ",lattice
      end if
      if(output_stdout_E) then
         write(6,*) "E: ",E
      end if
      if(output_stdout_M) then
         write(6,*) "M: ",M
      end if
      if(output_stdout_macro) then
         write(6,*) "macro: ",macro
      end if
      if(output_stdout_eta) then
         write(6,*) "eta: ",eta
      end if
      if(output_stdout_moves_lattice) then
         write(6,*) "moves_lattice: ",moves_lattice
      end if
      if(output_stdout_accepted_moves_lattice) then
         write(6,*) "accepted_moves_lattice: ",accepted_moves_lattice
      end if
      if(output_stdout_moves_part) then
         write(6,*) "moves_part: ",moves_part
      end if
      if(output_stdout_accepted_moves_part) then
         write(6,*) "accepted_moves_part: ",accepted_moves_part
      end if
      if(output_stdout_moves_vol) then
         write(6,*) "moves_vol: ",moves_vol
      end if
      if(output_stdout_accepted_moves_vol) then
         write(6,*) "accepted_moves_vol: ",accepted_moves_vol
      end if
      if(output_stdout_rejected_moves_M_OOB) then
         write(6,*) "rejected_moves_M_OOB: ",rejected_moves_M_OOB
      end if
      if(output_stdout_M_OOB_high) then
         write(6,*) "M_OOB_high: ",M_OOB_high
      end if
      if(output_stdout_M_OOB_low) then
         write(6,*) "M_OOB_low: ",M_OOB_low
      end if
      if(output_stdout_barrier_macro_low) then
         write(6,*) "barrier_macro_low: ",barrier_macro_low
      end if
      if(output_stdout_barrier_macro_high) then
         write(6,*) "barrier_macro_high: ",barrier_macro_high
      end if
      if(output_stdout_rejected_moves_M_barrier) then
         write(6,*) "rejected_moves_M_barrier: ",rejected_moves_M_barrier
      end if
      if(output_stdout_moves_since_lock) then
         write(6,*) "moves_since_lock: ",moves_since_lock
      end if
      if(output_stdout_melts) then
         write(6,*) "melts: ",melts
      end if
      if(output_stdout_equil_DeltaF) then
         write(6,*) "equil_DeltaF: ",equil_DeltaF
      end if
      if(output_stdout_sigma_equil_DeltaF) then
         write(6,*) "sigma_equil_DeltaF: ",sigma_equil_DeltaF
      end if
      if(output_stdout_equil_H_1) then
         write(6,*) "equil_H_1: ",equil_H_1
      end if
      if(output_stdout_sigma_equil_H_1) then
         write(6,*) "sigma_equil_H_1: ",sigma_equil_H_1
      end if
      if(output_stdout_equil_H_2) then
         write(6,*) "equil_H_2: ",equil_H_2
      end if
      if(output_stdout_sigma_equil_H_2) then
         write(6,*) "sigma_equil_H_2: ",sigma_equil_H_2
      end if
      if(output_stdout_equil_V_1) then
         write(6,*) "equil_V_1: ",equil_V_1
      end if
      if(output_stdout_sigma_equil_V_1) then
         write(6,*) "sigma_equil_V_1: ",sigma_equil_V_1
      end if
      if(output_stdout_equil_V_2) then
         write(6,*) "equil_V_2: ",equil_V_2
      end if
      if(output_stdout_sigma_equil_V_2) then
         write(6,*) "sigma_equil_V_2: ",sigma_equil_V_2
      end if
      if(output_stdout_equil_umsd_1) then
         write(6,*) "equil_umsd_1: ",equil_umsd_1
      end if
      if(output_stdout_sigma_equil_umsd_1) then
         write(6,*) "sigma_equil_umsd_1: ",sigma_equil_umsd_1
      end if
      if(output_stdout_equil_umsd_2) then
         write(6,*) "equil_umsd_2: ",equil_umsd_2
      end if
      if(output_stdout_sigma_equil_umsd_2) then
         write(6,*) "sigma_equil_umsd_2: ",sigma_equil_umsd_2
      end if
      if(output_stdout_equil_L_1) then
         write(6,*) "equil_L_1: ",equil_L_1
      end if
      if(output_stdout_sigma_equil_L_1) then
         write(6,*) "sigma_equil_L_1: ",sigma_equil_L_1
      end if
      if(output_stdout_equil_L_2) then
         write(6,*) "equil_L_2: ",equil_L_2
      end if
      if(output_stdout_sigma_equil_L_2) then
         write(6,*) "sigma_equil_L_2: ",sigma_equil_L_2
      end if
    end subroutine output_stdout

  end subroutine run




  !! <h3> Procedures to update weight functions </h3>




  !! <h4> <code> subroutine update_eta_VS() </code> </h4>
  !! <p>
  !! This subroutine updates <code>eta_grid</code> according to the information in the order
  !! parameter histograms, i.e. the <code>M_counts_1</code> and <code>M_counts_2</code> arrays. Specifically, an estimate
  !! <code>Pmac(M)</code> of the probability of the system being in macrostate <code>M</code> is
  !! created from the total number of visits so far to macrostate <code>M</code>, i.e.
  !! <code>M_counts_1(M)+M_counts_2(M)</code>, using the formula
  !! <code>Pmac(M)=(M_counts_1(M)+M_counts_2(M)+1)/sum_{M'}(M_counts_1(M')+M_counts_2(M')+1)</code>.
  !! Using this estimate, the weight function is updated as follows:
  !! <code>eta_grid(M)=eta_grid(M)-log(Pmac(M))+k</code>, where <code>k</code> is chosen such that
  !! the minimum value of <code>eta_grid(M)</code> over all <code>M</code> is 0.
  !! </p>
  subroutine update_eta_VS()
    ! Pmac(M) is the (inferred) probability of the system being in macrostate M
    real(rk), dimension(M_grid_size) :: Pmac
    Pmac=(M_counts_1+M_counts_2+1.0_rk)/sum(M_counts_1+M_counts_2+1.0_rk)
    eta_grid=eta_grid-log(Pmac)
    eta_grid=eta_grid-minval(eta_grid)
  end subroutine update_eta_VS




  !! <h4> <code> subroutine update_eta_shooting() </code> </h4>
  !! <p>
  !! This subroutine updates <code>eta_grid</code> according to the shooting method. The
  !! 'shooting' is started from the macrostate with the lowest order parameter.
  !! </p>
  subroutine update_eta_shooting()

    integer(ik) :: i
    
    eta_grid = 0.0_rk
    
    do i=1, M_grid_size-1

        eta_grid(i+1) = eta_grid(i) - log( &
	         ( (trans(3,i)+1.0_rk) / (trans(1,i+1)+1.0_rk)  ) * &
                 ( (sum(trans(:,i+1))+1.0_rk*M_grid_size) / (sum(trans(:,i))+1.0_rk*M_grid_size) )   &
                )
        
    end do

    eta_grid = eta_grid - minval(eta_grid)

    
  end subroutine update_eta_shooting




  !! <h3> Utility procedures </h3>




  !! <h4> <code> subroutine set_equil_properties() </code> </h4>
  !! <p>
  !! This subroutine sets all <code>equil_</code> and <code>sigma_</code> variables from their associated interblock sums and
  !! 'block counts'.
  !! </p>
  subroutine set_equil_properties()
    ! DeltaF
    equil_DeltaF=interblock_sum_DeltaF/block_counts_DeltaF
    sigma_equil_DeltaF=sqrt( interblock_sum_DeltaF_sqrd/block_counts_DeltaF - equil_DeltaF*equil_DeltaF ) &
         / sqrt(block_counts_DeltaF-1.0_rk)
    ! H_1
    equil_H_1=interblock_sum_H_1/block_counts_H_1
    sigma_equil_H_1=sqrt( interblock_sum_H_1_sqrd/block_counts_H_1 - equil_H_1*equil_H_1 ) &
         / sqrt(block_counts_H_1-1.0_rk)
    ! H_2
    equil_H_2=interblock_sum_H_2/block_counts_H_2
    sigma_equil_H_2=sqrt( interblock_sum_H_2_sqrd/block_counts_H_2 - equil_H_2*equil_H_2 ) &
         / sqrt(block_counts_H_2-1.0_rk)
    ! V_1
    equil_V_1=interblock_sum_V_1/block_counts_V_1
    sigma_equil_V_1=sqrt( interblock_sum_V_1_sqrd/block_counts_V_1 - equil_V_1*equil_V_1 ) &
         / sqrt(block_counts_V_1-1.0_rk)
    ! V_2
    equil_V_2=interblock_sum_V_2/block_counts_V_2
    sigma_equil_V_2=sqrt( interblock_sum_V_2_sqrd/block_counts_V_2 - equil_V_2*equil_V_2 ) &
         / sqrt(block_counts_V_2-1.0_rk)
    ! umsd_1
    equil_umsd_1=interblock_sum_umsd_1/block_counts_umsd_1
    sigma_equil_umsd_1=sqrt( interblock_sum_umsd_1_sqrd/block_counts_umsd_1 - equil_umsd_1*equil_umsd_1 ) &
         / sqrt(block_counts_umsd_1-1.0_rk)
    ! umsd_2
    equil_umsd_2=interblock_sum_umsd_2/block_counts_umsd_2
    sigma_equil_umsd_2=sqrt( interblock_sum_umsd_2_sqrd/block_counts_umsd_2 - equil_umsd_2*equil_umsd_2 ) &
         / sqrt(block_counts_umsd_2-1.0_rk)
    ! L_1
    equil_L_1=interblock_sum_L_1/block_counts_L_1
    sigma_equil_L_1=sqrt( interblock_sum_L_1_sqrd/block_counts_L_1 - equil_L_1*equil_L_1 ) &
         / sqrt(block_counts_L_1-1.0_rk)
    ! L_2
    equil_L_2=interblock_sum_L_2/block_counts_L_2
    sigma_equil_L_2=sqrt( interblock_sum_L_2_sqrd/block_counts_L_2 - equil_L_2*equil_L_2 ) &
         / sqrt(block_counts_L_2-1.0_rk)
  end subroutine set_equil_properties




  !! <h4> <code> function calc_M(E_1,E_2,V_1,V_2) </code> </h4>
  !! <p>
  !! This function returns the order parameter if the energies and volumes of the system in lattice types
  !! 1 and 2 are those of the arguments.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> E_1 </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   Energy corresponding to lattice type 1.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> E_2 </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   Energy corresponding to lattice type 2.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> V_1 </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   Volume corresponding to lattice type 1.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> V_2 </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   Volume corresponding to lattice type 2.
  !!   </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
  function calc_M(E_1,E_2,V_1,V_2)
    real(rk), intent(in) :: E_1
    real(rk), intent(in) :: E_2
    real(rk), intent(in) :: V_1
    real(rk), intent(in) :: V_2
    real(rk) :: calc_M
    calc_M = E_1-E_2 + P*(V_1-V_2) - M_const
  end function calc_M




  
  !! <h4> <code> function get_macro(M) </code> </h4>
  !! <p>
  !! This function returns the 'macrostate number' which <code>M</code> corresponds to. If M is within the supported
  !! range the function returns a value ><code>M_grid_size</code> or <1 then <code>M</code> as one would expect. If
  !! M is above the supported range then a value of -1 is returned; if M is below the supported range a value of -2
  !! is returned.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> M </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   The order parameter which the macrostate number is required for.
  !!   </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> integer(ik) </code> </p>
 function get_macro(M)
    real(rk), intent(in) :: M
    integer(ik) :: get_macro
    real(rk) :: dM

    ! Macrostate order parameter window size
    dM = M_grid(2)-M_grid(1)
    
    ! Check M is within the allowed range.
    ! Note that I do not simply use 'get_macro= floor( (M-M_grid(1)) / dM )+1' to calculate get_macro,
    ! and then check whether get_macro<1 or get_macro>M_grid_size to determine whether get_macro is
    ! within the supported range. Why? Because if M is very large then this expression can return 
    ! unexpected values for get_macro, possibly because of some two's compliment thing - I can't quite
    ! remember.
    if(M >= M_grid(M_grid_size)+dM) then
        ! M is above allowed range: return a value of -1
        get_macro = -1
    else if(M < M_grid(1)) then
        ! M is below allowed range: return a value of -2
        get_macro = -2
    else
        ! M is within allowed range
        get_macro= floor( (M-M_grid(1)) / dM )+1
    end if

  end function get_macro



  !! <h4> <code> function metropolis_prob(E_trial,macro_trial) </code> </h4>
  !! <p>
  !! This function gives the probability of a particle move with the specified trial energy and
  !! macrostate being accepted according to the Metropolis algorithm, and depends
  !! on the type of sampling being used (i.e. the value of the <code>enable_multicanonical</code> flag).
  !! The <code>trans</code> array is updated within
  !! this function if it is in use (i.e. if <code>update_trans=.true.</code>) and we have reached equilibration.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> E_trial </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   The trial energy.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> macro_trial </code> </td>
  !!   <td> <code> integer(ik), intent(in) </code> </td>
  !!   <td>
  !!   The trial macrostate number.
  !!   </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
  function metropolis_prob(E_trial,macro_trial)
    real(rk), intent(in) :: E_trial
    integer(ik), intent(in) :: macro_trial
    real(rk) :: metropolis_prob
    real(rk) :: prob_B
    if(enable_multicanonical) then
       metropolis_prob=metropolis_prob_MC(beta,E,E_trial,eta_grid(macro),eta_grid(macro_trial))
       if(update_trans .and. sweeps>=sweep_equil_reference+equil_sweeps) then
          ! Update 'trans' using the Boltzmann probability of the transition
          prob_B=metropolis_prob_B(beta,E,E_trial)
          call trans_update(macro,macro_trial,prob_B)
       end if
    else
       metropolis_prob=metropolis_prob_B(beta,E,E_trial)
       if(update_trans .and. sweeps>=sweep_equil_reference+equil_sweeps) then
          call trans_update(macro,macro_trial,metropolis_prob)
       end if
    end if
  end function metropolis_prob



  !! <h4> <code> subroutine trans_update(macro,macro_trial,prob) </code> </h4>
  !! <p>
  !! This procedure updates <code>trans</code> given that we are currently in macrostate <code>macro</code>, and
  !! have proposed a trial move to macrostate <code>macro_trial</code> whose Boltzman/canonical/unbiased
  !! probability of it being accepted is <code>prob</code>.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> macro </code> </td>
  !!   <td> <code> integer(ik), intent(in) </code> </td>
  !!   <td>
  !!   The current macrostate
  !!   </td>
  !!  </tr>
 !!  <tr>
  !!   <td> <code> macro_trial </code> </td>
  !!   <td> <code> integer(ik), intent(in) </code> </td>
  !!   <td>
  !!   The trial macrostate
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> prob </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   The probability of the transition being successful/accepted.
  !!   </td>
  !!  </tr>
  !! </table>
  subroutine trans_update(macro, macro_trial, prob)

      integer(ik), intent(in) :: macro
      integer(ik), intent(in) :: macro_trial
      real(rk), intent(in) :: prob

      integer(ik) :: macro_diff

      macro_diff = macro_trial - macro

      select case(macro_diff)
      case(-1,0,1)
         ! The trial macrostate is a nearest neighbour macro state, or the same state
         macro_diff = macro_diff + 2
         trans(macro_diff,macro) = trans(macro_diff,macro) + prob
         trans(2,macro) = trans(2,macro) + 1.0_rk - prob
      case default
         ! The trial macrostate is something else
         trans(4,macro) = trans(4,macro) + prob
         trans(2,macro) = trans(2,macro) + 1.0_rk - prob
      end select

  end subroutine trans_update
      
  


  !! <h4> <code> function metropolis_prob_vol(E_trial,macro_trial,V_trial) </code> </h4>
  !! <p>
  !! This function gives the probability of a volume move with the specified trial energy,
  !! macrostate and volume being accepted according to the Metropolis algorithm, and depends
  !! on the type of sampling being used (i.e. the value of the <code>enable_multicanonical</code> flag). The 
  !! <code>trans</code> array is updated within this function if it is in use (if <code>update_trans=.true.</code>)
  !! and we have reached equilibration.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> E_trial </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   The trial energy.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> macro_trial </code> </td>
  !!   <td> <code> integer(ik), intent(in) </code> </td>
  !!   <td>
  !!   The trial macrostate number.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> V_trial </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   The trial volume.
  !!   </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
  function metropolis_prob_vol(E_trial,macro_trial,V_trial)
    real(rk), intent(in) :: E_trial
    integer(ik), intent(in) :: macro_trial
    real(rk), intent(in) :: V_trial
    real(rk) :: metropolis_prob_vol
    real(rk) :: prob_B
    if(enable_multicanonical) then
        select case(vol_dynamics)
        case(vol_dynamics_FVM)
            metropolis_prob_vol = &
                metropolis_prob_MC_vol(beta,E,E_trial,eta_grid(macro),eta_grid(macro_trial),n_part,P,V,V_trial)
            if(update_trans .and. sweeps>=sweep_equil_reference+equil_sweeps) then
                ! Update 'trans' using the Boltzmann probability of the transition
                prob_B=metropolis_prob_B_vol(beta,E,E_trial,n_part,P,V,V_trial)
                call trans_update(macro,macro_trial,prob_B)
            end if
        case(vol_dynamics_UVM)
            metropolis_prob_vol = &
                metropolis_prob_MC_vol_uniaxial(beta,E,E_trial,eta_grid(macro),eta_grid(macro_trial),n_part,P,V,V_trial)
            if(update_trans .and. sweeps>=sweep_equil_reference+equil_sweeps) then
                ! Update 'trans' using the Boltzmann probability of the transition
                prob_B=metropolis_prob_B_vol_uniaxial(beta,E,E_trial,n_part,P,V,V_trial)
                call trans_update(macro,macro_trial,prob_B)
            end if
        case default
            write(0,*) "monteswitch_mod: Error. 'vol_dynamics' value is not recognised."
            stop 1
        end select
    else
        select case(vol_dynamics)
        case(vol_dynamics_FVM)
            metropolis_prob_vol=metropolis_prob_B_vol(beta,E,E_trial,n_part,P,V,V_trial)
        case(vol_dynamics_UVM)
            metropolis_prob_vol=metropolis_prob_B_vol_uniaxial(beta,E,E_trial,n_part,P,V,V_trial)
        case default
            write(0,*) "monteswitch_mod: Error. 'vol_dynamics' value is not recognised."
            stop 1
        end select
        if(update_trans .and. sweeps>=sweep_equil_reference+equil_sweeps) then
            call trans_update(macro,macro_trial,metropolis_prob_vol)
        end if
    end if
  end function metropolis_prob_vol




  !! <h4> <code> function metropolis_prob_lattice(E_trial,macro_trial,V_trial) </code> </h4>
  !! <p>
  !! This function gives the probability of a lattice move with the specified trial energy,
  !! macrostate and volume being accepted according to the Metropolis algorithm. Note that
  !! the acceptance probability is not the same as for a volume move because displacements are held
  !! constant during a lattice switch here. The probability returned depends
  !! on the type of sampling being used (i.e. the value of the <code>enable_multicanonical</code> flag). The 
  !! <code>trans</code> array is updated within this function if it is in use (if <code>update_trans=.true.</code>)
  !! and we have reached equilibration.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> E_trial </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   The trial energy.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> macro_trial </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   The trial macrostate number.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> V_trial </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   The trial volume.
  !!   </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
  function metropolis_prob_lattice(E_trial,macro_trial,V_trial)
      real(rk), intent(in) :: E_trial
      integer(ik), intent(in) :: macro_trial
      real(rk), intent(in) :: V_trial
      real(rk) :: metropolis_prob_lattice
      real(rk) :: prob_B
      if(enable_multicanonical) then
          metropolis_prob_lattice= &
              metropolis_prob_MC_vol_unscaled_pos(beta,E,E_trial,eta_grid(macro),eta_grid(macro_trial),n_part,P,V,V_trial)
          if(update_trans .and. sweeps>=sweep_equil_reference+equil_sweeps) then
              ! Update 'trans' using the Boltzmann probability of the transition
              prob_B=metropolis_prob_B_vol_unscaled_pos(beta,E,E_trial,n_part,P,V,V_trial)
              call trans_update(macro,macro_trial,prob_B)
          end if
      else
          metropolis_prob_lattice=metropolis_prob_B_vol_unscaled_pos(beta,E,E_trial,n_part,P,V,V_trial)
          if(update_trans .and. sweeps>=sweep_equil_reference+equil_sweeps) then
              call trans_update(macro,macro_trial,metropolis_prob_lattice)
          end if
      end if
  end function metropolis_prob_lattice




  !! <h4> <code> function amend_OOB_variables(M,key) </code> </h4>
  !! <p>
  !! This procedure updates <code>rejected_moves_M_OOB</code>, <code>M_OOB_high</code> and <code>M_OOB_low</code>
  !! to reflect the fact that M is outwith the supported range. <code>rejected_moves_M_OOB</code> is amended by 1,
  !! and <code>M_OOB_high</code> and <code>M_OOB_low</code> is amended to be <code>M</code> if <code>M</code> is
  !! greater than or less than <code>M_OOB_high</code> or <code>M_OOB_low</code> respectively.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> M </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   The out of bounds order parameter.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> key </code> </td>
  !!   <td> <code> integer(ik), intent(in) </code> </td>
  !!   <td>
  !!   Label specifying how the order parameter is out of bounds: -2 if it is below the supported range;
  !!   -1 if it is above the supported range.
  !!   </td>
  !!  </tr>
  !! </table>
  subroutine amend_OOB_variables(M,key)
    real(rk), intent(in) :: M
    integer(ik), intent(in) :: key

    select case(key)
    case(-2)
        ! M is below supported range
        if(M<M_OOB_low) then
            M_OOB_low=M
        end if
    case(-1)
        ! M is above the supported range
        if(M>M_OOB_high) then
            M_OOB_high=M
        end if
    case default
        write(0,*) "monteswitch_mod: Error. 'key' in 'amend_OOB_variables' does not equal -1 or -2."
        stop 1
    end select

    rejected_moves_M_OOB=rejected_moves_M_OOB+1

  end subroutine amend_OOB_variables



  !! <h4> <code> function is_macro_within_barriers(macro) </code> </h4>
  !! <p>
  !! This function returns a logical value which specifies whether the macrostate in the argument
  !! is within the order parameters. This function updates <code>rejected_moves_M_barrier</code> if it is not.
  !! This function tacitly assumes that the macrostate in the argument is within the supported range.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> macro </code> </td>
  !!   <td> <code> integer(ik), intent(in) </code> </td>
  !!   <td>
  !!   The macrostate to check.
  !!   </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> logical </code> </p>
  function is_macro_within_barriers(macro)
    integer(ik), intent(in) :: macro
    logical :: is_macro_within_barriers
    is_macro_within_barriers=.true.
    if(macro<barrier_macro_low .or. macro>barrier_macro_high) then 
       is_macro_within_barriers=.false.
       rejected_moves_M_barrier=rejected_moves_M_barrier+1
    end if
  end function is_macro_within_barriers




  !! <h4> <code> subroutine after_move() </code> </h4>
  !! <p>
  !! This subroutine performs general tasks which should be performed after every move:
  !! updating <code>moves</code>; updating <code>M_counts_1</code> and <code>M_counts_2</code> 
  !! after equilibration; updating the order parameter barriers if they are in use;
  !! updating sums used in evaluating equilibrium properties
  !! </p>
  subroutine after_move()
    real(rk) :: weight ! Used in updating the sums pertaining to the equilibrium properties
    integer(ik) :: n ! Particle number
    real(rk) :: usqrd ! The magnitude of u squared for a given particle

    ! Update interactions stuff
    call after_all_interactions(Lx(1), Ly(1), Lz(1), spec_1, pos_1, R_1, &
        Lx(2), Ly(2), Lz(2), spec_2, pos_2, R_2, u)

    ! Update the barrier stuff if need be
    if(enable_barriers) then
       call update_barriers()
    end if

    ! Update the order parameter histogram after equilibration
    if(sweeps>=sweep_equil_reference+equil_sweeps) then
       select case(lattice)
       case(1)
          call update_histogram(M_grid_size,M_grid,M_counts_1,(M_grid(2)-M_grid(1)),M)
       case(2)
          call update_histogram(M_grid_size,M_grid,M_counts_2,(M_grid(2)-M_grid(1)),M)
       case default
          write(0,*) "monteswitch_mod: Error. 'lattice' does not equal 1 or 2."
          stop 1
       end select
    end if

    ! Update sums used in evaluating equilibrium properties after equilibration: details of what's going
    ! on are provided in the 'run' subroutine.
    if(sweeps>=sweep_equil_reference+equil_sweeps .and. calc_equil_properties) then
       ! Set the weight which will be given to this state
       if(enable_multicanonical) then
          ! For multicanonical sampling
          weight=exp(-eta)
       else
          ! For Boltzmann sampling
          weight=1.0_rk
       end if
       
       ! Update the sums
       select case(lattice)
       case(1)
          intrablock_counts_1=intrablock_counts_1+weight
          ! Update the energy/enthalpy sum
          if(enable_vol_moves) then
             ! Enthalpy for NPT ensemble
             intrablock_sum_H_1=intrablock_sum_H_1+(E+P*V)*weight
          else
             ! Energy for NVT ensemble
             intrablock_sum_H_1=intrablock_sum_H_1+E*weight
          end if
          ! Update the volume sum
          intrablock_sum_V_1=intrablock_sum_V_1+V*weight
          ! Update the umsd sum
          do n=1,n_part
             usqrd=dot_product(u(:,n),u(:,n))
             intrablock_sum_umsd_1(n)=intrablock_sum_umsd_1(n)+weight*usqrd
          end do
          ! Update the L sum
          intrablock_sum_L_1(1) = intrablock_sum_L_1(1) + weight*Lx(1)
          intrablock_sum_L_1(2) = intrablock_sum_L_1(2) + weight*Ly(1)
          intrablock_sum_L_1(3) = intrablock_sum_L_1(3) + weight*Lz(1)

       case(2)
          intrablock_counts_2=intrablock_counts_2+weight
          ! Update the energy/enthalpy sum
          if(enable_vol_moves) then
             ! Enthalpy for NPT ensemble
             intrablock_sum_H_2=intrablock_sum_H_2+(E+P*V)*weight
          else
             ! Energy for NVT ensemble
             intrablock_sum_H_2=intrablock_sum_H_2+E*weight
          end if
          ! Update the volume sum
          intrablock_sum_V_2=intrablock_sum_V_2+V*weight
          ! Update the umsd sum
          do n=1,n_part
             usqrd=dot_product(u(:,n),u(:,n))
             intrablock_sum_umsd_2(n)=intrablock_sum_umsd_2(n)+weight*usqrd
          end do 
          ! Update the L sum
          intrablock_sum_L_2(1) = intrablock_sum_L_2(1) + weight*Lx(2)
          intrablock_sum_L_2(2) = intrablock_sum_L_2(2) + weight*Ly(2)
          intrablock_sum_L_2(3) = intrablock_sum_L_2(3) + weight*Lz(2)

       case default
          write(0,*) "monteswitch_mod: Error. 'lattice' does not equal 1 or 2."
          stop 1
       end select
    end if

    moves=moves+1
    
  contains
    
    ! This nested subroutine updates the order parameter barriers. 
    ! It assumes that the current microstate is between the order parameter barriers.
    subroutine update_barriers()      
      ! Ammend 'moves_since_lock'
      moves_since_lock=moves_since_lock+1
      
      if(moves_since_lock>lock_moves) then
         ! If we are waiting to lock into the new microstate, i.e. if moves_since_lock>lock_moves, then 
         ! check whether the order parameter is within the newly available macrostate. If so then lock the 
         ! system into the new macrostate.
         select case(barrier_to_lock)
         case(1)
            ! If 'barrier_to_lock'=1 then the lower barrier needs to be shifted upwards once the system
            ! enters the newly available macrostate 'barrier_macro_high' 
            if(macro==barrier_macro_high) then
               barrier_macro_low=barrier_macro_low+1
               moves_since_lock=0
            end if
         case(2)
            ! If 'barrier_to_lock'=2 then the higher barrier needs to be shifted downwards once the 
            ! system enters the newly available macrostate 'barrier_macro_low'
            if(macro==barrier_macro_low) then
               barrier_macro_high=barrier_macro_high-1
               moves_since_lock=0
            end if
         case default
            write(0,*) "monteswitch_mod: Error. 'barrier_to_lock' is not 1 or 2."
            stop 1
         end select
         
      else if(moves_since_lock==lock_moves) then
         ! Move the barriers to open up a new macrostate if moves_since_lock=lock_moves. Note that if 
         ! the upper barrier is as high as is possible then we automatically shift the lower barrier; and
         ! if the lower barrier is as low as is possible then we automatically shift the higher barrier.
         if(barrier_macro_low==1) then
            barrier_macro_high=barrier_macro_high+1
            ! We have just shifted the high barrier upwards. Set barrier_to_lock to 1 to reflect the
            ! fact that the low barrier is now to be shifted upwards once we move into the newly
            ! available macrostate
            barrier_to_lock=1
         else if(barrier_macro_high==M_grid_size) then
            barrier_macro_low=barrier_macro_low-1
            ! We have shifted the low barrier downwards. Set barrier_to_lock to 2 to reflect the
            ! fact that the high barrier is now to be shifted downwards once we move into the newly
            ! available macrostate
            barrier_to_lock=2
         else
            !
            ! If we are not in an 'edge' macrostate...
            !
            select case(barrier_dynamics)
            case(barrier_dynamics_random)
               ! Decide to shift either the high or low barrier at random
               if(get_random_number()>0.5) then
                  barrier_macro_high=barrier_macro_high+1
                  barrier_to_lock=1
               else
                  barrier_macro_low=barrier_macro_low-1
                  barrier_to_lock=2
               end if
            case(barrier_dynamics_pong_up, barrier_dynamics_pong_down)
               ! Code for 'pong' barrier_dynamics
               !
               ! If barrier_to_lock=1 then we previously locked 'upwards' and we will open a new macrostate
               ! upwards (and continue to do so until we reach the highest allowed macrostate). The opposite
               ! applies for barrier_to_lock=2
               !
               ! Note that the code below only comes into effect if we are not in the highest or lowest supported
               ! macrostate, in which case barrier_to_lock is set such that we will move away from the edge
               ! of the supported macrostate range.
               !
               select case(barrier_to_lock)
               case(1)
                  barrier_macro_high=barrier_macro_high+1
               case(2)
                  barrier_macro_low=barrier_macro_low-1
               case default
                  write(0,*) "monteswitch_mod: Error. 'barrier_to_lock' is not 1 or 2."
                  stop 1
               end select
            case default
               write(0,*) "monteswitch_mod: Error. 'barrier_dynamics' value is not recognised."
               stop 1
            end select
         end if
      end if
    end subroutine update_barriers

  end subroutine after_move




  !! <h3> <code> subroutine translate_positions(r,Lx,Ly,Lz) </code> </h3>
  !! <p>
  !! <code> translate_positions </code> translates the position vectors in the 
  !! array <code>r</code> so that they become their analogous positions in a cuboid 
  !! whose faces are x=0, x=<code>Lx</code>, y=0, y=<code>Ly</code>, z=0, and
  !! z=<code>Lz</code>. The format of <code>r</code> should be as follows: 
  !! (r(1,n),r(2,n),r(3,n)) is the nth position to be translated.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> r </code> </td>
  !!   <td> <code> real(rk), dimension(:,:), intent(inout) </code> </td>
  !!   <td> Positions to be translated. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> Lx </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Length of cuboid along the x-axis. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> Ly </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Length of cuboid along the y-axis. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> Lz </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Length of cuboid along the z-axis. </td>
  !!  </tr>
  !! </table>
  subroutine translate_positions(r,Lx,Ly,Lz)
    real(rk), dimension(:,:), intent(inout) :: r
    real(rk), intent(in) :: Lx, Ly, Lz
    integer(ik) :: i

    do i=1,n_part
       call translate_position(r(:,i),Lx,Ly,Lz)
    end do    
  end subroutine translate_positions




  !! <h3> <code> subroutine translate_position(r,Lx,Ly,Lz) </code> </h3>
  !! <p>
  !! <code> translate_positios </code> translates the position vector in the 
  !! array <code>r</code> so that it becomes the analogous positions in a cuboid 
  !! whose faces are x=0, x=<code>Lx</code>, y=0, y=<code>Ly</code>, z=0, and
  !! z=<code>Lz</code>.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> r </code> </td>
  !!   <td> <code> real(rk), dimension(3), intent(inout) </code> </td>
  !!   <td> Position to be translated. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> Lx </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Length of cuboid along the x-axis. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> Ly </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Length of cuboid along the y-axis. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> Lz </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Length of cuboid along the z-axis. </td>
  !!  </tr>
  !! </table>
  subroutine translate_position(r,Lx,Ly,Lz)
    real(rk), dimension(3), intent(inout) :: r
    real(rk), intent(in) :: Lx, Ly, Lz
    integer(ik) :: i

    r(1)=modulo(r(1),Lx)
    r(2)=modulo(r(2),Ly)
    r(3)=modulo(r(3),Lz)

  end subroutine translate_position




  !! <h3> <code>  subroutine update_histogram(size,bin_min,counts,bin_width,x) </code> </h3>
  !! <p>
  !! <code> update_histogram </code> updates a histogram defined by the arrays
  !! <code>bin_min</code> and <code>counts</code> - which each have <code>size</code> 
  !! elements - by adding 1 to the bin corresponding to the value <code>x</code>.
  !! The appropriate bin for <code>x</code> is determined as follows: if <code>x</code>
  !! is between <code>bin_min(n)</code> (including <code>bin_min(n)</code>) and
  !! <code>bin_min(n)+bin_width</code> (not including <code>bin_min(n)+bin_width</code>),
  !! then <code>x</code> is in bin n. Note that the valid range for <code>x</code> is
  !! <code>bin_min(1)</code> (including <code>bin_min(1)</code>) to
  !! <code>bin_min(nbins)+bin_width</code> (not including 
  !! <code>bin_min(nbins)+bin_width</code>). If <code>x</code> is not within this range,
  !! then the histogram is not updated. Note that the width of all bins must be the same;
  !! and all arrays must have the default range of indices, i.e. 1 to the dimension of
  !! the array.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> size </code> </td>
  !!   <td> <code> integer(ik), intent(in) </code> </td>
  !!   <td> The number of bins. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> bin_min </code> </td>
  !!   <td> <code> real(rk), dimension(size), intent(in)  </code> </td>
  !!   <td> Array containing the lower bounds of each bin. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> counts </code> </td>
  !!   <td> <code> integer(ik), dimension(size), intent(inout) </code> </td>
  !!   <td> Array containing the number of counts for each bin. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> bin_width </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> The width of each bin. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> x </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> The value whose corresponding bin is to be updated. </td>
  !!  </tr>
  !! </table>
  subroutine update_histogram(size,bin_min,counts,bin_width,x)
    integer(ik), intent(in) :: size
    real(rk), dimension(size), intent(in) :: bin_min
    integer(ik), dimension(size), intent(inout) :: counts
    real(rk), intent(in) :: bin_width
    real(rk), intent(in) :: x
    integer(ik) :: bin
    ! Determine which bin number 'x' corresponds to 
    bin= floor( (x-bin_min(1))/bin_width +1 )
    ! If 1<=bin<=size then 'x' is within the histogram's range and we update it; otherwise
    ! do nothing
    if(bin>=1 .and. bin<=size) then
       counts(bin)=counts(bin)+1
    end if
  end subroutine update_histogram




  !! <h3> Monte Carlo moves </h3>



  !! <h4> <code> subroutine move_lattice() </code> </h4>
  !! <p>
  !! This subroutine 'moves' the lattice and updates the relevant counters. Note that the current order parameter and
  !! macrostate are not affected by a lattice move.
  !! </p>
  subroutine move_lattice()
    ! The probability the move should be accepted
    real(rk) :: prob

    ! Calculate the probability of accepting a lattice switch
    select case(lattice)
    case(1)
       prob=metropolis_prob_lattice(E_2, macro, Lx(2)*Ly(2)*Lz(2) )
    case(2)
       prob=metropolis_prob_lattice(E_1, macro, Lx(1)*Ly(1)*Lz(1) )
    case default
       write(0,*) "monteswitch_mod: Error. 'lattice' is not 1 or 2."
       stop 1
    end select
    ! Update the lattice or not
    if(get_random_number()<prob) then
       select case(lattice)
       case(1)
          E=E_2
          V=Lx(2)*Ly(2)*Lz(2)
          lattice=2
       case(2)
          E=E_1
          V=Lx(1)*Ly(1)*Lz(1)
          lattice=1
       case default
          write(0,*) "monteswitch_mod: Error. 'lattice' is not 1 or 2."
          stop 1
       end select
       ! Ammend accepted counters
       accepted_moves_lattice=accepted_moves_lattice+1
       ! Perform 'interactions'-specific tasks after an accepted move
       call after_accepted_lattice_interactions(Lx(1), Ly(1), Lz(1), spec_1, pos_1, R_1, &
           Lx(2), Ly(2), Lz(2), spec_2, pos_2, R_2, u)
    end if
    ! Ammend counter
    moves_lattice=moves_lattice+1
    ! Perform other tasks
    call after_move()
  end subroutine move_lattice




  !! <h4> <code> subroutine move_particle(i) </code> </h4>
  !! <p>
  !! This subroutine moves particle 'i'.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> i </code> </td>
  !!   <td> <code>  integer(ik), intent(in) </code> </td>
  !!   <td>
  !!   The particle number to move.
  !!   </td>
  !!  </tr>
  !! </table>
  subroutine move_particle(i)
    integer(ik), intent(in) :: i
    ! 'u_trial' is the trial displacement vector for particle i
    real(rk), dimension(3) :: u_trial
    ! Trial position vectors for each lattice
    real(rk), dimension(3) :: pos_new_1
    real(rk), dimension(3) :: pos_new_2
    ! 'delta_E_1' and 'delta_E_2' are the energy differences between the trial move and the current energy
    ! lattice types 1 and 2.
    real(rk) :: delta_E_1
    real(rk) :: delta_E_2
    ! 'delta_E_trial' is the energy difference given the actual lattice type 'lattice'
    real(rk) :: delta_E
    ! The trial order parameter
    real(rk) :: M_trial
    ! The trial macrostate
    integer(ik) :: macro_trial
    ! The trial value of the weight function
    real(rk) :: eta_trial
    ! The probability the move should be accepted
    real(rk) :: prob
    ! This is a boolean which determines whether or not the move is within the required range
    ! of order parameters
    logical :: proceed
    ! The change in the particle displacement - which constitutes the move
    real(rk), dimension(3) :: delta_u
    ! For use in loops
    integer(ik) :: j

    ! Generate trial positions and displacement for particle 'i' according to a random walk
    delta_u(1)=top_hat_rand(part_step)
    delta_u(2)=top_hat_rand(part_step)
    delta_u(3)=top_hat_rand(part_step)
    u_trial=u(:,i)+delta_u
    pos_new_1=pos_1(:,i)+delta_u
    call translate_position(pos_new_1,Lx(1),Ly(1),Lz(1))
    pos_new_2=pos_2(:,i)+delta_u
    call translate_position(pos_new_2,Lx(2),Ly(2),Lz(2))

    ! Calculate the change in energies for each lattice, using the new positions for i
    delta_E_1=calc_energy_part_move(1,Lx(1),Ly(1),Lz(1),spec_1,pos_1,pos_new_1,R_1,u,u_trial,i)
    delta_E_2=calc_energy_part_move(2,Lx(2),Ly(2),Lz(2),spec_2,pos_2,pos_new_2,R_2,u,u_trial,i)

    ! Set delta_E depending on the lattice, and calculate M for the trial state
    select case(lattice)
    case(1)
       delta_E=delta_E_1
    case(2)
       delta_E=delta_E_2
    case default
       write(0,*) "monteswitch_mod: Error. 'lattice' is not 1 or 2."
       stop 1
    end select

    ! Calculate the trial order parameter and macrostate number
    M_trial=calc_M(E_1+delta_E_1, E_2+delta_E_2, Lx(1)*Ly(1)*Lz(1), Lx(2)*Ly(2)*Lz(2))
    macro_trial=get_macro(M_trial)
    ! get_macro(M_trial) is -2 or -1 if M_trial is outwith the supported range. In this case the
    ! move is rejected, and the OOB variables are updated.
    proceed=.true.
    if(macro_trial<1) then
        proceed=.false.
        call amend_OOB_variables(M_trial,macro_trial)
    end if

    if(proceed) then
       
       ! The trial move is within the supported order parameter range, and macro_trial is sensible...

       eta_trial=eta_grid(macro_trial)

       ! Calculate the probability of accepting the move. Note that this procedure updates
       ! 'trans' if it is in use. Hence we necessarily reject the move if macro_trial is outwith the order
       ! parameter barriers AFTER this procedure.
       prob=metropolis_prob(E+delta_E,macro_trial)
       ! Check if the move is within the barriers if the barriers are in use. If not then reject it on those grounds,
       ! i.e. set proceed to .false.
       if(enable_barriers) then
          proceed=is_macro_within_barriers(macro_trial)
       end if
       ! Finally, if proceed=.true. accept the move or not based on a random number
       if(proceed .and. get_random_number()<prob) then

           ! Accept the change to pos_1, pos_2 and u for particle i
           pos_1(:,i)=pos_new_1
           pos_2(:,i)=pos_new_2
           u(:,i)=u_trial
           
           ! Update the energies, M, macro and eta
           E_1=E_1+delta_E_1
           E_2=E_2+delta_E_2
           E=E+delta_E
           M=M_trial
           macro=macro_trial
           eta=eta_trial
           
           ! If we use the centre-of-mass frame...
           if(enable_COM_frame) then
               ! Ammend the particle displacements: shift all displacements by -delta_u/n_part
               ! to counter the fact that u(i,:) has shifted by delta_u
               u(1,:)=u(1,:)-delta_u(1)/n_part
               u(2,:)=u(2,:)-delta_u(2)/n_part
               u(3,:)=u(3,:)-delta_u(3)/n_part
               ! Update the positions to correspond to this shift - otherwise, while the separations between the
               ! particles according to the positions will be correct, the absolute positions of the particles will
               ! not be 'the lattice site positions plus the displacement'.
               ! For optimisation purposes, if one does not care about the positions not matching the displacements,
               ! and considers the arrays pos_1 and pos_2 simply as 'work arrays' to enable energies to be efficiently
               ! evaluated, then perhaps this can be removed?
               pos_1(1,:)=pos_1(1,:)-delta_u(1)/n_part
               pos_1(2,:)=pos_1(2,:)-delta_u(2)/n_part
               pos_1(3,:)=pos_1(3,:)-delta_u(3)/n_part
               call translate_positions(pos_1,Lx(1),Ly(1),Lz(1))
               pos_2(1,:)=pos_2(1,:)-delta_u(1)/n_part
               pos_2(2,:)=pos_2(2,:)-delta_u(2)/n_part
               pos_2(3,:)=pos_2(3,:)-delta_u(3)/n_part
               call translate_positions(pos_2,Lx(2),Ly(2),Lz(2))
           end if
           ! Ammend accepted counters
           accepted_moves_part=accepted_moves_part+1
           ! Perform 'interactions'-specific tasks after an accepted move
           call after_accepted_part_interactions(i, Lx(1), Ly(1), Lz(1), spec_1, pos_1, R_1, &
               Lx(2), Ly(2), Lz(2), spec_2, pos_2, R_2, u)

       end if
    end if
    ! Ammend counter
    moves_part=moves_part+1
    ! Perform other tasks
    call after_move()

  end subroutine move_particle




  !! <h4> <code> subroutine move_volume() </code> </h4>
  !! <p>
  !! This subroutine moves the volume. The flag <code>vol_dynamics</code> determines the
  !! nature of the volume moves.
  !! </p>
  subroutine move_volume()
    ! 'Lx_trial', 'Ly_trial' and 'Lz_trial' are the trial dimensions of the supercell for both lattice types
    real(rk), dimension(2) :: Lx_trial
    real(rk), dimension(2) :: Ly_trial
    real(rk), dimension(2) :: Lz_trial
    ! 'u_trial' are the trial displacements
    real(rk), dimension(3,n_part) :: u_trial
    ! 'R_1_trial' and 'R_2_trial' are the trial lattice vectors for lattice types 1 and 2 respectively
    real(rk), dimension(3,n_part) :: R_1_trial
    real(rk), dimension(3,n_part) :: R_2_trial
    ! Trial positions for both lattices
    real(rk), dimension(3,n_part) :: pos_trial_1
    real(rk), dimension(3,n_part) :: pos_trial_2
    ! 'E_1_trial' and 'E_2_trial' are the trial energies for lattice types 1 and 2.
    real(rk) :: E_1_trial
    real(rk) :: E_2_trial
    ! 'E_trial; is the trial energy, given the actual lattice type 'lattice'
    real(rk) :: E_trial
    ! The trial order parameter
    real(rk) :: M_trial
    ! The trial macrostate
    integer(ik) :: macro_trial
    ! The trial value of the weight function
    real(rk) :: eta_trial
    ! The trial volume
    real(rk) :: V_trial
    ! The probability the move should be accepted
    real(rk) :: prob
    ! This is a boolean which determines whether or not the move is within the required range
    ! of order parameters
    logical :: proceed

    ! Generate trial varaiables 'Lx_trial', 'Ly_trial', 'Lz_trial', 'V_trial',
    ! 'R_1_trial', 'R_2_trial', 'u_trial', 'pos_trial_1' and 'pos_trial_2'
    select case(vol_dynamics)
    case(vol_dynamics_FVM)
       call vol_trial_FVM()
    case(vol_dynamics_UVM)
       call vol_trial_UVM()
    case default
       write(0,*) "monteswitch_mod: Error. 'vol_dynamics' value is not recognised."
       stop 1
    end select

    ! Calculate the other trial variables
    E_1_trial=calc_energy_scratch(1,Lx_trial(1),Ly_trial(1),Lz_trial(1),spec_1,pos_trial_1,R_1_trial,u_trial)
    E_2_trial=calc_energy_scratch(2,Lx_trial(2),Ly_trial(2),Lz_trial(2),spec_2,pos_trial_2,R_2_trial,u_trial)

    select case(lattice)
    case(1)
       E_trial=E_1_trial
    case(2)
       E_trial=E_2_trial
    case default
       write(0,*) "monteswitch_mod: Error. 'lattice' is not 1 or 2."
       stop 1
    end select

    ! Set M_trial. Don't set macro_trial or eta_trial until after we verify that M_trial
    ! is within the allowed range later
    M_trial=calc_M(E_1_trial, E_2_trial, Lx_trial(1)*Ly_trial(1)*Lz_trial(1), Lx_trial(2)*Ly_trial(2)*Lz_trial(2))
    macro_trial=get_macro(M_trial)
    ! get_macro(M_trial) is -2 or -1 if M_trial is outwith the supported range. In this case the
    ! move is rejected, and the OOB variables are updated.
    proceed=.true.
    if(macro_trial<1) then
        proceed=.false.
        call amend_OOB_variables(M_trial,macro_trial)
    end if

    if(proceed) then

       ! The trial move is within the supported order parameter range, and macro_trial is sensible...

       eta_trial=eta_grid(macro_trial)

       ! Calculate the probability of accepting the move. Note that this procedure updates
       ! 'trans' if it is in use. Hence we necessarily reject the move if macro_trial is outwith the order
       ! parameter barriers AFTER this procedure.
       prob=metropolis_prob_vol(E_trial,macro_trial,V_trial)
       ! Check if the move is within the barriers if the barriers are in use. If not then reject it on those grounds,
       ! i.e. set proceed to .false.
       if(enable_barriers) then
          proceed=is_macro_within_barriers(macro_trial)
       end if
       ! Finally, if proceed=.true. accept the move or not based on a random number
       if(proceed .and. get_random_number()<prob) then
          Lx=Lx_trial
          Ly=Ly_trial
          Lz=Lz_trial
          V=V_trial
          u=u_trial
          R_1=R_1_trial
          R_2=R_2_trial
          pos_1=pos_trial_1
          pos_2=pos_trial_2
          E_1=E_1_trial
          E_2=E_2_trial
          E=E_trial
          M=M_trial
          macro=macro_trial
          eta=eta_trial
          ! Ammend accepted counters
          accepted_moves_vol=accepted_moves_vol+1
          ! Perform 'interactions'-specific tasks after an accepted move
          call after_accepted_vol_interactions(Lx(1), Ly(1), Lz(1), spec_1, pos_1, R_1, &
              Lx(2), Ly(2), Lz(2), spec_2, pos_2, R_2, u)
       end if
    end if
    ! Ammend counters
    moves_vol=moves_vol+1
    ! Perform other tasks
    call after_move()

  contains

    ! This nested subroutine sets the variables 'Lx_trial', 'Ly_trial', 'Lz_trial', 'V_trial',
    ! 'R_1_trial', 'R_2_trial' and 'u_trial' to reflect the following: ln(V) is increased/decreased
    ! by a random amount uniformly drawn from -vol_step to vol_step, where an equal scaling factor is
    ! applied to the x, y and z dimensions.
    ! (FVM=fixed aspect ratio volume moves)
    subroutine vol_trial_FVM()
      ! Scaling factor for the volume
      real(rk) :: S_vol
      ! Scaling factor for each dimension for each lattice
      real(rk) :: S
      ! Determine the scaling factor for the volume
      S_vol=exp(top_hat_rand(vol_step))
      ! Calculate the scaling factor for each dimension required to achieve this new volume
      S=S_vol**(1.0_rk/3.0_rk)

      ! Set the trial variables accordingly.
      V_trial=V*S_vol
      Lx_trial(1)=Lx(1)*S
      Ly_trial(1)=Ly(1)*S
      Lz_trial(1)=Lz(1)*S
      Lx_trial(2)=Lx(2)*S
      Ly_trial(2)=Ly(2)*S
      Lz_trial(2)=Lz(2)*S
      R_1_trial=R_1*S
      R_2_trial=R_2*S
      u_trial=u*S
      pos_trial_1=pos_1*S
      pos_trial_2=pos_2*S

    end subroutine vol_trial_FVM

    ! This nested subroutine sets the variables 'Lx_trial', 'Ly_trial', 'Lz_trial', 'V_trial',
    ! 'R_1_trial', 'R_2_trial' and 'u_trial' to reflect the following: a random dimension (x,y or z)
    ! is selected, and the supercell is stretched/compressed in that dimension such that
    ! ln(V) is increased/decreased by a random amount uniformly drawn from -vol_step to vol_step.
    ! (UVM=unconstrained aspect ratio volume moves)
    subroutine vol_trial_UVM()
      real(rk) :: rand
      ! Scaling factor for the volume
      real(rk) :: S

      ! Set the trial state vectors to the current state initially
      Lx_trial=Lx
      Ly_trial=Ly
      Lz_trial=Lz
      R_1_trial=R_1
      R_2_trial=R_2
      pos_trial_1=pos_1
      pos_trial_2=pos_2
      u_trial=u

      ! Determine the scaling factor for the volume
      S=exp(top_hat_rand(vol_step))

      ! Select a random dimension: x, y or z, and multiply that dimension by S
      rand=get_random_number()
      if(rand<1.0_rk/3.0_rk) then
         ! Move x
         V_trial=V*S
         Lx_trial(1)=Lx(1)*S
         Lx_trial(2)=Lx(2)*S
         R_1_trial(1,:)=R_1(1,:)*S
         R_2_trial(1,:)=R_2(1,:)*S
         u_trial(1,:)=u(1,:)*S
         pos_trial_1(1,:)=pos_1(1,:)*S
         pos_trial_2(1,:)=pos_2(1,:)*S

      else if(rand<2.0_rk/3.0_rk) then
         ! Move y
         V_trial=V*S
         Ly_trial(1)=Ly(1)*S
         Ly_trial(2)=Ly(2)*S
         R_1_trial(2,:)=R_1(2,:)*S
         R_2_trial(2,:)=R_2(2,:)*S
         u_trial(2,:)=u(2,:)*S
         pos_trial_1(2,:)=pos_1(2,:)*S
         pos_trial_2(2,:)=pos_2(2,:)*S

      else
         ! Move z
         V_trial=V*S
         Lz_trial(1)=Lz(1)*S
         Lz_trial(2)=Lz(2)*S
         R_1_trial(3,:)=R_1(3,:)*S
         R_2_trial(3,:)=R_2(3,:)*S
         u_trial(3,:)=u(3,:)*S
         pos_trial_1(3,:)=pos_1(3,:)*S
         pos_trial_2(3,:)=pos_2(3,:)*S

      end if

    end subroutine vol_trial_UVM

  end subroutine move_volume




  !! <h3> Procedures to simplify initialisation </h3>



  !! <h4> <code> subroutine initialise_M_variables(M_grid_size_in,M_grid_max) </code> </h4>
  !! <p>
  !! This procedure initialises <code>M_grid</code> such that the minimum order parameter considered is
  !! <code>M_grid_min</code>, the maximum order parameter considered is <code>M_grid_max</code>, and the
  !! number of macrostates is <code>M_grid_size_in</code>; <code>M_OOB_low</code> and <code>M_OOB_high</code> 
  !! to be <code>M_grid_min</code> and <code>M_grid_max</code> respectively; <code>trans</code> to be of 
  !! size <code>M_grid_size_in</code> along both dimensions, with all elements set to 0; <code>eta_grid</code>
  !! to be of size <code>M_grid_size_in</code> with all elements set to 0. Note that if <code>M_grid</code>, 
  !! <code>trans</code> or <code>eta_grid</code> is already allocated then this
  !! subroutine deallocates it before initialising it according to the arguments. Furthermore, this procedure
  !! performs some safety checks on <code>M_grid_size</code>, <code>M_grid_min</code> and 
  !! <code>M_grid_max</code>, and stops the program with an error if they are not sensible (in particular if
  !! <code>M_grid_size</code> < 2).
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> M_grid_size_in </code> </td>
  !!   <td> <code> integer(ik), intent(in) </code> </td>
  !!   <td> The number of macrostates which will be defined. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> M_grid_in </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> The minimum order parameter which will be considered. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> M_grid_max </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> The maximum order parameter which will be considered. </td>
  !!  </tr>
  !! </table>
  subroutine initialise_M_variables(M_grid_size_in,M_grid_min,M_grid_max)
    integer(ik), intent(in) :: M_grid_size_in
    real(rk), intent(in) :: M_grid_min
    real(rk), intent(in) :: M_grid_max
    integer(ik) :: i

    if(M_grid_min > M_grid_max) then
        write(0,*) "monteswitch_mod: Error. 'M_grid_min' > 'M_grid_max'."
        stop 1
    end if

    M_grid_size=M_grid_size_in
    if(M_grid_size < 2) then
        write(0,*) "monteswitch_mod: Error. 'M_grid_size' must be > 1."
        stop 1
    end if

    if(allocated(M_grid)) then
       deallocate(M_grid)
    end if
    allocate(M_grid(M_grid_size))
    do i=1,M_grid_size
       M_grid(i)=M_grid_min+(i-1)*(M_grid_max-M_grid_min)/M_grid_size
    end do

    M_OOB_low=M_grid_min
    M_OOB_high=M_grid_max

    if(allocated(trans)) then
       deallocate(trans)
    end if
    allocate(trans(4,M_grid_size))
    trans=0.0_rk

    if(allocated(eta_grid)) then
       deallocate(eta_grid)
    end if
    allocate(eta_grid(M_grid_size))
    eta_grid=0.0_rk

  end subroutine initialise_M_variables




  !! <h4> <code> subroutine initialise_counters() </code> </h4>
  !! <p>
  !! This procedure initialises all counters to be 0, except <code>trans</code> and anything to do with order parameter 
  !! barriers. 'Counters' include those all sums used to evaluate equilibrium properties, sweep numbers, etc.
  !! It does include the equilibrium values and their uncertainties, e.g. <code>equil_DeltaF</code>, <code>sigma_equil_DeltaF</code>; 
  !! however note that these values are meaningless until at least one block has been considered for these variables.
  !! It does include <code>sweep_equil_reference</code>, which is set to 0.
  !! If <code>M_counts_1</code> or <code>M_counts_2</code> is already allocated then this subroutine deallocates them and
  !! then allocates them to have a size corresponding to the current value of <code>M_grid_size</code>.
  !! </p>
  !! <p> <b> Dependencies: </b> <code>M_grid_size</code>and <code>n_part</code> must be set before this procedure is called. </p>
  subroutine initialise_counters()
    if(allocated(M_counts_1)) then
       deallocate(M_counts_1)
    end if
    if(allocated(M_counts_2)) then
       deallocate(M_counts_2)
    end if
    allocate(M_counts_1(M_grid_size))
    allocate(M_counts_2(M_grid_size))
    M_counts_1=0
    M_counts_2=0
    sweeps=0
    moves=0
    moves_lattice=0
    accepted_moves_lattice=0
    moves_part=0
    accepted_moves_part=0
    moves_vol=0
    accepted_moves_vol=0
    rejected_moves_M_OOB=0
    melts=0
    block_counts=0
    intrablock_counts_1=0.0_rk
    intrablock_counts_2=0.0_rk
    interblock_sum_DeltaF=0.0_rk
    interblock_sum_DeltaF_sqrd=0.0_rk
    block_counts_DeltaF=0 
    equil_DeltaF=0.0_rk
    sigma_equil_DeltaF=0.0_rk
    intrablock_sum_H_1=0.0_rk
    intrablock_sum_H_2=0.0_rk 
    interblock_sum_H_1=0.0_rk
    interblock_sum_H_2=0.0_rk
    interblock_sum_H_1_sqrd=0.0_rk
    interblock_sum_H_2_sqrd=0.0_rk
    block_counts_H_1=0 
    block_counts_H_2=0 
    equil_H_1=0.0_rk
    sigma_equil_H_1=0.0_rk
    equil_H_2=0.0_rk
    sigma_equil_H_2=0.0_rk
    intrablock_sum_V_1=0.0_rk
    intrablock_sum_V_2=0.0_rk
    interblock_sum_V_1=0.0_rk
    interblock_sum_V_2=0.0_rk
    interblock_sum_V_1_sqrd=0.0_rk
    interblock_sum_V_2_sqrd=0.0_rk
    block_counts_V_1=0 
    block_counts_V_2=0
    equil_V_1=0.0_rk
    sigma_equil_V_1=0.0_rk
    equil_V_2=0.0_rk
    sigma_equil_V_2=0.0_rk
    intrablock_sum_L_1=0.0_rk
    intrablock_sum_L_2=0.0_rk
    interblock_sum_L_1=0.0_rk
    interblock_sum_L_2=0.0_rk
    interblock_sum_L_1_sqrd=0.0_rk
    interblock_sum_L_2_sqrd=0.0_rk
    block_counts_L_1=0 
    block_counts_L_2=0
    equil_L_1=0.0_rk
    sigma_equil_L_1=0.0_rk
    equil_L_2=0.0_rk
    sigma_equil_L_2=0.0_rk

    ! umsd variables. This must be called after n_part is set
    if(allocated(equil_umsd_1)) then
       deallocate(equil_umsd_1)
    end if
    allocate(equil_umsd_1(n_part))
    if(allocated(equil_umsd_2)) then
       deallocate(equil_umsd_2)
    end if
    allocate(equil_umsd_2(n_part))
    if(allocated(sigma_equil_umsd_1)) then
       deallocate(sigma_equil_umsd_1)
    end if
    allocate(sigma_equil_umsd_1(n_part))
    if(allocated(sigma_equil_umsd_2)) then
       deallocate(sigma_equil_umsd_2)
    end if
    allocate(sigma_equil_umsd_2(n_part))
    if(allocated(intrablock_sum_umsd_1)) then
       deallocate(intrablock_sum_umsd_1)
    end if
    allocate(intrablock_sum_umsd_1(n_part))
    if(allocated(intrablock_sum_umsd_2)) then
       deallocate(intrablock_sum_umsd_2)
    end if
    allocate(intrablock_sum_umsd_2(n_part))
    if(allocated(interblock_sum_umsd_1)) then
       deallocate(interblock_sum_umsd_1)
    end if
    allocate(interblock_sum_umsd_1(n_part))
    if(allocated(interblock_sum_umsd_2)) then
       deallocate(interblock_sum_umsd_2)
    end if
    allocate(interblock_sum_umsd_2(n_part))
    if(allocated(interblock_sum_umsd_1_sqrd)) then
       deallocate(interblock_sum_umsd_1_sqrd)
    end if
    allocate(interblock_sum_umsd_1_sqrd(n_part))
    if(allocated(interblock_sum_umsd_2_sqrd)) then
       deallocate(interblock_sum_umsd_2_sqrd)
    end if
    allocate(interblock_sum_umsd_2_sqrd(n_part))
    intrablock_sum_umsd_1=0.0_rk
    intrablock_sum_umsd_2=0.0_rk
    interblock_sum_umsd_1=0.0_rk
    interblock_sum_umsd_2=0.0_rk
    interblock_sum_umsd_1_sqrd=0.0_rk
    interblock_sum_umsd_2_sqrd=0.0_rk
    block_counts_umsd_1=0
    block_counts_umsd_2=0
    equil_umsd_1=0.0_rk
    sigma_equil_umsd_1=0.0_rk
    equil_umsd_2=0.0_rk
    sigma_equil_umsd_2=0.0_rk
    ! End of umsd variables
    sweep_equil_reference=0

  end subroutine initialise_counters




  !! <h4> <code>  subroutine initialise_lattices(filename) </code> </h4>
  !! <p>
  !! This subroutine initialises <code>spec_1</code>, <code>spec_2</code>, <code>R_1</code>, <code>R_2</code>, <code>n_part</code>, 
  !! <code>Lx</code>, <code>Ly</code>, <code>Lz</code> by importing these variables from the specified file. This procedure also allocates 
  !! <code>u</code>, <code>pos_1</code> and <code>pos_2</code>, and initialises <code>switchscalex</code>, <code>switchscaley</code>, 
  !! <code>switchscalez</code> and <code>M_const</code> (assuming <code>beta</code> has already been initialised). 
  !! If any of the aforementioned arrays are already allocated, then this subroutine deallocates them before 
  !! initialising them (and the scalar variables) 'from scratch'. Be aware that this subroutine opens and closes unit 10.
  !! The format of the file <code>filename</code> must be as follows, where <code>n_part</code> is the number of particles in
  !! the supercells pertaining to both lattices 1 and 2; <code>Lx(1)</code> is the length of the supercell for lattice 1 along
  !! the x-direction and similarly for other such quantities; <code>spec_1(i)</code> is the species (integer) of particle i for
  !! lattice 1 and similarly for <code>spec_2(i)</code> for lattice 2; and <code>(R_1(i,1),R_1(i,2),R_1(i,3))</code> is the position of
  !! particle <code>i</code> in lattice 1 <i>in fractional coordinates</i> (note that this is the only situation in which
  !! these quantities are given in fractional coordinates) with respect to the supercell dimensions, and similarly for
  !! <code>(R_2(i,1),R_2(i,2),R_2(i,3))</code> for lattice 2:
  !! <pre><code>
  !! An optional comment of up to 120 characters describing the lattices; if this is unrequired this line can be left blank.
  !! n_part
  !! Lx(1)
  !! Ly(1)
  !! Lz(1)
  !! R_1(1,1) R_1(1,2) R_1(1,3) spec_1(1)
  !! R_1(2,1) R_1(2,2) R_1(2,3) spec_1(2)
  !! R_1(3,1) R_1(3,2) R_1(3,3) spec_1(3)
  !! ...
  !! R_1(n_part,1) R_1(n_part,2) R_1(n_part,3) spec_1(n_part)
  !! Lx(2)
  !! Ly(2)
  !! Lz(2)
  !! R_2(1,1) R_2(1,2) R_2(1,3) spec_2(1)
  !! R_2(2,1) R_2(2,2) R_2(2,3) spec_2(2)
  !! R_2(3,1) R_2(3,2) R_2(3,3) spec_2(3)
  !! ...
  !! R_2(n_part,1) R_2(n_part,2) R_2(n_part,3) spec_2(n_part)
  !! </code></pre>
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> filename </code> </td>
  !!   <td> <code> character(*), intent(in) </code> </td>
  !!   <td> The file from which the lattices will be imported. </td>
  !!  </tr>
  !! </table>
  subroutine initialise_lattices(filename)
    character(*), intent(in) :: filename
    integer(ik) :: error, i
    character(120) :: string
    real(rk) :: scalefactor

    if(allocated(spec_1)) then
       deallocate(spec_1)
    end if
    if(allocated(spec_2)) then
       deallocate(spec_2)
    end if
    if(allocated(pos_1)) then
       deallocate(pos_1)
    end if
    if(allocated(pos_2)) then
       deallocate(pos_2)
    end if
    if(allocated(R_1)) then
       deallocate(R_1)
    end if
    if(allocated(R_2)) then
       deallocate(R_2)
    end if
    if(allocated(u)) then
       deallocate(u)
    end if

    open(unit=10,file=filename,iostat=error,status="old")
    if(error/=0) then
       write(0,*) "monteswitch_mod: Error. Problem opening file '",trim(filename),"'"
       stop 1
    end if

    read(10,*) string
    read(10,*) n_part
    allocate(spec_1(n_part))
    allocate(spec_2(n_part))
    allocate(pos_1(3,n_part))
    allocate(pos_2(3,n_part))
    allocate(R_1(3,n_part))
    allocate(R_2(3,n_part))
    allocate(u(3,n_part))
    read(10,*) Lx(1)
    read(10,*) Ly(1)
    read(10,*) Lz(1)
    do i=1,n_part
       read(10,*) R_1(1,i),R_1(2,i),R_1(3,i), spec_1(i)
       R_1(1,i)=R_1(1,i)*Lx(1)
       R_1(2,i)=R_1(2,i)*Ly(1)
       R_1(3,i)=R_1(3,i)*Lz(1)
    end do
    read(10,*) Lx(2)
    read(10,*) Ly(2)
    read(10,*) Lz(2)
    do i=1,n_part
       read(10,*) R_2(1,i),R_2(2,i),R_2(3,i), spec_2(i)
       R_2(1,i)=R_2(1,i)*Lx(2)
       R_2(2,i)=R_2(2,i)*Ly(2)
       R_2(3,i)=R_2(3,i)*Lz(2)
    end do

    close(unit=10)

    switchscalex = Lx(2)/Lx(1)
    switchscaley = Ly(2)/Ly(1)
    switchscalez = Lz(2)/Lz(1)

    M_const = log( Lx(1)*Ly(1)*Lz(1) / (Lx(2)*Ly(2)*Lz(2)) )/beta

  end subroutine initialise_lattices




  !! <h4> <code> subroutine initialise_cold_microstate(lattice_in) </code> </h4>
  !! <p>
  !! This procedure initialises <code>u</code>, <code>lattice</code>, <code>E</code>, 
  !! <code>E_1</code>, <code>E_2</code>, <code>pos_1</code>, <code>pos_2</code>,
  !! <code>V</code>, <code>M</code>, <code>macro</code> and <code>eta</code> to reflect the microstate with 
  !! <code>u=0</code> for the specified lattice.
  !! </p>
  !! <p>
  !! <b> Dependencies: </b> <code>n_part</code>, <code>R_1</code>, <code>R_2</code>, <code>Lx</code>, <code>Ly</code>, 
  !! <code>Lz</code>, and the 'interactions variables' must be set before this
  !! procedure is called, and <code>u</code> must be allocated.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> lattice_in </code> </td>
  !!   <td> <code>  integer(ik), intent(in) </code> </td>
  !!   <td> The distance - at zero displacement - which defines whether or not two particles interact. </td>
  !!  </tr>
  !! </table>
  subroutine initialise_cold_microstate(lattice_in)
    integer(ik), intent(in) :: lattice_in
    select case(lattice_in)
    case(1)
       u=0.0_rk
       pos_1=R_1
       pos_2=R_2
       E=calc_energy_scratch(1,Lx(1),Ly(1),Lz(1),spec_1,pos_1,R_1,u)
       E_1=E
       E_2=calc_energy_scratch(2,Lx(2),Ly(2),Lz(2),spec_2,pos_2,R_2,u)
       lattice=1
       V=Lx(1)*Ly(1)*Lz(1)
    case(2)
       u=0.0_rk 
       pos_1=R_1
       pos_2=R_2
       E=calc_energy_scratch(2,Lx(2),Ly(2),Lz(2),spec_2,pos_2,R_2,u)
       E_2=E
       E_1=calc_energy_scratch(1,Lx(1),Ly(1),Lz(1),spec_1,pos_1,R_1,u)
       lattice=2
       V=Lx(2)*Ly(2)*Lz(2)
    case default
       write(0,*) "monteswitch_mod: Error. 'lattice_in' is not 1 or 2."
       stop 1
    end select
    M=calc_M(E_1, E_2, Lx(1)*Ly(1)*Lz(1), Lx(2)*Ly(2)*Lz(2))
    macro=get_macro(M)
    ! If M is outwith the supported range then get_macro will return -2 or -1.
    if(macro<1) then
       write(0,*) "monteswitch_mod: Error. Order parameter is outwith supported order parameter range: M = ",M
       stop 1
    end if
    eta=eta_grid(macro)
  end subroutine initialise_cold_microstate



 
  !! <h4> <code> subroutine initialise_barriers() </code> </h4>
  !! <p>
  !! This procedure sets the order parameter barriers such that we are locked into the current macrostate (the
  !! barriers do not span two macrostates). It also zeros counters pertaining to the order parameter barriers,
  !! and sets <code>barrier_to_lock</code> to be 1 or 2 if the <code>barrier_dynamics</code> flag is for <code>"pong_up"</code>
  !! or <code>"pong_down"</code> respectively.
  !! </p>
  !! <p>
  !! <b> Dependencies: </b> <code>M</code>, <code>macro</code> and <code>M_grid</code> should be set before this subroutine is called,
  !! since the barriers are set so that we are locked in the current macrostate. Furthermore <code>M</code> and <code>macro</code>
  !! must be within the supported range.
  !! </p>
  subroutine initialise_barriers()
    ! Set the barriers to encompass the current macrostate
    barrier_macro_low=macro
    barrier_macro_high=barrier_macro_low
    ! Set the correct initial value of barrier_to_lock if 'pong' dynamics are to be used
    select case(barrier_dynamics)
    case(barrier_dynamics_pong_up)
       barrier_to_lock=1
    case(barrier_dynamics_pong_down)
       barrier_to_lock=2
    end select
    ! Zero the relevant counters
    rejected_moves_M_barrier=0
    moves_since_lock=0
  end subroutine initialise_barriers




end module monteswitch_mod
!!
!! </body>
!! </html>
