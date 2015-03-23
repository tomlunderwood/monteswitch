!
! By extracting lines beginning with the regular expression '\!\! ?'
! (ignoring leading whitespace), and then removing matches to the
! regular expression, html documentation corresponding to this 
! source code will be created.
!
!! <html>
!! <head>
!!  <title> monteswitch_mpi Documentation </title>
!! </head>
!! <body>
!! 
!! <h1> <code>monteswitch_mpi</code> Documentation </h1>
!!
!! <h2> Author </h2>
!! <p> Tom Underwood </p>
!!
!! <h2> General Notes </h2>
!! <p>
!! The program is almost identical to <code>monteswitch</code>, but is parallelised using MPI. See the documentation for 
!! <code>monteswitch</code> for more details. Note that this program obviously must be compiled using MPI.
!! </p>
!! <p>
!! The nature of the parallelisation is as follows. Instead of <code>stop_sweeps</code> being performed by a single 
!! task - as is the case in <code>monteswitch</code> - the <code>stop_sweeps</code> are spread between the number of tasks
!! specified when running this program in conjunction with <code>mpiexec</code>. For each task, the initial state of the
!! system is that specified in the 'init' file (for a new simulation) or the 'state' file (for a resumed simulation). Each
!! task has a separate 'data' file to which information is output during its simulation: file 'data_n' corresponds to task 
!! 'n', where n=0,1,2,...  When all tasks are completed, the results are combined in a manner described below, and exported
!! to the file 'state'. Furthermore, if the second command line argument is "-explicit", then the state at the end of each 
!! task is output to a file 'state_n', where 'n' denotes the task number. If the second command line argument is absent,
!! then the files are created, but they are empty. The files can be quite large; hence the default being that they are empty.
!! </p>
!! <p>
!! Details of the parallelisation are as follows. All tasks are initialised to have an identical 'starting state'. For a new simulation this
!! is as determined in the 'init' file. For a resumed simulation this is the corresponding state in the 'state' file.
!! However, the 'counter' variables of all tasks are not identical. Only task '0' - the 'master' task - inherits the counter
!! variables of the starting state. All other tasks have counter variables set to zero. 
!! The simulations for all tasks are then run, during which the counter variables are ammended independently for each task. 
!! Once all simulations are completed then their states are exported to the 'state_n' files. The counter variables of all tasks 
!! other than task 0 are then summed and stored in task 0. Task 0 then has a state corresponding to an evolution of, if 
!! there are <code>T</code> tasks, approximately <code>stop_sweeps/T</code> sweeps from the starting state, but with counters 
!! drawn from <code>stop_sweeps</code> sweeps worth of information.
!! The state of task 0 is then updated to reflect its new counters, e.g., the weight function is recalculated using the information
!! from the new counters. This state is then exported to 'state'.
!! </p>
!! <p>
!! Note that 'data_n' are always overwritten by subsequent simulations since, given the nature of the parallelisation,
!! there is no continuity between the tasks in subsequent simulations.
!! </p>
!!
program monteswitch_mpi

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


  ! MPI header
  include 'mpif.h'

  ! MPI-related variables
  integer :: task_id ! The task ID
  integer :: num_tasks ! The number of tasks

  ! The total number of sweeps to be performed - the same as 'stop_sweeps' in the 'init' or 'state' file
  integer(ik) :: stop_sweeps_total
  ! The string of the 'data' file for the current task
  character(len=20) :: filename_data
  ! The string of the 'state' file for the current task
  character(len=20) :: filename_state
  ! The checkpoint_period flag 'for the whole simulation'
  integer(ik) :: checkpoint_period_overall
  ! The seed to be used for the random number generator. In the module 'rng_mod', if a specific seed is not given,
  ! then a random one is chosen - the value of the system clock. One cannot rely on the clock being different
  ! when the initalisation procedure in 'rng_mod' is called by each task; hence one should explicitly give the
  ! random number generators for each task different seeds. Otherwise it is possible that two or more tasks have
  ! the same stream of random numbers.
  integer :: seed
  
  ! Unimportant variables
  character(len=20) :: char
  character(len=20) :: char_int
  integer(ik) :: error

  ! THE PROGRAM
  
 
  ! Initialise the MPI environment
  call MPI_INIT(error)
  if(error/=0) then
     write(0,*) "monteswitch_mpi: Error. There was a problem initialising the MPI environment."
     stop 1
  end if

  ! Find out the task ID and the number of processes which were started
  call MPI_COMM_RANK(MPI_COMM_WORLD,task_id,error)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_tasks,error)

  ! Set the filenames for this task
  write(char_int,*) task_id
  filename_data='data_'//trim(adjustl(char_int))
  filename_state='state_'//trim(adjustl(char_int))
  ! Set the seed (different values for each task)
  call system_clock(count=seed)
  seed=seed+task_id

  ! Read the first command line argument 
  call getarg(1,char)
  if(char=="") then
     write(0,*) "monteswitch_mpi: Error. No command line argument detected."
     stop 1
  else if(trim(char)=="-new") then

     ! CODE FOR A NEW SIMULATION

     ! Import the appropriate variables from the 3 input files
     call initialise_from_files("params_in","lattices_in","interactions_in")
     
     ! Assign sweeps according to the task ID
     call assign_sweeps()
     ! Assign checkpoint_period according to the second command line argument
     call assign_checkpoint_period()
     ! Run
     call run(filename_data,filename_state,.false.,seed)
     ! Wait for all tasks to complete before proceeding
     call MPI_BARRIER(MPI_COMM_WORLD,error) 
     ! Combine the counter variables in all simulations into task 0, process, then export to 'state'
     call combine()

  else if(trim(char)=="-resume") then

     ! CODE FOR A RESUMED SIMULATION
     
     call import("state")

     ! Assign sweeps according to the task ID
     call assign_sweeps()
     ! Reset counters for all tasks other than task 0
     if(task_id/=0) then
         call initialise_counters()
     end if
     ! Assign checkpoint_period according to the second command line argument
     call assign_checkpoint_period()
     ! Run
     call run(filename_data,filename_state,.false.,seed)
     ! Wait for all tasks to complete before proceeding
     call MPI_BARRIER(MPI_COMM_WORLD,error) 
     ! Combine the counter variables in all simulations into task 0, process, then export to 'state'
     call combine()

  else if(trim(char)=="-reset") then
     
     ! CODE FOR SIMULATION INITIALISED FROM THE STATE FILE BUT WITH COUNTERS RESET
     
     call import("state")
     call initialise_counters()

     ! Assign sweeps according to the task ID
     call assign_sweeps()
     ! Assign checkpoint_period according to the second command line argument
     call assign_checkpoint_period()
     ! Run
     call run(filename_data,filename_state,.false.,seed)
     ! Combine the counter variables in all simulations into task 0, process, then export to 'state'
     call combine()

  else

     ! CODE FOR ANY OTHER COMMAND LINE ARGUMENTS

     write(0,*) "monteswitch_mpi: Error. Unrecognised first command line argument."
     stop 1

  end if



  ! Finalise the MPI environment
  call MPI_FINALIZE(error)
  if(error/=0) then
     write(0,*) "monteswitch_mpi: Error. There was a problem finalising the MPI environment."
     stop 1
  end if




contains

  
  ! This subroutine examines the second command argument. If it is
  ! "-explicit", then we are in explicit mode, and checkpoint_period takes the expected value; if it is absent then
  ! we want empty state files exported by all tasks, and hence we set checkpoint_period to be -1.
  subroutine assign_checkpoint_period()
    ! Get the original value for 'checkpoint_period'
    checkpoint_period_overall=checkpoint_period
    ! Read the second command line argument, and set 'explicit mode' accordingly.
    call getarg(2,char)
    if(trim(char)=="-explicit") then     
       ! Do nothing (checkpoint_period for all tasks is that specified in the imported 'state' or 'init')
    else if(char=="") then
       ! Have empty state files
       checkpoint_period=-1
    else 
       write(0,*) "monteswitch_mpi: Error. Unrecognised second argument detected."
       stop 1
    end if
  end subroutine assign_checkpoint_period



  ! This subroutine sets stop_sweeps for the current task: for each task we perform stop_sweeps_total/num_tasks sweeps; 
  ! with any extra sweeps required to bring us to a total of stop_sweeps_total assigned to the master task (task 0)
  subroutine assign_sweeps()
    stop_sweeps_total=stop_sweeps
    stop_sweeps=stop_sweeps_total/num_tasks
    if(task_id==0) then
       stop_sweeps = stop_sweeps + ( stop_sweeps_total-(stop_sweeps_total/num_tasks)*num_tasks )
    end if
  end subroutine assign_sweeps




  ! This subroutine combines 'counter' variables from all tasks with those in the master task. Note that
  ! 'intrablock' sum variables are not combined - since they shouldn't be. It then runs task 0 for 0 sweeps 
  ! in order to recalculate weight functions, and finally exports the final master task 
  ! variables to the file 'state'.
  ! Note that the data type used in MPI_REDUCE for the real variables is MPI_DOUBLE_PRECISION. Using MPI_REAL,
  ! with real(rk) amounting to double precision, caused errors.
  subroutine combine()
    integer(ik) :: sum_int
    real(rk) :: sum_real
    integer(ik), dimension(:), allocatable :: sum_int_array
    real(rk), dimension(:), allocatable :: sum_real_array
    real(rk), dimension(:,:), allocatable :: trans_sum

    ! Block the current task until all other tasks get to this point
    call MPI_BARRIER(MPI_COMM_WORLD,error)

    ! For each counter variable, sum over all variables in each task, and store them in the relevant 'sum' variable 
    ! in process 0.

    ! M_counts_1
    allocate(sum_int_array(size(M_counts_1)))
    sum_int_array=0
    call MPI_REDUCE(M_counts_1,sum_int_array,size(sum_int_array),MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       M_counts_1=sum_int_array
    end if
    deallocate(sum_int_array)
    ! M_counts_2
    allocate(sum_int_array(size(M_counts_2)))
    sum_int_array=0
    call MPI_REDUCE(M_counts_2,sum_int_array,size(sum_int_array),MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       M_counts_2=sum_int_array
    end if
    deallocate(sum_int_array)
    ! sweeps
    sum_int=0
    call MPI_REDUCE(sweeps,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       sweeps=sum_int
    end if
    ! moves
    sum_int=0
    call MPI_REDUCE(moves,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       moves=sum_int
    end if
    ! moves_lattice
    sum_int=0
    call MPI_REDUCE(moves_lattice,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       moves_lattice=sum_int
    end if
    ! accepted_moves_lattice
    sum_int=0
    call MPI_REDUCE(accepted_moves_lattice,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       accepted_moves_lattice=sum_int
    end if
    ! moves_part
    sum_int=0
    call MPI_REDUCE(moves_part,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       moves_part=sum_int
    end if
    ! accepted_moves_part
    sum_int=0
    call MPI_REDUCE(accepted_moves_part,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       accepted_moves_part=sum_int
    end if
    ! moves_vol
    sum_int=0
    call MPI_REDUCE(moves_vol,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       moves_vol=sum_int
    end if
    ! accepted_moves_vol
    sum_int=0
    call MPI_REDUCE(accepted_moves_vol,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       accepted_moves_vol=sum_int
    end if
    ! rejected_moves_M_OOB
    sum_int=0
    call MPI_REDUCE(rejected_moves_M_OOB,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       rejected_moves_M_OOB=sum_int
    end if
    ! melts
    sum_int=0
    call MPI_REDUCE(melts,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       melts=sum_int
    end if
    ! trans
    allocate(  trans_sum( size(trans,1) , size(trans,2) )  )
    trans_sum=0.0_rk
    call MPI_REDUCE(trans,trans_sum,size(trans_sum),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       trans=trans_sum
    end if
    deallocate(trans_sum)
    ! block_counts
    sum_int=0
    call MPI_REDUCE(block_counts,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       block_counts=sum_int
    end if
    ! interblock_sum_DeltaF
    sum_real=0.0_rk
    call MPI_REDUCE(interblock_sum_DeltaF,sum_real,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
      interblock_sum_DeltaF=sum_real
    end if
    ! interblock_sum_DeltaF_sqrd
    sum_real=0.0_rk
    call MPI_REDUCE(interblock_sum_DeltaF_sqrd,sum_real,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
      interblock_sum_DeltaF_sqrd=sum_real
    end if
    ! block_counts_DeltaF
    sum_int=0
    call MPI_REDUCE(block_counts_DeltaF,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       block_counts_DeltaF=sum_int
    end if
    ! interblock_sum_H_1
    sum_real=0.0_rk
    call MPI_REDUCE(interblock_sum_H_1,sum_real,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
      interblock_sum_H_1=sum_real
    end if
    ! interblock_sum_H_2
    sum_real=0.0_rk
    call MPI_REDUCE(interblock_sum_H_2,sum_real,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
      interblock_sum_H_2=sum_real
    end if
    ! interblock_sum_H_1_sqrd
    sum_real=0.0_rk
    call MPI_REDUCE(interblock_sum_H_1_sqrd,sum_real,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
      interblock_sum_H_1_sqrd=sum_real
    end if
    ! interblock_sum_H_2_sqrd
    sum_real=0.0_rk
    call MPI_REDUCE(interblock_sum_H_2_sqrd,sum_real,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
      interblock_sum_H_2_sqrd=sum_real
    end if
    ! block_counts_H_1
    sum_int=0
    call MPI_REDUCE(block_counts_H_1,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       block_counts_H_1=sum_int
    end if
    ! block_counts_H_2
    sum_int=0
    call MPI_REDUCE(block_counts_H_2,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       block_counts_H_2=sum_int
    end if
    ! interblock_sum_V_1
    sum_real=0.0_rk
    call MPI_REDUCE(interblock_sum_V_1,sum_real,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
      interblock_sum_V_1=sum_real
    end if
    ! interblock_sum_V_2
    sum_real=0.0_rk
    call MPI_REDUCE(interblock_sum_V_2,sum_real,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
      interblock_sum_V_2=sum_real
    end if
    ! interblock_sum_V_1_sqrd
    sum_real=0.0_rk
    call MPI_REDUCE(interblock_sum_V_1_sqrd,sum_real,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
      interblock_sum_V_1_sqrd=sum_real
    end if
    ! interblock_sum_V_2_sqrd
    sum_real=0.0_rk
    call MPI_REDUCE(interblock_sum_V_2_sqrd,sum_real,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
      interblock_sum_V_2_sqrd=sum_real
    end if
    ! block_counts_V_1
    sum_int=0
    call MPI_REDUCE(block_counts_V_1,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       block_counts_V_1=sum_int
    end if
    ! block_counts_V_2
    sum_int=0
    call MPI_REDUCE(block_counts_V_2,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       block_counts_V_2=sum_int
    end if
    ! interblock_sum_umsd_1
    allocate(sum_real_array(n_part))
    sum_real_array=0
    call MPI_REDUCE(interblock_sum_umsd_1,sum_real_array,size(sum_real_array),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       interblock_sum_umsd_1=sum_real_array
    end if
    deallocate(sum_real_array)
    ! interblock_sum_umsd_2
    allocate(sum_real_array(n_part))
    sum_real_array=0
    call MPI_REDUCE(interblock_sum_umsd_2,sum_real_array,size(sum_real_array),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       interblock_sum_umsd_2=sum_real_array
    end if
    deallocate(sum_real_array)
    ! interblock_sum_umsd_1_sqrd
    allocate(sum_real_array(n_part))
    sum_real_array=0
    call MPI_REDUCE(interblock_sum_umsd_1_sqrd,sum_real_array,size(sum_real_array), &
         MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       interblock_sum_umsd_1_sqrd=sum_real_array
    end if
    deallocate(sum_real_array)
    ! interblock_sum_umsd_2_sqrd
    allocate(sum_real_array(n_part))
    sum_real_array=0
    call MPI_REDUCE(interblock_sum_umsd_2_sqrd,sum_real_array,size(sum_real_array), &
         MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       interblock_sum_umsd_2_sqrd=sum_real_array
    end if
    deallocate(sum_real_array)
    ! block_counts_umsd_1
    sum_int=0
    call MPI_REDUCE(block_counts_umsd_1,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       block_counts_umsd_1=sum_int
    end if
    ! block_counts_umsd_2
    sum_int=0
    call MPI_REDUCE(block_counts_umsd_2,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,error)
    if(task_id==0) then
       block_counts_umsd_2=sum_int
    end if

    ! In task 0, perform the final things
    if(task_id==0) then
       ! Run the simulation once with a stop_sweeps of 0 to update the weight function and
       ! equilibrium properties.
       stop_sweeps=0
       call run(filename_data,"state",.true.)
       ! Set stop_sweeps back to the value before parallelisation, which is stored in stop_sweeps_total
       stop_sweeps=stop_sweeps_total
       ! Set checkpoint_period back to the value before parallelisation, which is stored 
       ! in checkpoint_period_overall. Recall that 'checkpoint_period' was modified to be -1, unless we
       ! are in explicit mode.
       checkpoint_period=checkpoint_period_overall
       ! Export the simulation to 'state'
       call export("state")
    end if

  end subroutine combine




end program monteswitch_mpi
!!
!! </body>
!! </html>
