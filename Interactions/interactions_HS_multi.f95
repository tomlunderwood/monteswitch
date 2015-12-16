! 'interactions.f95' file for monteswitch for multiple-component hard-sphere mixtures, in conjunction with 
! fixed neighbour lists.
! 
! Author: Tom L Underwood
!
! This file was adapted from the file 'interactions_TEMPLATE_pair.f95', but note that 'cutoff' is superfluous
! for hard-sphere interactions, as is the function 'pair_potential_trunc', and hence I have removed them.
! The variables to be specified in the 'interactions_in' file (see below) are as follows (in the following 
! order):
! 
! * epsilon (real(rk)) is the energy cost associated with overlappng spheres. This will normally be set to a 
!   'high' value.
!
! * n_species (integer(ik)) is the number of species to define hard-sphere diameters for (see the next 
!   variable).
!
! * sigma (real(rk), dimension(n_species)) are the hard-sphere diameters for each species: sigma(i) is the
!   diameter for species 'i'. Note that species labels range from 1 to 'n_species' (e.g., species '0' or '-2'
!   are not used).
!
! * list_cutoff (real(rk)) is the cut-off distance determining whether pairs of particles interact with each 
!   other throughout the simulation. (See below).
!
! * list_size (integer(ik)) determines the maximum number of particles any one particle is 'allowed' to 
!   interact with via the 'list' mechanism (determined by the 'list_cutoff' variable above). (See below).
!
!
!
!
! COMMENTS PRESENT IN THE TEMPLATE BEFORE I MODIFIED IT...
!
! This file corresponds to the Lennard-Jones potential until modified by the user. Search for 
! 'USER-DEFINED CODE' to find the relevant parts to modify. What follows is a description of functionality
! for the Lennard-Jones potential. The functionality will be similar if the 'USER-DEFINED CODE' blocks
! are altered, except the variables 'lj_epsilon' and 'lj_sigma' will be replaced by user-defined variables.
!
! Description of functionality of this file if it is unaltered by the user...
!
! This file implements a Lennard-Jones potential with truncated interactions, and a fixed neighbour list for
! each particle. The variables for this module are imported from a file 'interactions_in', and the format of 
! this file is as follows. On the first line there are two tokens. The first is a character(len=20) variable
! (I recommend: 'lj_epsilon='); the second is the value of 'lj_epsilon' (all variables are explained in a 
! moment). The second line is similar, but for 'lj_sigma'. The third line is for 'cutoff', the fourth is for 
! 'list_cutoff', and the fifth is for 'list_size'. The variables are as follows:
!
! * lj_epsilon (real(rk)) is the depth of the Lennard-Jones potential well.
!
! * lj_sigma (real(rk)) is the distance corresponding to 0 potential for the Lennard-Jones potential.
!
! * cutoff (real(rk)) is the cut-off distance for the Lennard-Jones potential.
!
! * list_cutoff (real(rk)) is the cut-off distance determining whether pairs of particles interact with each 
!   other throughout the simulation. Those within list_cutoff of each other at the start of the simulation, 
!   before any moves are made, will interact with each other forever more. Note that the set of interacting 
!   pairs does not change during the simulation, even if pairs of interacting particles later exceed 
!   'list_cutoff' in separation.
!
! * list_size (integer(ik)) determines the maximum number of particles any one particle is 'allowed' to 
!   interact with via the 'list' mechanism (determined by the 'list_cutoff' variable above). This should be
!   set to at least the number of neghbours one will interact with + 2. E.g., if 'list_cutoff' is set to (just over) 
!   the nearest neighbour distance, and there are 12 nearest neighbours, then 'list_size' should 
!   be set to at least 14. While there is nothing wrong in principle with setting 'list_size' to, say, 500, 
!   the associated arrays would be very large (500 integers per particle), and hence could slow down the 
!   simulation and/or use up too much memory.
!
! Regarding checkpointing, the variables in this module are stored in the 'state' file used to checkpoint all
! other monteswitch variables; see the comments for 'export_interactions' for more details.
!
module interactions_mod


    ! The module 'kinds_mod' (in 'kinds_mod.f95') contains the real and integer kinds
    ! which monteswitch uses: real(kind=rk) and integer(kind=ik). Please use these
    ! kinds throughout this module to allow portability.
    use kinds_mod
    
    implicit none

    ! <<<<<<<<<<<< USER-DEFINED CODE >>>>>>>>>>>> 
    ! Insert code defining the free parameters for the potential here.
    ! Example (for Lennard-Jones potential):
    !  real(rk) :: lj_epsilon
    !  real(rk) :: lj_sigma
    !
    real(rk) :: epsilon
    integer(ik) :: n_species
    real(rk), dimension(:), allocatable :: sigma
    ! <<<<<<<<<<<< END OF USER-DEFINED CODE >>>>>>>>>>>> 

    ! Variables for the potential cutoff and lists of interacting particles
    !real(rk) :: cutoff
    real(rk) :: list_cutoff
    integer(ik) :: list_size
    integer(ik), dimension(:,:), allocatable :: list_1
    integer(ik), dimension(:,:), allocatable :: list_2

    ! The below public procedures are called by 'monteswitch_mod' (in 'monteswitch_mod.f95'),
    ! and must be 'filled in' by the user.
    private
    public :: initialise_interactions, export_interactions, import_interactions, &
        after_accepted_part_interactions, after_accepted_vol_interactions, &
        after_accepted_lattice_interactions, after_all_interactions, &
        calc_energy_scratch, calc_energy_part_move


contains




! Initialises the  variables in this module. The initial configurations for both
! lattices are provided as arguments in case this information is required. Unless the interactions
! are 'hard-coded', the initialisation will involve reading from a file. This should be done here.
subroutine initialise_interactions(Lx1, Ly1, Lz1, species1, pos1, Lx2, Ly2, Lz2, species2, pos2)
    ! Dimensions of the initial (orthorhombic) supercell for lattices 1 and 2
    ! in each Cartesian dimension
    real(rk), intent(in) :: Lx1, Ly1, Lz1, Lx2, Ly2, Lz2
    ! Array containing the species of each particle for each lattice: e.g., species1(i) is the
    ! species of particle i in lattice 1
    integer(ik), intent(in), dimension(:) :: species1, species2
    ! Initial positions (Cartesian) of the particles for lattices 1 and 2: e.g., pos1(i,1) is
    ! the x-coordinate of particle i in lattice 1, pos1(i,2) is the y-coordinate, and pos1(i,3)
    ! is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2

    character(len=20) string
    integer(ik) :: error
    integer(ik) :: i,j,n,n_part

    n_part = size(pos1,1)

    ! Read the interactions variables from the file 'interactions_in'
    open(unit=10,file="interactions_in",iostat=error,status="old")
    if(error/=0) then
       write(0,*) "interactions: Error. Problem opening file 'interactions_in'"
       stop 1
    end if

    ! <<<<<<<<<<<< USER-DEFINED CODE >>>>>>>>>>>> 
    ! Insert code which reads the free parameters for the potential from unit 10, using list-directed
    ! formatting. Note that the way in which the free parameters are read must match the way in which
    ! they are written - in the USER-DEFINED CODE block below.
    ! Example (for Lennard-Jones potential, including an error message and exit code of 1 if there is
    ! an error):
    !  read(10,*,iostat=error) string, lj_epsilon
    !  if(error/=0) then
    !     write(0,*) "interactions: Error. Problem reading 'lj_epsilon' from file '",filename,"'"
    !     stop 1
    !  end if
    !  read(10,*,iostat=error) string, lj_sigma
    !  if(error/=0) then
    !     write(0,*) "interactions: Error. Problem reading 'lj_sigma' from file '",filename,"'"
    !     stop 1
    !  end if
    !
    read(10,*,iostat=error) string, epsilon
    if(error/=0) then
        write(0,*) "interactions: Error. Problem reading 'epsilon' from file 'interactions'"
        stop 1
    end if
    read(10,*,iostat=error) string, n_species
    if(error/=0) then
        write(0,*) "interactions: Error. Problem reading 'n_species' from file 'interactions'"
        stop 1
    end if
    if(allocated(sigma)) deallocate(sigma)
    allocate(sigma(n_species))
    read(10,*,iostat=error) string, sigma
    if(error/=0) then
        write(0,*) "interactions: Error. Problem reading 'sigma' from file 'interactions'"
        stop 1
    end if
    ! <<<<<<<<<<<< END OF USER-DEFINED CODE >>>>>>>>>>>> 

    !read(10,*,iostat=error) string, cutoff
    !if(error/=0) then
    !   write(0,*) "interactions: Error. Problem reading 'cutoff' from file 'interactions_in'"
    !   stop 1
    !end if
    read(10,*,iostat=error) string, list_cutoff
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'list_cutoff' from file 'interactions_in'"
       stop 1
    end if
    read(10,*,iostat=error) string, list_size
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'list_size' from file 'interactions_in'"
       stop 1
    end if
    close(10)

    ! Initialise lists
    if(allocated(list_1)) then
       deallocate(list_1)
    end if
    if(allocated(list_2)) then
       deallocate(list_2)
    end if
    allocate(list_1(n_part,list_size))
    allocate(list_2(n_part,list_size))

    list_1=0
    do i=1,n_part
       n=1
       do j=1,n_part
          if(min_image_distance(pos1(i,:),pos1(j,:),Lx1,Ly1,Lz1)<list_cutoff) then
             list_1(i,n)=j
             n=n+1
          end if
       end do
    end do
    list_2=0
    do i=1,n_part
       n=1
       do j=1,n_part
          if(min_image_distance(pos2(i,:),pos2(j,:),Lx2,Ly2,Lz2)<list_cutoff) then
             list_2(i,n)=j
             n=n+1
          end if
       end do
    end do

end subroutine initialise_interactions




! Export all of the variables in this module - for the purposes of checkpointing. There are 
! two options for this. If one outputs the variables to unit 10, without opening or closing that 
! unit in this procedure, then the variables will be output to, and stored within,
! the 'state' file which contains all other monteswitch variables, and is used for resuming
! old simulations as well as post-processing. If one wishes to use a separate file or files to 
! checkpoint the variables in this module, then one is free to do so in this procedure - 
! just don't use unit 10! Note that the format output to the checkpoint file(s) by this procedure
! must correspond to the format read in from these files by the procedure 'import_interactions()'.
subroutine export_interactions()

    ! <<<<<<<<<<<< USER-DEFINED CODE >>>>>>>>>>>> 
    ! Insert code which writes the free parameters for the potential to unit 10, using list-directed
    ! formatting. Note that the way in which the free parameters are read must match the way in which
    ! they are read - in the USER-DEFINED CODE block above.
    ! Example (for Lennard-Jones potential):
    !  write(10,*) "lj_epsilon= ",lj_epsilon
    !  write(10,*) "lj_sigma= ",lj_sigma
    !
    write(10,*) "epsilon= ",epsilon
    write(10,*) "n_species= ",n_species
    write(10,*) "sigma= ",sigma
    ! <<<<<<<<<<<< END OF USER-DEFINED CODE >>>>>>>>>>>> 

    !write(10,*) "cutoff= ",cutoff
    write(10,*) "list_cutoff= ",list_cutoff
    write(10,*) "list_size= ",list_size
    write(10,*) "list_1= ",list_1
    write(10,*) "list_2= ",list_2
    
end subroutine export_interactions




! Import all variables in this module from file(s) - for the purposes of checkpointing. The format 
! read from the file(s) should correspond to that output by 'export_interactions' above. If one is
! importing from within a 'state' file as described in that procedure, then use unit 10, but
! do not open or close that unit!
subroutine import_interactions(Lx1, Ly1, Lz1, species1, pos1, R1, Lx2, Ly2, Lz2, species2, pos2, R2, u)
    ! Dimensions of the (orthorhombic) supercell for lattices 1 and 2
    ! in each Cartesian dimension
    real(rk), intent(in) :: Lx1, Ly1, Lz1, Lx2, Ly2, Lz2
    ! Array containing the species of each particle for each lattice: e.g., species1(i) is the
    ! species of particle i in lattice 1
    integer(ik), intent(in), dimension(:) :: species1, species2
    ! Positions (Cartesian) of the particles for lattices 1 and 2: e.g., pos1(i,1) is
    ! the x-coordinate of particle i in lattice 1, pos1(i,2) is the y-coordinate, and pos1(i,3)
    ! is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2
    ! Positions (Cartesian) of the lattice sites for lattices 1 and 2: e.g., R1(i,1) is
    ! the x-coordinate of the lattice site for particle i in lattice 1, ..., R2(i,1) is
    ! the x-coordinate of the lattice site for particle i in lattice 2
    real(rk), intent(in), dimension(:,:) :: R1, R2
    ! Displacements of the particles: e.g., u(i,1) is the x-coordinate of particle i in lattice 1
    ! etc.
    real(rk), intent(in), dimension(:,:) :: u

    character(len=20) string
    integer(ik) :: error, n_part
    
    ! <<<<<<<<<<<< USER-DEFINED CODE >>>>>>>>>>>>
    ! Insert code which reads the free parameters for the potential from unit 10, using list-directed
    ! formatting. Note that the way in which the free parameters are read must match the way in which
    ! they are written - in the USER-DEFINED CODE block below.
    ! Example (for Lennard-Jones potential, including an error message and exit code of 1 if there is
    ! an error):
    !  read(10,*,iostat=error) string, lj_epsilon
    !  if(error/=0) then
    !     write(0,*) "interactions: Error. Problem reading 'lj_epsilon' from unit ",unit
    !     stop 1
    !  end if
    !  read(10,*,iostat=error) string, lj_sigma
    !  if(error/=0) then
    !     write(0,*) "interactions: Error. Problem reading 'lj_sigma' from unit ",unit
    !     stop 1
    !  end if
    !
    read(10,*,iostat=error) string, epsilon
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'epsilon' from unit 10"
       stop 1
    end if
    read(10,*,iostat=error) string, n_species
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'n_species' from unit 10"
       stop 1
    end if
    if(allocated(sigma)) deallocate(sigma)
    allocate(sigma(n_species))
    read(10,*,iostat=error) string, sigma
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'sigma' from unit 10"
       stop 1
    end if
    ! <<<<<<<<<<<< END OF USER-DEFINED CODE >>>>>>>>>>>> 

    !read(10,*,iostat=error) string, cutoff
    !if(error/=0) then
    !   write(0,*) "interactions: Error. Problem reading 'cutoff' from unit 10"
    !   stop 1
    !end if
    read(10,*,iostat=error) string, list_cutoff
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'list_cutoff' from unit 10"
       stop 1
    end if
    read(10,*,iostat=error) string, list_size
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'list_size' from unit 10"
       stop 1
    end if

    n_part = size(pos1,1)
    if(allocated(list_1)) then
       deallocate(list_1)
    end if
    if(allocated(list_2)) then
       deallocate(list_2)
    end if
    allocate(list_1(n_part,list_size))
    allocate(list_2(n_part,list_size))

    read(10,*,iostat=error) string, list_1
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'list_1' from unit 10"
       stop 1
    end if
    read(10,*,iostat=error) string, list_2
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'list_2' from unit 10"
       stop 1
    end if

end subroutine import_interactions




! Performs any tasks required for the variables in this module, after a particle move for particle i is accepted
subroutine after_accepted_part_interactions(i, Lx1, Ly1, Lz1, species1, pos1, R1, Lx2, Ly2, Lz2, species2, pos2, R2, u)
    ! The particle which has just been moved (where the move has been accepted)
    integer(ik), intent(in) :: i
    ! Current dimensions of the initial (orthorhombic) supercell for lattices 1 and 2
    ! in each Cartesian dimension (after the move has been accepted)
    real(rk), intent(in) :: Lx1, Ly1, Lz1, Lx2, Ly2, Lz2
    ! Array containing the species of each particle for each lattice: e.g., species1(i) is the
    ! species of particle i in lattice 1
    integer(ik), intent(in), dimension(:) :: species1, species2
    ! Current positions (Cartesian) of the particles for lattices 1 and 2 (after the move has been
    ! accepted): e.g., pos1(i,1) is the x-coordinate of particle i in lattice 1, pos1(i,2) is the 
    ! y-coordinate, and pos1(i,3) is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2
    ! Positions (Cartesian) of the lattice sites for lattices 1 and 2: e.g., R1(i,1) is
    ! the x-coordinate of the lattice site for particle i in lattice 1, ..., R2(i,1) is
    ! the x-coordinate of the lattice site for particle i in lattice 2
    real(rk), intent(in), dimension(:,:) :: R1, R2
    ! Current displacement vectors for the particles (after the move has been accepted); e.g., 
    ! u(i,1) is the x-displacement of particle 1 from its lattice site, etc.
    real(rk) ,intent(in), dimension(:,:) :: u

    return

end subroutine after_accepted_part_interactions




! Performs any tasks required for the variables in this module, after a volume move is accepted
subroutine after_accepted_vol_interactions(Lx1, Ly1, Lz1, species1, pos1, R1, Lx2, Ly2, Lz2, species2, pos2, R2, u)
    ! Current dimensions of the initial (orthorhombic) supercell for lattices 1 and 2
    ! in each Cartesian dimension (after the move has been accepted)
    real(rk), intent(in) :: Lx1, Ly1, Lz1, Lx2, Ly2, Lz2
    ! Array containing the species of each particle for each lattice: e.g., species1(i) is the
    ! species of particle i in lattice 1
    integer(ik), intent(in), dimension(:) :: species1, species2
    ! Current positions (Cartesian) of the particles for lattices 1 and 2 (after the move has been
    ! accepted): e.g., pos1(i,1) is the x-coordinate of particle i in lattice 1, pos1(i,2) is the 
    ! y-coordinate, and pos1(i,3) is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2
    ! Positions (Cartesian) of the lattice sites for lattices 1 and 2: e.g., R1(i,1) is
    ! the x-coordinate of the lattice site for particle i in lattice 1, ..., R2(i,1) is
    ! the x-coordinate of the lattice site for particle i in lattice 2
    real(rk), intent(in), dimension(:,:) :: R1, R2
    ! Current displacement vectors for the particles (after the move has been accepted); e.g., 
    ! u(i,1) is the x-displacement of particle 1 from its lattice site, etc.
    real(rk) ,intent(in), dimension(:,:) :: u


    return

end subroutine after_accepted_vol_interactions



! Performs any tasks required for the variables in this module, after a lattice move is accepted
subroutine after_accepted_lattice_interactions(Lx1, Ly1, Lz1, species1, pos1, R1, Lx2, Ly2, Lz2, species2, pos2, R2, u)
    ! Current dimensions of the initial (orthorhombic) supercell for lattices 1 and 2
    ! in each Cartesian dimension (after the move has been accepted)
    real(rk), intent(in) :: Lx1, Ly1, Lz1, Lx2, Ly2, Lz2
    ! Array containing the species of each particle for each lattice: e.g., species1(i) is the
    ! species of particle i in lattice 1
    integer(ik), intent(in), dimension(:) :: species1, species2
    ! Current positions (Cartesian) of the particles for lattices 1 and 2 (after the move has been
    ! accepted): e.g., pos1(i,1) is the x-coordinate of particle i in lattice 1, pos1(i,2) is the 
    ! y-coordinate, and pos1(i,3) is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2
    ! Positions (Cartesian) of the lattice sites for lattices 1 and 2: e.g., R1(i,1) is
    ! the x-coordinate of the lattice site for particle i in lattice 1, ..., R2(i,1) is
    ! the x-coordinate of the lattice site for particle i in lattice 2
    real(rk), intent(in), dimension(:,:) :: R1, R2
    ! Current displacement vectors for the particles (after the move has been accepted); e.g., 
    ! u(i,1) is the x-displacement of particle 1 from its lattice site, etc.
    real(rk) ,intent(in), dimension(:,:) :: u


    return

end subroutine after_accepted_lattice_interactions




! Performs any tasks required for the variables in this module after ALL moves (accepted or rejected)
subroutine after_all_interactions(Lx1, Ly1, Lz1, species1, pos1, R1, Lx2, Ly2, Lz2, species2, pos2, R2, u)
    ! Current dimensions of the initial (orthorhombic) supercell for lattices 1 and 2
    ! in each Cartesian dimension
    real(rk), intent(in) :: Lx1, Ly1, Lz1, Lx2, Ly2, Lz2
    ! Array containing the species of each particle for each lattice: e.g., species1(i) is the
    ! species of particle i in lattice 1
    integer(ik), intent(in), dimension(:) :: species1, species2
    ! Current positions (Cartesian) of the particles for lattices 1 and 2: e.g., pos1(i,1) is 
    ! the x-coordinate of particle i in lattice 1, pos1(i,2) is the y-coordinate, and pos1(i,3) 
    ! is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2
    ! Positions (Cartesian) of the lattice sites for lattices 1 and 2: e.g., R1(i,1) is
    ! the x-coordinate of the lattice site for particle i in lattice 1, ..., R2(i,1) is
    ! the x-coordinate of the lattice site for particle i in lattice 2
    real(rk), intent(in), dimension(:,:) :: R1, R2
    ! Current displacement vectors for the particles (after the move has been accepted); e.g., 
    ! u(i,1) is the x-displacement of particle 1 from its lattice site, etc.
    real(rk) ,intent(in), dimension(:,:) :: u


    return

end subroutine after_all_interactions




! Returns the energy for the specificed configuration of the specified lattice
function calc_energy_scratch(lattice, Lx, Ly, Lz, species, pos, R, u)
    ! The lattice (1 or 2)
    integer(ik), intent(in) :: lattice
    ! Dimensions of the (orthorhombic) supercell
    real(rk), intent(in) :: Lx, Ly, Lz
    ! Array containing the species of each particle for each lattice: e.g., species(i) is the
    ! species of particle i
    integer(ik), dimension(:), intent(in) :: species
    ! Positions (Cartesian) of the particles: e.g., pos(i,1) is the x-coordinate of particle 
    ! i, pos1(i,2) is the y-coordinate, and pos1(i,3) is the z-coordinate
    real(rk), dimension(:,:), intent(in) :: pos
    ! Positions (Cartesian) of the lattice sites for the configuration: e.g., R(i,1) is
    ! the x-coordinate of the lattice site for particle i, etc.
    real(rk), intent(in), dimension(:,:) :: R
    ! Displacement vectors for the particles; e.g., u(i,1) is the x-displacement of 
    ! particle 1 from its lattice site, etc.
    real(rk) ,intent(in), dimension(:,:) :: u

    real(rk) :: calc_energy_scratch

    select case(lattice)
    case(1)
       calc_energy_scratch = system_energy(species,Lx,Ly,Lz,list_1,pos)
    case(2)
       calc_energy_scratch = system_energy(species,Lx,Ly,Lz,list_2,pos)
    case default
       write(0,*) "interactions: Error. 'lattice' is not 1 or 2."
       stop 1
    end select

end function calc_energy_scratch




! Returns the energy for the specificed configuration of the specified lattice given
! that particle i has moved
function calc_energy_part_move(lattice, Lx, Ly, Lz, species, pos, pos_new, R, u, u_new, i)
    ! The lattice (1 or 2)
    integer(ik), intent(in) :: lattice
    ! The particle which has just been moved
    integer(ik), intent(in) :: i
    ! Dimensions of the (orthorhombic) supercell
    real(rk), intent(in) :: Lx, Ly, Lz
    ! Array containing the species of each particle for each lattice: e.g., species(i) is the
    ! species of particle i
    integer(ik), dimension(:), intent(in) :: species
    ! Positions (Cartesian) of the particles BEFORE particle i has been moved: e.g., pos(j,1) is 
    ! the x-coordinate of particle j, pos(j,2) is the y-coordinate, and pos(j,3) is the
    ! z-coordinate
    real(rk), dimension(:,:), intent(in) :: pos
    ! Position of particle i AFTER the particle has been moved
    real(rk), dimension(3), intent(in) :: pos_new
    ! Positions (Cartesian) of the lattice sites for the configuration: e.g., R(i,1) is
    ! the x-coordinate of the lattice site for particle i, etc.
    real(rk), intent(in), dimension(:,:) :: R
    ! Displacement vectors for the particles BEFORE particle i has been moved; e.g., u(j,1) is 
    ! the x-displacement of particle 1 from its lattice site, etc.
    real(rk) ,intent(in), dimension(:,:) :: u
    ! Displacement of particle i AFTER the particle has been moved
    real(rk), dimension(3), intent(in) :: u_new

    real(rk) :: calc_energy_part_move

    real(rk) :: sep, sep_new
    integer(ik) :: j,n

    calc_energy_part_move=0.0_rk

    n=1
    do
       select case(lattice)
       case(1)
          j=list_1(i,n)
       case(2)
          j=list_2(i,n)
       case default
          write(0,*) "interactions: Error. 'lattice' is not 1 or 2."
          stop 1
       end select

       if(j==0) then
          exit
       else
          if(j/=i) then
             sep = min_image_distance(pos(i,:),pos(j,:),Lx,Ly,Lz)
             sep_new = min_image_distance(pos_new,pos(j,:),Lx,Ly,Lz)
             calc_energy_part_move = calc_energy_part_move + &
                 !pair_potential_trunc(sep_new,species(i),species(j)) - &
                 !pair_potential_trunc(sep,species(i),species(j))
                 pair_potential(sep_new,species(i),species(j)) - &
                 pair_potential(sep,species(i),species(j))
          end if
       end if
       n=n+1
    end do

end function calc_energy_part_move




! Returns the separation between two positions within an orthorhombic
! cell with the specified dimensions
function min_image_distance(r_1, r_2, Lx, Ly, Lz)
    ! Positions of the two particles
    real(rk), dimension(3), intent(in) :: r_1, r_2
    ! Dimensions of the orthorhombic cell
    real(rk), intent(in) :: Lx, Ly, Lz

    real(rk) :: min_image_distance

    real(rk) :: xsep, ysep, zsep

    ! Calculate the x-sep
    xsep=abs(r_2(1)-r_1(1))
    xsep=xsep-Lx*floor(2.0_rk*xsep/Lx)
    ! Calculate the y-sep
    ysep=abs(r_2(2)-r_1(2))
    ysep=ysep-Ly*floor(2.0_rk*ysep/Ly)
    ! Calculate the z-sep
    zsep=abs(r_2(3)-r_1(3))
    zsep=zsep-Lz*floor(2.0_rk*zsep/Lz)
    ! Calculate the distance
    min_image_distance=sqrt(xsep*xsep+ysep*ysep+zsep*zsep)

end function min_image_distance




! The 'pure' pair potential, without truncation, between 2 particles belonging to
! the specified species and the specified separation
function pair_potential(r, species1, species2)
    ! Separation between the particles
    real(rk), intent(in) :: r
    ! Species of the two particles
    integer(ik), intent(in) :: species1, species2

    real(rk) :: pair_potential

    ! <<<<<<<<<<<< USER-DEFINED CODE >>>>>>>>>>>> 
    ! Insert code corresponding to the 'pure' pair potential, without truncation,
    ! and using the free parameters defined above.
    ! Example (for Lennard-Jones potential):
    !  pair_potential=4.0_rk*lj_epsilon*( (lj_sigma/r)**12-(lj_sigma/r)**6 )
    !
    if( r < 0.5_rk*(sigma(species1)+sigma(species2)) ) then
        pair_potential=epsilon
    else
        pair_potential=0.0_rk
    end if
    ! <<<<<<<<<<<< END OF USER-DEFINED CODE >>>>>>>>>>>> 

end function pair_potential




! Pair potential truncated at 'cutoff'. Note that there is no shifting of the potential;
! there is a discontinuity at the cut-off.
!function pair_potential_trunc(r, species1, species2)
!    ! Separation between the particles
!    real(rk), intent(in) :: r
!    ! Species of the two particles
!    integer(ik), intent(in) :: species1, species2
!
!    real(rk) :: pair_potential_trunc
!    
!    if(r<cutoff) then
!        pair_potential_trunc=pair_potential(r,species1,species2)
!    else
!        pair_potential_trunc=0.0_rk
!    end if
!
!end function pair_potential_trunc




! Returns the energy of the system calculated from scratch using the specified list
function system_energy(species, Lx, Ly, Lz, list, r)
    ! Species of all particles in the system
    integer(ik), dimension(:), intent(in) :: species
    ! Dimensions of the cell
    real(rk), intent(in) :: Lx
    real(rk), intent(in) :: Ly
    real(rk), intent(in) :: Lz
    ! List of interacting particles
    integer(ik), dimension(:,:), intent(in) :: list
    ! Positions of particles in the system
    real(rk), intent(in), dimension(:,:) :: r

    real(rk) :: system_energy

    integer(ik) :: i,j,n
    real(rk) :: sep

    system_energy=0.0_rk
    do i=1,ubound(r,1)
       n=1
       do
          j=list(i,n)
          if(j==0) then
             exit
          else
             if(j/=i) then
                sep=min_image_distance(r(i,:),r(j,:),Lx,Ly,Lz)
                !system_energy=system_energy+pair_potential_trunc(sep,species(i),species(j))
                system_energy=system_energy+pair_potential(sep,species(i),species(j))
             end if
          end if
          n=n+1
       end do
    end do
    system_energy=0.5*system_energy

end function system_energy




end module interactions_mod
