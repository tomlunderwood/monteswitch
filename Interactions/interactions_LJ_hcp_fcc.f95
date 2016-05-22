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
! 'interactions.f95' file for special LJ interactions for hcp-fcc system
!
! Author: Tom L Underwood

! This version of the file corresponds to Lennard-Jones interactions in hcp (lattice 1) and fcc (lattice 2).
! The interactions are not your typical Lennard-Jones interactions. Usually one truncates the interactions at
! a predetermined distance, or uses a predetermined 'list' of interacting pairs of particles. However, this
! gives errors in the ground state energies (see 'Structural Phase Behaviour Via Monte Carlo Techniques', 
! Andrew N. Jackson, PhD Thesis, University of Edinburgh, 2002). A better approach is to evaluate the 
! DIFFERENCE of the energy of the lattice under consideration relative to the ground state for the density
! under consideration, and apply the usual truncations to this difference. This is what is done here. Of 
! course, this approach necessitates a, say, fcc-specific Fortran procedure to treat the fcc lattice; while 
! in the conventional approach the Lennard-Jones procedure can be applied to any crystal. Hence this file
! should only be used in conjunction with hcp as lattice 1 and fcc with lattice 2 as the underlying
! lattices; the lattice sites should form a perfect hcp crysal for lattice 1 and a perfect fcc crystal for
! lattice 2 (though the displacements of the particles from their lattice sites are allowed to be non-zero).
!
! The energy here for particle positions {r} in lattice L at density V is: 
! E = Phi_{LJ,trunc}({r})-Phi_{LJ,trunc}({R_L})+ E_{GS}(L,\rho),
! where {R_L} are the lattice vectors corresponding to lattice L at density \rho, Phi_{LJ,trunc}({r}) is the 
! Lennard-Jones energy for set of particle positions {r} using truncated interactions, and E_{GS}(L,\rho) is
! the EXACT energy of lattice L at density \rho, i.e., what Phi_{LJ,trunc}({R_L}) would be if there were no
! truncations.
!
! The variables for this module are imported from a file 'interactions_in', and the format of this file is
! as follows. On the first line there are two tokens. The first is a character(len=20) variable (I recommend:
! 'lj_epsilon='); the second is the value of 'lj_epsilon' (all variables are explained in a moment). The second 
! line is similar, but for 'lj_sigma'. The third line is for 'cutoff', the fourth is for 'list_cutoff',
! and the fifth is for 'list_size'. The variables are as follows:
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

    ! Lennard-Jones parameters
    real(rk) :: lj_epsilon
    real(rk) :: lj_sigma
    ! Variables for the potential cutoff and lists of interacting particles
    real(rk) :: cutoff
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
    ! Initial positions (Cartesian) of the particles for lattices 1 and 2: e.g., pos1(1,i) is
    ! the x-coordinate of particle i in lattice 1, pos1(2,i) is the y-coordinate, and pos1(3,i)
    ! is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2

    character(len=20) string
    integer(ik) :: error
    integer(ik) :: i,j,n,n_part

    n_part = size(pos1,2)

    ! Read the interactions variables from the file 'interactions_in'
    open(unit=10,file="interactions_in",iostat=error,status="old")
    if(error/=0) then
       write(0,*) "interactions: Error. Problem opening file 'interactions_in'"
       stop 1
    end if

    read(10,*,iostat=error) string, lj_epsilon
    if(error/=0) then
        write(0,*) "interactions: Error. Problem reading 'lj_epsilon' from file 'interactions_in'"
        stop 1
    end if
    read(10,*,iostat=error) string, lj_sigma
    if(error/=0) then
        write(0,*) "interactions: Error. Problem reading 'lj_sigma' from file 'interactions_in'"
        stop 1
    end if

    read(10,*,iostat=error) string, cutoff
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'cutoff' from file 'interactions_in'"
       stop 1
    end if
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
    allocate(list_1(list_size,n_part))
    allocate(list_2(list_size,n_part))

    list_1=0
    do i=1,n_part
       n=1
       do j=1,n_part
          if(min_image_distance(pos1(:,i),pos1(:,j),Lx1,Ly1,Lz1)<list_cutoff) then
             list_1(n,i)=j
             n=n+1
          end if
       end do
    end do
    list_2=0
    do i=1,n_part
       n=1
       do j=1,n_part
          if(min_image_distance(pos2(:,i),pos2(:,j),Lx2,Ly2,Lz2)<list_cutoff) then
             list_2(n,i)=j
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

    write(10,*) "lj_epsilon= ",lj_epsilon
    write(10,*) "lj_sigma= ",lj_sigma

    write(10,*) "cutoff= ",cutoff
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
    ! Positions (Cartesian) of the particles for lattices 1 and 2: e.g., pos1(1,i) is
    ! the x-coordinate of particle i in lattice 1, pos1(2,i) is the y-coordinate, and pos1(3,i)
    ! is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2
    ! Positions (Cartesian) of the lattice sites for lattices 1 and 2: e.g., R1(1,i) is
    ! the x-coordinate of the lattice site for particle i in lattice 1, ..., R2(1,i) is
    ! the x-coordinate of the lattice site for particle i in lattice 2
    real(rk), intent(in), dimension(:,:) :: R1, R2
    ! Displacements of the particles: e.g., u(1,i) is the x-coordinate of particle i in lattice 1
    ! etc.
    real(rk), intent(in), dimension(:,:) :: u

    character(len=20) string
    integer(ik) :: error, n_part
    
    read(10,*,iostat=error) string, lj_epsilon
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'lj_epsilon' from unit 10"
       stop 1
    end if
    read(10,*,iostat=error) string, lj_sigma
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'lj_sigma' from unit 10"
       stop 1
    end if

    read(10,*,iostat=error) string, cutoff
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'cutoff' from unit 10"
       stop 1
    end if
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

    n_part = size(pos1,2)
    if(allocated(list_1)) then
       deallocate(list_1)
    end if
    if(allocated(list_2)) then
       deallocate(list_2)
    end if
    allocate(list_1(list_size,n_part))
    allocate(list_2(list_size,n_part))

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
! Leave empty if there are no such tasks
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
    ! accepted): e.g., pos1(1,i) is the x-coordinate of particle i in lattice 1, pos1(2,i) is the 
    ! y-coordinate, and pos1(3,i) is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2
    ! Positions (Cartesian) of the lattice sites for lattices 1 and 2: e.g., R1(1,i) is
    ! the x-coordinate of the lattice site for particle i in lattice 1, ..., R2(1,i) is
    ! the x-coordinate of the lattice site for particle i in lattice 2
    real(rk), intent(in), dimension(:,:) :: R1, R2
    ! Current displacement vectors for the particles (after the move has been accepted); e.g., 
    ! u(1,i) is the x-displacement of particle 1 from its lattice site, etc.
    real(rk) ,intent(in), dimension(:,:) :: u

    return

end subroutine after_accepted_part_interactions




! Performs any tasks required for the variables in this module, after a volume move is accepted
! Leave empty if there are no such tasks
subroutine after_accepted_vol_interactions(Lx1, Ly1, Lz1, species1, pos1, R1, Lx2, Ly2, Lz2, species2, pos2, R2, u)
    ! Current dimensions of the initial (orthorhombic) supercell for lattices 1 and 2
    ! in each Cartesian dimension (after the move has been accepted)
    real(rk), intent(in) :: Lx1, Ly1, Lz1, Lx2, Ly2, Lz2
    ! Array containing the species of each particle for each lattice: e.g., species1(i) is the
    ! species of particle i in lattice 1
    integer(ik), intent(in), dimension(:) :: species1, species2
    ! Current positions (Cartesian) of the particles for lattices 1 and 2 (after the move has been
    ! accepted): e.g., pos1(1,i) is the x-coordinate of particle i in lattice 1, pos1(2,i) is the 
    ! y-coordinate, and pos1(3,i) is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2
    ! Positions (Cartesian) of the lattice sites for lattices 1 and 2: e.g., R1(1,i) is
    ! the x-coordinate of the lattice site for particle i in lattice 1, ..., R2(1,i) is
    ! the x-coordinate of the lattice site for particle i in lattice 2
    real(rk), intent(in), dimension(:,:) :: R1, R2
    ! Current displacement vectors for the particles (after the move has been accepted); e.g., 
    ! u(1,i) is the x-displacement of particle 1 from its lattice site, etc.
    real(rk) ,intent(in), dimension(:,:) :: u


    return

end subroutine after_accepted_vol_interactions




! Performs any tasks required for the variables in this module, after a lattice move is accepted
! Leave empty if there are no such tasks
subroutine after_accepted_lattice_interactions(Lx1, Ly1, Lz1, species1, pos1, R1, Lx2, Ly2, Lz2, species2, pos2, R2, u)
    ! Current dimensions of the initial (orthorhombic) supercell for lattices 1 and 2
    ! in each Cartesian dimension (after the move has been accepted)
    real(rk), intent(in) :: Lx1, Ly1, Lz1, Lx2, Ly2, Lz2
    ! Array containing the species of each particle for each lattice: e.g., species1(i) is the
    ! species of particle i in lattice 1
    integer(ik), intent(in), dimension(:) :: species1, species2
    ! Current positions (Cartesian) of the particles for lattices 1 and 2 (after the move has been
    ! accepted): e.g., pos1(1,i) is the x-coordinate of particle i in lattice 1, pos1(2,i) is the 
    ! y-coordinate, and pos1(3,i) is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2
    ! Positions (Cartesian) of the lattice sites for lattices 1 and 2: e.g., R1(1,i) is
    ! the x-coordinate of the lattice site for particle i in lattice 1, ..., R2(1,i) is
    ! the x-coordinate of the lattice site for particle i in lattice 2
    real(rk), intent(in), dimension(:,:) :: R1, R2
    ! Current displacement vectors for the particles (after the move has been accepted); e.g., 
    ! u(1,i) is the x-displacement of particle 1 from its lattice site, etc.
    real(rk) ,intent(in), dimension(:,:) :: u


    return

end subroutine after_accepted_lattice_interactions




! Performs any tasks required for the variables in this module after ALL moves (accepted or rejected)
! Leave empty if there are no such tasks
subroutine after_all_interactions(Lx1, Ly1, Lz1, species1, pos1, R1, Lx2, Ly2, Lz2, species2, pos2, R2, u)
    ! Current dimensions of the initial (orthorhombic) supercell for lattices 1 and 2
    ! in each Cartesian dimension
    real(rk), intent(in) :: Lx1, Ly1, Lz1, Lx2, Ly2, Lz2
    ! Array containing the species of each particle for each lattice: e.g., species1(i) is the
    ! species of particle i in lattice 1
    integer(ik), intent(in), dimension(:) :: species1, species2
    ! Current positions (Cartesian) of the particles for lattices 1 and 2: e.g., pos1(1,i) is 
    ! the x-coordinate of particle i in lattice 1, pos1(2,i) is the y-coordinate, and pos1(3,i) 
    ! is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2
    ! Positions (Cartesian) of the lattice sites for lattices 1 and 2: e.g., R1(1,i) is
    ! the x-coordinate of the lattice site for particle i in lattice 1, ..., R2(1,i) is
    ! the x-coordinate of the lattice site for particle i in lattice 2
    real(rk), intent(in), dimension(:,:) :: R1, R2
    ! Current displacement vectors for the particles (after the move has been accepted); e.g., 
    ! u(1,i) is the x-displacement of particle 1 from its lattice site, etc.
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
    ! Positions (Cartesian) of the particles: e.g., pos(1,i) is the x-coordinate of particle 
    ! i, pos1(2,i) is the y-coordinate, and pos1(3,i) is the z-coordinate
    real(rk), dimension(:,:), intent(in) :: pos
    ! Positions (Cartesian) of the lattice sites for the configuration: e.g., R(1,i) is
    ! the x-coordinate of the lattice site for particle i, etc.
    real(rk), intent(in), dimension(:,:) :: R
    ! Displacement vectors for the particles; e.g., u(1,i) is the x-displacement of 
    ! particle 1 from its lattice site, etc.
    real(rk) ,intent(in), dimension(:,:) :: u

    real(rk) :: calc_energy_scratch

    ! Constants for calculating the exact ground state for the hcp and fcc lattices
    ! at a given density (see Table 6.1 of Andrew Jackson's thesis mentioned above)
    real(rk), parameter :: A12_hcp=12.132293768
    real(rk), parameter :: A6_hcp=14.454897093
    real(rk), parameter :: A12_fcc=12.131880196
    real(rk), parameter :: A6_fcc=14.453920885
    real(rk) :: root2=sqrt(2.0_rk)
    real(rk) :: E_gs
    integer(ik) :: n_part

    n_part = size(pos,2)

    select case(lattice)
    case(1)
       ! Calculate the ground state energy (for 0 particle displacements) for the volume of the input system (corresponding to Lx,
       ! Ly, Lz, r) (hcp)
       E_gs=2*n_part*lj_epsilon*( (lj_sigma**3*n_part/(Lx*Ly*Lz*root2))**4*A12_hcp &
            -(lj_sigma**3*n_part/(Lx*Ly*Lz*root2))**2*A6_hcp )
       ! Calculate the difference relative to this ground state using truncation (1st 2 terms), and add it to the
       ! ground state for the input volume. Note that the 2nd term is the energy of the ground state for the input volume
       ! with a truncated potential.
       calc_energy_scratch = lj_energy_trunc_list_2(lj_epsilon,lj_sigma,cutoff,Lx,Ly,Lz,list_1,pos) &
           -lj_energy_trunc_list_2(lj_epsilon,lj_sigma,cutoff,Lx,Ly,Lz,list_1,R) + E_gs
    case(2)
       ! As above, but for fcc
       E_gs=2*n_part*lj_epsilon*( (lj_sigma**3*n_part/(Lx*Ly*Lz*root2))**4*A12_fcc &
            -(lj_sigma**3*n_part/(Lx*Ly*Lz*root2))**2*A6_fcc )
       calc_energy_scratch = lj_energy_trunc_list_2(lj_epsilon,lj_sigma,cutoff,Lx,Ly,Lz,list_2,pos) &
           -lj_energy_trunc_list_2(lj_epsilon,lj_sigma,cutoff,Lx,Ly,Lz,list_2,R) + E_gs
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
    ! Positions (Cartesian) of the particles BEFORE particle i has been moved: e.g., pos(1,j) is 
    ! the x-coordinate of particle j, pos(2,j) is the y-coordinate, and pos(3,j) is the
    ! z-coordinate
    real(rk), dimension(:,:), intent(in) :: pos
    ! Position of particle i AFTER the particle has been moved
    real(rk), dimension(3), intent(in) :: pos_new
    ! Positions (Cartesian) of the lattice sites for the configuration: e.g., R(1,i) is
    ! the x-coordinate of the lattice site for particle i, etc.
    real(rk), intent(in), dimension(:,:) :: R
    ! Displacement vectors for the particles BEFORE particle i has been moved; e.g., u(1,j) is 
    ! the x-displacement of particle 1 from its lattice site, etc.
    real(rk) ,intent(in), dimension(:,:) :: u
    ! Displacement of particle i AFTER the particle has been moved
    real(rk), dimension(3), intent(in) :: u_new

    real(rk) :: calc_energy_part_move

    ! Squared separations between particles
    real(rk) :: sep2, sep_new2
    integer(ik) :: j,n

    calc_energy_part_move=0.0_rk
    n=1
    do
       select case(lattice)
       case(1)
          j=list_1(n,i)
       case(2)
          j=list_2(n,i)
       case default
          write(0,*) "interactions: Error. 'lattice' is not 1 or 2."
          stop 1
       end select

       if(j==0) then
          exit
       else
          if(j/=i) then
             sep2 = min_image_distance2(pos(:,i),pos(:,j),Lx,Ly,Lz)
             sep_new2 = min_image_distance2(pos_new,pos(:,j),Lx,Ly,Lz)
             calc_energy_part_move = calc_energy_part_move + lj_pot_trunc(lj_epsilon,lj_sigma,cutoff,sep_new2) &
                  - lj_pot_trunc(lj_epsilon,lj_sigma,cutoff,sep2)
          end if
       end if
       n=n+1
    end do

end function calc_energy_part_move




! Returns the separation between two positions within an orthorhombic
! cell with the specified dimensions
function min_image_distance(r_1,r_2,Lx,Ly,Lz)
    real(rk), dimension(3), intent(in) :: r_1, r_2
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




! Returns the squared separation between two positions within an orthorhombic
! cell with the specified dimensions
function min_image_distance2(r_1, r_2, Lx, Ly, Lz)
    ! Positions of the two particles
    real(rk), dimension(3), intent(in) :: r_1, r_2
    ! Dimensions of the orthorhombic cell
    real(rk), intent(in) :: Lx, Ly, Lz

    real(rk) :: min_image_distance2

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

    min_image_distance2=xsep*xsep+ysep*ysep+zsep*zsep

end function min_image_distance2




! This function returns the value of the Lennard-Jones potential at the specified squared
! distance, given the specified parametrisation of the potential.
function lj_pot(epsilon,sigma,r2)
    real(rk), intent(in) :: epsilon
    real(rk), intent(in) :: sigma
    real(rk), intent(in) :: r2
    real(rk) :: lj_pot

    lj_pot = ( sigma*sigma / r2 )**3
    lj_pot = 4.0_rk*epsilon * ( lj_pot - 1 ) * lj_pot

end function lj_pot




! This function returns the value of a truncated Lennard-Jones potential at the specified squared
! distance, given the specified parametrisation of the potential. The potential is truncated
! at distance 'cutoff'.
function lj_pot_trunc(epsilon,sigma,cutoff,r2)
    real(rk), intent(in) :: epsilon
    real(rk), intent(in) :: sigma
    real(rk), intent(in) :: cutoff
    real(rk), intent(in) :: r2
    real(rk) :: lj_pot_trunc
    if(r2<cutoff*cutoff) then
       lj_pot_trunc=lj_pot(epsilon,sigma,r2)
    else
       lj_pot_trunc=0.0_rk
    end if
end function lj_pot_trunc




! This function returns the Lennard-Jones energy for the specified system, list and parametrisation
! of the LJ potential (which is truncated)
function lj_energy_trunc_list_2(epsilon,sigma,cutoff,Lx,Ly,Lz,list,r)
    real(rk), intent(in) :: epsilon
    real(rk), intent(in) :: sigma
    real(rk), intent(in) :: cutoff
    real(rk), intent(in) :: Lx
    real(rk), intent(in) :: Ly
    real(rk), intent(in) :: Lz
    integer(ik), dimension(:,:), intent(in) :: list
    real(rk), intent(in), dimension(:,:) :: r
    real(rk) :: lj_energy_trunc_list_2
    integer(ik) :: i,j,n
    real(rk) :: sep2
    lj_energy_trunc_list_2=0.0_rk
    do i=1,ubound(r,2)
       n=1
       do
          j=list(n,i)
          if(j==0) then
             exit
          else
             if(j/=i) then
                sep2=min_image_distance2(r(:,i),r(:,j),Lx,Ly,Lz)
                lj_energy_trunc_list_2=lj_energy_trunc_list_2+lj_pot_trunc(epsilon,sigma,cutoff,sep2)
             end if
          end if
          n=n+1
       end do
    end do
    lj_energy_trunc_list_2=0.5*lj_energy_trunc_list_2
end function lj_energy_trunc_list_2




end module interactions_mod
