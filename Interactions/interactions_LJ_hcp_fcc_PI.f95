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
! 'interactions.f95' file for monteswitch for a path-integral lattice-switch Monte Carlo simulation of
! the fcc vs. hcp Lennard-Jones solids, in conjunction with a cut-off, fixed neighbour lists and 'tail corrections'
!
! Author: Tom L Underwood
!
! This file implements a system of disinguishable quantum particles (i.e. this file performs a path-integral
! lattice-switch Monte Carlo simulation) ineracting via Lennard-Jones interactions in hcp (lattice 1) and fcc
! (lattice 2). The interactions are not your typical Lennard-Jones interactions. Usually one truncates the
! interactions at a predetermined distance, or uses a predetermined 'list' of interacting pairs of particles.
! However, this gives errors in the ground state energies (see 'Structural Phase Behaviour Via Monte Carlo Techniques', 
! Andrew N. Jackson, PhD Thesis, University of Edinburgh, 2002). A better approach is to evaluate the 
! DIFFERENCE of the energy of the lattice under consideration relative to the ground state for the density
! under consideration, and apply the usual truncations to this difference. This is what is done here. Of 
! course, this approach necessitates a, say, fcc-specific Fortran procedure to treat the fcc lattice; while 
! in the conventional approach the Lennard-Jones procedure can be applied to any crystal. Hence this file
! should only be used in conjunction with hcp as lattice 1 and fcc with lattice 2 as the underlying
! lattices; the lattice sites should form a perfect hcp crysal for lattice 1 and a perfect fcc crystal for
! lattice 2 (though the displacements of the particles from their lattice sites are allowed to be non-zero).
!
! The potential energy here for particle positions {r} in lattice L at density V is: 
! E = Phi_{LJ,trunc}({r})-Phi_{LJ,trunc}({R_L})+ E_{GS}(L,\rho),
! where {R_L} are the lattice vectors corresponding to lattice L at density \rho, Phi_{LJ,trunc}({r}) is the 
! Lennard-Jones energy for set of particle positions {r} using truncated interactions, and E_{GS}(L,\rho) is
! the EXACT energy of lattice L at density \rho, i.e., what Phi_{LJ,trunc}({R_L}) would be if there were no
! truncations.
!
! The system treated here is comprised of N indistinguishable quantum particles over P time slices - a path
! integral representation of a system of M=N/P particles (see below). (See, e.g. 'A Guide to Monte Carlo 
! Simulations in Statistical Physics', Kurt Binder, section 8.2, for more information
! regarding path integral Monte Carlo). Particle i interacts with its corresponding particle in 'adjacent' time
! slices by a harmonic potential with spring constant 'kappa', where:
! kappa=m*P*(k_B*T)^2/hbar^2,
! with 'm' the mass of the particle, 'k_B' Boltzmann's constant and 'hbar' the reduced Plank's constant.
! Moreover particle i interacts with all particles in the same time slice via a Lennard-Jones potential in
! the expected manner.
! Here, 'P' and 'kappa' are read as input. The system of N particles specified in the 'lattices_in'
! file is then assumed to consist of P slices of M=N/P particles each, where M is the number of particles
! in the 'real' system we are interested in. Of the N particles, the 1st M correspond to slice 1, the 2nd M
! to slice 2, etc. In the code N, M and P are named 'n_part', 'n_part_slice' and 'slices', respectively.
!
! Note that all energies returned by the public procedures in this module correspond to M particles - 
! the number of particles in the 'real' system under consideration. In other words all energies are 'per slice'. By
! contrast there are N particles in the simulation, which is necessary to represent the real system of M particles
! quantum mechanically, energies returned by public procedures do NOT pertain to N particles.
!
! Fixed neighbour lists (that do not change throughout the simulation) are used for the truncated component
! of the LJ interactions for each particle.
!
! The variables for this module are imported from a file 'interactions_in', and the format of 
! this file is as follows. On the first line there are two tokens. The first is a character(len=20) variable
! (I recommend: 'lj_epsilon='); the second is the value of 'lj_epsilon' (all variables are explained in a 
! moment). The second and subsequent lines are similar, but for other input variables, which are listed
! below in the order in which they should appear in 'interactions_in' (one per line):
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
! * slices (integer(ik)) is the number of slices to consider in the calculation - refered to as 'P' above and here.
!   for brevity. As described above, the 'N' particles in the monteswitch 'lattices_in' input file are assumed to consist
!   of 'P' slices of M=N/P particles each. 'M' is the number of particles in the 'real' system to be simulated,
!   'N' being the number of particles necessary to capture the quantum aspect of this system. 'P' must be a
!   factor of 'N' otherwise an error is thrown. The 1st 'M' particles in the system comprise slice 1, the 2nd 'M'
!   comprise slice 2, etc.
!
! * kappa (real(rk)) is the harmonic constant for the inter-slice contribution to the energy - as described above.
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

    ! Variables for the potential cutoff and lists of interacting particles (intra-slice
    ! interactions only)
    real(rk) :: cutoff
    real(rk) :: list_cutoff
    integer(ik) :: list_size
    integer(ik), dimension(:,:), allocatable :: list_1
    integer(ik), dimension(:,:), allocatable :: list_2

    ! Number of time slices (refered to as 'P' in equations in comments).
    ! This must be a factor of the total number of particles in the
    ! system N: M=N/P is the number of particles which are actually simulated. It is assumed
    ! that particles 1 to M correspond to time slice 1; particles M+1 to 2M correspond to
    ! time slice 2; ...; particles (i-1)*M+1 to i*M are species i; ...; particles (P-1)*M to N
    ! correspond to time slice P
    integer(ik) :: slices
    ! Number of particles per time slice (see above) (referred to as 'M' in equations in
    ! comments)
    integer(ik) :: n_part_slice
    ! The inter-slice interaction energy between a particle and its image in another time slice
    ! is 0.5*kappa*sep2, where 'sep2' is the squared separation between the two particles
    real(rk) :: kappa
    
    ! Variables for the inter-slice interactions
    
    ! These are neighbour lists for inter-slice interactions.
    ! Each particle will have 2 neighborus: one in the time slice above and one in the
    ! time slice below that of the particle
    integer(ik), dimension(:,:), allocatable :: interslice_list_1
    integer(ik), dimension(:,:), allocatable :: interslice_list_2


    
    
    ! The below public procedures are called by 'monteswitch_mod' (in 'monteswitch_mod.f95')
    private
    public :: initialise_interactions, export_interactions, import_interactions, &
        after_accepted_part_interactions, after_accepted_vol_interactions, &
        after_accepted_lattice_interactions, after_all_interactions, &
        calc_energy_scratch, calc_energy_part_move


contains




! Initialises the variables in this module.
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
    integer(ik) :: jlow,jhigh
    
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
    ! Read the number of time slices and the harmonic constant
    read(10,*,iostat=error) string, slices
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'slices' from file 'interactions_in'"
       stop 1
    end if
    read(10,*,iostat=error) string, kappa
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'kappa' from file 'interactions_in'"
       stop 1
    end if
    close(10)

    ! Check that the number of time slices is divisible by the number of particles
    if(mod(n_part,slices)/=0) then
       write(0,*) "interactions: Error. 'slices' is not a factor of 'n_part'"
       stop 1
    end if

    ! Set n_part_slice
    n_part_slice=n_part/slices

    
    ! Initialise lists    

    if(allocated(list_1)) then
       deallocate(list_1)
    end if
    if(allocated(list_2)) then
       deallocate(list_2)
    end if
    allocate(list_1(list_size,n_part))
    allocate(list_2(list_size,n_part))

    allocate(interslice_list_1(2,n_part))
    allocate(interslice_list_2(2,n_part))
    interslice_list_1=0
    interslice_list_2=0
    
    list_1=0
    do i=1,n_part
        
       ! Intra-slice list...       
       n=1
       ! Only consider particles in the same timeslice for the neighbour list
       ! For a particle i, int((i-1)/M) gives the slice number the particle is in, with
       ! time slice labels from 0 to M-1, and M=N/P is the number of particles
       ! in each time slice and P is the number of slices. Particles in the same time slice as i have
       ! indices int((i-1)/M)*M+1 to int((i-1)/M)*M+M.
       jlow  = int((i-1)/n_part_slice)*n_part_slice+1
       jhigh = int((i-1)/n_part_slice)*n_part_slice+n_part_slice
       do j=jlow,jhigh
          if(min_image_distance(pos1(:,i),pos1(:,j),Lx1,Ly1,Lz1)<list_cutoff) then
             list_1(n,i)=j
             n=n+1
          end if
       end do

       ! Inter-slice list...
       ! Note that if there is only 1 slice then we have a classical system (no neighbours),
       ! and if there are 2 slices then we only have one neighbouring slice (i.e. the
       ! 'forward' slice). Only if there are 3 or more slices do we have both a 'forward'
       ! and a 'backward' neighbouring slice
       if(slices>1) then
           interslice_list_1(1,i) = mod(i+n_part_slice,n_part)
           if(interslice_list_1(1,i)<=0) interslice_list_1(1,i)=interslice_list_1(1,i)+n_part
           if(slices>2) then
               interslice_list_1(2,i) = mod(i-n_part_slice,n_part)
               if(interslice_list_1(2,i)<=0) interslice_list_1(2,i)=interslice_list_1(2,i)+n_part
           end if
       end if
       
    end do

    ! Do the same as above but for lattice 2
    list_2=0
    do i=1,n_part

       ! Intra-slice list...
       n=1
       jlow  = int((i-1)/n_part_slice)*n_part_slice+1
       jhigh = int((i-1)/n_part_slice)*n_part_slice+n_part_slice
       do j=jlow,jhigh
          if(min_image_distance(pos2(:,i),pos2(:,j),Lx2,Ly2,Lz2)<list_cutoff) then
             list_2(n,i)=j
             n=n+1
          end if
       end do

       ! Inter-slice list...
       if(slices>1) then
           interslice_list_2(1,i) = mod(i+n_part_slice,n_part)
           if(interslice_list_2(1,i)<=0) interslice_list_2(1,i)=interslice_list_2(1,i)+n_part
           if(slices>2) then
               interslice_list_2(2,i) = mod(i-n_part_slice,n_part)
               if(interslice_list_2(2,i)<=0) interslice_list_2(2,i)=interslice_list_2(2,i)+n_part
           end if
       end if
       
    end do
    
end subroutine initialise_interactions




! Export all of the variables in this module - for the purposes of checkpointing.
subroutine export_interactions()

    write(10,*) "lj_epsilon= ",lj_epsilon
    write(10,*) "lj_sigma= ",lj_sigma

    write(10,*) "cutoff= ",cutoff
    write(10,*) "list_cutoff= ",list_cutoff
    write(10,*) "list_size= ",list_size
    write(10,*) "slices= ",slices
    write(10,*) "kappa= ",kappa
    write(10,*) "n_part_slice= ",n_part_slice
    write(10,*) "list_1= ",list_1
    write(10,*) "list_2= ",list_2
    write(10,*) "interslice_list_1= ",interslice_list_1
    write(10,*) "interslice_list_2= ",interslice_list_2
    
end subroutine export_interactions




! Import all variables in this module from file(s) - for the purposes of checkpointing. The format 
! read from the file(s) should correspond to that output by 'export_interactions' above.
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
    read(10,*,iostat=error) string, slices
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'slices' from unit 10"
       stop 1
    end if
    read(10,*,iostat=error) string, kappa
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'kappa' from unit 10"
       stop 1
    end if
    read(10,*,iostat=error) string, n_part_slice
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'n_part_slice' from unit 10"
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

    if(allocated(interslice_list_1)) then
       deallocate(interslice_list_1)
    end if
    if(allocated(interslice_list_2)) then
       deallocate(interslice_list_2)
    end if
    allocate(interslice_list_1(2,n_part))
    allocate(interslice_list_2(2,n_part))

    read(10,*,iostat=error) string, interslice_list_1
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'interslice_list_1' from unit 10"
       stop 1
    end if
    read(10,*,iostat=error) string, interslice_list_2
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'interslice_list_2' from unit 10"
       stop 1
    end if

end subroutine import_interactions




! Performs any tasks required for the variables in this module, after a particle move for particle i is accepted
! Empty here - there are no such tasks
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
! Empty here - there are no such tasks
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
! Empty here - there are no such tasks
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
! Empty here - there are no such tasks
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




! Returns the energy for the specificed configuration of the specified lattice. 
! Note that the energy returned is 'per slice', i.e. the energy pertains to 'n_part_slice' particles.
! Moreover the contribution to this due to the intra-slice (Lennard-Jones) energy includes 'tail
! corrections'.
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
    
    select case(lattice)
    case(1)

       ! Calculate the ground state energy (for 0 particle displacements) for the volume of the input system (corresponding to Lx,
       ! Ly, Lz, r) (hcp)
       ! Note that this is the ground state energy per 'n_part_slice' particles, which is the number of particles in the 'real'
       ! system we are considering (see preamble to this module)
       E_gs = 2*n_part_slice*lj_epsilon*( (lj_sigma**3*n_part_slice/(Lx*Ly*Lz*root2))**4*A12_hcp &
              -(lj_sigma**3*n_part_slice/(Lx*Ly*Lz*root2))**2*A6_hcp ) 
       ! Calculate the difference relative to this ground state using truncation (1st 2 terms), and add it to the
       ! ground state for the input volume. Note that the 2nd term is the energy of the ground state for the input volume
       ! with a truncated potential. Note that 'system_potential_energy' returns the energy using a truncated potential
       ! PER TIME SLICE, i.e. the energy corresponds to 'n_part_slice' particles, which is the number in the 'real' system 
       ! under consideration). Note also that 'E_gs' amounts to 'tail corrections'
       calc_energy_scratch = system_potential_energy(species,Lx,Ly,Lz,list_1,pos) &
                           - system_potential_energy(species,Lx,Ly,Lz,list_1,R) + E_gs
        
       ! Add the interslice energy
       calc_energy_scratch = calc_energy_scratch + system_interslice_energy(Lx,Ly,Lz,interslice_list_1,pos)
       
    case(2)

       ! As above, but for fcc
       E_gs = 2*n_part_slice*lj_epsilon*( (lj_sigma**3*n_part_slice/(Lx*Ly*Lz*root2))**4*A12_fcc &
              -(lj_sigma**3*n_part_slice/(Lx*Ly*Lz*root2))**2*A6_fcc ) 
       calc_energy_scratch = system_potential_energy(species,Lx,Ly,Lz,list_2,pos) &
                           - system_potential_energy(species,Lx,Ly,Lz,list_2,R) + E_gs

       ! Add the interslice energy 
       calc_energy_scratch = calc_energy_scratch  + system_interslice_energy(Lx,Ly,Lz,interslice_list_2,pos)
       
    case default
       write(0,*) "interactions: Error. 'lattice' is not 1 or 2."
       stop 1
    end select

end function calc_energy_scratch




! Returns the energy for the specificed configuration of the specified lattice given
! that particle i has moved - includes intra- and inter-slice contributions.
! Note that the energy returned is 'per slice', i.e. the energy pertains to 'n_part_slice' particles
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

    ! Calculate the change in potential energy of the system (per 'n_part' particles)
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
             calc_energy_part_move = calc_energy_part_move + &
                 pair_potential_trunc(sep_new2,species(i),species(j)) - &
                 pair_potential_trunc(sep2,species(i),species(j))
          end if
       end if
       n=n+1
    end do

    ! Rescale to get the potential energy per 'n_part_slice' particles
    calc_energy_part_move = calc_energy_part_move / slices

    
    ! Calculate the change in the inter-slice energy... (only necessary if there are >1 slices)

    if(slices>1) then

        ! If we have >=2 slices then consider the interactions with the analogous particle in the 
        ! forward slice...

        ! Forwards slice
        select case(lattice)
        case(1)
            j=interslice_list_1(1,i)
        case(2)
            j=interslice_list_2(1,i)
        end select
        sep2=min_image_distance2(pos(:,i),pos(:,j),Lx,Ly,Lz)
        sep_new2 = min_image_distance2(pos_new,pos(:,j),Lx,Ly,Lz)
        calc_energy_part_move = calc_energy_part_move + &
            0.5_rk * kappa * (sep_new2 - sep2)

    end if
        
    if(slices>2) then         

        ! If we have >=3 slices then also consider the interactions with the analogous particle in the
        ! backwards slice...

        ! Backwards slice
        select case(lattice)
        case(1)
            j=interslice_list_1(2,i)
        case(2)
            j=interslice_list_2(2,i)
        end select
        sep2=min_image_distance2(pos(:,i),pos(:,j),Lx,Ly,Lz)
        sep_new2 = min_image_distance2(pos_new,pos(:,j),Lx,Ly,Lz)
        calc_energy_part_move = calc_energy_part_move + &
            0.5_rk * kappa * (sep_new2 - sep2)

    end if

    
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
    if(xsep>0.5_rk*Lx) xsep=Lx-xsep
    ! Calculate the y-sep
    ysep=abs(r_2(2)-r_1(2))
    if(ysep>0.5_rk*Ly) ysep=Ly-ysep
    ! Calculate the z-sep
    zsep=abs(r_2(3)-r_1(3))
    if(zsep>0.5_rk*Lz) zsep=Lz-zsep
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
    if(xsep>0.5_rk*Lx) xsep=Lx-xsep
    ! Calculate the y-sep
    ysep=abs(r_2(2)-r_1(2))
    if(ysep>0.5_rk*Ly) ysep=Ly-ysep
    ! Calculate the z-sep
    zsep=abs(r_2(3)-r_1(3))
    if(zsep>0.5_rk*Lz) zsep=Lz-zsep

    min_image_distance2=xsep*xsep+ysep*ysep+zsep*zsep

end function min_image_distance2




! The 'pure' Lennard-Jones pair potential, without truncation, between 2 particles belonging to
! the specified species and the specified squared separation
function pair_potential(r2, species1, species2)
    ! Squared separation between the particles
    real(rk), intent(in) :: r2
    ! Species of the two particles
    integer(ik), intent(in) :: species1, species2

    real(rk) :: pair_potential
    pair_potential=4.0_rk*lj_epsilon*( (lj_sigma**12/r2**6) - (lj_sigma**6/r2**3) )

end function pair_potential




! The Lennard-Jones potential truncated at 'cutoff'. Note that there is no shifting of the potential;
! there is a discontinuity at the cut-off.
function pair_potential_trunc(r2, species1, species2)
    ! Squared separation between the particles
    real(rk), intent(in) :: r2
    ! Species of the two particles
    integer(ik), intent(in) :: species1, species2

    real(rk) :: pair_potential_trunc
    
    if(r2<cutoff*cutoff) then
        pair_potential_trunc=pair_potential(r2,species1,species2)
    else
        pair_potential_trunc=0.0_rk
    end if

end function pair_potential_trunc




! Returns the intra-slice potential energy of the system PER TIME SLICE (i.e. 'n_part_slice' particles) calculated from 
! scratch using the specified list. Note that this energy does NOT include the tail corrections, which are added
! in 'calc_energy_scratch'. 
function system_potential_energy(species, Lx, Ly, Lz, list, r)
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

    real(rk) :: system_potential_energy
    
    integer(ik) :: i,j,n
    real(rk) :: sep2
    
    system_potential_energy=0.0_rk
    do i=1,ubound(r,2)
       n=1
       do
          j=list(n,i)
          if(j==0) then
             exit
          else
             if(j/=i) then
                sep2=min_image_distance2(r(:,i),r(:,j),Lx,Ly,Lz)
                system_potential_energy=system_potential_energy+pair_potential_trunc(sep2,species(i),species(j))
             end if
          end if
          n=n+1
       end do
    end do

    ! Multiply by 0.5 to eliminate double counting and divide by 'slices' to get the energy per 'n_part_slice' particles
    ! (the above is the potential energy per 'n_part' particles)
    system_potential_energy = 0.5_rk * system_potential_energy / slices

end function system_potential_energy




! Returns the inter-slice energy of the system calculated from scratch using the specified list
function system_interslice_energy(Lx, Ly, Lz, list, r)
    ! Dimensions of the cell
    real(rk), intent(in) :: Lx
    real(rk), intent(in) :: Ly
    real(rk), intent(in) :: Lz
    ! List of interacting particles (for inter-slice energy calculations) - dimension(2,n_part)
    integer(ik), dimension(:,:), intent(in) :: list
    ! Positions of particles in the system
    real(rk), intent(in), dimension(:,:) :: r

    real(rk) :: system_interslice_energy

    integer(ik) :: i,j,n
    real(rk) :: sep2

    system_interslice_energy=0.0_rk

    if(slices==1) then

        ! If we have 1 slice then we have no inter-slice energy
        
        return

    else if(slices==2) then

        ! If we have 2 slices then consider the interactions with the forward slice...

        do i=1,ubound(r,2)

            ! Forwards slice
            j=list(1,i)
            sep2=min_image_distance2(r(:,i),r(:,j),Lx,Ly,Lz)
            system_interslice_energy = system_interslice_energy + 0.5_rk * kappa * sep2

        end do

        
        ! A faster approach? Considers bonds between slices instead of double counting. A similar
        ! approach could be used for 3 or more slices below
        ! 
        !         ! If we have 2 slices then consider the interactions with each particle in slice 1 with its analogue
        !         ! in slice 2 (which is the 'forward slice' of slice 1). (Note that one could equivalently consider
        !         ! the interactions with each particle in slice 2, whose forward slice is slice 1).
        !         ! Note that particles with index 1 to n_part_slice belong to slice 1
        ! 
        !         do i=1,n_part_slice
        ! 
        !             ! j here is the index of the particle in slice 2
        !             j=list(1,i)
        !             sep2=min_image_distance2(r(:,i),r(:,j),Lx,Ly,Lz)
        !             system_interslice_energy = system_interslice_energy + kappa * sep2
        ! 
        !         end do



    else

        ! If we have 3 or more slices then consider the interactions with the forwards and backwards slices...
        
        do i=1,ubound(r,2)

            ! Forwards slice
            j=list(1,i)
            sep2=min_image_distance2(r(:,i),r(:,j),Lx,Ly,Lz)
            system_interslice_energy = system_interslice_energy + 0.5_rk * kappa * sep2

            ! Backwards slice
            j=list(2,i)
            sep2=min_image_distance2(r(:,i),r(:,j),Lx,Ly,Lz)
            system_interslice_energy = system_interslice_energy + 0.5_rk * kappa * sep2
            
       end do

    end if
   
    system_interslice_energy = 0.5_rk * system_interslice_energy

end function system_interslice_energy




end module interactions_mod
