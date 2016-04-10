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
! Minimal template for 'interactions.f95' file for monteswitch
! Author: Tom L Underwood

module interactions_mod


    ! The module 'kinds_mod' (in 'kinds_mod.f95') contains the real and integer kinds
    ! which monteswitch uses: real(kind=rk) and integer(kind=ik). Please use these
    ! kinds throughout this module to allow portability.
    use kinds_mod
    
    implicit none

    ! The below public procedures are called by 'monteswitch_mod' (in 'monteswitch_mod.f95'),
    ! and must be 'filled in' by the user.
    private
    public :: initialise_interactions, export_interactions, import_interactions, &
        after_accepted_part_interactions, after_accepted_vol_interactions, &
        after_accepted_lattice_interactions, after_all_interactions, &
        calc_energy_scratch, calc_energy_part_move


contains


! Initialises the variables in this module. The initial configurations for both
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

    return

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

    return
    
end subroutine export_interactions




! Import all variables in this module from file(s) - for the purposes of checkpointing. The format 
! read from the file(s) should correspond to that output by 'export_interactions' above. If one is
! importing from within a 'state' file as described in that procedure, then use unit 10, but
! do not open or close that unit!
subroutine import_interactions(Lx1, Ly1, Lz1, species1, pos1, R1, Lx2, Ly2, Lz2, species2, pos2, u, R2)
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

    return

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

    return

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

    return

end function calc_energy_part_move




end module interactions_mod
