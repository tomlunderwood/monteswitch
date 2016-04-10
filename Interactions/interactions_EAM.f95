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
! 'interactions.f95' file for EAM metal potentials (single-species '.eam.alloy' format) 
!
! Author: Tom L Underwood
!
! This file implements the embedded atom model (EAM) similarly to LAMMPS, as described at
! http://lammps.sandia.gov/doc/pair_eam.html. The input file which specifies the potential to be used by monteswitch
! must be named 'interactions_in', and must be in the DYNAMO/LAMMPS 'setfl' format; a file with this format is often
! indicated with the suffix '.eam.alloy'. A description of the 'setfl' format can be found at 
! http://lammps.sandia.gov/doc/pair_eam.html. However, note that this file ONLY WORKS FOR SINGLE-SPECIES
! SYSTEMS: the .eam.alloy file must not correspond to a multicomponent system - an error is returned
! if this is the case. Note also this file implements the EAM potential by using linear interpolation, using the
! tabulations of the various functions in the input file (the embedding function, the density function, 
! and the pair-potential) as a basis. 
!
! At initialisation this module creates three files: 'F.dat', 'rho.dat', and 'rphi.dat', which
! correspond to the embedding function, density function, and pair potential (multiplied by separation) read
! from the input file.
!
! NB: It seems that 'setfl' format files often have cut-offs which are very slightly
! higher than that suggested by the values of 'Nr' and 'dr. The aforementioned LAMMPS documentation specifies that
! elements 'i' of rho' and 'rphi' correspond to a separation of 'r=(i-1)*dr' - as opposed to
! 'r=i*dr'. Hence the maximum supported separation between particles is '(Nr-1)*dr', as opposed to 'Nr*dr'. 
! However, people often set the cut-off to the latter value, as opposed to the former. Here, to account for this, 
! if the cut-off is greater than '(Nr-1)*dr', then it is ammeded at initialisation to be '(Nr-1)*dr'.
!
! This file was adapted from the file 'interactions_TEMPLATE_minimal.f95' - some comments are inherited from that
! file.

module interactions_mod


    ! The module 'kinds_mod' (in 'kinds_mod.f95') contains the real and integer kinds
    ! which monteswitch uses: real(kind=rk) and integer(kind=ik). Please use these
    ! kinds throughout this module to allow portability.
    use kinds_mod
    
    implicit none

    ! Number of points on the density grid
    integer(ik) :: Nrho
    ! Spacing of points on the density grid
    real(rk) :: drho
    ! Number of points on the separation grid
    integer(ik) :: Nr
    ! Spacing of points on the separation grid
    real(rk) :: dr
    ! The interaction cut-off: pairs of particles beyond this separation do not interact.
    real(rk) :: cutoff

    ! Array defining the embedding function: 'F(i)' is the embedding function at density
    ! '(i-1)*drho'.
    real(rk), dimension(:), allocatable :: F
    ! Array defining the density function: 'rho(i)' is the density for inter-particle separation
    ! '(i-1)*dr'.
    real(rk), dimension(:), allocatable :: rho
    ! Array defining the pair-potential function: 'rphi(i)' is the pair-potential, multiplied by
    ! the separation, for inter-particle separation '(i-1)*dr'.
    real(rk), dimension(:), allocatable :: rphi
    ! Array defining containing the density - with regards to the embedding function - associated with
    ! each particle in the current microstate of the system for each latice. 
    ! 'part_rho(lattice,i)' is the density for particle 'i' in lattice <cdoe>lattice'. 
    ! This is used to speed up the calculation of the change in embedding energy as a result of a particle move. 
    real(rk), dimension(:,:), allocatable :: part_rho
    ! Like 'part_rho', but a 'buffer' value. Every particle or volume move, the energy of the trial
    ! microstate is evaluated by calling the 'calc_energy_scratch' or 'calc_energy_part_move'
    ! procedure. In each of these procedures, the densities for the trial state are stored in 'part_rho_buffer'.
    ! If the move is accepted, i.e., the trial microstate becomes the actual microstate of the system, 
    ! 'part_rho_buffer' is copied to 'part_rho'. If the move is rejected, then 
    ! 'part_rho_buffer' is not copied to 'part_rho'. Thus 'part_rho' always reflects the 
    ! actual microstate.
    real(rk), dimension(:,:), allocatable :: part_rho_buffer


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

    character(*),parameter :: filename="interactions_in"
    character(len=200) :: line
    integer(ik) :: error
    integer(ik) :: i,j,n

    integer(ik) :: indx, prev_indx
    integer(ik) :: start_of_word

    integer(ik) :: n_part

    n_part=size(pos1,1)

    ! Open the file and import the variables required for initialisation
    open(unit=10,file=filename,iostat=error,status="old")
    if(error/=0) then
       write(0,*) "interactions: Error. Problem opening file '",filename,"'"
       stop 1
    end if

    ! Read 3 comment lines
    read(10,*,iostat=error) line
    if(error/=0) then
       write(0,*) "interactions: Error reading file '",filename,"'" 
       stop 1
    end if
    read(10,*,iostat=error) line
    if(error/=0) then
       write(0,*) "interactions: Error reading file '",filename,"'"
       stop 1
    end if
    read(10,*,iostat=error) line
    if(error/=0) then
       write(0,*) "interactions: Error reading file '",filename,"'"
       stop 1
    end if

    ! Read the number of elements (number of elements must be 1), and ignoring the element name
    read(10,*,iostat=error) n
    if(error/=0) then
       write(0,*) "interactions: Error reading file '",filename,"'"
       stop 1
    end if
    if(n/=1) then
       write(0,*) "interactions: Error. File '",filename,"' seems to correspond to a ", &
            "multicomponent system; only unary systems supported."
       stop 1
    end if

    ! Read Nrho, drho, Nr, dr and cutoff
    read(10,*,iostat=error) Nrho,drho,Nr,dr,cutoff
    if(error/=0) then
       write(0,*) "interactions: Error reading file '",filename,"'"
       stop 1
    end if

    ! Read, but ignore, the line containing the atomic number, mass, lattice constant and lattice type
    read(10,*,iostat=error) line
    if(error/=0) then
       write(0,*) "interactions: Error reading file '",filename,"'"
       stop 1
    end if

    ! Initialise and read the embedding function
    allocate(F(Nrho))
    j=1
    do
       read(10,'(A)',iostat=error) line
       if(error/=0) then
          write(0,*) "interactions: Error reading file '",filename,"'"
          stop 1
       end if

       start_of_word=1
       prev_indx=0

       do i=1,len(line)
          ! Check if the current character in the line corresponds to a number (i.e., not a space).
          ! If so, then indx is NOT 0. If indx==0 then we have a blank
          indx = index('0123456789.-+Ee', line(i:i))

          ! We are at the start of a word if the current character is not a blank, but the previous
          ! was.
          if(indx/=0 .and. prev_indx==0) then
             start_of_word=i
          end if
          ! We are at the end of a word if the current character is a blank, but the previous
          ! was not.
          if(indx==0 .and. prev_indx/=0) then
             read(line(start_of_word:i-1), *) F(j)
             j=j+1
          end if

          if(j>Nrho) exit

          prev_indx=indx

       end do

       if(j>Nrho) exit

    end do

    ! Initialise and read the density function
    allocate(rho(Nr))
    j=1
    do
       read(10,'(A)',iostat=error) line
       if(error/=0) then
          write(0,*) "interactions: Error reading file '",filename,"'"
          stop 1
       end if

       start_of_word=1
       prev_indx=0

       do i=1,len(line)
          ! Check if the current character in the line corresponds to a number (i.e., not a space).
          ! If so, then indx is NOT 0. If indx==0 then we have a blank
          indx = index('0123456789.-+Ee', line(i:i))

          ! We are at the start of a word if the current character is not a blank, but the previous
          ! was.
          if(indx/=0 .and. prev_indx==0) then
             start_of_word=i
          end if
          ! We are at the end of a word if the current character is a blank, but the previous
          ! was not.
          if(indx==0 .and. prev_indx/=0) then
             read(line(start_of_word:i-1), *) rho(j)
             j=j+1
          end if

          if(j>Nr) exit

          prev_indx=indx

       end do

       if(j>Nr) exit

    end do

    ! Initialise and read the pair potential function. Note that the values in the file are actually r*phi(r),
    ! not phi(r). This is stored in the variable rphi.
    allocate(rphi(Nr))
    j=1
    do
       read(10,'(A)',iostat=error) line
       if(error/=0) then
          write(0,*) "interactions: Error reading file '",filename,"'"
          stop 1
       end if

       start_of_word=1
       prev_indx=0

       do i=1,len(line)
          ! Check if the current character in the line corresponds to a number (i.e., not a space).
          ! If so, then indx is NOT 0. If indx==0 then we have a blank
          indx = index('0123456789.-+Ee', line(i:i))

          ! We are at the start of a word if the current character is not a blank, but the previous
          ! was.
          if(indx/=0 .and. prev_indx==0) then
             start_of_word=i
          end if
          ! We are at the end of a word if the current character is a blank, but the previous
          ! was not.
          if(indx==0 .and. prev_indx/=0) then
             read(line(start_of_word:i-1), *) rphi(j)
             j=j+1
          end if

          if(j>Nr) exit

          prev_indx=indx

       end do

       if(j>Nr) exit

    end do

    ! Set cutoff to the maximum supported value if it is greater than (Nr-1)*dr
    if(cutoff>(Nr-1)*dr) cutoff=(Nr-1)*dr

    ! Initialise part_rho, and allocate - but not set any values to - part_rho_buffer
    allocate(part_rho(2,n_part))
    allocate(part_rho_buffer(2,n_part))
    call set_part_rho(1,Lx1,Ly1,Lz1,pos1)
    call set_part_rho(2,Lx2,Ly2,Lz2,pos2)

    ! Output the functions to files
    open(unit=10,file="F.dat")
    do i=1,Nrho
       write(10,*) (i-1)*drho,F(i)
    end do
    close(10)
    open(unit=10,file="rho.dat")
    do i=1,Nr
       write(10,*) (i-1)*dr,rho(i)
    end do
    close(10)
    open(unit=10,file="rphi.dat")
    do i=1,Nr
       write(10,*) (i-1)*dr,rphi(i)
    end do
    close(10)


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

    integer(ik), parameter :: unit=10

    write(unit,*) "Nrho= ",Nrho
    write(unit,*) "drho= ",drho
    write(unit,*) "Nr= ",Nr
    write(unit,*) "dr= ",dr
    write(unit,*) "cutoff= ",cutoff
    write(unit,*) "F= ",F
    write(unit,*) "rho= ",rho
    write(unit,*) "rphi= ",rphi
    write(unit,*) "part_rho= ",part_rho
    write(unit,*) "part_rho_buffer= ",part_rho_buffer
    
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

    integer(ik), parameter :: unit=10
    character(len=20) string
    integer(ik) :: error, n_part

    n_part=size(pos1,1)

    read(unit,*,iostat=error) string, Nrho
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'Nrho' from unit ",unit
       stop 1
    end if
    read(unit,*,iostat=error) string, drho
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'drho' from unit ",unit
       stop 1
    end if
    read(unit,*,iostat=error) string, Nr
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'Nr' from unit ",unit
       stop 1
    end if
    read(unit,*,iostat=error) string, dr
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'dr' from unit ",unit
       stop 1
    end if
    read(unit,*,iostat=error) string, cutoff
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'cutoff' from unit ",unit
       stop 1
    end if

    if(allocated(F)) then
       deallocate(F)
    end if
    if(allocated(rho)) then
       deallocate(rho)
    end if
    if(allocated(rphi)) then
       deallocate(rphi)
    end if
    allocate(F(Nrho))
    allocate(rho(Nr))
    allocate(rphi(Nr))
    allocate(part_rho(2,n_part))
    allocate(part_rho_buffer(2,n_part))

    read(unit,*,iostat=error) string, F
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'F' from unit ",unit
       stop 1
    end if
    read(unit,*,iostat=error) string, rho
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'rho' from unit ",unit
       stop 1
    end if
    read(unit,*,iostat=error) string, rphi
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'rphi' from unit ",unit
       stop 1
    end if
    read(unit,*,iostat=error) string, part_rho
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'part_rho' from unit ",unit
       stop 1
    end if
    read(unit,*,iostat=error) string, part_rho_buffer
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'part_rho_buffer' from unit ",unit
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


    part_rho=part_rho_buffer

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

    part_rho=part_rho_buffer

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

    integer(ik) :: i,j
    real(rk) :: sep
    real(rk) :: pair_energy
    real(rk) :: embedding_energy
    integer(ik) :: n_part

    n_part=size(pos,1)

    pair_energy=0.0_rk
    embedding_energy=0.0_rk
    do i=1,n_part
       ! Calculate the density associated with particle i, and store it in part_rho_buffer(lattice,i).
       ! At the same time we are resetting part_rho_buffer for 'lattice' to correspond to the 
       ! microstate in the argument
       part_rho_buffer(lattice,i)=0.0_rk
       do j=1,n_part
          if(j/=i) then
             sep=min_image_distance(pos(i,:),pos(j,:),Lx,Ly,Lz)
             if(sep<cutoff) then
                pair_energy=pair_energy+0.5_rk*phi_func(sep)
                part_rho_buffer(lattice,i)=part_rho_buffer(lattice,i)+rho_func(sep)
             end if
          end if
       end do
       embedding_energy=embedding_energy+F_func(part_rho_buffer(lattice,i))
    end do

    calc_energy_scratch = pair_energy + embedding_energy

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

    integer(ik) :: j
    real(rk) :: sep
    real(rk) :: sep_new
    real(rk) :: delta_pair_energy
    real(rk) :: delta_embedding_energy
    ! The change in density for microstate with r_new relative to the current microstate (r)
    real(rk), dimension(size(pos,1)) :: delta_rho    

    integer(ik) :: n_part

    n_part=size(pos,1)

    ! EAM code from an older version of monteswitch, which is optimised below...
    !
    !    ! Calculate the change in the pair energy
    !    delta_pair_energy=0.0_rk
    !    do j=1,n_part
    !       if(j/=i) then
    !          sep=min_image_distance(pos(j,:),pos(i,:),Lx,Ly,Lz)
    !          sep_new=min_image_distance(pos(j,:),pos_new,Lx,Ly,Lz)
    !          if(sep_new<cutoff) delta_pair_energy=delta_pair_energy+phi_func(sep_new)
    !          if(sep<cutoff) delta_pair_energy=delta_pair_energy-phi_func(sep)
    !       end if
    !    end do
    !
    !    ! Calculate the change in the embedding energy
    !
    !    ! Calculate the changes in all particles' densities for  r_new relative to r_ref
    !    delta_rho=0.0_rk
    !    do j=1,n_part
    !       if(j/=i) then
    !          sep=min_image_distance(pos(j,:),pos(i,:),Lx,Ly,Lz)
    !          sep_new=min_image_distance(pos(j,:),pos_new,Lx,Ly,Lz)
    !          if(sep_new<cutoff) then
    !             delta_rho(j)=delta_rho(j)+rho_func(sep_new)
    !             delta_rho(i)=delta_rho(i)+rho_func(sep_new)
    !          end if
    !          if(sep<cutoff) then
    !             delta_rho(j)=delta_rho(j)-rho_func(sep)
    !             delta_rho(i)=delta_rho(i)-rho_func(sep)
    !          end if
    !       end if
    !    end do
    
    ! Combine the above commented out code into one more efficient loop
    
    delta_pair_energy=0.0_rk
    delta_rho=0.0_rk
    do j=1,n_part
       if(j/=i) then
          sep=min_image_distance(pos(j,:),pos(i,:),Lx,Ly,Lz)
          sep_new=min_image_distance(pos(j,:),pos_new,Lx,Ly,Lz)
          if(sep_new<cutoff) then
              delta_pair_energy=delta_pair_energy+phi_func(sep_new)

              delta_rho(j)=delta_rho(j)+rho_func(sep_new)
              delta_rho(i)=delta_rho(i)+rho_func(sep_new)
          end if
          if(sep<cutoff) then
              delta_pair_energy=delta_pair_energy-phi_func(sep)
              
              delta_rho(j)=delta_rho(j)-rho_func(sep)
              delta_rho(i)=delta_rho(i)-rho_func(sep)
          end if
       end if
    end do

    ! Use the change to calculate the current densities, and store it in part_rho_buffer
    part_rho_buffer(lattice,:) = part_rho(lattice,:) + delta_rho

    ! Calculate the chage in embedding energy from the changes in densities
    delta_embedding_energy=0.0_rk
    do j=1,n_part
       delta_embedding_energy = delta_embedding_energy + &
            F_func(part_rho_buffer(lattice,j)) - F_func(part_rho(lattice,j))
    end do

    calc_energy_part_move = delta_pair_energy + delta_embedding_energy

end function calc_energy_part_move




! Set the part of part_rho for 'lattice' to correspond to the microstate
! in the argument
subroutine set_part_rho(lattice,Lx,Ly,Lz,r)
    integer(ik), intent(in) :: lattice
    real(rk), intent(in) :: Lx
    real(rk), intent(in) :: Ly
    real(rk), intent(in) :: Lz
    real(rk), intent(in), dimension(:,:) :: r

    integer(ik) :: i,j
    real(rk) :: sep
    integer(ik) :: n_part
    
    n_part=size(r,1)

    do i=1,n_part
       part_rho(lattice,i)=0.0_rk
       do j=1,n_part
          if(j/=i) then
             sep=min_image_distance(r(i,:),r(j,:),Lx,Ly,Lz)
             if(sep<cutoff) then
                part_rho(lattice,i)=part_rho(lattice,i)+rho_func(sep)
             end if
          end if
       end do
    end do
end subroutine set_part_rho




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




! Returns the embedding energy for a given density.
function F_func(rho)
    real(rk), intent(in) :: rho
    real(rk) :: F_func
    real(rk) :: bin_cont
    integer(ik) :: bin_below, bin_above
    ! F(bin) corresponds to the embedding function at rho=(bin-1)*drho.
    ! Hence bin=rho/drho+1 is the 'continuous' bin number
    bin_cont=rho/drho+1
    if(bin_cont<1 .or. bin_cont>Nrho) then
       write(0,*) "interactions_mod: Error in F_func: rho not supported. rho = ",rho
       stop 1
    end if
    ! Perform a linear interpolation to get the value of the embedding function at rho
    bin_below=floor(bin_cont)
    bin_above=ceiling(bin_cont)
    F_func=F(bin_below) + (F(bin_above)-F(bin_below))*(bin_cont-bin_below)
end function F_func




! Returns the density corresponding to a given separation.
function rho_func(r)
    real(rk), intent(in) :: r
    real(rk) :: rho_func
    real(rk) :: bin_cont
    integer(ik) :: bin_below, bin_above
    ! rho(bin) corresponds to the density function at r=(bin-1)*dr.
    ! Hence bin=r/dr+1 is the 'continuous' bin number
    bin_cont=r/dr+1
    if(bin_cont<1 .or. bin_cont>Nr) then
       write(0,*) "interactions_mod: Error in rho_func: r not supported. r = ",r
       stop 1
    end if
    ! Perform a linear interpolation to get the value of the density function at r
    bin_below=floor(bin_cont)
    bin_above=ceiling(bin_cont)
    rho_func=rho(bin_below) + (rho(bin_above)-rho(bin_below))*(bin_cont-bin_below)
end function rho_func




! Returns the value of the pair potential corresponding to a given separation.
function phi_func(r)
    real(rk), intent(in) :: r
    real(rk) :: phi_func
    real(rk) :: bin_cont
    integer(ik) :: bin_below, bin_above
    ! rphi(bin) corresponds to the value of r*phi at r=(bin-1)*dr
    ! Hence bin=r/dr+1 is the 'continuous' bin number
    bin_cont=r/dr+1
    if(bin_cont<1 .or. bin_cont>Nr) then
       write(0,*) "interactions_mod: Error in rho_func: r not supported. r = ",r
       stop 1
    end if
    ! Perform a linear interpolation to get the value of the phi function at r
    bin_below=floor(bin_cont)
    bin_above=ceiling(bin_cont)
    phi_func=rphi(bin_below) + (rphi(bin_above)-rphi(bin_below))*(bin_cont-bin_below)
    phi_func=phi_func/r
end function phi_func




end module interactions_mod
