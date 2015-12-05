! 'interactions.f95' file for monteswitch for Einstein Crystal
! Author: Tom L Underwood

! The Hamiltonian corresponding to this module is as follows. In lattice 1 the energy
! of each particle corresponds to a harmonic oscillator about the origin in which the
! spring constant is 'alpha_1'; the spring constant for lattice 2 is 'alpha_2'.
!
! 'alpha_1' and 'alpha_2' are imported at initialisation from a file 'interactions_in'.
! The format of the file is as follows. On the first line there are two tokens. 
! The first token is a character(len=20) variable (I recommend: 'alpha_1=');
! the second is the value of 'alpha_1'. The second line is similar, but for 'alpha_2'.
!
! Comments beginning with '!*' are inherited from the template which was used to create this
! file
!
module interactions_mod


    !* The module 'kinds_mod' (in 'kinds_mod.f95') contains the real and integer kinds
    !* which monteswitch uses: real(kind=rk) and integer(kind=ik). Please use these
    !* kinds throughout this module to allow portability.
    use kinds_mod
    
    implicit none


    real(rk) :: alpha_1
    real(rk) :: alpha_2


    !* The below public procedures are called by 'monteswitch_mod' (in 'monteswitch_mod.f95'),
    !* and must be 'filled in' by the user.
    private
    public :: initialise_interactions, export_interactions, import_interactions, &
        after_accepted_part_interactions, after_accepted_vol_interactions, &
        after_accepted_lattice_interactions, after_all_interactions, &
        calc_energy_scratch, calc_energy_part_move


contains


!* Initialises the  variables in this module. The initial configurations for both
!* lattices are provided as arguments in case this information is required. Unless the interactions
!* are 'hard-coded', the initialisation will involve reading from a file. This should be done here.
subroutine initialise_interactions(Lx1, Ly1, Lz1, species1, pos1, Lx2, Ly2, Lz2, species2, pos2)
    !* Dimensions of the initial (orthorhombic) supercell for lattices 1 and 2
    !* in each Cartesian dimension
    real(rk), intent(in) :: Lx1, Ly1, Lz1, Lx2, Ly2, Lz2
    !* Array containing the species of each particle for each lattice: e.g., species1(i) is the
    !* species of particle i in lattice 1
    integer(ik), intent(in), dimension(:) :: species1, species2
    !* Initial positions (Cartesian) of the particles for lattices 1 and 2: e.g., pos1(i,1) is
    !* the x-coordinate of particle i in lattice 1, pos1(i,2) is the y-coordinate, and pos1(i,3)
    !* is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2

    character(len=20) string
    integer(ik) :: error

    ! Open the file 'interactions_in' and import the module variables
    open(unit=10,file="interactions_in",iostat=error,status="old")
    if(error/=0) then
       write(0,*) "interactions: Error. Problem opening file 'interactions_in'"
       stop 1
    end if
    read(10,*,iostat=error) string, alpha_1
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'alpha_1' from file 'interactions_in'"
       stop 1
    end if
    read(10,*,iostat=error) string, alpha_2
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'alpha_2' from file 'interactions_in'"
       stop 1
    end if
    close(10)

end subroutine initialise_interactions



!* Export all of the variables in this module - for the purposes of checkpointing. There are 
!* two options for this. If one outputs the variables to unit 10, without opening or closing that 
!* unit in this procedure, then the variables will be output to, and stored within,
!* the 'state' file which contains all other monteswitch variables, and is used for resuming
!* old simulations as well as post-processing. If one wishes to use a separate file or files to 
!* checkpoint the variables in this module, then one is free to do so in this procedure - 
!* just don't use unit 10! Note that the format output to the checkpoint file(s) by this procedure
!* must correspond to the format read in from these files by the procedure 'import_interactions()'.
subroutine export_interactions()

    write(10,*) "alpha_1= ",alpha_1
    write(10,*) "alpha_2= ",alpha_2
    
end subroutine export_interactions




!* Import all variables in this module from file(s) - for the purposes of checkpointing. The format 
!* read from the file(s) should correspond to that output by 'export_interactions' above. If one is
!* importing from within a 'state' file as described in that procedure, then use unit 10, but
!* do not open or close that unit!
subroutine import_interactions()

    character(len=20) string
    integer(ik) :: error
    
    read(10,*,iostat=error) string, alpha_1
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'alpha_1' from unit 10 "
       stop 1
    end if
    read(10,*,iostat=error) string, alpha_2
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'alpha_2' from unit 10"
       stop 1
    end if

end subroutine import_interactions




!* Performs any tasks required for the variables in this module, after a particle move for particle i is accepted
subroutine after_accepted_part_interactions(i, Lx1, Ly1, Lz1, species1, pos1, Lx2, Ly2, Lz2, species2, pos2)
    !* The particle which has just been moved (where the move has been accepted)
    integer(ik), intent(in) :: i
    !* Current dimensions of the initial (orthorhombic) supercell for lattices 1 and 2
    !* in each Cartesian dimension (after the move has been accepted)
    real(rk), intent(in) :: Lx1, Ly1, Lz1, Lx2, Ly2, Lz2
    !* Array containing the species of each particle for each lattice: e.g., species1(i) is the
    !* species of particle i in lattice 1
    integer(ik), intent(in), dimension(:) :: species1, species2
    !* Current positions (Cartesian) of the particles for lattices 1 and 2 (after the move has been
    !* accepted): e.g., pos1(i,1) is the x-coordinate of particle i in lattice 1, pos1(i,2) is the 
    !* y-coordinate, and pos1(i,3) is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2

    return

end subroutine after_accepted_part_interactions




!* Performs any tasks required for the variables in this module, after a volume move is accepted
subroutine after_accepted_vol_interactions(Lx1, Ly1, Lz1, species1, pos1, Lx2, Ly2, Lz2, species2, pos2)
    !* Current dimensions of the initial (orthorhombic) supercell for lattices 1 and 2
    !* in each Cartesian dimension (after the move has been accepted)
    real(rk), intent(in) :: Lx1, Ly1, Lz1, Lx2, Ly2, Lz2
    !* Array containing the species of each particle for each lattice: e.g., species1(i) is the
    !* species of particle i in lattice 1
    integer(ik), intent(in), dimension(:) :: species1, species2
    !* Current positions (Cartesian) of the particles for lattices 1 and 2 (after the move has been
    !* accepted): e.g., pos1(i,1) is the x-coordinate of particle i in lattice 1, pos1(i,2) is the 
    !* y-coordinate, and pos1(i,3) is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2

    return

end subroutine after_accepted_vol_interactions



!* Performs any tasks required for the variables in this module, after a lattice move is accepted
subroutine after_accepted_lattice_interactions(Lx1, Ly1, Lz1, species1, pos1, Lx2, Ly2, Lz2, species2, pos2)
    !* Current dimensions of the initial (orthorhombic) supercell for lattices 1 and 2
    !* in each Cartesian dimension (after the move has been accepted)
    real(rk), intent(in) :: Lx1, Ly1, Lz1, Lx2, Ly2, Lz2
    !* Array containing the species of each particle for each lattice: e.g., species1(i) is the
    !* species of particle i in lattice 1
    integer(ik), intent(in), dimension(:) :: species1, species2
    !* Current positions (Cartesian) of the particles for lattices 1 and 2 (after the move has been
    !* accepted): e.g., pos1(i,1) is the x-coordinate of particle i in lattice 1, pos1(i,2) is the 
    !* y-coordinate, and pos1(i,3) is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2

    return

end subroutine after_accepted_lattice_interactions




!* Performs any tasks required for the variables in this module after ALL moves (accepted or rejected)
subroutine after_all_interactions(Lx1, Ly1, Lz1, species1, pos1, Lx2, Ly2, Lz2, species2, pos2)
    !* Current dimensions of the initial (orthorhombic) supercell for lattices 1 and 2
    !* in each Cartesian dimension
    real(rk), intent(in) :: Lx1, Ly1, Lz1, Lx2, Ly2, Lz2
    !* Array containing the species of each particle for each lattice: e.g., species1(i) is the
    !* species of particle i in lattice 1
    integer(ik), intent(in), dimension(:) :: species1, species2
    !* Current positions (Cartesian) of the particles for lattices 1 and 2: e.g., pos1(i,1) is 
    !* the x-coordinate of particle i in lattice 1, pos1(i,2) is the y-coordinate, and pos1(i,3) 
    !* is the z-coordinate
    real(rk), intent(in), dimension(:,:) :: pos1, pos2

    return

end subroutine after_all_interactions




!* Returns the energy for the specificed configuration of the specified lattice
function calc_energy_scratch(lattice, Lx, Ly, Lz, species, pos)
    !* The lattice (1 or 2)
    integer(ik), intent(in) :: lattice
    !* Dimensions of the (orthorhombic) supercell
    real(rk), intent(in) :: Lx, Ly, Lz
    !* Array containing the species of each particle for each lattice: e.g., species(i) is the
    !* species of particle i
    integer(ik), dimension(:), intent(in) :: species
    !* Positions (Cartesian) of the particles: e.g., pos(i,1) is the x-coordinate of particle 
    !* i, pos1(i,2) is the y-coordinate, and pos1(i,3) is the z-coordinate
    real(rk), dimension(:,:), intent(in) :: pos
    real(rk) :: calc_energy_scratch

    real(rk), dimension(3) :: zero = (/0.0_rk,0.0_rk,0.0_rk/)
    integer(ik) :: i

    calc_energy_scratch=0.0_rk

    select case(lattice)
    case(1)

        do i=1,size(pos,1)
     
            calc_energy_scratch = calc_energy_scratch + alpha_1 * min_image_distance(pos(i,:),zero,Lx,Ly,Lz)**2

        end do
        
    case(2)

        do i=1,size(pos,1)
     
            calc_energy_scratch = calc_energy_scratch + alpha_2 * min_image_distance(pos(i,:),zero,Lx,Ly,Lz)**2

        end do

    case default
    
        write(0,*) "Error: interactions: 'lattice' was not 1 or 2."
        stop 1

    end select

end function calc_energy_scratch




!* Returns the energy for the specificed configuration of the specified lattice given
!* that particle i has moved
function calc_energy_part_move(lattice, Lx, Ly, Lz, species, pos, pos_new, i)
    !* The lattice (1 or 2)
    integer(ik), intent(in) :: lattice
    !* The particle which has just been moved
    integer(ik), intent(in) :: i
    !* Dimensions of the (orthorhombic) supercell
    real(rk), intent(in) :: Lx, Ly, Lz
    !* Array containing the species of each particle for each lattice: e.g., species(i) is the
    !* species of particle i
    integer(ik), dimension(:), intent(in) :: species
    !* Positions (Cartesian) of the particles BEFORE particle i has been moved: e.g., pos(j,1) is 
    !* the x-coordinate of particle j, pos(j,2) is the y-coordinate, and pos(j,3) is the
    !* z-coordinate
    real(rk), dimension(:,:), intent(in) :: pos
    !* Position of particle i AFTER the particle has been moved
    real(rk), dimension(3), intent(in) :: pos_new
    real(rk) :: calc_energy_part_move

    real(rk), dimension(3) :: zero = (/0.0_rk,0.0_rk,0.0_rk/)

    select case(lattice)
    case(1)

       calc_energy_part_move = alpha_1 * min_image_distance(pos_new,zero,Lx,Ly,Lz)**2 &
           - alpha_1 * min_image_distance(pos(i,:),zero,Lx,Ly,Lz)**2

    case(2)

        calc_energy_part_move = alpha_2 * min_image_distance(pos_new,zero,Lx,Ly,Lz)**2 &
            - alpha_2 * min_image_distance(pos(i,:),zero,Lx,Ly,Lz)**2

    case default

       write(0,*) "Error: interactions: 'lattice' was not 1 or 2."
       stop 1

    end select

end function calc_energy_part_move




! Returns the distance between positions r_1 and r_2 within the orthorhombc
! cell with the specified dimensions according to the minimum image convention.
! Note that r_1 and r_2 must be such that 0<=x<Lx, 0<=y<Ly, and 0<=z<Lz.
function min_image_distance(r_1,r_2,Lx,Ly,Lz)
    ! Positions to consider
    real(rk), dimension(3), intent(in) :: r_1, r_2
    ! Dimensions of cell
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




end module interactions_mod
