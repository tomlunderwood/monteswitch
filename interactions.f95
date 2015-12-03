


module interactions_mod


    use kinds_mod
    

contains



! Initialises interactions, using the specified initial configurations for
! lattices 1 and 2. I recommend reading variables which define the interactions
! from a file 'interactions_in'
subroutine initialise_interactions(Lx1, Ly1, Lz1, species1, pos1, Lx2, Ly2, Lz2, species2, pos2)
    real(rk), intent(in) :: Lx1, Ly1, Lz1, Lx2, Ly2, Lz2
    integer(ik), intent(in), dimension(:) :: species1, species2
    real(rk), intent(in), dimension(:,:) :: pos1, pos2


end subroutine initialise_interactions



! Export all interactions variables - for checkpointing. There are two options
! for this. If one outputs the varaibles to unit 10, without opening or closing that 
! unit in this procedure, then the interactions variables will be stored within
! the 'state' file which contains all other monteswitch variables. If one wishes to
! use a separate file or files to checkpoint the interactions variables, then that
! is fine - just don't use unit 10!
subroutine export_interactions()


end subroutine export_interactions



! Import all interactions variables from a checkpoint file. The format read from the
! file should correspond to that output by 'export_interactions' above. If one is
! importing from within a 'state' file as described above, then use unit 10, but
! do not open or close that unit!
subroutine import_interactions()


end subroutine import_interactions



! WHAT TO DO WITH THIS?
subroutine after_accepted_interactions()


end subroutine after_accepted_interactions




! Calculate the energy for the specificed configuration of the specified lattice
function calc_energy_scratch(lattice,Lx,Ly,Lz,species,pos)
    integer(ik), intent(in) :: lattice
    real(rk), intent(in) :: Lx, Ly, Lz
    integer(ik), dimension(:), intent(in) :: species
    real(rk), dimension(:,:), intent(in) :: pos
    real(rk) :: calc_energy_scratch

end function calc_energy_scratch



! Calculate the energy for the specificed configuration of the specified lattice given
! that particle i has moved
function calc_energy_part_move(lattice,Lx,Ly,Lz,species,pos,pos_new,i)
    integer(ik), intent(in) :: lattice, i
    real(rk), intent(in) :: Lx, Ly, Lz
    integer(ik), dimension(:), intent(in) :: species
    real(rk), dimension(:,:), intent(in) :: pos, pos_new
    real(rk) :: calc_energy_part_move

end function calc_energy_part_move




end module interactions_mod




!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$  ! interactions.f95 for penetrable sphere potential: epsilon for r<sigma; 0 otherwise
!!$  ! Note that hard spheres are realised for epsilon*beta >> 1.
!!$  ! Author: Tom L Underwood
!!$  ! Adapted from interactions_PAIR_TEMPLATE.f95, but note that 'cutoff' is superfluous
!!$  ! here, as is the function 'pair_potential_trunc', and hence I have removed them.
!!$
!!$  ! <<<<<<<<<<<< USER-DEFINED CODE >>>>>1>>>>>>> 
!!$  ! Insert code defining the free parameters for the potential here.
!!$  ! Example (for Lennard-Jones potential):
!!$  !  real(rk) :: lj_epsilon
!!$  !  real(rk) :: lj_sigma
!!$  !
!!$  integer(ik) :: dspecies
!!$  real(rk) :: epsilon
!!$  real(rk), allocatable, dimension(:) :: sigma
!!$  ! <<<<<<<<<<<< END OF USER-DEFINED CODE >>>>>>>>>>>> 
!!$
!!$  ! Variables for the potential cutoff and lists of interacting particles
!!$  real(rk) :: list_cutoff
!!$  integer(ik) :: list_size
!!$  integer(ik), dimension(:,:), allocatable :: list_1
!!$  integer(ik), dimension(:,:), allocatable :: list_2
!!$
!!$
!!$contains
!!$
!!$
!!$  subroutine initialise_interactions(filename)
!!$    character(*), intent(in) :: filename
!!$    character(len=20) string
!!$    integer(ik) :: error
!!$    integer(ik) :: i,j,n
!!$
!!$    ! Read the interactions variables from the file 'filename'
!!$    open(unit=10,file=filename,iostat=error,status="old")
!!$    if(error/=0) then
!!$       write(0,*) "interactions: Error. Problem opening file '",filename,"'"
!!$       stop 1
!!$    end if
!!$
!!$    ! <<<<<<<<<<<< USER-DEFINED CODE >>>>>>>>>>>> 
!!$    ! Insert code which reads the free parameters for the potential from unit 10, using list-directed
!!$    ! formatting. Note that the way in which the free parameters are read must match the way in which
!!$    ! they are written - in the USER-DEFINED CODE block below.
!!$    ! Example (for Lennard-Jones potential, including an error message and exit code of 1 if there is
!!$    ! an error):
!!$    !  read(10,*,iostat=error) string, lj_epsilon
!!$    !  if(error/=0) then
!!$    !     write(0,*) "interactions: Error. Problem reading 'lj_epsilon' from file '",filename,"'"
!!$    !     stop 1
!!$    !  end if
!!$    !  read(10,*,iostat=error) string, lj_sigma
!!$    !  if(error/=0) then
!!$    !     write(0,*) "interactions: Error. Problem reading 'lj_sigma' from file '",filename,"'"
!!$    !     stop 1
!!$    !  end if
!!$    !
!!$    read(10,*,iostat=error) string, dspecies
!!$    if(error/=0) then
!!$        write(0,*) "interactions: Error. Problem reading 'dspecies' from file '",filename,"'"
!!$        stop 1
!!$    end if
!!$    read(10,*,iostat=error) string, epsilon
!!$    if(error/=0) then
!!$        write(0,*) "interactions: Error. Problem reading 'epsilon' from file '",filename,"'"
!!$        stop 1
!!$    end if
!!$    allocate(sigma(dspecies))
!!$    read(10,*,iostat=error) string, sigma
!!$    if(error/=0) then
!!$        write(0,*) "interactions: Error. Problem reading 'sigma' from file '",filename,"'"
!!$        stop 1
!!$    end if
!!$    ! <<<<<<<<<<<< END OF USER-DEFINED CODE >>>>>>>>>>>> 
!!$
!!$    read(10,*,iostat=error) string, list_cutoff
!!$    if(error/=0) then
!!$       write(0,*) "interactions: Error. Problem reading 'list_cutoff' from file '",filename,"'"
!!$       stop 1
!!$    end if
!!$    read(10,*,iostat=error) string, list_size
!!$    if(error/=0) then
!!$       write(0,*) "interactions: Error. Problem reading 'list_size' from file '",filename,"'"
!!$       stop 1
!!$    end if
!!$    close(10)
!!$
!!$    ! Initialise lists
!!$    if(allocated(list_1)) then
!!$       deallocate(list_1)
!!$    end if
!!$    if(allocated(list_2)) then
!!$       deallocate(list_2)
!!$    end if
!!$    allocate(list_1(n_part,list_size))
!!$    allocate(list_2(n_part,list_size))
!!$
!!$    list_1=0
!!$    do i=1,n_part
!!$       n=1
!!$       do j=1,n_part
!!$          if(min_image_distance(R_1(i,:),R_1(j,:),Lx(1),Ly(1),Lz(1))<list_cutoff) then
!!$             list_1(i,n)=j
!!$             n=n+1
!!$          end if
!!$       end do
!!$    end do
!!$    list_2=0
!!$    do i=1,n_part
!!$       n=1
!!$       do j=1,n_part
!!$          if(min_image_distance(R_2(i,:),R_2(j,:),Lx(2),Ly(2),Lz(2))<list_cutoff) then
!!$             list_2(i,n)=j
!!$             n=n+1
!!$          end if
!!$       end do
!!$    end do
!!$
!!$  end subroutine initialise_interactions
!!$
!!$
!!$
!!$
!!$  subroutine export_interactions_state(unit)
!!$    integer(ik), intent(in) :: unit
!!$
!!$    ! <<<<<<<<<<<< USER-DEFINED CODE >>>>>>>>>>>> 
!!$    ! Insert code which writes the free parameters for the potential to unit 'unit', using list-directed
!!$    ! formatting. Note that the way in which the free parameters are read must match the way in which
!!$    ! they are read - in the USER-DEFINED CODE block above.
!!$    ! Example (for Lennard-Jones potential):
!!$    !  write(unit,*) "lj_epsilon= ",lj_epsilon
!!$    !  write(unit,*) "lj_sigma= ",lj_sigma
!!$    !
!!$    write(unit,*) "dspecies= ",dspecies
!!$    write(unit,*) "epsilon= ",epsilon
!!$    write(unit,*) "sigma= ",sigma
!!$    ! <<<<<<<<<<<< END OF USER-DEFINED CODE >>>>>>>>>>>> 
!!$
!!$    write(unit,*) "list_cutoff= ",list_cutoff
!!$    write(unit,*) "list_size= ",list_size
!!$    write(unit,*) "list_1= ",list_1
!!$    write(unit,*) "list_2= ",list_2
!!$  end subroutine export_interactions_state
!!$
!!$
!!$
!!$
!!$  subroutine import_interactions_state(unit)
!!$    integer(ik), intent(in) :: unit
!!$    character(len=20) string
!!$    integer(ik) :: error
!!$    
!!$    ! <<<<<<<<<<<< USER-DEFINED CODE >>>>>>>>>>>>
!!$    ! Insert code which reads the free parameters for the potential from unit 'unit', using list-directed
!!$    ! formatting. Note that the way in which the free parameters are read must match the way in which
!!$    ! they are written - in the USER-DEFINED CODE block below.
!!$    ! Example (for Lennard-Jones potential, including an error message and exit code of 1 if there is
!!$    ! an error):
!!$    !  read(unit,*,iostat=error) string, lj_epsilon
!!$    !  if(error/=0) then
!!$    !     write(0,*) "interactions: Error. Problem reading 'lj_epsilon' from unit ",unit
!!$    !     stop 1
!!$    !  end if
!!$    !  read(unit,*,iostat=error) string, lj_sigma
!!$    !  if(error/=0) then
!!$    !     write(0,*) "interactions: Error. Problem reading 'lj_sigma' from unit ",unit
!!$    !     stop 1
!!$    !  end if
!!$    !
!!$    read(unit,*,iostat=error) string, dspecies
!!$    if(error/=0) then
!!$       write(0,*) "interactions: Error. Problem reading 'dspecies' from unit ",unit
!!$       stop 1
!!$    end if
!!$    read(unit,*,iostat=error) string, epsilon
!!$    if(error/=0) then
!!$       write(0,*) "interactions: Error. Problem reading 'epsilon' from unit ",unit
!!$       stop 1
!!$    end if
!!$    allocate(sigma(dspecies))
!!$    read(unit,*,iostat=error) string, sigma
!!$    if(error/=0) then
!!$       write(0,*) "interactions: Error. Problem reading 'sigma' from unit ",unit
!!$       stop 1
!!$    end if
!!$    ! <<<<<<<<<<<< END OF USER-DEFINED CODE >>>>>>>>>>>> 
!!$
!!$    read(unit,*,iostat=error) string, list_cutoff
!!$    if(error/=0) then
!!$       write(0,*) "interactions: Error. Problem reading 'list_cutoff' from unit ",unit
!!$       stop 1
!!$    end if
!!$    read(unit,*,iostat=error) string, list_size
!!$    if(error/=0) then
!!$       write(0,*) "interactions: Error. Problem reading 'list_size' from unit ",unit
!!$       stop 1
!!$    end if
!!$
!!$    if(allocated(list_1)) then
!!$       deallocate(list_1)
!!$    end if
!!$    if(allocated(list_2)) then
!!$       deallocate(list_2)
!!$    end if
!!$    allocate(list_1(n_part,list_size))
!!$    allocate(list_2(n_part,list_size))
!!$
!!$    read(unit,*,iostat=error) string, list_1
!!$    if(error/=0) then
!!$       write(0,*) "interactions: Error. Problem reading 'list_1' from unit ",unit
!!$       stop 1
!!$    end if
!!$    read(unit,*,iostat=error) string, list_2
!!$    if(error/=0) then
!!$       write(0,*) "interactions: Error. Problem reading 'list_2' from unit ",unit
!!$       stop 1
!!$    end if
!!$
!!$  end subroutine import_interactions_state
!!$
!!$
!!$
!!$
!!$  ! This does nothing for pair potentials
!!$  subroutine after_accepted_interactions()
!!$
!!$  end subroutine after_accepted_interactions
!!$
!!$
!!$
!!$  ! Returns the energy of the system calculated from scratch for the given lattice
!!$  function calc_energy_scratch(lattice,Lx,Ly,Lz,r)
!!$    integer(ik), intent(in) :: lattice
!!$    real(rk), intent(in) :: Lx
!!$    real(rk), intent(in) :: Ly
!!$    real(rk), intent(in) :: Lz
!!$    real(rk), intent(in), dimension(:,:) :: r
!!$    real(rk) :: calc_energy_scratch
!!$    integer(ik) :: i,j,n, speci,specj
!!$    real(rk) :: sep
!!$    
!!$    calc_energy_scratch=0.0_rk
!!$    do i=1,ubound(r,1)
!!$
!!$        select case(lattice)
!!$        case(1)
!!$            speci=spec_1(i)
!!$        case(2)
!!$            speci=spec_2(i)
!!$        case default
!!$            write(0,*) "interactions: Error. 'lattice' is not 1 or 2."
!!$            stop 1
!!$        end select
!!$        
!!$        n=1
!!$        do
!!$            select case(lattice)
!!$            case(1)
!!$                j=list_1(i,n)
!!$            case(2)
!!$                j=list_2(i,n)
!!$            case default
!!$                write(0,*) "interactions: Error. 'lattice' is not 1 or 2."
!!$                stop 1
!!$            end select
!!$           
!!$            if(j==0) then
!!$                exit
!!$            else
!!$                if(j/=i) then
!!$                    select case(lattice)
!!$                    case(1)
!!$                        specj=spec_1(j)
!!$                    case(2)
!!$                        specj=spec_2(j)
!!$                    end select
!!$                    sep=min_image_distance(r(i,:),r(j,:),Lx,Ly,Lz)
!!$                    calc_energy_scratch=calc_energy_scratch+pair_potential(sep,speci,specj)
!!$                end if
!!$            end if
!!$            n=n+1
!!$        end do
!!$    end do
!!$    calc_energy_scratch=0.5*calc_energy_scratch
!!$
!!$  end function calc_energy_scratch
!!$
!!$
!!$
!!$  ! Returns energy after a particle is moved for the given lattice
!!$  function calc_energy_part_move(lattice,Lx,Ly,Lz,r,r_new,i)
!!$    integer(ik), intent(in) :: lattice
!!$    real(rk), intent(in) :: Lx
!!$    real(rk), intent(in) :: Ly
!!$    real(rk), intent(in) :: Lz
!!$    real(rk), intent(in), dimension(:,:) :: r
!!$    real(rk), intent(in), dimension(:,:) :: r_new
!!$    integer(ik), intent(in) :: i
!!$    real(rk) :: calc_energy_part_move
!!$    real(rk) :: sep, sep_new
!!$    integer(ik) :: j,n, speci,specj
!!$
!!$    select case(lattice)
!!$    case(1)
!!$        speci=spec_1(i)
!!$    case(2)
!!$        speci=spec_2(i)
!!$    case default
!!$        write(0,*) "interactions: Error. 'lattice' is not 1 or 2."
!!$        stop 1
!!$    end select
!!$
!!$    calc_energy_part_move=0.0_rk
!!$    n=1
!!$    do
!!$       select case(lattice)
!!$       case(1)
!!$          j=list_1(i,n)
!!$       case(2)
!!$          j=list_2(i,n)
!!$       case default
!!$          write(0,*) "interactions: Error. 'lattice' is not 1 or 2."
!!$          stop 1
!!$       end select
!!$
!!$       if(j==0) then
!!$          exit
!!$       else
!!$          if(j/=i) then
!!$              select case(lattice)
!!$              case(1)
!!$                  specj=spec_1(j)
!!$              case(2)
!!$                  specj=spec_2(j)
!!$              end select
!!$              sep = min_image_distance(r(i,:),r(j,:),Lx,Ly,Lz)
!!$              sep_new = min_image_distance(r_new(i,:),r_new(j,:),Lx,Ly,Lz)
!!$              calc_energy_part_move = calc_energy_part_move &
!!$                  + pair_potential(sep_new, speci, specj) &
!!$                  - pair_potential(sep, speci, specj)
!!$          end if
!!$       end if
!!$       n=n+1
!!$    end do
!!$
!!$
!!$  end function calc_energy_part_move
!!$
!!$
!!$
!!$
!!$  function min_image_distance(r_1,r_2,Lx,Ly,Lz)
!!$    real(rk), dimension(3), intent(in) :: r_1, r_2
!!$    real(rk), intent(in) :: Lx, Ly, Lz
!!$    real(rk) :: min_image_distance
!!$    real(rk) :: xsep, ysep, zsep
!!$    ! Calculate the x-sep
!!$    xsep=abs(r_2(1)-r_1(1))
!!$    xsep=xsep-Lx*floor(2.0_rk*xsep/Lx)
!!$    ! Calculate the y-sep
!!$    ysep=abs(r_2(2)-r_1(2))
!!$    ysep=ysep-Ly*floor(2.0_rk*ysep/Ly)
!!$    ! Calculate the z-sep
!!$    zsep=abs(r_2(3)-r_1(3))
!!$    zsep=zsep-Lz*floor(2.0_rk*zsep/Lz)
!!$    ! Calculate the distance
!!$    min_image_distance=sqrt(xsep*xsep+ysep*ysep+zsep*zsep)
!!$  end function min_image_distance
!!$
!!$
!!$
!!$
!!$  ! The 'pure' pair potential, without truncation, for a pair with the specified species
!!$  ! separated by r
!!$  function pair_potential(r, speci, specj)
!!$    real(rk), intent(in) :: r
!!$    integer(ik), intent(in) :: speci, specj
!!$    real(rk) :: pair_potential
!!$    ! <<<<<<<<<<<< USER-DEFINED CODE >>>>>>>>>>>> 
!!$    ! Insert code corresponding to the 'pure' pair potential, without truncation,
!!$    ! and using the free parameters defined above.
!!$    ! Example (for Lennard-Jones potential):
!!$    !  pair_potential=4.0_rk*lj_epsilon*( (lj_sigma/r)**12-(lj_sigma/r)**6 )
!!$    !
!!$    if(r < (sigma(speci)+sigma(specj))/2 ) then
!!$        pair_potential=epsilon
!!$    else
!!$        pair_potential=0.0_rk
!!$    end if
!!$    ! <<<<<<<<<<<< END OF USER-DEFINED CODE >>>>>>>>>>>> 
!!$  end function pair_potential
!!$
