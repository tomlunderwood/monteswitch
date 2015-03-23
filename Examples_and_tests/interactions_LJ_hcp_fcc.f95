  !
  ! By extracting lines beginning with the regular expression '\!\! ?'
  ! (ignoring leading whitespace), and then removing matches to the
  ! regular expression, html documentation corresponding to this 
  ! source code will be created.
  !
  !! <html>
  !! <head>
  !!  <title> Lennard-Jones hcp-fcc interactions Documentation </title>
  !! </head>
  !! <body>
  !! 
  !! <h1> Lennard-Jones hcp-fcc <code>interactions</code> Documentation </h1>
  !!
  !! <h2> Author </h2>
  !! <p> Tom Underwood </p>
  !!
  !! <h2> General Notes </h2>
  !!
  !! <h3> Description </h3>
  !! <p>
  !! This file contains variables and procedures pertaining to how the particles interact. It is
  !! <code>include</code>d in the <code>monteswitch_mod</code> module.
  !! </p>
  !! <p>
  !! This version of the file corresponds to Lennard-Jones interactions in hcp (lattice 1) and fcc (lattice 2).
  !! The interactions are not your typical Lennard-Jones interactions. Usually one truncates the interactions at
  !! a predetermined distance, or uses a predetermined 'list' of interacting pairs of particles. However, this
  !! gives errors in the ground state energies (see 'Structural Phase Behaviour Via Monte Carlo Techniques', 
  !! Andrew N. Jackson, PhD Thesis, University of Edinburgh, 2002). A more accurate approach is to evaluate the 
  !! <i>difference</i> of the energy of the lattice under consideration relative to the ground state, and apply 
  !! the usual truncations to this difference. This is what is done here. Of course, this approach necessitates 
  !! a, say, fcc-specific Fortran procedure to treat the fcc lattice; while in the conventional approach the 
  !! Lennard-Jones procedure can be applied to any crystal.
  !! </p>
  !! <p>
  !! The format of the file to import the interactions parameters is as follows. On the first line
  !! there are two tokens. The first is a <code>character(len=20)</code> variable (I recommend: 'lj_epsilon=');
  !! the second is the value of <code>lj_epsilon</code> (all variables will be explained in a moment). The second 
  !! line is similar, but for <code>lj_sigma</jcode>. The third line is for <code>lj_cutoff</code>, and the 
  !! fourth is for <code>list_cutoff</code>. 
  !! </p>



  !! <h2> 'Interation variables' which are required for initialisation of the particle interactions. </h3>
  !! <p>
  !! The following variables are read from the <code>filename_params</code> file by 
  !! <code>initialise_from_files(filename_params,filename_lattice)</code> in <code>monteswitch_mod</code>.
  !! These varaibles are later used to initialise how the particles interact.
  !! </p>
  !! <table border="1">
  !! <tr>
  !!  <td> <b> Variable </b> </td>
  !!  <td> <b> Type </b> </td>
  !!  <td> <b> Description </b> </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>lj_epsilon<code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> Depth of the Lennard-Jones potential well. </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>lj_sigma<code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> Distance corresponding to 0 potential for the Lennard-Jones potential. </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>lj_cutoff<code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> Cut-off distance for the Lennard-Jones potential. </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>list_cutoff<code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> 
  !!  Cut-off distance determining whether pairs of particles interact with each other throughout the simulation.
  !!  Those within <code>list_cutoff<code> of each other at the start of the simulation, before any moves are made,
  !!  will interact with each other. Note that the set of interacting pairs does not change during the simulation, even if
  !!  pairs of interacting particles later exceed <code>list_cutoff<code> in separation.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>list_size<code> </td>
  !!  <td> <code>integer(ik)</code> </td>
  !!  <td> 
  !!  The maximum number of particles any one particle is 'allowed' to interact with via the 'list' mechanism (determined
  !!  by the <code>list_cutoff<code> variables above). E.g., if <code>list_cutoff<code> is set to (just over) the nearest
  !!  neighbour distance, and there are 12 nearest neighbours, then <code>list_size<code> should be set to at least 12.
  !!  While there is nothing wrong in principle with setting <code>list_size<code> to, say, 500, the associated arrays would
  !!  be very large (500 integers per particle), and hence would slow down the simulation and/or use up too much memory.
  !!  </td>
  !! </tr>
  !! </table>
  real(rk) :: lj_epsilon
  real(rk) :: lj_sigma
  real(rk) :: lj_cutoff
  real(rk) :: list_cutoff
  integer(ik) :: list_size


  !! <h2> Further 'interaction variables'</h3>
  !! <p>
  !! The following variables are initialised according to the values of the above variables.
  !! </p>
  !! <table border="1">
  !! <tr>
  !!  <td> <b> Variable </b> </td>
  !!  <td> <b> Type </b> </td>
  !!  <td> <b> Description </b> </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>list_1<code> </td>
  !!  <td> <code>integer(ik), dimension(:,:), allocatable</code> </td>
  !!  <td>
  !!  Array holding information regarding which pairs of particles interact with one
  !!  another when we are in lattice 1: <code>list_1(i,n)</code> is the nth particle which
  !!  particle 'i' interacts with. <code>list_1(i,n)=0</code> for 'n' greater than 'z', where 
  !!  'z' is the number of particles which 'i' interacts with.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>list_2<code> </td>
  !!  <td> <code>integer(ik), dimension(:,:), allocatable</code> </td>
  !!  <td> 
  !!  Array holding information regarding which pairs of particles interact with one
  !!  another when we are in lattice 2. This is similar to <code>list_1<code>.
  !!  </td>
  !! </tr>
  !! </table>
  integer(ik), dimension(:,:), allocatable :: list_1
  integer(ik), dimension(:,:), allocatable :: list_2


contains


  !! <h2> Key procedures used in <code>monteswitch_mod</code> </h2>


  !! <h3> <code> subroutine import_interactions_params(unit) </code> </h3>
  !! <p>
  !! This procedure imports interaction variables required to initialise the interactions from the specified file,
  !! and initialises all interaction variables. This procedure presumably requires all other simulation variables
  !! to be initialised; accordingly it is the 'last' initialisation procedure to be called.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> filename </code> </td>
  !!   <td> <code>  character(*), intent(in) </code> </td>
  !!   <td> The file from which the interactions will be imported.
  !!   </td>
  !!  </tr>
  !! </table>
  subroutine initialise_interactions(filename)
    character(*), intent(in) :: filename
    character(len=20) string
    integer(ik) :: error
    integer(ik) :: i,j,n

    ! Open the file and import the variables required for initialisation
    open(unit=10,file=filename,iostat=error,status="old")
    if(error/=0) then
       write(0,*) "interactions: Error. Problem opening file '",filename,"'"
       stop 1
    end if
    read(10,*,iostat=error) string, lj_epsilon
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'lj_epsilon' from file '",filename,"'"
       stop 1
    end if
    read(10,*,iostat=error) string, lj_sigma
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'lj_sigma' from file '",filename,"'"
       stop 1
    end if
    read(10,*,iostat=error) string, lj_cutoff
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'lj_cutoff' from file '",filename,"'"
       stop 1
    end if
    read(10,*,iostat=error) string, list_cutoff
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'list_cutoff' from file '",filename,"'"
       stop 1
    end if
    read(10,*,iostat=error) string, list_size
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'list_size' from file '",filename,"'"
       stop 1
    end if
    close(10)

    ! Initialise the interaction variables accordingly

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
          if(min_image_distance(R_1(i,:),R_1(j,:),Lx(1),Ly(1),Lz(1))<list_cutoff) then
             list_1(i,n)=j
             n=n+1
          end if
       end do
    end do
    list_2=0
    do i=1,n_part
       n=1
       do j=1,n_part
          if(min_image_distance(R_2(i,:),R_2(j,:),Lx(2),Ly(2),Lz(2))<list_cutoff) then
             list_2(i,n)=j
             n=n+1
          end if
       end do
    end do

  end subroutine initialise_interactions




  !! <h3> <code> subroutine export_interactions_params(unit) </code> </h3>
  !! <p>
  !! This procedure exports all interaction variables to the specified unit.
  !! </p>
  subroutine export_interactions_state(unit)
    integer(ik), intent(in) :: unit
    write(unit,*) "lj_epsilon= ",lj_epsilon
    write(unit,*) "lj_sigma= ",lj_sigma
    write(unit,*) "lj_cutoff= ",lj_cutoff
    write(unit,*) "list_cutoff= ",list_cutoff
    write(unit,*) "list_size= ",list_size
    write(unit,*) "list_1= ",list_1
    write(unit,*) "list_2= ",list_2
  end subroutine export_interactions_state




  !! <h3> <code> subroutine import_interactions_params(unit) </code> </h3>
  !! <p>
  !! This procedure imports all interaction variables from the specified unit.
  !! </p>
  subroutine import_interactions_state(unit)
    integer(ik), intent(in) :: unit
    character(len=20) string
    integer(ik) :: error
    
    read(unit,*,iostat=error) string, lj_epsilon
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'lj_epsilon' from unit ",unit
       stop 1
    end if
    read(unit,*,iostat=error) string, lj_sigma
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'lj_sigma' from unit ",unit
       stop 1
    end if
    read(unit,*,iostat=error) string, lj_cutoff
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'lj_cutoff' from unit ",unit
       stop 1
    end if
    read(unit,*,iostat=error) string, list_cutoff
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'list_cutoff' from unit ",unit
       stop 1
    end if
    read(unit,*,iostat=error) string, list_size
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'list_size' from unit ",unit
       stop 1
    end if

    if(allocated(list_1)) then
       deallocate(list_1)
    end if
    if(allocated(list_2)) then
       deallocate(list_2)
    end if
    allocate(list_1(n_part,list_size))
    allocate(list_2(n_part,list_size))

    read(unit,*,iostat=error) string, list_1
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'list_1' from unit ",unit
       stop 1
    end if
    read(unit,*,iostat=error) string, list_2
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'list_2' from unit ",unit
       stop 1
    end if

  end subroutine import_interactions_state




  !! <h3> <code> function calc_energy_scratch(lattice,Lx,Ly,Lz,r) </code> </h3>
  !! <p>
  !! This function calculates the energy of a given system 'from scratch'.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> lattice </code> </td>
  !!   <td> <code> integer(ik), intent(in) </code> </td>
  !!   <td>
  !!   Lattice type to calculate the energy for.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> Lx_in </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   Length of supercell in x direction.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> Ly_in </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   Length of supercell in y direction.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> Lz_in </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   Length of supercell in z direction.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> r </code> </td>
  !!   <td> <code> real(rk), intent(in), dimension(:,:) </code> </td>
  !!   <td>
  !!   Array containing the particle positions. These positions will be within the supercell 'box'.  
  !!   </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
  function calc_energy_scratch(lattice,Lx_in,Ly_in,Lz_in,r)
    integer(ik), intent(in) :: lattice
    real(rk), intent(in) :: Lx_in
    real(rk), intent(in) :: Ly_in
    real(rk), intent(in) :: Lz_in
    real(rk), intent(in), dimension(:,:) :: r
    real(rk) :: calc_energy_scratch
    ! Constants for calculating the exact ground state for the hcp and fcc lattices
    ! at a given density (see Table 6.1 of Andrew Jackson's thesis mentioned above)
    real(rk), parameter :: A12_hcp=12.132293768
    real(rk), parameter :: A6_hcp=14.454897093
    real(rk), parameter :: A12_fcc=12.131880196
    real(rk), parameter :: A6_fcc=14.453920885
    real(rk) :: E_gs
    ! A scalefactor (explained below)
    real(rk) :: scale

    select case(lattice)
    case(1)
       ! Calculate the ground state energy (for 0 particle displacements) for the volume of the input system (corresponding to Lx_in,
       ! Ly_in, Lz_in, r) (hcp)
       E_gs=2*n_part*lj_epsilon*( (lj_sigma**3*n_part/(Lx_in*Ly_in*Lz_in*sqrt(2.0_rk)))**4*A12_hcp &
            -(lj_sigma**3*n_part/(Lx_in*Ly_in*Lz_in*sqrt(2.0_rk)))**2*A6_hcp )
       ! Calculate the difference relative to this ground state using truncation (1st 2 terms), and add it to the
       ! ground state for the input volume. Note that the 2nd term is the energy of the ground state with a truncated
       ! potential, and that to get this we need to rescale R_1, Lx_in, Ly_in and Lz_in such that the underlying
       ! undistorted lattice corresponds to the INPUT volume to this function, not the CURRENT volume which R_1, Lx(1), Ly(1)
       ! and Lz(1) correspond to.
       scale=(Lx_in*Ly_in*Lz_in/(Lx(1)*Ly(1)*Lz(1)))**(1.0_rk/3.0_rk)
       calc_energy_scratch = lj_energy_trunc_list_2(lj_epsilon,lj_sigma,lj_cutoff,Lx_in,Ly_in,Lz_in,list_1,r) &
         -lj_energy_trunc_list_2(lj_epsilon,lj_sigma,lj_cutoff,Lx(1)*scale,Ly(1)*scale,Lz(1)*scale,list_1,R_1*scale) + E_gs
    case(2)
       ! As above, but for fcc
       E_gs=2*n_part*lj_epsilon*( (lj_sigma**3*n_part/(Lx_in*Ly_in*Lz_in*sqrt(2.0_rk)))**4*A12_fcc &
            -(lj_sigma**3*n_part/(Lx_in*Ly_in*Lz_in*sqrt(2.0_rk)))**2*A6_fcc )
       !
       scale=(Lx_in*Ly_in*Lz_in/(Lx(2)*Ly(2)*Lz(2)))**(1.0_rk/3.0_rk)
       calc_energy_scratch = lj_energy_trunc_list_2(lj_epsilon,lj_sigma,lj_cutoff,Lx_in,Ly_in,Lz_in,list_2,r) &
         -lj_energy_trunc_list_2(lj_epsilon,lj_sigma,lj_cutoff,Lx(2)*scale,Ly(2)*scale,Lz(2)*scale,list_2,R_2*scale) + E_gs
    case default
       write(0,*) "interactions: Error. 'lattice' is not 1 or 2."
       stop 1
    end select

  end function calc_energy_scratch




  !! <h3> <code> function calc_energy_part_move(lattice,Lx,Ly,Lz,r,r_new,i) </code> </h3>
  !! <p>
  !! This function calculates the energy <i>change</i> of the system if particle <code>i</code> is moved
  !! such that the positions, which were previously <code>r</code>, are now <code>r_new</code>.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> lattice </code> </td>
  !!   <td> <code> integer(ik), intent(in) </code> </td>
  !!   <td>
  !!   Lattice type to calculate the energy for.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> Lx </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   Length of supercell in x direction.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> Ly </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   Length of supercell in y direction.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> Lz </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td>
  !!   Length of supercell in z direction.
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> r </code> </td>
  !!   <td> <code> real(rk), intent(in), dimension(:,:) </code> </td>
  !!   <td>
  !!   Array containing the particle positions. These positions will be within the supercell 'box'.  
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> r_new </code> </td>
  !!   <td> <code> real(rk), intent(in), dimension(:,:) </code> </td>
  !!   <td>
  !!   Array containing the new particle positions with the positions of particle <code>i</code> moved. 
  !!   These positions will be within the supercell 'box'.  
  !!   </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> i </code> </td>
  !!   <td> <code> integer(ik) </code> </td>
  !!   <td>
  !!   The number of the particle which has been moved
  !!   </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
  function calc_energy_part_move(lattice,Lx,Ly,Lz,r,r_new,i)
    integer(ik), intent(in) :: lattice
    real(rk), intent(in) :: Lx
    real(rk), intent(in) :: Ly
    real(rk), intent(in) :: Lz
    real(rk), intent(in), dimension(:,:) :: r
    real(rk), intent(in), dimension(:,:) :: r_new
    integer(ik), intent(in) :: i
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
             sep = min_image_distance(r(i,:),r(j,:),Lx,Ly,Lz)
             sep_new = min_image_distance(r_new(i,:),r_new(j,:),Lx,Ly,Lz)
             calc_energy_part_move = calc_energy_part_move + lj_pot_trunc(lj_epsilon,lj_sigma,lj_cutoff,sep_new) &
                  - lj_pot_trunc(lj_epsilon,lj_sigma,lj_cutoff,sep)
          end if
       end if
       n=n+1
    end do


  end function calc_energy_part_move




  !! <h2>Procedures used internally </h2>



  !! <h3> <code> function min_image_distance(r_1,r_2,Lx,Ly,Lz) </code> </h3>
  !! <p>
  !! <code>  min_image_distance </code>  returns the distance between the 
  !! positions <code>r_1</code> and <code>r_2</code> according to the 
  !! minimum image convention for a periodic cuboid whose faces are x=0, 
  !! x=<code>Lx</code>, y=0, y=<code>Ly</code>, z=0, and z=<code>Lz</code>.
  !! <code>r_1(1)</code> is the x-component of <code>r_1</code>, 
  !! <code>r_1(2)</code> is the y-component, and <code>r_1(3)</code> is the 
  !! z-component; and similarly for <code>r_2</code>. Note that
  !! <code>r_1</code> and <code>r_2</code> must be such that
  !! 0<=x<<code>Lx</code>, 0<=y<<code>Ly</code>, and 0<=z<<code>Lz</code>.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> r_1 </code> </td>
  !!   <td> <code> real(rk), dimension(3), intent(in) </code> </td>
  !!   <td> Position within the cuboid. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> r_2 </code> </td>
  !!   <td> <code> real(rk), dimension(3), intent(in) </code> </td>
  !!   <td> Position within the cuboid. </td>
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
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
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




   !! <h3> <code> function lj_pot(epsilon,sigma,r) </code> </h3>
  !! <p>
  !! This function returns the value of the Lennard-Jones potential at the specified
  !! distance, given the specified parametrisation of the potential.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> epsilon </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Depth of the potential well. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> sigma </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Distance corresponding to 0 potential. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> r </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Distance at which the potential is to be evaluated. </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
  function lj_pot(epsilon,sigma,r)
    real(rk), intent(in) :: epsilon
    real(rk), intent(in) :: sigma
    real(rk), intent(in) :: r
    real(rk) :: lj_pot
    lj_pot=4.0_rk*epsilon*( (sigma/r)**12-(sigma/r)**6 )
  end function lj_pot




  !! <h3> <code> function lj_pot_trunc(epsilon,sigma,cutoff,r) </code> </h3>
  !! <p>
  !! This function returns the value of a truncated Lennard-Jones potential at the specified
  !! distance, given the specified parametrisation of the potential. The potential is truncated
  !! at distance <code>cutoff</code>, and the potential is shifted to avoid a discontinuity in the
  !! potential at this distance.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> epsilon </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Depth of the potential well. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> sigma </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Distance corresponding to 0 potential. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> cutoff </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> The cut-off distance for the potential. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> r </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Distance at which the potential is to be evaluated. </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
  function lj_pot_trunc(epsilon,sigma,cutoff,r)
    real(rk), intent(in) :: epsilon
    real(rk), intent(in) :: sigma
    real(rk), intent(in) :: cutoff
    real(rk), intent(in) :: r
    real(rk) :: lj_pot_trunc
    if(r<cutoff) then
       lj_pot_trunc=lj_pot(epsilon,sigma,r)-lj_pot(epsilon,sigma,cutoff)
    else
       lj_pot_trunc=0.0_rk
    end if
  end function lj_pot_trunc




  !! <h3> <code> function lj_energy_trunc_list_2(epsilon,sigma,cutoff,Lx,Ly,Lz,list,r) </code> </h3>
  !! <p>
  !! This function returns the Lennard-Jones energy for a system whose particles are contained
  !! in the array <code>r</code> - where (<code>r(i,1)</code>,<code>r(i,2)</code>,<code>r(i,3)</code>) 
  !! is the position of the 'i'th particle - where the supercell has dimensions <code>Lx</code>, 
  !! <code>Ly</code> and <code>Lz</code> respectively in the x-, y- and z- directions. The positions
  !! of the particles must conform to the rules for the <code>min_image_distance</code> function.
  !! <code>list</code> provides information regarding the pairs of particles which will interact. The 
  !! (i,n)th element of <code>list</code> is the nth particle which particle 'i' interacts with. If 'i'
  !! only interacts with 'z' particles, then element (i,n) for n>z should be 0.
  !! For this procedure, the lower bounds on the array indices should be 1.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> epsilon </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Depth of the potential well. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> sigma </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Distance corresponding to 0 potential. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> Lx </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Length of supercell along the x-axis. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> Ly </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Length of supercell along the y-axis. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> Lz </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Length of supercell along the z-axis. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> cutoff </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> The cut-off distance for the potential. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> list </code> </td>
  !!   <td> <code> integer(ik), dimension(:,:), intent(in) </code> </td>
  !!   <td> Array determining which particles interact. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> r </code> </td>
  !!   <td> <code> real(rk), dimension(:,:), intent(in) </code> </td>
  !!   <td> Particle positions. </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
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
    real(rk) :: sep
    lj_energy_trunc_list_2=0.0_rk
    do i=1,ubound(r,1)
       n=1
       do
          j=list(i,n)
          if(j==0) then
             exit
          else
             if(j/=i) then
                sep=min_image_distance(r(i,:),r(j,:),Lx,Ly,Lz)
                lj_energy_trunc_list_2=lj_energy_trunc_list_2+lj_pot_trunc(epsilon,sigma,cutoff,sep)
             end if
          end if
          n=n+1
       end do
    end do
    lj_energy_trunc_list_2=0.5*lj_energy_trunc_list_2
  end function lj_energy_trunc_list_2
