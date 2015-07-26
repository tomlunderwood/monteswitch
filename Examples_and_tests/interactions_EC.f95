  !
  ! By extracting lines beginning with the regular expression '\!\! ?'
  ! (ignoring leading whitespace), and then removing matches to the
  ! regular expression, html documentation corresponding to this 
  ! source code will be created.
  !
  !! <html>
  !! <head>
  !!  <title> Einstein crystal interactions Documentation </title>
  !! </head>
  !! <body>
  !! 
  !! <h1> Einstein crystal <code>interactions</code> Documentation </h1>
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
  !! This version of the file corresponds to a pair of Einstein crystals: in lattice 1 the energy
  !! of each particle corresponds to a harmonic oscillator about the lattice site in which the
  !! spring constant is <code>alpha_1</code>; the spring constant for lattice 2 is <code>alpha_2</code>.
  !! </p>
  !! <p>
  !! The format of the file to import the interactions parameters is as follows. On the first line
  !! there are two tokens. The first is a <code>character(len=20)</code> variable (I recommend: 'alpha_1=');
  !! the second is the value of <code>alpha_1</code>. The second line is similar, but for <code>alpha_2</code>.
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
  !!  <td> <code>alpha_1<code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> Strength of the harmonic potential well centred on all lattice sites in phase 1. </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>alpha_2<code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> Strength of the harmonic potential well centred on all lattice sites in phase 2. </td>
  !! </tr>
  !! </table>
  real(rk) :: alpha_1
  real(rk) :: alpha_2


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
  !! </tr>
  !! </table>
  
  ! NOTHING REQUIRED FOR THE EINSTEIN CRYSTAL


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

    ! Open the file and import the variables required for initialisation
    open(unit=10,file=filename,iostat=error,status="old")
    if(error/=0) then
       write(0,*) "interactions: Error. Problem opening file '",filename,"'"
       stop 1
    end if
    read(10,*,iostat=error) string, alpha_1
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'alpha_1' from file '",filename,"'"
       stop 1
    end if
    read(10,*,iostat=error) string, alpha_2
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'alpha_2' from file '",filename,"'"
       stop 1
    end if
    close(10)

    ! Initialise the interaction variables accordingly

    ! NOTHING REQUIRED FOR EINSTEIN CRYSTAL

  end subroutine initialise_interactions




  !! <h3> <code> subroutine export_interactions_params(unit) </code> </h3>
  !! <p>
  !! This procedure exports all interaction variables to the specified unit.
  !! </p>
  subroutine export_interactions_state(unit)
    integer(ik), intent(in) :: unit
    write(unit,*) "alpha_1= ",alpha_1
    write(unit,*) "alpha_2= ",alpha_2
  end subroutine export_interactions_state




  !! <h3> <code> subroutine import_interactions_params(unit) </code> </h3>
  !! <p>
  !! This procedure imports all interaction variables from the specified unit.
  !! </p>
  subroutine import_interactions_state(unit)
    integer(ik), intent(in) :: unit
    character(len=20) string
    integer(ik) :: error
    read(unit,*,iostat=error) string, alpha_1
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'alpha_1' from unit ",unit
       stop 1
    end if
    read(unit,*,iostat=error) string, alpha_2
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'alpha_2' from unit ",unit
       stop 1
    end if
  end subroutine import_interactions_state




  !! <h3> <code> subroutine after_accepted_interactions() </code> </h3>
  !! <p>
  !! This procedure is called after a particle or volume move (but not a lattice move) has been
  !! accepted. It can be used to update, say, neighbour lists, which require knowledge of the
  !! current microstate - or in this case, the current microstate pertaining to both lattices.
  !! </p>
  subroutine after_accepted_interactions()
   
  end subroutine after_accepted_interactions




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
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
  function calc_energy_scratch(lattice,Lx,Ly,Lz,r)
    integer(ik), intent(in) :: lattice
    real(rk), intent(in) :: Lx
    real(rk), intent(in) :: Ly
    real(rk), intent(in) :: Lz
    real(rk), intent(in), dimension(:,:) :: r
    real(rk) :: calc_energy_scratch
    real(rk), dimension(3) :: pos, site
    integer(ik) :: i
    calc_energy_scratch=0.0_rk
    do i=1,n_part
       select case(lattice)
       case(1)
          pos=r(i,:)
          site=R_1(i,:)
          calc_energy_scratch=calc_energy_scratch+alpha_1*min_image_distance(pos,site,Lx,Ly,Lz)**2
       case(2)
          pos=r(i,:)
          site=R_2(i,:)
          calc_energy_scratch=calc_energy_scratch+alpha_2*min_image_distance(pos,site,Lx,Ly,Lz)**2
       case default
          write(0,*) "Error: interactions: 'lattice' was not 1 or 2."
          stop 1
       end select
    end do
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
    real(rk), dimension(3) :: pos, pos_new, site
    select case(lattice)
    case(1)
       pos=r(i,:)
       pos_new=r_new(i,:)
       site=R_1(i,:)
       calc_energy_part_move=alpha_1*min_image_distance(pos_new,site,Lx,Ly,Lz)**2-alpha_1*min_image_distance(pos,site,Lx,Ly,Lz)**2
    case(2)
       pos=r(i,:)
       pos_new=r_new(i,:)
       site=R_2(i,:)
       calc_energy_part_move=alpha_2*min_image_distance(pos_new,site,Lx,Ly,Lz)**2-alpha_2*min_image_distance(pos,site,Lx,Ly,Lz)**2
    case default
       write(0,*) "Error: interactions: 'lattice' was not 1 or 2."
       stop 1
    end select
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
