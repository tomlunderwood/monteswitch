  !
  ! By extracting lines beginning with the regular expression '\!\! ?'
  ! (ignoring leading whitespace), and then removing matches to the
  ! regular expression, html documentation corresponding to this 
  ! source code will be created.
  !
  !! <html>
  !! <head>
  !!  <title> EAM interactions Documentation </title>
  !! </head>
  !! <body>
  !! 
  !! <h1> EAM <code>interactions</code> Documentation </h1>
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
  !! This version of the file corresponds to the embedded atom model.
  !! </p>
  !! <p>
  !! The format of the file to import the interactions parameters is as follows... (DYNAMO multi-element setfl format,
  !! but only one element must be present).
  !! </p>



  !! <h2> 'Interation variables' which are required for initialisation of the particle interactions. </h3>
  !! <p>
  !! The following variables are... 
  !! </p>
  integer(ik) :: Nrho
  real(rk) :: drho
  integer(ik) :: Nr
  real(rk) :: dr
  real(rk) :: cutoff


  !! <h2> Further 'interaction variables'</h3>
  !! <p>
  !! The following variables are initialised according to the values of the above variables...
  !! </p>
  ! Embedding function (function of density rho)
  real(rk), dimension(:), allocatable :: F
  ! Density function (function of separation r)
  real(rk), dimension(:), allocatable :: rho
  ! Pair potential function times r (function of separation r)
  real(rk), dimension(:), allocatable :: rphi


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
    character(len=20) :: string
    integer(ik) :: error
    integer(ik) :: i,j,n

    ! Open the file and import the variables required for initialisation
    open(unit=10,file=filename,iostat=error,status="old")
    if(error/=0) then
       write(0,*) "interactions: Error. Problem opening file '",filename,"'"
       stop 1
    end if
    
    ! Read 3 comment lines
    read(10,*) string
    read(10,*) string
    read(10,*) string
    ! Read the number of elements and the element name (number of elements must be 1)
    read(10,iostat=error) n,string
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading line 4 from file '",filename,"'"
       stop 1
    end if
    if(n/=1) then
       write(0,*) "interactions: Error. File '",filename,"' seems to correspond to a ", &
            "multicomponent system; only unary systems supported."
       stop 1
    end if
    ! Read Nrho, drho, Nr, dr and cutoff
    read(10,iostat=error) Nrho,drho,Nr,dr,cutoff
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading line 5 (Nrho, drho, Nr, dr, ", &
            "cutoff) from file '",filename,"'"
       stop 1
    end if
    ! Read, but ignore, the line containing the atomic number, mass, lattice constant and lattice type
    read(10,iostat=error) string
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading line 6 (atomic number, mass, ", &
            "lattice constant, lattice type) from file '",filename,"'"
       stop 1
    end if
    ! Initialise and read the embedding function
    allocate(F(Nrho))
    read(10,iostat=error) F
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading embedding function from file '",filename,"'"
       stop 1
    end if
    ! Initialise and read the density function
    allocate(rho(Nrho))
    read(10,iostat=error) rho
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading density function from file '",filename,"'"
       stop 1
    end if
    ! Initialise and read the pair potential function. Note that the values in the file are actually r*phi(r),
    ! not phi(r). This is stored in the variable rphi.
    allocate(rphi(Nrho))
    read(10,iostat=error) rphi
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading pair-potential function from file '",filename,"'"
       stop 1
    end if
    close(10)

    ! Initialise the interaction variables accordingly

  end subroutine initialise_interactions




  !! <h3> <code> subroutine export_interactions_params(unit) </code> </h3>
  !! <p>
  !! This procedure exports all interaction variables to the specified unit.
  !! </p>
  subroutine export_interactions_state(unit)
    integer(ik), intent(in) :: unit
    write(unit,*) "Nrho= ",Nrho
    write(unit,*) "drho= ",drho
    write(unit,*) "Nr= ",Nr
    write(unit,*) "dr= ",dr
    write(unit,*) "cutoff= ",cutoff
    write(unit,*) "F= ",F
    write(unit,*) "rho= ",rho
    write(unit,*) "rphi= ",rphi
  end subroutine export_interactions_state




  !! <h3> <code> subroutine import_interactions_params(unit) </code> </h3>
  !! <p>
  !! This procedure imports all interaction variables from the specified unit.
  !! </p>
  subroutine import_interactions_state(unit)
    integer(ik), intent(in) :: unit
    character(len=20) string
    integer(ik) :: error
    
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
    integer(ik) :: i
    calc_energy_scratch=0.0_rk
    do i=1,ubound(r,1)
       calc_energy_scratch=calc_energy_scratch+particle_energy(Lx_in,Ly_in,Lz_in,r,i)
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
    calc_energy_part_move=particle_energy(Lx,Ly,Lz,r_new,i)-particle_energy(Lx,Ly,Lz,r,i)
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




  ! The embedding function for a given rho, according to the array F, and the
  ! 'rho-grid' defined by Nrho and drho
  function F_func(rho)
    real(rk), intent(in) :: rho
    real(rk) :: F_func
    real(rk) :: bin_cont
    integer(ik) :: bin_below, bin_above
    ! F(bin) corresponds to the embedding function at rho=bin*drho
    ! Get the 'continuous' bin number
    bin_cont=rho/drho
    if(bin_cont<1 .or. bin_cont>Nrho) then
       write(0,*) "interactions_mod: Error in F_func: r not supported"
       stop 1
    end if
    ! Perform a linear interpolation to get the value of the embedding function at rho
    bin_below=floor(bin_cont)
    bin_above=ceiling(bin_cont)
    F_func=F(bin_below) + (F(bin_above)-F(bin_below))*(bin_cont-bin_below)
  end function F_func


  ! The density function for a given r, according to the array rho, and the
  ! 'r-grid' defined by Nr and dr
  function rho_func(r)
    real(rk), intent(in) :: r
    real(rk) :: rho_func
    real(rk) :: bin_cont
    integer(ik) :: bin_below, bin_above
    ! rho(bin) corresponds to the density function at r=bin*dr
    ! Get the 'continuous' bin number
    bin_cont=r/dr
    if(bin_cont<1 .or. bin_cont>Nr) then
       write(0,*) "interactions_mod: Error in rho_func: r not supported"
       stop 1
    end if
    ! Perform a linear interpolation to get the value of the density function at r
    bin_below=floor(bin_cont)
    bin_above=ceiling(bin_cont)
    rho_func=rho(bin_below) + (rho(bin_above)-rho(bin_below))*(bin_cont-bin_below)
  end function rho_func


  ! The phi function for a given r, according to the array rphi, and the
  ! 'r-grid' defined by Nr and dr
  function phi_func(r)
    real(rk), intent(in) :: r
    real(rk) :: phi_func
    real(rk) :: bin_cont
    integer(ik) :: bin_below, bin_above
    ! rphi(bin) corresponds to the value of r*phi at r=bin*dr
    ! Get the 'continuous' bin number
    bin_cont=r/dr
    if(bin_cont<1 .or. bin_cont>Nr) then
       write(0,*) "interactions_mod: Error in rho_func: r not supported"
       stop 1
    end if
    ! Perform a linear interpolation to get the value of the phi function at r
    bin_below=floor(bin_cont)
    bin_above=ceiling(bin_cont)
    phi_func=rphi(bin_below) + (rphi(bin_above)-rphi(bin_below))*(bin_cont-bin_below)
    phi_func=phi_func/r
  end function phi_func


  ! The total energy associated with atom i, given the particle positions and supercell
  ! dimensions
  function particle_energy(Lx,Ly,Lz,r,i)
    real(rk), intent(in) :: Lx
    real(rk), intent(in) :: Ly
    real(rk), intent(in) :: Lz
    real(rk), intent(in), dimension(:,:) :: r
    integer(ik), intent(in) :: i
    real(rk) :: particle_energy
    integer(ik) :: j
    real(rk) :: sep
    real(rk) :: rho_sum
    real(rk) :: phi_sum
    rho_sum=0.0_rk
    phi_sum=0.0_rk
    do j=1,ubound(r,1)
       if(j/=i) then
          sep=min_image_distance(r(i,:),r(j,:),Lx,Ly,Lz)
          if(sep<cutoff) then
             rho_sum=rho_sum+rho_func(sep)
             phi_sum=phi_sum+phi_func(sep)
          end if
       end if
    end do
    particle_energy=F_func(rho_sum)+0.5_rk*phi_sum    
  end function particle_energy
  
