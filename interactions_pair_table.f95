  !
  ! By extracting lines beginning with the regular expression '\!\! ?'
  ! (ignoring leading whitespace), and then removing matches to the
  ! regular expression, html documentation corresponding to this 
  ! source code will be created.
  !
  !! <html>
  !! <head>
  !!  <title> `Pair table' interactions Documentation </title>
  !! </head>
  !! <body>
  !! 
  !! <h1> `Pair table' <code>interactions</code> Documentation </h1>
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
  !! This version of the file implements a tabulated pair potential. The format of the file which contains the
  !! table (the 'interactions_in' file read by monteswitch) is as follows:
  !! The first line is a comment line.
  !! The second line contains <code>Nr</code>, then <code>dr</code>, then <code>cutoff</code>, which are, respectively, the 
  !! number of elements in the table, the separation between points on the `distance grid', and the cut-off distance for the
  !! potential.
  !! The next line, and all subsequent lines, contain the <code>Nr</code> values which defines the pair potential:
  !! the <code>i</code>th value corresponds to the value of the potential for inter-particle separation <code>(i-1)*dr</code>.
  !! Note that the potential is implemented here using linear interpolation between tabulated values.
  !! <p>
  !! At initialisation a file 'phi.dat' is created which is corresponds to the potential read from the
  !! 'interactions_in' file. Furthermore, the maximum supported separation is <code>(Nr-1)*dr</code>, as opposed to 
  !! <code>Nr*dr</code>. Hence if the cut-off is greater than <code>(Nr-1)*dr</code>, then it is ammeded at initialisation to 
  !! be <code>(Nr-1)*dr</code>.
  !! </p>
  
  !! <h2> 'Interation variables' which are required for initialisation of the particle interactions. </h3>
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
  !!  <td> <code>Nr<code> </td>
  !!  <td> <code>integer(ik)</code> </td>
  !!  <td> 
  !!  The number of points on the 'separation grid'.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>dr<code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> 
  !!  The spacing between points on the 'separation grid'.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>cutoff<code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> 
  !!  The interaction cut-off: pairs of particles beyond this separation do not interact.
  !!  </td>
  !! </tr>
  !! </table>
  integer(ik) :: Nr
  real(rk) :: dr
  real(rk) :: cutoff


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
  !!  <td> <code>phi<code> </td>
  !!  <td> <code>real(rk), dimension(:), allocatable</code> </td>
  !!  <td>
  !!  Array defining the pair-potential function: <code>phi(i)</code> is the pair-potential for inter-particle separation 
  !!  <code>(i-1)*dr</code>.
  !!  </td>
  !! </tr>
  !! </table>
  real(rk), dimension(:), allocatable :: phi

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
    character(len=200) :: line
    integer(ik) :: error
    integer(ik) :: i,j,n

    integer(ik) :: indx, prev_indx
    integer(ik) :: start_of_word

    ! Open the file and import the variables required for initialisation
    open(unit=10,file=filename,iostat=error,status="old")
    if(error/=0) then
       write(0,*) "interactions: Error. Problem opening file '",filename,"'"
       stop 1
    end if

    ! Read comment line
    read(10,*,iostat=error) line
    if(error/=0) then
       write(0,*) "interactions: Error reading file '",filename,"'"
       stop 1
    end if
    ! Read Nr, dr and cutoff
    read(10,*,iostat=error) Nr,dr,cutoff
    if(error/=0) then
       write(0,*) "interactions: Error reading file '",filename,"'"
       stop 1
    end if

    ! Initialise and read the pair potential function.
    allocate(phi(Nr))
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
             read(line(start_of_word:i-1), *) phi(j)
             j=j+1
          end if

          if(j>Nr) exit

          prev_indx=indx

       end do

       if(j>Nr) exit

    end do

    ! Set cutoff to the maximum supported value if it is greater than (Nr-1)*dr
    if(cutoff>(Nr-1)*dr) cutoff=(Nr-1)*dr

    ! Output the pair potential to a file
    open(unit=10,file="phi.dat")
    do i=1,Nr
       write(10,*) (i-1)*dr,phi(i)
    end do
    close(10)

  end subroutine initialise_interactions




  !! <h3> <code> subroutine export_interactions_params(unit) </code> </h3>
  !! <p>
  !! This procedure exports all interaction variables to the specified unit.
  !! </p>
  subroutine export_interactions_state(unit)
    integer(ik), intent(in) :: unit
    write(unit,*) "Nr= ",Nr
    write(unit,*) "dr= ",dr
    write(unit,*) "cutoff= ",cutoff
    write(unit,*) "phi= ",phi
  end subroutine export_interactions_state




  !! <h3> <code> subroutine import_interactions_params(unit) </code> </h3>
  !! <p>
  !! This procedure imports all interaction variables from the specified unit.
  !! </p>
  subroutine import_interactions_state(unit)
    integer(ik), intent(in) :: unit
    character(len=20) string
    integer(ik) :: error

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

    if(allocated(phi)) then
       deallocate(phi)
    end if
    allocate(phi(Nr))

    read(unit,*,iostat=error) string, phi
    if(error/=0) then
       write(0,*) "interactions: Error. Problem reading 'phi' from unit ",unit
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
    integer(ik) :: i,j
    real(rk) :: sep
    calc_energy_scratch=0.0_rk
    do i=1,n_part
       do j=i+1,n_part
          sep=min_image_distance(r(i,:),r(j,:),Lx_in,Ly_in,Lz_in)
          if(sep<cutoff) then
             calc_energy_scratch=calc_energy_scratch+phi_func(sep)
          end if
       end do
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
    integer(ik) :: j
    real(rk) :: sep
    real(rk) :: sep_new 

    ! Calculate the change in the pair energy
    calc_energy_part_move=0.0_rk
    do j=1,n_part
       if(j/=i) then
          sep=min_image_distance(r(j,:),r(i,:),Lx,Ly,Lz)
          sep_new=min_image_distance(r_new(j,:),r_new(i,:),Lx,Ly,Lz)
          if(sep_new<cutoff) calc_energy_part_move=calc_energy_part_move+phi_func(sep_new)
          if(sep<cutoff) calc_energy_part_move=calc_energy_part_move-phi_func(sep)
       end if
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




  !! <h3> <code> function phi_func(rho) </code> </h3>
  !! <p>
  !! <code>phi_func(r) </code> returns the value of the pair-potential corresponding to a given separation,
  !! given the array <code>phi</code>, and the 'separation grid' defined by <code>Nr/code>
  !! and <code>dr</code>. Note that linear interpolation is used.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code>r </code> </td>
  !!   <td> <code> real(rk) </code> </td>
  !!   <td> The distance between two particles. </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
  function phi_func(r)
    real(rk), intent(in) :: r
    real(rk) :: phi_func
    real(rk) :: bin_cont
    integer(ik) :: bin_below, bin_above
    ! phi(bin) corresponds to the value of phi at r=(bin-1)*dr
    ! Hence bin=r/dr+1 is the 'continuous' bin number
    bin_cont=r/dr+1
    if(bin_cont<1 .or. bin_cont>Nr) then
       write(0,*) "interactions_mod: Error in rho_func: r not supported. r = ",r
       stop 1
    end if
    ! Perform a linear interpolation to get the value of the phi function at r
    bin_below=floor(bin_cont)
    bin_above=ceiling(bin_cont)
    phi_func=phi(bin_below) + (phi(bin_above)-phi(bin_below))*(bin_cont-bin_below)
  end function phi_func

