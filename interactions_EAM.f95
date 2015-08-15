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
  !! This version of the file implements the embedded atom model (EAM) similarly to LAMMPS, as described at
  !! http://lammps.sandia.gov/doc/pair_eam.html. The input file which specifies the potential (the 'interactions_in'
  !! file read by monteswitch) must be in the DYNAMO/LAMMPS 'setfl' format; a file with this format is often
  !! indicated with the suffix '.eam.alloy'. A description of the 'setfl' format can be found at 
  !! http://lammps.sandia.gov/doc/pair_eam.html. However, note that, since monteswitch currently
  !! only supports unary systems, the file must not correspond to a multicomponent system - an error is returned
  !! if this is the case. Note also this file implements the EAM potential by using linear interpolation, using the
  !! tabulations of the various functions in the input file (the embedding function, the density function, 
  !! and the pair-potential) as a basis. At initialisation three files are created, 'F.dat', 'rho.dat', and 'rphi.dat', which
  !! correspond to the embedding function, density function, and pair potential (multiplied by separation) read
  !! from the input file.
  !! </p>
  !! <p>
  !! NB: It seems that 'setfl' format files often have cut-offs which are very slightly
  !! higher than that suggested by the values of <code>Nr</code> and
  !! <code>dr</code>. The aforementioned LAMMPS documentation specifies that
  !! elements <code>i</code> of <code>rho</code> and <code>rphi</code>
  !! correspond to a separation of <code>r=(i-1)*dr</code> - as opposed to
  !! <code>r=i*dr</code>. Hence the maximum supported separation is
  !! <code>(Nr-1)*dr</code>, as opposed to <code>Nr*dr</code>. However, people
  !! often set the cut-off to the latter value, as opposed to the former. Here,
  !! to account for this, if the cut-off is greater than <code>(Nr-1)*dr</code>, 
  !! then it is ammeded at initialisation to be <code>(Nr-1)*dr</code>.
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
  !!  <td> <code>Nrho<code> </td>
  !!  <td> <code>integer(ik)</code> </td>
  !!  <td> The number of points on the 'density grid'. </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>drho<code> </td>
  !!  <td> <code>real(rk)</code> </td>
  !!  <td> The spacing between points on the 'density grid'. </td>
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
  integer(ik) :: Nrho
  real(rk) :: drho
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
  !!  <td> <code>F<code> </td>
  !!  <td> <code>real(rk), dimension(:), allocatable</code> </td>
  !!  <td>
  !!  Array defining the embedding function: <code>F(i)</code> is the embedding function at density
  !!  <code>(i-1)*drho</code>.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>rho<code> </td>
  !!  <td> <code>real(rk), dimension(:), allocatable</code> </td>
  !!  <td>
  !!  Array defining the density function: <code>rho(i)</code> is the density for inter-particle separation
  !!  <code>(i-1)*dr</code>.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>rphi<code> </td>
  !!  <td> <code>real(rk), dimension(:), allocatable</code> </td>
  !!  <td>
  !!  Array defining the pair-potential function: <code>rphi(i)</code> is the pair-potential, multiplied by
  !!  the separation, for inter-particle separation <code>(i-1)*dr</code>.
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>part_rho<code> </td>
  !!  <td> <code>real(rk), dimension(:,:), allocatable</code> </td>
  !!  <td>
  !!  Array defining containing the density - with regards to the embedding function - associated with
  !!  each particle in the current microstate of the system for each latice. 
  !!  <code>part_rho(lattice,i)</code> is the density for particle <code>i</code> in lattice <cdoe>lattice</code>. 
  !!  This is used to speed up the calculation of the change in embedding energy as a result of a particle move. 
  !!  </td>
  !! </tr>
  !! <tr>
  !!  <td> <code>part_rho_buffer<code> </td>
  !!  <td> <code>real(rk), dimension(:,:), allocatable</code> </td>
  !!  <td>
  !!  Like <code>part_rho</code>, but a 'buffer' value. Every particle or volume move, the energy of the trial
  !!  microstate is evaluated by calling the <code>calc_energy_scratch</code> or <code>calc_energy_part_move</code>
  !!  procedure. In each of these procedures, the densities for the trial state are stored in <code>part_rho_buffer<code>.
  !!  If the move is accepted, i.e., the trial microstate becomes the actual microstate of the system, 
  !!  <code>part_rho_buffer<code> is copied to <code>part_rho</code>. If the move is rejected, then 
  !!  <code>part_rho_buffer<code> is not copied to <code>part_rho</code>.  Thus <code>part_rho</code> always reflects the 
  !!  actual microstate.
  !!  </td>
  !! </tr>
  !! </table>
  real(rk), dimension(:), allocatable :: F
  real(rk), dimension(:), allocatable :: rho
  real(rk), dimension(:), allocatable :: rphi
  real(rk), dimension(:,:), allocatable :: part_rho
  real(rk), dimension(:,:), allocatable :: part_rho_buffer

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
    call set_part_rho(1,Lx(1),Ly(1),Lz(1),R_1)
    call set_part_rho(2,Lx(2),Ly(2),Lz(2),R_2)

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
    write(unit,*) "part_rho= ",part_rho
    write(unit,*) "part_rho_buffer= ",part_rho_buffer
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

  end subroutine import_interactions_state



  
  !! <h3> <code> subroutine after_accepted_interactions() </code> </h3>
  !! <p>
  !! This procedure is called after a particle or volume move (but not a lattice move) has been
  !! accepted. It can be used to update, say, neighbour lists, which require knowledge of the
  !! current microstate - or in this case, the current microstate pertaining to both lattices.
  !! </p>
  subroutine after_accepted_interactions()
    part_rho=part_rho_buffer
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
    real(rk) :: pair_energy
    real(rk) :: embedding_energy
    pair_energy=0.0_rk
    embedding_energy=0.0_rk
    do i=1,n_part
       ! Calculate the density associated with particle i, and store it in part_rho_buffer(lattice,i).
       ! At the same time we are resetting part_rho_buffer for 'lattice' to correspond to the 
       ! microstate in the argument
       part_rho_buffer(lattice,i)=0.0_rk
       do j=1,n_part
          if(j/=i) then
             sep=min_image_distance(r(i,:),r(j,:),Lx_in,Ly_in,Lz_in)
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
    real(rk) :: delta_pair_energy
    real(rk) :: delta_embedding_energy
    ! The change in density for microstate with r_new relative to the current microstate (r)
    real(rk), dimension(n_part) :: delta_rho    

    ! Calculate the change in the pair energy
    delta_pair_energy=0.0_rk
    do j=1,n_part
       if(j/=i) then
          sep=min_image_distance(r(j,:),r(i,:),Lx,Ly,Lz)
          sep_new=min_image_distance(r_new(j,:),r_new(i,:),Lx,Ly,Lz)
          if(sep_new<cutoff) delta_pair_energy=delta_pair_energy+phi_func(sep_new)
          if(sep<cutoff) delta_pair_energy=delta_pair_energy-phi_func(sep)
       end if
    end do

    ! Calculate the change in the embedding energy

    ! Calculate the changes in all particles' densities for  r_new relative to r_ref
    delta_rho=0.0_rk
    do j=1,n_part
       if(j/=i) then
          sep=min_image_distance(r(j,:),r(i,:),Lx,Ly,Lz)
          sep_new=min_image_distance(r_new(j,:),r_new(i,:),Lx,Ly,Lz)
          if(sep_new<cutoff) then
             delta_rho(j)=delta_rho(j)+rho_func(sep_new)
             delta_rho(i)=delta_rho(i)+rho_func(sep_new)
          end if
          if(sep<cutoff) then
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




  !! <h2>Procedures used internally </h2>


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




  !! <h3> <code> function F_func(rho) </code> </h3>
  !! <p>
  !! <code>F_func(rho) </code> returns the embedding energy for a given density,
  !! given the array <code>F</code>, and the 'density grid' defined by <code>Nrho</code>
  !! and <code>drho</code>. Note that linear interpolation is used.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> rho </code> </td>
  !!   <td> <code> real(rk) </code> </td>
  !!   <td> The electron density. </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
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




  !! <h3> <code> function rho_func(rho) </code> </h3>
  !! <p>
  !! <code>rho_func(r) </code> returns the density corresponding to a given separation,
  !! given the array <code>rho</code>, and the 'separation grid' defined by <code>Nr/code>
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




  !! <h3> <code> function phi_func(rho) </code> </h3>
  !! <p>
  !! <code>phi_func(r) </code> returns the value of the pair-potential corresponding to a given separation,
  !! given the array <code>rphi</code>, and the 'separation grid' defined by <code>Nr/code>
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

