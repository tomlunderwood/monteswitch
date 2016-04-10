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
! By extracting lines beginning with the regular expression '\!\! ?'
! (ignoring leading whitespace), and then removing matches to the
! regular expression, html documentation corresponding to this 
! source code will be created.
!
!! <html>
!! <head>
!!  <title> monteswitch_post Documentation </title>
!! </head>
!! <body>
!! 
!! <h1> <code>monteswitch_post</code> Documentation </h1>
!!
!! <h2> Author </h2>
!! <p> Tom Underwood </p>
!!
!! <h2> General Notes </h2>
!! <p>
!! This program performs various post-processing operations on the 'state' file generated by <code>monteswitch</code> programs.
!! </p>
!! <p>
!! The command line options are as follows. Arguments in parenthesis are optional.
!! <table border="1">
!!  <tr>
!!   <td> <b> Argument </b> </td>
!!   <td> <b> Description </b> </td>
!!  </tr>
!!  <tr>
!!   <td> -extract_wf </td>
!!   <td>
!!   Extract the weight function from the 'state' file and output it to stdout. In the output the 
!!   first token on each line is the order parameter, and the second is the corresponding value of the weight
!!   function.
!!   </td>
!!  </tr>
!!  <tr>
!!   <td> -extract_M_counts </td>
!!   <td>
!!   Extract order parameter histograms from the 'state' file, and output them to stdout. In the output
!!   the first token on each line is the order parameter, the second is the corresponding number of
!!   counts for lattice type 1, and the third is the corresponding number of counts for lattice type 2.
!!   </td>
!!  </tr>
!!  <tr>
!!   <td> -extract_pos [<i>species</i>] </td>
!!   <td>
!!   Extract the positions of the particles, and output them to stdout. In the output
!!   the first, second and third tokens on each line are the x-, y- and z-coordinates respectively
!!   for a particle. <i>species</i> is an optional argument. If absent then particle positions for all species are output;
!!   if specified then only positions for particles belonging to species <i>species</i> are output.
!!   </td>
!!  </tr>
!!   <td> -extract_R_1 [<i>species</i>] </td>
!!   <td>
!!   Extract the positions of the lattice vectors for lattice type 1, and output them to stdout. In 
!!   the output the first, second and third tokens on each line are the x-, y- and z-coordinates respectively
!!   for a particle. <i>species</i> is an optional argument. If absent then vectors for particles of all species are output;
!!   if specified then only vectors for particles belonging to species <i>species</i> in lattice 1 are output.
!!   </td>
!!  </tr>
!!  </tr>
!!   <td> -extract_R_2 [<i>species</i>] </td>
!!   <td>
!!   Extract the positions of the lattice vectors for lattice type 2, and output them to stdout. In 
!!   the output the first, second and third tokens on each line are the x-, y- and z-coordinates respectively
!!   for a particle. <i>species</i> is an optional argument. If absent then vectors for particles of all species are output;
!!   if specified then only vectors for particles belonging to species <i>species</i> in lattice 2 are output.
!!   </td>
!!  </tr>
!!  <tr>
!!   <td> -extract_u [<i>species</i> <i>lattice</i>] </td>
!!   <td>
!!   Extract the displacements of the particles, and output them to stdout. In the output
!!   the first, second and third tokens on each line are the x-, y- and z-displacements respectively
!!   for a particle. <i>species</i> and <i>lattice</i> are optional arguments. If absent then displacements for particles of
!!   all species are output; if specified then only displacements for particles belonging to species 
!!   <i>species</i> in lattice <i>lattice</i> are output.
!!   </td>
!!  </tr>
!!  <tr>
!!   <td> -calc_rad_dist <i>bins</i> </td>
!!   <td>
!!   Calculate the instantaneous radial distribution function, which is output to stdout. In the output, each line corresponds
!!   to a distance, which is the first token; the second token contains the average number of particles at this distance
!!   from a particle. The output is like a histogram, with each line corresponding to a bin, the first token 
!!   corresponding to the minimum of the range covered by the bin, and the second token corresponding to the number of 
!!   counts for that bin. The range of the bin is inclusive at its minimum, and exclusive at its maximum. The range 
!!   chosen to evaluate the radial distribution function for is <code>min(Lx/2,Ly/2,Lz/2)</code>, where <code>Lx</code>, 
!!   <code>Ly</code> and <code>Lz</code> pertain to the lattice type of the system in 'state', and the number of bins 
!!   for the histogram is specified in the argument <i>bins</i>. Note that the range of distances corresponds to 
!!   the 'limit of periodicity' for the system.
!!   </td>
!!  </tr>
!!  <tr>
!!   <td> -merge_trans <i>state_in_1</i> <i>state_in_2</i> <i>state_out</i> </td>
!!   <td>
!!   Combine the <code>trans</code> matrices from the files <i>state_in_1</i> and <i>state_in_2</i>, and store the
!!   result in the file <i>state_out</i>, where all variables in <i>state_out</i> other than the matrix <code>trans</code>
!!   are inherited from <i>state_in_1</i>. Note that <code>M_grid</code> <i>must</i> be the same for both <i>state_in_1</i> and
!!   <i>state_in_2</i> (in which case the matrices <code>tran</code> are of the same size).
!!   </td>
!!  </tr>
!!  <tr>
!!   <td> -extract_lattices_in <i>vectors_in_1</i> <i>vectors_in_2</i></td>
!!   <td>
!!   Output the geometrical properties of the system in the format of a 'lattices_in' file. The arguments <i>vectors_in_1</i> 
!!   and <i>vectors_in_2</i> can be either <code>pos</code> or <code>R</code>. If <i>vectors_in_1</i> is <code>pos</code>, then
!!   the positions of the particles in lattice 1 of the 'lattices_in' file will be the positions of the particles in lattice 1
!!   in the 'state' file; if <i>vectors_in_1</i> is <code>R</code>, then the positions of the particles in lattice 1 of the 
!!   'lattices_in' file will be the lattice vectors (i.e., <code>R_1</code>) corresponding to lattice 1 in the 'state' file.
!!   Similar applies for <i>vectors_in_2</i> with lattice 2.
!!   </td>
!!  </tr>
!!  <tr>
!!   <td> -extract_pos_xyz </td>
!!   <td>
!!   Extract the positions of the particles, and output them to stdout in '.xyz' format. In the output
!!   the first line contains the number of particles, the second line is a comment line, and the subsequent lines contain the
!!   particle positions: the first token is the `element' ('A' for species 1, 'B' for species 2, ..., 'Z' for species 26, '?'
!!   otherwise), and the second, third, fourth and fifth tokens are the x-, y- and z-coordinates respectively.
!!   </td>
!!  </tr>
!!  <tr>
!!   <td> -set_wf <i>wf_file<i></td>
!!   <td>
!!   Alters the weight function in the 'state' file to correspond to that specified in the file <i>wf_file</i>. The format of the file
!!   must be analogous to the format of the weight function output by this program via the '-extract_wf' argument: the file
!!   must contain <code>M_grid_size</code> lines, each containing two tokens (extra lines and tokens are ignored). 
!!   Both tokens should be of type <code>real</code>.
!!   The first token on line <i>i</i> is ignored, and the second is the new value of the weight function for macrostate <i>i</i>.
!!   </td>
!!  </tr>
!! </table>
!! </p>
!! <p>
!! This program terminates with a non-zero exit status of 1 for errors. Note that exit statuses are not part of the 
!! Fortran standard, and may not work for all operating systems or compilers.
!! </p>
!!
program monteswitch_post

    !! <h2> Dependencies </h2>
    !! <p> 
    !! <ul>
    !!  <li> <code> kinds_mod </code> </li>
    !!  <li> <code> monteswitch_mod </code> </li>
    !! </ul>
    !! </p>
    use kinds_mod
    use monteswitch_mod

    implicit none


    character(len=20) :: char


    ! THE PROGRAM

    ! Read the command line argument 
    call getarg(1,char)

    if(char=="") then

        write(0,*) "monteswitch_post: Error: No command line argument detected."
        stop 1

    else if(trim(char)=="-extract_wf") then

        call extract_wf()

    else if(trim(char)=="-extract_M_counts") then

        call extract_M_counts()

    else if(trim(char)=="-extract_pos") then

        call extract_pos()

    else if(trim(char)=="-extract_R_1") then

        call extract_R_1()

    else if(trim(char)=="-extract_R_2") then

        call extract_R_2()

    else if(trim(char)=="-extract_u") then

        call extract_u()

    else if(trim(char)=="-calc_rad_dist") then

        call calc_rad_dist()

    else if(trim(char)=="-merge_trans") then

        call merge_trans()

    else if(trim(char)=="-extract_lattices_in") then

        call extract_lattices_in()

    else if(trim(char)=="-extract_pos_xyz") then

        call extract_pos_xyz()

    else if(trim(char)=="-set_wf") then

        call set_wf()

    else

        write(0,*) "monteswitch_post: Error: Unrecognised command line argument."
        stop 1

    end if




contains




    ! CODE FOR EXTRACTING THE WEIGHT FUNCTION FROM THE STATE FILE
    subroutine extract_wf()

        integer(ik) :: i

        call import("state")

        do i=1,M_grid_size
            write(*,*) M_grid(i), eta_grid(i)
        end do

    end subroutine extract_wf




    ! CODE FOR EXTRACTING THE ORDER PARAMETER HISTOGRAM FROM THE STATE FILE
    subroutine extract_M_counts()

        integer(ik) :: i

        call import("state")

        do i=1,M_grid_size
            write(*,*) M_grid(i), M_counts_1(i),M_counts_2(i)
        end do

    end subroutine extract_M_counts




    ! CODE FOR EXTRACTING THE PARTICLE POSITIONS
    subroutine extract_pos()

        real(rk), dimension(:,:), allocatable :: pos
        integer(ik), dimension(:), allocatable :: spec
        
        integer(ik) :: i, error, species

        call import("state")

        allocate(pos(n_part,3))    
        allocate(spec(n_part))

        select case(lattice)
        case(1)
            pos = pos_1
            spec = spec_1
        case(2)
            pos = pos_2
            spec = spec_2
        end select

        call getarg(2,char)
        if(char=="") then

            ! If there is no command line argument then print the positions of all particles
            do i=1,n_part
                write(*,*) pos(i,:)
            end do

        else

            ! Otherwise print only the positions of the particles belonging to the specified (integer) species
            read(char,*,iostat=error) species
            if(error/=0) then
                write(0,*) "monteswitch_post: Error. Problem reading integer after the command line argument '-extract_pos'"
                stop 1
            end if

            do i=1,n_part
                if(spec(i)==species) write(*,*) pos(i,:)
            end do

        end if

    end subroutine extract_pos




    ! CODE FOR EXTRACTING THE LATTICE VECTORS FOR LATTICE TYPE 1
    subroutine extract_R_1()

        integer(ik) :: i, species, error

        call import("state")

        call getarg(2,char)
        if(char=="") then

            ! If there is no command line argument then print the lattice vectors for all particles
            do i=1,n_part
                write(*,*) R_1(i,:)
            end do

        else

            ! Otherwise print only the lattice vectors of the particles belonging to the specified (integer) species
            read(char,*,iostat=error) species
            if(error/=0) then
                write(0,*) "monteswitch_post: Error. Problem reading integer after the command line argument '-extract_R_1'"
                stop 1
            end if

            do i=1,n_part
                if(spec_1(i)==species) write(*,*) R_1(i,:)
            end do

        end if

    end subroutine extract_R_1




    ! CODE FOR EXTRACTING THE LATTICE VECTORS FOR LATTICE TYPE 2
    subroutine extract_R_2()

        integer(ik) :: i, species, error

        call import("state")

        call getarg(2,char)
        if(char=="") then

            ! If there is no command line argument then print the lattice vectors for all particles
            do i=1,n_part
                write(*,*) R_2(i,:)
            end do

        else

            ! Otherwise print only the lattice vectors of the particles belonging to the specified (integer) species
            read(char,*,iostat=error) species
            if(error/=0) then
                write(0,*) "monteswitch_post: Error. Problem reading integer after the command line argument '-extract_R_2'"
                stop 1
            end if

            do i=1,n_part
                if(spec_2(i)==species) write(*,*) R_2(i,:)
            end do

        end if

    end subroutine extract_R_2




    ! CODE FOR EXTRACTING THE PARTICLE DISPLACEMENTS
    subroutine extract_u()

        integer(ik) :: i, species, l, error

        call import("state")

        call getarg(2,char)
        if(char=="") then

            ! If there is no command line argument then print the displacements for all particles
            do i=1,n_part
                write(*,*) u(i,:)
            end do

        else

            ! Otherwise print only the displacements of the particles belonging to the specified (integer) species
            ! in the specified (integer) lattice - which is the next argument
            read(char,*,iostat=error) species
            if(error/=0) then
                write(0,*) "monteswitch_post: Error. Problem reading integer after the command line argument '-extract_u'"
                stop 1
            end if

            call getarg(3,char)
            if(char=="") then
                write(0,*) "monteswitch_post: Error. No 3rd command line argument detected for '-extract_u'"
                stop 1
            end if
            read(char,*,iostat=error) l
            if(error/=0) then
                write(0,*) "monteswitch_post: Error. Problem reading 2nd integer after the command line argument '-extract_u'"
                stop 1
            end if
            if(l/=1 .and. l/=2) then
                write(0,*) "monteswitch_post: Error. Lattice/phase specified for '-extract_u' is not 1 or 2"
                stop 1
            end if
            
            do i=1,n_part
                select case(l)
                case(1)
                    if(spec_1(i)==species) write(*,*) u(i,:)
                case(2)
                    if(spec_2(i)==species) write(*,*) u(i,:)
                end select
            end do

        end if


    end subroutine extract_u




    ! CODE FOR CALCULATING THE RADIAL DISTRIBUTION FUNCTION
    subroutine calc_rad_dist()
    
        real(rk), dimension(:), allocatable :: rad_dist_seps
        integer(ik), dimension(:), allocatable :: rad_dist_counts
        real(rk), dimension(:), allocatable :: rad_dist_function
        character(len=20) :: char
        real(rk) :: sep
        integer(ik) :: bins, i,j

        call getarg(2,char)
        if(char=="") then
            write(0,*) "monteswitch_post: Error: No 2nd command line argument detected (the 'bins' file name for -calc_rad_dist)."
            stop 1
        end if

        read(char,*) bins

        if(bins<=0) then
            write(0,*) "monteswitch_post: Error: 'bins' is <=0."
            stop 1
        end if

        call import("state")

        ! Initialise histogram
        allocate(rad_dist_seps(bins))
        allocate(rad_dist_counts(bins))
        allocate(rad_dist_function(bins))
        do i=1,bins
            select case(lattice)
            case(1)
                rad_dist_seps(i)=(i-1)*min(Lx(1)/2.0_rk,Ly(1)/2.0_rk,Lz(1)/2.0_rk)/bins
            case(2)
                rad_dist_seps(i)=(i-1)*min(Lx(2)/2.0_rk,Ly(2)/2.0_rk,Lz(2)/2.0_rk)/bins
            end select
        end do
        rad_dist_counts=0

        ! Construct histogram
        do i=1,n_part
            do j=1,n_part
                select case(lattice)
                case(1)
                    sep=min_image_distance_fancy(R_1(i,:)+u(i,:),R_1(j,:)+u(j,:),Lx(1),Ly(1),Lz(1))
                case(2)
                    sep=min_image_distance_fancy(R_2(i,:)+u(i,:),R_2(j,:)+u(j,:),Lx(2),Ly(2),Lz(2))
                end select
                call update_histogram(bins,rad_dist_seps,rad_dist_counts,rad_dist_seps(2),sep)
            end do
        end do

        ! Output histogram
        rad_dist_function=(rad_dist_counts*1.0_rk)/n_part
        do i=1,bins
            write(*,*) rad_dist_seps(i),rad_dist_function(i)
        end do

    end subroutine calc_rad_dist




    ! CODE FOR MERGING TRANS MATRICES
    subroutine merge_trans()

        ! 'state' file names to merge and output
        character(len=20) :: state_in_1, state_in_2, state_out
        ! Matrix used for merging trans
        real(rk), dimension(:,:), allocatable :: trans_toadd

        character(len=20) :: char

        call getarg(2,state_in_1)
        if(state_in_1=="") then
            write(0,*) "monteswitch_post: Error: No 2nd command line argument detected for -merge_trans."
            stop 1
        end if
        call getarg(3,state_in_2)
        if(state_in_2=="") then
            write(0,*) "monteswitch_post: Error: No 3rd command line argument detected for -merge_trans."
            stop 1
        end if
        call getarg(4,state_out)
        if(state_out=="") then
            write(0,*) "monteswitch_post: Error: No 4th command line argument detected for -merge_trans."
            stop 1
        end if

        call import(trim(state_in_2))    

        allocate(trans_toadd(M_grid_size,M_grid_size))
        trans_toadd=trans

        call import(trim(state_in_1))
        trans=trans+trans_toadd

        call export(trim(state_out))

    end subroutine merge_trans




    ! CODE FOR EXPORTING IN 'lattices_in' FORMAT
    subroutine extract_lattices_in()

        character(len=20) :: vector_in_1, vector_in_2
        real(rk), dimension(:,:), allocatable :: pos_1_temp
        real(rk), dimension(:,:), allocatable :: pos_2_temp
        integer(ik) :: i

        call import("state")

        call getarg(2,vector_in_1)
        call getarg(3,vector_in_2)

        allocate(pos_1_temp(n_part,3))
        allocate(pos_2_temp(n_part,3))

        ! Set pos_1_temp
        if(trim(vector_in_1)=="pos") then
            pos_1_temp = pos_1
        else if(trim(vector_in_1)=="R") then
            pos_1_temp = R_1
        else if(vector_in_1=="") then
            write(0,*) "monteswitch_post: Error: No 2nd command line argument detected for -extract_lattices_in."
            stop 1
        else
            write(0,*) "monteswitch_post: Error: Unrecognised 2nd command line argument to -extract_lattices_in."
            stop 1
        end if

        ! Set pos_2_temp
        if(trim(vector_in_2)=="pos") then
            pos_2_temp = pos_2
        else if(trim(vector_in_2)=="R") then
            pos_2_temp = R_2
        else if(vector_in_2=="") then
            write(0,*) "monteswitch_post: Error: No 2nd command line argument detected for -extract_lattices_in."
            stop 1
        else
            write(0,*) "monteswitch_post: Error: Unrecognised 2nd command line argument for -extract_lattices_in."
            stop 1
        end if

        ! Output comment line and the number of particles
        write(*,*) "lattices_in file created by monteswitch_post"
        write(*,*) n_part
        ! Output lattice 1 dimensions and particle positions (in fractional coordinates)
        write(*,*) Lx(1)
        write(*,*) Ly(1)
        write(*,*) Lz(1)
        do i=1,n_part
            write(*,*) pos_1_temp(i,1)/Lx(1), pos_1_temp(i,2)/Ly(1), pos_1_temp(i,3)/Lz(1), spec_1(i)
        end do
        ! Output lattice 2 dimensions and particle positions (in fractional coordinates)
        write(*,*) Lx(2)
        write(*,*) Ly(2)
        write(*,*) Lz(2)
        do i=1,n_part
            write(*,*) pos_2_temp(i,1)/Lx(2), pos_2_temp(i,2)/Ly(2), pos_2_temp(i,3)/Lz(2), spec_2(i)
        end do

    end subroutine extract_lattices_in




    ! CODE FOR EXTRACTING THE PARTICLE POSITIONS IN XYZ FORMAT
    subroutine extract_pos_xyz()

        real(rk), dimension(:,:), allocatable :: pos
        integer(ik), dimension(:), allocatable :: spec
        integer(ik) :: i
        character(len=26), parameter :: element="ABCDEFGHIJKLMNOPQRSTUVWXYZ"

        call import("state")

        allocate(pos(n_part,3))
        allocate(spec(n_part))

        select case(lattice)
        case(1)
            pos = pos_1
            spec = spec_1
        case(2)
            pos = pos_2
            spec = spec_2
        end select

        write(*,*) n_part
        write(*,*) "xyz-format file created by monteswitch_post"

        do i=1,n_part

            if(spec(i)>0 .and. spec(i)<=26) then

                write(*,*) element(spec(i):spec(i)), pos(i,:)

            else

                write(*,*) "?", pos(i,:)

            end if

        end do

    end subroutine extract_pos_xyz




    ! CODE FOR SETTING THE WEIGHT FUNCTION IN THE 'STATE' FILE
    subroutine set_wf()

        integer(ik) :: i, error
        real(rk) :: x

        call import("state")

        ! Get the name of the file to read the weight function from
        call getarg(2,char)
        if(char=="") then
            write(0,*) "monteswitch_post: Error. No 2nd command line argument detected for '-set_wf'"
            stop 1
        end if

        ! Open the file
        open(unit=10,file=trim(char),iostat=error,status="old")
        if(error/=0) then
            write(0,*) "monteswitch_post: Error. Problem opening file '",trim(char),"'"
            stop 1
        end if

        ! Read the weight function from the file and set the weight function in 'state' accordingly        
        do i=1,M_grid_size
            read(10,*,iostat=error) x, eta_grid(i)
            if(error/=0) then
                write(0,*) "monteswitch_mod: Error. Problem reading weight function from line ",i," in file '",trim(char),"'"
                stop 1
            end if
        end do

        close(unit=10)

        ! Overwrite the old 'state' file with the new weight function
        call export("state")

    end subroutine set_wf



    !! <h3> <code> function min_image_distance_fancy(r_1,r_2,Lx,Ly,Lz) </code> </h3>
    !! <p>
    !! <code>  min_image_distance_fancy </code>  returns the distance between the 
    !! positions <code>r_1</code> and <code>r_2</code> according to the 
    !! minimum image convention for a periodic cuboid whose faces are x=0, 
    !! x=<code>Lx</code>, y=0, y=<code>Ly</code>, z=0, and z=<code>Lz</code>.
    !! <code>r_1(1)</code> is the x-component of <code>r_1</code>, 
    !! <code>r_1(2)</code> is the y-component, and <code>r_1(3)</code> is the 
    !! z-component; and similarly for <code>r_2</code>.
    !! This function has a wider applicability than the funciton
    !! <code>min_image_distance</code>: it allows <code>r_1</code> and
    !! <code>r_2</code> to be outwith the aforementioned cube.
    !! Note that <code>r_1</code> and <code>r_2</code> must be such that
    !! <code>-Lx</code><=x<<code>2*Lx</code>, <code>-Ly</code><=y<<code>2*Ly</code>
    !! and <code>-Lz</code><=z<<code>2*Lz</code>.
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
    function min_image_distance_fancy(r_1,r_2,Lx,Ly,Lz)
        real(rk), dimension(3), intent(in) :: r_1, r_2
        real(rk), intent(in) :: Lx, Ly, Lz
        real(rk) :: min_image_distance_fancy
        real(rk) :: xsep, ysep, zsep
        ! Calculate the x-sep
        xsep=abs(r_2(1)-r_1(1))
        xsep=xsep-Lx*floor(xsep/Lx)
        xsep=xsep-Lx*floor(2.0_rk*xsep/Lx)
        ! Calculate the y-sep
        ysep=abs(r_2(2)-r_1(2))
        ysep=ysep-Ly*floor(ysep/Ly)
        ysep=ysep-Ly*floor(2.0_rk*ysep/Ly)
        ! Calculate the z-sep
        zsep=abs(r_2(3)-r_1(3))
        zsep=zsep-Lz*floor(zsep/Lz)
        zsep=zsep-Lz*floor(2.0_rk*zsep/Lz)
        ! Calculate the distance
        min_image_distance_fancy=sqrt(xsep*xsep+ysep*ysep+zsep*zsep)
    end function min_image_distance_fancy


end program monteswitch_post
!!
!! </body>
!! </html>
