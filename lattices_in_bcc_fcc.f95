!
! By extracting lines beginning with the regular expression '\!\! ?'
! (ignoring leading whitespace), and then removing matches to the
! regular expression, html documentation corresponding to this 
! source code will be created.
!
!! <html>
!! <head>
!!  <title> lattices_in_bcc_fcc Documentation </title>
!! </head>
!! <body>
!! 
!! <h1> <code>lattices_in_bcc_fcc</code> Documentation </h1>
!!
!! <h2> Author </h2>
!! <p> Tom Underwood </p>
!!
!! <h2> Description </h2>
!! <p>
!! This program outputs a 'lattices_in' file for use by <code>monteswitch</code>
!! and <code>monteswitch_mpi</code> corresponding to lattice 1 being bcc and lattice
!! 2 being fcc (at the same density). The command line arguments are as follows:
!! </p>
!!
!! <table border="1">
!!  <tr>
!!   <td> <b> Argument </b> </td>
!!   <td> <b> Type </b> </td>
!!   <td> <b> Description </b> </td>
!!  </tr>
!!  <tr>
!!   <td> <code> rho </code> </td>
!!   <td> <code> real(rk) </code> </td>
!!   <td> Density (particles per unit volume) of both the bcc and the fcc phases.</td>
!!  </tr>
!!  <tr>
!!   <td> <code> nx </code> </td>
!!   <td> <code> integer(ik), intent(in) </code> </td>
!!   <td> Number of unit cells to create in x direction. </td>
!!  </tr>
!!  <tr>
!!   <td> <code> ny </code> </td>
!!   <td> <code> integer(ik), intent(in) </code> </td>
!!   <td> Number of unit cells to create in y direction. </td>
!!  </tr>
!!  <tr>
!!   <td> <code> nz </code> </td>
!!   <td> <code> integer(ik), intent(in) </code> </td>
!!   <td> Number of unit cells to create in z direction. </td>
!!  </tr>
!! </table>
!!
!! <p>
!! The unit cell here is the conventional 2-atom bct unit cell; for the bcc(fcc) lattice
!! the relative dimensions of the bct unit cell in each Cartesian direction correspond to the
!! bct representation of the bcc(fcc) lattice.
!! </p>
program lattices_in_bcc_fcc

  !! <h2> Dependencies </h2>
  !! <p> 
  !! <ul>
  !!  <li> <code> kinds_mod </code> </li>
  !! </ul>
  !! </p>
  use kinds_mod

  implicit none

  ! Command line arguments
  real(rk) :: rho
  integer(ik) :: nx,ny,nz

  ! fcc and bcc conventional lattice parameters corresponding to rho
  real(rk) :: dfcc, dbcc
  ! Dimensions of bct unit cell used to construct a bcc or fcc supercell
  real(rk) :: a, b, c

  ! Particle positions and supercell dimensions for the current phase (bcc or fcc) under consideration
  real(rk), dimension(:,:), allocatable ::  R
  real(rk) :: Lx
  real(rk) :: Ly
  real(rk) :: Lz

  ! Unimportant variables
  integer(ik) :: i
  character(len=20) :: char


  ! Read the command line arguments
  call getarg(1,char)
  if(char=="") then
     write(0,*) "lattices_in_bcc_fcc: Error. No 1st (rho) command line argument detected."
     stop 1
  end if
  read(char,*) rho
  call getarg(2,char)
  if(char=="") then
     write(0,*) "lattices_in_bcc_fcc: Error. No 2nd (nx) command line argument detected."
     stop 1
  end if
  read(char,*) nx
  call getarg(3,char)
  if(char=="") then
     write(0,*) "lattices_in_bcc_fcc: Error. No 3rd (ny) command line argument detected."
     stop 1
  end if
  read(char,*) ny
  call getarg(4,char)
  if(char=="") then
     write(0,*) "lattices_in_bcc_fcc: Error. No 4th (nz) command line argument detected."
     stop 1
  end if
  read(char,*) nz
  

  ! Output comment line and the number of particles
  write(*,*) "bcc-fcc, rho = ",rho,", nx,ny,nz = ",nx,ny,nz
  write(*,*) 2*nx*ny*nz  


  ! Construct and output the bcc supercell
  dbcc=2.0_rk**(1.0_rk/3.0_rk)*rho**(-1.0_rk/3.0_rk)
  a=dbcc
  b=dbcc
  c=dbcc
  call bct_lattice(a,b,c,nx,ny,nz,R,Lx,Ly,Lz)
  write(*,*) Lx
  write(*,*) Ly
  write(*,*) Lz
  do i=1,2*nx*ny*nz
     write(*,*)  R(i,1)/Lx, R(i,2)/Ly, R(i,3)/Lz, 1
  end do

  ! Construct and output the fcc supercell
  dfcc=4.0_rk**(1.0_rk/3.0_rk)*rho**(-1.0_rk/3.0_rk)
  a=0.5*sqrt(2.0_rk)*dfcc
  b=0.5*sqrt(2.0_rk)*dfcc
  c=dfcc
  call bct_lattice(a,b,c,nx,ny,nz,R,Lx,Ly,Lz)
  write(*,*) Lx
  write(*,*) Ly
  write(*,*) Lz
  do i=1,2*nx*ny*nz
     write(*,*)  R(i,1)/Lx, R(i,2)/Ly, R(i,3)/Lz, 1
  end do




contains




  ! This subroutine initialises the array R, Lx, Ly and Lz - the particle positions and supercell
  ! dimensionts - to correspond to a bct lattice comprised of nx, ny and nz unit cells in the x-, 
  ! y- and z-directions respectively, where a unit cell is a conventional 2-atom bct unit cell with 
  ! dimensions a, b anc c respectively in the x-, y- and z-directions.
  subroutine bct_lattice(a,b,c,nx,ny,nz,R,Lx,Ly,Lz)
    real(rk), intent(in) :: a,b,c
    integer(ik), intent(in) :: nx,ny,nz
    real(rk), dimension(:,:), allocatable, intent(out) :: R
    real(rk), intent(out) :: Lx,Ly,Lz

    ! The lattice vectors (which are used to form the supercell from a single unit
    ! cell)
    real(rk), dimension(3) :: lattice_x
    real(rk), dimension(3) :: lattice_y
    real(rk), dimension(3) :: lattice_z
    ! The position vectors of particles 1, and 2 within a unit cell
    real(rk), dimension(3) :: b_1
    real(rk), dimension(3) :: b_2
    ! Other variables used by this procedure
    integer(ik) :: n, ix, iy, iz

    ! Initialise the required variables
    lattice_x = (/ a, 0.0_rk, 0.0_rk /)
    lattice_y = (/ 0.0_rk, b, 0.0_rk /)
    lattice_z = (/ 0.0_rk, 0.0_rk, c /)
    b_1 = (/ 0.0_rk, 0.0_rk, 0.0_rk /)
    b_2 = (/ 0.5_rk*a, 0.5_rk*b, 0.5_rk*c /)

    ! Set 'Lx', 'Ly' and 'Lz'
    Lx=lattice_x(1)*nx
    Ly=lattice_y(2)*ny
    Lz=lattice_z(3)*nz

    ! Set 'R'. First allocate 'R' to correspond to the fact that the supercell contains
    ! 2*nx*ny*nz particles
    if(allocated(R)) then
       deallocate(R)
    end if
    allocate(R(2*nx*ny*nz,3))
    n=0
    do iz=0,nz-1
       do iy=0,ny-1
          do ix=0,nx-1
             n=n+1
             R(n,:) = ix*lattice_x + iy*lattice_y + iz*lattice_z +  b_1
             n=n+1
             R(n,:) = ix*lattice_x + iy*lattice_y + iz*lattice_z +  b_2 
          end do
       end do
    end do
  end subroutine bct_lattice




end program lattices_in_bcc_fcc
!!
!! </body>
!! </html>
