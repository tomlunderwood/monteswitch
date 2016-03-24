!
! By extracting lines beginning with the regular expression '\!\! ?'
! (ignoring leading whitespace), and then removing matches to the
! regular expression, html documentation corresponding to this 
! source code will be created.
!
!! <html>
!! <head>
!!  <title> lattices_init_hcp_fcc Documentation </title>
!! </head>
!! <body>
!! 
!! <h1> <code>lattices_init_hcp_fcc</code> Documentation </h1>
!!
!! <h2> Author </h2>
!! <p> Tom Underwood </p>
!!
!! <h2> Description </h2>
!! <p>
!! This program outputs a 'lattices_in' file for use by <code>monteswitch</code>
!! and <code>monteswitch_mpi</code> corresponding to lattice 1 being hcp and lattice
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
!!   <td> Density (particles per unit volume) of both the hcp and the fcc phases.</td>
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
!! The unit cell here contains 12 particles, 2 in each z-plane, with a stacking 
!! sequence of ABCABC for fcc and ABABAB for hcp.
!! </p>
program lattices_in_hcp_fcc

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

  ! The nearest neighbour distance
  real(rk) :: nn_dist
  ! The lattice vectors (which are used to form the supercell from a single unit
  ! cell)
  real(rk), dimension(3) :: lattice_x
  real(rk), dimension(3) :: lattice_y
  real(rk), dimension(3) :: lattice_z
  ! The position vectors of particles 1 and 2 in the z=0 plane
  real(rk), dimension(3) :: b_1
  real(rk), dimension(3) :: b_2
  real(rk), dimension(3) :: b_2_AB
  real(rk), dimension(3) :: b_2_C
  ! The z-offset between adjacent planes
  real(rk), dimension(3) :: z_offset
  ! The x/y offset for a B plane relative to an A plane
  real(rk), dimension(3) :: B_offset
  ! The x/y offset for a C plane relative to an A plane
  real(rk), dimension(3) :: C_offset
  ! The position of the particle under consideration
  real(rk), dimension(3) :: R
  ! The supercell dimensions
  real(rk) :: Lx
  real(rk) :: Ly
  real(rk) :: Lz
  ! Other variables used by this procedure
  integer(ik) :: n, ix, iy, iz, plane
  character(len=20) :: char
  real(rk), dimension(3) :: offset
  

  ! Read the command line arguments
  call getarg(1,char)
  if(char=="") then
     write(0,*) "lattices_in_hcp_fcc: Error. No 1st (rho) command line argument detected."
     stop 1
  end if
  read(char,*) rho
  call getarg(2,char)
  if(char=="") then
     write(0,*) "lattices_in_hcp_fcc: Error. No 2nd (nx) command line argument detected."
     stop 1
  end if
  read(char,*) nx
  call getarg(3,char)
  if(char=="") then
     write(0,*) "lattices_in_hcp_fcc: Error. No 3rd (ny) command line argument detected."
     stop 1
  end if
  read(char,*) ny
  call getarg(4,char)
  if(char=="") then
     write(0,*) "lattices_in_hcp_fcc: Error. No 4th (nz) command line argument detected."
     stop 1
  end if
  read(char,*) nz


  ! Initialise the required variables
  nn_dist = 2**(1.0_rk/6.0_rk)/rho**(1.0_rk/3.0_rk)
  lattice_x = (/ nn_dist, 0.0_rk, 0.0_rk /)
  lattice_y = (/ 0.0_rk, nn_dist*sqrt(3.0_rk), 0.0_rk /)
  lattice_z = (/ 0.0_rk, 0.0_rk, 6.0_rk*nn_dist*sqrt(6.0_rk)/3.0_rk /)
  b_1 = (/ 0.0_rk, 0.0_rk, 0.0_rk /)
  b_2_AB = (/ nn_dist/2.0_rk, nn_dist*sqrt(3.0_rk)/2.0_rk , 0.0_rk /)
  b_2_C = (/ -nn_dist/2.0_rk, nn_dist*sqrt(3.0_rk)/2.0_rk , 0.0_rk /)
  B_offset = (/ 0.0_rk, nn_dist/sqrt(3.0_rk), 0.0_rk /)
  C_offset = (/ nn_dist/2.0_rk, nn_dist*sqrt(3.0_rk)/6.0_rk, 0.0_rk /)
  z_offset = (/ 0.0_rk, 0.0_rk, nn_dist*sqrt(6.0_rk)/3.0_rk /)
  
  ! Set 'Lx', 'Ly' and 'Lz'
  Lx=lattice_x(1)*nx
  Ly=lattice_y(2)*ny
  Lz=lattice_z(3)*nz

  ! Output comment line and the number of particles
  write(*,*) "hcp-fcc, rho = ",rho,", nx,ny,nz = ",nx,ny,nz
  write(*,*) 12*nx*ny*nz

  ! Output the hcp particle positions (in fractional coordinates)
  write(*,*) Lx
  write(*,*) Ly
  write(*,*) Lz
  n=0
  do ix=0,nx-1
     do iy=0,ny-1
        do iz=0,nz-1
           do plane=0,5
              if(mod(plane,2)==0) then
                 offset = plane*z_offset
                 b_2=b_2_AB
              else if(mod(plane,2)==1) then
                 offset = plane*z_offset + B_offset
                 b_2=b_2_AB
              end if
              n=n+1
              R = ix*lattice_x + iy*lattice_y + iz*lattice_z +  b_1 + offset
              write(*,*) R(1)/Lx, R(2)/Ly, R(3)/Lz, 1
              n=n+1
              R = ix*lattice_x + iy*lattice_y + iz*lattice_z +  b_2 + offset
              write(*,*) R(1)/Lx, R(2)/Ly, R(3)/Lz, 1
           end do
        end do
     end do
  end do

  ! Output the fcc lattice positions
  write(*,*) Lx
  write(*,*) Ly
  write(*,*) Lz
  n=0
  do ix=0,nx-1
     do iy=0,ny-1
        do iz=0,nz-1
           do plane=0,5
              if(mod(plane,3)==0) then
                 offset = plane*z_offset
                 b_2=b_2_AB
              else if(mod(plane,3)==1) then
                 offset = plane*z_offset + B_offset
                 b_2=b_2_AB
              else if(mod(plane,3)==2) then
                 offset = plane*z_offset + C_offset
                 b_2=b_2_C
              end if
              n=n+1
              R = ix*lattice_x + iy*lattice_y + iz*lattice_z +  b_1 + offset
              write(*,*) R(1)/Lx, R(2)/Ly, R(3)/Lz, 1
              n=n+1
              R = ix*lattice_x + iy*lattice_y + iz*lattice_z +  b_2 + offset
              write(*,*) R(1)/Lx, R(2)/Ly, R(3)/Lz, 1
           end do
        end do
     end do
  end do




end program lattices_in_hcp_fcc
!!
!! </body>
!! </html>
