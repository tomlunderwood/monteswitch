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
!!  <title> lattices_in_bcc_hcp Documentation </title>
!! </head>
!! <body>
!! 
!! <h1> <code>lattices_in_bcc_hcp</code> Documentation </h1>
!!
!! <h2> Author </h2>
!! <p> Tom Underwood </p>
!!
!! <h2> Description </h2>
!! <p>
!! This program outputs a 'lattices_in' file for use by <code>monteswitch</code>
!! and <code>monteswitch_mpi</code> corresponding to lattice 1 being bcc and lattice
!! 2 being hcp (at the same density). The command line arguments are as follows:
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
!!   <td> Density (particles per unit volume) of both the bcc and the hcp phases.</td>
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
!! The unit cell here for both phases resembles a conventional 4-atom fcc unit cell stretched/compressed in all
!! three Cartesian directions, i.e., the unit cells resemble a conventional 4-atom face-centred tetragonal (fct) 
!! unit cell. For the bcc lattice the relative dimensions of the bct unit cell in each Cartesian direction 
!! correspond to the fct representation of the bcc lattice. The unit cell is a cuboid with dimensions 
!! <code>a=c=\sqrt(2)*b</code>, and the atomic positions within a unit cell are 
!! <code>(0,0,0), (a/2,b/2,0), (a/2,0,c/2), (0,b/2,c/2)</code>.
!! The hcp lattice cannot exactly be represented by the fct structure. However, it can if the basis is slightly
!! altered. For the hcp lattice the unit cell is a cuboid with dimensions <code>a=sqrt(3)*b</code> and <code>c=(2*sqrt(6)/3)*b</code>,
!! where <code>b</code> happens to be the nearest neighbour distance. The atomic positions within a unit cell are
!! <code>(0,0,0), (a/2,b/2,0), (a/3,0,c/2), (5*a/6,b/2,c/2)</code>; note the change of the 3rd and 4th atoms relative
!! to the exact fct (which is given above in the discussion for bcc).
!! </p>
program lattices_in_bcc_hcp

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

  ! Dimensions of fct unit cell used to construct a bcc or hcp supercell
  real(rk) :: a, b, c

  ! Supercell dimensions for the current phase (bcc or hcp) under consideration
  real(rk) :: Lx
  real(rk) :: Ly
  real(rk) :: Lz

  ! Origin of unit cell under consideration
  real(rk), dimension(3) :: cell_origin
  ! Unimportant variables
  integer(ik) :: i,j,k
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
  write(*,*) "bcc-hcp, rho = ",rho,", nx,ny,nz = ",nx,ny,nz
  write(*,*) 4*nx*ny*nz  


  ! Construct and output the bcc supercell
  a=2.0_rk**(5.0_rk/6.0_rk)/rho**(1.0_rk/3.0_rk)
  c=a
  b=a/sqrt(2.0_rk)
  Lx=nx*a
  Ly=ny*b
  Lz=nz*c
  write(*,*) Lx
  write(*,*) Ly
  write(*,*) Lz
  do i=1,nx
     do j=1,ny
        do k=1,nz
           cell_origin=(/(i-1)*a,(j-1)*b,(k-1)*c/)
           write(*,*)  cell_origin(1)/Lx,            cell_origin(2)/Ly,            cell_origin(3)/Lz,           1
           write(*,*) (cell_origin(1)+a/2.0_rk)/Lx, (cell_origin(2)+b/2.0_rk)/Ly,  cell_origin(3)/Lz,           1
           write(*,*) (cell_origin(1)+a/2.0_rk)/Lx,  cell_origin(2)/Ly,           (cell_origin(3)+c/2.0_rk)/Lz, 1
           write(*,*)  cell_origin(1)/Lx,           (cell_origin(2)+b/2.0_rk)/Ly, (cell_origin(3)+c/2.0_rk)/Lz, 1
        end do
     end do
  end do

  ! Construct and output the hcp supercell
  b=2.0_rk**(1.0_rk/6.0_rk)/rho**(1.0_rk/3.0_rk)
  a=sqrt(3.0_rk)*b
  c=(2.0_rk*sqrt(6.0_rk)/3.0_rk)*b
  Lx=nx*a
  Ly=ny*b
  Lz=nz*c
  write(*,*) Lx
  write(*,*) Ly
  write(*,*) Lz
  do i=1,nx
     do j=1,ny
        do k=1,nz
           cell_origin=(/(i-1)*a,(j-1)*b,(k-1)*c/)
           write(*,*)  cell_origin(1)/Lx,                   cell_origin(2)/Ly,            cell_origin(3)/Lz,           1
           write(*,*) (cell_origin(1)+a/2.0_rk)/Lx,        (cell_origin(2)+b/2.0_rk)/Ly,  cell_origin(3)/Lz,           1
           write(*,*) (cell_origin(1)+a/3.0_rk)/Lx,         cell_origin(2)/Ly,           (cell_origin(3)+c/2.0_rk)/Lz, 1
           write(*,*) (cell_origin(1)+5.0_rk*a/6.0_rk)/Lx, (cell_origin(2)+b/2.0_rk)/Ly, (cell_origin(3)+c/2.0_rk)/Lz, 1
        end do
     end do
  end do




end program lattices_in_bcc_hcp
!!
!! </body>
!! </html>
