!
! By extracting lines beginning with the regular expression '\!\! ?'
! (ignoring leading whitespace), and then removing matches to the
! regular expression, html documentation corresponding to this 
! source code will be created.
!
!! <html>
!! <head>
!!  <title> kinds_mod Documentation </title>
!! </head>
!! <body>
!! 
!! <h1> <code>kinds_mod</code> Documentation </h1>
!!
!! <h2> Author </h2>
!! <p> Tom Underwood </p>
!!
!! <h2> Description </h2>
!! <p>
!! <code>kinds_mod</code> is a module which contains definitions for the 'kinds' of 
!! integers and reals used by other modules and programs.
!! </p>
!!
module kinds_mod
  
  implicit none
  
  !! <h2> Variables </h2>
  !! <table border="1">
  !! <tr>
  !!  <td> <b> Variable </b> </td>
  !!  <td> <b> Type </b> </td>
  !!  <td> <b> Description </b> </td>
  !!  <td> <b> Default value </b> </td>
  !! </tr>

  !! <tr>
  !!  <td> <code> ik </code> </td>
  !!  <td> <code> integer, parameter </code> </td>
  !!  <td> Defines a kind for integer variables </td>
  !!  <td> <code> selected_int_kind(9) </code> </td>
  !! </tr>
  integer, parameter :: ik=selected_int_kind(9)
  !! <tr>
  !!  <td> <code> rk </code> </td>
  !!  <td> <code> integer, parameter </code> </td>
  !!  <td> Defines a kind for real variables </td>
  !!  <td> <code> selected_real_kind(15,307) </code> </td>
  !! </tr>
  integer, parameter :: rk=selected_real_kind(15,307)

  !! </table>

end module kinds_mod
!!
!! </body>
!! </html>
