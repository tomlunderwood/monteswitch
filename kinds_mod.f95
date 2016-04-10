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
