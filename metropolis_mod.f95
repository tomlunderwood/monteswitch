!
! By extracting lines beginning with the regular expression '\!\! ?'
! (ignoring leading whitespace), and then removing matches to the
! regular expression, html documentation corresponding to this 
! source code will be created.
!
!! <html>
!! <head>
!!  <title> metropolis_mod Documentation </title>
!! </head>
!! <body>
!! 
!! <h1> <code>metropolis_mod</code> Documentation </h1>
!!
!! <h2> Author </h2>
!! <p> Tom Underwood </p>
!!
!! <h2> Description </h2>
!! <p>
!! <code>metropolis_mod</code> contains functions for evaluating 
!! probabilities of accepting moves according to the Metropolis
!! algorithm.
!! </p>
!!
module metropolis_mod

  !! <h2> Dependencies </h2>
  !! <p>
  !! <code> kinds_mod </code>
  !! </p>
  use kinds_mod

  implicit none

contains

  !! <h2> Procedures </h2>

  !! <h3> <code> function metropolis_prob_B(beta,E,E_trial) </code> </h3>
  !! <p>
  !! <code> metropolis_prob_B </code> gives the probability of a move being 
  !! accepted via the Metropolis algorithm, given the current and trial 
  !! energies of the system and Boltzmann sampling. This function should
  !! not be used for volume moves.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> beta </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Thermodynamic beta. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> E </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Current energy of the system. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> E_trial </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Trial energy of the system. </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
  function metropolis_prob_B(beta,E,E_trial)
    real(rk), intent(in) :: beta
    real(rk), intent(in) :: E
    real(rk), intent(in) :: E_trial
    real(rk) :: metropolis_prob_B
    metropolis_prob_B=min(exp(-beta*(E_trial-E)),1.0_rk)
  end function metropolis_prob_B

  

  !! <h3> <code> function metropolis_prob_B_vol(beta,E,E_trial,N,P,V,V_trial) </code> </h3>
  !! <p>
  !! <code> metropolis_prob_B_vol </code> gives the probability of a move being 
  !! accepted via the Metropolis algorithm, given the current and trial 
  !! energies and volumes of the system and Boltzmann sampling, and volume moves
  !! corresponding to a random walk in 'ln(V)-space'.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> beta </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Thermodynamic beta. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> E </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Current energy of the system. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> E_trial </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Trial energy of the system. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> N </code> </td>
  !!   <td> <code> integer(rk), intent(in) </code> </td>
  !!   <td> Number of particles in the system. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> P </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Pressure of the system. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> V </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Current volume of the system. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> V_trial </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Trial volume of the system. </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
  function metropolis_prob_B_vol(beta,E,E_trial,N,P,V,V_trial)
    real(rk), intent(in) :: beta
    real(rk), intent(in) :: E
    real(rk), intent(in) :: E_trial
    integer(ik), intent(in) :: N
    real(rk), intent(in) :: P
    real(rk), intent(in) :: V
    real(rk), intent(in) :: V_trial
    real(rk) :: metropolis_prob_B_vol
    metropolis_prob_B_vol=min(exp(-beta*((E_trial-E)+P*(V_trial-V))+(N+1)*log(V_trial/V)),1.0_rk)
  end function metropolis_prob_B_vol



  !! <h3> <code> function metropolis_prob_MC(beta,E,E_trial,eta,eta_trial) </code> </h3>
  !! <p>
  !! <code> metropolis_prob_MC </code> gives the probability of a move being 
  !! accepted via the Metropolis algorithm, given the current and trial 
  !! energies of the system and weight functions and multicanonical sampling. 
  !! This function should not be used for volume moves. Note that the sign
  !! convention of the weight function is such that configurations of a 
  !! system with higher weight functions will be sampled more often. 
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> beta </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Thermodynamic beta. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> E </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Current energy of the system. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> E_trial </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Trial energy of the system. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> eta </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Current value of the weight function. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> eta_trial </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Trial value of the weight function. </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
  function metropolis_prob_MC(beta,E,E_trial,eta,eta_trial)
    real(rk), intent(in) :: beta
    real(rk), intent(in) :: E
    real(rk), intent(in) :: E_trial
    real(rk), intent(in) :: eta
    real(rk), intent(in) :: eta_trial
    real(rk) :: metropolis_prob_MC
    metropolis_prob_MC=min(exp(-beta*(E_trial-E)+eta_trial-eta),1.0_rk)
  end function metropolis_prob_MC

  

  !! <h3> <code> function metropolis_prob_MC_vol(beta,E,E_trial,eta,eta_trial,N,P,V,V_trial) </code> </h3>
  !! <p>
  !! <code> metropolis_prob_MC_vol </code> gives the probability of a move being 
  !! accepted via the Metropolis algorithm, given the current and trial 
  !! energies, volumes and values of the weight function of the system and 
  !! multicanonical sampling, and volume moves corresponding to a random walk in 
  !! 'ln(V)-space'.
  !! </p>
  !! <table border="1">
  !!  <tr>
  !!   <td> <b> Argument </b> </td>
  !!   <td> <b> Type </b> </td>
  !!   <td> <b> Description </b> </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> beta </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Thermodynamic beta. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> E </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Current energy of the system. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> E_trial </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Trial energy of the system. </td>
  !!  </tr>
  !!   <td> <code> eta </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Current value of the weight function. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> eta_trial </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Trial value of the weight function. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> N </code> </td>
  !!   <td> <code> integer(rk), intent(in) </code> </td>
  !!   <td> Number of particles in the system. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> P </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Pressure of the system. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> V </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Current volume of the system. </td>
  !!  </tr>
  !!  <tr>
  !!   <td> <code> V_trial </code> </td>
  !!   <td> <code> real(rk), intent(in) </code> </td>
  !!   <td> Trial volume of the system. </td>
  !!  </tr>
  !! </table>
  !! <p><b>Returns:</b> <code> real(rk) </code> </p>
  function metropolis_prob_MC_vol(beta,E,E_trial,eta,eta_trial,N,P,V,V_trial)
    real(rk), intent(in) :: beta
    real(rk), intent(in) :: E
    real(rk), intent(in) :: E_trial  
    real(rk), intent(in) :: eta
    real(rk), intent(in) :: eta_trial
    integer(ik), intent(in) :: N
    real(rk), intent(in) :: P
    real(rk), intent(in) :: V
    real(rk), intent(in) :: V_trial
    real(rk) :: metropolis_prob_MC_vol
    metropolis_prob_MC_vol=min(exp(-beta*((E_trial-E)+P*(V_trial-V))+eta_trial-eta+(N+1)*log(V_trial/V)),1.0_rk)
  end function metropolis_prob_MC_vol


end module metropolis_mod
!!
!! </body>
!! </html>
