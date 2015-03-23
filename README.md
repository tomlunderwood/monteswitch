# monteswitch
monteswitch is a program which uses the Monte Carlo method to evaluate the free energy difference between two 
user-defined phases with user-defined inter-particle interactions. It is special because it does this very efficiently, 
using, in addition to the conventional particle and volume moves, a 'switch' Monte Carlo move which jumps directly 
between the two phases, bypassing the free energy barrier separating them. The two 'phases' are user-defined, and could 
be two different crystalline solids, or two different Hamiltonians. 
