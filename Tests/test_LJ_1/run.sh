
../../monteswitch -new

echo
echo "Potential energy per particle (with tail corrections) should be -0.03102(6)"
echo "From the simulation (value, then uncertainty):"
grep ' equil_H_1' state | awk '{print $2/480-0.0009304166}'
grep ' sigma_equil_H_1' state | awk '{print $2/480}'
echo


