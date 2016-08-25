echo "Running test_EC_8_MPI..."
mpiexec -n 4 ../../monteswitch_mpi -new 
echo
echo "* equil_Delta_F should be -0.0415888:"
grep ' equil_DeltaF' state
grep ' sigma_equil_DeltaF' state
echo
echo "* equil_umsd_1 should be 0.015:"
grep ' equil_umsd_1' state
grep ' sigma_equil_umsd_1' state
echo
echo "* equil_umsd_2 should be 0.0075:"
grep ' equil_umsd_2' state
grep ' sigma_equil_umsd_2' state
echo
echo "* equil_H_1 should be 0.06:"
grep ' equil_H_1' state
grep ' sigma_equil_H_1' state
echo
echo "* equil_H_2 should be 0.06:"
grep ' equil_H_2' state
grep ' sigma_equil_H_2' state
echo
rm -f state_? data_?
