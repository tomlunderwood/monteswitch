echo "Running test_EC_9_MPI..."
echo
echo "Running simulation 1 of 10"
mpiexec -n 4 ../../monteswitch_mpi -new -explicit
for i in 2 3 4 5 6 7 8 9 10; do
  echo "Running simulation $i of 10"
  mpiexec -n 4 ../../monteswitch_mpi -resume -explicit
done
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
rm -f data state
