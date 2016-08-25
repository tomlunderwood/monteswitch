echo "Running simulation 1 of 2 to generate weight function"
mpiexec -n 2 ../../monteswitch_mpi -new
echo "Running simulation 2 of 2 to generate weight function"
mpiexec -n 2  ../../monteswitch_mpi -resume
../../monteswitch_post -extract_wf > wf.dat
