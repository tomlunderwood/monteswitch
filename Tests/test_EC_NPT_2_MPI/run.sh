echo "Running simulation 1 of 5"
mpiexec -n 2 ../../monteswitch_mpi -new 
for i in 2 3 4 5; do
  echo "Running simulation $i of 5"
  mpiexec -n 2 ../../monteswitch_mpi -resume 
done

