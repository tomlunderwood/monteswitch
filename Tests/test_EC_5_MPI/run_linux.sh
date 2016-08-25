echo "Running simulation to generate weight function"
mpiexec -n 2 ../../monteswitch_mpi -new
../../monteswitch_post -extract_wf > wf.dat
awk '$1~"barrier_macro_low:" {print $2,$3}' data_0 > macro_vs_t_0.dat
awk '$1~"M:" {print $2,$3}' data_0 > M_vs_t_0.dat
awk '$1~"barrier_macro_low:" {print $2,$3}' data_1 > macro_vs_t_1.dat
awk '$1~"M:" {print $2,$3}' data_1 > M_vs_t_1.dat
