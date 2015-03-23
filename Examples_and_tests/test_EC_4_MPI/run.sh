echo "Running simulation to generate weight function"
mpiexec -n 2 ../../monteswitch_mpi -new
# Alter the state file to prepare for the production run
sed -i 's/stop_sweeps=      1000000/stop_sweeps=    10000000/' state
sed -i 's/update_eta=  T/update_eta=  F/' state
sed -i 's/update_trans=  T/update_trans=  F/' state
sed -i 's/calc_equil_properties=  F/calc_equil_properties=  T/' state
echo "Running production simulation"
mpiexec -n 2 ../../monteswitch_mpi -reset
# Extract the weight function and the order parameter histogram from the final 'state' file
../../monteswitch_post -extract_wf > wf.dat
../../monteswitch_post -extract_M_counts | awk '{print $1,$2+$3}' > M_hist.dat
