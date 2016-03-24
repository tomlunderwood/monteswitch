mpiexec -n 4 ../../monteswitch_mpi -new -explicit
cp state state_AFTER_WF_GEN
cp state_0 state_0_AFTER_WF_GEN
cp state_1 state_1_AFTER_WF_GEN
cp state_2 state_2_AFTER_WF_GEN
cp state_3 state_3_AFTER_WF_GEN
cp data_0 data_0_AFTER_WF_GEN
cp data_1 data_1_AFTER_WF_GEN
cp data_2 data_2_AFTER_WF_GEN
cp data_3 data_3_AFTER_WF_GEN
sed -i '' 's/update_eta=[ \t]*T/update_eta= F/' state
sed -i '' 's/update_trans=[ \t]*T/update_trans= F/' state
sed -i '' 's/enable_barriers=[ \t]*T/enable_barriers= F/' state
sed -i '' 's/stop_sweeps=[ \t]*2000000/stop_sweeps= 125000000/' state
sed -i '' 's/equil_sweeps=[ \t]*0/equil_sweeps= 250000/' state
mpiexec -n 4 ../../monteswitch_mpi -reset -explicit

