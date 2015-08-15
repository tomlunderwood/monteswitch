mpiexec -n 4 ../../monteswitch_mpi -new
cp state state_AFTER_WF_GEN
sed -i 's/update_eta=[ \t]*T/update_eta= F/' state
sed -i 's/update_trans=[ \t]*T/update_trans= F/' state
sed -i 's/enable_barriers=[ \t]*T/enable_barriers= F/' state
sed -i 's/stop_sweeps=[ \t]*500000/stop_sweeps= 10000000/' state

