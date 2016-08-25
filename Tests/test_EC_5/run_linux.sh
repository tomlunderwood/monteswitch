echo "Running simulation to generate weight function"
../../monteswitch -new
../../monteswitch_post -extract_wf > wf.dat
awk '$1~"barrier_macro_low:" {print $2,$3}' data > macro_vs_t.dat
awk '$1~"M:" {print $2,$3}' data > M_vs_t.dat
