echo "Running simulation 1 of 10 to generate weight function"
../../monteswitch -new
for i in 2 3 4 5 6 7 8 9 10; do
  echo "Running simulation $i of 10 to generate weight function"
  ../../monteswitch -resume
done
../../monteswitch_post -extract_wf > wf.dat
awk '$1~"barrier_macro_low:" {print $2,$3}' data > macro_vs_t.dat
awk '$1~"M:" {print $2,$3}' data > M_vs_t.dat
