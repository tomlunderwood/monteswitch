echo "Running simulation 1 of 10"
../../monteswitch -new
for i in 2 3 4 5 6 7 8 9 10; do
  echo "Running simulation $i of 10"
  ../../monteswitch -resume
done
