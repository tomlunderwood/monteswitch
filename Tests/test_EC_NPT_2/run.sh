echo "Running simulation 1 of 5"
../../monteswitch -new
for i in 2 3 4 5; do
  echo "Running simulation $i of 5"
  ../../monteswitch -resume
done
