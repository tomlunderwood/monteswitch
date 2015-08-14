
rm -f *.dat
for i in  $(seq 1 100); do
  rho=`echo |awk -v i=$i 'END{print (1+2.0*i/100)/(2*2*2)}'` 
  ../../lattices_in_hcp_fcc $rho 7 4 2 > lattices_in
  echo "Running calculation $i of 100"
  ../../monteswitch -new
  awk -v rho=$rho '$1~"E_1=" {E1=$2} $1~"E_2=" {E2=$2} END{print rho*(2*2*2),(E2-E1)/(0.5*672)}' state >> Ediff_vs_rho.dat
done

