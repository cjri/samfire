for dd in do*.sh; do
  chmod 755 $dd
  ./$dd
  for ddd in D?; do
    cd $ddd
      ls Haps*out > Haps.list
      mv Haps.list ../
    cd ../
  done
done
