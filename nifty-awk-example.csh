#!/bin/csh -f

foreach n (`seq -w 0100 119`)

   echo starting $n

   setenv f /work/clas12/mcdata/generated/lund/sidis/clasdispr.00.e11.000.emn0.75tmn.09.xs65.61nb.dis.$n.dat
   
   foreach N (`seq 0 19`)
      setenv min `echo "$N * 2000" | bc` 
      setenv max `echo "($N + 1) * 2000" | bc` 

      echo "" > temp.txt # start with empty file
      awk -v min=$min -v max=$max '{ if(NF == 10) i++; if(i > min && i <= max) print $0 }' $f >> temp.txt
      tail -n +2 temp.txt > /volatile/clas12/first_experiment/sidis/production/pass4/lund/clasdispr.00.e11.000.emn0.75tmn.09.xs65.61nb.dis."$n"_$N.dat # remove first blank line and give real filename
      rm temp.txt

      echo "   ... $N" 
   end 
   
end
