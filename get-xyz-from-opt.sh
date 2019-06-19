#!/bin/bash

# Get .xyz files for each optimized geometry separately,
# to be used in subsequent identification of geometric center
# for NICS calculation:

for i in # file names
do
   for j in DMF gas o-DCB
   do
      for k in m062x b3lyp
      do
         dirs=$j'_'$k
         if [ -d "$dirs" ]
         then
            cd $j'_'$k
            if [ -d $i ]
            then
               cd $i
               filename=$i'.out'
               echo $filename $j $k
               if [ -f $filename ]
               then
                  echo "file exists"
                  # -n, the length of string is not zero
                  if [[ -n "$(grep "Normal termination" $filename | tail -1)" ]]
                  then
                     echo "job finished"
                     # print final geometries to .xyz file
                     sed -n 'H; /Standard orientation/h; ${g;p;}' $i.out | sed -n '/Standard orientation/,/Rotational/p' | sed -n '/1/,/----/p' | sed -n '/-------------/q;p' |  cut -c 17-19,32-95 | cat | wc -l | cat >> ../../$i.xyz
                     echo $i $j $k | cat >> ../../$i.xyz
                     sed -n 'H; /Standard orientation/h; ${g;p;}' $i.out | sed -n '/Standard orientation/,/Rotational/p' | sed -n '/1/,/----/p' | sed -n '/-------------/q;p' | cut -c 17-28,32-95 | sed -n 's/^17 / Cl /g;p' | sed -n 's/^ 6 / C /g;p' | sed -n 's/^ 1 / H /g;p' | sed -n 's/^ 7 / N /g;p' | sed -n 's/^ 8 / O /g;p' | sed -n 's/^16 / S /g;p' | cat >> ../../$i.xyz
                  else
                     echo "job did not finish"
                  fi
               else
                  echo "file does not exist"
                  exit
               fi
               cd ..
               cd ..
            else
               echo $i $j $k does not exist
               cd ..
            fi
#            sleep 2
         else
            echo $j $k does not exist
         fi
      done
   done
done
