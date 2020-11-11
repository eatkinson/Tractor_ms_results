#!/bin/bash

#wd=$(pwd)
#qsub simulate_power.sh -lwalltime=3:00:00 -d $wd -e errandout/ -o errandout/ -F "-n 7"

while getopts n: option
do
  case "${option}"
    in
      n) nodeuse=${OPTARG};;
    esac
done

 module load R
 
 #Write the start and stop points of the file
 jstart=1 #changed jstart to 11 to finish last job
 jstop=$nodeuse

 for j in $(seq -w $jstart 1 $jstop) 
 do
  Rscript /mnt/sdb/genetics/pgc_ptsd_grant_resub/lanc_simulation_v5_mafdif.txt $j &
 done

wait

