#!/bin/sh --login
#SBATCH -n 20      
#SBATCH -t 30:00:00  
#SBATCH -A aosc-hi 
#SBATCH -J maooam
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

ENS_NUM=21
MIN_NUM=21
TW_DA=1.D-1

while [ "${ENS_NUM}" -ge "${MIN_NUM}" ]
do 

   sed -i "9s/.*/  TW_DA = $TW_DA/" da_params.nml 
   sed -i "15s/.*/  ENS_NUM = $ENS_NUM/" da_params.nml # Don't forget to add "-i" to save
  ./etkf_maooam
  python compute_rmse.py ${ENS_NUM}
  ENS_NUM=`expr $ENS_NUM - 1`

done

#mv evol_field.dat evol_field_1_1.dat

#python gradual_coupling.py

exit 0
