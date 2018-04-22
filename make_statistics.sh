#!/bin/bash
DIRECTORY="keys"
#for (( cpus=1; cpus<=8; cpus++ ))
#do 
	#for (( gpus=1; gpus<=8; gpus++ ))
	#do 
		for (( keys=100; keys<=100; keys+=100 ))
		do 
			for (( key_size=512; key_size<=512; key_size++ ))
			do 
				for (( threads=512; threads<=512; threads+=8 ))
				do 
					echo "***********************************************"
					echo "cpus: $cpus\n gpus: $gpus\n keys: $keys\n key_size: $key_size\n blocks: $blocks\n threads: $threads\n" 
					#srun -p plgrid-gpu -t 00:30:00 -N 1 --cpus-per-task=$cpus --gres=gpu:$gpus  -A plgpkarbownik2017a ./GCD_RSA  $keys $key_size $threads $DIRECTORY$key_size
					./GCD_RSA $keys $key_size $threads $DIRECTORY$key_size

				done
			done
		done
	#done
#done
echo "**************DONE**************"

