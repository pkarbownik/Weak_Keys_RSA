#!/bin/bash
DIRECTORY="keys"
#for (( cpus=1; cpus<=8; cpus++ ))
#do 
	#for (( gpus=1; gpus<=8; gpus++ ))
	#do 
		for (( keys=2000; keys<=2000; keys+=100 ))
		do 
			for (( key_size=1024; key_size<=1024; key_size+=2 ))
			do 
				for (( threads=8; threads<=512; threads*=2 ))
				do 
					echo "***********************************************"
					echo "cpus: $cpus\n gpus: $gpus\n keys: $keys\n key_size: $key_size\n blocks: $blocks\n threads: $threads\n" 
					#srun -p plgrid-gpu -t 00:30:00 -N 1 --cpus-per-task=$cpus --gres=gpu:$gpus  -A plgpkarbownik2017a ./GCD_RSA  $keys $key_size $threads $DIRECTORY$key_size
					./GCD_RSA $keys $key_size $threads $DIRECTORY

				done
			done
		done
	#done
#done
echo "**************DONE**************"

