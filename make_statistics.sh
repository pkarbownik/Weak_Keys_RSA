#!/bin/bash
DIRECTORY="keys"
#for (( cpus=1; cpus<=8; cpus++ ))
#do 
	#for (( gpus=1; gpus<=8; gpus++ ))
	#do 
		for (( keys=10; keys<=10000; keys*=10 ))
		do 
			for (( key_size=1024; key_size<=1024; key_size+=2 ))
			do 
				for (( threads=8; threads<=512; threads*=2 ))
				do 
					echo "**************************************************************"
					echo "GCD_RSA: $cpus\n gpus: $gpus\n keys: $keys\n key_size: $key_size\n threads: $threads\n algorithm: euclid" 
					./GCD_RSA $keys $key_size $threads $DIRECTORY euclid
					echo "GCD_RSA: $cpus\n gpus: $gpus\n keys: $keys\n key_size: $key_size\n threads: $threads\n algorithm: binary" 
					./GCD_RSA $keys $key_size $threads $DIRECTORY binary
					echo "GCD_RSA: $cpus\n gpus: $gpus\n keys: $keys\n key_size: $key_size\n threads: $threads\n algorithm: fast" 
					./GCD_RSA $keys $key_size $threads $DIRECTORY fast
					#srun -p plgrid-gpu -t 00:30:00 -N 1 --cpus-per-task=$cpus --gres=gpu:$gpus  -A plgpkarbownik2017a ./GCD_RSA  $keys $key_size $threads $DIRECTORY$key_size
				done
			done
		done
	#done
#done
echo "**************DONE**************"

