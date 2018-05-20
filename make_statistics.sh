#!/bin/bash
DIRECTORY="keys"
for (( cpus=1; cpus<=1; cpus++ ))
do 
	for (( gpus=1; gpus<=1; gpus++ ))
	do 
		for (( keys=100; keys<=100; keys*=10 ))
		do 
			for (( key_size=1024; key_size<=1024; key_size+=2 ))
			do 
				for (( threads=512; threads<=512; threads*=2 ))
				do 
					echo "**************************************************************"
					echo "GCD_RSA: $cpus\n gpus: $gpus\n keys: $keys\n key_size: $key_size\n threads: $threads\n algorithm: euclid" 
					srun -p plgrid-gpu -t 01:00:00 -N 1 --cpus-per-task=$cpus --gres=gpu:$gpus  -A plgpkarbownik2017a ./GCD_RSA  $keys $key_size $threads $DIRECTORY$key_size euclid
					./GCD_RSA $keys $key_size $threads $DIRECTORY euclid
					echo "GCD_RSA: $cpus\n gpus: $gpus\n keys: $keys\n key_size: $key_size\n threads: $threads\n algorithm: binary" 
					srun -p plgrid-gpu -t 01:00:00 -N 1 --cpus-per-task=$cpus --gres=gpu:$gpus  -A plgpkarbownik2017a ./GCD_RSA  $keys $key_size $threads $DIRECTORY$key_size binary
					./GCD_RSA $keys $key_size $threads $DIRECTORY binary
					echo "GCD_RSA: $cpus\n gpus: $gpus\n keys: $keys\n key_size: $key_size\n threads: $threads\n algorithm: fast" 
					srun -p plgrid-gpu -t 01:00:00 -N 1 --cpus-per-task=$cpus --gres=gpu:$gpus  -A plgpkarbownik2017a ./GCD_RSA  $keys $key_size $threads $DIRECTORY$key_size fast
					./GCD_RSA $keys $key_size $threads $DIRECTORY fast
				done
			done
		done
	#done
#done
echo "**************DONE**************"

