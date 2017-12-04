#!/bin/bash

for (( c=1; c<=$1; c++ ))
do  
   openssl genrsa -out $c.pem 4096
   openssl rsa -pubout -in $c.pem -out $c.pem
   echo "$c keys generated"
done
echo "All $1 keys generated successfully"



