#!/bin/bash

COUNT=1

for (( i=1; i<=${COUNT}; i++ ))
do
	nohup ./Main.sage > ./res${i}.txt 2>&1 &
done
