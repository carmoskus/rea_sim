#!/bin/bash

i=1

while true
do
    nrun=`jobs|wc -l`
    if (( $nrun < 16 ))
    then
	echo 'Running line' $i
    else
	sleep 5
    fi
done
