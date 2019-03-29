#!/bin/bash

cmds="commands.txt"

i=1
n=`wc -l $cmds | awk '{print $1}'`

while (( $i <= $n ))
do
    nrun=`jobs|wc -l`
    if (( $nrun < 96 ))
    then
	echo 'Running line' $i
	`head -n $i $cmds | tail -n 1` &
	(( i = $i + 1 ))
    else
	sleep 2
    fi
done
