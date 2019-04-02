#!/bin/bash

cmds="commands.txt"

i=1
n=`wc -l $cmds | awk '{print $1}'`

while (( $i <= $n ))
do
    nrun=`jobs|wc -l`
    if (( $nrun < 96 ))
    then
	echo '--- Running line' $i of $n
	if (( "$i" % 100 == 0 ))
	then
	    n=`wc -l $cmds | awk '{print $1}'`
	fi
	l=`head -n $i $cmds | tail -n 1`
	echo --- $l
	$l &
	(( i = $i + 1 ))
    else
	sleep 2
    fi
done

echo '--- Finished'
