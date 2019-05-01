#!/bin/bash

cmds="commands.txt"

i=1
n=`wc -l $cmds | awk '{print $1}'`

echo '--- Reading from' $cmds

while (( $i <= $n ))
do
    (( nto = 96 - `jobs|wc -l` ))
    if (( $nto > 0 ))
    then
#	if (( nto == 1 ))
#	then
	    echo '--- Running line' $i of $n '-' $nto
	    l=`sed -n ${i}p "$cmds"`
	    echo --- $l
	    $l &
	    (( i = $i + 1 ))
#	fi
	
	if (( "$i" % 1024 == 0 ))
	then
	    n=`wc -l $cmds | awk '{print $1}'`
	fi
    else
	sleep 1
    fi
done

wait

echo '--- Finished'
