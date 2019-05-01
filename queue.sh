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
	if (( nto == 1 ))
	then
	    echo '--- Running line' $i of $n
	    l=`sed -n ${i}p "$cmds"`
	    echo --- $l
	    $l &
	    
	    if (( "$i" % 1024 == 0 ))
	    then
		n=`wc -l $cmds | awk '{print $1}'`
		echo Reloading $i
	    fi

	    (( i = $i + 1 ))
	else
	    # Run lines $i through $j
	    (( j = $nto - 1 ))
	    echo '--- Running' $j lines after $i of $n
	    code=`sed -n "${i},+${j}p" "$cmds"`
	    IFS=$'\n'
	    for l in $code
	    do
		IFS=$' \t\n'
		$l &
	    done

	    if (( $i % 1024 > ( $i + $j ) % 1024 ))
	    then
		n=`wc -l $cmds | awk '{print $1}'`
		echo Reloading $i
	    fi

	    (( i = $i + $j + 1 ))
	fi
	
    else
	sleep 1
    fi
done

wait

echo '--- Finished'
