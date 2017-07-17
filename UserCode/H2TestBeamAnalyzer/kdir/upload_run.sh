#!/bin/bash

if [ -z $1 ]; then
	echo "No run number specified"
	exit 1;
fi
filename=analysis_run_$1.root

if [ ! -f $filename ]; then
	echo "$filename does not exist"
	exit 2;
fi

hostname=`hostname`

if [ $hostname == cmskam05.cern.ch ]; then
	target_location="/data/import"
else
	target_location="/net/cmskam05/data/import"
fi

#exit 0;
cp $filename $target_location
sleep 1
firefox http://cmskam05.cern.ch/db_runs/jt_new_run_form.php?nrun=$1 &

