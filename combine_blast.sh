#!/bin/bash

BASEDIR=/bio/abchase/MG/elevation_gradient

count=0

while [ $count -lt 4 ]
do
	cd $BASEDIR/time${count}/

	ls *.sh | sed 's/.sh//g' > temp.txt

	while read line
	do
	
		cat ${line}*.blast.txt > ${line}.blast.total.txt
		wc -l ${line}.blast.total.txt

	done < temp.txt

	rm temp.txt
	# rm *.e*
	# rm *.o*

	count=`expr $count + 1`
done
