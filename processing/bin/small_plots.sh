#!/bin/bash


for file in $*
do
	if [[ $file != *_sm.png ]] ; then
	   #echo $file is big
	   newfile=`dirname $file`/`basename $file .png`_sm.png
	   convert $file -geometry  800x1200 $newfile
	#else
	#	echo $file is small
	fi
done
