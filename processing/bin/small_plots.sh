#!/bin/bash

for file in $*
do
	newfile=`dirname $file`/`basename $file .png`_sm.png
	convert $file -geometry 500x500  $newfile
done
