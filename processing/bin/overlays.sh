#!/usr/bin/env bash

# Create overlay plots atop map plot
thumbdir=./processing/Logos
plotdir=./plots

files="$plotdir/KIVU_GNSS.png  $plotdir/KIVU_GNSS_common.png"

# create thumbs on andy's laptop 
KIVU=/Users/anewman/Documents/Grants/Kivu2022-2024/Logo/Kivu-Logo.png 
NSF=/Users/anewman/Documents/figures/LOGOS/NSF_4-Color_vector_Logo.png
GT=/Users/anewman/Documents/figures/LOGOS/GT-LOGO/2022/GTVertical_Navy.png
TUL=/Users/anewman/Documents/figures/LOGOS/Tulane_logo_240x240.png

KIVUth=$thumbdir/Kivu_th.png
NSFth=$thumbdir/NSF_th.png
GTth=$thumbdir/GT_th.png
TULth=$thumbdir/TUL_th.png

# only do once
 #   inaction="-verbose -units "PixelsPerInch" -density 100 "
 #   outaction=" -resize 200x200 "
 #   outaction=" -resize 180x180 "
 #   convert $inaction $KIVU $outaction  $KIVUth
 #   convert $inaction $NSF  $outaction  $NSFth
 #   convert $inaction $GT   $outaction  $GTth
 #   convert $inaction $TUL   $outaction  $TULth


for file in $files
do 
  OUT=$plotdir/`basename $file .png`_ov.png
  convert $file \
   $KIVUth -geometry +160+70 -composite \
   $NSFth -geometry +160+1960 -composite \
   $GTth -geometry +360+1960 -composite \
   $TULth -geometry +560+1980 -composite \
   $OUT
  rm $file
  mv $OUT $file
done
