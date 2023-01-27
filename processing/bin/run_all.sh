#!/usr/bin/env bash

source ~/.bashrc

# set to run on linux network at GT
pushd /home/anewman/GPS/Kivu/KivuGNSS
      #conda run -n Maps ./processing/bin/get_UNAVCO_Slns.py
      #conda run -n Maps ./processing/bin/get_rates.py
      #conda run -n Maps ./processing/bin/UKF.py 
   # above processing within a single script now
   conda run -n Maps ./processing/bin/get_rates2.py
   #plot vector map
   conda run -n Maps ./processing/bin/Kivu_GNSS_Map.py
   # overlay images on map   
   ./processing/bin/overlays.sh
   # send it all to github
   ./processing/bin/gitpush.sh
   # create small copies of files
   processing/bin/small_plots.sh plots/*png
   # cp files to webdir
   cp ./plots/*.png ~/html/research/KivuGNSS/ 
popd
