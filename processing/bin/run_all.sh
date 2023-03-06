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
   # create small copies of files
   processing/bin/small_plots.sh plots/*png
   # organize plots
   mv plots/*TS.png plots/TS/
   mv plots/*TS_sm.png plots/TS/small
   mv plots/*UKF.png plots/UKF/
   mv plots/*UKF_sm.png plots/UKF/small


   # send it all to github
   ./processing/bin/gitpush.sh
   # cp files to webdir
   cp -a ./plots/* ~/html/research/KivuGNSS/plots 
popd
