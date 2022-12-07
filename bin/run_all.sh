#!/usr/bin/env bash

source ~/.bashrc

# set to run on linux network at GT
pushd /home/anewman/GPS/Kivu/KivuGNSS
   # get into environment
   #conda activate Maps
   # pull data from UNAVCO
   conda run -n Maps ./bin/get_UNAVCO_Slns.py
   # get rates and plot timeseries results
   conda run -n Maps ./bin/get_rates.py
   #plot vector map
   conda run -n Maps ./bin/Kivu_GNSS_Map.py
   # overlay images on map   
   ./bin/overlay.sh
   # send it all to github
   ./bin/gitpush.sh
popd
