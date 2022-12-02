#!/usr/bin/env bash

# get into environment
#conda activate Maps
# pull data from UNAVCO
./bin/get_UNAVCO_Slns.py
# get rates and plot timeseries results
./bin/get_rates.py
#plot vector map
conda run -n Maps ./bin/Kivu_GNSS_Map.py
