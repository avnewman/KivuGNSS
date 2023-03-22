#!/usr/bin/env python
"""
NOTE:  Items tested in here have been incorporated back into get_rates2.py...this version is only for testing!
Script to pull data, and process through velocity results and time series time series plots
"""
# 1. defaults
analysisCenter = 'cwu'  # or 'cwu'
#analysisCenter = 'unr'  # or 'cwu'
# what stations
#Stats=('BNZA','BYAH', 'KANZ', 'NYBA', 'RUBO', 'IWAW', 'KMBR','UGN3','MBAR','ASUM','MOIU','KYN7','KYN2')
# only MBAR has current data outside of our network
Stats=('BNZA','BYAH', 'KANZ', 'NYBA', 'RUBO', 'IWAW','KMBR') # above not processed at cwu
meanList=('BNZA','BYAH', 'KANZ', 'NYBA', 'KMBR') # get daily average of these stations for net-average
ratesfile='Kivu_rates.txt'
commonratesfile='Kivu_rates_common.txt'
mmpm=1000 # to get to mm displacements
UKFerr=[3,0.001,5] #[std measurement error, std process error, initial uncertainty in mm] 

### 2. Importing necessary libraries
import os
import sys
basedir=os.getcwd() # meant to be run in directory above procdir, if not, some otherthings need be done
# processing and base info
procdir=os.path.join(basedir,'processing')
bindir=os.path.join(procdir,'bin')
sys.path.append(os.path.join(bindir,'GNSSrates'))
from GNSSrates import timeSeriesPlotsCommon,commonMode

### 5. perform per-station processing
dfMaster,fit=commonMode(Stats,commonratesfile,analysisCenter,meanList)
timeSeriesPlotsCommon(dfMaster,fit,analysisCenter,meanList,yscale=mmpm)