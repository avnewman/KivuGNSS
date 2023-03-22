#!/usr/bin/env python
"""
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
from GNSSrates import ITRF14EulerPole,SariaEulerPole,getProcessedGNSS
from GNSSrates import linearRates,timeSeriesPlots,UKF
from GNSSrates import commonMode,timeSeriesPlotsCommon

### 3. get reference plate
iPlate='NUBI' # Nubia-ITRF14
NUB_ITRF = ITRF14EulerPole(iPlate)
# local plate of reference
lPlate='Victoria'   # read relative to Nubian plate
VIC_NUB = SariaEulerPole(lPlate)
# get local relative to it's refernce plate
VIC_ITRF = VIC_NUB + NUB_ITRF  
#print(NUB_ITRF,VIC_NUB,VIC_ITRF)

### 4. pull data from data center
getProcessedGNSS(Stats,analysisCenter=analysisCenter)

### 5. perform per-station processing
fitDict = linearRates(Stats,ratesfile,lPlate,VIC_ITRF,analysisCenter,mmpm)
# new common-mode processing
dfCommon,fitCommon=commonMode(Stats,commonratesfile,analysisCenter,meanList)

### 6. create timeseries plots
timeSeriesPlots(Stats,fitDict,analysisCenter,mmpm)
timeSeriesPlotsCommon(dfCommon,fitCommon,analysisCenter,meanList,yscale=mmpm)

### 7. Kalman Filtered PLots
UKF(Stats,fitDict,analysisCenter,mmpm,UKFerr)