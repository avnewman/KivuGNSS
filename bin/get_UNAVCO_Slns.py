#!/usr/bin/env python

# ### 1. Importing necessary libraries
import requests
import os

#creating url variable for each station
# CWU processing
ancenter='cwu'
#refFrame='nam14'
refFrame='igs14'
# UNR processing
#ancenter='unr'
#refFrame='igs08'
processing='Uncleaned'
refCoord='from_analysis_center'
import math

Stats=('BNZA','BYAH', 'KANZ', 'NYBA', 'RUBO', 'IWAW', 'KMBR')

nStats=len(Stats)
count=0
datadir='timeseries'
if not os.path.exists(datadir):
    os.mkdir(datadir)

for stat in Stats:
    URL='https://web-services.unavco.org/gps/data/position/' + stat +    '/v3?analysisCenter=' + ancenter +    '&referenceFrame=' + refFrame +    '&report=long&dataPostProcessing=' + processing +    '&refCoordOptio='+ refCoord
    
    #creating pandas readable csv from the url
    file=os.path.join(datadir,stat+'.csv')
    count += 1
    print('Requesting data for: '+stat+'  ['+str(count)+'/'+str(nStats)+']') # usable output
    req = requests.get(URL) # get data
    url_content = req.content
    csv_file = open(file, 'wb')
    csv_file.write(url_content)
    csv_file.close()
    
    
