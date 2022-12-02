#!/usr/bin/env python
# needs to run within the 'Maps' conda environment on Andy's Mac

# Rwanda-Kivu Rift region
# we will start by making a prelim plot to assure that this is the region we want

import pygmt
import elevation
import os
import pandas as pd
import math

plotdir='./plots'
if not os.path.exists(plotdir):
    os.mkdir(plotdir)
#xmin=28.8; xmax=30; ymin=-2.7;ymax=-1.2
xmin=29; xmax=29.5; ymin=-2.1;ymax=-1.5
region=[xmin,xmax,ymin,ymax]
proj="M0/0/15c"

libdir=os.path.join(os.getcwd(),'mapdata')
mapcpt=os.path.join(libdir,'map_gray.cpt')
DEM=os.path.join(libdir,'Rwanda_DEM.tif') 
VELS=os.path.join(os.getcwd(),'TS_rates.txt')
Plate=os.path.join(os.getcwd(),'SO-NNR_itrf2014.txt') # Somalian Plate wrt NNR in ITRF2014 using UNAVCO PMC

def round_up(n, rnd):
    '''
    round up number, n by rounding value, rnd.
    thus n=16.8, rnd=5 outputs 20
    '''
    return math.ceil(n / rnd) * rnd

def round_down(n, rnd):
    '''
    round down number, n by rounding value, rnd.
    thus n=16.8, rnd=5 outputs 15
    '''
    return math.floor(n / rnd) * rnd



vels=pd.read_csv(VELS, delim_whitespace=True)
plate=pd.read_csv(Plate,delim_whitespace=True)
# remove plate motion
vels["NPvel"]=vels["Nvel"]-plate["Nvel"]
vels["EPvel"]=vels["Evel"]-plate["Evel"]

vels["ENcor"]=0-vels["NEcor"]
data=vels[["Long","Lat","EPvel","NPvel","Eerr","Nerr","ENcor","Uvel","STAT"]] # EP NP are plate removed
# add legend velocity
vslat = float(-2.03)
vslong = float(29.35)
velscale = pd.Series([vslong,vslat, 10, 0, 5, 2, 0.4, 10, " "],index=data.columns)
vdf=velscale.to_frame().T
# Need to reset floats to all be floats again...pita
vdf[["Long","Lat","EPvel","NPvel","Eerr","Nerr","ENcor","Uvel"]] = vdf[["Long","Lat","EPvel","NPvel","Eerr","Nerr","ENcor","Uvel"]].apply(pd.to_numeric)
data = pd.concat([data,vdf],ignore_index=True)

### GET DEM ###
# gets proper path declaration independent of OS
# this will default to the global 90m resolution --way big enough for us
## should define product='SRTM3' = 90m data (3arcsec), otherwise will default to 1arcsec
#elevation.clean()
#elevation.clip(bounds=(xmin,ymin,xmax,ymax),output=DEM,product='SRTM3')
#elevation.clean()  # allows us to rerun get command
#print(pygmt.grdinfo(DEM)) # shows useful information about the grid file

fig1=pygmt.Figure()
pygmt.config(MAP_FRAME_TYPE="plain", # no alternating B&W frame
             FORMAT_GEO_MAP='ddd.xx') # decimal degrees

#pygmt.makecpt(series=[1000,4000,500],  #create a topo cpt just for this range (500-4500m)
#              continuous=True,
#              cmap='topo')

#fig1.grdimage(DEM,region=region,projection=proj,cmap=mapcpt,shading=True,dpi=300, transparency=60) # Use globe version (all high elev.)
DEMgrad = pygmt.grdgradient(grid=DEM, azimuth=[0,90], normalize='e.8')
pygmt.makecpt(cmap="gray",series=[-1,0.5,0.01])
fig1.grdimage(DEMgrad,region=region,projection=proj,cmap=True,dpi=600, transparency=60) # Use globe version (all high elev.)


# remap the national borders as dashed lines
fig1.coast(region=region,  # xmin,xmax,ymin,ymax
    projection=proj,
    frame=['p','WSen','xa0.1', 'ya.1'], 
    resolution='f', 
    borders='1/1.2p,150,-.-',
    transparency=10,
    )

# Legend box -- hardwired :( 
fig1.plot(x=29.25, y=-2.025, 
          style="R14/4/0.25", 
          color="255/255/235", 
          pen="2p,black",
         transparency=10,
         )
cinc=5  # CPT increment
pygmt.makecpt(cmap="turbo",  #reverse=True, 
              continuous=False,
              series=[round_down(vels.Uvel.min(),cinc),round_up(vels.Uvel.max(),cinc),cinc]
              #series=[0,60,5]
             )

fig1.plot(x=data.Long,y=data.Lat,
          style='d0.5c', 
          cmap=True,
          pen='0.5,0', 
          color=data.Uvel,
          transparency=0)


fig1.velo(data=data,
         region=region,
         pen="0.5p",
         #zvalue='u', # user-defined before STAT name column in data
         #cmap=True,
         spec="e0.2/0.65/10",
         vector="0.5c+p1p+e"
         )

fig1.colorbar(
    #position="jBC", 
    position="JBC+o0c/-1.1c+w11c/0.2c",
    #box="+gwhite+p2,black", # white fill
    #frame=['x+lvertical [mm/yr]']
    )

fig1.text(
    text=["Vertical velocity [mm/yr]", "Horizontal velocity 10 mm/yr (1 @~s@~)"],
    x=[29.25,29.25],
    y=[-2.05,vslat]
    )
fig1.text(text="Kivu GNSS Sites",x=29.05,y=-1.985, justify="LB", font="12p,Helvetica-Bold,black")
fig1.text(text="Observations are relative to stable Somalian Plate",x=29.05,y=-2.0, justify="LB", font="10p,Helvetica,black")
fig1.text(text=vels.Sdate[0] +" to "+vels.Edate[0],x=29.05,y=-2.015, justify="LB", font="10p,Helvetica,black")

fig1.basemap(region=region,  # xmin,xmax,ymin,ymax
    projection=proj,
    map_scale   = '29.40/-1.972/-1.9/10',
    )

fig1.savefig(os.path.join(plotdir,'TS_rates_SONNR.png'),  # types include png,jpg,pdf,bmp,tif,eps,kml
            transparent=False, # transp background for png only
            crop=True, # removes whitespace around fig
            anti_alias=True, # creates smoother plots
            show=True, # display externally too
            dpi=300 #this is default for png
            )
