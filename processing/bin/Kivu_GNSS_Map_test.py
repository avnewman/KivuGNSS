#!/usr/bin/env python
# needs to run within the 'Maps' conda environment on Andy's Mac

# Rwanda-Kivu Rift region
# we will start by making a prelim plot to assure that this is the region we want

import pygmt
import elevation
import os
import pandas as pd
import datetime
import math

plotdir='./plots'
if not os.path.exists(plotdir):
    os.mkdir(plotdir)
#xmin=28.8; xmax=30; ymin=-2.7;ymax=-1.2
xmin=29; xmax=29.5; ymin=-2.1;ymax=-1.5
region=[xmin,xmax,ymin,ymax]
proj="M0/0/15c"

libdir=os.path.join(os.getcwd(),'processing','mapdata')
mapcpt=os.path.join(libdir,'map_gray.cpt')
DEM=os.path.join(libdir,'Rwanda_DEM.tif') 
BATH=os.path.join(libdir,'kivu_mergedBIG.grd') 
ratesFile=os.path.join(os.getcwd(),'rates','Kivu_rates.txt')

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

rates=pd.read_csv(ratesFile, delim_whitespace=True)

rates["ENcor"]=0-rates["NEcor"]
data=rates[["Long","Lat","Evel-SO","Nvel-SO","Eerr","Nerr","ENcor","Uvel","STAT"]] # EP NP are plate removed
# add legend velocity
vslat = float(-2.19)
vslong = float(29.38)
velscale = pd.Series([vslong,vslat, 10, 0, 5, 2, 0.4, 10, " "],index=data.columns)
vdf=velscale.to_frame().T
# Need to reset floats to all be floats again...pita
vdf[["Long","Lat","Evel-SO","Nvel-SO","Eerr","Nerr","ENcor","Uvel"]] = vdf[["Long","Lat","Evel-SO","Nvel-SO","Eerr","Nerr","ENcor","Uvel"]].apply(pd.to_numeric)
data = pd.concat([data,vdf],ignore_index=True)

### GET DEM ###
# gets proper path declaration independent of OS
# this will default to the global 90m resolution --way big enough for us
## should define product='SRTM3' = 90m data (3arcsec), otherwise will default to 1arcsec
#   elevation.clean()
#   #elevation.clip(bounds=(xmin,ymin,xmax,ymax),output=DEM,product='SRTM3')
#   elevation.clip(bounds=(xmin,ymin,xmax,ymax),output=DEM)
#   elevation.clean()  # allows us to rerun get command
#   print(pygmt.grdinfo(DEM)) # shows useful information about the grid file

fig1=pygmt.Figure()
pygmt.config(MAP_FRAME_TYPE="plain", # no alternating B&W frame
    FORMAT_GEO_MAP='ddd.xx') # decimal degrees

#fig1.grdimage(DEM,region=region,projection=proj,cmap=mapcpt,shading=True,dpi=300, transparency=60) # Use globe version (all high elev.)
#DEMgrad = pygmt.grdgradient(grid=DEM, azimuth=[0,90], normalize='e.8')

# convert kml fault map into something usable (only needed 1x)
kml=os.path.join(libdir,'doc.kml')  # obtained by unzipping Cindy's kivufaults.kmz
faults=os.path.join(libdir,'faults.txt')  # obtained by unzipping Cindy's kivufaults.kmz
if not os.path.exists(faults):  
    cmd = "gmt kml2gmt "+kml+" > "+faults
    returned_value = os.system(cmd)  # returns the exit code in unix
    print('KML2GMT returned value:', returned_value)

# start building plot
#fig1.grdimage(DEMgrad,region=region,projection=proj,
hcoast=1458
pygmt.makecpt(series=[1220,hcoast,50],  #create a topo cpt just for this range (500-4500m)
   continuous=True,
   #reverse=True,
   cmap='wysiwyg')
BATHclip = pygmt.grdclip(BATH, above=[hcoast, 'NaN'])
BATHgrad = pygmt.grdgradient(grid=BATH, azimuth=[90], normalize='e100')
fig1.grdimage(BATH,region=region,projection=proj,
    shading=BATHgrad,
    cmap=True,
    dpi=300, 
    transparency=30)
fig1.coast(region=region,  # xmin,xmax,ymin,ymax
    projection=proj,
    resolution='f', 
    lakes=True,
    land='white',
    )

# Set all grid points < value to NaN
DEMclip = pygmt.grdclip(DEM, below=[hcoast, 'NaN'])
DEMgrad = pygmt.grdgradient(grid=DEMclip, azimuth=[0,90], normalize='e0.8')
pygmt.makecpt(cmap="gray",series=[-1,0.5,0.01])
fig1.grdimage(DEMgrad,region=region,projection=proj,
    shading=True,
    cmap=True,
    dpi=300, 
    transparency=60)


# Legend box -- hardwired :( 
fig1.plot(x=29.25, y=-2.188, 
    style="R15/4/0.25", 
    color="255/255/235", 
    pen="2p,black",
    transparency=50,
    no_clip=True,
    )

# add lake 
#fig1.coast(region=region,  # xmin,xmax,ymin,ymax
#    projection=proj,
#    resolution='f', 
#    lakes=True,
#    water='200/255/255',
#    transparency=60,
#    )
# remap the national borders as dashed lines
fig1.coast(region=region,  # xmin,xmax,ymin,ymax
    projection=proj,
    resolution='f', 
    borders='1/1.2p,120,-.-',
    transparency=10,
    )
# remap the national borders as dashed lines
fig1.text(text=["DRC","Rwanda"],
    x=[29.302,29.36],
    y=[-1.525,-1.53],
    font="14p,Helvetica-Bold,120",
    justify="LB", 
    transparency=10,
    )
fig1.basemap(
    frame=['WSeN','xa0.1i0.02','ya0.1i0.02' ,'+e'], 
    )
# faults from Cindy
fig1.plot(data=faults,
    pen="0.5p,darkred",
    transparency=10,
    )
cinc=5  # CPT increment
pygmt.makecpt(cmap="turbo",  #reverse=True, 
    continuous=False,
    series=[round_down(rates.Uvel.min(),cinc),round_up(rates.Uvel.max(),cinc),cinc]
    #series=[0,60,5]
    )
fig1.plot(x=data.Long,y=data.Lat,
    no_clip=True,
    style='d0.5c', 
    cmap=True,
    pen='0.5,0', 
    color=data.Uvel,
    transparency=0)
fig1.velo(data=data,
    no_clip=True,
    region=region,
    pen="0.5p",
    #zvalue='u', # user-defined before STAT name column in data
    #cmap=True,
    spec="e0.2/0.65/10",
    vector="0.5c+p1p+e"
    )
fig1.colorbar(
    #position="jBC", 
    position="JBC+o0c/3.6c+w11c/0.2c",
    #box="+gwhite+p2,black", # white fill
    #frame=['x+lvertical [mm/yr]']
    )
fig1.text(
    text=["Vertical velocity [mm/yr]", "Horizontal velocity 10 mm/yr (1 @~s@~)"],
    no_clip=True,
    x=[29.25,29.28],
    y=[-2.21,vslat]
    )
fig1.text(text="Kivu GNSS Sites",
    x=29.02,y=-2.145, 
    no_clip=True,
    justify="LB", 
    font="12p,Helvetica-Bold,black"
    )
now=datetime.datetime.utcnow().strftime("%d-%b-%Y %H:%M UTC")
print(now)
fig1.text(text=now,
    x=29.02,y=-2.157, 
    no_clip=True,
    justify="LB", 
    font="7p,Helvetica-Oblique"
    )
fig1.text(text="Observations are relative to stable Somalian Plate",
    x=29.20,y=-2.145, 
    no_clip=True,
    justify="LB", 
    font="10p,Helvetica,black"
    )
fig1.text(text=rates.Sdate[0] +" to "+rates.Edate[0],
    x=29.20,y=-2.16, 
    no_clip=True,
    justify="LB", 
    font="10p,Helvetica,black"
    )

pygmt.config(MAP_TICK_PEN_PRIMARY='2p')  # thicker line for map_scale
fig1.basemap(region=region,  # xmin,xmax,ymin,ymax
    projection=proj,
    map_scale   = '29.43/-2.070/-1.9/10',
    )

xmin=29; xmax=29.5; ymin=0;ymax=1
region=[xmin,xmax,ymin,ymax]
proj="X15c/4c"

# polygons and lines are always clipped :( need to offset plot
fig1.plot(x=[29.01, 29.04, 29.06], y=[0.41, 0.47, 0.47],
    region=region,  
    projection=proj,
    yshift=-4.5,
    pen="0.5p,darkred",
    transparency=10,
    )
fig1.text( text="Mapped Faults", x=29.08, y=0.45, justify="LB" )
fig1.text( text="(Wood et al., 2015, Smets et al., 2016)", x=29.01, y=0.34, justify="LB",font="8p,Helvetica-Oblique,black" )

fig1.savefig(os.path.join(plotdir,'KIVU_GNSS_test.png'),  # types include png,jpg,pdf,bmp,tif,eps,kml
    transparent=False, # transp background for png only
    crop=True, # removes whitespace around fig
    anti_alias=True, # creates smoother plots
    show=False, # display externally too
    dpi=300 #this is default for png
    )
