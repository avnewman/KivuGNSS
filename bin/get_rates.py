#!/usr/bin/env python

# ### 1. Importing necessary libraries
import pandas as pd
from matplotlib import pyplot as plt
import os
import numpy as np
import math

yscale=1000 # to get to mm displacements
xscale=1/(86400e9*365.25) #convert numeric nanosecond into fractional year

Stats=('BNZA','BYAH', 'KANZ', 'NYBA', 'RUBO', 'IWAW', 'KMBR')
#Stats=('BNZA',)
#GetStats = False  # turn off for testing once data has been downloaded

def wcorr(df):
    """
    df contains dN,E,U and wN,E,U for dataset to determine
    the weighed covarience matrix
    """
    Q=np.array(df[['dN','dE','dU']])
    W=np.array(df[['wN','wE','wU']])
    QW=Q*W
    C=QW.T.dot(QW)/W.T.dot(W)  # covariance
    CorrNE=C[0][1]/(C[0][0]*C[1][1])
    CorrNU=C[0][2]/(C[0][0]*C[2][2])
    CorrEU=C[1][2]/(C[1][1]*C[2][2])
    return CorrNE,CorrNU,CorrEU


def xyz2blh(x, y, z):
    """Convert XYZ coordinates to BLH,
    return tuple(latitude, longitude, height).
    
    from https://github.com/purpleskyfall/XYZ2BLH/blob/master/xyz2blh.py
    """
    A = 6378137.0
    B = 6356752.314245
    e = math.sqrt(1 - (B**2)/(A**2))
    # calculate longitude, in radians
    longitude = math.atan2(y, x)

    # calculate latitude, in radians
    xy_hypot = math.hypot(x, y)

    lat0 = 0
    latitude = math.atan(z / xy_hypot)

    while abs(latitude - lat0) > 1E-9:
        lat0 = latitude
        N = A / math.sqrt(1 - e**2 * math.sin(lat0)**2)
        latitude = math.atan((z + e**2 * N * math.sin(lat0)) / xy_hypot)

    # calculate height, in meters
    N = A / math.sqrt(1 - e**2 * math.sin(latitude)**2)
    if abs(latitude) < math.pi / 4:
        R, phi = math.hypot(xy_hypot, z), math.atan(z / xy_hypot)
        height = R * math.cos(phi) / math.cos(latitude) - N
    else:
        height = z / math.sin(latitude) - N * (1 - e**2)

    # convert angle unit to degrees
    longitude = math.degrees(longitude)
    latitude = math.degrees(latitude)

    return latitude, longitude, height

datadir='./timeseries'
plotdir='./plots'
if not os.path.exists(datadir):
    os.mkdir(datadir)
if not os.path.exists(plotdir):
    os.mkdir(plotdir)

ratesFile=open('TS_rates.txt', 'w')
ratesFile.write('STAT        Lat      Long   Height      Nvel     Evel     Uvel      Nerr     Eerr     Uerr    NEcor   NUcor   EUcor         Sdate       Edate\n') 
for stat in Stats:
    # errors
    csize=2
    elw=0.8
    ecol='k'
    # markers 
    mec='k' #edge color
    mew=0.8 #width
    mfmt='o' # shape 
    msz='4' # size
    # final
    falpha=1
    fcol='blue'
    # rapid
    ralpha=0.5
    rcol='red'
    # grid
    galpha=0.5
    gvwidth=0.5
    ghwidth=1
    gcol='gray'
    
    #creating pandas readable csv from the url
    file=os.path.join(datadir,stat+'.csv')
    df = pd.read_csv(file, header=8)
    df['Date'] =  pd.to_datetime(df['Datetime'])
    df['NumDate']=pd.to_numeric(df.Date)*xscale # converts to date in NANOseconds
    
    # get xyz coords from header (reread file)
    csv_file=open(file,'r')
    content = csv_file.readlines()
    csv_file.close()
    xyzline=content[6]
    words=xyzline.split()
    xyz=[float(words[5]),float(words[7]),float(words[9])] # xyz coordinates
    lat,lon,height=xyz2blh(xyz[0],xyz[1],xyz[2])

    # remove mean of data
    Comps=('N','E', 'U')
    for comp in Comps:
        ycolComp=str(' delta '+comp)
        df[ycolComp]=df[ycolComp]-df[ycolComp].mean()

    dfF1=df.loc[df[' Solution']=='final']          #final solutions
    dfF2=df.loc[df[' Solution']=='suppl']          #final solutions
    dfF3=df.loc[df[' Solution']=='suppf']          #final solutions
    dfFinal=pd.concat([dfF1,dfF2,dfF3]).sort_index()
    dfRapid=df.loc[df[' Solution']=='rapid'].reset_index()          #rapid solutions

    #creating Station (3 component) plots 
    f, ax = plt.subplots(3, 1, figsize=(14,8),sharey=False, sharex=True) 
    f.tight_layout(h_pad=0)
    
    x=df['Date']
    xND=df['NumDate']
    xF=dfFinal['Date']
    xFND=dfFinal['NumDate']
    xR=dfRapid['Date']
    xRND=dfRapid['NumDate']
    sp=0
    slopes=np.zeros(3) # store slopes and errors
    errs=np.zeros(3)
    dfFNEU=pd.DataFrame()

    for comp in Comps:
        ycolComp=str(' delta '+comp)
        wcolComp=str(' Std Dev '+comp)
        y=df[ycolComp]*yscale
        yF=dfFinal[ycolComp]*yscale
        yR=dfRapid[ycolComp]*yscale
        yeF=dfFinal[wcolComp]*yscale
        yeR=dfRapid[wcolComp]*yscale
        wghtF=1/dfFinal[wcolComp]/yscale
        
        # performing linear fits to data
        sdate=dfFinal.Date[0].strftime("%Y-%m-%d")
        edate=dfFinal.Date[len(dfFinal)-1].strftime("%Y-%m-%d")
        fit,var=np.polyfit(xFND,yF,1,w=wghtF, full=False, cov=True)
        fity=np.polyval(fit,xFND)
        errs[sp]=np.sqrt(var[0][0])
        slopes[sp]=fit[0]
        dfFNEU['d'+comp]=(yF-fity) # detrended sln for covariance determination
        dfFNEU['w'+comp]=(wghtF)   #  weights 
        
        #plot
        ax[sp].plot(xF, fity)
        ax[sp].errorbar(x=xF, y=yF, yerr=yeF, 
            fmt=mfmt, ms=msz, capsize=csize, label='Final', mfc=fcol, mec=mec, mew=mew, 
            ecolor=ecol, elinewidth=elw, alpha=falpha)
        ax[sp].errorbar(x=xR, y=yR, yerr=yeR, 
            fmt=mfmt, ms=msz, capsize=csize, label='Rapid', mfc=rcol, mec=mec, mew=mew, 
            ecolor=ecol, elinewidth=elw, alpha=ralpha)
        #ax[sp].text(right,bottom,, ha='bottom',va='right')
        
        #plot labels and legend
        ax[sp].grid(axis='x', linestyle='-', color=gcol, linewidth=gvwidth, alpha=galpha)
        ax[sp].axhline(0, linestyle='-', color=gcol,linewidth=ghwidth, alpha=galpha) 
        if sp == 0 :
            ax[sp].legend(loc='lower left',fancybox=True, shadow=True)
            ax[sp].set_ylabel('North [mm]')
        
        if sp == 1 :
            ax[sp].set_ylabel('East [mm]')
        elif sp == 2 :
            ax[sp].set_ylabel('Vertical [mm]')
            ax[sp].set_xlabel('Date')
            # get errors
            CorrNE,CorrNU,CorrEU=wcorr(dfFNEU)
            
            #errs=np.sqrt(Corr.diagonal())  not error  since slope is removed
            # do last but plot atop
            ax[0].set_title('Station: '+stat+'  Rates: N=%.1f±%.1f, E=%.1f±%.1f, V=%.1f±%.1f [mm/yr]' % (slopes[0], errs[0],slopes[1], errs[1],slopes[2], errs[2]))
            ratesFile.write('%s   %8.4f %9.4f %8.1f  %8.1f %8.1f %8.1f  %8.1f %8.1f %8.1f  %7.4f %7.4f %7.4f    %s  %s\n' % 
                            (stat, lat, lon, height, slopes[0], slopes[1],slopes[2], errs[0],  errs[1], errs[2], CorrNE, CorrNU, CorrEU, sdate, edate)) 
        sp+=1
    f.savefig(os.path.join(plotdir,stat+'_TS.png'), dpi=150, facecolor='white', bbox_inches='tight', pad_inches=0.5)
ratesFile.close()