#!/usr/bin/env python

from flask import Flask, jsonify, render_template, request
app = Flask(__name__)

from osgeo import ogr
from osgeo import osr
from osgeo import gdal
from osgeo.gdalnumeric import *
from osgeo.gdalconst import *
import json
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import os
from scipy import stats

#%%
@app.route('/_draft')
def draft():
    """Constructing rectangle..."""
    latO = request.args.get('latO', 0, type=float)
    lngO = request.args.get('lngO', 0, type=float)
    latA = request.args.get('latA', 0, type=float)
    lngA = request.args.get('lngA', 0, type=float)
    latB = request.args.get('latB', 0, type=float)
    lngB = request.args.get('lngB', 0, type=float)
    latC = request.args.get('latC', 0, type=float)
    lngC = request.args.get('lngC', 0, type=float)
#    sf = request.args.get('sf', 0, type=float)
    cs = request.args.get('cs', 0, type=float)
    print(lngO,latO)
    print(lngA,latA)
    print(lngB,latB)
    print(lngC,latC)
#    print(sf,cs)
    
    wgs84 = osr.SpatialReference()
    omerc = osr.SpatialReference()
    wgs84.SetWellKnownGeogCS("WGS84")
    
    omerc.SetHOM2PNO(clat=latO,dfLat1=latA,dfLong1=lngA,dfLat2=latB,dfLong2=lngB,
                     scale=1,fe=0,fn=0)
    wgs2om = osr.CoordinateTransformation(wgs84, omerc)
    om2wgs = osr.CoordinateTransformation(omerc, wgs84)

    om_OABC = ogr.Geometry(ogr.wkbMultiPoint)
    
    om_O = ogr.Geometry(ogr.wkbPoint)
    om_O.AddPoint_2D(lngO, latO)
    om_OABC.AddGeometry(om_O)
    
    om_A = ogr.Geometry(ogr.wkbPoint)
    om_A.AddPoint_2D(lngA, latA)
    om_OABC.AddGeometry(om_A)
    
    om_B = ogr.Geometry(ogr.wkbPoint)
    om_B.AddPoint_2D(lngB, latB)
    om_OABC.AddGeometry(om_B)
    
    om_C = ogr.Geometry(ogr.wkbPoint)
    om_C.AddPoint_2D(lngC, latC)
    om_OABC.AddGeometry(om_C)
    
    om_OABC.Transform(wgs2om)
    om_OABC = json.loads(om_OABC.ExportToJson())['coordinates']
    xOff = om_OABC[0][0]
    yOff = om_OABC[0][1]
    xA = om_OABC[1][0]
    yA = om_OABC[1][1]
    xB = om_OABC[2][0]
    yB = om_OABC[2][1]
    xC = om_OABC[3][0]
    yC = om_OABC[3][1]
    
    xA = xA - xOff
    xB = xB - xOff
    xC = xC - xOff
    yA = yA - yOff
    yB = yB - yOff
    yC = yC - yOff
    xO = 0
    yO = 0
    
    r = np.sqrt((xB-xA)**2+(yB-yA)**2)*.5
    dx = xC-xO
    dy = yC-yO
    dr = np.sqrt(dx*dx+dy*dy)
    D = xO*yC-xC*yO
    sqrtdel = np.sqrt(r*r*dr*dr-D*D)
    
    sgn = (-1)**(int(dy>=0)+1)
    
    x1 = (D*dy+sgn*dx*sqrtdel)/(dr*dr)
    y1 = (-D*dx+np.abs(dy)*sqrtdel)/(dr*dr)
    x2 = (D*dy-sgn*dx*sqrtdel)/(dr*dr)
    y2 = (-D*dx-np.abs(dy)*sqrtdel)/(dr*dr)
    
    d1sq = (x1-xC)**2+(y1-yC)**2
    d2sq = (x2-xC)**2+(y2-yC)**2
    
    if d1sq<d2sq:
        xM = x1
        yM = y1
        xN = x2
        yN = y2
    else:
        xM = x2
        yM = y2
        xN = x1
        yN = y1
    
    xO = xO+xOff
    yO = yO+yOff
    xA = xA+xOff
    xB = xB+xOff
    xC = xC+xOff
    yA = yA+yOff
    yB = yB+yOff
    yC = yC+yOff
    xM = xM+xOff
    xN = xN+xOff
    yM = yM+yOff
    yN = yN+yOff
    
    lenI = np.sqrt((xA-xM)**2+(yA-yM)**2)
    lenJ = np.sqrt((xB-xM)**2+(yB-yM)**2)
    
    if cs>0:
        cnI = int(np.ceil(lenI/cs))
        cnJ = int(np.ceil(lenJ/cs))
    else:
        cnI = 0
        cnJ = 0
    
    wgs_MN = ogr.Geometry(ogr.wkbMultiPoint)
    
    wgs_M = ogr.Geometry(ogr.wkbPoint)
    wgs_M.AddPoint_2D(xM,yM)
    wgs_MN.AddGeometry(wgs_M)
    
    wgs_N = ogr.Geometry(ogr.wkbPoint)
    wgs_N.AddPoint_2D(xN,yN)
    wgs_MN.AddGeometry(wgs_N)
    
    wgs_MN.Transform(om2wgs)
    wgs_MN = json.loads(wgs_MN.ExportToJson())['coordinates']
    lngM = wgs_MN[0][0]
    latM = wgs_MN[0][1]
    lngN = wgs_MN[1][0]
    latN = wgs_MN[1][1]
    
    return jsonify(lngM=lngM,latM=latM,
                   lngN=lngN,latN=latN,
                   lenI=lenI,lenJ=lenJ,
                   cnI=cnI,cnJ=cnJ)


#%%
@app.route('/_final')
def final():
    """Generating model grid..."""
    latO = request.args.get('latO', 0, type=float)
    lngO = request.args.get('lngO', 0, type=float)
    latA = request.args.get('latA', 0, type=float)
    lngA = request.args.get('lngA', 0, type=float)
    latB = request.args.get('latB', 0, type=float)
    lngB = request.args.get('lngB', 0, type=float)
    latC = request.args.get('latC', 0, type=float)
    lngC = request.args.get('lngC', 0, type=float)
    cs = request.args.get('cs', 0, type=float)
    
    status = ''
    
    demPath = request.args.get('demPath', 0, type=str)
    if not demPath:
        demPath = 'in/'+os.listdir('in/')[0]
    elif not os.path.exists(demPath):
        status = 'Wrong DEM'
        return jsonify(status=status)
    outPath = 'out/'
    print('DEM path: '+demPath)
    print('Output path: '+outPath)
    
    wgs84 = osr.SpatialReference()
    omerc = osr.SpatialReference()
    wgs84.SetWellKnownGeogCS("WGS84")
    
    omerc.SetHOM2PNO(clat=latO,dfLat1=latA,dfLong1=lngA,dfLat2=latB,dfLong2=lngB,
                     scale=1,fe=0,fn=0)
    
    wgs2om = osr.CoordinateTransformation(wgs84, omerc)
    om2wgs = osr.CoordinateTransformation(omerc, wgs84)
    
    wgs2om = osr.CoordinateTransformation(wgs84, omerc)
    om2wgs = osr.CoordinateTransformation(omerc, wgs84)
    
    om_OABC = ogr.Geometry(ogr.wkbMultiPoint)
    
    om_O = ogr.Geometry(ogr.wkbPoint)
    om_O.AddPoint_2D(lngO, latO)
    om_OABC.AddGeometry(om_O)
    
    om_A = ogr.Geometry(ogr.wkbPoint)
    om_A.AddPoint_2D(lngA, latA)
    om_OABC.AddGeometry(om_A)
    
    om_B = ogr.Geometry(ogr.wkbPoint)
    om_B.AddPoint_2D(lngB, latB)
    om_OABC.AddGeometry(om_B)
    
    om_C = ogr.Geometry(ogr.wkbPoint)
    om_C.AddPoint_2D(lngC, latC)
    om_OABC.AddGeometry(om_C)
    
    om_OABC.Transform(wgs2om)
    om_OABC = json.loads(om_OABC.ExportToJson())['coordinates']
    xOff = om_OABC[0][0]
    yOff = om_OABC[0][1]
    xA = om_OABC[1][0]
    yA = om_OABC[1][1]
    xB = om_OABC[2][0]
    yB = om_OABC[2][1]
    xC = om_OABC[3][0]
    yC = om_OABC[3][1]
    
    xA = xA - xOff
    xB = xB - xOff
    xC = xC - xOff
    yA = yA - yOff
    yB = yB - yOff
    yC = yC - yOff
    xO = 0
    yO = 0
    
    r = np.sqrt((xB-xA)**2+(yB-yA)**2)*.5
    dx = xC-xO
    dy = yC-yO
    dr = np.sqrt(dx*dx+dy*dy)
    D = xO*yC-xC*yO
    sqrtdel = np.sqrt(r*r*dr*dr-D*D)
    
    sgn = (-1)**(int(dy>=0)+1)
    
    x1 = (D*dy+sgn*dx*sqrtdel)/(dr*dr)
    y1 = (-D*dx+np.abs(dy)*sqrtdel)/(dr*dr)
    x2 = (D*dy-sgn*dx*sqrtdel)/(dr*dr)
    y2 = (-D*dx-np.abs(dy)*sqrtdel)/(dr*dr)
    
    d1sq = (x1-xC)**2+(y1-yC)**2
    d2sq = (x2-xC)**2+(y2-yC)**2
    
    if d1sq<d2sq:
        xM = x1
        yM = y1
        xN = x2
        yN = y2
    else:
        xM = x2
        yM = y2
        xN = x1
        yN = y1
    
    xO = xO+xOff
    yO = yO+yOff
    xA = xA+xOff
    xB = xB+xOff
    xC = xC+xOff
    yA = yA+yOff
    yB = yB+yOff
    yC = yC+yOff
    xM = xM+xOff
    xN = xN+xOff
    yM = yM+yOff
    yN = yN+yOff
    
    lenI = np.sqrt((xA-xM)**2+(yA-yM)**2)
    lenJ = np.sqrt((xB-xM)**2+(yB-yM)**2)
    
    if cs>0:
        cnI = int(np.ceil(lenI/cs))
        cnJ = int(np.ceil(lenJ/cs))
    else:
        cnI = 0
        cnJ = 0
    
    wgs_MN = ogr.Geometry(ogr.wkbMultiPoint)
    
    wgs_M = ogr.Geometry(ogr.wkbPoint)
    wgs_M.AddPoint_2D(xM,yM)
    wgs_MN.AddGeometry(wgs_M)
    
    wgs_N = ogr.Geometry(ogr.wkbPoint)
    wgs_N.AddPoint_2D(xN,yN)
    wgs_MN.AddGeometry(wgs_N)
    
    wgs_MN.Transform(om2wgs)
    wgs_MN = json.loads(wgs_MN.ExportToJson())['coordinates']
    lngM = wgs_MN[0][0]
    latM = wgs_MN[0][1]
    lngN = wgs_MN[1][0]
    latN = wgs_MN[1][1]
    
    #%%
    Jm,Im = np.meshgrid(range(cnJ),range(cnI))
    Jm = Jm+1
    Im = Im+1
    
    Idx = cs*(xM-xA)/lenI
    Idy = cs*(yM-yA)/lenI
    
    xI0 = np.linspace(xA,xA+Idx*cnI,num=cnI,dtype=float)
    yI0 = np.linspace(yA,yA+Idy*cnI,num=cnI,dtype=float)
    
    Jdx = cs*(xN-xA)/lenJ
    Jdy = cs*(yN-yA)/lenJ
    
    Xm = np.array([np.linspace(v,v+Jdx*cnJ,num=cnJ,dtype=float) for v in xI0])
    Ym = np.array([np.linspace(v,v+Jdy*cnJ,num=cnJ,dtype=float) for v in yI0])
    
    #%%
    def center2corners(CENTER):
        CORNERS = np.zeros([CENTER.shape[0]+1,CENTER.shape[1]+1])
        I_diff_half = np.diff(CENTER,axis=0)*.5
        J_diff_half = np.diff(CENTER,axis=1)*.5
        
        I_interim = CENTER[:-1,:]+I_diff_half
        J_interim_diff_half = np.diff(I_interim,axis=1)*.5
        CORNERS[1:-1,1:-1] = I_interim[:,:-1]+J_interim_diff_half
        
        # Sides
        I_W_interim = CENTER[0,:]-I_diff_half[0,:]
        J_W_diff_half = np.diff(I_W_interim)*.5
        CORNERS[0,1:-1] = I_W_interim[:-1]+J_W_diff_half
        
        I_E_interim = CENTER[-1,:]+I_diff_half[-1,:]
        J_E_diff_half = np.diff(I_E_interim)*.5
        CORNERS[-1,1:-1] = I_E_interim[:-1]+J_E_diff_half
        
        I_S_interim = CENTER[:,0]-J_diff_half[:,0]
        J_S_diff_half = np.diff(I_S_interim)*.5
        CORNERS[1:-1,0] = I_S_interim[:-1]+J_S_diff_half
        
        I_N_interim = CENTER[:,-1]+J_diff_half[:,-1]
        J_N_diff_half = np.diff(I_N_interim)*.5
        CORNERS[1:-1,-1] = I_N_interim[:-1]+J_N_diff_half
        
        # Corners
        CORNERS[0,0] = CENTER[0,0]-I_diff_half[0,0]-J_diff_half[0,0]
        CORNERS[-1,0] = CENTER[-1,0]+I_diff_half[-1,0]-J_diff_half[-1,0]
        CORNERS[0,-1] = CENTER[0,-1]-I_diff_half[0,-1]+J_diff_half[0,-1]
        CORNERS[-1,-1] = CENTER[-1,-1]+I_diff_half[-1,-1]+J_diff_half[-1,-1]
        
        return CORNERS
    
    cn_Xm = center2corners(Xm)
    cn_Ym = center2corners(Ym)
    
    
    #%% Centers
    wgs_grid = ogr.Geometry(ogr.wkbMultiPoint)
    
    for iJ in range(cnJ):
        for iI in range(cnI):
            wgs_node = ogr.Geometry(ogr.wkbPoint)
            wgs_node.AddPoint_2D(Xm[iI,iJ],Ym[iI,iJ])
            wgs_grid.AddGeometry(wgs_node)
    wgs_grid.Transform(om2wgs)
    wgs_grid = json.loads(wgs_grid.ExportToJson())['coordinates']
    
    lonm = np.zeros(Xm.shape)
    latm = np.zeros(Ym.shape)
    count = 0
    for iJ in range(cnJ):
        for iI in range(cnI):
            lonm[iI,iJ] = wgs_grid[count][0]
            latm[iI,iJ] = wgs_grid[count][1]
            count+=1
    wgs_grid = []
    
    #%% Corners
    wgs_grid = ogr.Geometry(ogr.wkbMultiPoint)
    
    for iJ in range(cnJ+1):
        for iI in range(cnI+1):
            wgs_node = ogr.Geometry(ogr.wkbPoint)
            wgs_node.AddPoint_2D(cn_Xm[iI,iJ],cn_Ym[iI,iJ])
            wgs_grid.AddGeometry(wgs_node)
    wgs_grid.Transform(om2wgs)
    wgs_grid = json.loads(wgs_grid.ExportToJson())['coordinates']
    
    cn_lonm = np.zeros(cn_Xm.shape)
    cn_latm = np.zeros(cn_Ym.shape)
    count = 0
    for iJ in range(cnJ+1):
        for iI in range(cnI+1):
            cn_lonm[iI,iJ] = wgs_grid[count][0]
            cn_latm[iI,iJ] = wgs_grid[count][1]
            count+=1
    wgs_grid = []
    
    
    #%% H1, H2, ANG
    I_interim_lonm = cn_lonm[:-1,:]+np.diff(cn_lonm,axis=0)*.5
    I_interim_latm = cn_latm[:-1,:]+np.diff(cn_latm,axis=0)*.5
    J_interim_lonm = cn_lonm[:,:-1]+np.diff(cn_lonm,axis=1)*.5
    J_interim_latm = cn_latm[:,:-1]+np.diff(cn_latm,axis=1)*.5
    
    def dist_greatcircle(lat1,lon1,lat2,lon2):
        R = 6371000   # m
        latrad1 = np.deg2rad(lat1)
        latrad2 = np.deg2rad(lat2)
        dLat = latrad2-latrad1
        dLon = np.deg2rad(lon2-lon1)
        a = (np.sin(dLat/2) * np.sin(dLat/2) +
            np.cos(latrad1) * np.cos(latrad2) *
            np.sin(dLon/2) * np.sin(dLon/2))
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
        return R * c
    
    H1m = dist_greatcircle(J_interim_latm[:-1,:],J_interim_lonm[:-1,:],
                           J_interim_latm[1:,:],J_interim_lonm[1:,:])
    H2m = dist_greatcircle(I_interim_latm[:,:-1],I_interim_lonm[:,:-1],
                           I_interim_latm[:,1:],I_interim_lonm[:,1:])
    
    # ANG
    def bear_greatcircle(lat1,lon1,lat2,lon2):
        latrad1 = np.deg2rad(lat1)
        latrad2 = np.deg2rad(lat2)
        lonrad1 = np.deg2rad(lon1)
        lonrad2 = np.deg2rad(lon2)
        y = np.sin(lonrad2-lonrad1) * np.cos(latrad2)
        x = (np.cos(latrad1)*np.sin(latrad2)-
            np.sin(latrad1)*np.cos(latrad2)*np.cos(lonrad2-lonrad1))
        return np.rad2deg(np.arctan2(y, x))
    
    bearm = bear_greatcircle(latm,lonm,J_interim_latm[1:,:],J_interim_lonm[1:,:])
    
    degQ4 = bearm>=270
    ANGm = np.zeros(bearm.shape)
    ANGm[degQ4] = 360+90-bearm[degQ4]
    ANGm[~degQ4] = 90-bearm[~degQ4]
    
    print('H1', stats.describe(H1m.ravel()))
    print('H2', stats.describe(H2m.ravel()))
    print('ANG', stats.describe(ANGm.ravel()))
    
    #%% Depth
    #%% H1, H2 Distance Matrix
    def latlonlen(latdeg):
        lat = np.deg2rad(latdeg)
        m1 = 111132.92;
        m2 = -559.82;
        m3 = 1.175;
        m4 = -0.0023;
        p1 = 111412.84;
        p2 = -93.5;
        p3 = 0.118;
        latlen = m1+(m2*np.cos(2*lat))+(m3*np.cos(4*lat))+(m4*np.cos(6*lat));
        lonlen = (p1*np.cos(lat))+(p2*np.cos(3*lat))+(p3*np.cos(5*lat));
        return (latlen,lonlen) # m
    
    
    # Open the dataset
    bandNum1 = 1
    DEM = gdal.Open(demPath, GA_ReadOnly )
    band1 = DEM.GetRasterBand(bandNum1)
    
    geotransform = DEM.GetGeoTransform()
    x_ul = geotransform[0]
    y_ul = geotransform[3]
    x_size = geotransform[1]
    y_size = geotransform[5]
    
    data_raw = BandReadAsArray(band1)
    (y_cell,x_cell) = data_raw.shape
    xv, yv = meshgrid(range(x_cell), range(y_cell), indexing='xy')
    x_coor = xv * x_size + x_ul + (x_size*.5)
    y_coor = yv * y_size + y_ul + (y_size*.5)
    
    nodata = np.nanmin(data_raw)
    mask_domain = ~(data_raw==nodata)
    data = np.copy(data_raw)
    data[~mask_domain] = np.nan
    
    # Close the datasets
    band1 = None
    DEM = None
    
    distm_lat,distm_lon = latlonlen(latm)
    distm_y,distm_x = latlonlen(y_coor)
    
    depthm = griddata((np.ravel(x_coor*distm_x),np.ravel(y_coor*distm_y)), np.ravel(data), (np.ravel(lonm*distm_lon),np.ravel(latm*distm_lat)), method='linear')
    depthm = -depthm.reshape(lonm.shape)
    
    #%%
    datumm = np.zeros(lonm.shape)
    
    #%%
    # Write to csv
    MG = pd.DataFrame({'I':Im.ravel(),
                       'J':Jm.ravel(),
                       'H1':H1m.ravel(),
                       'H2':H2m.ravel(),
                       'depth':depthm.ravel(),
                       'ang':ANGm.ravel(),
                       'lat':latm.ravel(),
                       'lon':lonm.ravel(),
                       'datum':datumm.ravel()})[['I','J','H1','H2','depth','ang','lat','lon','datum']]
    
    MG.to_csv(outPath+'model_grid.csv', index=False)
    
    # Write to model_grid_hor
    s = []
    s += "Horizontal Segmentations\n"
    s += "{nI:5d}{nJ:5d}".format(nI=cnI,nJ=cnJ)
    for j in range(1,cnJ-1):
        for i in range(1,cnI-1):
            if ~isnan(depthm[i][j]):
                s += "\n{I:5d}{J:5d}{H1:10.2f}{H2:10.2f}{depth:10.2f}{ang:10.2f}{lat:10.6f}{lon:10.6f}{datum:5.1f}".format(I=Im[i][j],J=Jm[i][j],H1=H1m[i][j],H2=H2m[i][j],depth=depthm[i][j],ang=ANGm[i][j],lat=latm[i][j],lon=lonm[i][j],datum=datumm[i][j])  
    
    with open(outPath+'model_grid_hor', 'w') as f:
        f.writelines(s)

    # Write to corner_loc
    s = []
    for j in range(cnJ+1):
        for i in range(cnI+1):
            s += "{I:5d}{J:5d}{lon:12.6f}{lat:12.6f}{mask:5d}\n".format(I=i+1,J=j+1,lat=cn_latm[i][j],lon=cn_lonm[i][j],mask=1)  
    
    with open(outPath+'corner_loc', 'w') as f:
        f.writelines(s)
    
    # 
    status = 'Job completed'
    
    return jsonify(lngM=lngM,latM=latM,
                   lngN=lngN,latN=latN,
                   lenI=lenI,lenJ=lenJ,
                   cnI=cnI,cnJ=cnJ,
                   status=status)


#%%
@app.route('/')
def index():
    return render_template('index.html')

if __name__ == '__main__':
#    app.run(host= '0.0.0.0',debug=True)
    app.run(debug=True)
