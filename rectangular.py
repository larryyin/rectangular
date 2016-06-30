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
import re
import shutil
from subprocess import call
from scipy import stats
import time

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

def bearing(lat1,lon1,lat2,lon2):
    latrad1 = np.deg2rad(lat1)
    latrad2 = np.deg2rad(lat2)
    lonrad1 = np.deg2rad(lon1)
    lonrad2 = np.deg2rad(lon2)
    y = np.sin(lonrad2-lonrad1) * np.cos(latrad2)
    x = (np.cos(latrad1)*np.sin(latrad2)-
        np.sin(latrad1)*np.cos(latrad2)*np.cos(lonrad2-lonrad1))
    return np.rad2deg(np.arctan2(y, x))

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

def _det(xvert, yvert):
    '''Compute twice the area of the triangle defined by points with using
    determinant formula.

    Input parameters:

    xvert -- A vector of nodal x-coords (array-like).
    yvert -- A vector of nodal y-coords (array-like).

    Output parameters:

    Twice the area of the triangle defined by the points.

    Notes:

    _det is positive if points define polygon in anticlockwise order.
    _det is negative if points define polygon in clockwise order.
    _det is zero if at least two of the points are concident or if
        all points are collinear.

    '''
    xvert = np.asfarray(xvert)
    yvert = np.asfarray(yvert)
    x_prev = np.concatenate(([xvert[-1]], xvert[:-1]))
    y_prev = np.concatenate(([yvert[-1]], yvert[:-1]))
    return np.sum(yvert * x_prev - xvert * y_prev, axis=0)

class Polygon:
    '''Polygon object.

    Input parameters:

    x -- A sequence of nodal x-coords.
    y -- A sequence of nodal y-coords.

    '''

    def __init__(self, x, y):
        if len(x) != len(y):
            raise IndexError('x and y must be equally sized.')
        self.x = np.asfarray(x)
        self.y = np.asfarray(y)
        # Closes the polygon if were open
        x1, y1 = x[0], y[0]
        xn, yn = x[-1], y[-1]
        if x1 != xn or y1 != yn:
            self.x = np.concatenate((self.x, [x1]))
            self.y = np.concatenate((self.y, [y1]))
        # Anti-clockwise coordinates
        if _det(self.x, self.y) < 0:
            self.x = self.x[::-1]
            self.y = self.y[::-1]

    def is_inside(self, xpoint, ypoint, smalld=1e-12):
        '''Check if point is inside a general polygon.

        Input parameters:

        xpoint -- The x-coord of the point to be tested.
        ypoint -- The y-coords of the point to be tested.
        smalld -- A small float number.

        xpoint and ypoint could be scalars or array-like sequences.

        Output parameters:

        mindst -- The distance from the point to the nearest point of the
                  polygon.
                  If mindst < 0 then point is outside the polygon.
                  If mindst = 0 then point in on a side of the polygon.
                  If mindst > 0 then point is inside the polygon.

        Notes:

        An improved version of the algorithm of Nordbeck and Rydstedt.

        REF: SLOAN, S.W. (1985): A point-in-polygon program. Adv. Eng.
             Software, Vol 7, No. 1, pp 45-47.

        '''
        xpoint = np.asfarray(xpoint)
        ypoint = np.asfarray(ypoint)
        # Scalar to array
        if xpoint.shape is tuple():
            xpoint = np.array([xpoint], dtype=float)
            ypoint = np.array([ypoint], dtype=float)
            scalar = True
        else:
            scalar = False
        # Check consistency
        if xpoint.shape != ypoint.shape:
            raise IndexError('x and y has different shapes')
        # If snear = True: Dist to nearest side < nearest vertex
        # If snear = False: Dist to nearest vertex < nearest side
        snear = np.ma.masked_all(xpoint.shape, dtype=bool)
        # Initialize arrays
        mindst = np.ones_like(xpoint, dtype=float) * np.inf
        j = np.ma.masked_all(xpoint.shape, dtype=int)
        x = self.x
        y = self.y
        n = len(x) - 1  # Number of sides/vertices defining the polygon
        # Loop over each side defining polygon
        for i in range(n):
            d = np.ones_like(xpoint, dtype=float) * np.inf
            # Start of side has coords (x1, y1)
            # End of side has coords (x2, y2)
            # Point has coords (xpoint, ypoint)
            x1 = x[i]
            y1 = y[i]
            x21 = x[i + 1] - x1
            y21 = y[i + 1] - y1
            x1p = x1 - xpoint
            y1p = y1 - ypoint
            # Points on infinite line defined by
            #     x = x1 + t * (x1 - x2)
            #     y = y1 + t * (y1 - y2)
            # where
            #     t = 0    at (x1, y1)
            #     t = 1    at (x2, y2)
            # Find where normal passing through (xpoint, ypoint) intersects
            # infinite line
            t = -(x1p * x21 + y1p * y21) / (x21 ** 2 + y21 ** 2)
            tlt0 = t < 0
            tle1 = (0 <= t) & (t <= 1)
            # Normal intersects side
            d[tle1] = ((x1p[tle1] + t[tle1] * x21) ** 2 +
                       (y1p[tle1] + t[tle1] * y21) ** 2)
            # Normal does not intersects side
            # Point is closest to vertex (x1, y1)
            # Compute square of distance to this vertex
            d[tlt0] = x1p[tlt0] ** 2 + y1p[tlt0] ** 2
            # Store distances
            mask = d < mindst
            mindst[mask] = d[mask]
            j[mask] = i
            # Point is closer to (x1, y1) than any other vertex or side
            snear[mask & tlt0] = False
            # Point is closer to this side than to any other side or vertex
            snear[mask & tle1] = True
        if np.ma.count(snear) != snear.size:
            raise IndexError('Error computing distances')
        mindst **= 0.5
        # Point is closer to its nearest vertex than its nearest side, check if
        # nearest vertex is concave.
        # If the nearest vertex is concave then point is inside the polygon,
        # else the point is outside the polygon.
        jo = j.copy()
        jo[j == 0] -= 1
        area = _det([x[j + 1], x[j], x[jo - 1]], [y[j + 1], y[j], y[jo - 1]])
        mindst[~snear] = np.copysign(mindst, area)[~snear]
        # Point is closer to its nearest side than to its nearest vertex, check
        # if point is to left or right of this side.
        # If point is to left of side it is inside polygon, else point is
        # outside polygon.
        area = _det([x[j], x[j + 1], xpoint], [y[j], y[j + 1], ypoint])
        mindst[snear] = np.copysign(mindst, area)[snear]
        # Point is on side of polygon
        mindst[np.fabs(mindst) < smalld] = 0
        # If input values were scalar then the output should be too
        if scalar:
            mindst = float(mindst)
        return mindst

def kml2polygon(kml,extentW=None,extentS=None,extentE=None,extentN=None):
    with open(kml,'r') as f_poly:
        text_all = f_poly.read().replace('\n', '')
    
    item = text_all.split("</outerBoundaryIs>")[0]
    if "<outerBoundaryIs>" in item:
        testStr = item[item.find("<outerBoundaryIs>")+len("<outerBoundaryIs>"):]
        if ',0.' not in testStr:
            isNear0 = 0
            if ',0' in testStr:
                is3D = 1
            else:
                is3D = 0
        else:
            isNear0 = 1
            if ',0' in testStr:
                is3D = 1
            else:
                is3D = 0
    
    outer_block = []
    
    if (isNear0==0) and (is3D==1):
        stripper = re.compile(r'[^\d.,-;]+')
        for item in text_all.split("</outerBoundaryIs>"):
            if "<outerBoundaryIs>" in item:
                block = (stripper.sub('', item[item.find("<outerBoundaryIs>")+
                         len("<outerBoundaryIs>"):].replace(',0',';'))).rstrip('/').rstrip(';')
                outer_block.append(block)
    elif (isNear0==1) and (is3D==1):
        stripper = re.compile(r'[^\d.,-;]+')
        for item in text_all.split("</outerBoundaryIs>"):
            if "<outerBoundaryIs>" in item:
                block = (stripper.sub('', item[item.find("<outerBoundaryIs>")+
                         len("<outerBoundaryIs>"):].replace(',0.',',999.').replace(',0',';'))).rstrip('/').rstrip(';').replace(',999.',',0.')
                outer_block.append(block)
    elif (is3D==0):
        stripper = re.compile(r'[^\d.,-;]+')
        for item in text_all.split("</outerBoundaryIs>"):
            if "<outerBoundaryIs>" in item:
                block = (stripper.sub('', item[item.find("<outerBoundaryIs>")+
                         len("<outerBoundaryIs>"):].replace(' ',';'))).lstrip(';').rstrip('/').rstrip(';').rstrip('/').rstrip(';')
                outer_block.append(block)
    text_all = None
    
    outer = np.array([np.array([[float(v6) for v6 in v5] for v5 in v4]) for v4 in 
                     [[v3.split(',') for v3 in v2] for v2 in 
                      [v.split(';') for v in outer_block]]])
    outer_block = None
    
    if np.array([extentW,extentS,extentE,extentN]).all():
        extentWS = np.array([extentW,extentS])
        extentEN = np.array([extentE,extentN])
        
        WS = np.array([np.min(v1,axis=0) for v1 in outer])
        EN = np.array([np.max(v1,axis=0) for v1 in outer])
        
        isExtent = np.hstack((WS>extentWS,EN<extentEN)).all(axis=1)
        outer = np.extract(isExtent, outer)
    
    polygons = [Polygon(v[:,0], v[:,1]) for v in outer]
    return polygons

cs_prec = 0.004 # m

#%%
@app.route('/_draft')
def draft():
    """Constructing rectangle..."""
    if not os.path.exists('abc/'):
        os.makedirs('abc/')
        
    abcPath = request.args.get('abcPath', 0, type=str)
    if not abcPath:
        if not os.listdir('abc/'):
            isABC = 0
        else:
            isABC = 1
            abcPath = 'abc/'+os.listdir('abc/')[0]
            print('ABC path: '+abcPath)
    elif not os.path.exists(abcPath):
        status = 'Invalid ABC'
        return jsonify(status=status)
    else:
        isABC = 1
        print('ABC path: '+abcPath)
    
    if isABC:
        ABC = pd.read_csv(abcPath)
        lngA,latA = ABC.iloc[0]
        lngB,latB = ABC.iloc[1]
        lngC,latC = ABC.iloc[2]
        lngO,latO = ABC.iloc[3]
    else:    
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
    
    outPath = 'out/'
    if not os.path.exists(outPath):
        os.makedirs(outPath)
    
    # Save ABC
    s = []
    s += 'lng,lat'
    s += '\n{:f},{:f}'.format(lngA,latA)
    s += '\n{:f},{:f}'.format(lngB,latB)
    s += '\n{:f},{:f}'.format(lngC,latC)
    s += '\n{:f},{:f}'.format(lngO,latO)
    with open(outPath+'ABC', 'w') as f:
        f.writelines(s)
    
    
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
    
    return jsonify(lngA=lngA,latA=latA,
                   lngB=lngB,latB=latB,
                   lngC=lngC,latC=latC,
                   lngO=lngO,latO=latO,
                   lngM=lngM,latM=latM,
                   lngN=lngN,latN=latN,
                   lenI=lenI,lenJ=lenJ,
                   cnI=cnI,cnJ=cnJ)


#%%
@app.route('/_final')
def final():
    """Generating model grid..."""
    run_start = time.time()
    
    if not os.path.exists('abc/'):
        os.makedirs('abc/')
    
    abcPath = request.args.get('abcPath', 0, type=str)
    if not abcPath:
        if not os.listdir('abc/'):
            isABC = 0
        else:
            isABC = 1
            abcPath = 'abc/'+os.listdir('abc/')[0]
            print('ABC path: '+abcPath)
    elif not os.path.exists(abcPath):
        status = 'Invalid ABC'
        return jsonify(status=status)
    else:
        isABC = 1
        print('ABC path: '+abcPath)
    
    if isABC:
        ABC = pd.read_csv(abcPath)
        lngA,latA = ABC.iloc[0]
        lngB,latB = ABC.iloc[1]
        lngC,latC = ABC.iloc[2]
        lngO,latO = ABC.iloc[3]
    else:    
        latO = request.args.get('latO', 0, type=float)
        lngO = request.args.get('lngO', 0, type=float)
        latA = request.args.get('latA', 0, type=float)
        lngA = request.args.get('lngA', 0, type=float)
        latB = request.args.get('latB', 0, type=float)
        lngB = request.args.get('lngB', 0, type=float)
        latC = request.args.get('latC', 0, type=float)
        lngC = request.args.get('lngC', 0, type=float)
    cs = request.args.get('cs', 0, type=float)
    
    #%% Paths
    status = ''
    
    tmpPath = 'tmp/'
    if not os.path.exists(tmpPath):
        os.makedirs(tmpPath)
    else:
        shutil.rmtree(tmpPath)
        os.makedirs(tmpPath)
    
    demPath = request.args.get('demPath', 0, type=str)
    if not demPath:
        demPath = 'dem/'+os.listdir('dem/')[0]
    elif not os.path.exists(demPath):
        status = 'Invalid DEM'
        return jsonify(status=status)
    print('DEM path: '+demPath)
    
    demclippedPath = tmpPath+'dem_clipped.tif'
    
    buildingsPathList = request.args.get('buildingsPath', 0, type=str)
    if not buildingsPathList:
        if not os.listdir('buildings/'):
            isBuilding = 0
        else:
            isBuilding = 1
            buildingsPathList = [('buildings/'+v) for v in os.listdir('buildings/')]
            for v in buildingsPathList:
                print('Buildings path: '+v)
    elif not os.path.exists(buildingsPathList.split(',')[0]):
        status = 'Invalid buildings'
        return jsonify(status=status)
    else:
        isBuilding = 1
        buildingsPathList = buildingsPathList.split(',')
        for v in buildingsPathList:
            print('Buildings path: '+v)
    
    nlcdPath = request.args.get('nlcdPath', 0, type=str)
    if not nlcdPath:
        if not os.listdir('nlcd/'):
            isNLCD = 0
        else:
            isNLCD = 1
            nlcdPath = 'nlcd/'+os.listdir('nlcd/')[0]
            print('NLCD path: '+nlcdPath)
    elif not os.path.exists(nlcdPath):
        status = 'Invalid NLCD'
        return jsonify(status=status)
    else:
        isNLCD = 1
        print('NLCD path: '+nlcdPath)
    
    nlcdclippedPath = tmpPath+'nlcd_clipped.tif'
    
    outPath = 'out/'
    if not os.path.exists(outPath):
        os.makedirs(outPath)
    else:
        shutil.rmtree(outPath)
        os.makedirs(outPath)
    print('Output path: '+outPath)

    # Save ABC
    s = []
    s += 'lng,lat'
    s += '\n{:f},{:f}'.format(lngA,latA)
    s += '\n{:f},{:f}'.format(lngB,latB)
    s += '\n{:f},{:f}'.format(lngC,latC)
    s += '\n{:f},{:f}'.format(lngO,latO)
    with open(outPath+'ABC', 'w') as f:
        f.writelines(s)
    
    #%%
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
    
    lenI = np.sqrt((xA-xM)**2+(yA-yM)**2)
    lenJ = np.sqrt((xB-xM)**2+(yB-yM)**2)

#%% Squarization
    csI = cs
    csJ = cs
    med_H1 = 0
    med_H2 = 0
    loop = 0
    
    loopCount = 0
    while (abs(cs-med_H1)>=cs_prec) and (loopCount<100):
        loop+=1
        loopCount+=1
        if med_H1*med_H2>0:
            csI = csI-(med_H1-cs)/2
        else:
            csI = cs
        
        print("Iteration "+str(loop)+': '+str(med_H1)+' x '+str(med_H2))
        
        if csI>cs_prec*3 and csJ>cs_prec*3:
            cnI = int(np.ceil(lenI/csI))
            cnJ = int(np.ceil(lenJ/csJ))
        else:
            cnI = 0
            cnJ = 0
        
        #
        Jm,Im = np.meshgrid(range(cnJ),range(cnI))
        Jm = Jm+1
        Im = Im+1
        
        Idx = csI*(xM-xA)/lenI
        Idy = csI*(yM-yA)/lenI
        
        xI0 = np.linspace(xA,xA+Idx*cnI,num=cnI,dtype=float)
        yI0 = np.linspace(yA,yA+Idy*cnI,num=cnI,dtype=float)
        
        Jdx = csJ*(xN-xA)/lenJ
        Jdy = csJ*(yN-yA)/lenJ
        
        Xm = np.array([np.linspace(v,v+Jdx*cnJ,num=cnJ,dtype=float) for v in xI0])
        Ym = np.array([np.linspace(v,v+Jdy*cnJ,num=cnJ,dtype=float) for v in yI0])
        
        #
        cn_Xm = center2corners(Xm)
        cn_Ym = center2corners(Ym)
        
        
        # Centers
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
        
        # Corners
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
        
        
        # H1, H2
        I_interim_lonm = cn_lonm[:-1,:]+np.diff(cn_lonm,axis=0)*.5
        I_interim_latm = cn_latm[:-1,:]+np.diff(cn_latm,axis=0)*.5
        J_interim_lonm = cn_lonm[:,:-1]+np.diff(cn_lonm,axis=1)*.5
        J_interim_latm = cn_latm[:,:-1]+np.diff(cn_latm,axis=1)*.5
        
        H1m = dist_greatcircle(J_interim_latm[:-1,:],J_interim_lonm[:-1,:],
                               J_interim_latm[1:,:],J_interim_lonm[1:,:])
        H2m = dist_greatcircle(I_interim_latm[:,:-1],I_interim_lonm[:,:-1],
                               I_interim_latm[:,1:],I_interim_lonm[:,1:])
        
        med_H1 = np.median(H1m.ravel())
        med_H2 = np.median(H2m.ravel())
    
    #
    loopCount = 0
    while (abs(cs-med_H2)>=cs_prec) and (loopCount<100):
        loop+=1
        loopCount+=1
        if med_H2>0:
            csJ = csJ-(med_H2-cs)/2
        else:
            csJ = cs
        
        print("Iteration "+str(loop)+': '+str(med_H1)+' x '+str(med_H2))
        
        if csI>cs_prec*3 and csJ>cs_prec*3:
            cnI = int(np.ceil(lenI/csI))
            cnJ = int(np.ceil(lenJ/csJ))
        else:
            cnI = 0
            cnJ = 0
        
        #
        Jm,Im = np.meshgrid(range(cnJ),range(cnI))
        Jm = Jm+1
        Im = Im+1
        
        Idx = csI*(xM-xA)/lenI
        Idy = csI*(yM-yA)/lenI
        
        xI0 = np.linspace(xA,xA+Idx*cnI,num=cnI,dtype=float)
        yI0 = np.linspace(yA,yA+Idy*cnI,num=cnI,dtype=float)
        
        Jdx = csJ*(xN-xA)/lenJ
        Jdy = csJ*(yN-yA)/lenJ
        
        Xm = np.array([np.linspace(v,v+Jdx*cnJ,num=cnJ,dtype=float) for v in xI0])
        Ym = np.array([np.linspace(v,v+Jdy*cnJ,num=cnJ,dtype=float) for v in yI0])
        
        #
        cn_Xm = center2corners(Xm)
        cn_Ym = center2corners(Ym)
        
        
        # Centers
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
        
        # Corners
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
        
        
        # H1, H2
        I_interim_lonm = cn_lonm[:-1,:]+np.diff(cn_lonm,axis=0)*.5
        I_interim_latm = cn_latm[:-1,:]+np.diff(cn_latm,axis=0)*.5
        J_interim_lonm = cn_lonm[:,:-1]+np.diff(cn_lonm,axis=1)*.5
        J_interim_latm = cn_latm[:,:-1]+np.diff(cn_latm,axis=1)*.5
        
        H1m = dist_greatcircle(J_interim_latm[:-1,:],J_interim_lonm[:-1,:],
                               J_interim_latm[1:,:],J_interim_lonm[1:,:])
        H2m = dist_greatcircle(I_interim_latm[:,:-1],I_interim_lonm[:,:-1],
                               I_interim_latm[:,1:],I_interim_lonm[:,1:])
        
        med_H1 = np.median(H1m.ravel())
        med_H2 = np.median(H2m.ravel())
    
    # ANG
    bearm = bearing(latm,lonm,J_interim_latm[1:,:],J_interim_lonm[1:,:])
    
    degQ4 = bearm>=270
    ANGm = np.zeros(bearm.shape)
    ANGm[degQ4] = 360+90-bearm[degQ4]
    ANGm[~degQ4] = 90-bearm[~degQ4]
    
    print('H1', stats.describe(H1m.ravel()))
    print('H2', stats.describe(H2m.ravel()))
    print('ANG', stats.describe(ANGm.ravel()))
    
    
    #%% DEM
    clipMargin = cs*.00005 # 5 times cell size in meters
    clipW = np.nanmin(cn_lonm.ravel())-clipMargin
    clipS = np.nanmin(cn_latm.ravel())-clipMargin
    clipE = np.nanmax(cn_lonm.ravel())+clipMargin
    clipN = np.nanmax(cn_latm.ravel())+clipMargin
#    clipRes = cs*.00001*.2
    # gdalwarp -te <x_min> <y_min> <x_max> <y_max> input.tif clipped_output.tif
    # gdalwarp -tr 30 30 -r average equals2.tif equals2-averaged_30m.tif
    call(['gdalwarp', 
          '-te', '{:f}'.format(clipW), '{:f}'.format(clipS), 
          '{:f}'.format(clipE), '{:f}'.format(clipN), 
#          '-tr', '{:f}'.format(clipRes), '{:f}'.format(clipRes), '-r', 'bilinear',
          demPath, demclippedPath])
    
    #%% Depth
    print('Extracting depths from DEM...')
    #%% H1, H2 Distance Matrix
    bandNum1 = 1
    DEM = gdal.Open(demclippedPath, GA_ReadOnly )
    band1 = DEM.GetRasterBand(bandNum1)
    
    geotransform = DEM.GetGeoTransform()
    x_ul = geotransform[0]
    y_ul = geotransform[3]
    x_size = geotransform[1]
    y_size = geotransform[5]
    print('DEM cellsize: {xSize:f} x {ySize:f}'.format(xSize=x_size,ySize=y_size))
    
    data_raw = BandReadAsArray(band1)
    (y_cell,x_cell) = data_raw.shape
    xv, yv = meshgrid(range(x_cell), range(y_cell), indexing='xy')
    x_coor = xv * x_size + x_ul + (x_size*.5)
    y_coor = yv * y_size + y_ul + (y_size*.5)
    
    mask_domain = (data_raw>-9999)*(data_raw<9999)
    data = np.copy(data_raw).astype(float)
    data[~mask_domain] = np.nan
    
    band1 = None
    DEM = None
    
#==============================================================================
#     # Cell average scheme
#     depthm = np.empty(lonm.shape)
#     depthm.fill(np.nan)
#     sampleCountm = np.zeros(lonm.shape)
#     
#     for j in range(cnJ):
#         for i in range(cnI):
#             cell_vertices = np.array([[cn_lonm[i,j],cn_latm[i,j]],
#                                       [cn_lonm[i+1,j],cn_latm[i+1,j]],
#                                       [cn_lonm[i+1,j+1],cn_latm[i+1,j+1]],
#                                       [cn_lonm[i,j+1],cn_latm[i,j+1]],
#                                       [cn_lonm[i,j],cn_latm[i,j]]])
#             cell_polygon = Polygon(cell_vertices[:,0], cell_vertices[:,1])
# #            depth_polygon = data[cell_polygon.is_inside(x_coor,y_coor)>0]
#             depth_polygon = data[(cell_polygon.is_inside(x_coor,y_coor)>0)*mask_domain]
#             if len(depth_polygon)>0:
#                 depthm[i,j] = np.mean(depth_polygon)
#                 sampleCountm[i,j] = depth_polygon.size
#             print(depthm[i,j],sampleCountm[i,j])
#     print(depthm)
#     print(sampleCountm)
#==============================================================================
    # Griddata scheme
    distm_lat,distm_lon = latlonlen(latm)
    distm_y,distm_x = latlonlen(y_coor)
    
    depthm = griddata((np.ravel(x_coor*distm_x),np.ravel(y_coor*distm_y)), 
                      np.ravel(data), (np.ravel(lonm*distm_lon),
                      np.ravel(latm*distm_lat)), method='linear')
    depthm = -depthm.reshape(lonm.shape)
    
    #%%
    datumm = np.zeros(lonm.shape)
    
    #%% Buildings
    if isBuilding:
        print("Importing buildings...")
        for buildingsPath in buildingsPathList:
            polygons = kml2polygon(buildingsPath,extentW=clipW,extentS=clipS,
                                   extentE=clipE,extentN=clipN)
            bad_building = 0
            for v in polygons:
                building_polygon = v.is_inside(lonm,latm)>0
                if (building_polygon.sum()/lonm.size)<.2:
                    depthm[building_polygon] = np.nan
                else:
                    bad_building+=1
                    print('Bad building shape #{nBad:d}...'.format(nBad=bad_building))
            print("Imported from {buildingsPath:s}: {nBuilding:d} buildings.".format(buildingsPath=buildingsPath,nBuilding=len(polygons)))
    
    #%% NLCD
    if isNLCD:
        call(['gdalwarp', '-q', 
              '-te', '{:f}'.format(clipW), '{:f}'.format(clipS), 
              '{:f}'.format(clipE), '{:f}'.format(clipN), 
              nlcdPath, nlcdclippedPath])
        
        print('Extracting NLCD classes from raster...')
        
        bandNum1 = 1
        DEM = gdal.Open(nlcdclippedPath, GA_ReadOnly )
        band1 = DEM.GetRasterBand(bandNum1)
        
        geotransform = DEM.GetGeoTransform()
        x_ul = geotransform[0]
        y_ul = geotransform[3]
        x_size = geotransform[1]
        y_size = geotransform[5]
        print('NLCD raster cellsize: {xSize:f} x {ySize:f}'.format(xSize=x_size,ySize=y_size))
        
        data_raw = BandReadAsArray(band1)
        (y_cell,x_cell) = data_raw.shape
        xv, yv = meshgrid(range(x_cell), range(y_cell), indexing='xy')
        x_coor = xv * x_size + x_ul + (x_size*.5)
        y_coor = yv * y_size + y_ul + (y_size*.5)
        
        data = np.copy(data_raw).astype(float)
        
        band1 = None
        DEM = None
        
        # Griddata scheme
        distm_y,distm_x = latlonlen(y_coor)
        
        nlcdm = griddata((np.ravel(x_coor*distm_x),np.ravel(y_coor*distm_y)), 
                          np.ravel(data), (np.ravel(lonm*distm_lon),
                          np.ravel(latm*distm_lat)), method='nearest')
        nlcdm = nlcdm.reshape(lonm.shape)
        nlcdm[np.isnan(depthm)] = np.nan
        
        # NLCD to Manning's
        LC = pd.read_csv('templates/nlcd_table.csv')
        LC_match = list(zip(LC.NLCD.values,LC.Manning.values))
        LC_dict = dict(zip(LC.NLCD.values,LC.Name.values))
        
        manm = np.ones(nlcdm.shape)*.02 # Conservative base value
        for v in LC_match:
            manm[nlcdm==v[0]] = round(v[1],3)
        
        BFRIC_base = 0.0025
    
    #%% Output
    print('Write to output...')
    
    # Write to model_grid_hor
    s = []
    s += "Horizontal Segmentations\n"
    s += "{nI:5d}{nJ:5d}".format(nI=cnI,nJ=cnJ)
    for j in range(1,cnJ-1):
        for i in range(1,cnI-1):
            if ~np.isnan(depthm[i][j]):
                s += "\n{I:5d}{J:5d}{H1:10.2f}{H2:10.2f}{depth:10.3f}{ang:10.2f}{lat:10.6f}{lon:10.6f}{datum:5.1f}".format(I=Im[i][j],J=Jm[i][j],H1=H1m[i][j],H2=H2m[i][j],depth=depthm[i][j],ang=ANGm[i][j],lat=latm[i][j],lon=lonm[i][j],datum=datumm[i][j])  
    with open(outPath+'model_grid_hor', 'w') as f:
        f.writelines(s)

    # Write to corner_loc
    s = []
    for j in range(cnJ+1):
        for i in range(cnI+1):
            s += "{I:5d}{J:5d}{lon:12.6f}{lat:12.6f}{mask:5d}\n".format(I=i+1,J=j+1,lat=cn_latm[i][j],lon=cn_lonm[i][j],mask=1)  
    with open(outPath+'corner_loc', 'w') as f:
        f.writelines(s)
    
    if not isNLCD:
        # Write model_grid to csv
        s = []
        s += "I,J,H1,H2,depth,ang,lat,lon,datum"
        for j in range(1,cnJ-1):
            for i in range(1,cnI-1):
                if ~np.isnan(depthm[i][j]):
                    s += "\n{I:d},{J:d},{H1:.2f},{H2:.2f},{depth:.3f},{ang:.2f},{lat:.6f},{lon:.6f},{datum:.1f}".format(I=Im[i][j],J=Jm[i][j],H1=H1m[i][j],H2=H2m[i][j],depth=depthm[i][j],ang=ANGm[i][j],lat=latm[i][j],lon=lonm[i][j],datum=datumm[i][j])  
        with open(outPath+'model_grid.csv', 'w') as f:
            f.writelines(s)
    else:    
        # Write to bfric2d.inp
        s = []
        s += "NVARBF    BFRIC\n"
        s += "   -1{base:10.5f}\n".format(base=BFRIC_base)
        s += "    I    J     VARBF"
        for j in range(cnJ):
            for i in range(cnI):
                s += "\n{I:5d}{J:5d}{Man:10.5f}".format(I=i+1,J=j+1,Man=manm[i][j])  
        with open(outPath+'bfric2d.inp', 'w') as f:
            f.writelines(s)
        
        # Write model_grid and NLCD/Manning's to csv
        s = []
        s += "I,J,H1,H2,depth,ang,lat,lon,datum,NLCD,Mannings,Land"
        for j in range(cnJ):
            for i in range(cnI):
                if ~np.isnan(depthm[i][j]):
                    s += "\n{I:d},{J:d},{H1:.2f},{H2:.2f},{depth:.3f},{ang:.2f},{lat:.6f},{lon:.6f},{datum:.1f},{NLCD:.0f},{Mannings:.3f},{Land:s}".format(I=Im[i][j],J=Jm[i][j],H1=H1m[i][j],H2=H2m[i][j],depth=depthm[i][j],ang=ANGm[i][j],lat=latm[i][j],lon=lonm[i][j],datum=datumm[i][j],NLCD=nlcdm[i][j],Mannings=manm[i][j],Land=LC_dict[int(nlcdm[i][j])])  
        with open(outPath+'model_grid.csv', 'w') as f:
            f.writelines(s)
    
    # Save to binary
    np.savez('out/bin.npz', 
             cnI=cnI,cnJ=cnJ,Im=Im,Jm=Jm,H1m=H1m,H2m=H2m,depthm=depthm,
             ANGm=ANGm,latm=latm,lonm=lonm,datumm=datumm,nlcdm=nlcdm,manm=manm,
             cn_latm=cn_latm,cn_lonm=cn_lonm)
    
    #%% Stats
    stats_node = cnI*cnJ
    stats_area = np.nansum((H1m*H2m).ravel())
    
    stats_dt1 = np.nanmin((0.5*H1m/np.sqrt(9.80665*(depthm+3))).ravel())
    stats_dt2 = np.nanmin((0.5*H2m/np.sqrt(9.80665*(depthm+3))).ravel())

    stats_lon_max = np.nanmax(lonm.ravel())
    stats_lon_min = np.nanmin(lonm.ravel())
    stats_lat_max = np.nanmax(latm.ravel())
    stats_lat_min = np.nanmin(latm.ravel())
    
    stats_H1_mean = np.nanmean(H1m.ravel())
    stats_H1_median = np.nanmedian(H1m.ravel())
    stats_H1_max = np.nanmax(H1m.ravel())
    stats_H1_min = np.nanmin(H1m.ravel())
    stats_H2_mean = np.nanmean(H2m.ravel())
    stats_H2_median = np.nanmedian(H2m.ravel())
    stats_H2_max = np.nanmax(H2m.ravel())
    stats_H2_min = np.nanmin(H2m.ravel())
    stats_ANG_mean = np.nanmean(ANGm.ravel())
    stats_ANG_median = np.nanmedian(ANGm.ravel())
    stats_ANG_max = np.nanmax(ANGm.ravel())
    stats_ANG_min = np.nanmin(ANGm.ravel())
    stats_depth_mean = np.nanmean(depthm.ravel())
    stats_depth_median = np.nanmedian(depthm.ravel())
    stats_depth_max = np.nanmax(depthm.ravel())
    stats_depth_min = np.nanmin(depthm.ravel())
    
    # Write to stats.txt
    s=[]
    s+='Stats\n'
    s+='Nodes: {:d} x {:d} = {:d}\n'.format(cnI,cnJ,stats_node)
    s+='Extent: {:.6f}, {:.6f}, {:.6f}, {:.6f}\n'.format(stats_lon_min,stats_lat_min,stats_lon_max,stats_lat_max)
    s+='Area: {:.2f} m^2\n'.format(stats_area)
    s+='H1: mean({:.2f}), median({:.2f}), min({:.2f}), max({:.2f})\n'.format(stats_H1_mean,stats_H1_median,stats_H1_min,stats_H1_max)
    s+='H2: mean({:.2f}), median({:.2f}), min({:.2f}), max({:.2f})\n'.format(stats_H2_mean,stats_H2_median,stats_H2_min,stats_H2_max)
    s+='ANG: mean({:.2f}), median({:.2f}), min({:.2f}), max({:.2f})\n'.format(stats_ANG_mean,stats_ANG_median,stats_ANG_min,stats_ANG_max)
    s+='Depth: mean({:.3f}), median({:.3f}), min({:.3f}), max({:.3f})\n'.format(stats_depth_mean,stats_depth_median,stats_depth_min,stats_depth_max)
    s+='Min time step along I: {:.3f} s\n'.format(stats_dt1)
    s+='Min time step along J: {:.3f} s\n'.format(stats_dt2)
    with open(outPath+'stats.txt', 'w') as f:
        f.writelines(s)
    
    #%%
    shutil.rmtree(tmpPath)
    status = 'Job completed'
    
    run_end = time.time()
    print(run_end-run_start)
    print('Job completed successfully.\n')
    
    return jsonify(lngA=lngA,latA=latA,
                   lngB=lngB,latB=latB,
                   lngC=lngC,latC=latC,
                   lngO=lngO,latO=latO,
                   lngM=lngM,latM=latM,
                   lngN=lngN,latN=latN,
                   lenI=lenI,lenJ=lenJ,
                   cnI=cnI,cnJ=cnJ,
                   status=status)

#%%
@app.route('/')
def index():
    return render_template('index.html')

if __name__ == '__main__':
#    app.run(host= '0.0.0.0',port=7110,debug=True)
    app.run(port=7110,debug=True)
