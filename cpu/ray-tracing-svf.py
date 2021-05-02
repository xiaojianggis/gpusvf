### This script is used to compute the SVF based on the ray-tracing algorithm
### by Xiaojiang Li, Temple University, May 2, 2021

import os, os.path
import rasterio as rio
from osgeo import gdal
from osgeo.gdalconst import *
from matplotlib import pyplot as plt
import numpy as np
import rasterio
import time
import math


demfile = r'/drive2/thermal-env/data/miami/ground_dsm/row5-col7.tif'

rast = demfile

#open raster layer
src_ds=gdal.Open(rast) 
gt=src_ds.GetGeoTransform()
rb=src_ds.GetRasterBand(1)
arrayDEM = rb.ReadAsArray()
gdal.UseExceptions() #so it doesn't print to screen everytime point is outside grid

# get the size of the raster
rows = arrayDEM.shape[0]
cols = arrayDEM.shape[1]
print ('The size rows is %s, the cols is %s :', rows,cols)


## because the ray-tracing algorithm is very slow, therefore, here only compute the SVF for certain pixels
px = 120
py = 400


for px in range(100, 1000, 10):
    try: #in case raster isnt full extent
        structval=rb.ReadRaster(px,py,1,1,buf_type=gdal.GDT_Float32) #Assumes 32 bit int- 'float'
        intval = struct.unpack('f' , structval) #assume float
        val=intval[0]
    except:
        val=-9999 #or some value to indicate a fail
    
    t0 = time.time()
    
    ## CALCULATE SKY VIEW FACTOR FOR ALL GSV SITES
    SVF_res = 0
    
    # search all 360 angles for each pano site
    for thetaN in range(360):
        rangeDist = 200 # 500m
        #rangeFeet = int(rangeDist*3.28084)  # convert the meter to feet
        
        # create points along each angle
        radiusRange = range(5, rangeDist,1)
        theta = np.pi*thetaN/180
        
        # create an empty beta list to store the betas for each spike
        betaLst = []
        
        # create points along the ray line in 200m or 656 feet, one pixel is one foot
        for radius in radiusRange:
            rayX = int(px + radius*math.cos(theta))
            rayY = int(py + radius*math.sin(theta))
            
            # the corresponding building height is, consider the search region could out of the image
            if rayX >= cols or rayX < 0 or rayY >= rows or rayY < 0:
                continue
            
            # because the ground value is not zero, therefore, need to use the relative height
            buildH = arrayDEM[rayY,rayX] - arrayDEM[py, px]
            
            # if the pixel has its height lower than 2.5m, do not consider anymore
            #if buildH < 6.5:
            #    continue
            
            # considering the GSV pano is captured at height of 2.5m
            #beta = math.atan((buildH - 2.5)*3.28084/radius)
            beta = math.atan(buildH/radius)
            betaLst.append(beta)
        
        if len(betaLst)>0:
            maxBeta = max(betaLst)
        else:
            maxBeta = 0
        SVF_res = SVF_res + math.cos(maxBeta)**2
        
    SVF = SVF_res*1.0/360
    
    print('The svf and the time need is:', SVF, time.time() - t0)
    


