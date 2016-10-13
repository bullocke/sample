##Strata_Map=raster
##No_Data_Value=number 255


import gdal, ogr, osr, os
import numpy as np
from osgeo import gdal
import sys
import math
from qgis.core import *
from qgis.utils import iface
from PyQt4.QtCore import *
from PyQt4.QtGui import *




path = Strata_Map

ndv = No_Data_Value


gdalData = gdal.Open(path)
if gdalData is None:
  sys.exit( "ERROR: can't open raster" )

ycount=gdalData.RasterYSize
xcount=gdalData.RasterXSize

band_i = gdalData.GetRasterBand(1)
count = {}
total = 0
percent=0
for i in range(0,ycount,ycount/10):
    progress.setPercentage(percent)
    percent+=10
    if (i + ycount/10) < ycount:
        raster_full = band_i.ReadAsArray(0,i,xcount,ycount/10).astype(np.byte)
    else:
        ynew = ycount - i
        raster_full = band_i.ReadAsArray(0,i,xcount,ynew).astype(np.byte)
    classes = np.unique(raster_full)
    classes = classes[classes != ndv]
    classes = classes[classes > 0]
    if np.shape(classes)[0] == 0:
        continue
    for pix in classes:
        try:
            count[pix] = count[pix] + len(raster_full[raster_full == pix])
        except:
            count[pix] = len(raster_full[raster_full == pix])
        total += count[pix]

  # print results sorted by cell_value
for key in sorted(count.iterkeys()):
    classtotal = count[key] / float(total)
    pixels = count[key]
    progress.setText('# Pixels in class {n}: {l}'.format(n=key,l=pixels))
