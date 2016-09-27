##Strata_Map=raster
##No_Data_Values=string '0;255'


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

ndv = []
ndvs = No_Data_Values.split(';')
for i in ndvs:
    ndv.append(int(i))


gdalData = gdal.Open(path)
if gdalData is None:
  sys.exit( "ERROR: can't open raster" )

# get number of bands
bands = gdalData.RasterCount

band_i = gdalData.GetRasterBand(1)
raster_full = band_i.ReadAsArray()
raster_flat = raster_full.flatten()
del raster_full
raster_full = None
raster = raster_flat[~np.in1d(raster_flat, ndv)]
del raster_flat
raster_flat = None
  # create dictionary for unique values count
count = {}
classes = np.unique(raster)
total = 0
percent = 0
progress.setPercentage(percent)
total_num = len(classes)
for pix in classes:
    progress.setPercentage(percent)
    count[pix] = len(raster[raster == pix])
    total += count[pix]
    percent += 10

  # print results sorted by cell_value
for key in sorted(count.iterkeys()):
  classtotal = count[key] / float(total)
  pixels = count[key]
  progress.setText('# Pixels in class {n}: {l}'.format(n=key,l=pixels))
  progress.setText('Percent total area class {n}: {l}'.format(n=key,l=classtotal))
