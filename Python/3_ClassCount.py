# -*- coding: UTF-8 -*-

"""
Get count of classes in raster.
Modified from script from:
https://ssrebelious.wordpress.com/

Usage:
  ClassCount.py <filename> <ndv>

Example:
  3_ClassCount.py ../Material/peru/from_henning/Uso2011.tif '0;255'

"""



from docopt import docopt
import gdal, ogr, osr, os
import numpy as np
from osgeo import gdal
import sys
import math


if __name__ == '__main__':
    args = docopt(__doc__, version='0.1.0')


path = args['<filename>']

ndv = []
ndvs = args['<ndv>'].split(';')
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
for pix in classes:
    count[pix] = len(raster[raster == pix])
    total += count[pix]

  # print results sorted by cell_value
for key in sorted(count.iterkeys()):
  classtotal = count[key] / float(total)
  print "# pixels class %s: %s" %(key, count[key])
  print "Percent total area class %s: %s" %(key, classtotal)
