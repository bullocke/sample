##Landcover_Map=raster
##Change_Map=raster
##Input_Forest_Class=number 2
##Output_Forest_Class=number 5
##Output_NonForest_Class=number 6
##Input_Other_Classes=string '0;0'
##Output_Other_Classes=number 0
##No_Data_ChangeMap=number 255
##No_Data_LCMap=number 255
##Output=output raster

from osgeo import ogr, osr
import numpy as np
import gdal
import os
#import sys
from qgis.core import *
from qgis.utils import iface
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from osgeo import ogr, osr
import numpy as np
import gdal
import os
import sys


ogr.UseExceptions()
ogr.RegisterAll()
gdal.PushErrorHandler('CPLQuietErrorHandler')




def create_strata(changemap, lcmap, ndv, output, forest, nonforest, inforest, lc_ndv, inother, outother):
    """
    Create final change strata with the classes:
    1..n: n number of change classes
    <f>: specified class for stable forest
    <o>: specified class for stable non-forest
    """

    cm_open = gdal.Open(changemap)

    try:
        cm_ar = cm_open.GetRasterBand(1).ReadAsArray()
    except:
        progress.setText('Change map corrupted or not raster file')
        sys.exit(1)
    cm_0 = np.shape(cm_ar)[0]
    cm_1 = np.shape(cm_ar)[1]

    lc_open = gdal.Open(lcmap)

    try:
        lc_ar = lc_open.GetRasterBand(1).ReadAsArray()
    except:
        progress.setText('Change map corrupted or not raster file')
        sys.exit(1)

    lc_0 = np.shape(lc_ar)[0]
    lc_1 = np.shape(lc_ar)[1]

    if np.logical_or(lc_0 != cm_0, lc_1 != cm_1):
        progress.setText('Landcover and change maps must be aligned!')
        sys.exit(1)

    out_ar = np.zeros((lc_0, lc_1))

    #Create some dummy variables for logging
    itera = 0
    percent = 0
    ten_perc = lc_0 / 10
    for i in range(lc_0):

        if itera == ten_perc:
            percent += 10
            progress.setText('{n}% Complete'.format(n=percent))
            itera = 0

        forest_ind = np.logical_and(cm_ar[i,:] == ndv, lc_ar[i,:] == inforest)
        out_ar[i,:][forest_ind] = forest

        nonforest_ind = np.logical_and(cm_ar[i,:] == ndv, lc_ar[i,:] != inforest)
        out_ar[i,:][nonforest_ind] = nonforest

        change_ind = np.where(cm_ar[i,:] != ndv)
        out_ar[i,:][change_ind] = cm_ar[i,:][change_ind]
        change_ind = None

        nodata_ind = np.where(lc_ar[i,:] == lc_ndv)
        out_ar[i,:][nodata_ind] = ndv
        itera += 1

        if inother:
            for o in inother:
                other_ind = np.where(cm_ar[i,:] == o)
                out_ar[i,:][other_ind] = outother

    progress.setText('Writing results')

    #write output
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()
    cols = cm_open.RasterXSize
    rows = cm_open.RasterYSize
    bands = 1
    outDataset = driver.Create(output, cols, rows, bands, gdal.GetDataTypeByName('Int16'))
    geoTransform = cm_open.GetGeoTransform()
    outDataset.SetGeoTransform(geoTransform )
    proj = cm_open.GetProjection()
    outDataset.SetProjection(proj)
    outBand = outDataset.GetRasterBand(1)
    outBand.WriteArray(out_ar, 0, 0)

    cm_open = None
    lc_open = None
    outBand = None
    outDataset = None


changemap = Change_Map
lcmap = Landcover_Map
output = Output
forest = Output_Forest_Class
nonforest = Output_NonForest_Class
inforest = Input_Forest_Class
lc_ndv = No_Data_LCMap
ndv = No_Data_ChangeMap
inothers = Input_Other_Classes.split(';')
inother = []
for i in inothers:
    inother.append(int(i))
if np.all(np.array(inother) == 0):
    inother = None
outother = Output_Other_Classes

create_strata(changemap, lcmap, ndv, output, forest, nonforest, inforest, lc_ndv, inother, outother)

