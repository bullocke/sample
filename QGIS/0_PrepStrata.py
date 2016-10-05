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
    cm_band = cm_open.GetRasterBand(1)

    cm_0=cm_open.RasterYSize
    cm_1=cm_open.RasterXSize

    lc_open = gdal.Open(lcmap)
    lc_0=lc_open.RasterYSize
    lc_1=lc_open.RasterXSize
    lc_band = lc_open.GetRasterBand(1)

    if np.logical_or(lc_0 != cm_0, lc_1 != cm_1):
        progress.setText('Landcover and change maps must be aligned!')

    out_ar = np.zeros((lc_0, lc_1), dtype=np.byte)

    #Create some dummy variables for logging
    percent = 0
    progress.setPercentage(percent)
    percent += 10
    block = 100
    itera = 0
    for i in range(0,lc_0,block):
        if (i + block) < lc_0:
            cm_ar = cm_band.ReadAsArray(0,i,cm_1,(block))
            lc_ar = lc_band.ReadAsArray(0,i,cm_1,(block))
            yend = (i + block)
        else:
            ynew = lc_0 - i
            yend = lc_0
            cm_ar = cm_band.ReadAsArray(0,i,cm_1,ynew)
            lc_ar = lc_band.ReadAsArray(0,i,cm_1,ynew)

        if itera >= lc_0 / 10:
            progress.setPercentage(percent)
            percent += 10
            itera = 0
        itera += block
        forest_ind = np.logical_and(cm_ar == ndv, lc_ar == inforest)
        out_ar[i:yend,:][forest_ind] = forest

        nonforest_ind = np.logical_and(cm_ar == ndv, lc_ar != inforest)
        out_ar[i:yend,:][nonforest_ind] = nonforest

        change_ind = np.where(cm_ar != ndv)
        out_ar[i:yend,:][change_ind] = cm_ar[change_ind]

        nodata_ind = np.where(lc_ar == lc_ndv)
        out_ar[i:yend,:][nodata_ind] = ndv

        if inother:
            for o in inother:
                other_ind = np.where(cm_ar == o)
                out_ar[i:yend,:][other_ind] = outother

    progress.setText('Writing Results')

    #write output
    driver = gdal.GetDriverByName('GTiff')
    cols = cm_1
    rows = cm_0
    bands = 1
    outDataset = driver.Create(output, cols, rows, bands, gdal.GDT_UInt16, options = [ 'COMPRESS=PACKBITS' ])
    geoTransform = cm_open.GetGeoTransform()
    outDataset.SetGeoTransform(geoTransform )
    proj = cm_open.GetProjection()
    outDataset.SetProjection(proj)
    outBand = outDataset.GetRasterBand(1)
    outBand.WriteArray(out_ar, 0, 0)
    del cm_open, lc_open, outBand, outDataset

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

