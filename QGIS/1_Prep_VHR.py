##RapidEye_Tile=vector
##Threshold=number .4
##Landcover_Map=raster
##Change_Map=raster
##Strata_Values=string 3;4;5
##No_Data_Values=string 0;255
##Output=output vector

from osgeo import ogr, osr
import numpy as np
import gdal
import os
#import sys
from qgis.core import *
from qgis.utils import iface
from PyQt4.QtCore import *
from PyQt4.QtGui import *



ogr.UseExceptions()
ogr.RegisterAll()
gdal.PushErrorHandler('CPLQuietErrorHandler')


def prep_vhr(changemap, rapideye, output,lcmap, thresh, ndv, strata):
    """ Prepare the VHR tile vector based on a corresponding change map"""

    #Open the change map and vector tiles
    changemap_open, _ = open_raster(changemap)
    lc_open, _ = open_raster(lcmap)
    rapideye_open, _ = open_shapefile(rapideye)

    #Open layer on VHR vector
    rapideyelayer = rapideye_open.GetLayer()

    #Create output strata
    outShapefile = output
    outDriver = ogr.GetDriverByName("ESRI Shapefile")

    #Deleting file if it already exist
    if os.path.exists(outShapefile):
        outDriver.DeleteDataSource(outShapefile)

    outDataSource = outDriver.CreateDataSource(outShapefile)
    srs = rapideyelayer.GetSpatialRef()

    #Create the layer
    outLayer = outDataSource.CreateLayer("strata", srs, geom_type=ogr.wkbPolygon)

    #Copy attributes from Rapid Eye Tile
    inLayerDefn = rapideyelayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)

    #Create new fields:
    ##area: Area of change within individual tile
    ##proportion: Proportion of tile that contains change
    ##ch_pix: Total # of change pixels within tile
    ##noch_pix: Total # of non-change pixels within tile
    area_field = ogr.FieldDefn("area", ogr.OFTInteger)
    prop_field = ogr.FieldDefn("proportion", ogr.OFTReal)
    pixel_field = ogr.FieldDefn("ch_pix", ogr.OFTInteger)
    total_field = ogr.FieldDefn("noch_pix", ogr.OFTInteger)
    outLayer.CreateField(area_field)
    outLayer.CreateField(prop_field)
    outLayer.CreateField(pixel_field)
    outLayer.CreateField(total_field)
    outLayerDefn = outLayer.GetLayerDefn()

    #Total number of tiles in vector file
    totalfeats = len(rapideyelayer)

    itera = 0
    percent = 0
    ten_perc = totalfeats / 10

    #Iterate over features, retrieving zonal statistics
    for i in range(totalfeats):

        if itera == ten_perc:
            percent += 10
            progress.setPercentage(percent)
            itera = 0

        feat = rapideyelayer.GetFeature(i)

        try:
            area, proportion, pix, totalpix = zonal_stats(feat, changemap_open, rapideyelayer, ndv, strata)
            if proportion < thresh:
                itera += 1
                continue
        except:
            itera += 1
            continue
        outFeature = ogr.Feature(outLayerDefn)

         # Add field values from input Layer
        for i in range(0, inLayerDefn.GetFieldCount()):
             outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), feat.GetField(i))

        #Fill zonal statistic fields in output file
        outFeature.SetField('area',area)
        outFeature.SetField('proportion',proportion)
        outFeature.SetField('ch_pix',pix)
        outFeature.SetField('noch_pix', totalpix)

        # Set geometry as centroid
        geom = feat.GetGeometryRef()
        outFeature.SetGeometry(geom)

        # Add new feature to output Layer
        outLayer.CreateFeature(outFeature)

        itera += 1

    #Close and destroy the data source
    changemap_open = None
    rapideye_open.Destroy()
    outDataSource.Destroy()


def zonal_stats(feat, raster, layer, ndv, strata):
    """Perform zonal statistics of vector feature
       on change map"""

    #Get extent information
    transform = raster.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]

    #Pixel size
    pixelWidth = transform[1]
    pixelHeight = transform[5]
    geom = feat.GetGeometryRef()
    if (geom.GetGeometryName() == 'MULTIPOLYGON'):
        count = 0
        pointsX = []; pointsY = []
        for polygon in geom:
            geomInner = geom.GetGeometryRef(count)
            ring = geomInner.GetGeometryRef(0)
            numpoints = ring.GetPointCount()
            for p in range(numpoints):
                lon, lat, z = ring.GetPoint(p)
                pointsX.append(lon)
                pointsY.append(lat)
            count += 1
    elif (geom.GetGeometryName() == 'POLYGON'):
        ring = geom.GetGeometryRef(0)
        numpoints = ring.GetPointCount()
        pointsX = []; pointsY = []
        for p in range(numpoints):
                lon, lat, z = ring.GetPoint(p)
                pointsX.append(lon)
                pointsY.append(lat)

    #Extent of vector feature
    xmin = min(pointsX)
    xmax = max(pointsX)
    ymin = min(pointsY)
    ymax = max(pointsY)

    # Specify offset and rows and columns to read
    xoff = int((xmin - xOrigin)/pixelWidth)
    yoff = int((yOrigin - ymax)/pixelWidth)
    xcount = int((xmax - xmin)/pixelWidth)+1
    ycount = int((ymax - ymin)/pixelWidth)+1

    # Create memory target raster
    target_ds = gdal.GetDriverByName('MEM').Create('', xcount, ycount, gdal.GDT_Byte)
    target_ds.SetGeoTransform((
        xmin, pixelWidth, 0,
        ymax, 0, pixelHeight,
    ))

    # Create for target raster the same projection as for the value raster
    raster_srs = osr.SpatialReference()
    raster_srs.ImportFromWkt(raster.GetProjectionRef())
    target_ds.SetProjection(raster_srs.ExportToWkt())

    # Rasterize zone polygon to raster
    gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[1])

    # Read raster as arrays
    banddataraster = raster.GetRasterBand(1)
    dataraster = banddataraster.ReadAsArray(xoff, yoff, xcount, ycount).astype(np.float)
    bandmask = target_ds.GetRasterBand(1)
    datamask = bandmask.ReadAsArray(0, 0, xcount, ycount).astype(np.float)

    # Mask zone of raster
    zonemask = np.ma.masked_array(dataraster,  np.logical_not(datamask))
    zone_raster_full = np.ma.compressed(zonemask)
    zone_masked = zone_raster_full[np.in1d(zone_raster_full, strata)]

    #Area of change. 1 pixel = 30 X 30 m = 900m^2
    area = len(zone_masked) * 900

    #Proportion of change
    proportion = float(len(zone_masked)) / len(zone_raster_full)
    return area, proportion, len(zone_masked), len(zone_raster_full)


def open_raster(raster):
    """ Open raster file """

    raster_open = gdal.Open(raster)
    if raster_open:
        success = True
    else:
        success = False
    return raster_open, success

def open_shapefile(shapefile):
    """ Open vector file """

    success = False
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(shapefile, 0)
    if dataSource:
        success = True
    return dataSource, success





rapideye = RapidEye_Tile

#Add a Sample ID field incase not all tiles are kept
driver = ogr.GetDriverByName('ESRI Shapefile')
dataSource = driver.Open(rapideye, 1) #1 is read/write

#Create new field for keeping track of sample ID
fldDef = ogr.FieldDefn('SampID', ogr.OFTInteger)

#get layer and add the field:
layer = dataSource.GetLayer()

attributes=[]
inFieldDefn = layer.GetLayerDefn()
for i in range(inFieldDefn.GetFieldCount()):
    attributes.append(inFieldDefn.GetFieldDefn(i).GetNameRef())

if 'SampID' not in attributes:
    layer.CreateField(fldDef)

    sid=0
    for feat in layer:
        feat.SetField('SampID',sid)
        layer.SetFeature(feat)
        sid+=1
dataSource=None

ndv = []
ndvs = No_Data_Values.split(';')
for i in ndvs:
    ndv.append(int(i))

strata = []
stratas = Strata_Values.split(';')
for i in stratas:
    strata.append(int(i))


threshold = float(Threshold)


prep_vhr(Change_Map, RapidEye_Tile, Output, Landcover_Map, threshold,  ndv, strata)
