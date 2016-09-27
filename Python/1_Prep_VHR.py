#!/usr/bin/env python
"""
    Prep Very-High Resolution (VHR) vector tile data for sampling in the first
    stage of a two-stage cluster design. The VHR tile is stratified based on an
    underlying change map. This script calculates change statistics for each
    vector feature.

    Author: Eric Bullock
    Date: September 2016
    Version: 1.1

Usage:
    1_Prep_VHR.py [options] <lcmap> <changemap> <rapideye> <output>

Options:
    -v --verbose                Show verbose debugging messages
    --thresh <t>                Threshold for proportion area for including in data [default: .4]
    -h --help                   Show help
    -n --ndv <n>                Colon seperated list of no data values (see note below) [default [0;255]]
    <lcmap>                     Input landcover map for determining study area (raster)
    <map>                       Input change map (raster)
    <output>                    Output strata file (vector)

Notes:
    The '--ndv' flag indicates values of the landcover map that are to be treated as
    outside the study region. This should include no data values, and in a multi-class
    landcover product any classes that should be ignored. An example would be a 7-class
    landcover class in which '1' indicates forest. If the units are to be sampled only
    from tiles containing over 40% forest, then the flag should be "--ndv '0;2;3;4;5;6;7;255'"
    assuming 0 and 255 are no data values in the image. The '--thresh' flag should
    then be set to '.4'.

    If a final strata including both change and stable classes has been created, the file
    can be used for both the <lcmap> and <changemap> inputs.

Example:

    > 1_Prep_VHR.py --ndv '0; 255' -v --thresh .4 changemap.shp rapideye.shp strata.shp

"""
from docopt import docopt
from osgeo import ogr, osr
import numpy as np
import gdal
import os
import sys

#Logging
import logging
VERBOSE = False


ogr.UseExceptions()
ogr.RegisterAll()
gdal.PushErrorHandler('CPLQuietErrorHandler')

logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s',
                    level=logging.INFO,
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


def prep_vhr(changemap, rapideye, output,lcmap, thresh, ndv):
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
            logger.debug('{n}% Complete'.format(n=percent))
            itera = 0

        feat = rapideyelayer.GetFeature(i)

        #if int(feat.GetField('TILE_ID')) == 1837116:
            #import pdb; pdb.set_trace()

        try:
            _1, proportion_lc, _2, _3 = zonal_stats(feat, lc_open, rapideyelayer, ndv)
            if proportion_lc < thresh:
                itera += 1
                continue
            area, proportion, pix, totalpix = zonal_stats(feat, changemap_open, rapideyelayer, ndv)
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


def zonal_stats(feat, raster, layer, ndv):
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
    zone_masked = zone_raster_full[~np.in1d(zone_raster_full, ndv)]

    #Area of change. 1 pixel = 30 X 30 m = 900m^2
    area = len(zone_masked) * 900

    #Proportion of change
    proportion = float(len(zone_masked)) / len(zone_raster_full)
    return area, proportion, len(zone_raster_full), len(zone_raster_full)


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

def main():
    """ Read and test inputs """


    method = 'raster'

    #Input map
    changemap = args['<changemap>']
    if not os.path.isfile(changemap):
        logger.error(
            'Specified <map> file {f} does not exist'.format(f=changemap))
        sys.exit(1)
    logger.debug('Using map image {f}'.format(f=changemap))

    changemap = args['<changemap>']
    if not os.path.isfile(changemap):
        logger.error(
            'specified <map> file {f} does not exist'.format(f=changemap))
        sys.exit(1)
    logger.debug('using map image {f}'.format(f=changemap))

    lcmap = args['<lcmap>']
    if not os.path.isfile(changemap):
        logger.error(
            'specified <map> file {f} does not exist'.format(f=changemap))
        sys.exit(1)
    logger.debug('using land cover map image {f}'.format(f=changemap))

    if method == 'raster':
        _, success = open_raster(changemap)
        if not success:
            logger.error('Specified <map> file {f} is not a raster vector file'.format(f=changemap))
            sys.exit(1)

    output = args['<output>']
    rapideye = args['<rapideye>']
    if not os.path.isfile(rapideye):
        logger.error(
            'Specified <map> file {f} does not exist'.format(f=rapideye))
        sys.exit(1)

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

    if args['--ndv']:
        ndv = []
        ndvs = args['--ndv'].split(';')
        for i in ndvs:
            ndv.append(int(i))
    else:
        ndv = [0, 255]

    if args['--thresh']:
        threshold = float(args['--thresh'])
    else:
        threshold = .4

    prep_vhr(changemap, rapideye, output, lcmap, threshold, ndv)

if __name__ == '__main__':
    args = docopt(__doc__,)
    if args['--verbose']:
        VERBOSE = True
        logger.setLevel(logging.DEBUG)

    main()
