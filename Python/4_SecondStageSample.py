#!/usr/bin/env python
""" Generate random sample of high-resolution scenes based on change map

    Author: Eric Bullock
    Date: September 2016
    Version: 1.1

Usage:
    4_SecondStageSample.py [options] (Random | Stratified) <shapefile>  <changemap> <output>

Options:
    --size <size>               Sample size for random allocation [default: 100]
    -v --verbose                Show verbose debugging messages
    -h --help                   Show help
    --allocation <allocation>   Comma-seperated list of sample allocations


Examples:

    > 3_SecondStageSample.py --size 100 Random VHR_sample.shp changemap.tif point_sample.shp

    > 3_SecondStageSample.py --size 200 --allocation '20; 60; 20; 100' Stratified VHR_sample.shp changemap.tif point_sample.shp
"""

from docopt import docopt
from osgeo import ogr, osr
import numpy as np
import gdal
import sys

#Logging
import logging
VERBOSE = False

gdal.UseExceptions()
gdal.AllRegister()
ogr.UseExceptions()
ogr.RegisterAll()

logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s',
                    level=logging.INFO,
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)

def do_point_sample(method, size,shapefile, changemap, output, strata):
    """ Main function for doing point sampling """

    #Open shapefile from first stage sampling
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(shapefile, 0)
    layer = dataSource.GetLayer()

    #Open changemap
    map_ds = gdal.Open(changemap)
    map_ar = map_ds.GetRasterBand(1).ReadAsArray()

    #Create new vector file with samples
    map_sr = osr.SpatialReference()
    map_sr.ImportFromWkt(map_ds.GetProjectionRef())

    # Get OGR driver
    driver = ogr.GetDriverByName('ESRI Shapefile')
    # Create OGR dataset and layer
    sample_ds = driver.CreateDataSource(output)
    out_layer = sample_ds.CreateLayer('sample', map_sr, geom_type=ogr.wkbPolygon)

    # Add fields for layer
    # ID: Sample ID
    out_layer.CreateField(ogr.FieldDefn('ID', ogr.OFTInteger))

    ## Tile: VHR tile
    out_layer.CreateField(ogr.FieldDefn('Tile', ogr.OFTInteger))

    ## Strata1: First stage strata
    out_layer.CreateField(ogr.FieldDefn('Strata1', ogr.OFTInteger))

    ## Strata: Change strata
    out_layer.CreateField(ogr.FieldDefn('Strata', ogr.OFTInteger))

    ## Reference: Reference class
    out_layer.CreateField(ogr.FieldDefn('Reference', ogr.OFTInteger))

    ## TotalPix: Total pixels in strata
    out_layer.CreateField(ogr.FieldDefn('TotalPix', ogr.OFTInteger))

    ## Inclu_1: First-stage inclusion probability
    out_layer.CreateField(ogr.FieldDefn('Inclu_1', ogr.OFTReal))

    ## Inclu_2: Second-stage inclusion probability
    out_layer.CreateField(ogr.FieldDefn('Inclu_2', ogr.OFTReal))

    ## Inclu_Fin: Final inclusion probability
    out_layer.CreateField(ogr.FieldDefn('Inclu_Fin', ogr.OFTReal))

    ## TilePop: First-stage strata population
    out_layer.CreateField(ogr.FieldDefn('TilePop', ogr.OFTInteger))
    #Sample ID
    pix_id = 1

    if method == 'stratified':

                #Return changemap in chosen vector tiles
                changemap_mask, tile_ar = extract_alltiles(layer, changemap)

                #Sample the selected scenes
                strata, sample_y, sample_x, total, inclu2, classes = sample_stratified(size, changemap_mask, strata)

                #Return first-stage inclusion probabilities
                inclu1 = get_first_inclusion(layer)

                #Write output
                _, _ = write_vector_output(sample_x, sample_y,
                                    map_ds,out_layer, map_ar,
                                    output, classes, pix_id, tile_ar, total, inclu2, inclu1,
                                    ogr_frmt='ESRI Shapefile')
    elif method == 'random':

        #Loop over features
        for feat in range(layer.GetFeatureCount()):

            feature = layer.GetFeature(feat)
            tile = feature.GetField("SampID")
            inclu1 = feature.GetField("inclu_1")
            if feature.GetField("selection") == 0:
                continue

            #Extract indices of pixels within sampled vector tiles
            _, barray, xoff, yoff, xcount, ycount = extract_tile(feature, changemap, layer)

            #Sample the raster indices
            sample_y, sample_x, total, inclu2 = sample_random(size, xoff, yoff, xcount, ycount)

            #Write output
            pix_id, out_layer = write_vector_output(sample_x, sample_y,
                                                       map_ds,out_layer, map_ar,
                                                    output,None, pix_id, tile,
                                                    total, inclu2, inclu1,
                                                ogr_frmt='ESRI Shapefile')



def get_first_inclusion(layer):
    """ Return first-stage inclusion probability for each vector tile """

    inclu1 = []
    tiles = []
    tilepop = []
    combined = []
    strata1 = []
    for i in range(layer.GetFeatureCount()):
        feat = layer.GetFeature(i)
        inclu1.append(feat.GetField('inclu_1'))
        tiles.append(feat.GetField('SampId'))
        tilepop.append(feat.GetField('pop_stage1'))
        strata1.append(feat.GetField('strata'))
    inclu1 = np.array(inclu1)
    tiles = np.array(tiles)
    population = np.array(tilepop)
    first_strata = np.array(strata1)
    combined.append(tiles)
    combined.append(inclu1)
    combined.append(population)
    combined.append(first_strata)

    return combined

def sample_stratified(size, changemap_mask, strata):
    """ perform random stratified sample of map
        Modified from code by Chris Holden
        https://github.com/ceholden/misc """

    # Find map classes within image
    classes = np.sort(np.unique(changemap_mask))
    # Exclude masked values
    mask = [0, 255]
    classes = classes[~np.in1d(classes, mask)]

    counts = np.array(strata)
    class_counts = []
    inclu2 = []
    for a, i in enumerate(classes):
        px = np.sum(changemap_mask == i)
        class_counts.append(px)
        inclu2.append(counts[a]/float(px))
    inclu2 = np.array(inclu2)
    logger.debug('Found {n} classes'.format(n=classes.size))

    if classes.size != counts.size:
        raise ValueError(
            'Sample counts must be given for each unmasked class in map')

    # Initialize outputs
    strata = np.array([])
    rows = np.array([])
    cols = np.array([])

    logger.debug('Performing sampling')

    for c, n in zip(classes, counts):
        logger.debug('Sampling class {c}'.format(c=c))

        # Find pixels containing class c
        row, col = np.where(changemap_mask == c)

        # Check for sample size > population size
        if n > col.size:
            logger.warning(
                'Class {0} sample size larger than population'.format(c))
            logger.warning('Reducing sample count to size of population')

            n = col.size

        samples = np.random.choice(col.size, n, replace=False)

        logger.debug('    collected samples')

        strata = np.append(strata, np.repeat(c, n))
        rows = np.append(rows, row[samples])
        cols = np.append(cols, col[samples])

    return strata, rows, cols, class_counts, inclu2, classes


def extract_alltiles(layer, changemap):
    """ Extract changemap for areas sampled in first-stage of sampling """

    ch_open = gdal.Open(changemap)
    ch_ar = ch_open.GetRasterBand(1).ReadAsArray()
    mask_ar = np.zeros_like(ch_ar)

    #Create matching tile array to save corresponding tile ID
    tile_ar = np.zeros((np.shape(mask_ar)[0],np.shape(mask_ar)[1]))

    for feat in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(feat)
        if feature.GetField("selection") == 1:
            #Extract subset for feature
            _, barray, xoff, yoff, xcount, ycount = extract_tile(feature, changemap, layer)

            #File in new array with subset
            mask_ar[yoff:(yoff+ycount), xoff:(xoff+xcount)] = barray

            #File in matching tile array in sample ID
            tile_ar[yoff:(yoff+ycount), xoff:(xoff+xcount)] = feature.GetField("SampID")

    return mask_ar, tile_ar

def extract_tile(feat, raster_file, layer):
    """Subset changemap for vector feature
       Modified from the Python GDAL/Ogr Cookbook
       https://pcjericks.github.io/py-gdalogr-cookbook/ """

    raster = gdal.Open(raster_file)
    transform = raster.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
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

    #Get extent of feature
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
    zonemask = np.ma.masked_array(dataraster, np.logical_not(datamask))
    zonal_data = zonemask.data
    return target_ds, zonal_data, xoff, yoff, xcount, ycount


def sample_random(size, xoff, yoff, xcount, ycount):
    """ Perform random sample """

    #Buffer 50 pixels to avoid edges
    xoff += 25
    yoff += 25
    xcount -= 50
    ycount -= 50
    sample_x = np.random.choice(np.arange(xoff,(xoff+xcount)), size, replace=False)
    sample_y = np.random.choice(np.arange(yoff,(yoff+ycount)), size, replace=False)
    total = xcount * ycount
    inclu2 = float(size) / total
    return sample_y, sample_x, total, inclu2

def write_vector_output(cols, rows, map_file,layer,changemap, output,
                        classes, pix_id, tiles, _total, _inclu2, _inclu1, ogr_frmt='ESRI Shapefile'):
    """ Write output row and column indices to point samples
        Modified from code by Chris Holden
        https://github.com/ceholden/misc """

    # Corners of pixel in pixel coordinates
    corners = [(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]

    gt = map_file.GetGeoTransform()
    # Loop through samples adding to layer
    for i, (col, row) in enumerate(zip(cols, rows)):
        # Feature
        if type(tiles) == int:
            tile = tiles
            strata = int(changemap[row,col])
            inclu2 = _inclu2
            inclu1 = _inclu1
            total = _total
            #total = _TilePop #TODO
        else:
            tile = int(tiles[row, col])
            strata = int(changemap[row, col])
            inclu2 = _inclu2[classes == strata][0]
            inclu2 = float(inclu2)
            total = np.array(_total)[classes == strata][0]
            total = int(total)
            inc_id = np.where(tile == _inclu1[0])[0]
            inclu1 = float(_inclu1[1][inc_id])
            tilepop = int(_inclu1[2][inc_id])
            strata1 = int(_inclu1[3][inc_id])

        final_inclusion = inclu2 * inclu1
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetField('ID', pix_id)
        feature.SetField('Tile', tile)
        feature.SetField('TilePop', tilepop)
        feature.SetField('Strata1', strata1)
        feature.SetField('Strata', strata)
        feature.SetField('TotalPix', total)
        feature.SetField('Inclu_1', inclu1)
        feature.SetField('Inclu_2', inclu2)
        feature.SetField('Inclu_Fin', final_inclusion)

        # Geometry
        ring = ogr.Geometry(type=ogr.wkbLinearRing)

        for corner in corners:
            ring.AddPoint(
                gt[0] + (col + corner[0]) * gt[1] + (row + corner[1]) * gt[2],
                gt[3] + (col + corner[1]) * gt[4] + (row + corner[1]) * gt[5])
        square = ogr.Geometry(type=ogr.wkbPolygon)
        square.AddGeometry(ring)

        feature.SetGeometry(square)

        layer.CreateFeature(feature)

        feature.Destroy()

        pix_id += 1

    tiles = None
    return pix_id, layer


def main():
    """ Read in arguments and test them """

    # Sample size
    try:
        size = int(args['--size'])
    except:
        logger.error('Sample size must be an integer')
        sys.exit(1)

    # Sampling method
    if args['Random']:
        method = 'random'
        strata = None
    elif args['Stratified']:
        method = 'stratified'
        if not args['--allocation']:
            logger.error('Must specify sample allocation for stratified sample')
            sys.exit(1)
        strata = []
        allocation = args['--allocation'].split(';')
        for i in allocation:
            strata.append(int(i))
        if sum(strata) != size:
            logger.error('Sample size must equal strata allocation')
            sys.exit(1)


    changemap = args['<changemap>']

    shapefile = args['<shapefile>']
    output = args['<output>']

    do_point_sample(method, size, shapefile, changemap, output, strata)


if __name__ == '__main__':
    args = docopt(__doc__,)

    if args['--verbose']:
        VERBOSE = True
        logger.setLevel(logging.DEBUG)

    main()
