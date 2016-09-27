#!/usr/bin/env python
""" Generate random sample of high-resolution scenes based on change map

    Author: Eric Bullock
    Date: September 2016
    Version: 1.1

Usage:
    2_FirstStageSamply.py [options] (Neyman | Specified) <shapefile> <output>

Options:
    --allocation <allocation>   Sample allocation [default: 'calc']
    --size <n>                  Sample size for allocation [default: 100]
    -v --verbose                Show verbose debugging messages
    -h --help                   Show help

Allocation (--allocation) "<allocation>" options:
    <specified>                 Comma or space separated list of integers
    calc                        Calculate based on Neyman sampling rules


Notes:

    Output stratified sample of high-resolution scenes based on overlapping
        change map. Scenes are grouped into two strata: The highest concentrated
        change scenes that make up 50% of total change and the lower 50% concentration of
        change.

    If the sample allocation is to be entered manually it must be in the format of
        'low; high', with low corresponding of the proportion of total sample scenes
        to be chosen from low-change scenes and high corresponding to high change
        scenes. For example, if 25 of the high-resolution scenes should correspond with
        Landsat Path Rows with low-density change, and 75 with high-density of change,
        the correct format would be '25; 75'.

Example:

    > ./2_SampleHighRes.py --allocation 'calc' --size 100 Neyman changes.shp sample.shp

"""

from docopt import docopt
from osgeo import ogr
import numpy as np
import os
import sys

#Logging
import logging
VERBOSE = False


ogr.UseExceptions()
ogr.RegisterAll()

logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s',
                    level=logging.INFO,
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


def do_firststage_sample(method, size, allocation, shapefile):
    """Main function for sampling vector tiles """

    #Get total number of change pixels
    total_area = 0
    percent_total_change = []
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(shapefile, 0)
    layer = dataSource.GetLayer()
    num_feat = layer.GetFeatureCount()
    for feature in layer:
        total_area+=feature.GetField("area")

    logger.debug('There are a total of {n} change pixels in the study area (Nh)'.format(n=total_area))
    logger.debug('Total samples that you have specified is {n} (n)'.format(n=size))

    #Get total area
    for feature in range(num_feat):
        count = layer[feature].GetField("area")
        percent = float(count) / total_area
        percent_total_change.append(percent)

    percent_total_change = np.array(percent_total_change)
    sort = np.argsort(percent_total_change)[::-1]
    percent_not_sorted = np.copy(percent_total_change)

    #Sorted by amount of change
    percent_total_change = percent_total_change[sort]

    #Cumulative sum
    cumulative_sum=np.cumsum(percent_total_change)

    #Get indices of the two strata
    strata_1 = sort[np.where( cumulative_sum < .50 )]
    strata_2 = sort[np.where( cumulative_sum >= .50 )]

    #Standard deviation
    sd_strata_1 = np.std(percent_not_sorted[strata_1]) * 100
    sd_strata_2 = np.std(percent_not_sorted[strata_2]) * 100

    #Do allocation
    if method == 'neyman':
        highchange = do_neyman(size, strata_1, strata_2, sd_strata_1, sd_strata_2)
        lowchange = do_neyman(size, strata_2, strata_1, sd_strata_2, sd_strata_1)
    elif allocation == 'equal':
        highchange = size / 2
        lowchange = size / 2

    #Do sampling
    strata_1_selection = np.random.choice(strata_1,size=highchange,replace=False)
    strata_2_selection = np.random.choice(strata_2,size=lowchange,replace=False)


    return strata_1, strata_2, strata_1_selection, strata_2_selection, percent_not_sorted

def write_output(shapefile, output, strata_1, strata_2, high, low, percent_not_sorted):

    #Open shapefile #TODO This is redundant
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(shapefile, 0)
    layer = dataSource.GetLayer()

    #Create memory copy of layer
    outShapefile = output
    outDriver = ogr.GetDriverByName("ESRI Shapefile")

    #TODO: overwrite or no?
    if os.path.exists(outShapefile):
        outDriver.DeleteDataSource(outShapefile)
    outDataSource = outDriver.CreateDataSource(outShapefile)
    srs = layer.GetSpatialRef()
    outLayer = outDataSource.CreateLayer("strata", srs, geom_type=ogr.wkbPolygon)

    inLayerDefn = layer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)

    #Get inclusion probabilities
    inclu1_high = float(high.shape[0]) / strata_1.shape[0]
    inclu1_low = float(low.shape[0]) / strata_2.shape[0]

    strata_field = ogr.FieldDefn("strata", ogr.OFTInteger)
    outLayer.CreateField(strata_field)
    percent_field = ogr.FieldDefn("percent", ogr.OFTReal)
    outLayer.CreateField(percent_field)
    selection_field = ogr.FieldDefn("selection", ogr.OFTInteger)
    outLayer.CreateField(selection_field)
    selection_field = ogr.FieldDefn("pop_stage1", ogr.OFTInteger)
    outLayer.CreateField(selection_field)
    selection_field = ogr.FieldDefn("inclu_1", ogr.OFTReal)
    outLayer.CreateField(selection_field)

    # Get the output Layer's Feature Definition
    outLayerDefn = outLayer.GetLayerDefn()
    # Add features to the ouput Layer
    for i in range(0, layer.GetFeatureCount()):
        # Get the input Feature
        if i in strata_1:
            feature_strata = 1
            inclusion = inclu1_high
            pop = len(strata_1)
        elif i in strata_2:
            feature_strata = 2
            inclusion = inclu1_low
            pop = len(strata_2)
        if (i in high) or (i in low):
            selection = 1
        else:
            continue
        feature_percent = percent_not_sorted[i] * 100



        inFeature = layer.GetFeature(i)
        # Create output Feature
        outFeature = ogr.Feature(outLayerDefn)
        # Add field values from input Layer
        for i in range(0, inLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
        outFeature.SetField('strata',feature_strata)
        outFeature.SetField('percent',feature_percent)
        outFeature.SetField('selection',selection)
        outFeature.SetField('inclu_1',inclusion)
        outFeature.SetField('pop_stage1',pop)
        # Set geometry as centroid
        geom = inFeature.GetGeometryRef()
        outFeature.SetGeometry(geom)
        # Add new feature to output Layer
        outLayer.CreateFeature(outFeature)
        # Close DataSources
    dataSource.Destroy()
    outDataSource.Destroy()

def do_neyman(n, strata1, strata2, sd1, sd2):
    """
    Neyman sampling aims at providing the most precision. #TODO: Elaborate

    Based on Neyman sampling, the idal sample size for stratum h is:

    nh = n * ( Nh * Oh ) / [ Sum ( Ni * Oi ) ]

    Where:
    n = total sample size
    Nh = sample size of stratum h
    Oh = standard deviation of stratum h
    """
    nh = n * (len(strata1) * sd1) / ((len(strata1) * sd1) + (len(strata2) * sd2))
    nh_int = int(round(nh))
    return nh_int

def main():
    """ Read in arguments and test them """
    # Sampling method
    if args['Neyman']:
        method = 'neyman'
    elif args['Equal']:
        method = 'equal'

    # Sample size
    try:
        size = int(args['--size'])
    except:
        logger.error('Sample size must be an integer')
        sys.exit(1)


    # Test if allocation is built-in; if not then it needs to be list of ints
    allocation = args['--allocation']
    if allocation is None:
        if method != 'neyman':
            logger.error('Must specify allocation for designs other than\
                simple random sampling')
            sys.exit(1)
    elif (args['--allocation'].lower() != 'calc') and (method == 'specified'):
        try:
            allocation = np.array([int(i) for i in
                                   allocation.replace(',', ' ').split(' ') if
                                   i != ''])
        except:
            logger.error(
                'Allocation strategy must be built-in method or user must'
                ' specify sequence of integers separated by commas or spaces')
            sys.exit(1)

    shapefile = args['<shapefile>']
    if not os.path.isfile(shapefile):
        logger.error(
            'Specified <shapefile> file {f} does not exist'.format(f=shapefile))
        sys.exit(1)
    logger.debug('Using strata shapefile image {f}'.format(f=shapefile))

    output = args['<output>']

    strata1, strata2, high, low, percent_not_sorted = do_firststage_sample(method, size, allocation, shapefile)

    write_output(shapefile, output, strata1, strata2, high, low, percent_not_sorted)


if __name__ == '__main__':
    args = docopt(__doc__,)

    if args['--verbose']:
        VERBOSE = True
        logger.setLevel(logging.DEBUG)

    main()
