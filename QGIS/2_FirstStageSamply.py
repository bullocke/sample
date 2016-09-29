##RapidEye_Tile=vector
##Sample_Size=number 100
##Threshold=number .5
##Allocation_Method=selection 'Neyman';'Equal'
##Output=output vector

from osgeo import ogr, osr
import numpy as np
import gdal
import os
import sys
from qgis.core import *
from qgis.utils import iface
from PyQt4.QtCore import *
from PyQt4.QtGui import *


ogr.UseExceptions()
ogr.RegisterAll()



def do_firststage_sample(method, size, shapefile):
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

    progress.setText('There are a total of {n} change pixels in the study area (Nh)'.format(n=total_area))
    progress.setText('Total samples that you have specified is {n} (n)'.format(n=size))

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
    strata_1 = sort[np.where( cumulative_sum < high_threshold )]
    strata_2 = sort[np.where( cumulative_sum >= high_threshold )]

    #Standard deviation
    sd_strata_1 = np.std(percent_not_sorted[strata_1]) * 100
    sd_strata_2 = np.std(percent_not_sorted[strata_2]) * 100
    #Do allocation
    if str(method) == '0':
        highchange = do_neyman(size, strata_1, strata_2, sd_strata_1, sd_strata_2)
        lowchange = do_neyman(size, strata_2, strata_1, sd_strata_2, sd_strata_1)
    elif str(method) == '1':
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


high_threshold = Threshold

strata1, strata2, high, low, percent_not_sorted = do_firststage_sample(Allocation_Method, Sample_Size, RapidEye_Tile, high_threshold)


write_output(RapidEye_Tile, Output, strata1, strata2, high, low, percent_not_sorted)


