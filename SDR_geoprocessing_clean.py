# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 14:42:32 2021
# author: Eleonore Pierrat
# citation : Global modelling of water consumption impacts on freshwater fish biodiversity in LCA, submitted


this code georprocesses the river basins characteristics (aggregation at the river basin scale)





@author: easpi
"""



import os
from netCDF4 import Dataset
from pcraster import readmap, pcr2numpy
import numpy as np
from numpy import  ma, load,isnan
from matplotlib import pyplot as plt
import tifffile
import pandas as pd
from scipy import stats
import math
#from tabulate import tabulate
# import time

os.chdir('C:/Users/easpi/Documents/PhD Water Footprint/GLAM/code')#where did you put the modules
import Input
import module

# input directory
os.chdir(Input.inputDir)

#%%spatial unit definition


#raster of river basins units
spatial_unit = tifffile.imread(Input.inputDir + '/' + 'basin_5min_pcrglob_adjusted_resized.tif' )#pcrglob basin delineation
#spatial_unit = tifffile.imread(Input.inputDir + '/' + 'gsim_BAF_MD.tif' )#gsim basin delineation - IT HAS A PROBLEM NOW
spatial_unit = ma.masked_values(spatial_unit, -3.40282e+38)#some values are weird, they were filtered from the index array


#array containing the coordinates of the grid cells included in each river basin
pointer_array = np.load(Input.inputDir +'/'+ 'filter_array_catchments_valerio_adjusted.npy', allow_pickle = 1) 
idlist= np.unique(spatial_unit)#does not work, select idlist from the fitler_Array
idlist= idlist[3:,] 

#for gsim sensitivity
#pointer_array = np.load(Input.inputDir +'/'+ 'filter_gsim_basins.npy', allow_pickle = 1)
#idlist = pointer_array[1,]
#idlist = ma.masked_values(idlist,-1.7976931348623157e+308)
#pointer_array = np.asarray(pointer_array[0,], dtype = 'object')
#pointer_array[0][1]


#%%load variables

area = pcr2numpy(readmap(Input.inputDir +'/'+ 'Global_CellArea_m2_05min.map'), mv = 1e20)

#catchment area
area_avg = module.aggregate_sum(area/1e6, pointer_array)#km2

#DEM=elevation m
dem = pcr2numpy(readmap(Input.inputDir +'/'+ 'DEM_05min.map'), mv = 1e20)
dem_agg = module.aggregate_average(dem, pointer_array, area)
#m

#slope
slope = tifffile.imread(Input.inputDir + '/' + 'Slope_DM_WGS84_VB_resized.tif' )
slope= ma.masked_values(slope, -3.40282e+38)
slope_agg = module.aggregate_average(slope, pointer_array, area)#degrees
h= module.plot_percentiles(slope_agg)

#topography index
tan_a = np.tan(slope*math.pi/180)#radians
tan_a = ma.masked_values(tan_a, 0,rtol = 1e-10, atol = 1e-10)
h=module.plot_percentiles(tan_a)

ti_map =np.log(area*1e6/tan_a)
ti_agg = module.aggregate_average(ti_map, pointer_array, area)#careful that the area is in km2!!!
h=module.plot_percentiles(ti_agg)

#discharge

#simulated data including water consumption
d  = Dataset(Input.outputDir + '/'+ 'discharge_outlet_human_1960_2010_pcrglob5.nc')
lat = d.variables['latitude'][:]
lon = d.variables['longitude'][:]
q = d.variables['discharge flow'][:]#m3/s monthly
q = np.append(q, np.full((1,q.shape[1]), 9.96921e+36), axis = 0)
q = ma.masked_values(q, 9.96921e+36)

q_yr = module.convert_month_to_year_sum(q*3600*24*365/12)#in m3/yr. average is wrong!! it has to bethe sum of monthly values to get annual values, this is the monthly average not yearly  average. rerun the whole thing
q_avg = np.mean(q_yr[:,-40:-10], axis = 1)#1970-2000 in m3/yr
q_avg = q_avg/(365*24*3600) #m3/s annual average
h=module.plot_percentiles(q_avg)

#simulated data excluding water consumption (Q natural or Q0)
d  = Dataset(Input.outputDir + '/'+ 'discharge_outlet_natural_1960-2010_pcrglob5.nc')
qn = d.variables['discharge flow'][:]#m3/s monthly
qn = np.append(qn, np.full((1,qn.shape[1]), 9.96921e+36), axis = 0)
qn = ma.masked_values(qn, 9.96921e+36)

qn_yr = module.convert_month_to_year_sum(qn*3600*24*365/12)#in m3/yr. average is wrong!! it has to bethe sum of monthly values to get annual values, this is the monthly average not yearly  average. rerun the whole thing
qn_avg = np.mean(qn_yr[:,-40:-10], axis = 1)#1970-2000 in m3/yr
qn_avg = qn_avg/(365*24*3600) #m3/s annual average
h=module.plot_percentiles(qn_avg)


#precipitation in the catchment 

prec = tifffile.imread(Input.inputDir + '/' + 'wc2.1_5m_bio_12.tif' )#from world clim 2.1
prec = ma.masked_values(prec, -3.4e+38)

prec_pres_agg = module.aggregate_average(prec, pointer_array, area)
#prec_pres_agg = module.aggregate_mean(prec, pointer_array)
h=module.plot_percentiles(prec_pres_agg)



#temperature in the catchment
temp = tifffile.imread(Input.inputDir + '/' + 'wc2.1_5m_bio_1.tif' )
#temp = tifffile.imread('D:/SDR/CHELSA_cur_V1_2B_r5m/5min/bio_1_resized.tif' )
temp = ma.masked_values(temp, -3.4e+38)#some values are weird, they were filtered from the index array
temp = temp + 273.15
temp_pres_agg = module.aggregate_average(temp, pointer_array, area)
h = module.plot_percentiles(temp_pres_agg)



#precipitation in last glacial maximum

prec_21_a = tifffile.imread(Input.inputDir + '/' + 'cclgmbi12.tif' )
prec_21_a = np.append(prec_21_a, np.full((2160-prec_21_a.shape[0],4320), -3.40282e+38), axis = 0)
prec_21_a = ma.masked_values(prec_21_a, -3.40282e+38)#some values are weird, they were filtered from the index array

prec_21_b = tifffile.imread(Input.inputDir + '/' + 'melgmbi12.tif' )
prec_21_b = np.append(prec_21_b, np.full((2160-prec_21_b.shape[0],4320), -3.40282e+38), axis = 0)
prec_21_b = ma.masked_values(prec_21_b, -3.40282e+38)
#some values are weird, they were filtered from the index array

prec_21_c = tifffile.imread(Input.inputDir + '/' + 'mrlgmbi12.tif' )
prec_21_c = np.append(prec_21_c, np.full((2160-prec_21_c.shape[0],4320), -3.40282e+38), axis = 0)
prec_21_c = ma.masked_values(prec_21_c, -3.40282e+38)#some values are weird, they were filtered from the index array

prec_21  = (prec_21_a + prec_21_b + prec_21_c)/3
prec_21 = ma.masked_values(prec_21,0)
prec_21_agg = module.aggregate_average(prec_21, pointer_array, area)
prec_21_agg = ma.masked_invalid(prec_21_agg)
h=module.plot_percentiles(prec_21_agg)

prec_var = (prec - prec_21)/prec_21

prec_21_agg = module.aggregate_average(prec_var, pointer_array, area)
h = module.plot_percentiles(prec_21_agg)

bug = ma.masked_where(prec_21_agg>1, prec_21_agg)#showing only where prec_21_agg is correct
bug_map = module.make_map(bug, pointer_array, idlist, spatial_unit)
module.new_map_netcdf(Input.outputDir + '/'+ 'prec_delta_bug_map', bug_map, 'pred_delta', '1', lat, lon)
plt.imshow(bug_map)#seems ok
plt.imshow(spatial_unit)



#the raster are projected in WGS84 while the filter is expressed in geogrpahical coordinates


#temperature in the last glacial maximum
 
temp_21_a = tifffile.imread(Input.inputDir + '/' + 'cclgmbi1.tif' )
temp_21_a = np.append(temp_21_a, np.full((2160-temp_21_a.shape[0],4320), -3.40282e+38), axis = 0)
temp_21_a = ma.masked_values(temp_21_a, -3.40282e+38)#some values are weird, they were filtered from the index array
temp_21_a = temp_21_a/10 + 273.15

temp_21_b = tifffile.imread(Input.inputDir + '/' + 'melgmbi1.tif' )
temp_21_b = np.append(temp_21_b, np.full((2160-temp_21_b.shape[0],4320), -3.40282e+38), axis = 0)
temp_21_b = ma.masked_values(temp_21_b, -3.40282e+38)#some values are weird, they were filtered from the index array
temp_21_b = temp_21_b/10 + 273.15

temp_21_c = tifffile.imread(Input.inputDir + '/' + 'mrlgmbi1.tif' )
temp_21_c = np.append(temp_21_c, np.full((2160-temp_21_c.shape[0],4320), -3.40282e+38), axis = 0)
temp_21_c = ma.masked_values(temp_21_c, -3.40282e+38)#some values are weird, they were filtered from the index array
temp_21_c = temp_21_c/10 + 273.15


temp_21  = (temp_21_a  + temp_21_b + temp_21_c)/3#variation
temp_21 = ma.masked_values(temp_21, 0)

temp_var = (temp-temp_21)/temp_21
temp_21_agg = module.aggregate_average(temp_var, pointer_array, area)
h = module.plot_percentiles(temp_21_agg)
#kelvins


#climate classification
kg =  tifffile.imread(Input.inputDir + '/' + 'Beck_KG_V1_present_0p083.tif' )
kg = ma.masked_values(kg, 0)#some values are weird, they were filtered from the index array
kg_agg = module.attribute_label(kg, pointer_array)#check alignment, it does not fit perfectly.


out = np.full(kg_agg.shape,0, dtype ='int')

for i in range(len(kg_agg)):
    if kg_agg[i]<= 3:
        out[i]= 1 #'A'
    else:
        if kg_agg[i] <=7 :
            out[i] =2# 'B'
        else:
            if kg_agg[i]<= 16:
                out[i]= 3#'C'
            else:
                if kg_agg[i]<=28:
                    out[i]= 4#'D'
                else:
                    if kg_agg[i]<=30:
                        out[i]= 5#'E'


climate_5 = ma.masked_equal(out,0)


#freshwater ecoregion, habitat, and realm 

fer = tifffile.imread(Input.inputDir + '/' + 'FER.tif' )
fer = ma.masked_values(fer, 0)
ecoregion = module.attribute_label(fer, pointer_array)

#upload attribute table (extracted manually from  https://www.feow.org/)
table = pd.read_excel(Input.inputDir + '/'+'freshwater_ecoregions.xlsx')

out = np.full((ecoregion.shape[0],3), 1e20)

for i in range(len(ecoregion)):
        temp = table.loc[table['ID'] == ecoregion[i]] 
        if temp.empty == False:
            ecoregion_ID = temp.iloc[0,0]
            realm = temp.iloc[0,3]
            hab  = temp.iloc[0,5]
            out[i,0] = ecoregion_ID
            out[i,1]  =realm
            out[i,2] = hab


habitat = ma.masked_values(out, 1e20)


#%%matrix construction q in m3/s,area in km2
data = np.stack((q_avg,qn_avg, area_avg, prec_pres_agg, temp_pres_agg, prec_21_agg, temp_21_agg, dem_agg, ti_agg, slope_agg, habitat[:,0], habitat[:,1],habitat[:,2], climate_5, idlist), axis = 1)

#*365*3600*24/1e9

#no q if using gsim data
#data = np.stack((area_avg, prec_pres_agg, temp_pres_agg, prec_21_agg, temp_21_agg, dem_agg, ti_agg, slope_agg, habitat[:,0], habitat[:,1],habitat[:,2], climate_5, kg_agg, idlist), axis = 1)

#additionnal data for valerio's sensitivity analysis
#data = np.stack((area_avg, prec_pres_agg, temp_pres_agg, prec_21_agg, temp_21_agg, dem, slope, habitat[:,0], habitat[:,1],habitat[:,2], climate_5, kg_agg, idlist), axis = 1)

data = ma.masked_values(data, value = 1e20)
data = ma.masked_invalid(data)
data = data[:-1,]#the last line corresponds to the are not included in the study

# Create the dataframe with the other data-------------------------------------
df = pd.DataFrame(data,columns=[
                        'q',
                        "qnat",
                       'area',
                       'prec',
                       'temp',
                       'prec_delta',
                       'temp_delta',
                       'elevation',
                       'ti',
                       'slope',
                       'ecoregion',
                       'realm',
                       'habitat',
                       'climate5',
                       'id_basin_pcrglob'])
# df = pd.DataFrame(data,columns=[ 
#                       'area',
#                       'prec',
#                       'temp',
#                       'prec_delta',
#                       'temp_delta',
#                       'elevation',
#                       'ti',
#                       'slope',
#                       'ecoregion',
#                       'realm',
#                       'habitat',
#                       'climate5',
#                       'climate30',
#                       'id_basin_gsim'])





#%% #add species richness 
# #not needed for sensitivity
table_sr = pd.read_excel(Input.inputDir + '/'+'pcrglobwb_5min_basins_SR.xlsx')
table_sr.columns =['id_basin_pcrglob',"SR_tot","SR_exo"]

df_join = df.set_index('id_basin_pcrglob').join(table_sr.set_index('id_basin_pcrglob'))
df_join.shape
df_join.head
df_join.columns
#df_join = ma.masked_invalid(df_join)
#df_join = df
df_join.to_csv(Input.outputDir + '/'+ 'dataset_SDR_no_filter_20220316.csv', index = True)



#%% generate map scope




d1=df_join.copy()
threshold_area = module.plot_percentiles(area_avg)[80]
threshold_area

d1.mask((d1['area']<threshold_area), inplace=(True))
d1.mask((d1['climate5']==0), inplace=True)
d1.mask((d1['SR_tot']<1), inplace=True)
d1.mask((d1['q']<=0), inplace=True)
#d1.to_csv(Input.outputDir + '/'+ 'dataset_SDR_filter.csv', index = True)
#df.to_csv(Input.outputDir + '/'+ 'dataset_SDR_GSIM_filter.csv', index = False)


df_join[df_join['q']>0].shape#19579 the rest is null
np.min(df_join['q'])
df_join[df_join['qnat']>0].shape#19577 the rest is null
np.min(df_join['qnat'])
df_join[df_join['SR_tot']>=0].shape#15433 with data, the rest is NA
np.min(df_join['SR_tot'])
df_join[df_join['climate5']>=1].shape#20270 the rest is na
np.min(df_join['climate5'])

df_f = df_join[df_join['area']>threshold_area]#(4064) the rest is small basin
df_f.shape
np.sum(df_f['area'])/np.sum(df_join['area'])

df_f =df_f[df_f['q']>0]#(3931)
df_f=df_f[df_f['SR_tot']>=0]
df_f.shape#3134
df_f=df_f[df_f['climate5']>=1]
df_f.shape
np.sum(df_f['area'])/np.sum(df_join['area'])#88% landcover has valid data

#stats.scoreatpercentile(df['q'],10)

h=module.plot_percentiles(df_f['SR_tot'])#km3/yr
h=module.plot_percentiles(df_f['q'])#m3/s annual average
h=module.plot_percentiles(df_f['qnat'])#m3/s annual average
h=module.plot_percentiles(df_f['area'])
h=module.plot_percentiles(df_f['prec'])
h=module.plot_percentiles(df_f['temp'])
h=module.plot_percentiles(df_f['prec_delta'])
h=module.plot_percentiles(df_f['temp_delta'])
h=module.plot_percentiles(df_f['ti'])
h=module.plot_percentiles(df_f['elevation'])


mask = np.zeros(idlist.shape, dtype = 'float')
len(df_f.index)

for i in df_f.index:#id basin pcrglob inlcuded in the scope
    index  = np.where(idlist == i)[0][0]#in the position where idlistall basins
    mask[index]=1#mask is true 

mask = ma.masked_values(mask, 0)
np.sum(mask)#ok

map_scope = module.make_map(mask, pointer_array, idlist, spatial_unit)#broken
plt.imshow(map_scope)
#export map scope
module.new_map_netcdf(Input.outputDir + '/'+ 'map_scope', map_scope, 'id', '1', lat, lon)


