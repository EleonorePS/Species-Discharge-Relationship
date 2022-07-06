# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 10:09:30 2021

@author: easpi
"""

#import os
from netCDF4 import Dataset
from pcraster import readmap, pcr2numpy
import numpy as np
from numpy import isnan, ma, trapz
from matplotlib import pyplot as plt
from scipy import stats
import tifffile
import Input


def weight(pointer_array, idlist, spatial_unit, area):
    n_spatial_unit = pointer_array.shape[0]
    out = np.full((n_spatial_unit), 1e20)
    for i in range(n_spatial_unit):#select catchment in the same order as the pointer array sot hat the idlist is still valid
        out[i] = np.sum(area[pointer_array[i][0], pointer_array[i][1]])
    out = ma.masked_where(isnan(out), out)

    return out

def exclude(s, mask):
    out = np.full(s.shape, 1e20)
    for  i in range(s.shape[0]):
        if mask[i]==0:
            out[i,:] = s[i,:]
    out = ma.masked_values(out, 1e20)
    return out

def index_filter(spatial_unit):
        """
        Parameters
        ----------
        spatial_unit : TYPE array (lat,lon)
            DESCRIPTION. array containing the ID of the catchments starting with ID = 0
    
        Returns
        -------
        index_filter : TYPE tuple ([spatial unit id, (lat list,lon list)], [id list])
        idlsit : TYPE array where the ID of the catchment is stored in the same position  as index_filter row number
            DESCRIPTION. the list contains for each catchmetn ID, the list of latitude index in position 0 and the list of longitude indexes for position 1
    
        """
        index_filter = []
    # select ID of spatial units present in the map
        idlist = np.unique(spatial_unit)
        for n in idlist:
            a = np.where(spatial_unit == n)
            index_filter.append(a)
        pointer_array = np.array(index_filter[0], dtype=list)
        return pointer_array, idlist

 
def convert_month_to_year_avg(s):
    """
    Parameters
    ----------
    s : TYPE scalar array of shape (spatial unit id, time)
        DESCRIPTION.it represents the scalar variable timseries at monthly timestep

    Returns
    -------
    out : TYPE scalar array (spatial unit id, time)
        DESCRIPTION. converted array to yearly timestep (average 12 months)

    """

    out = np.zeros( (s.shape[0], s.shape[1]//12) )
    for i in range(s.shape[0]):#for all catchments
        for j in range(out.shape[1]):#for all year timesteps
            out[i,j] = np.mean(s[i, j*12:j*12+12])
    #out_mask = ma.masked_where(isnan(s), out)
    out_mask = ma.masked_where(isnan(out) == 1, out)
    return out_mask
                   
def convert_month_to_year_sum(s):
    """
    Parameters
    ----------
    s : TYPE scalar array of shape (spatial unit id, time)
        DESCRIPTION.it represents the scalar variable timseries at monthly timestep

    Returns
    -------
    out : TYPE scalar array (spatial unit id, time)
        DESCRIPTION. converted array to yearly timestep (average 12 months)

    """

    out = np.zeros( (s.shape[0], s.shape[1]//12) )
    for i in range(s.shape[0]):#for all catchments
        for j in range(out.shape[1]):#for all year timesteps
            out[i,j] = np.sum(s[i, j*12:j*12+12])
    #out_mask = ma.masked_where(isnan(s), out)
    out_mask = ma.masked_where(isnan(out) == 1, out)
    return out_mask

def moving_average(a, n):
    #https://stackoverflow.com/questions/14313510/how-to-calculate-rolling-moving-average-using-numpy-scipy
    cumsum = np.cumsum(a, axis = 1)
    out = np.full((a.shape[0], a.shape[1]-n+1),1e20)
    for i in range(a.shape[0]):
        ret = cumsum[i,:]
        ret[n:] =ret[n:]-ret[:-n]
        out[i,:] = ret[n-1:]/n
    out = ma.masked_values(out, 1e20)
    return out


def make_map(stressor_aggregated, index_filter, idlist, spatial_unit):#correct that taking into acount the correct location of the catchment 
    """array inputs! spatial unit with mask"""
    """indexfilter has to be calculated based on spatial_unit beforehand"""
    l1=idlist
    # l1 = l1[np.where(isnan(l1) == 0)]
    # l = ma.getdata(l1).tolist()
    map_temp = np.full(spatial_unit.shape, 1e20)
       
   #  for i in l1:#for all catchment id 
   #      map_temp[index_filter[l1==i][0][0], index_filter[l1==i][0][1]] = stressor_aggregated[l1==i]
   #      #fill the map with the value associated with each catchment.
    for i in l1:#for all catchment id 
        a = np.where(l1 == i)[0][0]#verify position of the catchment in the id list
        if stressor_aggregated.mask[a]==0:
            map_temp[index_filter[a][0], index_filter[a][1]] = stressor_aggregated[a]
        #fill the map with the value associated with each catchment.
   
    #map_out = ma.masked_where(isnan(map_temp) == 1, map_temp)
    #map_out = ma.masked_values(map_temp, 1e20)
    map_out = ma.masked_where(isnan(spatial_unit) == 1, map_temp)
    map_out = ma.masked_values(map_out, 1e20)
    return map_out

    
def new_map_netcdf(filename, map_out, name_scalar, unit_scalar, lat, lon):
    """
    Parameters
    ----------
   
    Returns
    -------
    None. Netcdf file is written based on the inputs. and the dataset is 
    returned to main program.

    """
    dataset_out = Dataset(filename +'.nc','w',format = 'NETCDF4')
    dataset_out.createDimension("latitude",lat.shape[0])
    dataset_out.createDimension("longitude", lon.shape[0])
    

    latitude = dataset_out.createVariable("latitude","f8",("latitude",), zlib = True)
    longitude = dataset_out.createVariable("longitude","f8",("longitude",), zlib = True)
    scalar = dataset_out.createVariable(name_scalar,"f8",("latitude","longitude"), zlib = True)
    scalar.units = unit_scalar

    #fill NETCDF with results
    latitude[:] = lat[:]
    longitude[:] = lon[:]
    scalar[:] = map_out[:]

    dataset_out.sync()#write into the  saved file
    print(dataset_out)
    dataset_out.close()
    return "netcdf created, check folder"

def new_stressor_out_netcdf(filename, stressor_out, ID, t, name_stressor, unit_time,unit_stressor):
    """
    Parameters
    ----------
   
    Returns
    -------
    None. Netcdf file is written based on the inputs. and the dataset is 
    returned to main program.

    """
    dataset_out = Dataset(filename +'.nc','w',format = 'NETCDF4')
    dataset_out.createDimension("spatial_unit_id",ID.shape[0])
    dataset_out.createDimension("time", t.shape[0])
    

    time = dataset_out.createVariable("time","f4",("time",), zlib = True)
    time.units = unit_time
    spatial_unit_id = dataset_out.createVariable("spatial_unit_id","f4",("spatial_unit_id",), zlib = True)
    stressor_aggregated_timeserie = dataset_out.createVariable(name_stressor,"f4",("spatial_unit_id","time"), zlib = True)
    stressor_aggregated_timeserie.units = unit_stressor

    #fill NETCDF with results
    time[:] = t[:]
    spatial_unit_id[:] = ID[:]
    stressor_aggregated_timeserie [:] = stressor_out[:]

    dataset_out.sync()#write into the  saved file
    print(dataset_out)
    dataset_out.close()
    return "netcdf created, check folder"




def aggregate_Q(idlist,pointer_array):


#the flow accumulation file is used to identify the outlet point
    flowacc = tifffile.imread(Input.inputDir + '/' + 'flowAcc_5min.tif')
    flowacc = ma.masked_equal(flowacc, -2147483647)

    d = Dataset(Input.inputDir + '/' + Input.fn_discharge)
    #area = pcr2numpy(readmap(Input.inputDir + '/' + Input.fn_area_map), mv = 1e20)

    discharge = d.variables["discharge"]#m3/s
    times = d.variables["time"][:]#month
        
#aggrgeate 
    ntime, nlat, nlon = discharge.shape
    n_spatial_unit = pointer_array.shape[0]
    s_aggregated = np.full((n_spatial_unit, ntime), 1e20)
    
    for t in range(ntime):
        temp = discharge[t,:,:]
        for k in range(n_spatial_unit):
            coord = np.argmax(flowacc[pointer_array[k][0],pointer_array[k][1]], axis=None)#returns the index of the max value of the flattenned array.
            s_aggregated[k,t] = np.ravel(temp[pointer_array[k][0],pointer_array[k][1]])[coord]#select the value at the coordinate point
            #s_aggregated[k,t] = np.max(temp[pointer_array[k][0],pointer_array[k][1]])
    
    s_aggregated = ma.masked_where(isnan(s_aggregated), s_aggregated)
    
    ID = ma.getdata(idlist[:-1])

    new_stressor_out_netcdf(Input.outputDir + '/'+ 'discharge_outlet_human_'  + Input.name_timeperiod, s_aggregated[:-1,:], ID, times, 'discharge flow', 'month', 'm3/s') 

    world_Q = np.sum(s_aggregated, axis = 0)
    fig, ax = plt.subplots()  # Create a figure and an axes.
    ax.plot(world_Q, label='q human')  # Plot some data on the axes.
    ax.set_xlabel('year')  # Add an x-label to the axes.
    ax.set_ylabel('m3/s')  # Add a y-label to the axes.
    ax.set_title("World discharge at outlet 1960/2010")  # Add a title to the axes.
    ax.legend()  # Add a legend.
    
    return s_aggregated


def aggregate_max_Q(idlist,pointer_array):


    d = Dataset(Input.inputDir + '/' + Input.fn_discharge)
    #area = pcr2numpy(readmap(Input.inputDir + '/' + Input.fn_area_map), mv = 1e20)

    discharge = d.variables["discharge"]#m3/s
    times = d.variables["time"][:]#month
        
#aggrgeate 
    ntime, nlat, nlon = discharge.shape
    n_spatial_unit = pointer_array.shape[0]
    s_aggregated = np.full((n_spatial_unit, ntime), 1e20)
    
    for t in range(ntime):
        temp = discharge[t,:,:]
        for k in range(n_spatial_unit):
            s_aggregated[k,t] = np.max(temp[pointer_array[k][0],pointer_array[k][1]])#select the value at the coordinate point
            #s_aggregated[k,t] = np.max(temp[pointer_array[k][0],pointer_array[k][1]])
    
    s_aggregated = ma.masked_where(isnan(s_aggregated), s_aggregated)
    
    ID = ma.getdata(idlist[:-1])

    new_stressor_out_netcdf(Input.outputDir + '/'+ 'discharge_max_human_'  + Input.name_timeperiod, s_aggregated[:-1,:], ID, times, 'discharge flow', 'month', 'm3/s') 

    world_Q = np.sum(s_aggregated, axis = 0)
    fig, ax = plt.subplots()  # Create a figure and an axes.
    ax.plot(world_Q, label='q human')  # Plot some data on the axes.
    ax.set_xlabel('year')  # Add an x-label to the axes.
    ax.set_ylabel('m3/s')  # Add a y-label to the axes.
    ax.set_title("World discharge at outlet 1960/2010")  # Add a title to the axes.
    ax.legend()  # Add a legend.
    
    return s_aggregated


def aggregate_sum(stressor, pointer_array):
    nlat, nlon = stressor.shape
    n_spatial_unit = pointer_array.shape[0]
    out = np.full((n_spatial_unit), 1e20)
    
    
    for k in range(n_spatial_unit):
        s = stressor[pointer_array[k][0],pointer_array[k][1]]
        #s = stressor[pointer_array[k,0][0],pointer_array[k,0][1]]
        out[k] = np.sum(s)
    
    out = ma.masked_where(isnan(out), out)
    
    #ID = ma.getdata(idlist[:-1])

    #new_stressor_out_netcdf(Input.outputDir + '/'+ 'discharge_outlet_human_'  + Input.name_timeperiod, s_aggregated[:-1,:], ID, times, 'discharge flow', 'month', 'm3/s') 
    
    return out

def attribute_label(label,pointer_array):

    n_spatial_unit = pointer_array.shape[0]
    out = np.full((n_spatial_unit), 1e20)
    
    for k in range(n_spatial_unit):
        unique = np.unique(label[pointer_array[k][0],pointer_array[k][1]], return_counts=1)
        #unique = np.unique(label[pointer_array[k,0][0],pointer_array[k,0][1]], return_counts=1)
        unique_count = unique[1][unique[0].mask == 0]
        if len(unique_count)!=0:
            position = np.argmax(unique[1][unique[0].mask == 0])#exclude from the count the masked values in the selection of label value
            value_label = unique[0][position]
            out[k] = value_label
        #out[k]=np.sum(s*a)/np.sum(a)
        
        
    out = ma.masked_values(out, 1e20)
    
    #ID = ma.getdata(idlist[:-1])

    #new_stressor_out_netcdf(Input.outputDir + '/'+ 'discharge_outlet_human_'  + Input.name_timeperiod, s_aggregated[:-1,:], ID, times, 'discharge flow', 'month', 'm3/s') 
    
    return out

def aggregate_average(stressor,pointer_array, area):
    nlat, nlon = stressor.shape
    n_spatial_unit = pointer_array.shape[0]
    out = np.full((n_spatial_unit), 1e20)
    
    for k in range(n_spatial_unit):
        s = stressor[pointer_array[k][0],pointer_array[k][1]]
        a = area[pointer_array[k][0], pointer_array[k][1]]
        #s = stressor[pointer_array[k,0][0],pointer_array[k,0][1]]
        #a = area[pointer_array[k,0][0], pointer_array[k,0][1]]
        out[k] = np.average(s, weights = a)
        #out[k]=np.sum(s*a)/np.sum(a)
        
        
    out = ma.masked_where(isnan(out), out)
    
    #ID = ma.getdata(idlist[:-1])

    #new_stressor_out_netcdf(Input.outputDir + '/'+ 'discharge_outlet_human_'  + Input.name_timeperiod, s_aggregated[:-1,:], ID, times, 'discharge flow', 'month', 'm3/s') 
    
    return out


def aggregate_mean(stressor,pointer_array):
    nlat, nlon = stressor.shape
    n_spatial_unit = pointer_array.shape[0]
    out = np.full((n_spatial_unit), 1e20)
    
    for k in range(n_spatial_unit):
        s = stressor[pointer_array[k][0],pointer_array[k][1]]
        #a = area[pointer_array[k][0], pointer_array[k][1]]
        #s = stressor[pointer_array[k,0][0],pointer_array[k,0][1]]
        #a = area[pointer_array[k,0][0], pointer_array[k,0][1]]
        out[k] = np.mean(s)
        #out[k]=np.sum(s*a)/np.sum(a)
        
        
    out = ma.masked_where(isnan(out), out)
    
    #ID = ma.getdata(idlist[:-1])

    #new_stressor_out_netcdf(Input.outputDir + '/'+ 'discharge_outlet_human_'  + Input.name_timeperiod, s_aggregated[:-1,:], ID, times, 'discharge flow', 'month', 'm3/s') 
    
    return out

def plot_percentiles(s):
    l=[]
    m = np.arange(0,100,1)
    for i in m:
       l.append(stats.scoreatpercentile(s,i)) 
    plt.plot(m,l)
    print(min(s))
    print(l[50])
    print(max(s))
    return l


