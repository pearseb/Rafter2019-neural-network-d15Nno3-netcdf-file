# -*- coding: utf-8 -*-
"""
Created on Thursday Sept 20 16:00:00 2019

@author: pearseb
"""

#%% imports

from __future__ import unicode_literals

import os
import numpy as np
import netCDF4 as nc


#%% get Rafter's dataset

os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\data\\nitrogen isotopes\\convert_Rafter2019_netcdf')
d15N = np.genfromtxt('Rafter2019_inversemodel_d15N_no3.txt',skip_header=1, usecols=(0,1,2,3,4))

# change the longitudes to positives
d15N[:,1][d15N[:,1]<0] += 360.0

# check dimensions
print(np.shape(d15N))

# check for nans in important columns
print(np.any(np.isnan(d15N[:,0])))
print(np.any(np.isnan(d15N[:,1])))
print(np.any(np.isnan(d15N[:,2])))
print(np.any(np.isnan(d15N[:,3])))
print(np.any(np.isnan(d15N[:,4])))

# check for how many values are -999
print(np.any(d15N[:,0]==-999))
print(np.any(d15N[:,1]==-999))
print(np.any(d15N[:,2]==-999))
print(np.any(d15N[:,3]==-999))
print(np.any(d15N[:,4]==-999))
print(len(d15N[d15N[:,4]==-999]))

print(len(d15N[:,0]), "data points")

print("unweighted mean delta15NO3 = ", np.mean(d15N[:,3]), " plus/minus ", np.std(d15N[:,3]) ) 
print("min delta15NO3 = ", np.min(d15N[:,3]) )
print("max delta15NO3 = ", np.max(d15N[:,3]) )

print(np.min(d15N[:,2]), np.max(d15N[:,2]))
print(np.unique(d15N[:,2]))


#%% load WOA grid

data = nc.Dataset('../../nutrients/woa13_all_n00_01.nc', 'r')
lon = data.variables['lon'][...]
lat = data.variables['lat'][...]

lon_bnds = data.variables['lon_bnds'][...]
lat_bnds = data.variables['lat_bnds'][...]

# fix longitudes
lon[lon<0.0] += 360.0
lon_bnds[lon_bnds<0.0] += 360.0

dep = np.array(np.unique(d15N[:,2]))
dep_bnds = np.zeros((len(dep),2))
dep_dif = dep[1::]-dep[0:32]
dep_bnds[0,0] = 0.0; dep_bnds[0,1] = 5.0
dep_bnds[1::,0] = dep[1::] - dep_dif*0.5 
dep_bnds[1::,1] = dep[1::] + dep_dif*0.5 

for i,dep1 in enumerate(dep_bnds[0:32,1]):
    if dep1 != dep_bnds[i+1,0]:
        ave = (dep_bnds[i+1,0] + dep1)*0.5
        dep_bnds[i+1,0] = ave
        dep_bnds[i,1] = ave

print(dep_bnds)


#%% average data at each grid point of the model


d15n_obs_grid = np.zeros((33,180,360))
std_obs_grid = np.zeros((33,180,360))

# for each grid cell
for i,lon1 in enumerate(lon):
    print("index %i, longitude %.1f"%(i, lon1))
    for j,lat1 in enumerate(lat):
        for k,dep1 in enumerate(dep):
            
            d15N_a = d15N[d15N[:,0]==lat1]
            d15N_b = d15N_a[d15N_a[:,1]==lon1]
            d15N_c = d15N_b[d15N_b[:,2]==dep1]
            
            if len(d15N_c) == 0:
                d15n_obs_grid[k,j,i] = np.nan
                std_obs_grid[k,j,i] = np.nan
            elif len(d15N_c) == 1:
                d15n_obs_grid[k,j,i] = d15N_c[0,3]
                std_obs_grid[k,j,i] = d15N_c[0,4]
            else:
                print("More than one value found!")
                exit

np.savez('Rafter_ann_d15N_no3_gridded.npz', d15n_obs_grid=d15n_obs_grid, std_obs_grid=std_obs_grid)


#%% save as npz file (gridded)

data = np.load('Rafter_ann_d15N_no3_gridded.npz')
data.files
d15n_obs_grid = data['d15n_obs_grid']
std_obs_grid = data['std_obs_grid']


#%% change longitudes back to being negative (-180 --> 0)

lon[lon>180.0] -= 360.0

#%% save as netcdf file

## activate line below if overwriting file
os.remove('Rafter2019_ann_d15N_no3_gridded.nc')
data = nc.Dataset('Rafter2019_ann_d15N_no3_gridded.nc', 'w', format='NETCDF4_CLASSIC')

lond = data.createDimension('lon', 360)
latd = data.createDimension('lat', 180)
depd = data.createDimension('dep', 33)
bndd = data.createDimension('nbounds', 2)

lonv = data.createVariable('lon', np.float64, ('lon',))
latv = data.createVariable('lat', np.float64, ('lat',))
depv = data.createVariable('dep', np.float64, ('dep',))
lonbndv = data.createVariable('lon_bnds', np.float64, ('lon','nbounds'))
latbndv = data.createVariable('lat_bnds', np.float64, ('lat','nbounds'))
depbndv = data.createVariable('dep_bnds', np.float64, ('dep','nbounds'))
d15nv = data.createVariable('d15N_no3', np.float64, ('dep','lat','lon'))
stdv = data.createVariable('d15N_stdev', np.float64, ('dep','lat','lon'))

data.description = 'Global d15N of NO3 data produced by Rafter et al. (2019) Biogeosciences regridded onto World Ocean Atlas coordinates'
data.history = "Created by Pearse J. Buchanan on 20th September 2019"
data.source = "Rafter et al. (2019) Biogeosciences dataset - https://hdl.handle.net/1912/24277"

lonv.units = "degrees_east"
latv.units = "degrees_north"
depv.units = "m"
d15nv.units = "per mil"
stdv.units = "per mil"

lonv.standard_name = "longitude"
latv.standard_name = "latitude"
depv.standard_name = "depth"

lonv.long_name = "longitude"
latv.long_name = "latitude"
depv.long_name = "depth"
lonbndv.long_name = "longitude bounds"
latbndv.long_name = "latitude bounds"
depbndv.long_name = "depth bounds"
d15nv.long_name = "d15N of Nitrate (NO3)"
stdv.long_name = "standard deviations of d15N of NO3 data"

lonv.bounds = "lon_bnds"
latv.bounds = "lat_bnds"
depv.bounds = "dep_bnds"

lonv.axis = "X"
latv.axis = "Y"
depv.axis = "Z"

depv.positive = "down"

d15nv.coordinates = "dep lat lon"
stdv.coordinates = "dep lat lon"


lonv[:] = lon
latv[:] = lat
depv[:] = dep
lonbndv[:,:] = lon_bnds
latbndv[:,:] = lat_bnds
depbndv[:,:] = dep_bnds
d15nv[:,:,:] = d15n_obs_grid
stdv[:,:,:] = std_obs_grid


data.close()

