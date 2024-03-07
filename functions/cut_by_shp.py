# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 20:39:57 2024

@author: fetoffah
"""

import geopandas as gpd
import xarray as xr

def cut_by_shp(shapefile_path, netcdf_path):
    # read in file path for shapefile
    fp_shp = shapefile_path
    # read in netcdf file path
    ncs = netcdf_path
    
    if not fp_shp.lower().endswith('.shp'):
        print(shapefile_path + 'not a shapefile')
        return
    if not ncs.lower().endswith('.nc'):
        print(netcdf_path + 'not a netCDF4 file')
        return
    
    # Read in NETCDF as a pandas dataframe
    ds = xr.open_dataset(ncs, decode_times=False)
    edgar = ds.to_dataframe()
    
    # Reset the index in the data frame
    edgar = edgar.reset_index()
    
    # Read shapefile
    shp = gpd.read_file(fp_shp)
    
    # use geopandas points_from_xy() to transform Longitude and Latitude into a list of shapely.Point objects and set it as a geometry while creating the GeoDataFrame
    edgar_gdf = gpd.GeoDataFrame(edgar, geometry=gpd.points_from_xy(edgar.lon, edgar.lat))
    
    
    # set coordinates equal to each other
    edgar_gdf.crs = shp.crs
    
    # # Clip points, lines, or polygon geometries to the mask extent.
    mask = gpd.clip(edgar_gdf, shp)
    mask.to_csv(ncs.split('.')[0] + '.csv')
    print('Done!!!')
    return
