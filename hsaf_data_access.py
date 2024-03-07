# -*- coding: utf-8 -*-
"""

This is a class containing methods for accessing the H-SAF data products on the ftp server.

"""

import os
import gzip
import shutil
from datetime import timedelta


import numpy as np
import netCDF4 as nc
import geopandas as gpd
import xarray as xr
import pandas as pd
import cartopy as cart
from ftplib import FTP
from matplotlib import path


class HSAFDataAccess:
    """ A class containing methods for accessing data from the H-SAF ftp server """
    

    def __init__(self):
        pass
        
    def extract_gz_file(directory):
        """ 
            This method looks for the zipped files in the specified directory, extracts all is contents in and delete the zipped file.
        
            Parameters
                directory: is the path in which to look for the zipped files
        """
        
        
        if not os.path.isdir(directory):
            print('Specified path is not a directory')
        else:
            os.chdir(directory)
            extension = ".gz"

            for item in os.listdir(directory):  # loop through items in dir
                if item.endswith(extension):  # check for ".gz" extension
                    gz_name = os.path.abspath(item)  # get full path of files
                    file_name = (os.path.basename(gz_name)).rsplit('.', 1)[0]  # get file name for file within
                    with gzip.open(gz_name, "rb") as f_in:
                        with open(file_name, "wb") as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(gz_name)  # delete zipped file
                else:
                    print('No zipped file found in the specified directory')
        return
    
    def add_border(ax):
        """
            This method adds ocean and land boundaries to the plot containing the ‘ax’ object
            
            Parameters
            ax: is the axes on which the boundary is to be shown
        
        """
        axi = ax
        param = ['ocean', 'land']
        for i in param:
            if i == 'ocean':
                facecolor='water'
            else:
                facecolor=i
            feature = cart.feature.NaturalEarthFeature(category='physical', name=i, scale='10m',
                                        edgecolor='none', facecolor=cart.feature.COLORS[facecolor], zorder=-1)
            axi.add_feature(feature)
        
        # coastlines
        axi.coastlines(resolution='10m',zorder=-1)
        return axi
    
    def to_netCDF(param_name, param_unit, outname, rr, lat1, lon1, datestart, dateend):
    
        """
            Parameters
                param_name: The name for the main variable being studied
                
                param_unit: The unit of the main variable
                
                outname:  is the name for the netCDF4 file to be created
        
        """
        
        sz = rr.shape
        date_list = []
        while datestart <= dateend:
            date_list.append(datestart)
            datestart += timedelta(days=1)
        
        ncfile = nc.Dataset(outname, mode='w', format='NETCDF4', mmap=False)
        ncfile.createDimension('time', sz[0])  # time axis
        ncfile.createDimension('x', sz[1])  # longitude axis
        
        lat = ncfile.createVariable('lat', np.float32, ('x'), zlib=True)
        lat.units = 'degrees north'
        lat.long_name = 'Latitude'
        lat[:] = lat1
        
        lon = ncfile.createVariable('lon', np.float32, ('x'), zlib=True)
        lon.units = 'degrees east'
        lon.long_name = 'Longitude'
        lon[:] = lon1
        
        time = ncfile.createVariable('time', np.float64, ('time'), zlib=True)
        time.units = 'days since 1970-01-01'
        time.long_name = 'Date'
        time.calendar = 'standard'    
        time[:] = nc.date2num(date_list, units=time.units, calendar=time.calendar)
       
        main_param = ncfile.createVariable(param_name, np.float64, ('time', 'x'), zlib=True)
        main_param.units = param_unit
        main_param.standard_name = param_name 
        main_param[:] = rr
        ncfile.close()
        
        return

    def create_polygon(xq, yq, xv, yv):
        shape = xq.shape
        xq = xq.reshape(-1)
        yq = yq.reshape(-1)
        xv = xv.reshape(-1)
        yv = yv.reshape(-1)
        q = [(xq[i], yq[i]) for i in range(xq.shape[0])]
        p = path.Path([(xv[i], yv[i]) for i in range(xv.shape[0])])
        return p.contains_points(q).reshape(shape)

    def cut_by_shp(shapefile_path, netcdf_path):
        """
           Cut the extent of a netCDF4 file to the region specified by the boundary of the given shapefile and export the result as a CSV file.
           
        """
        
        if not shapefile_path.lower().endswith('.shp'):
            print(shapefile_path + 'not a shapefile')
            return
        if not netcdf_path.lower().endswith('.nc'):
            print(netcdf_path + 'not a netCDF4 file')
            return
        
        ds = xr.open_dataset(netcdf_path, decode_times=True)
        edgar = ds.to_dataframe()
        
        shp = gpd.read_file(shapefile_path)
        
        # use geopandas points_from_xy() to transform Longitude and Latitude into a list of shapely.Point objects and set it as a geometry while creating the GeoDataFrame
        edgar_gdf = gpd.GeoDataFrame(edgar, geometry=gpd.points_from_xy(edgar.lon, edgar.lat))
        edgar_gdf.crs = shp.crs
        
        mask = gpd.clip(edgar_gdf, shp)
        mask.to_csv(netcdf_path.split('.')[0] + '.csv')
        return

    def download_h64(UserName, PassWord, datestart, dateend, storedir, product_category):
        datelist = pd.date_range(f'{datestart}', f'{dateend}', freq='d')
        prod_type = ''
        prod_ext = ''
        spt = ''
        
        with FTP('ftphsaf.meteoam.it', user=f'{UserName}', passwd=f'{PassWord}') as ftp:
            if product_category == 'h64':
                prod_type = product_category
            else:
                print('Use product category h64')
                return
            
            spt = prod_type + '/' + prod_type + '_cur_mon_data/'
                   
            ftp.cwd(spt)
            
            print('Preparing data download ...')
            
            for ii in range(0, len(datelist)):
                filename = prod_type+ '_' + f'{datelist[ii].year}{str(datelist[ii].month).zfill(2)}{str(datelist[ii].day).zfill(2)}' + prod_ext + '_0000_24_hea.nc.gz'
                local_filename = os.path.join(f'{storedir}', filename)                
                with open(local_filename, 'wb') as loc_file:
                    if filename in ftp.nlst():
                        ftp.retrbinary('RETR ' + filename, loc_file.write)
                        print(filename + ' downloaded')
                    else:
                        print('there is no data for ' + '%s/%s/%s' % (datelist[ii].day, datelist[ii].month, datelist[ii].year))
            print('Done!!!')

    def download_h60(UserName, PassWord, datestart, dateend, storedir, product_category):        
        datelist = pd.date_range(f'{datestart}', f'{dateend}', freq='d')
        prod_type = ''
        
        with FTP('ftphsaf.meteoam.it', user=f'{UserName}', passwd=f'{PassWord}') as ftp:
            if product_category == 'h60':
                prod_type = product_category
            else:
                print('Use product category h60')
                return
                   
            spt = prod_type + '/' + prod_type + '_cur_mon_data/'
            ftp.cwd(spt)
            
            print('Preparing data download ...')
            
            for ii in range(0, len(datelist)):
                if len(ftp.nlst()) ==0:
                    print('there is no data for ' + '%s/%s/%s' % (datelist[ii].day, datelist[ii].month, datelist[ii].year))
                else:
                    for file in ftp.nlst():
                        filename = prod_type+ '_' + f'{datelist[ii].year}{str(datelist[ii].month).zfill(2)}{str(datelist[ii].day).zfill(2)}'
                        if filename in file:
                            local_filename = os.path.join(f'{storedir}', file)
                            with open(local_filename, 'wb') as loc_file:
                                ftp.retrbinary('RETR ' + file, loc_file.write)
                                print(file + ' downloaded')
            
            print('Done!!!')
        
    def download_h26(UserName, PassWord, datestart, dateend, storedir, product_category):
        
        datelist = pd.date_range(f'{datestart}', f'{dateend}', freq='d')
        prod_type = ''
        
        with FTP('ftphsaf.meteoam.it', user=f'{UserName}', passwd=f'{PassWord}') as ftp:
            if product_category == 'h26':
                prod_type = product_category
            else:
                print('Use product category h26')
                return
            spt = prod_type + '/' + prod_type + '_cur_mon_nc/'
                       
            ftp.cwd(spt)
            
            print('Preparing data download ...')
            
            for ii in range(0, len(datelist)):
                filename = prod_type+ '_' + f'{datelist[ii].year}{str(datelist[ii].month).zfill(2)}{str(datelist[ii].day).zfill(2)}' + '00_R01.nc'
                local_filename = os.path.join(f'{storedir}', filename)                
                with open(local_filename, 'wb') as loc_file:
                    if filename in ftp.nlst():
                        ftp.retrbinary('RETR ' + filename, loc_file.write)
                        print(filename + ' downloaded')
                    else:
                        print('there is no data for ' + '%s/%s/%s' % (datelist[ii].day, datelist[ii].month, datelist[ii].year))
            print('Done!!!')


