# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 13:50:26 2023

@author: Felix Enyimah Toffah
"""
import numpy as np
import netCDF4 as nc


def ncfileWrite(outname, rr, lat1, lon1, datestart, dateend):
    dateend=dateend.toordinal()+1
    datelist = np.arange(datestart.toordinal(), dateend)
    sz = rr.shape
    ncfile = nc.Dataset(outname, mode='w', format='NETCDF4', mmap=False)
    ncfile.createDimension('time', sz[0])  # time axis
    ncfile.createDimension('x', sz[1])  # longitude axis
    lat = ncfile.createVariable('lat', np.float64, ('x'), zlib=True)
    lat.units = 'degrees_north'
    lat.long_name = 'latitude'
    lat[:] = lat1
    lon = ncfile.createVariable('lon', np.float64, ('x'), zlib=True)
    lon.units = 'degrees_east'
    lon.long_name = 'longitude'
    lon[:] = lon1
    time = ncfile.createVariable('time', np.float32, ('time'), zlib=True)
    time.units = 'Days since 1-Jan-0001'
    time.long_name = 'Rainfall event extent'
    time[:] = datelist
    Rainfall = ncfile.createVariable('Rainfall', np.float64, ('time', 'x'))  # note: unlimited dimension is leftmost
    Rainfall.units = 'mm'  # was in kg/m^2/s  *3600-> mm/h
    Rainfall.standard_name = 'Daily rainfall'  # this is a CF standard name
    Rainfall[:] = rr
    ncfile.close()
