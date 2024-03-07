# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 13:49:28 2023

@author: Felix Enyimah Toffah
"""

from ftplib import FTP
import os
import pandas as pd

# from datetime import date

def Download_h64(UserName, PassWord, datestart, dateend, storedir, product_category):
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


def Download_h60(UserName, PassWord, datestart, dateend, storedir, product_category):
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


def Download_h26(UserName, PassWord, datestart, dateend, storedir, product_category):
    datelist = pd.date_range(f'{datestart}', f'{dateend}', freq='d')
    
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
