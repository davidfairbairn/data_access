# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 13:49:29 2023

@author: Felix Enyimah Toffah
"""

import os
import gzip
import shutil

def gz_extract(directory):
    extension = ".gz"
    os.chdir(directory)
    for item in os.listdir(directory):  # loop through items in dir
        if item.endswith(extension):  # check for ".gz" extension
            gz_name = os.path.abspath(item)  # get full path of files
            file_name = (os.path.basename(gz_name)).rsplit('.', 1)[0]  # get file name for file within
            with gzip.open(gz_name, "rb") as f_in:
                with open(file_name, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(gz_name)  # delete zipped file