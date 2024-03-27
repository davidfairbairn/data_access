#!/bin/bash

# Create the virtual environment
python -m venv hsaf

# Activate the virtual environment
source hsaf/bin/activate

# Install required libraries
pip install numpy netCDF4 geopandas xarray pandas matplotlib cartopy ipywidgets

# Run Jupyter Notebook server
jupyter notebook

# Deactivate the virtual environment when done
deactivate
