#!/bin/bash

# Create the Conda environment
conda create -n hsaf python=3.8

# Activate the Conda environment
conda activate hsaf

# Install required libraries
conda install numpy netCDF4 geopandas xarray pandas matplotlib cartopy ipywidgets

# Install Jupyter Notebook
conda install jupyter

# Run Jupyter Notebook server
jupyter notebook

# Deactivate the Conda environment when done
conda deactivate