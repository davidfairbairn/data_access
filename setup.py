from setuptools import setup, find_packages

with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name='hsaf_toolkit',
    version='0.1',
    author='Felix Enyimah Toffah',
    author_email='gmfetoffah@gmail.com',
    description='Toolkit for accessing H64 data products and analyzing them using Jupyter notebooks',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/H-SAF/data_access',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'netCDF4',
        'geopandas',
        'xarray',
        'pandas',
        'matplotlib',
        'cartopy',
        'ipywidgets',
        'ftplib',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    # Include the additional files within the package
    package_data={
        'hsaf_data_access': ['hsaf_data_access.py'],
        'notebooks': ['template_h60.ipynb', 'template_h64.ipynb'],
    },
)
