# -*- coding: utf-8 -*-


import os
import gzip
import shutil

import numpy as np
import netCDF4 as nc
import geopandas as gpd
import xarray as xr
import pandas as pd
from ftplib import FTP
from matplotlib.path import Path
import cartopy.feature as cfeature


class HSAFDataAccess:
    """ A class containing methods for accessing data from the H-SAF ftp server """
    

    def __init__(self):
        pass

    @staticmethod        
    def extract_gz_files(directory):
        """ 
        Extracts all gzipped files in the specified directory and deletes the original gzipped files.

        Parameters:
            directory (str): The path of the directory containing the gzipped files.

        Returns:
            bool: True if extraction was successful, False otherwise.
        """
        if not os.path.isdir(directory):
            raise ValueError("Specified path is not a directory")
            
            
        # Get list of files in the directory
        file_list = os.listdir(directory)
        file_list.sort()

        # Check if any files are available
        if not file_list:
            print("No data files found in the specified directory.")
            return
        
        

        try:
            extracted = False
             
            # Check if files are of expected type (gz)
            gz_files = [file for file in file_list if file.endswith('.gz')]
            if not gz_files:
                print("No .gz files found in the specified directory.")
                return
            # Extract and clean the data
            for file in gz_files:
                gz_name = os.path.join(directory, file)
                file_name = os.path.splitext(gz_name)[0]
                try:
                    with gzip.open(gz_name, 'rb') as f_in:
                        with open(file_name, "wb") as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(gz_name)
                    extracted = True

                except Exception as e:
                    print(f"Error extracting data from {file}: {e}")
                    continue

        except Exception as e:
            raise RuntimeError(f"Error while extracting gzipped files: {e}")

        return extracted
    




    @staticmethod
    def add_border(ax):
        """
        Add ocean and land boundaries to the given axes object.

        Parameters:
            ax (matplotlib.axes.Axes): The axes object to which boundaries are to be added.

        Returns:
            matplotlib.axes.Axes: The modified axes object.
        """

        # Define boundary types and their corresponding face colors
        boundaries = {'ocean': (173/255, 216/255, 230/255), 'land': 'none'}  # Light blue color for ocean

        try:
            for boundary_type, face_color in boundaries.items():
                feature = cfeature.NaturalEarthFeature(category='physical', name=boundary_type, scale='10m',
                                                        edgecolor='none', facecolor=face_color)
                ax.add_feature(feature)

            # Add coastlines
            ax.coastlines(resolution='10m', zorder=-1)

            return ax

        except Exception as e:
            print(f"An error occurred while adding boundaries: {e}")
            return None




    
    @staticmethod
    def create_netCDF_from_data(param_name, param_unit, outname, rr, lat1, lon1, datestart, dateend):
        """
        Create a NetCDF file from provided data.

        Parameters:
            param_name (str): The name for the main variable being studied.
            param_unit (str): The unit of the main variable.
            outname (str): The name for the NetCDF4 file to be created.
            rr (numpy.ndarray): The main variable data array.
            lat1 (float): The latitude coordinate.
            lon1 (float): The longitude coordinate.
            datestart (datetime.datetime): Start date of the time dimension.
            dateend (datetime.datetime): End date of the time dimension.
        """
        
        try:
            sz = rr.shape
            date_list = pd.date_range(datestart, dateend)

            ncfile = nc.Dataset(outname, mode='w', format='NETCDF4', mmap=False)
            ncfile.createDimension('time', len(date_list))  # time axis
            ncfile.createDimension('x', sz[1])  # longitude axis

            lat = ncfile.createVariable('lat', np.float32, ('x'), zlib=True)
            lat.units = 'degrees_north'
            lat.long_name = 'Latitude'
            lat[:] = lat1

            lon = ncfile.createVariable('lon', np.float32, ('x'), zlib=True)
            lon.units = 'degrees_east'
            lon.long_name = 'Longitude'
            lon[:] = lon1

            time = ncfile.createVariable('time', np.float64, ('time',), zlib=True)
            time.units = 'days since 1970-01-01'
            time.long_name = 'Time'
            time.calendar = 'standard'
            time[:] = nc.date2num(date_list.to_pydatetime(), units=time.units, calendar=time.calendar)

            main_param = ncfile.createVariable(param_name, np.float64, ('time', 'x'), zlib=True)
            main_param.units = param_unit
            main_param.standard_name = param_name 
            main_param[:] = rr

            ncfile.close()
            return "netCDF created successfully"
        except Exception as e:
            print(f'An error occurred: {str(e)}')
            return 'could not create netCDF from data'


    @staticmethod
    def points_in_polygon(query_x, query_y, vertex_x, vertex_y):
        """
        Determine whether points are inside a polygon.

        Parameters:
            query_x (numpy.ndarray): X coordinates of the query points.
            query_y (numpy.ndarray): Y coordinates of the query points.
            vertex_x (numpy.ndarray): X coordinates of the polygon vertices.
            vertex_y (numpy.ndarray): Y coordinates of the polygon vertices.

        Returns:
            numpy.ndarray: Boolean array indicating whether each query point is inside the polygon.
        """
        # Ensure inputs are numpy arrays
        query_x = np.array(query_x)
        query_y = np.array(query_y)
        vertex_x = np.array(vertex_x)
        vertex_y = np.array(vertex_y)

        shape = query_x.shape
        query_x = query_x.reshape(-1)
        query_y = query_y.reshape(-1)
        vertex_x = vertex_x.reshape(-1)
        vertex_y = vertex_y.reshape(-1)

        # Convert vertices and query points to list of tuples
        query_points = [(query_x[i], query_y[i]) for i in range(query_x.shape[0])]
        polygon_vertices = [(vertex_x[i], vertex_y[i]) for i in range(vertex_x.shape[0])]

        # Create a Path object representing the polygon
        polygon_path = Path(polygon_vertices)

        # Determine whether each query point is inside the polygon
        points_inside_polygon = polygon_path.contains_points(query_points).reshape(shape)

        return points_inside_polygon


    @staticmethod
    def points_in_shapefile(lon, lat, shapefile_path):
        """
        Checks if points are within the boundaries defined by a shapefile.

        Args:
        - lon (array-like): Array-like object containing longitudes of the points.
        - lat (array-like): Array-like object containing latitudes of the points.
        - shapefile_path (str): Path to the shapefile.

        Returns:
        - list: List containing boolean values indicating whether each point is inside the shapefile boundary.
        """

        try:
            shapefile = gpd.read_file(shapefile_path)
            points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(lon, lat))
            points_in_shapefile = gpd.sjoin(points, shapefile, how="inner", op="within")
            points_indices_inside_shapefile = points_in_shapefile.index.tolist()
            
            # Construct a boolean mask indicating whether each point is inside the shapefile
            points_in_shapefile_mask = [True if i in points_indices_inside_shapefile else False for i in range(len(lon))]
            
            return points_in_shapefile_mask
        
        except Exception as e:
            print(f'Error checking points in shapefile: {str(e)}')


    @staticmethod
    def filter_data_by_bounding_box(lon, lat, P_h60, bounding_box):
        """
        Filter data points by a bounding box.
    
        This method filters data points based on whether they fall within the specified bounding box.
        It uses the `points_in_polygon` method to determine which points are inside the bounding box.
    
        Parameters:
            lon (array-like): Array-like object containing longitudes of the data points.
            lat (array-like): Array-like object containing latitudes of the data points.
            P_h60 (array-like): Array-like object containing the data values associated with each point.
            bounding_box (tuple): A tuple containing the coordinates (x, y) defining the bounding box.
    
        Returns:
            tuple: A tuple containing three arrays: filtered_lon, filtered_lat, and filtered_P_h60.
                   These arrays contain the filtered data points that fall within the bounding box.
        """
        try:
            x, y = bounding_box
            lon_h60, lat_h60 = np.meshgrid(lat, lon, sparse=True)
            lat_h60 = np.ravel(lat_h60)
            lon_h60 = np.ravel(lon_h60)
            P_h60 = np.ravel(P_h60)
            IN = HSAFDataAccess.points_in_polygon(lon_h60, lat_h60, x, y)
            filtered_lon = lon_h60[IN]
            filtered_lat = lat_h60[IN]
            filtered_P_h60 = P_h60[IN]
            return filtered_lon, filtered_lat, filtered_P_h60
        except Exception as e:
            print(f'Error filtering data by bounding box: {str(e)}')
                                            

    @staticmethod
    def filter_data_by_shapefile(lon, lat, P_h60, shapefile_path):
        """
        Filter data points by a shapefile boundary.
    
        This method filters data points based on whether they fall within the boundaries
        defined by a shapefile. It uses the GeoPandas library to read the shapefile and
        perform spatial join operations to determine which points are inside the shapefile.
    
        Parameters:
            lon (array-like): Array-like object containing longitudes of the data points.
            lat (array-like): Array-like object containing latitudes of the data points.
            P_h60 (array-like): Array-like object containing the data values associated with each point.
            shapefile_path (str): Path to the shapefile defining the boundary.
    
        Returns:
            tuple: A tuple containing three arrays: filtered_lon, filtered_lat, and filtered_P_h60.
                   These arrays contain the filtered data points that fall within the shapefile boundary.
        """
        try:
            lon_h60, lat_h60 = np.meshgrid(lat, lon, sparse=True)
            lat_h60 = np.ravel(lat_h60)
            lon_h60 = np.ravel(lon_h60)
            P_h60 = np.ravel(P_h60)
            IN = HSAFDataAccess.points_in_shapefile(lon_h60, lat_h60, shapefile_path)
            filtered_lon = lon_h60[IN]
            filtered_lat = lat_h60[IN]
            filtered_P_h60 = P_h60[IN]
            return filtered_lon, filtered_lat, filtered_P_h60
        except Exception as e:
            print(f'Error filtering data by shapefile: {str(e)}')



    @staticmethod
    def cut_netcdf_by_shapefile(shapefile_path, netcdf_path):
        """
        Cut the extent of a NetCDF file to the region specified by the boundary of the given shapefile 
        and export the result as a CSV file.

        Parameters:
            shapefile_path (str): Path to the shapefile defining the region of interest.
            netcdf_path (str): Path to the NetCDF file to be cut.

        Returns:
            bool: True if the operation succeeds, False otherwise.
        """
        # Input validation
        if not shapefile_path.lower().endswith('.shp'):
            raise ValueError(shapefile_path + ' is not a shapefile.')
        if not netcdf_path.lower().endswith('.nc'):
            raise ValueError(netcdf_path + ' is not a NetCDF4 file.')

        try:
            # Open NetCDF dataset and convert to DataFrame
            ds = xr.open_dataset(netcdf_path, decode_times=True)
            df = ds.to_dataframe()

            # Read shapefile
            shp = gpd.read_file(shapefile_path)

            # Convert DataFrame to GeoDataFrame
            gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat))
            gdf.crs = shp.crs

            # Clip GeoDataFrame by shapefile boundary
            clipped_gdf = gpd.clip(gdf, shp)

            # Export clipped GeoDataFrame to CSV
            output_filename = os.path.splitext(netcdf_path)[0] + '_clipped.csv'
            clipped_gdf.to_csv(output_filename)

            return True
        except Exception as e:
            print(f'An error occurred: {str(e)}')
            return False

    @staticmethod
    def download_h64(UserName, PassWord, datestart, dateend, storedir, product_category):
        """
        Download H64 data from the FTP server.

        Parameters:
            UserName (str): FTP username.
            PassWord (str): FTP password.
            datestart (str): Start date in YYYY-MM-DD format.
            dateend (str): End date in YYYY-MM-DD format.
            storedir (str): Directory to store downloaded files.
            product_category (str, optional): Product category. Defaults to 'h64'.
        """

        # Validate product category
        if product_category != 'h64':
            print('Error: Invalid product category. Use product category h64.')
            return

        # Validate date range
        try:
            pd.to_datetime(datestart)
            pd.to_datetime(dateend)
        except ValueError:
            print('Error: Invalid date format. Use YYYY-MM-DD.')
            return

        if datestart > dateend:
            print('Error: Invalid date range. Ensure datestart is less than or equal to dateend.')
            return

        # Validate directory path
        if not os.path.isdir(storedir):
            print('Error: Invalid directory path. Ensure storedir is a valid directory.')
            return

        # Establish FTP connection
        try:
            with FTP('ftphsaf.meteoam.it') as ftp:
                ftp.login(user=UserName, passwd=PassWord)

                # Change directory to product category
                spt = f'{product_category}/{product_category}_cur_mon_data/'
                ftp.cwd(spt)

                print('Info: Preparing data download...')

                # Generate list of dates
                datelist = pd.date_range(datestart, dateend, freq='d')

                # Download files for each date
                for date in datelist:
                    filename = f'{product_category}_{date.strftime("%Y%m%d")}_0000_24_hea.nc.gz'

                    # Check if filename is a substring of any file on the server
                    files = ftp.nlst()
                    if any(filename in file for file in files):
                        local_filename = os.path.join(storedir, filename)
                        with open(local_filename, 'wb') as loc_file:
                            ftp.retrbinary(f'RETR {filename}', loc_file.write)
                        print(f'Info: {filename} downloaded')
                    else:
                        print(f'Warning: No data available for {date.strftime("%d/%m/%Y")}')

                print('Infor: Download completed successfully')

        except Exception as e:
            print(f'Error: An error occurred: {str(e)}')

    @staticmethod
    def download_h60(UserName, PassWord, datestart, dateend, storedir, product_category):
        """
        Download H60 data from the FTP server.

        Parameters:
            UserName (str): FTP username.
            PassWord (str): FTP password.
            datestart (str): Start date in YYYY-MM-DD format.
            dateend (str): End date in YYYY-MM-DD format.
            storedir (str): Directory to store downloaded files.
            product_category (str, optional): Product category. Defaults to 'h60'.
        """
    
        # Validate product category
        if product_category != 'h60':
            print('Error: Invalid product category. Use product category h60.')
            return

        # Validate date format
        try:
            pd.to_datetime(datestart)
            pd.to_datetime(dateend)
        except ValueError:
            print('Error: Invalid date format. Use YYYY-MM-DD.')
            return

        # Validate date range
        if datestart > dateend:
            print('Error: Invalid date range. Ensure datestart is less than or equal to dateend.')
            return

        # Validate directory path
        if not os.path.isdir(storedir):
            print('Error: Invalid directory path. Ensure storedir is a valid directory.')
            return

        # Generate list of dates
        datelist = pd.date_range(datestart, dateend, freq='d')

        # Establish FTP connection
        try:
            with FTP('ftphsaf.meteoam.it') as ftp:
                ftp.login(user=UserName, passwd=PassWord)

                # Change directory to product category
                spt = f'{product_category}/{product_category}_cur_mon_data/'
                ftp.cwd(spt)

                print('Preparing data download...')

                # Download files for each date
                files = ftp.nlst()

                for date in datelist:
                    filename = f'{product_category}_{date.strftime("%Y%m%d")}'
                    for file in files:
                        if filename in file:
                            local_filename = os.path.join(storedir, file)
                            with open(local_filename, 'wb') as loc_file:
                                ftp.retrbinary(f'RETR {file}', loc_file.write)
                            print(f'Info: {file} downloaded')
                    else:
                        print(f'Warning: No data available for {date.strftime("%d/%m/%Y")}')


                print('Info: Download completed successfully')

        except Exception as e:
            print(f'An error occurred: {str(e)}')

            
    @staticmethod
    def download_h26(UserName, PassWord, datestart, dateend, storedir, product_category):
        """
        Download H26 data from the FTP server.

        Parameters:
            UserName (str): FTP username.
            PassWord (str): FTP password.
            datestart (str): Start date in YYYY-MM-DD format.
            dateend (str): End date in YYYY-MM-DD format.
            storedir (str): Directory to store downloaded files.
            product_category (str, optional): Product category. Defaults to 'h26'.
        """
        # Validate product category
        if product_category != 'h26':
            print('Error: Invalid product category. Use product category h26.')
            return

        # Validate date format
        try:
            pd.to_datetime(datestart)
            pd.to_datetime(dateend)
        except ValueError:
            print('Error: Invalid date format. Use YYYY-MM-DD.')
            return

        # Validate date range
        if datestart > dateend:
            print('Error: Invalid date range. Ensure datestart is less than or equal to dateend.')
            return

        # Validate directory path
        if not os.path.isdir(storedir):
            print('Error: Invalid directory path. Ensure storedir is a valid directory.')
            return

        try:
            # Establish FTP connection
            with FTP('ftphsaf.meteoam.it') as ftp:
                ftp.login(user=UserName, passwd=PassWord)

                # Change directory to product category
                spt = f'{product_category}/{product_category}_cur_mon_nc/'
                ftp.cwd(spt)

                print('Info: Preparing data download...')

                # Generate list of dates
                datelist = pd.date_range(datestart, dateend, freq='d')

                # Download files for each date
                files = ftp.nlst()

                for date in datelist:
                    filename = f'{product_category}_{date.strftime("%Y%m%d")}00_R01.nc'
                    file_found = False
                    for file in files:
                        if filename in file:
                            file_found = True
                            local_filename = os.path.join(storedir, file)
                            with open(local_filename, 'wb') as loc_file:
                                ftp.retrbinary(f'RETR {file}', loc_file.write)
                            print(f'Info: {file} downloaded')
                            break  # Exit the inner loop once the file is found
                    if not file_found:
                        print(f'Warning: No data available for {date.strftime("%d/%m/%Y")}')

                print('Info: Download completed successfully')

        except Exception as e:
            print(f'Error: An error occurred: {str(e)}')

        
    @staticmethod            
    def get_lat_lon(username, psw, filename='lat_lon_0.nc.gz', local_filename='latlon.nc.gz'):
        """
        Download latitude and longitude data from the FTP server.

        Parameters:
            username (str): FTP username.
            psw (str): FTP password.
            filename (str, optional): Name of the file to download. Defaults to 'lat_lon_0.nc.gz'.
            local_filename (str, optional): Name to save the downloaded file locally. Defaults to 'latlon.nc.gz'.

        Returns:
            str: The path to the downloaded file.
        """
        try:
            # Connect to FTP server
            with FTP('ftphsaf.meteoam.it') as ftp:
                ftp.login(user=username, passwd=psw)
                ftp.cwd('utilities/matlab_code')

                # Download file
                with open(local_filename, 'wb') as loc_file:
                    if filename in ftp.nlst():
                        ftp.retrbinary(f'RETR {filename}', loc_file.write)
                        print(f'Info: {filename} downloaded')
                    else:
                        print('Warning: No positional data available')
                        return None

            # Extract gz file
            with gzip.open(local_filename, 'rb') as f_in:
                with open('lat_lon.nc', 'wb') as f_out:
                    f_out.write(f_in.read())
                
    
            print('Info: Download completed successfully')
            return os.path.join('./', 'lat_lon.nc')

        except Exception as e:
            print(f'Error: An error occurred: {str(e)}')
            return None
    
    @staticmethod
    def create_folders():
        """
        Creates necessary folders if they don't exist.
    
        The method checks if the 'output' and 'data' folders exist. If they do not exist,
        it creates them. If they exist, it deletes all files inside the 'data' folder.
    
        Parameters:
            None
    
        Returns:
            None
        """
        if not os.path.exists('output'):
            os.makedirs('output')
        if not os.path.exists('data'):
            os.makedirs('data')
        else:
            for f in os.listdir('data'):
                os.remove(os.path.join('data', f))


    @staticmethod
    def extract_and_clean_data(store_dir):
        """
        Extracts all gzipped files in the specified directory and removes any empty files after extraction.
    
        Parameters:
            store_dir (str): The directory containing the gzipped files to be extracted and cleaned.
    
        Returns:
            None
        """
        
        try:
            HSAFDataAccess.extract_gz_files(store_dir)
    
            # Remove empty files after extraction
            for file_name in os.listdir(store_dir):
                file_path = os.path.join(store_dir, file_name)
                if os.stat(file_path).st_size == 0:
                    os.remove(file_path)
        except Exception as e:
            print(f'Error extracting and cleaning data: {str(e)}')

# def get_selected_area():
#     """
#     Get the selected area based on the radio button option.

#     Returns:
#         dict or None: A dictionary containing the selected area's information. 
#         If the bounding box option is selected, it returns the coordinates of the bounding box. 
#         If the shapefile option is selected and a shapefile is uploaded, it returns the content of the shapefile. 
#         If neither option is selected or no shapefile is uploaded, it returns None.
#     """
#     if radio_button.value == 'Bounding Box':
#         lower_left_x = bounding_box_widgets.children[1].value
#         lower_left_y = bounding_box_widgets.children[2].value
#         upper_right_x = bounding_box_widgets.children[3].value
#         upper_right_y = bounding_box_widgets.children[4].value
#         return {
#             'area_type': 'Bounding Box',
#             'lower_left_x': lower_left_x,
#             'lower_left_y': lower_left_y,
#             'upper_right_x': upper_right_x,
#             'upper_right_y': upper_right_y
#         }
#     elif radio_button.value == 'Shapefile':
#         if shapefile_widget.value:
#             file_info = list(shapefile_widget.value.values())[0]
#             shapefile_content = file_info['content']
#             return {
#                 'area_type': 'Shapefile',
#                 'shapefile_content': shapefile_content
#             }
#         else:
#             return None

# # Example usage
# selected_area = get_selected_area()
# print(selected_area)




# # Accessing bounding box widget values
# lower_left_x = bounding_box_widgets.children[1].value
# lower_left_y = bounding_box_widgets.children[2].value
# upper_right_x = bounding_box_widgets.children[3].value
# upper_right_y = bounding_box_widgets.children[4].value  
    
# # Accessing shapefile widget value
# uploaded_files = aoi_shp.value
# if uploaded_files:
#     # Iterate over the uploaded files
#     for file_info in uploaded_files:
#         # Get the content of the uploaded file
#         shapefile_content = file_info['content']
#         print(shapefile_content)
#         # Process the shapefile content as needed
        
        
# # Function to create bounding box widgets
# def create_bounding_box_widgets():
#     return widgets.VBox([
#         widgets.Label('Bounding Box Coordinates:'),
#         widgets.FloatText(value=0.0, description='Lower Left X:'),
#         widgets.FloatText(value=0.0, description='Lower Left Y:'),
#         widgets.FloatText(value=0.0, description='Upper Right X:'),
#         widgets.FloatText(value=0.0, description='Upper Right Y:')
#     ])

# # Area of study selection
# area_of_study = widgets.RadioButtons(
#     options=['Bounding Box', 'Shapefile'],
#     description='Area of Study:',
#     disabled=False
# )

# # Bounding box widgets (initially hidden)
# bounding_box_widgets = create_bounding_box_widgets()
# bounding_box_widgets.layout.visibility = 'hidden'

# # Shapefile upload widget (initially hidden)
# aoi_shp = widgets.FileUpload(description='Upload Shapefile')
# aoi_shp.layout.visibility = 'hidden'

# # Container for bounding box and shapefile widgets
# area_widgets_container = widgets.VBox([
#     bounding_box_widgets,
#     aoi_shp
# ])

# # Function to handle area of study selection
# def handle_area_of_study_selection(change):
#     if change.new == 'Bounding Box':
#         bounding_box_widgets.layout.visibility = 'visible'
#         aoi_shp.layout.visibility = 'hidden'
#     elif change.new == 'Shapefile':
#         bounding_box_widgets.layout.visibility = 'hidden'
#         aoi_shp.layout.visibility = 'visible'

# area_of_study.observe(handle_area_of_study_selection, names='value')

# # Display widgets
# display(widgets.HBox([area_of_study, area_widgets_container]))
