Certainly! Here's a combined version of the test cases for the three methods:

```python
import unittest
import os
import matplotlib.pyplot as plt
import numpy as np
import datetime

from hsaf_data_access import HSAFDataAccess


class TestHSAFDataAccess(unittest.TestCase):
    def test_extract_gz_files(self):
        # Create a temporary directory for testing
        temp_dir = 'test_dir'
        os.makedirs(temp_dir, exist_ok=True)

        # Create a dummy gzipped file
        with open(os.path.join(temp_dir, 'dummy.gz'), 'w') as f:
            f.write('Dummy content')

        # Call the method to extract the gzipped file
        success = HSAFDataAccess.extract_gz_files(temp_dir)

        # Check if the extraction was successful
        self.assertTrue(success)
        # Check if the original gzipped file is deleted
        self.assertFalse(os.path.exists(os.path.join(temp_dir, 'dummy.gz')))

        # Clean up the temporary directory
        os.rmdir(temp_dir)

    def test_add_border(self):
        # Create a dummy plot
        fig, ax = plt.subplots()
        # Call the method to add borders
        ax = HSAFDataAccess.add_border(ax)
        # Check if the axes object is modified
        self.assertIsNotNone(ax)

    def test_create_netCDF_from_data(self):
        # Generate sample data
        param_name = 'temperature'
        param_unit = 'Celsius'
        outname = 'test.nc'
        rr = np.random.rand(5, 10)  # Sample random data array
        lat1 = 0.0
        lon1 = 0.0
        datestart = datetime.datetime(2023, 1, 1)
        dateend = datetime.datetime(2023, 1, 5)

        # Call the method to create netCDF file
        result = HSAFDataAccess.create_netCDF_from_data(param_name, param_unit, outname, rr, lat1, lon1, datestart,
                                                        dateend)

        # Check if the netCDF file is created successfully
        self.assertEqual(result, "netCDF created successfully")


if __name__ == '__main__':
    unittest.main()
```

