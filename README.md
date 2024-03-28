# H SAF Data Access

This toolkit provides a collection of utilities for accessing and analyzing data products from the Hydrology SAF (H SAF), with a focus on datasets related to snow, precipitation, and soil moisture clusters. It includes functions for data retrieval, processing, visualization, and analysis, making it easier to work with H SAF data in Python environments.


## Usage

### Data Access

To access data products, import the `HSAFDataAccess` class from the `hsaf_data_access` module:

```python
from hsaf_data_access import HSAFDataAccess as data_access
```

### Data Processing

You can download, extract, and preprocess data products using the methods provided by the `HSAFDataAccess` class.

### Data Analysis

Utilize template Jupyter notebooks (`template_h60.ipynb` and `template_h64.ipynb`) included in the `notebooks` directory to perform data analysis and visualization tasks.

### Saving Results

Save processed data into NetCDF4 and CSV file formats using the provided functions.

## Prerequisites

Before using the H SAF data access library and the template tools, ensure that you have the following prerequisites installed:

- Python 3.x
- Jupyter Notebook
- Required Python libraries (listed in `requirements.txt`)

## Data Sources

You can access the data products through the H-SAF server by obtaining valid credentials from the [H SAF website](https://hsaf.meteoam.it/).

## Documentation

For detailed usage instructions and examples, refer to the guideline documentation.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

## Acknowledgements

The H SAF Toolkit is developed as part of the Hydrology SAF (H SAF) initiative.
