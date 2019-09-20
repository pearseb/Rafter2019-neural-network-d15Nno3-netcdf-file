
The python code here takes the .txt file and converts the 1D ASCII data into a 3D gridded netcdf file.

For this to work, you must download the woa13_all_n00_01.nc file from https://www.nodc.noaa.gov/cgi-bin/OC5/woa13/woa13oxnu.pl

Also, some dependencies in the python (version 3.7) environment are required:
	- numpy
	- netCDF4


FINAL OUTPUT  =  Rafter2019_ann_d15n_no3_gridded.nc

	(Rafter2019 is for the data source)
	(ann is for Artificial Neural Network)
	(d15n_no3 is for the variable type)
	(gridded represents the data format, no longer ASCII)

	Variables 	=	d15n_no3
			=	d15n_stdev


Free licence for use. Please contact me at pearse.buchanan@liverpool.ac.uk if you have any questions.


