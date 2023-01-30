#!/usr/bin/env python

# # Script screated by Davide Marchegiani ar ACCESS-NRI - davide.marchegiani@anu.edu.au

# # Interpolate an ancillary file onto another grid

import argparse

# Parse arguments
description='''Regrid UM ancillary file onto another ancillary file or netCDF file grid.'''
parser = argparse.ArgumentParser(description=description, allow_abbrev=False)
parser.add_argument('-i', '--input', dest='um_input_file', required=True, type=str,
                    help='UM ancillary input file. (Required)')
parser.add_argument('-g', "--grid", dest='gridfile', required=True, type=str,
                    help='UM ancillary file or netCDF file on the new grid. (Required)')
parser.add_argument('-o', '--output', dest='um_output_file', required=False, type=str,
                    help='UM ancillary output file.')
parser.add_argument('--lat', '--latitude', dest='ncfile_latitude_name', required=False, type=str,
                    help='Name of the latitude dimension in the netCDF file.')
parser.add_argument('--lon', '--longitude', dest='ncfile_longitude_name', required=False, type=str,
                    help='Name of the longitude dimension in the netCDF file.')
args = parser.parse_args()
inputFilename=args.um_input_file
gridFilename=args.gridfile
outputFilename=args.um_output_file
latcoord=args.ncfile_latitude_name
loncoord=args.ncfile_longitude_name

print("====== Reading files... ======", end="\r")

import sys
import numpy as np
import cdms2
import regrid2
import mule
from mule.validators import ValidateError
import warnings
import xarray as xr
warnings.filterwarnings("ignore")


# READ FILES
inputFile = mule.AncilFile.from_file(inputFilename)
try:
    gridFile = mule.AncilFile.from_file(gridFilename)
    umgrid=True
except ValueError:
    try:
        gridFile = xr.load_dataset(gridFilename)
        umgrid=False
    except ValueError:
        sys.exit(f"GRIDFILE must be either a UM ancillary file or a netCDF file.")
print("====== Reading files OK! ======")

print("====== Consistency check... ======", end="\r")
# Check that both ancillary file and grid file are valid.
try:
    inputFile.validate()
except ValidateError as e:
    sys.exit(f"Validation failed for the input ancillary file '{inputFilename}'.\nValidator spawned the following error message:\n\n"
            "{}".format('\n'.join(str(e).split('\n')[1:])))

if umgrid:
    try:
        gridFile.validate()
    except ValidateError as e:
        sys.exit(f"Validation failed for the grid file '{gridFilename}'.\nValidator spawned the following error message:\n\n"
            "{}".format('\n'.join(str(e).split('\n')[1:])))
else:
    # Check latitude
    if latcoord is None: #If user has not defined any latitude name
        if "latitude" in gridFile.dims:
            latcoord="latitude"
        elif "lat" in gridFile.dims:
            latcoord="lat"
        else:
            sys.exit("No latitude dimension found in the netCDF file."
                "\nTo specify the name of the latitude dimension in the netCDF file use the option --latitude <name>.")
    elif latcoord not in gridFile.dims:
        sys.exit(f"Specified latitude dimension '{latcoord}' not found in netCDF file.")
    # Check longitude
    if loncoord is None: #If user has not defined any latitude name
        if "longitude" in gridFile.dims:
            loncoord="longitude"
        elif "lon" in gridFile.dims:
            loncoord="lon"
        else:
            sys.exit("No longitude dimension found in the netCDF file."
                "\nTo specify the name of the longitude dimension in the netCDF file use the option --longitude <name>.")
    elif loncoord not in gridFile.dims:
        sys.exit(f"Specified latitude dimension '{loncoord}' not found in netCDF file.")
print("====== Consistency check OK! ======")

print("====== Regridding and writing new UM ancil file... ======", end="\r")
# Get output grid properties
if umgrid:
    nlat_target = gridFile.integer_constants.num_rows
    firstlat_target = gridFile.real_constants.start_lat
    dlat_target = gridFile.real_constants.row_spacing
    lastlat_target = gridFile.real_constants.north_pole_lat
    nlon_target = gridFile.integer_constants.num_cols
    firstlon_target = gridFile.real_constants.start_lon
    dlon_target = gridFile.real_constants.col_spacing
    lastlon_target = gridFile.real_constants.north_pole_lon
    lat0_target = firstlat_target - dlat_target
    lon0_target = firstlon_target - dlon_target
else:
    nlat_target = len(gridFile[latcoord])
    firstlat_target = np.float32(gridFile[latcoord][0].values)
    dlat_target = np.float32((gridFile[latcoord][1]-gridFile[latcoord][0]).values)
    lastlat_target = np.float32(gridFile[latcoord][-1].values)
    nlon_target = len(gridFile[loncoord])
    firstlon_target = np.float32(gridFile[loncoord][0].values)
    dlon_target = np.float32((gridFile[loncoord][1]-gridFile[loncoord][0]).values)
    lastlon_target = np.float32(gridFile[loncoord][-1].values)
    lat0_target = firstlat_target - dlat_target
    lon0_target = firstlon_target - dlon_target

# Output grid
outgrid = cdms2.createUniformGrid(
        startLat=firstlat_target, 
        nlat=nlat_target, 
        deltaLat=dlat_target, 
        startLon=firstlon_target, 
        nlon=nlon_target, 
        deltaLon=dlon_target,
        )

# Create Output file from ancil file
outputFile = inputFile.copy(include_fields=True)
# Change the grid values in the output header to match the chosen origin and size
outputFile.integer_constants.num_rows = nlat_target
outputFile.integer_constants.num_cols = nlon_target
outputFile.real_constants.start_lat = firstlat_target
outputFile.real_constants.start_lon = firstlon_target
outputFile.real_constants.north_pole_lat = lastlat_target
outputFile.real_constants.north_pole_lon = lastlon_target
outputFile.real_constants.row_spacing = dlat_target
outputFile.real_constants.col_spacing = dlon_target

# Loop over all the fields
for f in outputFile.fields:
    # Get input grid properties
    nlat = f.lbrow
    nlon = f.lbnpt
    lat0 = f.bzy
    lon0 = f.bzx
    dlat = f.bdy 
    dlon = f.bdx
    firstlat = lat0 + dlat
    firstlon = lon0 + dlon
    
    # Set output grid properties
    f.lblrec = nlon_target*nlat_target
    f.lbrow = nlat_target
    f.lbnpt = nlon_target
    f.bdy = dlat_target
    f.bdx = dlon_target
    f.bzy = lat0_target
    f.bzx = lon0_target
    
    # Get input grid
    ingrid = cdms2.createUniformGrid(
            startLat=firstlat,
            nlat=nlat, 
            deltaLat=dlat, 
            startLon=firstlon, 
            nlon=nlon, 
            deltaLon=dlon
            )
    
    # Regrid
    regridfunc = regrid2.Horizontal(ingrid, outgrid)
    newdata = regridfunc(f.get_data())

    # If the grid is a global one, force polar values to be the zonal means
    if (inputFile.fixed_length_header.horiz_grid_type == 0 and 
        np.allclose(firstlat, -90) and
        np.allclose(firstlon, 0)):
        newdata[0,:] = newdata[0,:].mean()
        newdata[-1,:] = newdata[-1,:].mean()
    
    f.set_data_provider(mule.ArrayDataProvider(newdata))

# Write output file
if outputFilename is None:
    import os
    outputFilename=inputFilename+"_regridded"
    k=0
    while os.path.isfile(outputFilename):
        k+=1
        outputFilename = outputFilename+str(k)

outputFile.to_file(outputFilename)
print("====== Regridding and writing new UM ancil file OK! ======")
