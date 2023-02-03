#!/usr/bin/env python

# # Script screated by Davide Marchegiani ar ACCESS-NRI - davide.marchegiani@anu.edu.au

# # Interpolate an ancillary file onto another grid

import argparse
import os
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
parser.add_argument('--lev', '--level', dest='ncfile_level_name', required=False, type=str,
                    help='Name of the level dimension in the netCDF file.')
parser.add_argument('--fix', dest='fix_validation', required=False, action='store_true',
                    help='Try to fix any ancillary validation error.')

args = parser.parse_args()
inputFilename=os.path.abspath(args.um_input_file)
gridFilename=os.path.abspath(args.gridfile)
outputFilename=args.um_output_file
latcoord=args.ncfile_latitude_name
loncoord=args.ncfile_longitude_name
levcoord=args.ncfile_level_name
fix=args.fix_validation

print(f"====== Reading '{inputFilename}' ancillary file... ======", end="\r")

import sys
import numpy as np
import cdms2
import regrid2
import mule
from mule.validators import ValidateError
import xarray as xr
from fix_validation_error import validate_and_fix, FixValidateError
import warnings
warnings.filterwarnings("ignore")

inputFilename='/g/data3/tm70/dm5220/ancil/abhik/ancil-from-uk/ozone/sparc/1994-2005/qrclim.ozone_L85_O85'
gridFilename='/g/data/access/projects/access/data/ancil/n48_hadgem1/RSPARC_O3_1990.anc'
outputFilename='/g/data3/tm70/dm5220/ancil/abhik/newancil/soil_temp/amip/qrclim.slt'
fix=True

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
print(f"====== Reading '{inputFilename}' ancillary file OK! ======")

print("====== Consistency check... ======", end="\r")
# Check that both ancillary file and grid file are valid.
try:
    if fix:
        inputFile=validate_and_fix(inputFile)
    else:
        inputFile.validate()
except ValidateError as e:
    sys.exit(f"Validation failed for the input ancillary file '{inputFilename}'.\nValidator spawned the following error message:\n"+\
            "'{}'\n".format('\n'.join(str(e).split('\n')[1:]))+\
            "If you want to try and fix the error, use the '--fix' option.")
except FixValidateError as e:
    sys.exit(f"Validation and fix failed for the input ancillary file '{inputFilename}'.\n\n{str(e)}")

# Check that dimensions are consistent
nlat_input = inputFile.integer_constants.num_rows
nlon_input = inputFile.integer_constants.num_cols
nlev_input = inputFile.integer_constants.num_levels
if umgrid: #If the grid is a um ancil file
    try:
        if fix:
            gridFile=validate_and_fix(gridFile)
        else:
            gridFile.validate()
    except ValidateError as e:
        sys.exit(f"Validation failed for the grid file '{gridFilename}'.\nValidator spawned the following error message:\n"+\
                "'{}'\n".format('\n'.join(str(e).split('\n')[1:]))+\
                "If you want to try and fix the error, use the '--fix' option.")
    except FixValidateError as e:
        sys.exit(f"Validation and fix failed for the grid file '{gridFilename}'.\n\n{str(e)}")
    # Check that both InputFile and gridFile are consistent in case latitude, longitude or levels are of length 1
    # latitude
    nlat_target = gridFile.integer_constants.num_rows
    firstlat_target = gridFile.real_constants.start_lat
    dlat_target = gridFile.real_constants.row_spacing
    lastlat_target = gridFile.real_constants.north_pole_lat
    lat0_target = firstlat_target - dlat_target
    # longitude
    nlon_target = gridFile.integer_constants.num_cols
    firstlon_target = gridFile.real_constants.start_lon
    dlon_target = gridFile.real_constants.col_spacing
    lastlon_target = gridFile.real_constants.north_pole_lon
    lon0_target = firstlon_target - dlon_target
    # level
    nlev_target = gridFile.integer_constants.num_levels
else: #If the grid is a netCDF file
    # Check latitude
    if latcoord is None: #If user has not defined any latitude name
        if "latitude" in gridFile.dims:
            latcoord="latitude"
        elif "lat" in gridFile.dims:
            latcoord="lat"
        else:
            sys.exit("No latitude dimension found in the netCDF file."
                "\nTo specify the name of the latitude dimension in the netCDF file use the '--latitude <name>' option.")
    elif latcoord not in gridFile.dims:
        sys.exit(f"Specified latitude dimension '{latcoord}' not found in netCDF file.")
    nlat_target = len(gridFile[latcoord])
    firstlat_target = np.float32(gridFile[latcoord][0].values)
    dlat_target = np.float32((gridFile[latcoord][1]-gridFile[latcoord][0]).values)
    lastlat_target = np.float32(gridFile[latcoord][-1].values)
    lat0_target = firstlat_target - dlat_target
    # Check longitude
    if loncoord is None: #If user has not defined any latitude name
        if "longitude" in gridFile.dims:
            loncoord="longitude"
        elif "lon" in gridFile.dims:
            loncoord="lon"
        else:
            sys.exit("No longitude dimension found in the netCDF file."
                "\nTo specify the name of the longitude dimension in the netCDF file use the '--longitude <name>' option.")
    elif loncoord not in gridFile.dims:
        sys.exit(f"Specified latitude dimension '{loncoord}' not found in netCDF file.")
    nlon_target = len(gridFile[loncoord])
    firstlon_target = np.float32(gridFile[loncoord][0].values)
    dlon_target = np.float32((gridFile[loncoord][1]-gridFile[loncoord][0]).values)
    lastlon_target = np.float32(gridFile[loncoord][-1].values)    
    lon0_target = firstlon_target - dlon_target
    # Check level
    if levcoord is None: #If user has not defined any level name
        vert_levs=["hybrid","sigma","pressure","depth","surface"]
        # Check if ancillary file has level variable
        for vlev in vert_levs:
            if sum([vlev in s for s in gridFile.dims]) == 1:
                levcoord=gridFile.dims[[vlev in s for s in gridFile.dims].index(True)]
                break
            sys.exit(f"No vertical level dimension found in the netCDF file."
                "\nTo specify the name of the vertical level dimension in the netCDF file use the '--level <name>' option.")
    elif levcoord not in gridFile.dims:
        sys.exit(f"Specified vertical level dimension '{levcoord}' not found in netCDF file.")
    nlev = len(gridFile.dims[levcoord])

# Check Latitude
if (nlat_input == 1 and nlat_target != 1):
    sys.exit(f"Latitude is inconsistent! Ancillary file '{inputFilename}' has no latitude dimension "+\
            f"but the grid file '{gridFilename}' has a latitude length of {nlat_target}.")
elif (nlat_input != 1 and nlat_target == 1):
    sys.exit(f"Latitude is inconsistent! Grid file '{gridFilename}' has no latitude dimension "+\
            f"but the ancillary file '{inputFilename}' has a latitude length of {nlat_input}.")
# Check Longitude
if (nlon_input == 1 and nlon_target != 1):
    sys.exit(f"Longitude is inconsistent! Ancillary file '{inputFilename}' has no longitude dimension "+\
            f"but the grid file '{gridFilename}' has a longitude length of {nlon_target}.")
elif (nlon_input != 1 and nlon_target == 1):
    sys.exit(f"Longitude is inconsistent! Grid file '{gridFilename}' has no longitude dimension "+\
            f"but the ancillary file '{inputFilename}' has a longitude length of {nlon_input}.")
# Check Level
if (nlev_input == 1 and nlev_target != 1):
    sys.exit(f"Levels are inconsistent! Ancillary file '{inputFilename}' has no vertical level dimension "+\
            f"but the grid file '{gridFilename}' has {nlev_target} vertical levels.")
elif (nlev_input != 1 and nlev_target == 1):
    sys.exit(f"Levels are inconsistent! Grid file '{gridFilename}' has no vertical level dimension "+\
            f"but the ancillary file '{inputFilename}' has {nlev_input} vertical levels.")
elif nlev_input > 1:
    vertical_regrid=True
else:
    vertical_regrid=False
print("====== Consistency check OK! ======")


if outputFilename is None:
    outputFilename=inputFilename+"_regridded"
    k=0
    while os.path.isfile(outputFilename):
        k+=1
        outputFilename = outputFilename+str(k)
else:
    outputFilename = os.path.abspath(outputFilename)

print(f"====== Regridding and writing '{outputFilename}' ancillary file... ======", end="\r")
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

outputFile.to_file(outputFilename)
print(f"====== Regridding and writing '{outputFilename}' ancillary file OK! ======")
