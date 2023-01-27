#!/usr/bin/env python

# Script created by Davide Marchegiani (davide.marchegiani@anu.edu.au) at ACCESS-NRI.

# Modify a UM ancillary file using data from a netCDF file

import argparse

# Parse arguments
description = '''Script to modify a UM ancillary file using data from a netCDF file.'''
parser = argparse.ArgumentParser(description=description)
parser.add_argument('-i', '--input', dest='um_input_file', required=True, type=str,
                    help='UM ancillary input file. (Required)')
parser.add_argument('--nc', '--ncfile', dest='ncfile', required=True, type=str,
                    help='NetCDF file to turn into UM ancillary file. (Required)')
parser.add_argument('-o', '--output', dest='um_output_file', required=False, type=str,
                    help='UM ancillary output file.')                    
parser.add_argument('--lat', '--latitude', dest='latitude_name', required=False, type=str,
                    help='Name of the latitude dimension in the netCDF file.')
parser.add_argument('--lon', '--longitude', dest='longitude_name', required=False, type=str,
                    help='Name of the longitude dimension in the netCDF file.')
parser.add_argument('-t', '--time', dest='time_name', required=False, type=str,
                    help='Name of the time dimension in the netCDF file.')
parser.add_argument('--lev', '--level', dest='level_name', required=False, type=str,
                    help='Name of the level dimension in the netCDF file.') 
parser.add_argument('--ps', '--pseudo', dest='pseudo_level_name', required=False, type=str,
                    help='Name of the pseudo-level dimension in the netCDF file.')
parser.add_argument('--nan', dest='nanval', required=False, type=float,
                    help='Value for NaNs in the netCDF file.')

args = parser.parse_args()
inputFilename=args.um_input_file
ncFilename=args.ncfile
outputFilename=args.um_output_file
latcoord=args.latitude_name
loncoord=args.longitude_name
tcoord=args.time_name
levcoord=args.level_name
pseudocoord=args.pseudo_level_name
nanval=args.nanval

import mule
from mule.validators import ValidateError
import xarray as xr
import warnings
import sys
warnings.filterwarnings("ignore")

# Set UM nanval
UM_NANVAL=-1073741824.0 #(-2.0**30)

# READ FILES
print("====== Reading files... ======", end="\r")
ncFile = xr.load_dataset(ncFilename).squeeze()
inputFile = mule.AncilFile.from_file(inputFilename)
print("====== Reading files OK! ======")


# CONSISTENCY CHECK
print("====== Consistency check... ======", end="\r")
# Check that ancillary file is valid.
try:
    inputFile.validate()
except ValidateError as e:
    sys.exit(f"Validation failed for the input ancillary file '{inputFilename}'.\nValidator spawned the following error message:\n\n"
            "{}".format('\n'.join(str(e).split('\n')[1:])))

# Check that longitude, latitude and time coords are present in the .nc file and they have consistent lenghts with the ancillary file
dims=list(ncFile.dims)

# Latitude
if latcoord is None: #If user has not defined any latitude name
    if "latitude" in dims:
        latcoord="latitude"
    elif "lat" in dims:
        latcoord="lat"
    else:
        sys.exit("No latitude dimension found in the netCDF file."
            "\nTo specify the name of the latitude dimension in the netCDF file use the option --latitude <name>.")
elif latcoord not in dims:
    sys.exit(f"Specified latitude dimension '{latcoord}' not found in netCDF file.")
lat=len(ncFile[latcoord])
# Check that latitude dimension is consistent
if lat != inputFile.integer_constants.num_rows:
    sys.exit(f"Latitude dimension not consistent!\nLength of netCDF file's latitude: {lat}."
        f" Length of ancillary file's latitude: {inputFile.integer_constants.num_rows}.")
dims.remove(latcoord)

# Longitude
if loncoord is None: #If user has not defined any longitude name
    if "longitude" in dims:
        loncoord="longitude"
    elif "lon" in dims:
        loncoord="lon"
    else:
        sys.exit("No longitude dimension found in the netCDF file."
            "\nTo specify the name of the longitude dimension in the netCDF file use the option --longitude <name>.")
elif loncoord not in dims:
    sys.exit(f"Specified longitude dimension '{loncoord}' not found in netCDF file.")
lon=len(ncFile[loncoord])
# Check that longitude dimension is consistent
if lon != inputFile.integer_constants.num_cols:
    sys.exit(f"Longitude dimension not consistent!\nLength of netCDF file's longitude: {lon}."
        f" Length of ancillary file's longitude: {inputFile.integer_constants.num_cols}.")
dims.remove(loncoord)

# Time
if tcoord is None: #If user has not defined any time name
    if "time" in dims:
        tcoord="time"
    elif "t" in dims:
        tcoord="t"
    else:
        sys.exit("No time dimension found in the netCDF file."
            "\nTo specify the name of the time dimension in the netCDF file use the option --time <name>.")
elif tcoord not in dims:
    sys.exit(f"Specified time dimension '{tcoord}' not found in netCDF file.")
time=len(ncFile[tcoord])
# Check that level time is consistent
if time != inputFile.integer_constants.num_times:
    sys.exit(f"Time dimension not consistent!\nLength of netCDF file's time: {time}."
        f" Length of ancillary file's time: {inputFile.integer_constants.num_times}.")
dims.remove(tcoord)

# Level
if inputFile.integer_constants.num_levels > 1:
    if levcoord is None: #If user has not defined any level name
        if len(dims) >= 1:
            vert_levs=["hybrid","sigma","pressure","depth"]
            # Check if ancillary file has level variable
            for vlev in vert_levs:
                if sum([vlev in s for s in dims]) == 1:
                    levcoord=dims[[vlev in s for s in dims].index(True)]
                    break
                sys.exit(f"Not able to identify vertical level dimension name."
                    "\nTo specify the name of the vertical level dimension in the netCDF file use the option --level <name>.")  
        else:
            sys.exit(f"Vertical levels found in the ancillary file, but no vertical level dimension found in the netCDF file.")
    elif levcoord not in dims:
        sys.exit(f"Specified vertical level dimension '{levcoord}' not found in netCDF file.")
    lev = len(ncFile[levcoord])
    # Check that level dimension is consistent
    if lev != inputFile.integer_constants.num_levels:
        sys.exit(f"Vertical level dimension not consistent!\nLength of netCDF file's vertical level: {lev}."
            f" Length of ancillary file's vertical level: {inputFile.integer_constants.num_levels}.")
    dims.remove(levcoord)
else:
    if levcoord is not None:
        sys.exit(f"Vertical level dimension '{levcoord}' specified, but no vertical level dimension found in the ancillary file.")

# Pseudo-levels
if sum([f.lbuser5 for f in inputFile.fields]) != 0:
    if pseudocoord is None: #If user has not defined any pseudo-level name
        # Check if ancillary file has pseudo-levels
            if len(dims) == 0:
                sys.exit(f"Pseudo-levels found in the ancillary file, but no pseudo-level dimension found in the netCDF file.")
            elif len(dims) > 1:
                sys.exit(f"Not able to identify pseudo-level dimension name."
                            "\nTo specify the name of the pseudo-level dimension in the netCDF file use the option --pseudo <name>.")
            else:
                pseudocoord=dims.pop()
    elif pseudocoord not in dims:
        sys.exit(f"Specified pseudo-level dimension '{levcoord}' not found in netCDF file.")
    pseudo = len(ncFile[pseudocoord])
    # Check that pseudo-levels are consistent
    pseudoLevs=[f.lbuser5 for f in inputFile.fields[:(len(inputFile.fields)//time)]]
    k=0
    for i,var in enumerate(ncFile.data_vars):
        if pseudo in ncFile[var].dims:
            if len(set(pseudoLevs[k:k+pseudo])) != pseudo:
                sys.exit(f"Pseudo dimension not consistent!\n"
                f"Num of pseudo-levels in the ancillary field n. {i}: {len(set(pseudoLevs[k:k+pseudo]))}."
                f" Length of the netCDF variable n. {i} ('{ncFile[var].name}'): {pseudo}.")
            k+=pseudo
        else:
            k+=1
else:
    if pseudocoord is not None:
        sys.exit(f"Pseudo-level dimension '{pseudocoord}' specified, but no pseudo-level dimension found in the ancillary file.")

print("====== Consistency check OK! ======")

# WRITE FILE
print("====== Writing new ancillary file... ======", end="\r")
# Get the correct shape and substitute netCDF nan with the UM nan value
def get_2D_data(data):
    if nanval is None:
        return data.where(data.notnull(),UM_NANVAL).transpose(latcoord,loncoord).values
    else:
        return data.where(data != nanval, UM_NANVAL).transpose(latcoord,loncoord).values

# Create a copy of the ancillary file to modify
if len(dims) > 0:
    newAncilFile = inputFile.copy(include_fields=True).drop(*dims)
else:
    newAncilFile = inputFile.copy(include_fields=True)

fldind = iter(newAncilFile.fields)
for time in range(time):
    if pseudocoord is None and levcoord is None:
        for var in ncFile.data_vars:
            data_2d = get_2D_data(ncFile[var].isel({tcoord:time},drop=True))
            next(fldind).set_data_provider(mule.ArrayDataProvider(data_2d))
    elif pseudocoord is None:
        for var in ncFile.data_vars:
            for lv in range(lev):
                data_2d = get_2D_data(ncFile[var].isel({levcoord:lv,tcoord:time},drop=True))
                next(fldind).set_data_provider(mule.ArrayDataProvider(data_2d))
    elif levcoord is None:
        for var in ncFile.data_vars:
            if pseudocoord in ncFile[var].dims:
                for ps in range(pseudo):
                    data_2d = get_2D_data(ncFile[var].isel({pseudocoord:ps,tcoord:time},drop=True))
                    next(fldind).set_data_provider(mule.ArrayDataProvider(data_2d))
    else:
        for var in ncFile.data_vars:
            for lv in range(lev):
                if pseudocoord in ncFile[var].dims:
                    for ps in range(pseudo):
                        data_2d = get_2D_data(ncFile[var].isel({levcoord:lv,pseudocoord:ps,tcoord:time},drop=True))
                        next(fldind).set_data_provider(mule.ArrayDataProvider(data_2d))

if outputFilename is None:
    import os
    outputFilename=inputFilename+"_modified"
    k=0
    while os.path.isfile(outputFilename):
        k+=1
        outputFilename = outputFilename+str(k)
newAncilFile.to_file(outputFilename)
print("====== Writing new ancillary file OK! ======")
