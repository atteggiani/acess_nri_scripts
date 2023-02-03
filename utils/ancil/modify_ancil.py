#!/usr/bin/env python

# Script created by Davide Marchegiani (davide.marchegiani@anu.edu.au) at ACCESS-NRI.

# Modify a UM ancillary file using data from a netCDF file

import argparse
import os

# Parse arguments
description = '''Script to modify a UM ancillary file using data from a netCDF file.'''
parser = argparse.ArgumentParser(description=description, allow_abbrev=False)
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
parser.add_argument('--regrid', dest='regrid', action="store_true",
                    help='Regrid the input UM ancillary file onto the netCDF file grid.')
parser.add_argument('--ignore-levels', dest='ignore_levels', action="store_true",
                    help="Ignore the level dimension found in the ancillary file and \
                          consider it as having only one level.\
                          (Necessary for a mismatch between the actual level dimension \
                          and the 'integer_constans.num_levels' header in some ancillary files.)")
parser.add_argument('--fix', dest='fix_validation', required=False, action='store_true',
                    help='Try to fix any ancillary validation error.')

args = parser.parse_args()
inputFilename=os.path.abspath(args.um_input_file)
ncFilename=os.path.abspath(args.ncfile)
outputFilename=args.um_output_file
latcoord=args.latitude_name
loncoord=args.longitude_name
tcoord=args.time_name
levcoord=args.level_name
pseudocoord=args.pseudo_level_name
nanval=args.nanval
regrid=args.regrid
ignore_levels=args.ignore_levels
fix=args.fix_validation

print(f"====== Reading '{inputFilename}' ancillary file... ======", end="\r")

import mule
from mule.validators import ValidateError
import xarray as xr
import sys
import numpy as np
from fix_validation_error import validate_and_fix, FixValidateError
import warnings
warnings.filterwarnings("ignore")

# Set UM nanval
UM_NANVAL=-1073741824.0 #(-2.0**30)

# READ FILES
inputFile = mule.AncilFile.from_file(inputFilename)
ncFile = xr.load_dataset(ncFilename, decode_times=False).squeeze()
print(f"====== Reading '{inputFilename}' ancillary file OK! ======")


del inputFile.fields[-1]
for f in inputFile.fields:
    f.lbhem=0


# CONSISTENCY CHECK
print("====== Consistency check... ======", end="\r")
# Check that ancillary file is valid.
try:
    if fix:
        inputFile=validate_and_fix(inputFile)
    else:
        inputFile.validate()
except ValidateError as e:
    sys.exit(f"Validation failed for the input ancillary file '{inputFilename}'.\nValidator spawned the following error message:\n"+\
            "{}".format('\n'.join(str(e).split('\n')[1:]))+\
            "If you want to try and fix the error, use the '--fix' option.")
except FixValidateError as e:
    sys.exit(f"Validation and fix failed for the input ancillary file '{inputFilename}'.\n\n{str(e)}")

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
            "\nTo specify the name of the latitude dimension in the netCDF file use the '--latitude <name>' option.")
elif latcoord not in dims:
    sys.exit(f"Specified latitude dimension '{latcoord}' not found in netCDF file.")
nlat=len(ncFile[latcoord])
# Check that latitude dimension is consistent
if (not regrid) and (nlat != inputFile.integer_constants.num_rows):
    sys.exit(f"Latitude dimension not consistent!\nLength of netCDF file's latitude: {nlat}."
        f" Length of ancillary file's latitude: {inputFile.integer_constants.num_rows}."
        "\nIf you want to regrid the ancillary file onto the netCDF grid, please use the '--regrid' option.")
dims.remove(latcoord)

# Longitude
if loncoord is None: #If user has not defined any longitude name
    if "longitude" in dims:
        loncoord="longitude"
    elif "lon" in dims:
        loncoord="lon"
    else:
        sys.exit("No longitude dimension found in the netCDF file."
            "\nTo specify the name of the longitude dimension in the netCDF file use the '--longitude <name>' option.")
elif loncoord not in dims:
    sys.exit(f"Specified longitude dimension '{loncoord}' not found in netCDF file.")
nlon=len(ncFile[loncoord])
# Check that longitude dimension is consistent
if (not regrid) and (nlon != inputFile.integer_constants.num_cols):
    sys.exit(f"Longitude dimension not consistent!\nLength of netCDF file's longitude: {nlon}."
        f" Length of ancillary file's longitude: {inputFile.integer_constants.num_cols}."
        "\nIf you want to regrid the ancillary file onto the netCDF grid, please use the '--regrid' option.")
dims.remove(loncoord)

# Time
if inputFile.integer_constants.num_times > 1:
    if tcoord is None: #If user has not defined any time name
        if "time" in dims:
            tcoord="time"
        elif "t" in dims:
            tcoord="t"
        else:
            sys.exit("No time dimension found in the netCDF file."
                "\nTo specify the name of the time dimension in the netCDF file use the '--time <name>' option.")
    elif tcoord not in dims:
        sys.exit(f"Specified time dimension '{tcoord}' not found in netCDF file.")
    ntime=len(ncFile[tcoord])
    dims.remove(tcoord)
else:
    ntime = 1
# Check that level time is consistent
if ntime != inputFile.integer_constants.num_times:
    sys.exit(f"Time dimension not consistent!\nLength of netCDF file's time: {ntime}."
        f" Length of ancillary file's time: {inputFile.integer_constants.num_times}.")

# Level
if inputFile.integer_constants.num_levels > 1 and not ignore_levels:
    if levcoord is None: #If user has not defined any level name
        if len(dims) >= 1:
            vert_levs=["hybrid","sigma","pressure","depth","surface"]
            # Check if ancillary file has level variable
            for vlev in vert_levs:
                if sum([vlev in s for s in dims]) == 1:
                    levcoord=dims[[vlev in s for s in dims].index(True)]
                    break
                sys.exit(f"Vertical levels found in the ancillary file, but not able to identify vertical level dimension name in the netCDF file."
                    "\nTo specify the name of the vertical level dimension in the netCDF file use the '--level <name>' option.")  
        else:
            sys.exit(f"Vertical levels found in the ancillary file, but no vertical level dimension found in the netCDF file.")
    elif levcoord not in dims:
        sys.exit(f"Specified vertical level dimension '{levcoord}' not found in netCDF file.")
    nlev = len(ncFile[levcoord])
    dims.remove(levcoord)
else:
    if levcoord is not None:
        sys.exit(f"Vertical level dimension '{levcoord}' specified, but no vertical level dimension found in the ancillary file.")
    nlev = 1
# Check that level dimension is consistent
if nlev != inputFile.integer_constants.num_levels and not ignore_levels:
    sys.exit(f"Vertical level dimension not consistent!\nLength of netCDF file's vertical level: {nlev}."
        f" Length of ancillary file's vertical level: {inputFile.integer_constants.num_levels}.")
    
# Pseudo-levels
if sum([f.lbuser5 for f in inputFile.fields]) != 0:
    if pseudocoord is None: #If user has not defined any pseudo-level name
        # Check if ancillary file has pseudo-levels
            if len(dims) == 0:
                sys.exit(f"Pseudo-levels found in the ancillary file, but no pseudo-level dimension found in the netCDF file.")
            elif len(dims) > 1:
                sys.exit(f"Pseudo-levels found in the ancillary file, but not able to identify pseudo-level dimension name in the netCDF file."
                            "\nTo specify the name of the pseudo-level dimension in the netCDF file use the '--pseudo <name>' option.")
            else:
                pseudocoord=dims.pop()
    elif pseudocoord not in dims:
        sys.exit(f"Specified pseudo-level dimension '{pseudocoord}' not found in netCDF file.")
    npseudo = len(ncFile[pseudocoord])
    # Check that pseudo-levels are consistent
    pseudoLevs=[f.lbuser5 for f in inputFile.fields[:(len(inputFile.fields)//ntime)]]
    k=0
    for i,var in enumerate(ncFile.data_vars):
        if pseudocoord in ncFile[var].dims:
            if len(set(pseudoLevs[k:k+npseudo])) != npseudo:
                sys.exit(f"Pseudo dimension not consistent!\n"
                f"Number of pseudo-levels in the ancillary variable n. {i}: {len(set(pseudoLevs[k:k+npseudo]))}."
                f" Length of the '{pseudocoord}' dimension in the netCDF variable n. {i} ('{ncFile[var].name}'): {npseudo}.")
            k+=npseudo
        else:
            k+=1
else:
    if pseudocoord is not None:
        sys.exit(f"Pseudo-level dimension '{pseudocoord}' specified, but no pseudo-level dimension found in the ancillary file.")
    npseudo=1
    
# Check that the number of variables is consistent
nvarsnc = len(ncFile.data_vars)
nvarsancil = len(set([f.lbuser4 for f in inputFile.fields]))
if nvarsnc != nvarsancil:
    sys.exit(f"Number of data variables not consistent!\n"
             f"Number of data variables in the ancillary file: {nvarsancil}. "
             f"Number of data variables in the netCDF file: {nvarsnc}.")

print("====== Consistency check OK! ======")

# WRITE FILE
if outputFilename is None:
    outputFilename=inputFilename+"_modified"
    k=0
    while os.path.isfile(outputFilename):
        k+=1
        outputFilename = outputFilename+str(k)
else:
    outputFilename = os.path.abspath(outputFilename)
print(f"====== Writing '{outputFilename}' ancillary file... ======", end="\r")
# Get the correct shape and substitute netCDF nan with the UM nan value

def substitute_nanval(data):
    if nanval is None:
        return data.where(data.notnull(),UM_NANVAL).transpose(latcoord,loncoord).values
    else:
        return data.where(data != nanval, UM_NANVAL).transpose(latcoord,loncoord).values

def select_time_level(data,tind,lvind):
    if tcoord is not None:
        data = data.isel({tcoord:tind}, drop=True)
    if levcoord is not None:
        data = data.isel({levcoord:lvind}, drop=True)
    return data

# Create a copy of the ancillary file to modify
newAncilFile = inputFile.copy(include_fields=True)

# Exclude extra dimensions
if len(dims) > 0:
    ncFile=ncFile.drop(*dims)

# Regrid parameters
if regrid:
    # Get output grid properties
    nlat_target = nlat
    firstlat_target = np.float32(ncFile[latcoord][0].values)
    dlat_target = np.float32((ncFile[latcoord][1]-ncFile[latcoord][0]).values)
    lastlat_target = np.float32(ncFile[latcoord][-1].values)
    nlon_target = nlon
    firstlon_target = np.float32(ncFile[loncoord][0].values)
    dlon_target = np.float32((ncFile[loncoord][1]-ncFile[loncoord][0]).values)
    lastlon_target = np.float32(ncFile[loncoord][-1].values)
    lat0_target = firstlat_target - dlat_target
    lon0_target = firstlon_target - dlon_target

    # Change the grid values in the output header to match the chosen origin and size
    newAncilFile.integer_constants.num_rows = nlat_target
    newAncilFile.integer_constants.num_cols = nlon_target
    newAncilFile.real_constants.start_lat = firstlat_target
    newAncilFile.real_constants.start_lon = firstlon_target
    newAncilFile.real_constants.north_pole_lat = lastlat_target
    newAncilFile.real_constants.north_pole_lon = lastlon_target
    newAncilFile.real_constants.row_spacing = dlat_target
    newAncilFile.real_constants.col_spacing = dlon_target

fldind = iter(newAncilFile.fields)
for t in range(ntime):
    f=next(fldind)
    if regrid:
        # Set field grid properties
        f.lblrec = nlon_target*nlat_target
        f.lbrow = nlat_target
        f.lbnpt = nlon_target
        f.bdy = dlat_target
        f.bdx = dlon_target
        f.bzy = lat0_target
        f.bzx = lon0_target

    for var in ncFile.data_vars:
        for lv in range(nlev):
            if pseudocoord in ncFile[var].dims:
                for ps in range(npseudo):
                    data_2d = substitute_nanval(select_time_level(ncFile[var],t,lv).isel({pseudocoord:ps},drop=True))
                    f.set_data_provider(mule.ArrayDataProvider(data_2d))
            else:
                data_2d = substitute_nanval(select_time_level(ncFile[var],t,lv))
                f.set_data_provider(mule.ArrayDataProvider(data_2d))

newAncilFile.to_file(outputFilename)
print(f"====== Writing '{outputFilename}' ancillary file OK! ======")
