#!/usr/bin/env python

import argparse
# import matplotlib.pyplot as plt
# import myfuncs as my
# import pandas as pd

# Parse arguments
description = '''Modify UM ancillary file.'''
parser = argparse.ArgumentParser(description="Modify UM ancillary file on different grid")
parser.add_argument('--input','-i','-bubbi', dest='um_input_file', required=True, type=str,
                    help='UM ancillary input file')
parser.add_argument('--output','-o', dest='um_output_file', required=True, type=str,
                    help='UM ancillary output file')
parser.add_argument('--ncfile','--nc', dest='ncfile', required=True, type=str,
                    help='NetCDF file to turn into UM ancillary file')
args = parser.parse_args()
ancilFilename=args.um_input_file
# ancilFilename = "/g/data/tm70/dm5220/scripts/abhik/ancil/LANDFRAC"
fileout=args.um_output_file
# fileout = "/g/data/tm70/dm5220/scripts/abhik/ancil_shift/LANDFRAC_shift"
ncFilename=args.ncfile
# ncFilename = "/g/data/tm70/dm5220/scripts/abhik/ancil_shift/nc/MASK_shift.nc"

import mule
import xarray as xr
import numpy as np
import itertools

print("====== Reading files... ======", end="\r")
# Read files
ncFile = xr.load_dataset(ncFilename)
ancilFile = mule.AncilFile.from_file(ancilFilename)
print("====== Reading files OK! ======")

nanval=-1073741824.0
quit()
print("====== Consistency check... ======", end="\r")
# CONSISTENCY CHECK
# Check that both longitude and latitude coords are present in the .nc file and they have consistent lenghts with the ancil file
if "latitude" in ncFile.dims:
    lat=len(ncFile["latitude"])
    if lat != ancilFile.integer_constants.num_rows:
        raise Exception(f"Latitude coord not consistent!\nLength of netCDF file's latitude: {lat}. Length of ancil file's latitude: {ancilFile.integer_constants.num_rows}.")
else:
    raise Exception("NetCDF file does not contain the 'latitude' coordinate.")

if "longitude" in ncFile.dims:
    lon=len(ncFile["longitude"])
    if lon != ancilFile.integer_constants.num_cols:
        raise Exception(f"Longitude coord not consistent!\nLength of netCDF file's longitude: {lon}. Length of ancil file's longitude: {ancilFile.integer_constants.num_cols}.")
else:
    raise Exception("NetCDF file does not contain the 'longitude' coordinate.")

# Check that the time coordinate is present in the .nc file  (if not properly specified) and it has consistent lenght with the ancil file
#========================== ADD SPECIFICATION LINE!!!! ==========================

if ("t" in ncFile.dims):
    tcoord="t"
    time = len(ncFile["t"])    
elif ("time" in ncFile.dims):
    tcoord="time"
    time = len(ncFile["time"])
else:
    raise Exception(f"Time coordinate not present or not understood.\nPlease specify the name of the time coordinate with the option '-t'.")
if time != ancilFile.integer_constants.num_times:
        raise Exception(f"Time coord not consistent!\nLength of netCDF file's time: {time}. Length of ancil file's time: {ancilFile.integer_constants.num_times}.")    

# Check if ancil file has level variable
#========================== ADD LINE!!!! ==========================

# Check if ancil file has pseudo-level
# Check that dimensions are consistent
if sum([f.lbuser5 for f in ancilFile.fields]) == 0: # If no pseudo-levels
    nfields = len(ancilFile.fields) 
    nrecords = len(ncFile.data_vars)*len(ncFile[tcoord])
    if nfields != nrecords:
        raise Exception(f"UM file and netCDF file not consistent.\nUM file has {nfields} fields, but netCDF file has {nrecords} records.")

print("====== Consistency check OK! ======")

print("====== Writing new ancil file... ======", end="\r")
# Change data in the UM ancillary file
newAncilFile = ancilFile.copy(include_fields=True)
ncFile = ncFile.squeeze()

vind = itertools.cycle(ncFile.data_vars)
tind = np.repeat(len(ncFile.data_vars),time)

for i,f in enumerate(newAncilFile.fields):
    var=next(vind)
    d = ncFile[var].where(ncFile[var].notnull(),nanval).values
    # if d.shape != f.shape:
    #     sh1=" x ".join([str(j) for j in d.shape])
    #     sh2=" x ".join([str(j) for j in f.shape])
    #     raise Exception(f"Shape mismatch at field {i}!\nNetCDF shape: {sh1}. AncilFile shape: {sh2}")
    newAncilFile.fields[i].set_data_provider(mule.ArrayDataProvider(d))    

newAncilFile.to_file(fileout)
print("====== Writing new ancil file OK! ======")

