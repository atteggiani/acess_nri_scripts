#!/usr/bin/env python

# Interpolate an ancillary file to another grid
# Output word size and endianness match input.

# Converts all variables in the file (assumed to all be on the same grid)

# Original version by Martin Dix - martin.dix@csiro.au
# Modified by Davide Marchegiani - davide.marchegiani@anu.edu.au

import sys
sys.path.append("/g/data/access/projects/access/apps/pythonlib/umfile_utils")
import numpy as np
import argparse
import umfile
from um_fileheaders import *
import cdms2
import regrid2

# Parse arguments
parser = argparse.ArgumentParser(description="Interpolate UM ancillary file on different grid")
parser.add_argument('-i', dest='ifile', required=True, type=str,
                    help='UM input file')
parser.add_argument('-o', dest='ofile', required=True, type=str,
                    help='UM output file')
parser.add_argument('-m', dest='gridfile', required=True, type=str,
                    help='UM file on the new grid')
args = parser.parse_args()
ifile=args.ifile
ofile=args.ofile
gridfile=args.gridfile

m = umfile.UMFile(gridfile)
f = umfile.UMFile(ifile)

# Target grid properties
nlon_target = m.inthead[IC_XLen]
nlat_target = m.inthead[IC_YLen]

# Create new axes for the output grid
outgrid = cdms2.createUniformGrid(m.realhead[RC_FirstLat], nlat_target, m.realhead[RC_LatSpacing], m.realhead[RC_FirstLong], nlon_target, m.realhead[RC_LongSpacing])

g = umfile.UMFile(ofile, "w")
g.copyheader(f)

# Change the grid values in the output header to match the chosen origin and
# size
g.inthead[IC_XLen] = nlon_target
g.inthead[IC_YLen] = nlat_target
g.realhead[RC_FirstLat] = m.realhead[RC_FirstLat]
g.realhead[RC_FirstLong] = m.realhead[RC_FirstLong]
g.realhead[RC_PoleLong] = m.realhead[RC_PoleLong]
g.realhead[RC_PoleLat] = m.realhead[RC_PoleLat]
g.realhead[RC_LatSpacing] = m.realhead[RC_LatSpacing]
g.realhead[RC_LongSpacing] = m.realhead[RC_LongSpacing]
lat0 = g.realhead[RC_FirstLat]  - g.realhead[RC_LatSpacing]
lon0 = g.realhead[RC_FirstLong] - g.realhead[RC_LongSpacing]

# Loop over all the fields
for k in range(f.fixhd[FH_LookupSize2]):
    ilookup = f.ilookup[k]
    rlookup = f.rlookup[k]
    lbegin = ilookup[LBEGIN] # lbegin is offset from start
    if lbegin == -99:
        break
    npts = ilookup[LBNPT]
    nrows = ilookup[LBROW]

    # Set modified output grid for this field
    g.ilookup[k,LBLREC] = nlon_target*nlat_target
    g.ilookup[k,LBROW] = nlat_target
    g.ilookup[k,LBNPT] = nlon_target
    # Need to hold float values in an integer array
    g.rlookup[k,BDY] = g.realhead[RC_LatSpacing]
    g.rlookup[k,BDX] = g.realhead[RC_LongSpacing]
    g.rlookup[k,BZY] = lat0
    g.rlookup[k,BZX] = lon0

    data = f.readfld(k)

    # May be different to the overall file settings if it's not a proper
    # ancillary file
    lat1 = rlookup[BZY] + rlookup[BDY]
    lon1 = rlookup[BZX] + rlookup[BDX]
    ingrid = cdms2.createUniformGrid(lat1, nrows, rlookup[BDY], lon1, npts, rlookup[BDX])

    regridfunc = regrid2.Regridder(ingrid, outgrid)

    newdata = regridfunc(data)
    if newdata.dtype == np.float32:
        # Output from regrid is float32, so set packing to match
        g.ilookup[k,LBPACK] = (f.ilookup[k,LBPACK]//10)*10 + 2
    else:
        # Don't expect this case but handle it anyway
        g.ilookup[k,LBPACK] = (f.ilookup[k,LBPACK]//10)*10
    
    # If this is a global grid force polar values to be the zonal means
    if ( f.fixhd[FH_HorizGrid] == 0 and
         np.allclose(rlookup[BZY] + rlookup[BDY], -90.) and
         np.allclose(rlookup[BZX] + rlookup[BDX], 0.) ):
        newdata[0,:] = newdata[0,:].mean()
        newdata[-1,:] = newdata[-1,:].mean()

    g.writefld(newdata,k)

g.close()
