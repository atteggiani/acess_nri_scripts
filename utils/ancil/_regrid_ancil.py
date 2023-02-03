def _regrid_ancil(ancilFile,lat_out=None,lon_out=None,nlev_out=None):
    '''
    Regrids a UM ancilFile over latitude, longitude and UM vertical levels, using scipy.interpolate.interpn function.

    PARAMETERS
    - ancilFile is a mule.ancilFile
    - lat_out is an array-like variable with the output latitude coordinate. If set to None, no regridding will
      be performed over latitude.
    - lon_out is an array-like variable with the output longitude coordinate. If set to None, no regridding will
      be performed over longitude.
    - nlev_out is an integer with the the number of vertical levels in output. If set to None, no regridding 
      will be performed over vertical levels. Pseudo-levels are not counted as vertical levels and regridding
      will not be performed for pseudo-levels. If ancilFile has pseudo-levels, nlev_out needs to be 1.
    '''
    import mule
    import itertools
    import numpy as np
    from scipy.interpolate import interpn
    
    # Parse input file
    if not isinstance(ancilFile,mule.AncilFile):
        raise TypeError("'ancilFile' needs to be a mule.ancilFile object.")
    # Get the input coordinates from the first field of the ancilFile
    f=ancilFile.fields[0].copy()
    lat_in = np.linspace(f.bzy+f.bdy,
                        f.bzy+f.bdy+f.bdy*(ancilFile.integer_constants.num_rows-1),
                        ancilFile.integer_constants.num_rows)
    lon_in = np.linspace(f.bzx+f.bdx,
                        f.bzx+f.bdx+f.bdx*(ancilFile.integer_constants.num_cols-1),
                        ancilFile.integer_constants.num_cols)
    lev_in = np.arange(1,ancilFile.integer_constants.num_levels+1)
    ntimes = ancilFile.integer_constants.num_times
    lbegin = f.lbegin
    # Check if ancil file has pseudolevs
    pseudoLevs = [f.lbuser5 for f in ancilFile.fields[:(len(ancilFile.fields)//ntimes)]]
    pseudoLevs = [val for i,val in enumerate(pseudoLevs[:-1]) if val == 0 or pseudoLevs[i+1]<pseudoLevs[i]]+[pseudoLevs[-1]]
    if sum(pseudoLevs) == 0:
        pseudoLevs = [1]
    elif nlev_out != 1:
        raise ValueError("Pseudo-levels found in the ancilFile, but 'nlev_out' is not 1.")
    else:
        pseudoLevs = [l+1 if l == 0 else l for l in pseudoLevs]
    # Parse output coords 
    # Latitude
    if lat_out is None:
        lat_out = lat_in
    elif isinstance(lat_out,(int, float)) and not isinstance(lat_out, bool):
        lat_out = [lat_out]
    elif not hasattr(lat_out,'__iter__'):
        raise TypeError("'lat_out' needs to be an iterable.")
    # Longitude
    if lon_out is None:
        lon_out = lon_in
    elif isinstance(lon_out,(int, float)) and not isinstance(lon_out, bool):
        lon_out = [lon_out]
    elif not hasattr(lon_out,'__iter__'):
        raise TypeError("'lon_out' needs to be an iterable.")
    # Vertical levels
    if nlev_out is None:
        lev_out = lev_in        
    elif not isinstance(nlev_out,int):
        raise TypeError("'nlev_out' needs to be an integer.")
    else:
        lev_out = np.linspace(1,ancilFile.integer_constants.num_levels,nlev_out)
    
    nlat_out = len(lat_out)
    nlon_out = len(lon_out)
    if len(lat_out) > 1:
        dlat_out = lat_out[1] - lat_out[0]
    else:
        dlat_out = 180.
    if len(lon_out) > 1:
        dlon_out = lon_out[1] - lon_out[0]
    else:
        dlon_out = 360.
    
    outpoints = list(itertools.product(lat_out,lon_out,lev_out))
    
    # Create new ancil file 
    regriddedFile = ancilFile.copy(include_fields=False)
    # Change regridded file header
    regriddedFile.integer_constants.num_rows = nlat_out
    regriddedFile.integer_constants.num_cols = nlon_out
    regriddedFile.integer_constants.num_levels = nlev_out
    regriddedFile.real_constants.start_lat = lat_out[0]
    regriddedFile.real_constants.start_lon = lon_out[0]
    regriddedFile.real_constants.north_pole_lat = lat_out[-1]
    regriddedFile.real_constants.north_pole_lon = lon_out[-1]
    regriddedFile.real_constants.row_spacing = dlat_out
    regriddedFile.real_constants.col_spacing = dlon_out

    f = 0
    interp=[]
    while f < len(ancilFile.fields):
        data = []
        for _ in range(ancilFile.integer_constants.num_levels):
            data.append(ancilFile.fields[f].copy().get_data())
            f += 1
        values = np.stack(data,axis=2)
        interp.append(interpn((lat_in,lon_in,lev_in), values, outpoints, bounds_error=False, fill_value=None).reshape(nlat_out,nlon_out,nlev_out))
    newfields = np.stack(interp,axis=3).reshape(nlat_target,nlon_target,-1)
    # If the grid is a global one, force polar values to be the zonal means
    if (ancilFile.fixed_length_header.horiz_grid_type == 0 and 
        np.allclose(lat_out[0], -90) and
        np.allclose(lon_out[0], 0)):
        newfields[0,...]=newfields[0,...].mean(axis=0)
        newfields[-1,...]=newfields[-1,...].mean(axis=0)

    k=0
    fldind = [0]+np.cumsum(pseudoLevs).tolist()[:-1]
    for _ in range(ntimes):
        for _ in range(nlev_out):
            for ips,nl in enumerate(pseudoLevs):
                for l in range(1,nl+1):
                    regriddedFile.fields.append(ancilFile.fields[fldind[ips]].copy())
                    f = regriddedFile.fields[k]
                    # Change field pseudo-level
                    if len(pseudoLevs) > 1:
                        f.lbuser5 = l
                    else:
                        f.lbuser5 = 0
                    # Change field headers relative to coordinates
                    f.lblrec = nlat_out*nlon_out
                    f.lbuser2 = 1+(f.lblrec*k)
                    f.lbegin = lbegin + (f.lbnrec*k)
                    f.lbrow = nlat_out
                    f.lbnpt = nlon_out
                    f.bdy = dlat_out
                    f.bdx = dlon_out
                    f.bzy = lat_out[0] - dlat_out
                    f.bzx = lon_out[0] - dlon_out
                    f.set_data_provider(mule.ArrayDataProvider(newfields[:,:,k]))
                    k+=1
    return regriddedFile

d=_regrid_ancil(ancilFile,lat_out,lon_out,nlev_out)