#!/usr/bin/env python

import mule
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import itertools

veg_file="/g/data/tm70/dm5220/scripts/abhik/ancil/VEGINIT"
veg=mule.AncilFile.from_file(veg_file)

# Change veg file so it accounts for multiple functional types + minor other changes
new_veg=veg.copy(include_fields=True)
# Change model version
new_veg.fixed_length_header.model_version = 703
# Remove total_prognostic_fields (not used)
new_veg.fixed_length_header.total_prognostic_fields = -32768
# Change number of fields to 3
new_veg.integer_constants.num_field_types = 3

functype=itertools.cycle(range(1,6))
# Modify for every field
# field=new_veg.fields[-1].copy()
# new_veg.fields.append(field)
for f in new_veg.fields:
    # Set 5 Pseudo-levels for functional types but to the canopy (stash 213)
    if f.lbuser4 != 213:
        f.lbuser5 = next(functype)
    else:
        f.lbuser5 = 0

new_veg.to_file(veg_file+"_modif")