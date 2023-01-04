import numpy as np
import xarray as xr
import os
import mule
from matplotlib import pyplot as plt

def plot(x,**kwargs):
    plt.figure()
    plt.contourf(x.get_data(),**kwargs)
    plt.axis("tight")
    plt.colorbar()

file="veg_modified.nc"
ancil_file = "/g/data/tm70/dm5220/scripts/abhik/veg_original"

ancil=mule.AncilFile.from_file(ancil_file)
ancil_modif=ancil.copy(include_fields=True)

nanval=ancil.fields[-12].get_data()[0,0]

d=xr.open_dataset(file).squeeze()
k=0
for v in d.data_vars:
    # cond=d[v]
    cond=d[v].where(~d[v].isnull(),nanval)
    data=cond.values
    
    for i in range(data.shape[0]):
        ancil_modif.fields[k+i].set_data_provider(mule.ArrayDataProvider(data[i,...]))
    k+=data.shape[0]

newfile=os.path.splitext(file)[0]
ancil_modif.to_file(newfile)

# test
newancil=mule.AncilFile.from_file(newfile)

