#!/usr/bin/env python

# Script created by Davide Marchegiani (davide.marchegiani@anu.edu.au) at ACCESS-NRI.

# Validate and try to fix a UM ancillary file

import mule
from mule.validators import ValidateError
# /g/data3/hh5/public/apps/miniconda3/envs/analysis3-22.07/lib/python3.9/site-packages/mule/validators.py
class FixValidateError(Exception):
    pass

def fix_error(error,ancilFile):
    errorString = '\n'.join(str(error).split('\n')[1:])
    newAncilFile = ancilFile.copy(include_fields=True)
    if errorString.startswith("Ancillary file contains header components other than the"):
        newAncilFile.level_dependent_constants = None
    elif errorString.startswith("Unsupported grid_staggering"):
        newAncilFile.fixed_length_header.grid_staggering = 3
    else:
        raise FixValidateError("The following validation error could not be fixed or is currently not supported:\n"
        f"'{errorString}'")
    return newAncilFile

def validate_and_fix(ancilFile):
    if not isinstance(ancilFile,mule.AncilFile):
        raise ValueError("The file provided needs to be a mule.AncilFile.")
    try:
        ancilFile.validate()
    except ValidateError as e:
        newAncilFile = fix_error(e,ancilFile)
        return newAncilFile
    else:
        return ancilFile

