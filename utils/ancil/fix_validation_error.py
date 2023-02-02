#!/usr/bin/env python

# Script created by Davide Marchegiani (davide.marchegiani@anu.edu.au) at ACCESS-NRI.

# Validate and try to fix a UM ancillary file

import mule
from mule.validators import ValidateError

class FixValidateError(Exception):
    pass

def fix_error(error,ancilFile):
    errorString = '\n'.join(str(error).split('\n')[1:])
    newAncilFile = ancilFile.copy(include_fields=True)
    if errorString == \
            "Ancillary file contains header components other than the "\
            "row/column dependent constants - these should be set to "\
            "\"None\" for Ancillary files":
        newAncilFile.level_dependent_constants = None
        return newAncilFile
    else:
        # sys.exit("The following validation error could not be fixed or is currently not supported:\n\n"
        # f"{errorString}")
        raise FixValidateError("The following validation error could not be fixed or is currently not supported:\n"
        f"{errorString}")

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

