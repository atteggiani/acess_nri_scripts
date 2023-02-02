#!/usr/bin/env python

# Script created by Davide Marchegiani (davide.marchegiani@anu.edu.au) at ACCESS-NRI.

# Validate a UM ancillary file

def process(file):
    import mule
    import warnings
    import os
    import sys
    from mule.validators import ValidateError
    warnings.filterwarnings("ignore")

    file=os.path.abspath(file)
    data=mule.load_umfile(file)
    if not isinstance(data,mule.ancil.AncilFile):
        sys.exit(f"{file} does not appear to be a UM Ancillary file.")
    try:
        data.validate()
    except ValidateError as e:
        print(f"Failed to validate '{file}'.\n"
        "{}".format('\n'.join(str(e).split('\n')[1:])))
    else:
        print(f"'{file}' is a valid ancillary file!")

if __name__ == '__main__':
    import argparse

    # Parse arguments
    parser = argparse.ArgumentParser(description="Validate UM ancillary file.")
    parser.add_argument('ancilfile', type=str, help='UM ancil file')
    args = parser.parse_args()
    ancilFile=args.ancilfile
    
    process(ancilFile)
