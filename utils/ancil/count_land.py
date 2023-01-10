#!/usr/bin/env python

# Count land points in land_mask UM ancil file
# Created by Davide Marchegiani - davide.marchegiani@anu.edu.au
def process(maskfile):
    import mule
    mask=mule.AncilFile.from_file(maskfile)
    print(mask.fields[0].get_data().sum())

if __name__ == '__main__':
    import argparse

    # Parse arguments
    parser = argparse.ArgumentParser(description="Count land points in land mask UM ancillary file")
    parser.add_argument('maskfile', type=str, help='UM mask file')
    args = parser.parse_args()
    maskfile=args.maskfile
    
    process(maskfile)

