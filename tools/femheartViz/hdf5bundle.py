
import numpy as np
import h5py
import re

def main():
    import argparse
    import os
    parser = argparse.ArgumentParser()
    parser.add_argument("outdir", help="Directory used for femheart output")
    args = parser.parse_args()
    
    outdir = args.outdir
    outdir = os.path.realpath(outdir)
    print outdir
    #Make HDF5 file
    h5file = h5py.File(outdir+".h5", 'w')
    #Scan output for timesteps
    for possibleTime in os.listdir(outdir):
        #for each timestep
        timedir = os.path.join(outdir,possibleTime)
        if (re.match(r"^tm[0-9]+", possibleTime)
            and
            os.path.isdir(timedir)):
            #scan directory for npy files
            for filename in os.listdir(timedir):
                m = re.match(r'^(.*)\.npy',filename)
                if m:
                    arrayName = m.group(1)
                    arrayFile = os.path.join(timedir,filename)
                    if possibleTime not in h5file:
                        h5file.create_group(possibleTime)
                    AAA = np.load(arrayFile)
                    h5file[possibleTime][arrayName] = AAA

if __name__=='__main__':
    main()
