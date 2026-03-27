#Compliments of Will Lenthe 
import h5py
import numpy as np
import struct

input_file = "/Volumes/Drive2/playarea/Oxford/Al-large.h5"
input_dataset = "1/EBSD/Data/Processed Patterns"

# open input file and read dataset
pats = h5py.File(input_file)[input_dataset]
output_file = '.'.join(input_file.split('.')[:-1]) + ".up"

if pats.dtype == np.dtype("u1"):
                output_file += '1'
elif pats.dtype == np.dtype("u2"):
                output_file += '2'
else:
                raise TypeError("unsupported datatype for up1/up2 output file")

 

# write data to output file

with open(output_file, 'wb') as f:
# may need to switch pats.shape[1] and pats.shape[2] if the patterns are from a non-square detector
                header = [1, pats.shape[2], pats.shape[1], 16] # version, pattern width, pattern height, offset to data start
                print (header)
                f.write(struct.pack('4i', *header))
                f.write(np.array(pats).tobytes())
