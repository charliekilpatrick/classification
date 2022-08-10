#!/usr/bin/env python3
from astropy.table import Table
import glob
import os
import classification

datadir = 'data'
objtable = Table.read('typeII_sample.csv', format='csv')

seps = []
for row in objtable:

    filename = glob.glob(os.path.join(datadir, row['Name']+'*.cat'))
    if len(filename)>0:
        file = os.path.basename(filename[0])
        distance = float(row['Distance'].split()[0])
        sep = classification.do_classification(file,
            distance, row['RA'], row['Dec'])

        seps.append((row['Name'], sep))

for sep in seps:
    print(*sep)
