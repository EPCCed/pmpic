#!/usr/bin/env python
from __future__ import print_function
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys


if (len(sys.argv) != 3):
    print("Usage: ./visualise [frame_number] [num_processes]")
    sys.exit()


def read_filename(filename):
    f=open(filename,"rb")
    time=np.fromfile(f,dtype=np.float64,count=1)
    n=np.fromfile(f,dtype=np.int64,count=1)

    n=n[0]

    print("Number of parcels = %d"%n)

    x=np.fromfile(f,dtype=np.float64,count=n)
    y=np.fromfile(f,dtype=np.float64,count=n)
    z=np.fromfile(f,dtype=np.float64,count=n)
    tag=np.fromfile(f,dtype=np.float64,count=n)

    f.close()

    return x, y, z, tag, n


nt=int(sys.argv[1])
num_p=int(sys.argv[2])

xs=[]
ys=[]
zs=[]
tags=[]

for p in range(num_p):
    filename="parcels_%03d_%04d.dat"%(p,nt)
#filename="parcels_000_0100.dat"

    print(filename)
    x, y, z,tag,n = read_filename(filename)



    nperdir=np.cbrt(n)
    #print("Number of parcels per direction is %d"%nperdir)

    #res= The number of parcels per direction we want to plot
    res=16.

    stride=int(n/(res*res*res))

    print("stride= %d"%stride)

    stride=1



    for i in range(n):
        if i%stride == 0:
            xs.append(x[i])
            ys.append(y[i])
            zs.append(z[i])
            tags.append(tag[i])

print(len(xs), np.cbrt(len(xs)))

fig = plt.figure()
ax = fig.add_subplot(111)#, projection='3d')
ax.scatter(xs, ys, c=tags,s=4)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')


plt.show()
