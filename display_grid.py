#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

file="grids_000_0163.dat"

def readframe(num):
    file=str("grids_000_%04d.dat"%num)
    print(file)
    f=open(file,"rb")


    t=np.fromfile(f,dtype=np.float64,count=1)
    rg=np.fromfile(f,dtype=np.float64,count=6)
    nx=np.fromfile(f,dtype=np.int32,count=1)
    ny=np.fromfile(f,dtype=np.int32,count=1)
    nz=np.fromfile(f,dtype=np.int32,count=1)
    print("Time= %f"%t[0])

    nx=int(nx[0])
    ny=int(ny[0])
    nz=int(nz[0])

    #We read in eah variable. Note: in fortran they are stored [z,y,x] but as numpy is row major we
    # read them in as [x,y,z]
    u=np.fromfile(f,dtype=np.float64,count=nx*ny*nz)
    v=np.fromfile(f,dtype=np.float64,count=nx*ny*nz)
    w=np.fromfile(f,dtype=np.float64,count=nx*ny*nz)
    p=np.fromfile(f,dtype=np.float64,count=nx*ny*nz)
    q=np.fromfile(f,dtype=np.float64,count=nx*ny*nz)
    r=np.fromfile(f,dtype=np.float64,count=nx*ny*nz)
    b=np.fromfile(f,dtype=np.float64,count=nx*ny*nz)
    hg=np.fromfile(f,dtype=np.float64,count=nx*ny*nz)
    hgliq=np.fromfile(f,dtype=np.float64,count=nx*ny*nz)
    vol=np.fromfile(f,dtype=np.float64,count=nx*ny*nz)
    f.close()


    b=b.reshape((nx,ny,nz))
    hgliq=hgliq.reshape((nx,ny,nz))


    sliceb=b[nx/2,:,:]
    plt.imshow(sliceb.T,origin='lower',extent=[0,1,0,1])#,cmap="rainbow")
    plt.colorbar()
    plt.title(str("t= %2.2f"%t[0]))
    plt.savefig("b%03d.png"%num)
    plt.clf()

    slicehgliq=hgliq[nx/2,:,:]
    plt.imshow(slicehgliq.T,origin='lower',extent=[0,1,0,1],vmin=0.,vmax=0.08)#,cmap="rainbow")
    plt.colorbar()
    plt.title(str("t= %2.2f"%t[0]))
    plt.savefig("hgliq%03d.png"%num)
    plt.clf()


for i in range(228):
    readframe(i+1)
