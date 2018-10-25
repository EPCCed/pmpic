#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

file="grids_00000_00163.dat"

def readfile(fname):
    print("reading '%s'"%fname)
    f=open(fname,"rb")


    t=np.fromfile(f,dtype=np.float64,count=1)
    rg=np.fromfile(f,dtype=np.float64,count=6)
    nx=np.fromfile(f,dtype=np.int32,count=1)
    ny=np.fromfile(f,dtype=np.int32,count=1)
    nz=np.fromfile(f,dtype=np.int32,count=1)
    print("Time= %f"%t[0])

    nx=int(nx[0])
    ny=int(ny[0])
    nz=int(nz[0])

    #We read in each variable. Note: in fortran they are stored [z,y,x] but as numpy is row major we
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

    u=u.reshape((nx,ny,nz))
    v=v.reshape((nx,ny,nz))
    w=w.reshape((nx,ny,nz))
    p=p.reshape((nx,ny,nz))
    q=q.reshape((nx,ny,nz))
    r=r.reshape((nx,ny,nz))
    b=b.reshape((nx,ny,nz))
    hg=hg.reshape((nx,ny,nz))
    hgliq=hgliq.reshape((nx,ny,nz))


    # b=b.reshape((nx,ny,nz))
    # hgliq=hgliq.reshape((nx,ny,nz))
    #
    #
    # sliceb=b[nx/2,:,:]
    # plt.imshow(sliceb.T,origin='lower',extent=[0,1,0,1])#,cmap="rainbow")
    # plt.colorbar()
    # plt.title(str("t= %2.2f"%t[0]))
    # plt.savefig("b%03d.png"%num)
    # plt.clf()
    #
    # slicehgliq=hgliq[nx/2,:,:]
    # plt.imshow(slicehgliq.T,origin='lower',extent=[0,1,0,1],vmin=0.,vmax=0.08)#,cmap="rainbow")
    # plt.colorbar()
    # plt.title(str("t= %2.2f"%t[0]))
    # plt.savefig("hgliq%03d.png"%num)
    # plt.clf()

    return u, v, w, p, q, r, b,hg, hgliq


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python display_grid.py [frame number] [number of processes]")

    xcut = np.pi

    fnumber = int(sys.argv[1])
    nprocs=int(sys.argv[2])

    #first check that the files exist
    fail=False
    for i in range(nprocs):
        fname = "grids_%05d_%05d.dat"%(i,fnumber)
        if not os.path.isfile(fname):
            print("Error: cannot find '%s'"%fname)
            fail=True
        else:
            print("Found '%s'"%fname)
    if fail:
        print("Aborting")
        sys.exit()

    xrange=[]
    yrange=[]
    zrange=[]

    nx=[]
    ny=[]
    nz=[]
    xmin=[]
    ymin=[]
    zmin=[]
    xmax=[]
    ymax=[]
    zmax=[]

    #now get the global x/y/z ranges
    for i in range(nprocs):
        fname = "grids_%05d_%05d.dat"%(i,fnumber)
        f=open(fname,"rb")
        t=np.fromfile(f,dtype=np.float64,count=1)
        t=t[0]
        xmn=float(np.fromfile(f,dtype=np.float64,count=1))
        xmx=float(np.fromfile(f,dtype=np.float64,count=1))
        ymn=float(np.fromfile(f,dtype=np.float64,count=1))
        ymx=float(np.fromfile(f,dtype=np.float64,count=1))
        zmn=float(np.fromfile(f,dtype=np.float64,count=1))
        zmx=float(np.fromfile(f,dtype=np.float64,count=1))
        nx.append(float(np.fromfile(f,dtype=np.int32,count=1)))
        ny.append(float(np.fromfile(f,dtype=np.int32,count=1)))
        nz.append(float(np.fromfile(f,dtype=np.int32,count=1)))
        xmin.append(xmn)
        ymin.append(ymn)
        zmin.append(zmn)
        xmax.append(xmx)
        ymax.append(ymx)
        zmax.append(zmx)

        if i==0:
            xrange.append(xmn)
            xrange.append(xmx)
            yrange.append(ymn)
            yrange.append(ymx)
            zrange.append(zmn)
            zrange.append(zmx)
        else:
            if xmn<xrange[0]: xrange[0]=xmn
            if xmx>xrange[1]: xrange[1]=xmx
            if ymn<yrange[0]: yrange[0]=ymn
            if ymx>yrange[1]: yrange[1]=ymx
            if zmn<zrange[0]: zrange[0]=zmn
            if zmx>zrange[1]: zrange[1]=zmx
        f.close()
    # print("Global xrange=",xrange)
    # print("Global yrange=",yrange)
    # print("Global zrange=",zrange)
    # print("Local nx=",nx)
    # print("Local ny=",ny)
    # print("Local nz=",nz)
    # print("Local xmin",xmin)
    # print("Local ymin",ymin)
    # print("Local zmin",zmin)
    # print("Local xmax",xmax)
    # print("Local ymax",ymax)
    # print("Local zmax",zmax)

    #figure out the global size. We assume that the grid spacing is constant
    dx = (xmax[0]-xmin[0])/(nx[0]-1)
    dy = (ymax[0]-ymin[0])/(ny[0]-1)
    dz = (zmax[0]-zmin[0])/(nz[0]-1)

    print("dx, dy, dz = ", dx, dy, dz)

    nxg = int(round((xrange[1]-xrange[0])/dx +1))
    nyg = int(round((yrange[1]-yrange[0])/dy +1))
    nzg = int(round((zrange[1]-zrange[0])/dz +1))

    print("nxg, nyg, nzg = ", nxg, nyg, nzg)

    xc = np.arange(xrange[0],xrange[1]+dx,dx)
    yc = np.arange(yrange[0],yrange[1]+dy,dy)
    zc = np.arange(zrange[0],zrange[1]+dz,dz)

    #now determine the global start and end indices of each file
    irange=[]
    jrange=[]
    krange=[]

    for p in range(nprocs):
        istart = int(round((xmin[p]-xrange[0])/dx))
        istop = int(round((xmax[p]-xrange[0])/dx))
        irange.append([istart,istop])

        jstart = int(round((ymin[p]-yrange[0])/dy))
        jstop = int(round((ymax[p]-yrange[0])/dy))
        jrange.append([jstart,jstop])

        kstart = int(round((zmin[p]-zrange[0])/dz))
        kstop = int(round((zmax[p]-zrange[0])/dz))
        krange.append([kstart,kstop])

    print("irange = ",irange)
    print("jrange = ",jrange)
    print("krange = ",krange)

    img = np.zeros((nyg,nzg))

    # We now loop through each file
    for i in range(nprocs):

        #determine if we need this file
        if (xcut >= xmin[i] and xcut < xmax[i]):
            fname = "grids_%05d_%05d.dat"%(i,fnumber)
            u,v,w,p,q,r,b,hg,hgliq = readfile(fname)

            #determine x_coord of cut
            xindex = int(round((xcut-xrange[0])/dx))
            xindex = xindex-irange[i][0]

            img[jrange[i][0]:jrange[i][1]+1][krange[i][0]:krange[i][1]+1] = b[xindex][:][:]

    plt.imshow(img.T,origin='lower',extent=[0,1,0,1])
    plt.title("Time = %f"%t)
    plt.xlabel("y/L_y")
    plt.ylabel("z/L_z")
    plt.colorbar()
    plt.show()
