#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import os

file="parcels_000_0100.dat"
gwidth=4
root="vort"

def read_file(fname):
    print("File: %s:"%fname)
    f=open(fname,"rb")
    t=np.fromfile(f,dtype=np.float64,count=1)
    t=t[0]
    xrange=np.fromfile(f,dtype=np.float64,count=2)
    yrange=np.fromfile(f,dtype=np.float64,count=2)
    zrange=np.fromfile(f,dtype=np.float64,count=2)
    n=np.fromfile(f,dtype=np.int64,count=1)
    n=n[0]

    x=np.fromfile(f,dtype=np.float64,count=n)
    y=np.fromfile(f,dtype=np.float64,count=n)
    z=np.fromfile(f,dtype=np.float64,count=n)

    p=np.fromfile(f,dtype=np.float64,count=n)
    q=np.fromfile(f,dtype=np.float64,count=n)
    r=np.fromfile(f,dtype=np.float64,count=n)

    dxdt=np.fromfile(f,dtype=np.float64,count=n)
    dydt=np.fromfile(f,dtype=np.float64,count=n)
    dzdt=np.fromfile(f,dtype=np.float64,count=n)

    dpdt=np.fromfile(f,dtype=np.float64,count=n)
    dqdt=np.fromfile(f,dtype=np.float64,count=n)
    drdt=np.fromfile(f,dtype=np.float64,count=n)

    h=np.fromfile(f,dtype=np.float64,count=n)
    b=np.fromfile(f,dtype=np.float64,count=n)
    vol=np.fromfile(f,dtype=np.float64,count=n)

    stretch=np.fromfile(f,dtype=np.float64,count=n)
    tag=np.fromfile(f,dtype=np.int32,count=n)




    # print("time= %f"%t)
    # print("xrange= %f to %f"%(xrange[0],xrange[1]))
    # print("yrange= %f to %f"%(yrange[0],yrange[1]))
    # print("zrange= %f to %f"%(zrange[0],zrange[1]))
    print("Num parcels= %d"%n)



    f.close()

    # f = np.sqrt(p**2 + q**2 + r**2)
    # hliq=np.zeros(n)
    # for i in range(len(hliq)):
    #     hliq[i] = max(0.,h[i] - np.exp(-z[i]))
    #
    # render_slice(x,y,z,b, vol,xrange,yrange,zrange,"x",(xrange[1]-xrange[0])/2,200,kernel="gaussian",number=number,root="buoyancy",time=t)
    # render_slice(x,y,z,hliq, vol,xrange,yrange,zrange,"x",(xrange[1]-xrange[0])/2,200,kernel="gaussian",number=number,root="hgliq",time=t)
    # render_slice(x,y,z,f, vol,xrange,yrange,zrange,"x",(xrange[1]-xrange[0])/2,200,kernel="gaussian",number=number,root="vort",time=t)

    return x,y,z,p,q,r,dxdt,dydt,dzdt,dpdt,dqdt,drdt,h,b,vol,stretch,tag


def render_slice(x,y,z,var,vol,xrange,yrange,zrange,plane,loc,resolution,kernel="gaussian",number=0, root=root, time=0.):
    dx = (xrange[1]-xrange[0])/(resolution)
    dy = (yrange[1]-yrange[0])/(resolution)
    dz = (zrange[1]-zrange[0])/(resolution)

    print('rendering %s slice at %s=%f using %s kernel'%(root,plane,loc,kernel))

    #create new xp, yp coordinates corresponding to x, y of image, and a z for depth.
    xp=[]
    yp=[]
    zp=[]
    volp=[]
    varp=[]

    if plane == "x":
        if kernel == "gaussian":
            xmin = loc-gwidth*dx
            xmax = loc+gwidth*dx
        else:
            xmin = loc-dx
            xmax = loc+dx
        dxp=dy
        dyp=dz
        dzp=dx
        for i in range(len(x)):
            if x[i] > xmin and x[i] < xmax:
                xp.append(y[i])
                yp.append(z[i])
                zp.append(x[i])
                volp.append(vol[i])
                varp.append(var[i])
        xg=np.arange(yrange[0],yrange[1]+dy,dy)
        yg=np.arange(zrange[0],zrange[1]+dz,dz)
    elif plane == "y":
        if kernel == "gaussian":
            ymin = loc-gwidth*dy
            ymax = loc+gwidth*dy
        else:
            ymin = loc-dy
            ymax = loc+dy
        dxp=dx
        dyp=dz
        dzp=dy
        for i in range(len(y)):
            if y[i] > ymin and y[i] < ymax:
                xp.append(x[i])
                yp.append(z[i])
                zp.append(y[i])
                volp.append(vol[i])
                varp.append(var[i])
            xg=np.arange(xrange[0],xrange[1]+dx,dx)
            yg=np.arange(zrange[0],zrange[1]+dz,dz)
    elif plane == "z":
        if kernel == "gaussian":
            zmin = loc-gwidth*dz
            zmax = loc+gwidth*dz
        else:
            zmin = loc-dz
            zmax = loc+dz
        dxp=dx
        dyp=dy
        dzp=dz
        for i in range(len(z)):
            if z[i] > zmin and z[i] < zmax:
                xp.append(x[i])
                yp.append(y[i])
                zp.append(z[i])
                volp.append(vol[i])
                varp.append(var[i])
        xg=np.arange(xrange[0],xrange[1]+dx,dx)
        yg=np.arange(yrange[0],yrange[1]+dy,dy)
    else:
        print("Invalid plane. Aborting")
        print(plane)
        sys.exit()

    print("Identified %d parcels in plane"%len(xp))

    img=np.zeros((resolution+1,resolution+1))
    wgt=np.zeros((resolution+1,resolution+1))

    for num in range(len(xp)):
        #identify the cell we're in (index of lower left corner)
        i=int(math.floor(xp[num]/dxp))
        j=int(math.floor(yp[num]/dyp))

        if kernel == "linear":

            for ii in range(i,i+2):
                for jj in range(j,j+2):
                    w = (1.-np.abs(xp[num]-xg[ii])/dxp)*(1.-np.abs(yp[num]-yg[jj])/dyp)*(1-np.abs(zp[num]-loc)/dzp)*volp[num]
                    img[ii,jj] += varp[num]*w
                    wgt[ii,jj] += w


        elif kernel == "gaussian":


            for ii in range(i-gwidth,i+gwidth):
                if ii >= 0 and ii <= resolution:
                    for jj in range(j-gwidth,j+gwidth):
                        if jj >= 0 and jj <= resolution:
                            w = np.exp( -(xp[num] - xg[ii])**2/dxp**2 - (yp[num] - yg[jj])**2/dyp**2 -(zp[num] - loc)**2/dzp**2) * volp[num]
                            img[ii,jj] += varp[num]*w
                            wgt[ii,jj] += w


        else:
            print("Error, invalid kernel '%s'. Aborting."%kernel)
            print(sys.abort())


    #Remove any zero weights
    counter=0
    for i in range(resolution+1):
        for j in range(resolution+1):
            if wgt[i,j] == 0.:
                wgt[i,j]=1.
                counter+=1
    print("Removed %d zeroed pixels"%counter)

    img=img/wgt

    img=img.T
    if root == "vort":
        plt.imshow(img,origin='lower',extent=[0,1,0,1],vmin=0.,vmax=max(10.,np.amax(img)))
    elif root == "hgliq":
        plt.imshow(img,origin='lower',extent=[0,1,0,1],vmin=0.,vmax=0.08)
    else:
        plt.imshow(img,origin='lower',extent=[0,1,0,1])
    plt.colorbar()
    plt.title("Time = %f"%time)

    #plt.savefig("%s_%04d.png"%(root,number))
    #plt.clf()
    plt.show()



if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: 'display_parcels.py [framenumber] [optional: number of files]'")
        sys.exit()

    filenum = int(sys.argv[1])
    if len(sys.argv) == 3:
        numfiles=int(sys.argv[2])
    else:
        numfiles=1

    #first check that the files exist
    fail=False
    for i in range(numfiles):
        fname = "parcels_%03d_%04d.dat"%(i,filenum)
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

    #now get the global x/y/z ranges
    for i in range(numfiles):
        fname = "parcels_%03d_%04d.dat"%(i,filenum)
        f=open(fname,"rb")
        t=np.fromfile(f,dtype=np.float64,count=1)
        t=t[0]
        xmin=float(np.fromfile(f,dtype=np.float64,count=1))
        xmax=float(np.fromfile(f,dtype=np.float64,count=1))
        ymin=float(np.fromfile(f,dtype=np.float64,count=1))
        ymax=float(np.fromfile(f,dtype=np.float64,count=1))
        zmin=float(np.fromfile(f,dtype=np.float64,count=1))
        zmax=float(np.fromfile(f,dtype=np.float64,count=1))
        if i==0:
            xrange.append(xmin)
            xrange.append(xmax)
            yrange.append(ymin)
            yrange.append(ymax)
            zrange.append(zmin)
            zrange.append(zmax)
        else:
            if xmin<xrange[0]: xrange[0]=xmin
            if xmax>xrange[1]: xrange[1]=xmax
            if ymin<yrange[0]: yrange[0]=ymin
            if ymax>yrange[1]: yrange[1]=ymax
            if zmin<zrange[0]: zrange[0]=zmin
            if zmax>zrange[1]: zrange[1]=zmax
        f.close()
    print("xrange=",xrange)
    print("yrange=",yrange)
    print("zrange=",zrange)

    #now read in the files
    x=np.zeros(0)
    y=np.zeros(0)
    z=np.zeros(0)
    p=np.zeros(0)
    q=np.zeros(0)
    r=np.zeros(0)
    dxdt=np.zeros(0)
    dydt=np.zeros(0)
    dzdt=np.zeros(0)
    dpdt=np.zeros(0)
    dqdt=np.zeros(0)
    drdt=np.zeros(0)
    h=np.zeros(0)
    b=np.zeros(0)
    vol=np.zeros(0)
    stretch=np.zeros(0)
    tag=np.zeros(0)

    for i in range(numfiles):
        fname = "parcels_%03d_%04d.dat"%(i,filenum)
        xp,yp,zp,pp,qp,rp,dxdtp,dydtp,dzdtp,dpdtp,dqdtp,drdtp,hp,bp,volp,sp,tp = read_file(fname)
        x=np.append(x,xp)
        y=np.append(y,yp)
        z=np.append(z,zp)
        p=np.append(p,pp)
        q=np.append(q,qp)
        r=np.append(r,rp)
        dxdt=np.append(dxdt,dxdtp)
        dydt=np.append(dydt,dydtp)
        dzdt=np.append(dzdt,dzdtp)
        dpdt=np.append(dpdt,dpdtp)
        dqdt=np.append(dqdt,dqdtp)
        drdt=np.append(drdt,drdtp)
        h=np.append(h,hp)
        b=np.append(b,bp)
        vol=np.append(vol,volp)
        stretch=np.append(stretch,sp)
        tag=np.append(tag,tp)
        print(len(x))

    print("read everything properly")
    render_slice(x,y,z,b, vol,xrange,yrange,zrange,"x",(xrange[1]-xrange[0])/2,200,kernel="gaussian",number=filenum,root="buoyancy",time=t)
