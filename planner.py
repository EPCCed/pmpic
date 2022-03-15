#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np

#This code calculates the aproximate memory footprint of PMPIC (total and per process)
#It can thus be used to choose a value of maxparcels to maximise memory usage

#number of grid/parcel variables
n_grids=32
n_parcels=38

#number of bytes in a variable
sizeof=8

#default values (can be overridden)
ppn=128 #processes per node
mpn=256 #memory per node (GB)

#parse an input value of the form "quantity=number" to return the number
def getvalue(arg):
    loc=arg.find("=")
    if loc == -1:
        print("Warning: Cannot parse arguent '%s'"%arg)
        return None

    num=arg[loc+1:]
    try:
        num=int(num)
        return num
    except:
        print("Warning: cound not convert %s value '%s' to a string"%(arg[0:loc],num))
        return None


def getdims(nprocs):
    #first we want to factorise nprocs
    n=nprocs
    factors=[]
    i=2
    while n>=i:
        if n%i==0:
            factors.append(i)
            n /= i
        else:
            i+=1

    #print("Factors are", factors)

    dims=[1,1]

    while len(factors) > 0:
        val = factors.pop()
        dims.sort()
        dims[0] *= val

    dims.sort()

    return dims







#command line arguments:
#nprocs - number of processes
nprocs=None
#ngrid - number of grid cells (nx=ny=nz=ngrid)
ngrid=None
#nx - number of cells in x
nx=None
#ny - number of cells in y
ny=None
#nz - number of cells in z
nz=None
#ppn - processes per node
#mpn - memory per node


#get the command line arguments and remove the first one (script name)
args=sys.argv
args.pop(0)

if len(args) == 0:
    print('Error: no command line arguments supplied')
    print("\nUsage: 'planner.py [args]'")
    print("\nPossible arguments:")
    print("nprocs - number of processes used (REQUIRED)")
    print("ngrid  - number of grid cells in each direction*")
    print("nx     - number of cells in x direction*")
    print("ny     - number of cells in y direction*")
    print("nz     - number of cells in z direction*")
    print("ppn    - number of processes per node (default=24)")
    print("mpn    - memory per node (default 64GB)")
    print("* either n or nx, ny and nx must be specified")
    print("\nExample usage: 'planner.py nprocs=4 ngrid=128'")
    sys.exit()

#parse arguments:
print("\nThe following options were registered:")
for arg in args:
    if "nprocs=" in arg:
        print(arg)
        nprocs=getvalue(arg)
        args.remove(arg)

for arg in args:
    if "ngrid=" in arg:
        print(arg)
        ngrid=getvalue(arg)
        nx=ngrid
        ny=ngrid
        nz=ngrid
        args.remove(arg)

for arg in args:
    if "nx=" in arg:
        print(arg)
        nx=getvalue(arg)
        args.remove(arg)

for arg in args:
    if "ny=" in arg:
        print(arg)
        ny=getvalue(arg)
        args.remove(arg)

for arg in args:
    if "nz=" in arg:
        print(arg)
        nz=getvalue(arg)
        args.remove(arg)

for arg in args:
    if "ppn=" in arg:
        print(arg)
        ppn=getvalue(arg)
        args.remove(arg)

for arg in args:
    if "mpn=" in arg:
        print(arg)
        mpn=getvalue(arg)
        args.remove(arg)

if len(args) > 0:
    print("\nWarning: There are un-parsed arguments:")
    for arg in args:
        print("'%s'"%arg)

if nprocs==None:
    print('Error: nprocs must defined')
    sys.exit()

if ngrid == None:
    if nx == None or ny==None or nz==None:
        print("Error: either ngrid or nx, ny and nz must be defined")
        sys.exit()


#determine nodes used:
nnodes=int(np.ceil(float(nprocs)/ppn))
#determine the decomposition of the grid between processes (should return the same result as MPI_Dims_create)
dims = getdims(nprocs)

print("\n################################################################################")
#print("01234567890123456789012345678901234567890123456789012345678901234567890123456789")
print("#                            Parallel Decomposition                            #")
print("################################################################################")
print("#                                                                              #")
print("#   Number of Processes         = %5d                                        #"%nprocs)
print("#   Number of nodes             = %4d                                         #"%nnodes)
print("#   Processes in x, y and z     = %3d x %3d x %3d                              #"%(dims[0],dims[1],1))
#print("################################################################################")
#print("\n################################################################################")
#print("01234567890123456789012345678901234567890123456789012345678901234567890123456789")
#print("#                                  Grids                                       #")
#print("################################################################################")

if nx%dims[0] !=0 or ny%dims[1] !=0:
    print("Error: MPI cannot divide the grids evenly between processes!")
    sys.exit()
nxlocal = nx/dims[0]
nylocal = ny/dims[1]
print("#   Global gridsizes (nx,ny,nz) = (%4d,%4d,%4d)                             #"%(nx,ny,nz))
print("#   Local gridsizes             = (%3d,%3d,%3d)                                #"%(nxlocal,nylocal,nz))
print("#   Local gridsizes (inc halos) = (%3d,%3d,%3d)                                #"%(nxlocal+4,nylocal+4,nz+1))

print("#                                                                              #")

print("################################################################################")
#print("01234567890123456789012345678901234567890123456789012345678901234567890123456789")
print("#                                 Memory                                       #")
print("################################################################################")
print("#                                                                              #")
print("#   Memory per node    = %6.2f GB                                             #"%mpn)
if nprocs<=ppn:
    mperproc=mpn/float(nprocs)
else:
    mperproc = mpn/float(ppn)
print("#   Memory per process = %5.2f GB                                              #"%mperproc)
print("#                                                                              #")

print("#   Assuming 10% of the memory is required by PMPIC and temporary variables:   #")
mperproc *= 0.9
print("#   Memory for parcels and grids = %5.2f GB                                    #"%mperproc)

gridsize=(nxlocal+4)*(nylocal+4)*(nz+1)*sizeof
grids=gridsize*n_grids
print("#                                                                              #")
print("#   Grids:                                                                     #")
print("#   Number of gridded variables = %2d                                           #"%n_grids)
print("#   A single gridded variable takes up %7.2f MB                              #"%(gridsize/1000000.))
print("#   All gridded variables take up %7.2f MB                                   #"%(grids/1000000.))

# we now calculate OMP_Stacksize requirements
stacksize=int(np.ceil(gridsize*14./1E6))
print("#   Recommended minimum OMP_STACKSIZE =%4dM                                   #"%(stacksize))


if (grids/1E9 > mperproc) :
    print('Error: insufficient memory for grids')
    sys.exit()

mperproc = mperproc - grids/1E9


numparcels = mperproc*1E9/float(n_parcels)/float(sizeof)
print("#                                                                              #")
print("#   Parcels:                                                                   #")
print("#   Available memory for parcels = %5.2f GB                                    #"%(mperproc))
print("#   Number of parcel variables = %2d                                            #"%n_parcels)
print("#   Therefore a single parcel takes up %3d Bytes                               #"%(sizeof*n_parcels))
print("#   Therefore each process can have up to %10d parcels                   #"%(numparcels))
print("#   Maximum parcel count for simulation = %12d                         #"%(numparcels*nprocs))

totalparcels=numparcels*nprocs
npercell = totalparcels/nx/ny/nz
print("#   Maximum mean number of parcels per cell = %5.1f  (%6.2f^3)                #"%(npercell,np.power(npercell,1./3.)))
print("#                                                                              #")
print("################################################################################\n")
