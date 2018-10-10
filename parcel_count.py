#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

f=open("parcels.dat","r")

data=[]
total=[]
n=0
for line in f:
    vals=line.split()
    i=0
    total.append(0.)
    if n == 0:
        for val in vals:
            data.append([])
            n+=1
    for val in vals:
        data[i].append(float(val))
        total[-1]+=float(val)
        i+=1
print("Total number of processes = %d"%n)
print("Total number of timesteps = %d"%(len(total)))

for i in range(n):
    plt.plot(data[i])

plt.show()

plt.plot(total)
plt.show()
