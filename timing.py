#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys


if len(sys.argv) != 2:
    print("Usage: ./timing.py [filename]")
    sys.exit()

fname=sys.argv[1]

f=open(fname,"r")

t=[]
iteration=[]
ttot=[]
runningtot=[]
n=0
dt=0.

for line in f:
    if line[0]=="#":
        continue

    if line[0:14] == "TIMER_START  ,":
        tstart=float(line[15:])
    elif line[0:14] == "TIMER_PAUSE  ,":
        dt+= float(line[15:])-tstart
    elif line[0:14] == "TIMER_RESUME ,":
        tstart=float(line[15:])
    elif line[0:14] == "TIMER_STOP   ,":
        tstop=float(line[15:])

        time=dt + tstop-tstart
        t.append(time)
        iteration.append(n)
        ttot.append(tstop)
        if n==0:
            runningtot.append(tstop)
        else:
            runningtot.append(runningtot[-1]+time)
        n+=1
        dt=0.
    else:
        print("Error: unknown event 's'"%line[0:14])


plt.plot(iteration,t)
plt.xlabel("Call number")
plt.ylabel("Time in routine (s)")
plt.show()

plt.plot(ttot,t)
plt.xlabel("Runtime (s)")
plt.ylabel("Time in routine (s)")
plt.show()

plt.plot(ttot,runningtot)
plt.xlabel("Runtime (s)")
plt.ylabel("Cumulative time (s)")
plt.show()

f.close()
