#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys

def read_index():
    fname="timing_data/INDEX"
    f=open(fname,"r")
    routines=[]

    print("Reading '%s':"%fname)

    for line in f:
        if line[0]=="#":
            continue
        elif line[0:8]=="nprocs =":
            n=int(line[8:])
            print("Number of processes= %5d"%n)
        elif line[0:8]=="routine=":
            routines.append(line[8:-1])
            print("Routine: '%s'"%routines[-1])
        else:
            print("Error unknown option '%s'"%line[0:8])
            sys.exit()

    f.close()

    return (n,routines)


def read_file(fname):
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

    f.close()

    return (t,iteration,ttot,runningtot)


#get info on number of processes and routines traced
(n,routines) = read_index()

#create dictonary to hold data
#final dictionary structure will be data[routine name][process number][name of array]
data={}

print("\nReading files:")
for routine in routines:

    data[routine] = []

    for proc in range(n):
        fname="timing_data/"+routine+"_"+"%05d"%proc+".log"
        print(fname)
        timings={}
        (t,iteration,ttot,runningtot) = read_file(fname)
        timings["t"]=t
        timings["iteration"]=iteration
        timings["ttot"] = ttot
        timings["runningtot"] = runningtot
        data[routine].append(timings)

print("Done")

for routine in routines:
    x=data[routine][0]["ttot"]
    y=data[routine][0]["t"]
    plt.plot(x,y,marker="x",linestyle="",label=routine)

plt.title("Summary of rank 0's routines runtime")
plt.xlabel("runtime (s)")
plt.ylabel("Time in routine (s)")
plt.legend()
plt.show()

for routine in routines:
    x=data[routine][0]["ttot"]
    y=data[routine][0]["runningtot"]
    plt.plot(x,y,marker="x",linestyle="",label=routine)
plt.title("Summary of rank 0's routines cumulative runtime")
plt.xlabel("runtime (s)")
plt.ylabel("Cumulative time in routine (s)")
plt.legend()
plt.show()

for routine in routines:
    for proc in range(n):
        x=data[routine][proc]["ttot"]
        y=data[routine][proc]["t"]
        plt.plot(x,y,marker="x",linestyle="",label="%5d"%proc)

    plt.title(routine)
    plt.xlabel("runtime (s)")
    plt.ylabel("Time in routine (s)")
    plt.legend()
    plt.show()

# plt.plot(iteration,t,marker="x",linestyle="")
# plt.xlabel("Call number")
# plt.ylabel("Time in routine (s)")
# plt.show()
#
# plt.plot(ttot,t,marker="x",linestyle="")
# plt.xlabel("Runtime (s)")
# plt.ylabel("Time in routine (s)")
# plt.show()
#
# plt.plot(ttot,runningtot,marker="x",linestyle="")
# plt.xlabel("Runtime (s)")
# plt.ylabel("Cumulative time (s)")
# plt.show()
