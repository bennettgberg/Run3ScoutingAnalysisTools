import os, sys
import multiprocessing as mp
from math import ceil
from time import sleep

#function to do the hadding from a command
def hadder(comm):
    print(comm)

    os.system(comm)
    sleep(1)
    return

outname = sys.argv[1]
#total number of files that need to be hadded
ntot = int(sys.argv[2]) 

#what file number to start with
startn = 0
if len(sys.argv) > 3:
    startn = int(sys.argv[3])

#number of processors we should use
nproc = 6
#how many files per hadd call
nfile = 10

nloops = int(ceil(ntot / nfile))
#how many loops must be gone through
commands = []
for mainnum in range(nloops):
    cmd = "hadd -f hists_" + outname + "_mainnum" + str(mainnum) + ".root  "

    #have to do + 1 bc stupid files don't start at 0
    startnum = mainnum*nfile + startn
    endnum = startnum + nfile + 1
    nreal = 0
    for num in range(startnum, endnum):
        fname = "hists_" + outname + "_" + str(num) + ".root"
        if os.path.exists(fname):
            cmd += fname
            cmd += " "
            nreal += 1

    if nreal > 0:
        commands.append(cmd)

if len(commands) > 0:
    pool  = mp.Pool(nproc)
    pool.map(hadder, commands)
    pool.close()
    pool.join()

##now we have a bunch of mainnum files, and need to hadd those together
#newnloops = int(ceil(nloops / nfile))
#commands = []
#for mainnum in range(newnloops):
#    cmd = "hadd hists_" + outname + "_finalnum" + str(mainnum) + ".root  "
#    startnum = mainnum*nfile
#    endnum = startnum + nfile 
#    nreal = 0
#    for num in range(startnum, endnum):
#        fname = "hists_" + outname + "_mainnum" + str(num) + ".root"
#        if os.path.exists(fname):
#            cmd += fname
#            cmd += " "
#            nreal += 1
#
#    if nreal > 0:
#        commands.append(cmd)
#
#if len(commands) > 0:
#    pool  = mp.Pool(nproc)
#    pool.map(hadder, commands)
#    pool.close()
#    pool.join()

#now just hadd all the 'final' files, shouldn't be that many
#cmd = "hadd -f hists_" + outname + ".root hists_" + outname + "_finalnum*.root" 
cmd = "hadd -f hists_" + outname + ".root hists_" + outname + "_mainnum*.root" 
print(cmd)

os.system(cmd)
