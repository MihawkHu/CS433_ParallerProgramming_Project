# File:     hadoop reducer file
# Function: 2-dimension N-body simulator, compute one time interval, using 
#           hadoop to compute and get result, the input file and output file 
#           is as the same format, so it can do loop simulate
# Implementation: 
#           Get the input infomation from mapper, updating pos and vel, which
#           used Newton's second law and Law of universal gravitation. In the 
#           reducer part, it updates pos and vel, and print updated information
#           to output

import sys

time_interval = 0.1 # time interval

# dist to record infomation of each particle
par_vel = {}

for line in sys.stdin:
    # get data from mapper 
    # the vel change of each particle
    dvel = [0, 0]
    vel = [0, 0]
    pos = [0, 0]
    
    id, pos[0], pos[1], vel[0], vel[1], m, dvel[0], dvel[1] = \
        [float(i) for i in line.strip().split()]
    id = int(id)
    
    # update vel 
    p = par_vel.get(id, [pos[0], pos[1], vel[0], vel[1], m])
    p[2] = p[2] + dvel[0]
    p[3] = p[3] + dvel[1]
    par_vel[id] = p

for id, p in par_vel.items():
    # update pos 
    p[0] = p[0] + p[2] * time_interval
    p[1] = p[1] + p[3] * time_interval
    
    print id, p[0], p[1], p[2], p[3], p[4] # output updated infomation

