# File:     hadoop mapper file
# Function: 2-dimension N-body simulator, compute one time interval, using 
#           hadoop to compute and get result, the input file and output file 
#           is as the same format, so it can do loop simulate
# Implementation: 
#           Get the input infomation from init.file, computing vel change, which
#           used Newton's second law and Law of universal gravitation. In the 
#           mapper part, it computes the vel change in each map and send it to 
#           reducer part.

import sys
import math

G = 6.67e-11
time_interval = 0.1 # time interval

class Particle(object):
    def __init__(self, line):
        self.pos = [0, 0]
        self.vel = [0, 0]
        
        self.id, self.pos[0], self.pos[1], self.vel[0], self.vel[1], \
            self.m = [float(i) for i in line.strip().split()]
        self.id = int(self.id)
        
    def compute(self, q):
        if q.id == self.id:
            return [0, 0]
        
        # Compute the x and y distances and total distance d between
        # bodies self and q
        diff_x = self.pos[0] - q.pos[0]
        diff_y = self.pos[1] - q.pos[1]
        dist = math.sqrt(math.pow(diff_x, 2) + math.pow(diff_y, 2))
        
        if dist < 5:
            d = 5
            
        dist_cubed = math.pow(dist, 3)
        a_t =  G * q.m / dist_cubed # the mass of self was reduced
        
        return [a_t * diff_x, a_t * diff_y]


for line in sys.stdin:
    p = Particle(line)
    for qline in open('nbody/init.txt'):
        q = Particle(qline)
        acc = p.compute(q) # acc is the acceleration
        dvel = [i * time_interval for i in acc] # vel change of each particle
        
        print p.id, p.pos[0], p.pos[1], p.vel[0], p.vel[1], p.m, \
            dvel[0], dvel[1]

