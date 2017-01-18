# File:     generate.py
# Function: generate random initial state of particles, 
#           used for 2-dimension N-body simulator
# Particle Info:
#           index, pos_x, pos_y, vel_x, vel_y, mass
#           Particle number is N, the default is 200, index from 0 to N-1, 
#           initial pos from 0 to 400, initial vel from -10 to 10, mass from 
#           1 to 100

from random import uniform

N = 200

for a in xrange(N):
    print a, uniform(0, 400), uniform(0, 400), uniform(-10, 10), \
        uniform(-10, 10), uniform(1, 100)
