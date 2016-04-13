from __future__ import division
import numpy as np
from math import *
from visual import *

scene.background = (1,1,1)

M = 1 #This doesn't matter, included for completeness
PI = 3.14159265359
G = 6.67408e-11 #Grav. constant in kgs
p_initial = .0035 * 6.767991e-20 #kg per m^3
R_0 = 5000 * 3.086e16 #this is used in the potential function
r_initial = 6000 * 3.086e16 #m
r_max = 50000 * 3.086e16 #m
psi_initial = 0
V_r_initial = -43500 #m/s
V_psi_initial = 100000 #m/s

def cart2pol(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return(r, theta)

def pol2cart(r, theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return(x, y)

def potential(r_in):
    return -4*PI*G*p_initial*(R_0**2)*(1 + math.log((r_max / r_in)))

#the center of galaxy
center = sphere(pos = (0,0), color = color.magenta, radius = 5e18)
#the star
orbiter = sphere(pos = pol2cart(r_initial, psi_initial), color = color.black, make_trail = True)

E = potential(r_initial) + .5*(V_r_initial**2) + .5*(V_psi_initial**2) #energy is a constant
L = r_initial*V_psi_initial # angular momentum is a constant

t = 0
dt = 5e5 * 3.154e7 # first number is in yrs
t_final = 2e10 * 3.154e7 #first number is in yrs

r = r_initial
psi = psi_initial
V_r = V_r_initial
V_psi = V_psi_initial

r_new = 0
psi_new = 0
V_r_new = 0
V_psi_new = 0

V_psi_max = 0
r_at_V_psi_max = 0

V_r_min = np.abs(V_r_initial)
r_at_V_r_min = 0

t_2pi = 0
t_peri = 0
peri_or_apo = 0 #if we're at a aop-galacticon, this is 0, if at peri, it's one.
sum_of_peri_times = 0
num_of_peri_times = -1 #starts at -1 so we can throw out the first number, not a full cycle
while t < t_final: 
    rate(10000)
    r_new = r + V_r*dt
    psi_new = psi + math.atan(dt*V_psi/r_new)
    V_psi_new = L / r_new  

    #we want to keep the sign of V_r
    if (2*E - 2*potential(r_new) - V_psi_new**2) < 0:
        V_r_new = copysign(np.sqrt(np.abs(2*E - 2*potential(r_new) - V_psi_new**2)), V_r_new*(-1))
        if (peri_or_apo == 1):
            #this if throws out the first number, not a full cycle yet
            if (num_of_peri_times == -1):
                num_of_peri_times += 1
            #so we can get the average peri to peri time
            else:
                sum_of_peri_times += t_peri
                num_of_peri_times += 1
                t_peri = 0
                peri_or_apo -= 1
        else:
            peri_or_apo += 1
    else:
        V_r_new = copysign(np.sqrt(np.abs(2*E - 2*potential(r_new) - V_psi_new**2)), V_r_new)
        t_peri = t_peri + dt
    #update    
    r = r_new
    psi = psi_new
    V_psi = V_psi_new
    V_r = V_r_new

    #for finding the min V_r, where the max r is.
    if (np.abs(V_r) < V_r_min):
        V_r_min = V_r
        r_at_V_r_min = r
    #for finding the max V_psi, where the min r is.
    if (V_psi > V_psi_max):
        V_psi_max = V_psi
        r_at_V_psi_max = r
    #how long to go thru 2pi, to return to starting psi.
    if ((2*PI)-.001 < psi < (2*PI)+.001):
        t_2pi = t / 3.154e7 # to get back to yrs
                
    orbiter.pos = pol2cart(r, psi)

    t = t + dt
print "a)"
print "    r min (pc) = "
print r_at_V_psi_max / 3.086e16
print "    r max (pc) = "
print r_at_V_r_min / 3.086e16
print ""
print "b)"
print "    time to travel through 2pi (years) ="
print t_2pi
print "    time from one peri to next (years) ="
print sum_of_peri_times / num_of_peri_times/ 3.154e7 
