#!/bin/python


from math import *


area = 1212

elist = [0.3, 0.6, 0.8, 0.9]

g = open("oblate_params.dat", "w")
g.write("e\ta\tc\n")
# Oblate:
for e in elist:
    a = area/(1 + (1-e*e)/e*atanh(e))/(2*pi)
    a = sqrt(a)
    c = a*sqrt(1-e*e)
    g.write(str(e)+"\t"+str(a)+"\t"+str(c)+"\n")

g.close()

# Prolate
g = open("prolate_params.dat", "w")
g.write("e\ta\tc\n")
for e in elist:
#    a = area/(1 + sqrt(1/(1-e*e))/e*asin(e))/(2*pi)
    a = area/(1 + asin(e)/sqrt(1-e*e)/e)/2/pi
    a = sqrt(a)
    c = a/sqrt(1-e*e)
    g.write(str(e)+"\t"+str(a)+"\t"+str(c)+"\n")

g.close()
    
g = open("sphere_params.dat", "w")
g.write("e\ta\tc\n")
r = area/4/pi
r = sqrt(r)
g.write("0\t"+str(r)+"\t"+str(r)+"\n")
g.close()
