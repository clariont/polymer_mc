#!/usr/bin/python

############################################################################################################
##############################
##  Date: 16 Sept 2012
##  Description: Reads in a file with points on a sphere (specifying the location of the ligands), then maps it onto an ellipse.
##		 A LAMMPS .dat file is made for running MD on the points.
##############################
############################################################################################################

#import clarionmath
import math
import random
from os import sys
from decimal import *


#System Parameters:
nlig1 = 10
nlig2 = 10
lig1len = 4
lig2len = 8
xlo = -50
xhi = 50
ylo = -50
yhi = 50
zlo = -50
zhi = 50
rsphere = 1 
elp_a = float(sys.argv[1])
elp_c = float(sys.argv[2])
elp_b = elp_a;


# Main Program
# Read in points data
readname="ptsonsphere.dat"
f=open(readname, "r")
lines=f.readlines()
print "length ", len(lines)
allpts=[]
for line in lines:
    line=line.split()
    myX=float(line[0])
#    a = Decimal(myX)
    myX = myX*1.0
    myY=float(line[1])
    myZ=float(line[2])
#    print "pos: ",myX, ", ",myY, ", ", myZ
    pos=[myX,myY,myZ]
    allpts.append(pos)

# Rescale Points to the real sphere
#polepts=[]
#equatorpts=[]
#poledef = rsphere/2.0
#print "poledef is: ", poledef
#xgreater = 0
#ygreater = 0
#zgreater = 0
#zless = 0
#rotatefactor=0.105
#myvec = [0,rotatefactor,math.sqrt(rsphere*rsphere/4.0 - rotatefactor*rotatefactor)]
elppts = []
for points in allpts:
    myX = points[0]*rsphere*elp_a
    myY = points[1]*rsphere*elp_b
    myZ = points[2]*rsphere*elp_c
    pos=[myX,myY,myZ]
    elppts.append(pos)
#    print pos, "\t", poledef
#    dp1 = myX*myvec[0] + myY*myvec[1] + myZ*myvec[2]
#    if (dp1 > rsphere*rsphere/4.0 or dp1 < (-rsphere*rsphere/4.0)):
#	polepts.append(pos)
#    else:
#	equatorpts.append(pos)
#    if (myX > 0):
#	xgreater = xgreater + 1
#    if (myY > 0):
#	ygreater = ygreater + 1
#    if (myZ > 0):
#	zgreater = zgreater + 1
#    if (myZ < 0):
#	zless = zless + 1
#print "polepts len: ", len(polepts)
#print "equatorpts len: ", len(equatorpts)
#print "xgreater: ", xgreater
#print "ygreater: ", ygreater
#print "zgreater: ", zgreater
#print "zless: ", zless




#nlig1 = len(allpts)/2
#nlig2 = len(allpts)/2

# Write LAMMPS configuration file


atoms = 2*len(elppts)
bonds = len(elppts) 
#if (lig1len < 2 and lig2len < 2):
#    angles = 0
#elif (lig1len < 2):
#    angles = nlig2*(lig2len-2)
#elif (lig2len < 2):
#    angles = nlig1*(lig1len-2)
#else:
#    angles = nlig1*(lig1len-2) + nlig2*(lig2len-2)
angles = 0
dihedrals = 0
impropers = 0
#if (nlig2 == 0):
#    atomtypes = 3
#else:
#    atomtypes = 5
atomtypes = 2
bondtypes = 1
#if (lig1len > 2 or lig2len > 2):
#	angletypes = 1
#else:
#	angletypes = 0
angletypes = 0
dihedraltypes = 0
impropertypes = 0

g=open("elppts.dat", "w")
g.write("LAMMPS Description\n\n")
g.write("\t"+str(atoms)+"\tatoms\n")
g.write("\t"+str(bonds)+"\tbonds\n")
g.write("\t"+str(angles)+"\tangles\n")
g.write("\t"+str(dihedrals)+"\tdihedrals\n")
g.write("\t"+str(impropers)+"\timpropers\n\n")
g.write("\t"+str(atomtypes)+"\tatom types\n")
g.write("\t"+str(bondtypes)+"\tbond types\n")
g.write("\t"+str(angletypes)+"\tangle types\n")
g.write("\t"+str(dihedraltypes)+"\tdihedral types\n")
g.write("\t"+str(impropertypes)+"\timproper types\n\n\n")

#Write Box Lengths:
g.write(str(xlo)+" "+str(xhi)+" xlo xhi\n")
g.write(str(ylo)+" "+str(yhi)+" ylo yhi\n")
g.write(str(zlo)+" "+str(zhi)+" zlo zhi\n\n\n")

#Write Atomic Masses:
# Legend: 1 - fixed atom at origin.  2 - ligand1 head. 3 - ligand2 head. 4 - ligand tail.
g.write("Masses\n\n")
g.write("\t1 1.0\n")
g.write("\t2 1.0\n")
#if (nlig2 > 0):
#    g.write("\t3 1.0\n")
#    if (nlig1 > 1 or nlig2 > 1):
#	g.write("\t4 1.0\n")
#	g.write("\t5 1.0\n")
#if (nlig1 > 1 and nlig2 < 0):
#    g.write("\t4 1.0\n")
#    g.write("\t5 1.0\n")
#
#Write Atoms:
## Ligands are constrained to move on a sphere by (stiffly) bonding them to a fixed atom at the origin.
g.write("\nAtoms\n\n")
num = 1
for i in range(len(elppts)):
    g.write("\t"+str(num)+" "+str(num)+" 1 0 ")
    g.write("0.0 0.0 0.0 0 0 0\n")
    num = num + 1

for i in range(len(elppts)):
    g.write("\t"+str(num)+" "+str(num)+" 2 0 ")
    pos = elppts[i]
    g.write(str(pos[0])+" "+str(pos[1])+" "+str(pos[2])+" 0 0 0\n")
    num = num + 1

#
#num = 1   # counts the total number of atoms
#lig1pos = []
#lig2pos = []
#allligs = []
#overlap = 0
#
## Ligand 1:
## A ligand 1 has head atom "2" and tail atoms "4".
#for i in range(nlig1):
#    g.write("\t"+str(num)+" "+str(num)+" 1 0 ")
#    g.write("0.0 0.0 0.0 0 0 0\n")
#    num = num + 1
#
##    pos = randPtEllipse(rsphere, elp_a, elp_b, elp_c) 
##   Ligand 1 is the pole position.
#    pos = polepts[i]
##    if ((pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]) != rsphere*rsphere):
##	print "pointsq: ", (pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2])
#    allligs.append(pos)
#    g.write("\t"+str(num)+" "+str(num)+" 2 0 ")
#    g.write(str(pos[0])+" "+str(pos[1])+" "+str(pos[2])+" 0 0 0\n")
#    num = num + 1
#    for j in range(lig1len-1):
#	g.write("\t"+str(num)+" "+str(num)+" 4 0 ")
#	magnitude = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]
#	magnitude = math.sqrt(magnitude)
#	g.write(str(pos[0]*(1+(j+1)/magnitude))+" "+str(pos[1]*(1+(j+1)/magnitude))+" "+str(pos[2]*(1+(j+1)/magnitude))+" 0 0 0\n")
#	num = num + 1
#   
## Ligand 2 
## A ligand 2 atom has head atom "3" and tail atoms "5".
#for i in range(nlig2):
#    g.write("\t"+str(num)+" "+str(num)+" 1 0 ")
#    g.write("0.0 0.0 0.0 0 0 0\n")
#    num = num + 1
#
##    pos = randPtEllipse(rsphere, elp_a, elp_b, elp_c) 
##   Ligand 2 is the equator position.
#    pos = equatorpts[i]
#    overlap = 0
#    allligs.append(pos)
#    g.write("\t"+str(num)+" "+str(num)+" 3 0 ")
#    g.write(str(pos[0])+" "+str(pos[1])+" "+str(pos[2])+" 0 0 0\n")
#    num = num + 1
#    for j in range(lig2len-1):
#	g.write("\t"+str(num)+" "+str(num)+" 5 0 ")
#	magnitude = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]
#	magnitude = math.sqrt(magnitude)
#	g.write(str(pos[0]*(1+(j+1)/magnitude))+" "+str(pos[1]*(1+(j+1)/magnitude))+" "+str(pos[2]*(1+(j+1)/magnitude))+" 0 0 0\n")
#	num = num + 1
#
#Write Bonds:
g.write("\nBonds\n\n")
bondnum = 1
for i in range(len(elppts)):
    g.write("\t"+str(bondnum)+" 1 "+str(i+1)+" "+str(i+1+bonds)+"\n")
    bondnum = bondnum + 1
# Ligand 1:
#for i in range(nlig1):
#    g.write("\t"+str(bondnum)+" 1 "+str(i*(lig1len+1)+1)+" "+str(i*(lig1len+1)+2)+"\n")
#    bondnum = bondnum + 1  
#    # Bind head to first tail atom:
#    if (lig1len > 1):
#	g.write("\t"+str(bondnum)+" 2 " +str(i*(lig1len+1)+2)+" "+str(i*(lig1len+1)+3)+"\n")
#	bondnum = bondnum + 1
#	    # Bind remaining tail atoms:
#	for j in range(lig1len-2):
#		g.write("\t"+str(bondnum)+" 2 " +str(i*(lig1len+1)+j+3)+" "+str(i*(lig1len+1)+j+4)+"\n")
#		bondnum = bondnum + 1
#
## Ligand 2 
#if (lig2len > 1):
#	for i in range(nlig2):
#	    g.write("\t"+str(bondnum)+" 1 "+str(nlig1*(lig1len+1)+i*(lig2len+1)+1)+" "+str(nlig1*(lig1len+1)+i*(lig2len+1)+2)+"\n")
#	    bondnum = bondnum + 1  
#	    # Bind head to first tail atom:
#	    if (lig2len > 1):
#		g.write("\t"+str(bondnum)+" 2 " +str(nlig1*(lig1len+1)+i*(lig2len+1)+2)+" "+str(nlig1*(lig1len+1)+i*(lig2len+1)+3)+"\n")
#		bondnum = bondnum + 1
#	    # Bind remaining tail atoms:
#		for j in range(lig2len-2):
#			g.write("\t"+str(bondnum)+" 2 " +str(nlig1*(lig1len+1)+i*(lig2len+1)+j+3)+" "+str(nlig1*(lig1len+1)+i*(lig2len+1)+j+4)+"\n")
#			bondnum = bondnum + 1
#
#
##Write Angles:
#anglenum = 1
#if (lig1len > 2 or lig2len > 2):
#    g.write("\nAngles\n\n")
## Ligand 1:
#    if (lig1len > 2):
#	for i in range(nlig1):
#	    for j in range(lig1len-2):
#		g.write("\t"+str(anglenum)+" 1 "+str(i*(lig1len+1)+2+j)+" "+str(i*(lig1len+1)+3+j)+" "+str(i*(lig1len+1)+4+j)+"\n")
#		anglenum = anglenum + 1
## Ligand 2:
#    if (lig2len > 2):
#	for i in range(nlig2):
#	    for j in range(lig2len-2):
#		g.write("\t"+str(anglenum)+" 1 "+str(nlig1*(lig1len+1)+i*(lig2len+1)+2+j)+" "+str(nlig1*(lig1len+1)+i*(lig2len+1)+3+j)+" "+str(nlig1*(lig1len+1)+i*(lig2len+1)+4+j)+"\n")
#		anglenum = anglenum + 1
#    	
#
#
#Write LAMMPS in.file:
#g1=open("in.sphere", "w")

     
    
     
