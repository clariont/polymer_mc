#!/bin/bash

#rm equil_all.lammpstrj
mpiexec -np 4 ~/Programs/lmp_protein/src/lmp_openmpi -in in.equil

natoms=$( head -n 4 poly.lammpstrj | tail -n 1 )
nlines=$( echo "$natoms + 9" | bc )

tail -n $nlines poly.lammpstrj > equil.lammpstrj
cat equil.lammpstrj >> equil_all.lammpstrj






