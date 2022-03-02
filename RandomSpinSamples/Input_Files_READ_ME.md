# This folder has the input files for specifying the inital spin displacements as an initial configuration of the lattice


## Square spins

row of 4 spins -1/1 for each of the Lx*Ly squares in the lattices

-1,1,-1,1 \n
1,-1,1,-1, \n
...
...
...
1,-1,1,-1 \n

## Spins between back to back triangles 

1/2 = up/down   3/4 = right left   there are Lx*Ly triangle pairs 

1 \n
4 \n 
2
3
1
4
1
3
....


## I have included the starting file for the 4 x 4 PBC lattice in the paper

0 = the mechanical ground state

2,6,7 = special ordered higher energy states (not used in the paper) 

177, 5001-5005 These are unpaired defect states based on the magnetic Shakti Ground states

8001-8005  These are randomly oriented squares , central spins chosen randomly

12001-12005 THere are randomly chosen spin values 
