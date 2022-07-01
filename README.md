# Shakti-Chaco-Metamaterial-PBC-Anneal

Code for my project on spin-ice inspired metamaterials in the lab of [Yair Shokef](https://shokef.tau.ac.il/index.php/research) at Tel Aviv University. 

Please see our papers on the Arxiv for background : 

[Square spin-ice metamaterial and dynamic hysteresis](https://doi.org/10.48550/arXiv.2010.05227)

[Chaco metamaterial inspired by frustrated Shakti artificial spin ice](
https://doi.org/10.48550/arXiv.2204.04000)


C Code for Shakti/Chaco Metamaterial Simulations, Gradient Descent, PBC, "Anneal " Version k2 increases slowly 

# Main files:

Shakti_PBC_Small_Anneal_k2.c  main program, implements Periodic Boundary Conditions
can start from a desired k2 value and increase k2 in increments...

ran2.h function for uniform random numbers.  
ran2.c 

# Input Files: 

RandomSpinSamples/Square_Spins_Lx_4_Ly_4_StartStateNUMBER_0.txt
specify starting values for square edges s_i(0)

use StartStateNUMBER to specify different initial configurations 

in/out as 1 or -1, rows of 4 spins for each square, comma separated
list for each of Lx*Ly squares
1,-1,1,-1
-1,1,-1,1
-1,1,-1,1
etc...

RandomSpinSamples/Diamond_Spins_Lx_4_Ly_4_StartStateNUMBER_0.txt
edge/spin s_i(0) between back to back triangles... 

up/down = 1/2 right/left = 3/4, list of values for Lx*Ly spins
1
3
2
4
2
3
etc...

# Compiling the excutable: 
I typically use gcc to compile the c files to get an executable for the program, can then run on personal PC or on TAU cluster...

gcc -o shakti_pbc_anneal.exe Shakti_PBC_Small_Aneal_k2.c ran2.h ran2.c -lm -O3

-o executable_name.exe give a name to the compiled executable...

flags: -lm is needed so some math functions like sqrt() work,  -03 optimizes the code for efficiency...

# Running the executable: 

./shakti_pbc_anneal.exe Lx Ly gamma R1 R2 R3 T_relax data_set_number start_state_number T_step total_steps

Lx, Ly number of squares for total Lx*Ly squares under PBC, there are also 2*Lx*Ly triangles between the squares

gamma I usually use 0.001, this is the time step or the learning rate for the gradient descent algorithm. 

Spring Constants

k1 --> framework springs, edges of the squares
k2 --> interaction springs inside of the squares
k3 --> framework springs between pairs of back to back triangles
k4 --> interaction springs inside of the triangles

k1 = 10.0 by default

R1 ---> k2 = R1*k1  Ratio of interaction springs to k1 springs, I generally used 0.01 < R1 < 1.0
R2 ---> k4 = R2*k2  R2 allows the ineraction springs inside the triangles to be different from those in the squares, to change energy of triangles relative to squares.  
			default is R2 = 1 so k4 = k2, simplest choice
R3 ----> k3 = R3*k1  R3 in case you want to change the stiffness of framework springs between the 

alpha --> compression factor, alpha = 0.9 is hard coded into the .c file.  Can change if you want... 

T_relax  time to allow the energy to relax after the initial condition is prepared 

data_set_number -- just a number to label different runs if you want to save them

start_state_number -- number that goes with the input files in order to specify all the spin value/edge displacements in the initial condition

T_step --- time to wait between each increment in k2.  

right now, the increments in k2 are log spaced, have to go into the code to change this, modify the step size. 

total_steps --- number of k2 increase increments to carry out. Total elasped time for the run will be T_relax + T_step*total_steps . Small systems take a few minutes. 

