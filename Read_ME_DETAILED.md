
# DETAILED OVERVIEW OF .C PROGRAM FOR THE SHAKTI/CHACO METAMATERIAL SPRING NETWORK SIMULATIONS...

## Here is a detailed overview of the main C program for the Shakti/Chaco Mechanical Metamaterial
Simulation C FILE. 

"Shakti_PBC_Small_Anneal_k2.c"  

The Code is very long linear file, with all the steps/sections to carry the simulation in order. 
I will try to describe how each section works. 

Probably if you are a grad student or postdoc reading this, you will want to just get an idea of 
how the simulation/algorimths works and write your own code...

Sorry, I know reading someone elses code is never fun... 


if you have questions , please email me carl.merrigan@fulbrightmail.org 

### Good LUCK! 






include statements at the beginning for needed libraries...

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# function declatations

	"create_square_cell" and "create_diamdond_cell" these add the nodes/bonds to construct the lattice

	"make_data_directory" creates the folders for saving text data files  NetworkData/Lx_?_Ly_?/k2_?_k3_?_k4_?/g_?/Run_?_State_?/

	"update_k2_spring_constant"  this increments the increases in the k2 value... 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# main( ) program...


parameter defintions 

	like Lx, Ly, gamma, Nv - number of nodes, Nb number of bonds, etc... 

	parameter inputs from command line : Input: Lx Ly gamma R1 R2 R3 T_relax data_set_number start_state_number T_step total_steps\n

	if the program doesn't get a list of the right size with appropriate input arguments, it should terminate with an error.  Can modify this piece of code to add more desired input parameters.   

	A lot of other parameters like alpha,  the k2_step size, are hard coded in .c file.    

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# file names and opening files for data input and output 

	node_coords_file =  each line is a lisy x1,y1,x2,y2,x3,y3...xNv,yNv,current_step     of all node cooridates at each output time step. last entry is current time step

	energy_file = output the total energy and the average energy per square / average energy per triangle 

	bond_pairs_file = node i, node j, l_relax, k_spring    list of details for all the bonds in the lattice, 

	node_bonds_file = each row is a node i, and the entries are bond indices for bonds attached to that node (bonds are number in order of their creation)

	unit_cell_nodes_file = each row is list of nodes associated with each square in the lattices

	triangleA_nodes_file = each row is nodes associated with half the triangles upper/right triangles in each pair

	triangleB_nodes file = each row is nodes associated with half the triangles lower/left triangle in each pair 

	shakti_ground_state_file = input file to specify initial displacement spin values -1/1 for spins on edges of the squares

	shakti_ground_state_islands_file = input file to specify initial spin values for edges/spins between back to back triangles 


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# declare, allocate memory, and initialize variables for all lists and arrays in the main program

in c everything has to be declared, defined very precisely.  

	Examples . Matrix/Array containing the current (x,y) values for each node "coords" 

	double **coords           (** means pointer to a pointer to a chunk of memory that is size of double ) 

	// coordinate array holds x,y index for each NODE
   	 coords = (double **)malloc(2 * sizeof(double *));
    	if(coords != NULL){
    	    for(i=0;i<2;i++){
        	    coords[i] = (double *)malloc(Nv_max * sizeof(double));
            	if(coords[i] != NULL){
                	for(j=0;j<Nv_max;j++) coords[i][j] = -1.0;
            	}else{
                	printf("\nMalloc Error.\n"); exit(1);
            	}
        }
    	}else{
        	printf("\nMalloc Error. \n"); exit(1);
    	}

	hopeful the names of various lists/arrays are fairly understandable, should be able to figure out what they mean by where they are used

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# input files shakti_ground_state_file and shakti_ground_state_islands_file specify square/ centeral spin values for the initial state

	square_spin_directions[i][j] holds spin values for square spins i = 0,1,2,3 for the spin and j = 0,1,2... Nv-1    (note C indices always start from 0 ) 


	diamond_spin_values[i]   specified as 1/2 or 3/4 for up/down  and right/left for the direction of spins between back to back triangles ....

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# CREATING the lattice

## these loops sweep over all the squares and triangles to add all the nodes and bonds to the appropriate lists/arrays

loop over squares for x = 0...Lx for y = 0...Ly  creates squares (0,0) , (0,1) ...(0,Ly) ... (Lx-1,Ly-1). 

loop over triangles addes the k3 bonds and k4 for bonds defining the internal structure of the back to back triangle pairs 


look inside "create_square_cell" and "create_diamond_center" functions for details

the first arguments (x,y) are the center of the square, or (x,y) is the center of the pair of back to back triangles

important! the lattice constant is a = (1 + 2*cos(30))*alpha 

important! this version is for periodic Boundary Conditions, so all coordinates are limited with a box [0,Lx*y] and [0,Ly*a]

	

the lattice constant is for the prestressed compressed lattice, where all distances 1 --> alpha<1 is the compression factor for pinning the corners of the squares...

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# Prepare the initial condition 

## The initial lattice is created in the undeformed state, with everything compressed by factor alpha to implement the prestress....

spin values for squares / triangles read in from files are applied to the appropriate edges


this is the section of the code where randomization of the intial node positions could also be applied. 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

# MAIN PROGRAM LOOP         

## finally after all the preparation, here is the main simulation loop implementing the gradient descent energy minimization algorithm 

all updates take place within a large while loop 

	while(current_step < n_steps) {

		.....

	}

calculating node displacement updates

	1. loop over nodes i =  0 ... Nv-1 

	2. loop over bonds j connected to node i 

	3. compute the energy gradient (force ) with respect to x_i and y_i 

	4. displace node i += - 1  * gamma * gradient_ i     in gradient descent, each position increment is small time_step/learning rate times gradient in the energy 


NOTE: updates to the node positions are calculated for all nodes based on the current node positions x,y in "coords"   

Loop over all nodes to update the positions, after the updates have all been calculated

	1. loop over all nodes, update the value in "coords" 

	2. output the "coords" values to the coordinates value at desired output intervals  "output_step_size" 


### Energy Output step

	at desired intervals, compute and output the total energy (= sum all bonds ) and the average energy of squares/triangles

	NOTE: the loops for defining the average energy for squares and triangles count the framework bonds k1, k3 by a factor 1/2

	      the k1 bonds are all shared between two squares
	
	      the k3 bonds are all shared between back to back triangles

	       k2, k4 bonds are unique, internal to their respect squares and triangles... 

	This the convention ensures that summing the energy of all squares and all triangles give back the total energy = sum of all bond energies 


### ANNEALING PROTOCOL:  (slowly increasing k2/k1 and letting the energy relax between increments ...) 

	inside the MAIN LOOP is an if statement for incrementing the value of k2 at desired steps.....

	if(current_step>n_relax){  // allows for the initial condition to relax first

		if((current_step-n_relax)%n_step==0){  // increments k2 after waiting n_step  iterations

		     printf("\n step = %d k2_last = %f k2_next = %f",current_step,k2_last,k2_next);
	

		     k2_last = update_k2_spring_constants(Nb,bond_spring_constant,k2_last,k2_next);  all spring constants with value "k2_last" replaced by "k2_next"

		     k2_next = k2_last + k2_last*growth_factor;  steps chosen so that increments are spaced logarithmically , not linearly... 

    		     //getchar();
		}

	}


### COMMENT out this section to have a simulation with fixed k2, k4 values for all time..                        

#### CAUTION: I coded the annealing protocol assuming k2 = k4 , so things might not work if different values are used for the triangle internal bonds k4. 


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# code after the int main( ... ) program loop gives all the needed function defitions. 
