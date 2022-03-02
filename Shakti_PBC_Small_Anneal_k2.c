
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "ran2.h";

int max ( int a, int b );
int max ( int a, int b );

void print_coordinates(double **coords,int **square_vertex_index,int Nv,int step);
void print_bonds(int **bond_pairs,int Nb,int step);
void print_node_bonds(int **node_bonds,int Nv,int step);

int create_square_cell(double x0,double y0,double theta,double **coords,int **square_vertex_index,int **node_bonds,int **bond_pairs,int **squares_bonds,double *bond_equilibrium_length,double *bond_spring_constant,int *corner_node_check,int *Nv,int *Nb,int *N_corner,int square_index,double alpha,double k_internal,double k_small_diagonal,FILE * unit_cell_file,int Lx,int Ly,double a);

int create_diamond_center(double x0,double y0,double theta,double **coords,int **square_vertex_index,int **node_bonds,int **bond_pairs,int **diamonds_bonds,double *bond_equilibrium_length,double *bond_spring_constant,int *Nv,int *Nb,int diamond_center_index,double k3_tri,double k4_tri,double alpha,double k_prestress,FILE *triangleA_file,FILE *triangleB_file,int Lx,int Ly,double a);

double bond_gradient_x(double x1,double y1,double x2,double y2,double l_eq,double k,int Lx,int Ly,double a,double SIGMA_X,double SIGMA_Y);
double bond_gradient_y(double x1,double y1,double x2,double y2,double l_eq,double k,int Lx,int Ly,double a,double SIGMA_X,double SIGMA_Y);
double bond_energy(double x1,double y1,double x2,double y2,double l_eq,double k,int Lx,int Ly,double a,double SIGMA_X,double SIGMA_Y);

double bond_pressure_X(double x1,double y1,double x2,double y2,double l_eq,double k,int Lx,int Ly,double a,double SIGMA_X,double SIGMA_Y);
double bond_pressure_Y(double x1,double y1,double x2,double y2,double l_eq,double k,int Lx,int Ly,double a,double SIGMA_X,double SIGMA_Y);

int make_data_directory(int Lx,int Ly,double gamma,double k2_square,double k3_tri,double k4_tri,int data_set_number,int start_state_number);

double update_k2_spring_constants(int Nb,double *bond_spring_constant,double k2_old,double k2_new);

int main(int argc,char **argv)
{
    long seed = 41482342;

    int Lx,Ly,Nv_max,Nb_max,N_corner_total; // Lx by Ly square cells

    int max_bonds = 10;  // list of up to max_bonds=10 for each node 

    int Nv=0,Nb=0,N_corner=0; // number vertices, number bonds, number of pinned corner nodes of the squares

    double T_relax;
    int n_steps;
    int n_relax;
    double gamma;
    double epsilon = 0; // scale for uniform random noise to initial node coordinates 

    double alpha_square=0.900; // compression factor, side length of a square is reduced to 2*alpha to give spontaneous displacments 
    double k_internal_square=0.0; // "prestress" springs from old version, instead of pinning the corners 

    double delta = sqrt((1-(alpha_square*alpha_square))); // amplitude for spin displacements that leaves the k1 framework bonds unstrained 
r
    double k2_square = 10.0; // k1 square the framework are all set to 10.0 by default, k2_square are interaction bonds inside the square. 
    double k3_tri = 10.0; // default for k3 = k1
    double k4_tri = 10.0; //k4 are interactions inside each triangle 

    double R1=1.0,R2=1.0,R3=1.0;

    double k2_last,k2_next;

    double growth_factor = pow(10.0,0.01)-1.0;

    //double k2_step = 0.05;

    double T_step;
    int n_step,total_steps; // careful, n_step is number steps in each k2 internal of time, n_steps is total number of steps entire run. 

    int data_set_number;

    int output_step_size = 1;

    int energy_output_stepsize=100;

    int start_state_number = 0;  // for reading from a file with initial spin configurations specified 

    if((argc-1)== 11){
        Lx = 2*atoi(argv[1]);
        Ly = 2*atoi(argv[2]);
        gamma = atof(argv[3]);

	energy_output_stepsize = floor(0.1/gamma);

	R1 = atof(argv[4]);
	R2 = atof(argv[5]);
	R3 = atof(argv[6]);

	k2_square = R1*10.0;  // fix squares themselves to have k1 = 10

	k4_tri = R2*k2_square;  // R2 = 1 gives triangles and squares with same k2 internal bonds

	k3_tri = R3*10.0;  // R3 if you want bonds between back to back triangle bonds different from k1 bonds for square framework...

        T_relax = atof(argv[7]);

        n_relax = ceil(T_relax/gamma);

        data_set_number = atoi(argv[8]);

	start_state_number = atoi(argv[9]);   // fix the folder names, the input eventually.... 

	T_step = atof(argv[10]);

	total_steps = atoi(argv[11]);

	n_step = ceil(T_step/gamma);

	n_steps = n_relax + n_step*total_steps;

	output_step_size = (int) floor(n_steps/2000);

	if(output_step_size<1){
		output_step_size = 1;
	}

        seed += data_set_number*8118;

        Nv_max = 2*8*Lx*Ly;
        Nb_max = 2*14*Lx*Ly;

	N_corner_total = 2*Lx*Ly; 

	k2_last = k2_square;

	k2_next = k2_last + k2_last*growth_factor;

    }else{
        printf("Input: Lx Ly gamma R1 R2 R3 T_relax data_set_number start_state_number T_step total_steps\n");
        exit(1);
    }

/////       CREATE DATA DIRECTORY FOR OUTPUT TEXT FILES   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Create a Directory for the Files from this Run
    char directoryname[200];

    sprintf(directoryname,"NetworkData/Lx_%d_Ly_%d/k2_%.6f_k3_%.6f_k4_%.6f/g_%.4f/Run_%d_State_%d",Lx,Ly,k2_square,k3_tri,k4_tri,10000*gamma,data_set_number,start_state_number);

    int status;
    status = make_data_directory(Lx,Ly,gamma,k2_square,k3_tri,k4_tri,data_set_number,start_state_number);

    if (!status)
        printf("Directory created %d  \n",status);
    else {
        printf("Unable to create directory %d  \n",status);
        exit(1);
    }

    FILE * node_coords_file;
    FILE * bond_pairs_file;
    FILE * node_bonds_file;
    FILE * energy_file;
    FILE * unit_cell_nodes_file;
    FILE * triangleA_file;
    FILE * triangleB_file;

    FILE * shakti_ground_state_file; // files to read in GS spin configuration of a Shakti lattice
    FILE * shakti_ground_state_islands_file;

    //Names of all the files for data output 
    char node_coords_filename[200];
    sprintf(node_coords_filename,"%s/NodeCoords_Lx_%d_Ly_%d_gamma_%.3f_k_%.1f_%d.txt",directoryname,Lx,Ly,10000*gamma,k_internal_square,data_set_number);

    char bond_pairs_filename[200];
    sprintf(bond_pairs_filename,"%s/BondPairs_Lx_%d_Ly_%d_gamma_%.3f_k_%.1f_%d.txt",directoryname,Lx,Ly,10000*gamma,k_internal_square,data_set_number);

    char node_bonds_filename[200];
    sprintf(node_bonds_filename,"%s/NodeBonds_Lx_%d_Ly_%d_gamma_%.3f_k_%.1f_%d.txt",directoryname,Lx,Ly,10000*gamma,k_internal_square,data_set_number);
    char energy_filename[200];
    sprintf(energy_filename,"%s/Energy_Lx_%d_Ly_%d_gamma_%.3f_k_%.1f_%d.txt",directoryname,Lx,Ly,10000*gamma,k_internal_square,data_set_number);

    char unit_cell_nodes_filename[200];
    sprintf(unit_cell_nodes_filename,"%s/UnitCellNodes_Lx_%d_Ly_%d_gamma_%.3f_k_%.1f_%d.txt",directoryname,Lx,Ly,10000*gamma,k_internal_square,data_set_number);

    char triangleA_filename[200];
    sprintf(triangleA_filename,"%s/TriangleANodes_Lx_%d_Ly_%d_gamma_%.3f_k_%.1f_%d.txt",directoryname,Lx,Ly,10000*gamma,k_internal_square,data_set_number);

    char triangleB_filename[200];
    sprintf(triangleB_filename,"%s/TriangleBNodes_Lx_%d_Ly_%d_gamma_%.3f_k_%.1f_%d.txt",directoryname,Lx,Ly,10000*gamma,k_internal_square,data_set_number);

    char shakti_ground_state_filename[200];
    sprintf(shakti_ground_state_filename,"RandomSpinSamples/Square_Spins_Lx_%d_Ly_%d_StartStateNUMBER_%d.txt",Lx,Ly,start_state_number);

    char shakti_ground_state_islands_filename[200];
    sprintf(shakti_ground_state_islands_filename,"RandomSpinSamples/Diamond_Spins_Lx_%d_Ly_%d_StartStateNUMBER_%d.txt",Lx,Ly,start_state_number);
	
    //open the files for writing data 
    if(1){
        unit_cell_nodes_file = fopen(unit_cell_nodes_filename,"w");
        if(!unit_cell_nodes_file){
        printf("Can't open file <%s> for writing.\n",unit_cell_nodes_filename);
        exit(1);
        }
    }

    if(1){
        triangleA_file = fopen(triangleA_filename,"w");
        if(!triangleA_file){
        printf("Can't open file <%s> for writing.\n",triangleA_filename);
        exit(1);
        }
    }

    if(1){
        triangleB_file = fopen(triangleB_filename,"w");
        if(!triangleB_file){
        printf("Can't open file <%s> for writing.\n",triangleB_filename);
        exit(1);
        }
    }

    if(1){
        node_coords_file = fopen(node_coords_filename,"w");
        if(!node_coords_file){
        printf("Can't open file <%s> for writing.\n",node_coords_filename);
        exit(1);
        }
    }

    if(1){
        bond_pairs_file = fopen(bond_pairs_filename,"w");
        if(!bond_pairs_file){
        printf("Can't open file <%s> for writing.\n",bond_pairs_filename);
        exit(1);
        }
    }

    if(1){
        node_bonds_file = fopen(node_bonds_filename,"w");
        if(!node_bonds_file){
        printf("Can't open file <%s> for writing.\n",node_bonds_filename);
        exit(1);
        }
    }

    if(1){
        energy_file = fopen(energy_filename,"w");
        if(!energy_file){
        printf("Can't open file <%s> for writing.\n",energy_filename);
        exit(1);
        }
    }

    // open the files to read in the initial spin configuration
    if(1){
        shakti_ground_state_file = fopen(shakti_ground_state_filename,"r");
        if(!shakti_ground_state_file){
        printf("Can't open file <%s> for reading.\n",shakti_ground_state_filename);
        exit(1);
        }
    }

    if(1){
        shakti_ground_state_islands_file = fopen(shakti_ground_state_islands_filename,"r");
        if(!shakti_ground_state_islands_file){
        printf("Can't open file <%s> for reading.\n",shakti_ground_state_islands_filename);
        exit(1);
        }
    }


///////////              ALLOCATE MEMORY / INITIAL ALL MATRICSE AND LISTS NEEDED FOR THE MAIN SIMULATION    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    int i,j;

    double a = (1.0 + 2.0*cos(M_PI/6))*alpha_square;  // a is the lattice constant , size of PBC system is Lx*a by Ly*a

    double b = 1.0/cos(M_PI/6)*alpha_square;

    double c = sin(M_PI/6)/cos(M_PI/6)*alpha_square;

    double rand_uniform;

    int node1,node2;
    double l_eq,k_eq;

    double x,y,theta;
    int square_index_current,spin_index_current,spin_inout_current;
    int bond_count;

    int s1,s2,s3,s4;

    double x_square_center,y_square_center;

    double **coords;
    int **node_bonds;
    int **bond_pairs;
    int **squares_bonds;
    int **diamonds_bonds;

    //double SIGMA = 1; // sigma is variable parameter that tunes the units of length, based on overall pressure inside the periodic system.  
    double SIGMA_X = 1; // idea is sigma is a "metric" that changes to allow expansion and contraction of the periodic system, like surface of a torus that inflates/contracts...  
    double SIGMA_Y = 1; //fixing the corner nodes, so not using the sigma dynamics in this version, but might want again another time...

    double *bond_equilibrium_length;

    double *bond_spring_constant;

    double *node_displacement_x;
    double *node_displacement_y;

    int *corner_node_check;

    int **square_vertex_index;
    int square_cell_count = 0,diamond_center_count=0;

    int **square_spin_directions;
    int *diamond_spin_values; 

    //allocating memory for all the various arrays needed....

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

    // list of relaxed lengths for each bond
    bond_equilibrium_length = (double *)malloc(Nb_max*sizeof(double));
    if(bond_equilibrium_length!=NULL){
        for(i=0;i<Nb_max;i++){
            bond_equilibrium_length[i]=-1.0;
        }
    }else{
        printf("\nMalloc Error. \n"); exit(1);
    }

    // list of spring constants for each bond 
    bond_spring_constant = (double *)malloc(Nb_max*sizeof(double));
    if(bond_spring_constant!=NULL){
        for(i=0;i<Nb_max;i++){
            bond_spring_constant[i]=-1.0;
        }
    }else{
        printf("\nMalloc Error. \n"); exit(1);
    }

    // holds the next increment in x for each node 
    node_displacement_x = (double *)malloc(Nv_max*sizeof(double));
    if(node_displacement_x!=NULL){
        for(i=0;i<Nv_max;i++){
            node_displacement_x[i]=0.0;
        }
    }else{
        printf("\nMalloc Error. \n"); exit(1);
    }

    //holds the next increment in y for each node
    node_displacement_y = (double *)malloc(Nv_max*sizeof(double));
    if(node_displacement_y!=NULL){
        for(i=0;i<Nv_max;i++){
            node_displacement_y[i]=0.0;
        }
    }else{
        printf("\nMalloc Error. \n"); exit(1);
    }

    // 1 to label all the nodes which are corners of the squares. These are pinned in place and skipped over so they don't move during the simulation...
    corner_node_check = (int *)malloc(Nv_max*sizeof(int));
    if(corner_node_check!=NULL){
        for(i=0;i<Nv_max;i++){
            corner_node_check[i]=0;
        }
    }else{
        printf("\nMalloc Error. \n"); exit(1);
    }

    square_vertex_index = (int **)malloc(2 * sizeof(int *));
    if(square_vertex_index != NULL){
        for(i=0;i<2;i++){
            square_vertex_index[i] = (int *)malloc(Nv_max * sizeof(int));
            if(square_vertex_index[i] != NULL){
                for(j=0;j<Nv_max;j++) square_vertex_index[i][j] = -1;
            }else{
                printf("\nMalloc Error.\n"); exit(1);
            }
        }
    }else{
        printf("\nMalloc Error. \n"); exit(1);
    }

    node_bonds = (int **)malloc(max_bonds * sizeof(int *));
    if(node_bonds != NULL){
        for(i=0;i<max_bonds;i++){
            node_bonds[i] = (int *)malloc(Nv_max * sizeof(int));
            if(node_bonds[i] != NULL){
                for(j=0;j<Nv_max;j++) node_bonds[i][j] = -1;
            }else{
                printf("\nMalloc Error.\n"); exit(1);
            }
        }
    }else{
        printf("\nMalloc Error.\n"); exit(1);
    }

    squares_bonds = (int **)malloc(Lx*Ly * sizeof(int *));
    if(squares_bonds != NULL){
        for(i=0;i<Lx*Ly;i++){
            squares_bonds[i] = (int *)malloc(14 * sizeof(int));
            if(squares_bonds[i] != NULL){
                for(j=0;j<14;j++) squares_bonds[i][j] = -1;
            }else{
                printf("\nMalloc Error.\n"); exit(1);
            }
        }
    }else{
        printf("\nMalloc Error.\n"); exit(1);
    }

    diamonds_bonds = (int **)malloc(Lx*Ly * sizeof(int *));
    if(diamonds_bonds != NULL){
        for(i=0;i<Lx*Ly;i++){
            diamonds_bonds[i] = (int *)malloc(15 * sizeof(int));
            if(diamonds_bonds[i] != NULL){
                for(j=0;j<7;j++) diamonds_bonds[i][j] = -1;
            }else{
                printf("\nMalloc Error.\n"); exit(1);
            }
        }
    }else{
        printf("\nMalloc Error.\n"); exit(1);
    }

    bond_pairs = (int **)malloc(2 * sizeof(int *));
    if(bond_pairs != NULL){
        for(i=0;i<2;i++){
            bond_pairs[i] = (int *)malloc(Nb_max * sizeof(int));
            if(bond_pairs[i] != NULL){
                for(j=0;j<Nb_max;j++) bond_pairs[i][j] = -1;
            }else{
                printf("\nMalloc Error.\n"); exit(1);
            }
        }
    }else{
        printf("\nMalloc Error. \n"); exit(1);
    }

    square_spin_directions = (int **)malloc(4 * sizeof(int *));
    if(square_spin_directions != NULL){
        for(i=0;i<4;i++){
            square_spin_directions[i] = (int *)malloc(Lx*Ly * sizeof(int));
            if(square_spin_directions[i] != NULL){
                for(j=0;j<Lx*Ly;j++) square_spin_directions[i][j] = 0;
            }else{
                printf("\nMalloc Error.\n"); exit(1);
            }
        }
    }else{
        printf("\nMalloc Error. \n"); exit(1);
    }

    diamond_spin_values = (int *)malloc(Lx*Ly*sizeof(int));
    if(diamond_spin_values!=NULL){
        for(i=0;i<Lx*Ly;i++){
            diamond_spin_values[i]=-1.0;
        }
    }else{
        printf("\nMalloc Error. \n"); exit(1);
    }

//////////////////////////////////////////////  READ IN SPIN VALUES / INITIAL EDGE DISPLACEMENTS FOR PREPARTING THE DESIRED INITIAL STATE /////////////////////////////////////////////////////////////////////////////////////////////    


    //reading in the values of the square cell spins. 
    for(i=0;i<Lx*Ly;i++){

    fscanf(shakti_ground_state_file,"%d,%d,%d,%d",&s1,&s2,&s3,&s4);

    printf("\t%d\t%d\t%d\t%d\n",s1,s2,s3,s4);

    square_spin_directions[0][i]=s1;
    square_spin_directions[1][i]=s2;
    square_spin_directions[2][i]=s3;
    square_spin_directions[3][i]=s4;

    }

    fclose(shakti_ground_state_file);

    for(i=0;i<(Lx)*(Ly);i++){
    fscanf(shakti_ground_state_islands_file,"%d",&s1); 
    
    printf("\t%d\n",s1);

    //getchar();
    
    diamond_spin_values[i]=s1;
    }

    fclose(shakti_ground_state_islands_file);


 /////////////////////////////////////  CREATE ALL THE NODES AND BONDS IN THE NETWORK /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 

    // important a is the lattice constant . first all squares are created, and then next bonds defining internal structure of the back to back triangles and added on top ... 

    printf("\n Creating Network.... Lx = %d, Ly = %d",Lx,Ly);


    //loop to create all squares....
    for(i=0;i<Lx;i++){
        x = i*(a);
        for(j=0;j<Ly;j++){
            y = j*(a);

            if(((i+j)%2)==0){
                theta = M_PI/6;
            }else if(((i+j)%2)==1){
                theta = -M_PI/6;
            }

            //printf("\n x = %f y = %f theta = %f",x,y,theta);

	    create_square_cell(x,y,theta,coords,square_vertex_index,node_bonds,bond_pairs,squares_bonds,bond_equilibrium_length,bond_spring_constant,corner_node_check,&Nv,&Nb,&N_corner,square_cell_count,alpha_square,k_internal_square,k2_square,unit_cell_nodes_file,Lx,Ly,a);

	    square_cell_count++;

            //print_coordinates(coords,square_vertex_index,Nv,0);
	    //getchar();
	    //print_node_bonds(node_bonds,Nv,0); 
            //getchar();
        }
    }

    if(N_corner!=N_corner_total){
	printf("\nWRONG NUMBER OF CORNER NODES !!! \n\n"); exit(0);
    }

    // loop to add bonds for the triangles   (by DIAMOND I Mean pairs of back to back triangles ). 
    for(i=0;i<(Lx);i++){

            x = 0.5*a + i*a;

        for(j=0;j<(Ly);j++){

            y  = 0.5*a + j*a;

            if(((i+j)%2)==0){
                theta = 0.0;
            }else if(((i+j)%2)==1){
                theta = M_PI/2;
            }

            create_diamond_center(x,y,theta,coords,square_vertex_index,node_bonds,bond_pairs,diamonds_bonds,bond_equilibrium_length,bond_spring_constant,&Nv,&Nb,diamond_center_count,k3_tri,k4_tri,alpha_square,k_internal_square,triangleA_file,triangleB_file,Lx,Ly,a);

            diamond_center_count++;

        }

    }

    print_coordinates(coords,square_vertex_index,Nv,0);

    fclose(triangleA_file);
    fclose(triangleB_file);
    fclose(unit_cell_nodes_file);

    //getchar();

    //print_bonds(bond_pairs,Nb,0);

    //print_node_bonds(node_bonds,Nv,0);
    


    //print the initial undeformed lattice coordinates as first line in the coordinates file 
    for(i=0;i<Nv;i++){

        x = coords[0][i]; y = coords[1][i];

        fprintf(node_coords_file,"%.9f,%.9f,",x,y);
        if(i==Nv-1){
            //fprintf(node_coords_file,"%f,%f,%d,\n",SIGMA_X,SIGMA_Y,0);  // last entries are sigma_x sigma_y and the current step.
	    fprintf(node_coords_file,"%d,\n",0);
        }

    }


    //print matrix row for each node, and list of bond indices associated with that node ... 
    for(i=0;i<Nv;i++){

	for(j=0;j<max_bonds;j++){
        
        fprintf(node_bonds_file," %d,",node_bonds[j][i]);
        if(j==max_bonds-1){
            fprintf(node_bonds_file,"\n");
        }
        }
    }

    fclose(node_bonds_file);

    // print matrix with all bond parameters to a file 
    for(i=0;i<Nb;i++){
        node1 = bond_pairs[0][i]; node2 = bond_pairs[1][i];

        l_eq = bond_equilibrium_length[i];

	k_eq = bond_spring_constant[i];

        fprintf(bond_pairs_file,"%d,%d,%f,%f\n",node1,node2,l_eq,k_eq);
    }

    fclose(bond_pairs_file);


//////////////////////////////////////// IMPOSE THE INITIAL CONDITION / INITIAL SPIN VALUES ON THE SQUARES AND TRIANGLES ///////////////////////////////////////////////////////////////////////////////////////    

    for(i=0;i<Nv;i++){

        x = coords[0][i]; y = coords[1][i];

        square_index_current = square_vertex_index[0][i]; // 0 to Lx*Ly-1, total number of square unit cells, z = 4 vertices for SHAKTI
        spin_index_current = square_vertex_index[1][i];

        bond_count = 0;

        for(j=0;j<10;j++){
            if(node_bonds[j][i]!=-1){
                bond_count++;
            }
        }

        if((bond_count==5)){

        x_square_center = floor(square_index_current/Ly)*a;
        y_square_center = (square_index_current - Ly*(x_square_center/a))*a;

	if((y-y_square_center)>1){
		y_square_center = y_square_center + Ly*a;
	}else if((y-y_square_center)<-1){
		y_square_center = y_square_center - Ly*a;
	}

	if((x-x_square_center)>1){
                x_square_center = x_square_center + Lx*a;
        }else if((x-x_square_center)<-1){
                x_square_center = x_square_center - Lx*a;
        }

        if(spin_index_current==0){
                spin_inout_current = square_spin_directions[1][square_index_current];
        }else if(spin_index_current==1){
                spin_inout_current = square_spin_directions[0][square_index_current];
        }else if(spin_index_current==2){
                spin_inout_current = square_spin_directions[3][square_index_current];
        }else if(spin_index_current==3){
                spin_inout_current = square_spin_directions[2][square_index_current];
        }else{
                printf("\n\nError. This is not a face node spin!!!\n\n");exit(1);
        }

        coords[0][i]=x+spin_inout_current*delta*(x - x_square_center)/alpha_square;
        coords[1][i]=y+spin_inout_current*delta*(y- y_square_center)/alpha_square;

        }else if(bond_count==6){

		printf("%d\n",square_index_current); 

        spin_inout_current = diamond_spin_values[square_index_current];  // being lazy, this is a value now, 1 up , 2 down 3 right 4 left

        if(spin_inout_current==1){
                coords[0][i]=x;
                coords[1][i]=y+delta;
        }else if(spin_inout_current==2){
                coords[0][i]=x;
                coords[1][i]=y-delta;
        }else if(spin_inout_current==3){
                coords[0][i]=x+delta;
                coords[1][i]=y;
        }else if(spin_inout_current==4){
                coords[0][i]=x-delta;
                coords[1][i]=y;
        }else{
                printf("\n\nError. This is not an island spin value!!!\n\n"); exit(1);
        }



	}

	coords[0][i]=coords[0][i] - floor( (double) coords[0][i]/(Lx*a))*Lx*a; // keep coordinates within periodic box (0,Lx*a) (0,Ly*a)
        coords[1][i]=coords[1][i] - floor( (double) coords[1][i]/(Ly*a))*Ly*a;

	fprintf(node_coords_file,"%.9f,%.9f,",coords[0][i],coords[1][i]);
        if(i==Nv-1){
		fprintf(node_coords_file,"%d,\n",0);
        }

    }

    double x1,y1,x2,y2;

    double E_current,E_square_avg_current,E_diamond_avg_current;

    double PX_current;
    double PY_current;

    E_current = 0;


    // compute the total lattice energy = sum over all bonds in the lattice 
    for(j=0;j<Nb;j++){

        node1=bond_pairs[0][j]; node2 = bond_pairs[1][j];

        x1 = coords[0][node1]; y1 = coords[1][node1];
        x2 = coords[0][node2]; y2 = coords[1][node2];

        E_current += bond_energy(x1,y1,x2,y2,bond_equilibrium_length[j],bond_spring_constant[j],Lx,Ly,a,SIGMA_X,SIGMA_Y);

    }

    //fprintf(energy_file,"%f,\n",(double) E_current);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////// MAIN GRADIENT DESCENT LOOP FOR MINIMIZING THE ENERGY ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

    int current_step=0;

    int current_bond;

    double grad_xj,grad_yj;

    while(current_step<n_steps){

	if(current_step%5000==0){    
         printf("\n step = %d SX = %f SY = %f ",current_step,SIGMA_X,SIGMA_Y);
	}

        if(current_step>n_relax){

		if((current_step-n_relax)%n_step==0){

		     printf("\n step = %d k2_last = %f k2_next = %f",current_step,k2_last,k2_next);
	

		     k2_last = update_k2_spring_constants(Nb,bond_spring_constant,k2_last,k2_next);

		     k2_next = k2_last + k2_last*growth_factor;

    		     //getchar();
		}

	}

        //getchar();

        for(i=0;i<Nv;i++){

            x1 = coords[0][i];
            y1 = coords[1][i];

            node_displacement_x[i]=0.0;
            node_displacement_y[i]=0.0;

            for(j=0;j<max_bonds;j++){

                current_bond = node_bonds[j][i];

                if(current_bond!=-1){

                    node1 = bond_pairs[0][current_bond];
                    node2 = bond_pairs[1][current_bond];

                    l_eq = bond_equilibrium_length[current_bond];

                    k_eq = bond_spring_constant[current_bond];

                    if(node2==i){
                        node2 = node1;
                    }

                    x2 = coords[0][node2]; y2 = coords[1][node2];

		    grad_xj = bond_gradient_x(x1,y1,x2,y2,l_eq,k_eq,Lx,Ly,a,SIGMA_X,SIGMA_Y);
	            grad_yj = bond_gradient_y(x1,y1,x2,y2,l_eq,k_eq,Lx,Ly,a,SIGMA_X,SIGMA_Y);

	            node_displacement_x[i] += -1*gamma*grad_xj;
                    node_displacement_y[i] += -1*gamma*grad_yj;

                }

            }

        }
	

        for(i=0;i<Nv;i++){

            x = coords[0][i]; y = coords[1][i];

	    if(corner_node_check[i]==0){

            coords[0][i]=x+node_displacement_x[i];
            coords[1][i]=y+node_displacement_y[i];

            coords[0][i]=coords[0][i] - floor( (double) coords[0][i]/(Lx*a))*Lx*a; // keep coordinates within periodic box (0,Lx*a) (0,Ly*a)
            coords[1][i]=coords[1][i] - floor( (double) coords[1][i]/(Ly*a))*Ly*a;
	    
	    }
	    /////////
            if((current_step)%output_step_size==0){
            fprintf(node_coords_file,"%.9f,%.9f,",coords[0][i],coords[1][i]);
            if(i==Nv-1){
                //fprintf(node_coords_file,"%.9f,%.9f,%d,\n",SIGMA_X,SIGMA_Y,current_step);
		fprintf(node_coords_file,"%d,\n",current_step);
            }
            }

        }


	if(current_step%energy_output_stepsize==0){
        E_current = 0;
	E_square_avg_current = 0;
	E_diamond_avg_current = 0;

        for(j=0;j<Nb;j++){

            node1=bond_pairs[0][j]; node2 = bond_pairs[1][j];

            x1 = coords[0][node1]; y1 = coords[1][node1];
            x2 = coords[0][node2]; y2 = coords[1][node2];

            E_current += bond_energy(x1,y1,x2,y2,bond_equilibrium_length[j],bond_spring_constant[j],Lx,Ly,a,SIGMA_X,SIGMA_Y);
        }

	for(i=0;i<Lx*Ly;i++){
		for(j=0;j<12;j++){
			
			current_bond = squares_bonds[i][j];

			node1=bond_pairs[0][current_bond]; node2 = bond_pairs[1][current_bond];

		        x1 = coords[0][node1]; y1 = coords[1][node1];
            		x2 = coords[0][node2]; y2 = coords[1][node2];
			
			if(j<8){
				// add 0.5*k1 to energy of each unit, then total energy = sum all bonds is recovered when all units summed . 
				E_square_avg_current += 0.5*bond_energy(x1,y1,x2,y2,bond_equilibrium_length[current_bond],bond_spring_constant[current_bond],Lx,Ly,a,SIGMA_X,SIGMA_Y);
			}else{
				E_square_avg_current += bond_energy(x1,y1,x2,y2,bond_equilibrium_length[current_bond],bond_spring_constant[current_bond],Lx,Ly,a,SIGMA_X,SIGMA_Y);
			}

		}
	}

	for(i=0;i<Lx*Ly;i++){
		for(j=0;j<15;j++){
			current_bond = diamonds_bonds[i][j];

			node1=bond_pairs[0][current_bond]; node2 = bond_pairs[1][current_bond];

		       	x1 = coords[0][node1]; y1 = coords[1][node1];
            		x2 = coords[0][node2]; y2 = coords[1][node2];

			if(j==2||j==5){
				// these are shared between two triangles so counted twice = 0.5*2 
				E_diamond_avg_current += bond_energy(x1,y1,x2,y2,bond_equilibrium_length[current_bond],bond_spring_constant[current_bond],Lx,Ly,a,SIGMA_X,SIGMA_Y);
			}else if(j>6){
				// these are the edges of the triangles so weighted by 0.5*k1 
				E_diamond_avg_current += 0.5*bond_energy(x1,y1,x2,y2,bond_equilibrium_length[current_bond],bond_spring_constant[current_bond],Lx,Ly,a,SIGMA_X,SIGMA_Y);
			}else{
				//internal bonds to each triangle wheighted 1*k2 
				E_diamond_avg_current += bond_energy(x1,y1,x2,y2,bond_equilibrium_length[current_bond],bond_spring_constant[current_bond],Lx,Ly,a,SIGMA_X,SIGMA_Y);
			}
		}
	}

	E_square_avg_current = E_square_avg_current/(Lx*Ly);
	E_diamond_avg_current = E_diamond_avg_current/(2*Lx*Ly);

        //fprintf(energy_file,"%d,%f,%f,%f,%f,%f,\n",current_step,E_current,(1000000)*SIGMA_X,(1000000)*SIGMA_Y,PX_current,PY_current);
	fprintf(energy_file,"%d,%f,%f,%f,%f,\n",current_step,E_current,k2_last,E_square_avg_current,E_diamond_avg_current);
	}

        current_step++;
    }

    fclose(energy_file);
    fclose(node_coords_file);

    exit(0);
}

int min ( int a, int b ) { return a < b ? a : b; }
int max ( int a, int b ) { return a > b ? a : b; }

void print_coordinates(double **coords,int **square_vertex_index,int Nv,int step){

    int i;
    double x,y;

    int square_index,spin_index;

    printf("\n Step = %d NODE LIST Nv = %d\n",step,Nv);

    for(i=0;i<Nv;i++){
        x = coords[0][i]; y = coords[1][i];
	square_index = square_vertex_index[0][i];
        spin_index = square_vertex_index[1][i];
        printf("\t(%f,%f)\t [%d,%d] \t %d\n",x,y,square_index,spin_index,i);
	//getchar();
    }
}

void print_bonds(int **bond_pairs,int Nb,int step){

    int i;
    int node1,node2;

    printf("\n Step = %d BOND LIST -> Nb = %d\n",step,Nb);

    for(i=0;i<Nb;i++){
        node1 = bond_pairs[0][i]; node2 = bond_pairs[1][i];
        printf("\t(%d,%d) %d\n",node1,node2,i);
    }
}

void print_node_bonds(int **node_bonds,int Nv,int step){

    int i,j,bond;

    printf("\n Step = %d NODE BOND LIST -> Nv = %d\n",step,Nv);

    for(i=0;i<Nv;i++){

        for(j=0;j<8;j++){
        bond = node_bonds[j][i];
        printf(" %d,",bond);
        if(j==7){
            printf("\n");
        }
        }
    }
}

int create_diamond_center(double x0,double y0,double theta,double **coords,int **square_vertex_index,int **node_bonds,int **bond_pairs,int **diamonds_bonds,double *bond_equilibrium_length,double *bond_spring_constant,int *Nv,int *Nb,int diamond_center_index,double k3_tri,double k4_tri,double alpha,double k_prestress,FILE *triangleA_file,FILE *triangleB_file,int Lx,int Ly,double a){

    int i,j;
    double x,y;
    double delta = 0.2;

    double k_reduced = 2*sqrt(2)*k_prestress*cos((5*M_PI/12));

    double x_node[9];
    double y_node[9];

    int node_list_number[9];

    x_node[0]=0.0; y_node[0]=0.0;
    x_node[1]=cos(1*M_PI/3); y_node[1]=sin(1*M_PI/3);
    x_node[2]=cos(2*M_PI/3); y_node[2]=sin(2*M_PI/3);
    x_node[3]=cos(3*M_PI/3); y_node[3]=sin(3*M_PI/3);
    x_node[4]=cos(4*M_PI/3); y_node[4]=sin(4*M_PI/3);
    x_node[5]=cos(5*M_PI/3); y_node[5]=sin(5*M_PI/3);
    x_node[6]=cos(6*M_PI/3); y_node[6]=sin(6*M_PI/3);
    x_node[7]=0.0; y_node[7]=2*sin(M_PI/3);
    x_node[8]=0.0; y_node[8]=-2*sin(M_PI/3);

    for(i=0;i<9;i++){
        x = x_node[i]; y = y_node[i];

        x_node[i]= floor(10000000*(cos(theta)*x - sin(theta)*y))/10000000;
        y_node[i]= floor(10000000*(sin(theta)*x + cos(theta)*y))/10000000;

        //printf("\n x = %f y = %f",x_node[i],y_node[i]);
    }

    int node1,node2,node_min,node_max;
    int test_min,test_max;
    int bond_found;
    int bond_place,place_test;

    int first_bond_node[15];
    int second_bond_node[15];

    first_bond_node[0]=0;second_bond_node[0]=1;
    first_bond_node[1]=0;second_bond_node[1]=2;
    first_bond_node[2]=0;second_bond_node[2]=3;
    first_bond_node[3]=0;second_bond_node[3]=4;
    first_bond_node[4]=0;second_bond_node[4]=5;
    first_bond_node[5]=0;second_bond_node[5]=6;
    first_bond_node[6]=3;second_bond_node[6]=6;

    //to find bonds on the edges of the triangles... 
    first_bond_node[7]=1;second_bond_node[7]=6;
    first_bond_node[8]=1;second_bond_node[8]=7;
    first_bond_node[9]=2;second_bond_node[9]=7;
    first_bond_node[10]=2;second_bond_node[10]=3;
    first_bond_node[11]=3;second_bond_node[11]=4;
    first_bond_node[12]=4;second_bond_node[12]=8;
    first_bond_node[13]=5;second_bond_node[13]=8;
    first_bond_node[14]=5;second_bond_node[14]=6;

    //create needed nodes, and add to global node coordinate list
    for(i=0;i<9;i++){
        x = x0+x_node[i];
        y = y0+y_node[i];

        x = x - floor((double) x/(Lx*a))*Lx*a; //coords variables always in box (0,Lx*a)  (0,Ly*a)
        y = y - floor((double) y/(Ly*a))*Ly*a;

        node_list_number[i]=-1;

        for(j=0;j<*Nv;j++){
            if(((x-coords[0][j])*(x-coords[0][j])<delta)&&((y-coords[1][j])*(y-coords[1][j])<delta)){
                node_list_number[i]=j;

                //printf("\n %i %j %f %f ",i,j,(x-coords[0][j]),(y-coords[1][j]));
            }
        }

        if(node_list_number[i]==-1){
            *Nv=*Nv+1;
            coords[0][*Nv-1]=x;coords[1][*Nv-1]=y;
            node_list_number[i]=*Nv-1;
	    if(i==0){
            square_vertex_index[0][*Nv-1]=diamond_center_index;
	    }
            //printf("\n CREATING CENTER NODE \n");

        }

        if(i==0||i==1||i==2||i==3||i==6){
            fprintf(triangleA_file,"%d,",node_list_number[i]);
        }else if(i==7){
	    fprintf(triangleA_file,"%d,\n",node_list_number[i]);
	}

        if(i==0||i==3||i==4||i==5||i==6){
            fprintf(triangleB_file,"%d,",node_list_number[i]);
        }else if(i==8){
	    fprintf(triangleB_file,"%d,\n",node_list_number[i]);
	}

    }

    //write down bond pairs between nodes of this square
    for(i=0;i<15;i++){

        node1 = node_list_number[first_bond_node[i]];
        node2 = node_list_number[second_bond_node[i]];


	if(i==6){
        //printf("\nbond = %i n1 = %d n2 = %d",i,node1,node2);
	}

        node_min=min(node1,node2);
        node_max=max(node1,node2);
   
	if(i==6){
        //printf("\nbond = %i min = %d max = %d",i,node_min,node_max);
        }

        bond_found = 0;
        j = 0;

        while(bond_found==0&&j<*Nb){
            test_min = bond_pairs[0][j];
            test_max = bond_pairs[1][j];

            if((node_min==test_min)&&(node_max==test_max)){
                bond_found=1;
            }else{
                j++;
            }
        }

	diamonds_bonds[diamond_center_index][i]=j;
	
	if(i==6){
	//	printf(" bound found? =  %d  j = %d  ",bond_found,j);
	}

        if(bond_found==0){
            *Nb=*Nb+1;
            bond_pairs[0][*Nb-1]=node_min; bond_pairs[1][*Nb-1]=node_max;

            if(i==2||i==5){
                bond_equilibrium_length[*Nb-1]=1.0;
                bond_spring_constant[*Nb-1]=k3_tri;   // k3, the back to back triangle bonds.... 
            }else if(i==6){
	        bond_equilibrium_length[*Nb-1]=2*alpha;	    
		bond_spring_constant[*Nb-1]=k_reduced;  // prestress bond, set to zero in this version   
	    }else{
		bond_equilibrium_length[*Nb-1]=1.0;
		bond_spring_constant[*Nb-1]=k4_tri;   //k4, or interaction bonds inside each triangle...
	    }

            bond_place = 0;
            place_test = node_bonds[0][node_min];
            while(place_test!=-1){
                bond_place++;
                place_test = node_bonds[bond_place][node_min];
            }

            node_bonds[bond_place][node_min]=*Nb-1;

            bond_place = 0;
            place_test = node_bonds[0][node_max];
            while(place_test!=-1){
                bond_place++;
                place_test = node_bonds[bond_place][node_max];
            }

            node_bonds[bond_place][node_max]=*Nb-1;

        }

    }

    return 0;

}

double bond_gradient_x(double x1,double y1,double x2,double y2,double l_eq,double k,int Lx,int Ly,double a,double SIGMA_X,double SIGMA_Y){

    double l_current;

    double x_diff = (x1-x2);
    double y_diff = (y1-y2);

    if(fabs(x_diff)>fabs(Lx*a-fabs(x_diff))){

	if(x_diff>0){    
	    x_diff = -1*fabs(Lx*a-fabs(x_diff));
	}else{
	    x_diff = fabs(Lx*a-fabs(x_diff));
	}

    }

    if(fabs(y_diff)>fabs(Ly*a-fabs(y_diff))){

        if(y_diff>0){
            y_diff = -1*fabs(Ly*a-fabs(y_diff));
        }else{
            y_diff = fabs(Ly*a-fabs(y_diff));
        }

    }

    l_current = x_diff*x_diff*SIGMA_X*SIGMA_X + y_diff*y_diff*SIGMA_Y*SIGMA_Y;

    l_current = sqrt(l_current);

    double grad_x1 = k*(l_current - l_eq)*(x_diff*SIGMA_X*SIGMA_X)/l_current;

    return grad_x1;
}

double bond_gradient_y(double x1,double y1,double x2,double y2,double l_eq,double k,int Lx,int Ly,double a,double SIGMA_X,double SIGMA_Y){

    double l_current;

    double x_diff = (x1-x2);
    double y_diff = (y1-y2);

    if(fabs(x_diff)>fabs(Lx*a-fabs(x_diff))){

        if(x_diff>0){
            x_diff = -1*fabs(Lx*a-fabs(x_diff));
        }else{
            x_diff = fabs(Lx*a-fabs(x_diff));
        }

    }
    
    if(fabs(y_diff)>fabs(Ly*a-fabs(y_diff))){

        if(y_diff>0){
            y_diff = -1*fabs(Ly*a-fabs(y_diff));
        }else{
            y_diff = fabs(Ly*a-fabs(y_diff));
        }

    }
    
    l_current = x_diff*x_diff*SIGMA_X*SIGMA_X + y_diff*y_diff*SIGMA_Y*SIGMA_Y;

    l_current = sqrt(l_current);

    double grad_y1 = k*(l_current - l_eq)*(y_diff*SIGMA_Y*SIGMA_Y)/l_current;

    return grad_y1;
}

double bond_energy(double x1,double y1,double x2,double y2,double l_eq,double k,int Lx,int Ly,double a,double SIGMA_X,double SIGMA_Y){

    double E;

    double l_current;

    double x_diff = (x1-x2);
    double y_diff = (y1-y2);

    if(fabs(x_diff)>fabs(Lx*a-fabs(x_diff))){
            x_diff = fabs(Lx*a-fabs(x_diff));
    }

    if(fabs(y_diff)>fabs(Ly*a-fabs(y_diff))){
            y_diff = fabs(Ly*a-fabs(y_diff));
    }

    l_current = x_diff*x_diff*SIGMA_X*SIGMA_X + y_diff*y_diff*SIGMA_Y*SIGMA_Y;
    l_current = sqrt(l_current);

    E = 0.5*k*(l_current-l_eq)*(l_current-l_eq);

    return E;
}

double bond_pressure_X(double x1,double y1,double x2,double y2,double l_eq,double k,int Lx,int Ly,double a,double SIGMA_X,double SIGMA_Y){

    double P;

    double l_current;

    double x_diff = (x1-x2);
    double y_diff = (y1-y2);

    if(fabs(x_diff)>fabs(Lx*a-fabs(x_diff))){
            x_diff = fabs(Lx*a-fabs(x_diff));
    }

    if(fabs(y_diff)>fabs(Ly*a-fabs(y_diff))){
            y_diff = fabs(Ly*a-fabs(y_diff));
    }

    l_current = x_diff*x_diff*SIGMA_X*SIGMA_X + y_diff*y_diff*SIGMA_Y*SIGMA_Y;
    l_current = sqrt(l_current);

    P = k*(l_current-l_eq)*(x_diff*x_diff*SIGMA_X)/l_current;

    return P;
}

double bond_pressure_Y(double x1,double y1,double x2,double y2,double l_eq,double k,int Lx,int Ly,double a,double SIGMA_X,double SIGMA_Y){

    double P;

    double l_current;

    double x_diff = (x1-x2);
    double y_diff = (y1-y2);

    if(fabs(x_diff)>fabs(Lx*a-fabs(x_diff))){
            x_diff = fabs(Lx*a-fabs(x_diff));
    }

    if(fabs(y_diff)>fabs(Ly*a-fabs(y_diff))){
            y_diff = fabs(Ly*a-fabs(y_diff));
    }

    l_current = x_diff*x_diff*SIGMA_X*SIGMA_X + y_diff*y_diff*SIGMA_Y*SIGMA_Y;
    l_current = sqrt(l_current);

    P = k*(l_current-l_eq)*(y_diff*y_diff*SIGMA_Y)/l_current;

    return P;
}

int create_square_cell(double x0,double y0,double theta,double **coords,int **square_vertex_index,int **node_bonds,int **bond_pairs,int **squares_bonds,double *bond_equilibrium_length,double *bond_spring_constant,int *corner_node_check,int *Nv,int *Nb,int *N_corner,int square_index,double alpha,double k_internal,double k_small_diagonal,FILE * unit_cell_file,int Lx,int Ly,double a){

    int i,j; 
    double x,y; 
    double delta = 0.2; 

    double x_node[8];
    double y_node[8];

    int node_list_number[8];

    x_node[0]=1.0; y_node[0]=0.0;
    x_node[1]=0.0; y_node[1]=1.0;
    x_node[2]=-1.0; y_node[2]=0.0;
    x_node[3]=0.0; y_node[3]=-1.0;
    x_node[4]=1.0; y_node[4]=1.0;
    x_node[5]=-1.0; y_node[5]=1.0;
    x_node[6]=-1.0; y_node[6]=-1.0;
    x_node[7]=1.0; y_node[7]=-1.0;

    for(i=0;i<8;i++){

        x = alpha*x_node[i]; y = alpha*y_node[i];

        //x_node[i] = cos(theta)*x - sin(theta*y);
	//y_node[i] = sin(theta)*x + cos(theta*y);

        x_node[i]= floor(1000000*(cos(theta)*x - sin(theta)*y))/1000000;
        y_node[i]= floor(1000000*(sin(theta)*x + cos(theta)*y))/1000000;

        //printf("\n x = %f y = %f",x_node[i],y_node[i]);

    }    

    int node1,node2,node_min,node_max;
    int test_min,test_max;
    int bond_found;
    int bond_place,place_test;

    int first_bond_node[14];
    int second_bond_node[14];

    first_bond_node[0]=0;second_bond_node[0]=4;
    first_bond_node[1]=1;second_bond_node[1]=4;
    first_bond_node[2]=1;second_bond_node[2]=5;
    first_bond_node[3]=2;second_bond_node[3]=5;
    first_bond_node[4]=2;second_bond_node[4]=6;
    first_bond_node[5]=3;second_bond_node[5]=6;
    first_bond_node[6]=3;second_bond_node[6]=7;
    first_bond_node[7]=0;second_bond_node[7]=7;
    first_bond_node[8]=0;second_bond_node[8]=3;
    first_bond_node[9]=0;second_bond_node[9]=1;
    first_bond_node[10]=1;second_bond_node[10]=2;
    first_bond_node[11]=2;second_bond_node[11]=3;
    first_bond_node[12]=4;second_bond_node[12]=6;
    first_bond_node[13]=5;second_bond_node[13]=7;

    //create needed nodes, and add to global node coordinate list
    for(i=0;i<8;i++){
        x = x0+x_node[i];
        y = y0+y_node[i];

	x = x - floor((double) x/(Lx*a))*Lx*a; //coords variables always in box (0,Lx*a)  (0,Ly*a)
	y = y - floor((double) y/(Ly*a))*Ly*a;

        node_list_number[i]=-1;

        for(j=0;j<*Nv;j++){
            if(((x-coords[0][j])*(x-coords[0][j])<delta)&&((y-coords[1][j])*(y-coords[1][j])<delta)){
                node_list_number[i]=j;
            }
        }

        if(node_list_number[i]==-1){
            *Nv=*Nv+1;
            coords[0][*Nv-1]=x;coords[1][*Nv-1]=y;
            node_list_number[i]=*Nv-1;
	    square_vertex_index[0][*Nv-1]=square_index;
	    square_vertex_index[1][*Nv-1]=i;

	    if(i>3){
		*N_corner = *N_corner + 1; 
		corner_node_check[*Nv-1]=1;
	    }
        }

        fprintf(unit_cell_file,"%d,",node_list_number[i]);

        if(i==7){
            fprintf(unit_cell_file,"0,\n");
        }

    }    

    //write down bond pairs between nodes of this square
    for(i=0;i<14;i++){

        node1 = node_list_number[first_bond_node[i]];
        node2 = node_list_number[second_bond_node[i]];

        //printf("\nbond = %i n1 = %d n2 = %d",i,node1,node2);

        node_min=min(node1,node2);
        node_max=max(node1,node2);

        //printf("\nbond = %i min = %d max = %d",i,node_min,node_max);

        bond_found = 0; 
        j = 0; 

        while(bond_found==0&&j<*Nb){
            test_min = bond_pairs[0][j];
            test_max = bond_pairs[1][j];

            if((node_min==test_min)&&(node_max==test_max)){
                bond_found=1;
            }else{
                j++;
            }
        }

	//j should have the abs index for the bond associated with this square
	squares_bonds[square_index][i]=j;

        if(bond_found==0){
            *Nb=*Nb+1;
            bond_pairs[0][*Nb-1]=node_min; bond_pairs[1][*Nb-1]=node_max;

            if(i<8){
                bond_equilibrium_length[*Nb-1]=1.0;
                bond_spring_constant[*Nb-1]=10.0;
            }else if(i<12){
                bond_equilibrium_length[*Nb-1]=sqrt(2);
                bond_spring_constant[*Nb-1]=k_small_diagonal;   // k2 bonds, interactions inside each square... 
            }else{
                bond_equilibrium_length[*Nb-1]=alpha*2*sqrt(2);
                bond_spring_constant[*Nb-1]=k_internal;  // old versions used long diagonal springs to generate prestress, now corners are pinned 
            }

            bond_place = 0; 
            place_test = node_bonds[0][node_min];
            while(place_test!=-1){
                bond_place++;
                place_test = node_bonds[bond_place][node_min];
            }

            node_bonds[bond_place][node_min]=*Nb-1;

            bond_place = 0; 
            place_test = node_bonds[0][node_max];
            while(place_test!=-1){
                bond_place++;
                place_test = node_bonds[bond_place][node_max];
            }

            node_bonds[bond_place][node_max]=*Nb-1;

        }

    }    

    return 0;

}

int make_data_directory(int Lx,int Ly,double gamma,double k2_square,double k3_tri,double k4_tri,int data_set_number,int start_state_number){

    int status;

    char directoryname1[200];
    char directoryname2[200];
    char directoryname3[200];
    char directoryname4[200];
    char directoryname5[200];


    sprintf(directoryname1,"NetworkData");
    sprintf(directoryname2,"NetworkData/Lx_%d_Ly_%d",Lx,Ly);
    sprintf(directoryname3,"NetworkData/Lx_%d_Ly_%d/k2_%.6f_k3_%.6f_k4_%.6f",Lx,Ly,k2_square,k3_tri,k4_tri);
    sprintf(directoryname4,"NetworkData/Lx_%d_Ly_%d/k2_%.6f_k3_%.6f_k4_%.6f/g_%.4f",Lx,Ly,k2_square,k3_tri,k4_tri,10000*gamma);
    sprintf(directoryname5,"NetworkData/Lx_%d_Ly_%d/k2_%.6f_k3_%.6f_k4_%.6f/g_%.4f/Run_%d_State_%d",Lx,Ly,k2_square,k3_tri,k4_tri,10000*gamma,data_set_number,start_state_number);

    status = mkdir(directoryname1,0700);
    status = mkdir(directoryname2,0700);
    status = mkdir(directoryname3,0700);
    status = mkdir(directoryname4,0700);
    status = mkdir(directoryname5,0700);


    return status;
}

double update_k2_spring_constants(int Nb,double *bond_spring_constant,double k2_old,double k2_new){

     double k2_current;

     int count = 0;

     for(int i=0;i<Nb;i++){

	k2_current = bond_spring_constant[i];

	if(fabs(k2_current-k2_old)<0.0001){

		bond_spring_constant[i]=k2_new;
		count++;
		//printf("\n i=%d count = %d k2_old = %f k2_new = %f",i,count,k2_current,k2_new);

	}


     }

     return k2_new;
}

