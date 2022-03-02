clearvars; close all; clc;

%state_number = 0 % figure 2(a)
state_number = 177 % figure 2(b)
%state_number = 8001 % figure 2(c)
%state_number = 12001 % figure2(d)

% uncomment the right state number to plot generate the .fig file 

arrow_color = [255 215 0]/255  % this will change arrow color. 

arrow_color_edge = arrow_color;
shape = [0.5,0.5,0.5,0.05];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx =  4; Ly = 4; 

lw = 1;

patch_alpha = 1;

data_set_number = 1;

n_steps = 	-1; step_size = 1;  %% 390 gives k2/k1 0.015 , small value where all states are stable

%% remember I used n_steps = 500 for the paper figure (a) (b) (c) (d) !!!! 

n_steps = 1500; step_size = 100;

n_steps = 600; step_size = 100; 

alpha_square = 0.9;

a = (1 + 2*cos(pi/6))*alpha_square;

gamma = [10]*10^-4; % gamma is the gradient descent step size. 

T_relax = 300;

n_relax = ceil(T_relax/gamma);

k1_square = 10.0;

R1 = 0.01

R2 = 1

R3 = 1.0;

k2_square = R1*k1_square;

k4_tri = R2*k2_square;

k3_tri = R3*k1_square;



%6    26    32    36    40    41    42    45    71    73    75    77    93

k3 = 0
k_reduced = 0

directory_name = ['NetworkData/Lx_',num2str(Lx),'_Ly_',num2str(Ly),'/k2_',num2str(k2_square,'%.6f'),'_k3_',num2str(k3_tri,'%.6f'),'_k4_',num2str(k4_tri,'%.6f'),'/g_',num2str(10000*gamma,'%.4f'),'/Run_',num2str(data_set_number),'_State_',num2str(state_number),'/']

node_filename = [directory_name,'NodeCoords_Lx_',num2str(Lx),'_Ly_',num2str(Ly),'_gamma_',num2str(10000*gamma,'%.3f'),'_k_',num2str(k3,'%.1f'),'_',num2str(data_set_number),'.txt'];
bond_filename = [directory_name,'BondPairs_Lx_',num2str(Lx),'_Ly_',num2str(Ly),'_gamma_',num2str(10000*gamma,'%.3f'),'_k_',num2str(k3,'%.1f'),'_',num2str(data_set_number),'.txt'];

nodebond_filename = [directory_name,'NodeBonds_Lx_',num2str(Lx),'_Ly_',num2str(Ly),'_gamma_',num2str(10000*gamma,'%.3f'),'_k_',num2str(k3,'%.1f'),'_',num2str(data_set_number),'.txt'];

cell_filename = [directory_name,'UnitCellNodes_Lx_',num2str(Lx),'_Ly_',num2str(Ly),'_gamma_',num2str(10000*gamma,'%.3f'),'_k_',num2str(k3,'%.1f'),'_',num2str(data_set_number),'.txt'];

triangleA_filename = [directory_name,'TriangleANodes_Lx_',num2str(Lx),'_Ly_',num2str(Ly),'_gamma_',num2str(10000*gamma,'%.3f'),'_k_',num2str(k3,'%.1f'),'_',num2str(data_set_number),'.txt'];

triangleB_filename = [directory_name,'TriangleBNodes_Lx_',num2str(Lx),'_Ly_',num2str(Ly),'_gamma_',num2str(10000*gamma,'%.3f'),'_k_',num2str(k3,'%.1f'),'_',num2str(data_set_number),'.txt'];

fileID = fopen(node_filename);

firstline = fgetl(fileID);

Nv = (sum(firstline==',')-1)/2

energy_filename = [directory_name,'Energy_Lx_',num2str(Lx),'_Ly_',num2str(Ly),'_gamma_',num2str(10000*gamma,'%.3f'),'_k_',num2str(k3,'%.1f'),'_',num2str(data_set_number),'.txt'];

energy_data = importdata(energy_filename)

time_list_energy_file = energy_data(:,1);

energy_list = energy_data(:,2);

k2_current_energy_file = energy_data(:,3);



fclose all;

bond_pairs = importdata(bond_filename);

%diagonal_list = find(bond_pairs(:,4)==10.0);
%bond_pairs(diagonal_list,:)=[];

Nb = length(bond_pairs);

unit_cell_data = importdata(cell_filename);

unit_cell_nodes_matrix = unit_cell_data(:,1:8)+1;

triangle_A_data = importdata(triangleA_filename);

triangle_A_nodes_matrix = triangle_A_data(:,1:6)+1;

triangle_B_data = importdata(triangleB_filename);

triangle_B_nodes_matrix = triangle_B_data(:,1:6)+1;

clear unit_cell_data triangle_A_data;

bond_count_list = zeros(1,Nv);

figure(1);

set(gcf,'Units','normalized');
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);

for i = 1:Nv
    
    bond_count_list(i) = sum(bond_pairs(:,1)==(i-1)) + sum(bond_pairs(:,2)==(i-1));
    
end

fileID = fopen(node_filename);

energy_filename = [directory_name,'Energy_Lx_',num2str(Lx),'_Ly_',num2str(Ly),'_gamma_',num2str(10000*gamma,'%.3f'),'_k_',num2str(k3,'%.1f'),'_',num2str(data_set_number),'.txt'];

for current_step = -1:n_steps

data = textscan(fileID,"%f",2*Nv+1,'Delimiter',',');

display(current_step);

if(mod(current_step,step_size)==0&&(current_step>-2))

node_coords = data{1,1}';

step_current_file = node_coords(end);
time_list(current_step+2)=step_current_file;

x_list = node_coords(1,1:2:2*Nv);
y_list = node_coords(1,2:2:2*Nv);

x_list = mod(x_list+1*a,Lx*a);
y_list = mod(y_list+1*a,Ly*a); 

[~,myIndex]=min(abs(time_list_energy_file-step_current_file));

k2_current_VALUE = k2_current_energy_file(myIndex);

%% figure out patches for adding background color to differentiate squares and triangles
%% pray this code works....

x_cell = x_list(unit_cell_nodes_matrix);
y_cell = y_list(unit_cell_nodes_matrix);

x_cell = x_cell(:,[1 5 2 6 3 7 4 8]);
y_cell = y_cell(:,[1 5 2 6 3 7 4 8]);

x_cell_new =x_cell;
y_cell_new =y_cell;

for cell = 1:(Lx*Ly)

   new_diff_x = zeros(1,7); 
   new_diff_y = zeros(1,7);

   test1=ismember(1,abs(diff(x_cell(cell,:)))>(Lx/2));
   test2=ismember(1,abs(diff(y_cell(cell,:)))>(Ly/2));
   
   if(test1||test2)
   
   x_cell_new(cell,1)=x_cell(cell,1);
   y_cell_new(cell,1)=y_cell(cell,1);
   
for q = 1:7
    
    delta_test = x_cell(cell,q+1)-x_cell(cell,q);
    
    if((Lx*a-abs(delta_test))<abs(delta_test))
        delta_test = (Lx*a - abs(delta_test))*(-1)*sign(delta_test);
    end
   
    new_diff_x(q) = delta_test;
    
    x_cell_new(cell,1+q) = x_cell_new(cell,1)+sum(new_diff_x(1:q));   
    
    delta_test = y_cell(cell,q+1)-y_cell(cell,q);
    
    if((Ly*a-abs(delta_test))<abs(delta_test))
        delta_test = (Lx*a - abs(delta_test))*(-1)*sign(delta_test);
    end
   
    new_diff_y(q) = delta_test;
    
    y_cell_new(cell,1+q) = y_cell_new(cell,1)+sum(new_diff_y(1:q)); 
    
end

   end

end

x_A = x_list(triangle_A_nodes_matrix);
y_A = y_list(triangle_A_nodes_matrix);

x_A = x_A(:,[1 5 2 6 3 4]);
y_A = y_A(:,[1 5 2 6 3 4]);


x_A_new =x_A;
y_A_new =y_A;

for cell = 1:(Lx*Ly)

   new_diff_x = zeros(1,7); 
   new_diff_y = zeros(1,7);

   test1=ismember(1,abs(diff(x_A(cell,:)))>(Lx/2));
   test2=ismember(1,abs(diff(y_A(cell,:)))>(Ly/2));
   
   if(test1||test2)
   
   x_A_new(cell,1)=x_A(cell,1);
   y_A_new(cell,1)=y_A(cell,1);
   
for q = 1:5
    
    delta_test = x_A(cell,q+1)-x_A(cell,q);
    
    if((Lx*a-abs(delta_test))<abs(delta_test))
        delta_test = (Lx*a - abs(delta_test))*(-1)*sign(delta_test);
    end
   
    new_diff_x(q) = delta_test;
    
    x_A_new(cell,1+q) = x_A_new(cell,1)+sum(new_diff_x(1:q));   
    
    delta_test = y_A(cell,q+1)-y_A(cell,q);
    
    if((Ly*a-abs(delta_test))<abs(delta_test))
        delta_test = (Lx*a - abs(delta_test))*(-1)*sign(delta_test);
    end
   
    new_diff_y(q) = delta_test;
    
    y_A_new(cell,1+q) = y_A_new(cell,1)+sum(new_diff_y(1:q)); 
    
end

   end

end

x_B = x_list(triangle_B_nodes_matrix);
y_B = y_list(triangle_B_nodes_matrix);

x_B = x_B(:,[1 2 3 6 4 5]);
y_B = y_B(:,[1 2 3 6 4 5]);


x_B_new =x_B;
y_B_new =y_B;

for cell = 1:(Lx*Ly)

   new_diff_x = zeros(1,7); 
   new_diff_y = zeros(1,7);

   test1=ismember(1,abs(diff(x_B(cell,:)))>(Lx/2));
   test2=ismember(1,abs(diff(y_B(cell,:)))>(Ly/2));
   
   if(test1||test2)
   
   x_B_new(cell,1)=x_B(cell,1);
   y_B_new(cell,1)=y_B(cell,1);
   
for q = 1:5
    
    delta_test = x_B(cell,q+1)-x_B(cell,q);
    
    if((Lx*a-abs(delta_test))<abs(delta_test))
        delta_test = (Lx*a - abs(delta_test))*(-1)*sign(delta_test);
    end
   
    new_diff_x(q) = delta_test;
    
    x_B_new(cell,1+q) = x_B_new(cell,1)+sum(new_diff_x(1:q));   
    
    delta_test = y_B(cell,q+1)-y_B(cell,q);
    
    if((Ly*a-abs(delta_test))<abs(delta_test))
        delta_test = (Lx*a - abs(delta_test))*(-1)*sign(delta_test);
    end
   
    new_diff_y(q) = delta_test;
    
    y_B_new(cell,1+q) = y_B_new(cell,1)+sum(new_diff_y(1:q)); 
    
end

   end

end

SIGMA_X = 1; SIGMA_Y = 1; k2 = k2_square; k1=k1_square;


delta_x_new = SIGMA_X.*SIGMA_X.*(circshift(x_cell_new,-1,2)-x_cell_new).^2;

delta_y_new = SIGMA_Y.*SIGMA_Y.*(circshift(y_cell_new,-1,2)-y_cell_new).^2;

l_new = (delta_x_new+delta_y_new).^(1/2);

%energy_squares_k1 = 0.5*k1*sum((l_new-ones(size(l_new))).^2,2);

energy_squares_k1 = 0.5*0.5*k1*sum((l_new-ones(size(l_new))).^2,2);   %% EXTRA FACTOR HALF SO k1 SPRINGS NOT DOUBLE COUNTED. 


delta_x_new = SIGMA_X.*SIGMA_X.*(circshift(x_cell_new(:,[1:2:end]),-1,2)-x_cell_new(:,[1:2:end])).^2;

delta_y_new = SIGMA_Y.*SIGMA_Y.*(circshift(y_cell_new(:,[1:2:end]),-1,2)-y_cell_new(:,[1:2:end])).^2;

l_new = (delta_x_new+delta_y_new).^(1/2);

energy_squares_k2 = 0.5*k2_current_VALUE*sum((l_new-(2^(1/2)).*ones(size(l_new))).^2,2);

PX_squares_k2 = k1*sum((l_new-ones(size(l_new))).*(delta_x_new./SIGMA_X)./l_new,2);
PY_squares_k2 = k1*sum((l_new-ones(size(l_new))).*(delta_y_new./SIGMA_Y)./l_new,2);


%% triangleA energy 

delta_x_new = SIGMA_X.*SIGMA_X.*(circshift(x_A_new,-1,2)-x_A_new).^2;

delta_y_new = SIGMA_Y.*SIGMA_Y.*(circshift(y_A_new,-1,2)-y_A_new).^2;

l_new = (delta_x_new+delta_y_new).^(1/2);

%%energy_trianglesA_k1 = 0.5*k1*sum((l_new-ones(size(l_new))).^2,2);
energy_trianglesA_k1 = 0.5*0.5*k1*sum((l_new-ones(size(l_new))).^2,2);   %% EXTRA FACTOR HALF SO UNIT ENERGIES ADD TO TOTAL ENERGY ... 

PX_A_k1 = k1*sum((l_new-ones(size(l_new))).*(delta_x_new./SIGMA_X)./l_new,2);
PY_A_k1 = k1*sum((l_new-ones(size(l_new))).*(delta_y_new./SIGMA_Y)./l_new,2);


delta_x_new = SIGMA_X.*SIGMA_X.*(x_A_new(:,[3,5])-repmat(x_A_new(:,1),1,2)).^2;

delta_y_new = SIGMA_Y.*SIGMA_Y.*(y_A_new(:,[3,5])-repmat(y_A_new(:,1),1,2)).^2;

l_new = (delta_x_new+delta_y_new).^(1/2);

energy_trianglesA_k2 = 0.5*k2_current_VALUE*sum((l_new-ones(size(l_new))).^2,2);

PX_A_k2 = k1*sum((l_new-ones(size(l_new))).*(delta_x_new./SIGMA_X)./l_new,2);
PY_A_k2 = k1*sum((l_new-ones(size(l_new))).*(delta_y_new./SIGMA_Y)./l_new,2);


%% triangleB energy 

delta_x_new = SIGMA_X.*SIGMA_X.*(circshift(x_B_new,-1,2)-x_B_new).^2;

delta_y_new = SIGMA_Y.*SIGMA_Y.*(circshift(y_B_new,-1,2)-y_B_new).^2;

l_new = (delta_x_new+delta_y_new).^(1/2);

%%energy_trianglesB_k1 = 0.5*k1*sum((l_new-ones(size(l_new))).^2,2);
energy_trianglesB_k1 = 0.5*0.5*k1*sum((l_new-ones(size(l_new))).^2,2);   %% EXTRA FACTOR HALF SO UNIT ENERGIES ADD TO TOTAL ENERGY.... 


PX_B_k1 = k1*sum((l_new-ones(size(l_new))).*(delta_x_new./SIGMA_X)./l_new,2);
PY_B_k1 = k1*sum((l_new-ones(size(l_new))).*(delta_y_new./SIGMA_Y)./l_new,2);

delta_x_new = SIGMA_X.*SIGMA_X.*(x_B_new(:,[3,5])-repmat(x_B_new(:,1),1,2)).^2;

delta_y_new = SIGMA_Y.*SIGMA_Y.*(y_B_new(:,[3,5])-repmat(y_B_new(:,1),1,2)).^2;

l_new = (delta_x_new+delta_y_new).^(1/2);

energy_trianglesB_k2 = 0.5*k2_current_VALUE*sum((l_new-ones(size(l_new))).^2,2);

PX_B_k2 = k1*sum((l_new-ones(size(l_new))).*(delta_x_new./SIGMA_X)./l_new,2);
PY_B_k2 = k1*sum((l_new-ones(size(l_new))).*(delta_y_new./SIGMA_Y)./l_new,2);


%%


figure(1); cla;

caxis([log10(0.001) log10(1)])
colorbar

patch(x_cell_new',y_cell_new',log10(energy_squares_k1'+energy_squares_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;

patch(x_cell_new'+Lx*a,y_cell_new',log10(energy_squares_k1'+energy_squares_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;
patch(x_cell_new'-Lx*a,y_cell_new',log10(energy_squares_k1'+energy_squares_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;
patch(x_cell_new',y_cell_new'+Ly*a,log10(energy_squares_k1'+energy_squares_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;
patch(x_cell_new',y_cell_new'-Ly*a,log10(energy_squares_k1'+energy_squares_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;

patch(x_cell_new'+Lx*a,y_cell_new'+Ly*a,log10(energy_squares_k1'+energy_squares_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;


patch(x_A_new',y_A_new',log10(energy_trianglesA_k1'+energy_trianglesA_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;

patch(x_A_new'+Lx*a,y_A_new',log10(energy_trianglesA_k1'+energy_trianglesA_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;
patch(x_A_new'-Lx*a,y_A_new',log10(energy_trianglesA_k1'+energy_trianglesA_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;
patch(x_A_new',y_A_new'+Ly*a,log10(energy_trianglesA_k1'+energy_trianglesA_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;
patch(x_A_new',y_A_new'-Ly*a,log10(energy_trianglesA_k1'+energy_trianglesA_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;

patch(x_B_new',y_B_new',log10(energy_trianglesB_k1'+energy_trianglesB_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;

patch(x_B_new'+Lx*a,y_B_new',log10(energy_trianglesB_k1'+energy_trianglesB_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;
patch(x_B_new'-Lx*a,y_B_new',log10(energy_trianglesB_k1'+energy_trianglesB_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;
patch(x_B_new',y_B_new'+Ly*a,log10(energy_trianglesB_k1'+energy_trianglesB_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;
patch(x_B_new',y_B_new'-Ly*a,log10(energy_trianglesB_k1'+energy_trianglesB_k2'),'FaceAlpha',patch_alpha,'EdgeColor',[0 0 0]); hold on;


%% and now put the bond strains on top... 

%spin 1 
node1=5;node2=8; spin_index = 1;

x_diff = x_list(unit_cell_nodes_matrix(:,node1)) - x_list(unit_cell_nodes_matrix(:,node2));
y_diff = y_list(unit_cell_nodes_matrix(:,node1)) - y_list(unit_cell_nodes_matrix(:,node2));

x_diff(x_diff<-2)=x_diff(x_diff<-2)+Lx*a;
y_diff(y_diff<-2)=y_diff(y_diff<-2)+Ly*a;

x_diff(x_diff>2)=x_diff(x_diff>2)-Lx*a;
y_diff(y_diff>2)=y_diff(y_diff>2)-Ly*a;

x1_list = x_list(unit_cell_nodes_matrix(:,node2)) + x_diff*0.5;
y1_list = y_list(unit_cell_nodes_matrix(:,node2)) + y_diff*0.5;

x1_list = x1_list - floor(x1_list./(Lx*a))*(Lx*a);
y1_list = y1_list - floor(y1_list./(Ly*a))*(Ly*a);

u1_list = x_list(unit_cell_nodes_matrix(:,spin_index))-x1_list;
v1_list = y_list(unit_cell_nodes_matrix(:,spin_index))-y1_list;

%spin 2 
node1=5;node2=6; spin_index=2;

x_diff = x_list(unit_cell_nodes_matrix(:,node1)) - x_list(unit_cell_nodes_matrix(:,node2));
y_diff = y_list(unit_cell_nodes_matrix(:,node1)) - y_list(unit_cell_nodes_matrix(:,node2));

x_diff(x_diff<-2)=x_diff(x_diff<-2)+Lx*a;
y_diff(y_diff<-2)=y_diff(y_diff<-2)+Ly*a;

x_diff(x_diff>2)=x_diff(x_diff>2)-Lx*a;
y_diff(y_diff>2)=y_diff(y_diff>2)-Ly*a;

x2_list = x_list(unit_cell_nodes_matrix(:,node2)) + x_diff*0.5;
y2_list = y_list(unit_cell_nodes_matrix(:,node2)) + y_diff*0.5;

x2_list = x2_list - floor(x2_list./(Lx*a))*(Lx*a);
y2_list = y2_list - floor(y2_list./(Ly*a))*(Ly*a);

u2_list = x_list(unit_cell_nodes_matrix(:,spin_index))-x2_list;
v2_list = y_list(unit_cell_nodes_matrix(:,spin_index))-y2_list;

%spin 3
node1=6;node2=7; spin_index=3;

x_diff = x_list(unit_cell_nodes_matrix(:,node1)) - x_list(unit_cell_nodes_matrix(:,node2));
y_diff = y_list(unit_cell_nodes_matrix(:,node1)) - y_list(unit_cell_nodes_matrix(:,node2));

x_diff(x_diff<-2)=x_diff(x_diff<-2)+Lx*a;
y_diff(y_diff<-2)=y_diff(y_diff<-2)+Ly*a;

x_diff(x_diff>2)=x_diff(x_diff>2)-Lx*a;
y_diff(y_diff>2)=y_diff(y_diff>2)-Ly*a;

x3_list = x_list(unit_cell_nodes_matrix(:,node2)) + x_diff*0.5;
y3_list = y_list(unit_cell_nodes_matrix(:,node2)) + y_diff*0.5;

x3_list = x3_list - floor(x3_list./(Lx*a))*(Lx*a);
y3_list = y3_list - floor(y3_list./(Ly*a))*(Ly*a);

u3_list = x_list(unit_cell_nodes_matrix(:,spin_index))-x3_list;
v3_list = y_list(unit_cell_nodes_matrix(:,spin_index))-y3_list;

%spin 4
node1=7;node2=8; spin_index=4;

x_diff = x_list(unit_cell_nodes_matrix(:,node1)) - x_list(unit_cell_nodes_matrix(:,node2));
y_diff = y_list(unit_cell_nodes_matrix(:,node1)) - y_list(unit_cell_nodes_matrix(:,node2));

x_diff(x_diff<-2)=x_diff(x_diff<-2)+Lx*a;
y_diff(y_diff<-2)=y_diff(y_diff<-2)+Ly*a;

x_diff(x_diff>2)=x_diff(x_diff>2)-Lx*a;
y_diff(y_diff>2)=y_diff(y_diff>2)-Ly*a;

x4_list = x_list(unit_cell_nodes_matrix(:,node2)) + x_diff*0.5;
y4_list = y_list(unit_cell_nodes_matrix(:,node2)) + y_diff*0.5;

x4_list = x4_list - floor(x4_list./(Lx*a))*(Lx*a);
y4_list = y4_list - floor(y4_list./(Ly*a))*(Ly*a);

u4_list = x_list(unit_cell_nodes_matrix(:,spin_index))-x4_list;
v4_list = y_list(unit_cell_nodes_matrix(:,spin_index))-y4_list;



x_islands = mean(x_list(triangle_A_nodes_matrix(:,[4,5])),2)';
y_islands = mean(y_list(triangle_A_nodes_matrix(:,[4,5])),2)';

u_islands = x_list(triangle_A_nodes_matrix(:,1))-x_islands;
v_islands = y_list(triangle_A_nodes_matrix(:,1))-y_islands;






 

x_pairs = x_list((bond_pairs(:,1:2)+1)');
y_pairs = y_list((bond_pairs(:,1:2)+1)');

delta_x = (x_pairs(1,:)-x_pairs(2,:));
delta_y = (y_pairs(1,:)-y_pairs(2,:));

delta_x = min(abs(delta_x),Lx*a-abs(delta_x));
delta_y = min(abs(delta_y),Ly*a-abs(delta_y));

delta_l = (1*(delta_x.^2 + delta_y.^2).^(1/2)) - (bond_pairs(:,3)');




strain_list = delta_l./(bond_pairs(:,3)');

delta_l = strain_list;

%jetmap = colormap('jet');

delta_l = (delta_l + 0.3)/(2*0.3);

color_index = floor(delta_l*256);

color_index(color_index<1)=1;
color_index(color_index>256)=256;

%springcolors = jetmap(color_index,:);

n = 100; % the number of levels in said colormap
t = linspace(0,1,n)';
w = [1 1 1];
p = [75 0 130]/255; % rescale 8 bit colors into [0,1].
CM = (1-t)*w + t*p;
colormap(CM)

figure(1); 

axis([-0.05 (2*a+0.05) -0.05 (2*a+0.05)])

list1 = bond_count_list==5;
list2 = bond_count_list==7;
list3 = bond_count_list==6;




x_box = linspace(0,(Lx)*a,100);
y_box = zeros(1,100);

plot(x_box,y_box,'k-','Linewidth',1);
plot(y_box,x_box,'k-','Linewidth',1);
plot(x_box,Ly*a*ones(1,100),'k-','Linewidth',1);
plot(Ly*a*ones(1,100),x_box,'k-','Linewidth',1);

x_box = linspace(4,(Lx)*a-4,100);

scatter(x_list(bond_count_list==8),y_list(bond_count_list==8),20,[0 0 0],'filled');

set(gca,'Xtick',[]);
set(gca,'Ytick',[]);

set(gca,'Visible','off');

axis square;

x1 = a - 1.4*a; 
x2 = Lx*a - x1; 

y1 = x1; 
y2 = x2; 

axis([x1 x2 y1 y2])

%% HERE IS WHERE THE ARROWS ARE DRAWN... 

if(current_step==n_steps)

    %arrow_color = [0.4940 0.1840 0.5560

    %arrow_color = [255 203 37]./256;  %% see power point 

    %arrow_color = [1 1 1];

    %shape = [0.5,0.5,0.5,0.05];

arrows(x1_list,y1_list,u1_list,v1_list,'Cartesian',shape,'FaceColor',arrow_color,'EdgeColor',arrow_color_edge,'LineWidth',1);
arrows(x2_list,y2_list,u2_list,v2_list,'Cartesian',shape,'FaceColor',arrow_color,'EdgeColor',arrow_color_edge,'LineWidth',1);

arrows(x3_list,y3_list,u3_list,v3_list,'Cartesian',shape,'FaceColor',arrow_color,'EdgeColor',arrow_color_edge,'LineWidth',1);

arrows(x4_list,y4_list,u4_list,v4_list,'Cartesian',shape,'FaceColor',arrow_color,'EdgeColor',arrow_color_edge,'LineWidth',1);
arrows(x_islands,y_islands,u_islands,v_islands,'Cartesian',shape,'FaceColor',arrow_color,'EdgeColor',arrow_color_edge,'LineWidth',1);

%quiver(x1_list,y1_list,u1_list,v1_list,'AutoScale','off','Color',[1 0 0],'MaxHeadSize',1); 
%quiver(x2_list,y2_list,u2_list,v2_list,'AutoScale','off','Color',[0 0 0]);
%quiver(x3_list,y3_list,u3_list,v3_list,'AutoScale','off','Color',[0 0 0]);
%quiver(x4_list,y4_list,u4_list,v4_list,'AutoScale','off','Color',[0 0 0]);

%quiver(x_islands,y_islands,u_islands,v_islands,'Autoscale','off','Color',[0 0 0]);

end














set(gca,'Xtick',[]);
set(gca,'Ytick',[]);

set(gca,'Visible','on');

axis square;

x1 = a - 1.4*a; 
x2 = Lx*a - x1; 

y1 = x1; 
y2 = x2; 

axis([x1 x2 y1 y2])



v = [x1 y1; x2 y1;x2 0;x1 0];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',0.5,'EdgeColor','none')

v = [x1 y2; x2 y2;x2 Ly*a;x1 Ly*a];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',0.5,'EdgeColor','none')

v = [x1 0;0 0;0 Ly*a;x1 Ly*a];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',0.5,'EdgeColor','none')

v = [Ly*a 0;x2 0;x2 Ly*a;Ly*a Ly*a];
f = [1 2 3 4];
patch('Faces',f,'Vertices',v,'FaceColor','black','FaceAlpha',0.5,'EdgeColor','none')

test = 0;

x1 = a - 1.35*a; 
x2 = Lx*a - x1

y1 = x1; 
y2 = x2; 

axis([x1 x2 y1 y2])

title(['state number = ',num2str(state_number),' k_2/k_1 = ',num2str(k2_current_VALUE/10),' '])

end


end

saveas(gcf,['FIGURE2_FINAL_FIGURES/ENERGYARROWS_Sample_time_',num2str(step_current_file*gamma),'_k2_current_',num2str(k2_current_VALUE),'_STATE_',num2str(state_number),'_RATIO_1000_Zoom.fig']);
