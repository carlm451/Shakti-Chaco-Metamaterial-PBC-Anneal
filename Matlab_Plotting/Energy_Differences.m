clear; close all; clc;

Lx = 4; 
Ly = 4; 

Nb = 288;

R1 = 0.01

R2 = 1

R3 = 1

k1 = 10

k2_square = R1*k1;

k3 = 0
k3_tri = 10;
k4_tri = R2*k2_square;

gamma = 0.001

data_set_number = 1

%[0,177,2,6]

%[0,177,8001,12001]

colors = [0 0 0;0 1 0;0 0 1;1 0 0]

%colors = lines(4)

state_number_list = [0,177,8001,12001]

for i = 1:4
    
    state_number = state_number_list(i);


directory_name = ['NetworkData/Lx_',num2str(Lx),'_Ly_',num2str(Ly),'/k2_',num2str(k2_square,'%.6f'),'_k3_',num2str(k3_tri,'%.6f'),'_k4_',num2str(k4_tri,'%.6f'),'/g_',num2str(10000*gamma,'%.4f'),'/Run_',num2str(data_set_number),'_State_',num2str(state_number),'/'];

node_filename = [directory_name,'NodeCoords_Lx_',num2str(Lx),'_Ly_',num2str(Ly),'_gamma_',num2str(10000*gamma,'%.3f'),'_k_',num2str(k3,'%.1f'),'_',num2str(data_set_number),'.txt'];
bond_filename = [directory_name,'BondPairs_Lx_',num2str(Lx),'_Ly_',num2str(Ly),'_gamma_',num2str(10000*gamma,'%.3f'),'_k_',num2str(k3,'%.1f'),'_',num2str(data_set_number),'.txt'];

energy_filename = [directory_name,'Energy_Lx_',num2str(Lx),'_Ly_',num2str(Ly),'_gamma_',num2str(10000*gamma,'%.3f'),'_k_',num2str(k3,'%.1f'),'_',num2str(data_set_number),'.txt'];

energy_data = importdata(energy_filename)

time_list = energy_data(:,1);

energy_list = energy_data(:,2);

k2_current = energy_data(:,3);

energy_square_list = energy_data(:,4);

energy_diamond_list = energy_data(:,5);


T_relax = 250

T_step = 20

total_steps = 200

n_relax = ceil(T_relax/gamma)

n_step = ceil(T_step/gamma)

half_step = ceil(n_step/2);

sample_steps = n_relax+half_step:2*half_step:(n_relax + n_step*total_steps - half_step)

index_list = find(ismember(time_list,sample_steps))

figure(1)

subplot(2,1,1);
plot(time_list*gamma,energy_list,'r-','LineWidth',0.5); hold on;

plot(time_list(index_list)*gamma,energy_list(index_list),'^','MarkerSize',5);

set(gca,'Yscale','log');

subplot(2,1,2);
plot(time_list*gamma,k2_current/k1,'b-','LineWidth',2); hold on;

figure(2);

if(i==1)
    y_list_s = energy_square_list(index_list)
    y_list_t = energy_diamond_list(index_list)
    
    y_list_avg = y_list_s + 2*y_list_t;
end


plot(k2_current(index_list)/k1,energy_square_list(index_list)-y_list_s,'-','LineWidth',2,'Color',colors(i,:)); hold on;
plot(k2_current(index_list)/k1,energy_diamond_list(index_list)-y_list_t,':','LineWidth',2,'Color',colors(i,:)); hold on;

set(gca,'Xscale','log');
%set(gca,'Yscale','log');

figure(3);

plot(k2_current(index_list)/k1,(energy_square_list(index_list))./y_list_avg,'-','LineWidth',2,'Color',colors(i,:)); hold on;
plot(k2_current(index_list)/k1,(energy_diamond_list(index_list))./y_list_avg,':','LineWidth',2,'Color',colors(i,:)); hold on;



title('Lx = 4, Ly = 4 Annealed Energy Protocol')

set(gca,'Xscale','log');
set(gca,'Yscale','lin');

figure(4); 

if(i==1)
    plot(k2_current(index_list)/k1,y_list_avg,'-','LineWidth',2,'Color',colors(i,:)); hold on;
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
else
    %plot(k2_current(index_list)/k1,(energy_square_list(index_list))+2*(energy_diamond_list(index_list)),'LineWidth',2,'Color',colors(i,:));
end

end

%plot(0.0252*ones(1,100),linspace(-0.1,0.6,100),'k--'); hold on; 

figure(3)

%ylab = ylabel('$ E^{\framebox(4,4)} / E_0 , E^{\Delta} / E_0 $');
%xlabel('k_2/k_1');

axis([0.01 1 0 2])

%legend('(a) Squares','(a) Triangles','(b) Squares','(b) Triangles','(c) Squares','(c) Triangles','(d) Squares','(d) Triangles')
