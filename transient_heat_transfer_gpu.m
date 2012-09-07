%transient_heat_transfer.m

clear;
clc;
close('all');

Num_ts = 1500000;
dt = 1.2e-2;
time_step_rep_freq = 1000;
plot_frequency = 1000;

lumped = true;


gpu_accel=4;
%1 = jacket single precision
%2 = jacket double precision
%3 = gpuArray
%4 = nothing

time_int = 1;
% 1 = first order explicit Euler
% 2 = 2nd order Explicit RK
% 3 = 3rd order Explicit RK

% geometric definition
x_left = 0; x_right = 5;
y_bot = 0; y_top = 3;

% circular hole in the middle of the rectangle.
x_c = 2.5; y_c = 1.5;
r = 0.5;

% hot_spots
x_chs1=1; y_chs1=1.5;
r_hs1=.1;

x_chs2=4; y_chs2=1.5;
r_hs2=.1;

T_outside = 400;
T_hole = 300;
T_init = T_outside;

k_mat = 50e0; % W/m-C <---- thermal conductivity for internal conduction
h_boundary_outside = 12e5; % W/m^2-K <----heat transfer coefficient for convective boundaries.
h_boundary_inside = 1e3;
G = 2.2e5; % W/m^3 <---- internal heat generation rate
heat_capacity = 2.42e6; % W/m^3-C


% region 12 = bottom
% region 13 = right
% region 14 = top
% region 15 = left

% regions 16-19 = surface of the hole

% region 20 = whole domain...

mesh_file_name = 'prob_dom2.msh';
%[gcoord,nodes_tri,nodes_lin]=getGMSH_mesh_r2(mesh_file_name);

msh = load_gmsh4(mesh_file_name,[8 9 15]);
gcoord = msh.POS;
nodes_lin3=msh.LINES3;
nodes_tri6=msh.TRIANGLES6;


[nnodes,~]=size(gcoord);
[nel_tri,~]=size(nodes_tri6);
% now I know the nodes and element info: assemble the system matrices.

% get the centroids of the triangular elements
x_c_tri=zeros(nel_tri,1);
y_c_tri=zeros(nel_tri,1);
for el=1:nel_tri
    el_nodes=nodes_tri6(el,1:6);
   x_c_tri(el)=sum(gcoord(el_nodes,1))/6;
   y_c_tri(el)=sum(gcoord(el_nodes,2))/6;
end


% find nodes that are within the hot-spot regions
nds_hs1=find(((x_c_tri-x_chs1).^2 + (y_c_tri-y_chs1).^2)<=r_hs1^2);
nds_hs2=find(((x_c_tri-x_chs2).^2 + (y_c_tri-y_chs2).^2)<=r_hs2^2);



% form the B matrix

% all of the triangular nodes that are part of region 20 are

surface_elements = find(nodes_tri6(:,7)==20);
nel=length(surface_elements);

k_vec = k_mat*ones(nel,1);
k_vec(nds_hs1)=0.1*k_mat;
k_vec(nds_hs2)=0.1*k_mat;
K = makeK_2D_scalar_tri_isot6(gcoord,nodes_tri6(surface_elements,1:6),k_vec);



% get lists of elements on each boundary
bottom_edge = find(nodes_lin3(:,4)==12);
right_edge = find(nodes_lin3(:,4)==13);
top_edge = find(nodes_lin3(:,4)==14);
left_edge = find(nodes_lin3(:,4)==15);

circle = find((nodes_lin3(:,4)==16)|(nodes_lin3(:,4)==17)|...
    (nodes_lin3(:,4)==18)|(nodes_lin3(:,4)==19));

% need to compute the contribution to K from each boundary 
nel = length(bottom_edge);
Kh_bottom=convectionBoundary1Disol3(gcoord,nodes_lin3(bottom_edge,1:3),h_boundary_outside*ones(nel,1));
nel = length(right_edge);
Kh_right=convectionBoundary1Disol3(gcoord,nodes_lin3(right_edge,1:3),h_boundary_outside*ones(nel,1));
nel = length(top_edge);
Kh_top = convectionBoundary1Disol3(gcoord,nodes_lin3(top_edge,1:3),h_boundary_outside*ones(nel,1));
nel = length(left_edge);
Kh_left = convectionBoundary1Disol3(gcoord,nodes_lin3(left_edge,1:3),h_boundary_outside*ones(nel,1));
nel = length(circle);
Kh_circle = convectionBoundary1Disol3(gcoord,nodes_lin3(circle,1:3),h_boundary_inside*ones(nel,1));

K = K + Kh_bottom + Kh_right + Kh_top + Kh_left + Kh_circle;

% get contributions to the RHS from the convection boundaries
% first, lump all the edges
all_outside_edges = [bottom_edge;right_edge;top_edge;left_edge];
nel = length(all_outside_edges);
f_h_outside = convectionForcing1Disol3(gcoord,nodes_lin3(all_outside_edges,1:3),...
    h_boundary_outside*ones(nel,1),T_outside*ones(nel,1));
nel = length(circle);
f_h_inside = convectionForcing1Disol3(gcoord,nodes_lin3(circle,1:3),...
    h_boundary_inside*ones(nel,1),T_hole*ones(nel,1));

% get contributions to the RHS from the internal Source term...
G_vec = G*ones(nel_tri,1);
G_vec(nds_hs1)=1.5*G;
G_vec(nds_hs2)=1.5*G;


f_h = sourceForcing2Dtri_isot6(gcoord,nodes_tri6(:,1:6),G_vec);

f = f_h + f_h_inside + f_h_outside;
% get the heat mass-matrix
C = makeC_2D_tri_r2_isot6(gcoord,nodes_tri6(:,1:6),heat_capacity*ones(nel_tri,1),lumped);

% C*dT/dt + K*T = f

% all boundary conditions and forcing terms applied....let's iterate!

T = T_init*ones(nnodes,1);

switch lumped
    case true
        K = diag(C)\K;
        f = diag(C)\f;
    case false
        K = C\K;
        f = C\f;
end


switch gpu_accel
    case 1
        T = gsingle(T);
        K = gsingle(K);
        f = gsingle(f);
        T_new = gsingle(T);
    case 2
        T = gdouble(T);
        K = gdouble(K);
        f = gdouble(f);
        T_new = gdouble(T);
        
    case 3
        T = gpuArray(T);
        K = gpuArray(full(K));
        f = gpuArray(f);
        T_new = gpuArray(T);
    case 4
end


tic
for ts = 1:Num_ts
    
    if((ts ~=1)&&(mod(ts,time_step_rep_freq)==0))
        fprintf('Executing time step %g. \n', ts);
    end
    switch time_int
        
        case 1
            T_new = T-dt*(K*T)+dt*f;
            T = T_new;
        case 2
            
        case 3
            q0=T;
            q1=q0+dt*(f-K*q0);
            q2=(3/4)*q0+(1/4)*q1+(dt/4)*(f-K*q1);
            q3=(1/3)*q0+(2/3)*q2+((2*dt)/3)*(f-K*q2);
            T=q3;
    end
    
    if(mod(ts,plot_frequency) == 0)
        % do some sort of visualization
        scatter3(gcoord(:,1),gcoord(:,2),T,25,T,'filled')
        %imagesc(gcoord(:,1),gcoord(:,2),T)
        drawnow
    end
end
iteration_time = toc;
fprintf('Node updates per second = %g.\n',nnodes*Num_ts/iteration_time);


figure(2)
switch gpu_accel
    
    case 1
        T = double(T);
        
    case 2
        T = double(T);
        
    case 3
        T = gather(T);
        
end

scatter3(gcoord(:,1),gcoord(:,2),T,25,T,'filled')
drawnow


