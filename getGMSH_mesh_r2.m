function [gcoord,nodes_tri,nodes_lin]=getGMSH_mesh_r2(mesh_file_name)

mesh = load_gmsh(mesh_file_name);

gcoord = mesh.POS(:,1:2); 
nodes_tri = mesh.TRIANGLES(1:mesh.nbTriangles,:);
nodes_lin = mesh.LINES(1:mesh.nbLines,:);