function C = makeC_2D_tri_r2(gcoord,nodes,RhoC,lumped)

[nnode,~] = size(gcoord);
[nel,nnel]=size(nodes);
ndof = 1;
edof = ndof*nnel;
xi_coord = [1/6 2/3 1/6];
eta_coord = [1/6 1/6 2/3];
weight2 = [1/6 1/6 1/6];

nzmax=nnode*20;
C_vec = zeros(nzmax*3,3);

num_entries = 0;

for iel = 1:nel
    nd = nodes(iel,:); % get nodes for iel element
    xcoord = zeros(nnel,1); ycoord = zeros(nnel,1);
    index = feeldof(nd,nnel,ndof); % get global DOFs for nodes
    
    
    for i = 1:nnel
        xcoord(i) = gcoord(nd(i),1); ycoord(i)=gcoord(nd(i),2);
    end % get nodal coordinates
    
    for int_p = 1:3
        xi = xi_coord(int_p);
        eta = eta_coord(int_p);
        wt = weight2(int_p);
        [shape,dhdr,dhds]=feisot3(xi,eta);
        jacob2 = fejacob2(nnel,dhdr,dhds,xcoord,ycoord);
        detjacob = det(jacob2);
        %invjacob = inv(jacob2);
        %[dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob);
        for i = 1:edof
            ii = index(i);
            for  j = 1:edof
                jj = index(j);
                num_entries = num_entries+1;
                C_vec(num_entries,1:3) = [ii,jj,RhoC(iel)*shape(i)*shape(j)*wt*detjacob];
              
            end %j
        end % i
    end % int_p
      
end % iel

C_vec = C_vec(1:num_entries,:);
C = sparse(C_vec(:,1),C_vec(:,2),C_vec(:,3),nnode,nnode);

if lumped==1
    alpha = sum(diag(C));
    tot_mass = sum(sum(C));
    C = diag(C).*(tot_mass/alpha);
end
