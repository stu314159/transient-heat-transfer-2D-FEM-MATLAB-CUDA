function [K]= makeK_2D_scalar_tri_isot6(gcoord,nodes,k_vec)

% gcoord - global coordinate array
% node - local to global node numbering for elements in domain
% B = temperature kinematic array 2 x nnode


[nnode,~]=size(gcoord);
[nel,nnel]=size(nodes);

ndof = 1;
edof = ndof*nnel;
xi_coord = [1/6 2/3 1/6];
eta_coord = [1/6 1/6 2/3];
weight2 = [1/6 1/6 1/6];

nzmax = nnode*20;
K_vec = zeros(nzmax*3,3);

B = zeros(2,nnel);
num_entries = 0;

for iel = 1:nel
    nd = nodes(iel,:); % get nodes for iel element
    xcoord = zeros(nnel,1); ycoord = zeros(nnel,1);
    index = feeldof(nd,nnel,ndof); % get global DOFs for nodes
    
    D = eye(2,2)*k_vec(iel);
    
    for i = 1:nnel
        xcoord(i)=gcoord(nd(i),1);ycoord(i)=gcoord(nd(i),2);
    end
    
    for int_p = 1:3
        
        xi = xi_coord(int_p);
        eta = eta_coord(int_p);
        wt = weight2(int_p);
        [~,dhdr,dhds]=feisot6(xi,eta);
        jacob2=fejacob2(nnel,dhdr,dhds,xcoord,ycoord);
        detjacob=det(jacob2);
        invjacob=inv(jacob2);
        [dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob);
        
        % form elemental B matrix
        for i = 1:edof
            B(1,i)=dhdx(i);
            B(2,i)=dhdy(i);
        end
        
        % get elemental K matrix for this integration point
        K_ele = B'*D*B*wt*detjacob;
        
        % assemble into the global K
        for i = 1:edof
            ii = index(i);
            for j = 1:edof
                num_entries = num_entries+1;
                jj=index(j);
                K_vec(num_entries,1:3)=[ii,jj,K_ele(i,j)];
            end
        end
        
    end % int_p
end % iel

K_vec = K_vec(1:num_entries,:);

% this operation will add up all non-zero entries associated
% with each dof.
K = sparse(K_vec(:,1),K_vec(:,2),K_vec(:,3),nnode,nnode);

