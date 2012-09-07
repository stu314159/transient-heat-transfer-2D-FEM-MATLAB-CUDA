function [B]= makeB_2D_scalar_tri6(gcoord,nodes)

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

B = zeros(2,nnode);


for iel = 1:nel
    nd = nodes(iel,:); % get nodes for iel element
    xcoord = zeros(nnel,1); ycoord = zeros(nnel,1);
    index = feeldof(nd,nnel,ndof); % get global DOFs for nodes
    
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
                
        for i = 1:edof
            ii = index(i);
            B(1,ii)=dhdx(i)*wt*detjacob;
            B(2,ii)=dhdy(i)*wt*detjacob;
            
           
        end
        
    end % int_p
end % iel

