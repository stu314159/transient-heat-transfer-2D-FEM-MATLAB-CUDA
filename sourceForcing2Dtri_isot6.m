function f_h = sourceForcing2Dtri_isot6(gcoord,nodes,G)

[nnode,~]=size(gcoord);
[nel,nnel]=size(nodes);
ndof=1;
edof = ndof*nnel;

xi_coord = [1/6 2/3 1/6];
eta_coord = [1/6 1/6 2/3];
weight2 = [1/6 1/6 1/6];

f_h = zeros(nnode,1);

for iel = 1:nel
    nd = nodes(iel,:);
    xcoord = zeros(nnel,1);ycoord=zeros(nnel,1);
    index = feeldof(nd,nnel,ndof);
     for i = 1:nnel
        xcoord(i) = gcoord(nd(i),1); ycoord(i)=gcoord(nd(i),2);
    end % get nodal coordinates
    
    for int_p = 1:3
        xi = xi_coord(int_p);
        eta = eta_coord(int_p);
        wtx = weight2(int_p);
        [shape,dhdr,dhds]=feisot6(xi,eta);
        jacob2=fejacob2(nnel,dhdr,dhds,xcoord,ycoord);
        detjacob=det(jacob2);
        %invjacob=inv(jacob2);
        for i = 1:edof
            ii = index(i);
            f_h(ii) = f_h(ii) + G(iel)*shape(i)*wtx*detjacob;
            %f_h(ii)=f_h(ii)+G(iel)*detjacob*shape(i)*wtx*detjacob;
        end %i
        
    end % int_p
        
    
    
end %iel