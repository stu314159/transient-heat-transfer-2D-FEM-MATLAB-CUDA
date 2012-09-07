function Kh = convectionBoundary1Disol3(gcoord,nodes,h_vec)

[nnode,~]=size(gcoord);
[nel,nnel]=size(nodes);
ndof = 1;
edof = ndof*nnel;
ngl = 2;
[point1,weight1]=feglqd1(ngl);
nzmax = nel*15;
M_vec = zeros(nzmax*2,3);

num_entries = 0;
for iel=1:nel
    nd = nodes(iel,:); % get nodes for iel element
    xcoord = zeros(nnel,1);
    ycoord = zeros(nnel,1);
    index = feeldof(nd,nnel,ndof);
    
    for i = 1:nnel
        xcoord(i)=gcoord(nd(i),1); ycoord(i)=gcoord(nd(i),2);
    end
    length_ele = sqrt((xcoord(end)-xcoord(1))^2+...
        (ycoord(end)-ycoord(1))^2);
    fake_length = zeros(3,1);
    fake_length(1) = -length_ele/2; fake_length(3)=length_ele/2;
    
    for intx = 1:ngl
        x = point1(intx);
        wtx = weight1(intx);
        [shape,dhdr]=feisol3(x);
        jacob1=fejacob1(nnel,dhdr,fake_length);
        detjacob = jacob1;
        %invjacob=1/jacob1;
        %dhdx=invjacob*dhdr;
        for i = 1:edof
            ii=index(i);
            for j = 1:edof
                jj = index(j);
                num_entries=num_entries+1;
                M_vec(num_entries,1:3)=[ii,jj,h_vec(iel)*shape(i)*shape(j)*wtx*detjacob];
            end
        end
      
    end % intx
    
end % iel

M_vec=M_vec(1:num_entries,:);
Kh = sparse(M_vec(:,1),M_vec(:,2),M_vec(:,3),nnode,nnode);