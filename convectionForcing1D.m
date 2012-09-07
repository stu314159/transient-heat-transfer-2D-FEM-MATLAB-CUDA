% 
function f_h = convectionForcing1D(gcoord,nodes,h_vec,T_infty)

[nnode,~]=size(gcoord);
[nel,nnel]=size(nodes);
ndof = 1;
edof = ndof*nnel;
ngl = 2;
[point1,weight1]=feglqd1(ngl);

f_h = zeros(nnode,1);

for iel = 1:nel
    nd = nodes(iel,:);
    xcoord=zeros(nnel,1);
    ycoord = zeros(nnel,1);
    index = feeldof(nd,nnel,ndof);
    for i = 1:nnel
        xcoord(i)=gcoord(nd(i),1);ycoord(i)=gcoord(nd(i),2);
    end
    length_ele=sqrt((xcoord(end)-xcoord(1))^2+...
        (ycoord(end)-ycoord(1))^2);
    fake_length=zeros(2,1);
    fake_length(1)=-length_ele/2; fake_length(2)=length_ele/2;
    
    for intx=1:ngl
        x=point1(intx);
        wtx=weight1(intx);
        [shape,dhdr]=feisol2(x);
        jacob1=fejacob1(nnel,dhdr,fake_length);
        detjacob=jacob1;
        for i = 1:edof
            ii = index(i);
            f_h(ii)= f_h(ii)+T_infty(iel)*h_vec(iel)*shape(i)*wtx*detjacob;
        end
    end
end