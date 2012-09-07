% written by: Stu Blair
% Date: 11/3/11
% Purpose: implement 6-node triangular isoperimetric shape functions
% based on T.J.R. Hughes p167,168.

function [shapet6,dhdrt6,dhdst6]=feisot6(r,s)


% shape functions
shapet6(1)=2*r*r-r;
shapet6(2)=2*s*s-s;
shapet6(3)=1-3*r-3*s+4*r*s+2*r*r+2*s*s;
shapet6(4)=4*r*s;
shapet6(5)=4*s-4*r*s-4*s*s;
shapet6(6)=4*r-4*r*r-4*r*s;

% derivatives
dhdrt6(1)=4*r-1;
dhdrt6(2)=0;
dhdrt6(3)=-3+4*s+4*r;
dhdrt6(4)=4*s;
dhdrt6(5)=-4*s;
dhdrt6(6)=4-8*r-4*s;

dhdst6(1)=0;
dhdst6(2)=4*s-1;
dhdst6(3)=-3+4*r+4*s;
dhdst6(4)=4*r;
dhdst6(5)=4-4*r-8*s;
dhdst6(6)=-4*r;