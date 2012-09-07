% written by: Stu Blair
% Date: 11/3/11
% Purpose: isoparametric 1D 2nd order line element based on Kwon pg 162.

function [shape,dhdr]=feisol3(r)

shape(1)=0.5*(r*r-r);
shape(2)=1-r*r;
shape(3)=0.5*(r*r+r);

dhdr(1)=r-1;
dhdr(2)=-2*r;
dhdr(3)=r+1;