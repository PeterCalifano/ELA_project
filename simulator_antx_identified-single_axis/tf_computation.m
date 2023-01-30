% Start by clearing everything
clear 
clc


syms Xu Xq Mu Mq Xd Md s


A=[Xu Xq -9.81; Mu Mq 0; 0 1 0];
B=[Xd; Md; 0;];
C=[0 1 0; Xu Xq 0];
D=[0;Xd];


Phi=inv(s*eye(3)-A);
 

H=C*Phi*B+D;
 
 
%Display
pretty(H)