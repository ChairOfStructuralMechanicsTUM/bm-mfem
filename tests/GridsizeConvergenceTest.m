clc

omega=10;
load=10;
Lx = .05;
Ly = 0.01;

u=0;
for i = 1: 1 : 20
    u=u+1;
    nx = 5*i;
    ny = i;
    DisplacementPorousEndNodeOmega100_F_10(u,:) = ModelTestForGridCovergence(omega,load,Lx,Ly,nx,ny);  

end
fprintf("end")