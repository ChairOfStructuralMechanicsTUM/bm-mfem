clear all
clc

Lx = .05;
Ly = Lx/5;
u=15;

for i = 16: 1 : 20
    u=u+1;
    nx = 5*i;
    ny = i;
    DisplacementPorousEndNode(u,:) = ModelTestForGridCovergence(Lx,Ly,nx,ny);
    
end
fprintf("end")