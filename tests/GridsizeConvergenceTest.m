clc


omega=100;
load=10;
Lx = .05;
Ly = 0.01;

u=0;
for i = 1: 1 : 20
    u=u+1;
    nx = 5*i;
    ny = i;
    [DisplacementPorousNodeA(u,:), DisplacementPorousNodeB(u,:), DisplacementPorousNodeC(u,:), ...
        DisplacementPorousNodeD(u,:), DisplacementPorousTime(u,:)] ...
        = ModelTestForGridCovergence(omega,load,Lx,Ly,nx,ny);  
end
fprintf("end")