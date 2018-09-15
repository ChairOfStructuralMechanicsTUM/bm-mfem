clc

omega=100;
load=10;
Lx = .05;
Ly = 0.01;

u=15;
for i = 16: 1 : 20
    u=u+1;
    nx = 10*i;
    ny = i;
    [DisplacementPorousNode1(u,:), DisplacementPorousNode2(u,:), DisplacementPorousNode3(u,:), ...
        DisplacementPorousNode4(u,:), DisplacementPorousTime(u,:), DisplacementPorousRatio(u,:)] ...
        = ModelTestForGridCovergence(omega,load,Lx,Ly,nx,ny);  

end
fprintf("end")