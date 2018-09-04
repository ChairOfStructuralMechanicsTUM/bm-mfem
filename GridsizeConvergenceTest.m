clear all
clc

Lx = .1;
Ly = 0.02;
u=0;

for i = 1: 1 : 15
    u=u+1;
    nx = Lx/Ly*i;
    ny = i;
    minLambda(u,1) = ModelTestForGridCovergence(Lx,Ly,nx,ny);
    minLambda(u,2) = real(minLambda(u,1))/max([Lx/nx,Ly/ny]);
end