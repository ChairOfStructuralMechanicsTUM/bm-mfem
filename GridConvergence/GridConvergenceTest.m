Lx = 0.05;
Ly = 0.01;
u=0;
for i = 1: 1 : 20
   u=u+1;
   nx = 5*i;
   ny = i;
   [DisplacementPorousNode1(u,:), DisplacementPorousNode2(u,:), DisplacementPorousNode3(u,:), DisplacementPorousNode4(u,:), DisplacementPorousTime(u,:)] = ModelTestForGridCovergence(Lx,Ly,nx,ny);  
end
fprintf("end")