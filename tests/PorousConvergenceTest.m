clc

u=0;
for i = 0.96: -0.1 : 0.06
    u=u+1;
    [DisplacementPorousNode1(u,:), DisplacementPorousNode2(u,:), DisplacementPorousNode3(u,:), ...
        DisplacementPorousNode4(u,:)] ...
        = ModelTestForPorousConvergence(i);  

end
fprintf("end")