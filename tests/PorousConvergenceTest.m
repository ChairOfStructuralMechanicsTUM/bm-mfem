clc

u=0;
for i = 0.06: -0.01 : 0.01
    u=u+1;
    [DisplacementPorousNode1(u,:), DisplacementPorousNode2(u,:), DisplacementPorousNode3(u,:), ...
        DisplacementPorousNode4(u,:)] ...
        = ModelTestForPorousConvergence(i);  

end
fprintf("end")