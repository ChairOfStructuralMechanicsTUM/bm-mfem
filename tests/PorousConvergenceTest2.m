clc

u=0;
for i = 0.009: -0.001 : 0.001
    u=u+1;
    [DisplacementPorousNode1_009_001(u,:), DisplacementPorousNode2_009_001(u,:), DisplacementPorousNode3_009_001(u,:), ...
        DisplacementPorousNode4_009_001(u,:)] ...
        = ModelTestForPorousConvergence(i);  

end
fprintf("end")