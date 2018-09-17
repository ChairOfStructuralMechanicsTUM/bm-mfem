clc

u=0;
for i = 0.0000000009: -0.0000000001 : 0.0000000001
    u=u+1;
    [DisplacementPorousNode1_9e10_1e10(u,:), DisplacementPorousNode2_9e10_1e10(u,:), DisplacementPorousNode3_9e10_1e10(u,:), ...
        DisplacementPorousNode4_9e10_1e10(u,:)] ...
        = ModelTestForPorousConvergence(i);  

end
fprintf("end")