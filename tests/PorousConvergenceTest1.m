clc

u=0;
for i = 0.96: -0.1 : 0.06
    u=u+1;
    [DisplacementPorousNode1_96_006(u,:), DisplacementPorousNode2_96_006(u,:), DisplacementPorousNode3_96_006(u,:), ...
        DisplacementPorousNode4_96_006(u,:)] ...
        = ModelTestForPorousConvergence(i);  

end
fprintf("end")