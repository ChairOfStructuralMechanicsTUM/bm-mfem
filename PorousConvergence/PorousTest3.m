u=0;
for i = 0.96:-0.1:0.06
    u=u+1;
    EndNode_096_006(u,:) = ModelTestForPorousConvergence(i);
end
fprintf("end")

%%
u=0;
for i = 0.06:-0.01:0.01
    u=u+1;
    EndNode_006_001(u,:) = ModelTestForPorousConvergence(i);
end

%%

u=0;
for i = 0.009:-0.001:0.001
    u=u+1;
    EndNode_0009_0001(u,:) = ModelTestForPorousConvergence(i);
end