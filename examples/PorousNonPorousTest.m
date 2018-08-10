
u=0
for i=0.95:-0.1:0.05
    u=u+1
    Diff(u,:) = ModelTest(i)
end

