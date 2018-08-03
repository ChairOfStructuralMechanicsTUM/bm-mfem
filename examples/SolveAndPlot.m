function SolveAndPlot(SolvingMethod, model)

if SolvingMethod == 1
    
    solver = SimpleHarmonicSolvingStrategy(model,100);
    x = solver.solve();
    
    step = 1;
    
    VerschiebungDofs_biot = model.getDofArray.getValue(step);
    nodalForces_biot = solver.getNodalForces(step);
    
    v = Visualization(model);
    v.setScaling(1e3);
    v.plotUndeformed
    v.plotDeformed
end


if SolvingMethod == 2
    
    v = Visualization(model);
    v.setScaling(1e3);
    v = Visualization(model);
    v.plotUndeformed
    dt = .05;
    time = 0;
    endTime = 1.5;
    ls  = linspace(time,endTime,endTime/dt+1);
    solver = NewmarkSolvingStrategy(model, dt);
    
    x = solver.solve();
    step = 1;
    VerschiebungDofs_biot = model.getDofArray.getValue(step);
    
    
    for j = 1:30
        solver.solve();
        time = j*0.05;
        hold on
        v.plotDeformed();
        %view(3)
        axis([-1.1 10 -1.1 1])
        F(j) = getframe(gcf);
    end
    fig = figure;
    movie(fig,F,1,2)
    

end

end