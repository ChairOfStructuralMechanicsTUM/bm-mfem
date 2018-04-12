clear;

model = FemModel();

for i=1:11
    model.addNewNode(i,(i-1)*0.2,0);
end

for i=12:22
    model.addNewNode(i,(i-1)*0.2,2);
end

model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

for i=1:10
    model.addNewElement('QuadrilateralElement2d4n',i,[i i+1 i+12 i+11]);
end

model.getNode(1).fixDof('DISPLACEMENT_X');
model.getNode(1).fixDof('DISPLACEMENT_Y');
model.getNode(12).fixDof('DISPLACEMENT_X');
model.getNode(12).fixDof('DISPLACEMENT_Y');

model.getAllElements.setPropertyValue('YOUNGS_MODULUS',96);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.3);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',7860);

addPointLoad(model.getNode(22),1,[0 -1]);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
            
massMatrix = assembling.assembleGlobalMassMatrix(model);

solver = SimpleSolvingStrategy(model);
x = solver.solve();

step = 1;

VerschiebungDofs = model.getDofArray.getValue(step);

nodalForces = solver.getNodalForces(step);

fig=figure;
subplot(2,2,1);
plot(0:0.2:2,VerschiebungDofs(2:2:22));
title('w_unten', 'FontSize', 8);

subplot(2,2,2);
plot(0:0.2:2,VerschiebungDofs(24:2:end));
title('w_oben', 'FontSize', 8);

subplot(2,2,3);
plot(0:0.2:2,VerschiebungDofs(1:2:21));
title('u_unten', 'FontSize', 8);

subplot(2,2,4);
plot(0:0.2:2,VerschiebungDofs(23:2:end));
title('u_oben', 'FontSize', 8);
    
