clear;

model = FemModel();

for i=1:3
    model.addNewNode(i,(i-1)*0.5,0);
end

for i=4:6
    model.addNewNode(i,(i-4)*0.5,0.5);
end

for i=7:9
    model.addNewNode(i,(i-7)*0.5,1);
end

model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

for i=1:2
    model.addNewElement('QuadrilateralElement2d4n',i,[i i+1 i+4 i+3]);
end

for i=3:4
    model.addNewElement('QuadrilateralElement2d4n',i,[i+1 i+2 i+5 i+4]);
end

model.getNode(1).fixDof('DISPLACEMENT_X');
model.getNode(1).fixDof('DISPLACEMENT_Y');
model.getNode(7).fixDof('DISPLACEMENT_X');
model.getNode(7).fixDof('DISPLACEMENT_Y');
% model.getNode(21).fixDof('DISPLACEMENT_X');
% model.getNode(21).fixDof('DISPLACEMENT_Y');
% model.getNode(42).fixDof('DISPLACEMENT_X');
% model.getNode(42).fixDof('DISPLACEMENT_Y');

model.getAllElements.setPropertyValue('YOUNGS_MODULUS',210000000000);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.3);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',3);
model.getAllElements.setPropertyValue('DENSITY',7860);

addPointLoad(model.getNode(9),1,[0 -1]);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
            
massMatrix = assembling.assembleGlobalMassMatrix(model);

solver = SimpleSolvingStrategy(model);
x = solver.solve();

step = 1;

VerschiebungDofs = model.getDofArray.getValue(step);

nodalForces = solver.getNodalForces(step);

v = Visualization(model);
v.setScaling(1000000000);
v.plotUndeformed
v.plotDeformed
    
