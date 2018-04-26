clear;

model = FemModel();

for i=1:21
    model.addNewNode(i,(i-1)*0.5,0);
end

for i=22:42
    model.addNewNode(i,(i-22)*0.5,0.25);
end

model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

for i=1:20
    model.addNewElement('QuadrilateralElement2d4n',i,[i i+1 i+22 i+21]);
end

model.getNode(1).fixDof('DISPLACEMENT_X');
model.getNode(1).fixDof('DISPLACEMENT_Y');
model.getNode(22).fixDof('DISPLACEMENT_X');
model.getNode(22).fixDof('DISPLACEMENT_Y');

model.getAllElements.setPropertyValue('YOUNGS_MODULUS',96);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.3);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',7860);

addPointLoad(model.getNode(42),1,[0 -1]);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
            
massMatrix = assembling.assembleGlobalMassMatrix(model);

solver = SimpleSolvingStrategy(model);
x = solver.solve();

step = 1;

VerschiebungDofs = model.getDofArray.getValue(step);

nodalForces = solver.getNodalForces(step);

v = Visualization(model);
v.plotUndeformed
v.plotDeformed

% fig=figure;
% subplot(2,2,1);
% plot(0:0.5:10,VerschiebungDofs(2:2:42));
% title('w_unten', 'FontSize', 8);
% 
% subplot(2,2,2);
% plot(0:0.5:10,VerschiebungDofs(44:2:end));
% title('w_oben', 'FontSize', 8);
% 
% subplot(2,2,3);
% plot(0:0.5:10,VerschiebungDofs(1:2:41));
% title('u_unten', 'FontSize', 8);
% 
% subplot(2,2,4);
% plot(0:0.5:10,VerschiebungDofs(43:2:end));
% title('u_oben', 'FontSize', 8);
    
