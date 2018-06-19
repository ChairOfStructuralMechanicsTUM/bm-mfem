%Quadril.Element2d6n
clear;

model = FemModel();

for i=1:3
    model.addNewNode(i,(i-1)*0.5,0); 
end
%Erstellen von oberer Knotenreihe mit y=0 und x-Abstand 0,5
for i=4:6
    model.addNewNode(i,(i-4)*0.5,1);
end

model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
%Benennen der x und y Verschiebungen;

for i=1:2
    model.addNewElement('QuadrilateralElement2d4n',i,[i i+1 i+4 i+3]);
end %Elemente mit ID i, bestehend aus 4 Knoten (obere und mittlere Reihe);innere Knoten jew Doppelt
%2 Elemente 

model.getNode(1).fixDof('DISPLACEMENT_X'); %?
model.getNode(1).fixDof('DISPLACEMENT_Y');
model.getNode(4).fixDof('DISPLACEMENT_X');
model.getNode(4).fixDof('DISPLACEMENT_Y');


model.getAllElements.setPropertyValue('YOUNGS_MODULUS',210000000000);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.3);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',7860);

addPointLoad(model.getNode(6),1,[0 -1]);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
            
massMatrix = assembling.assembleGlobalMassMatrix(model);

solver = SimpleSolvingStrategy(model);
x = solver.solve();

step = 1;

VerschiebungDofs = model.getDofArray.getValue(step);

nodalForces = solver.getNodalForces(step);

v = Visualization(model);
v.setScaling(10000000);
v.plotUndeformed
v.plotDeformed






% elementArray.setPropertyValue('YOUNGS_MODULUS',96);
% elementArray.setPropertyValue('POISSON_RATIO',1/3);
% elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);
% elementArray.setPropertyValue('DENSITY',7860);
% 
% node01.fixDof('DISPLACEMENT_X');
% node01.fixDof('DISPLACEMENT_Y');
% node04.fixDof('DISPLACEMENT_X');
% node04.fixDof('DISPLACEMENT_Y');
% 
% addPointLoad(node03,1,[0 -1]);
% 
% model = FemModel(nodeArray, elementArray);
% 
% assembling = SimpleAssembler(model);
% stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
%             
% massMatrix = assembling.assembleGlobalMassMatrix(model);
% 
% solver = SimpleSolvingStrategy(model);
% x = solver.solve();
% 
% step = 1;
% 
% VerschiebungDofs = model.getDofArray.getValue(step);
% 
% nodalForces = solver.getNodalForces(step);
% 
% v = Visualization(model);
% v.setScaling(1);
% v.plotUndeformed
% v.plotDeformed