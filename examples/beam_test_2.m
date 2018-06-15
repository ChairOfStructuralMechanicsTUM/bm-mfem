clear;

model = FemModel();

for i=1:21
    model.addNewNode(i,(i-1)*0.5,0); 
end
%Erstellen von oberer Knotenreihe mit y=0 und x-Abstand 0,5
for i=22:42
    model.addNewNode(i,(i-22)*0.5,1);
end
%Erstellen von mittlerer Knotenreihe mit y=1 
for i=43:63
    model.addNewNode(i,(i-43)*0.5,2);
end
%Erstellen von oberer Knotenreihe mit y=2 
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

for i=1:20
    model.addNewElement('QuadrilateralElement2d4n',i,[i i+1 i+22 i+21]);
end %Elemente mit ID i, bestehend aus 4 Knoten (obere und mittlere Reihe);innere Knoten jew Doppelt

for i=21:40
    model.addNewElement('QuadrilateralElement2d4n',i,[i+1 i+2 i+23 i+22]);
end %Elemente mit ID i, bestehend aus 4 Knoten (mittlere und untere Reihe);innere Knoten bis zu 4 mal 



model.getNode(1).fixDof('DISPLACEMENT_X');
model.getNode(1).fixDof('DISPLACEMENT_Y');
model.getNode(22).fixDof('DISPLACEMENT_X');
model.getNode(22).fixDof('DISPLACEMENT_Y');
model.getNode(43).fixDof('DISPLACEMENT_X');
model.getNode(43).fixDof('DISPLACEMENT_Y');

model.getAllElements.setPropertyValue('YOUNGS_MODULUS',210000000000);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.3);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',7860);

addPointLoad(model.getNode(63),1,[0 -1]);

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
    
