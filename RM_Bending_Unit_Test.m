% example taken from Introduction to FEM (Felippa) 23-5
a=5; 
b=2*a;
node01 = Node(1,0,0,0);
node02 = Node(2,a,0,0);
node03 = Node(3,a,b,0);
node04 = Node(4,0,b,0);

nodeArray = [node01 node02 node03 node04]; 
nodeArray.addDof({'DISPLACEMENT_Z', 'ROTATION_X', 'ROTATION_Y'});
ele = ReissnerMindlinElement3d4n(1,nodeArray); 

elementArray = ele;
elementArray.setPropertyValue('THICKNESS', 1);
elementArray.setPropertyValue('YOUNGS_MODULUS', 8000);
elementArray.setPropertyValue('SHEAR_MODULUS', 3000);
elementArray.setPropertyValue('POISSON_RATIO', 1/3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT', 2);
elementArray.setPropertyValue('DENSITY', 1);

model = FemModel(nodeArray,elementArray);

[K,~] = SimpleAssembler.assembleGlobalStiffnessMatrix(model);
M = SimpleAssembler.assembleGlobalMassMatrix(model);