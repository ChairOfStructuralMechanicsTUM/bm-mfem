%%%Test Beam

clear;

node01 = Node(1,0,0);
node02 = Node(2,0.5,0);
node03 = Node(3,0.5,0.25);
node04 = Node(4,0,0.25);

nodeArray = [node01 node02 node03 node04];

nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

ele01 = QuadrilateralElement2d4n(1,[node01 node02 node03 node04]);

elementArray = [ele01];

elementArray.setPropertyValue('YOUNGS_MODULUS',96);
elementArray.setPropertyValue('POISSON_RATIO',1/3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);
elementArray.setPropertyValue('DENSITY',7860);

elementIds = elementArray.getId;

node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');
node04.fixDof('DISPLACEMENT_X');
node04.fixDof('DISPLACEMENT_Y');

addPointLoad(node03,10,[0 -1]);

model = FemModel(nodeArray, elementArray);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
            
massMatrix = assembling.assembleGlobalMassMatrix(model);

solver = SimpleSolvingStrategy(model);
x = solver.solve();

step = 1;

VerschiebungDofs = model.getDofArray.getValue(step);

nodalForces = solver.getNodalForces(step);

v = Visualization(model);
v.setScaling(0.01);
v.plotUndeformed
v.plotDeformed
    
    
