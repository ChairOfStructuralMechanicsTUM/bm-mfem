clear;

node01 = Node(1,-1,-1);
node02 = Node(2,1,-1);
node03 = Node(3,1,1);
node04 = Node(4,-1,1);



nodeArray = [node01 node02 node03 node04];

nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

ele01 = QuadrilateralElement2d4n(1,[node01 node02 node03 node04]);

elementArray = [ele01];

elementArray.setPropertyValue('YOUNGS_MODULUS',32);
elementArray.setPropertyValue('POISSON_RATIO',1/3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);
elementArray.setPropertyValue('DENSITY',7860);

elementIds = elementArray.getId;

model = FemModel(nodeArray, elementArray);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);

massMatrix = assembling.assembleGlobalMassMatrix(model);
    
    
