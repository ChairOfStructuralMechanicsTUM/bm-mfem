node01 = Node(1,0,0);
node02 = Node(2,2,0);
node03 = Node(3,3,0);
node04 = Node(4,4,0);
node05 = Node(5,6,0);
node06 = Node(6,2,0.5);
node07 = Node(7,3,0.5);
node08 = Node(8,4,0.5);
node09 = Node(9,0,1);
node10 = Node(10,2,1);
node11 = Node(11,3,1);
node12 = Node(12,4,1);
node13 = Node(13,6,1);
node14 = Node(14,0,2);
node15 = Node(15,2,2);
node16 = Node(16,4,2);
node17 = Node(17,6,2);

nodeArray = [node01 node02 node03 node04 node05 node06 node07 node08 ...
    node09 node10 node11 node12 node13 node14 node15 node16 node17];

nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

ele01 = QuadrilateralElement2d4n(1,[node01 node02 node10 node09]);
ele02 = QuadrilateralElement2d4n(2,[node02 node03 node07 node06]);
ele03 = QuadrilateralElement2d4n(3,[node03 node04 node08 node07]);
ele04 = QuadrilateralElement2d4n(4,[node06 node07 node11 node10]);
ele05 = QuadrilateralElement2d4n(5,[node07 node08 node12 node11]);
ele06 = QuadrilateralElement2d4n(6,[node04 node05 node13 node12]);
ele07 = QuadrilateralElement2d4n(7,[node09 node10 node15 node14]);
ele08 = QuadrilateralElement2d4n(8,[node10 node12 node16 node15]);
ele09 = QuadrilateralElement2d4n(9,[node12 node13 node17 node16]);



elementArray = [ele01 ele02 ele03 ele04 ele05 ele06 ele07 ele08 ele09];

elementArray.setPropertyValue('YOUNGS_MODULUS',96);
elementArray.setPropertyValue('POISSON_RATIO',1/3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);
elementArray.setPropertyValue('DENSITY',7860);


model = FemModel(nodeArray, elementArray);

obj = BlochInverse1D(model);
assembling = SimpleAssembler(model);

stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
            
massMatrix = assembling.assembleGlobalMassMatrix(model);


initialize(obj);
