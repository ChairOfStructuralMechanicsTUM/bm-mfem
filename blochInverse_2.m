clear;

node01 = Node(1,0,0);
node02 = Node(2,2,0);
node03 = Node(3,4,0);
node04 = Node(4,6,0);
node05 = Node(5,0,1);
node06 = Node(6,2,1);
node07 = Node(7,4,1);
node08 = Node(8,6,1);
node09 = Node(9,0,2);
node10 = Node(10,2,2);
node11 = Node(11,4,2);
node12 = Node(12,6,2);

nodeArray = [node01 node02 node03 node04 node05 ...
    node06 node07 node08 node09 node10 node11 node12];

nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});


ele01 = QuadrilateralElement2d4n(1,[node01 node02 node06 node05]);
ele02 = QuadrilateralElement2d4n(2,[node02 node03 node07 node06]);
ele03 = QuadrilateralElement2d4n(3,[node03 node04 node08 node07]);
ele04 = QuadrilateralElement2d4n(4,[node05 node06 node10 node09]);
ele05 = QuadrilateralElement2d4n(5,[node06 node07 node11 node10]);
ele06 = QuadrilateralElement2d4n(6,[node07 node08 node12 node11]);




elementArray = [ele01 ele02 ele03 ele04 ele05 ele06];

elementArray.setPropertyValue('YOUNGS_MODULUS',96);
elementArray.setPropertyValue('POISSON_RATIO',1/3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);
elementArray.setPropertyValue('DENSITY',7860);


model = FemModel(nodeArray, elementArray);

obj = BlochInverse1D(model);
assembling = SimpleAssembler(model);

stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
            
massMatrix = assembling.assembleGlobalMassMatrix(model);


initialize(obj)


%solver = SimpleSolvingStrategy(model);
%x = solver.solve();
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
