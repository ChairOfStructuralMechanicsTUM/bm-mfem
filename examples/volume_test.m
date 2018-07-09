clear;

node01 = Node(1,-1,-1,-1);
node02 = Node(2,1,-1,-1);
node03 = Node(3,1,1,-1);
node04 = Node(4,-1,1,-1);
node05 = Node(5,-1,-1,1);
node06 = Node(6,1,-1,1);
node07 = Node(7,1,1,1);
node08 = Node(8,-1,1,1);
node09 = Node(9,-1,-1,2);
node10 = Node(10,1,-1,2);
node11 = Node(11,1,1,2);
node12 = Node(12,-1,1,2);



nodeArray = [node01 node02 node03 node04 node05 node06 node07 node08 node09 node10 node11 node12];

nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});

ele01 = HexahedronElement3d8n(1,[node01 node02 node03 node04 node05 node06 node07 node08]);
ele02 = HexahedronElement3d8n(2,[node05 node06 node07 node08 node09 node10 node11 node12]);

elementArray = [ele01 ele02];

elementArray.setPropertyValue('YOUNGS_MODULUS',32);
elementArray.setPropertyValue('POISSON_RATIO',1/3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',3);
elementArray.setPropertyValue('DENSITY',7860);


node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');
node01.fixDof('DISPLACEMENT_Z');
node04.fixDof('DISPLACEMENT_X');
node04.fixDof('DISPLACEMENT_Y');
node04.fixDof('DISPLACEMENT_Z');
node05.fixDof('DISPLACEMENT_X');
node05.fixDof('DISPLACEMENT_Y');
node05.fixDof('DISPLACEMENT_Z');
node08.fixDof('DISPLACEMENT_X');
node08.fixDof('DISPLACEMENT_Y');
node08.fixDof('DISPLACEMENT_Z');

addPointLoad(node03,1,[0 -1]);
addPointLoad(node09,1,[1 0]);

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
v.setScaling(1);
v.plotUndeformed
v.plotDeformed

% a=[16,6,6,-8,2,2,-6,-6,1,4,-2,3,4,3,-2,-6,1,-6,-4,-3,-3,0,-1,-1;
%     6,16,6,-2,4,3,-6,-6,1,2,-8,2,3,4,-2,-1,0,-1,-3,-4,-3,1,-6,-6;
%     6,6,16,-2,3,4,-1,-1,0,3,-2,4,2,2,-8,-6,1,-6,-3,-3,-4,1,-6,-6;
%     -8,-2,-2,16,-6,-6,4,2,-3,-6,6,-1,-6,-1,6,4,-3,2,0,1,1,-4,3,3;
%     2,4,3,-6,16,6,-2,-8,2,6,-6,1,1,0,-1,-3,4,-2,-1,-6,-6,3,-4,-3;
%     2,3,4,-6,6,16,-3,-2,4,1,-1,0,6,1,-6,-2,2,-8,-1,-6,-6,3,-3,-4;
%     -6,-6,-1,4,-2,-3,16,6,-6,-8,2,-2,-4,-3,3,0,-1,1,4,3,2,-6,1,6;
%     -6,-6,-1,2,-8,-2,6,16,-6,-2,4,-3,-3,-4,3,1,-6,6,3,4,2,-1,0,1;
%     1,1,0,-3,2,4,-6,-6,16,2,-3,4,3,3,-4,-1,6,-6,-2,-2,-8,6,-1,-6;
%     4,2,3,-6,6,1,-8,-2,2,16,-6,6,0,1,-1,-4,3,-3,-6,-1,-6,4,-3,-2;
%     -2,-8,-2,6,-6,-1,2,4,-3,-6,16,-6,-1,-6,6,3,-4,3,1,0,1,-3,4,2;
%     3,2,4,-1,1,0,-2,-3,4,6,-6,16,1,6,-6,-3,3,-4,-6,-1,-6,2,-2,-8;
%     4,3,2,-6,1,6,-4,-3,3,0,-1,1,16,6,-6,-8,2,-2,-6,-6,-1,4,-2,-3;
%     3,4,2,-1,0,1,-3,-4,3,1,-6,6,6,16,-6,-2,4,-3,-6,-6,-1,2,-8,-2;
%     -2,-2,-8,6,-1,-6,3,3,-4,-1,6,-6,-6,-6,16,2,-3,4,1,1,0,-3,2,4;
%     -6,-1,-6,4,-3,-2,0,1,-1,-4,3,-3,-8,-2,2,16,-6,6,4,2,3,-6,6,1;
%     1,0,1,-3,4,2,-1,-6,6,3,-4,3,2,4,-3,-6,16,-6,-2,-8,-2,6,-6,-1;
%     -6,-1,-6,2,-2,-8,1,6,-6,-3,3,-4,-2,-3,4,6,-6,16,3,2,4,-1,1,0;
%     -4,-3,-3,0,-1,-1,4,3,-2,-6,1,-6,-6,-6,1,4,-2,3,16,6,6,-8,2,2;
%     -3,-4,-3,1,-6,-6,3,4,-2,-1,0,-1,-6,-6,1,2,-8,2,6,16,6,-2,4,3;
%     -3,-3,-4,1,-6,-6,2,2,-8,-6,1,-6,-1,-1,0,3,-2,4,6,6,16,-2,3,4;
%     0,1,1,-4,3,3,-6,-1,6,4,-3,2,4,2,-3,-6,6,-1,-8,-2,-2,16,-6,-6;
%     -1,-6,-6,3,-4,-3,1,0,-1,-3,4,-2,-2,-8,2,6,-6,1,2,4,3,-6,16,6;
%     -1,-6,-6,3,-3,-4,6,1,-6,-2,2,-8,-3,-2,4,1,-1,0,2,3,4,-6,6,16];
% 
% Fehler=stiffnessMatrix-a;
%     
