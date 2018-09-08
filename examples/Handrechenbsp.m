clear;

node01=Node(1,0,0);
node02=Node(2,0,-5);
node03=Node(3,0,-10);
%
node04=Node(4,5,0);
node05=Node(5,5,-5);
node06=Node(6,5,-10);
%
node07=Node(7,10,0);
node08=Node(8,10,-5);
node09=Node(9,10,-10);

nodeArray = [node01 node02 node03 node04 node05 node06 node07 node08 node09];

nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});


ele01 = BarElement2d2n(1,[node01 node02]);
ele02 = BarElement2d2n(2,[node02 node03]);
ele03 = BarElement2d2n(3,[node01 node04]);
ele04 = BarElement2d2n(4,[node02 node04]);
ele05 = BarElement2d2n(5,[node02 node05]);
ele06 = BarElement2d2n(6,[node03 node05]);
ele07 = BarElement2d2n(7,[node03 node06]);
ele08 = BarElement2d2n(8,[node04 node05]);
ele09 = BarElement2d2n(9,[node05 node06]);
ele10 = BarElement2d2n(10,[node04 node07]);
ele11 = BarElement2d2n(11,[node05 node07]);
ele12 = BarElement2d2n(12,[node05 node08]);
ele13 = BarElement2d2n(13,[node06 node08]);
ele14 = BarElement2d2n(14,[node06 node09]);
ele15 = BarElement2d2n(15,[node07 node08]);
ele16 = BarElement2d2n(16,[node08 node09]);



elementArray = [ele01 ele02 ele03 ele04 ele05 ele06 ele07 ele08 ele09 ...
    ele10 ele11 ele12 ele13 ele14 ele15 ele16];


elementArray.setPropertyValue('CROSS_SECTION',1);
elementArray.setPropertyValue('YOUNGS_MODULUS',1000);

elementIds = elementArray.getId;

node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');

node03.fixDof('DISPLACEMENT_X');


addPointLoad(node07,100,[0 1]);

model = FemModel(nodeArray, elementArray);
assembling = SimpleAssembler2(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
[forceVector,reducedForceVector] = assembling.applyExternalForces(model);



solver = SimpleSolvingStrategy(model);
x = solver.solve();

step = 1;

% stressVector = computeElementStress(elementArray, step)';

VerschiebungDofs = model.getDofArray.getValue(step);

nodalForces = solver.getNodalForces(step);


v = Visualization(model);
%v.setScaling(10000000);

v.plotUndeformed
v.plotDeformed
