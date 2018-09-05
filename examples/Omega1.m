clear;

node01=Node(1,0,0);
node02=Node(2,0,-5);
node03=Node(3,0,-10);
%
node04=Node(4,5,0);
node05=Node(5,5,-5);
node06=Node(6,5,-10);


nodeArray = [node01 node02 node03 node04 node05 node06];

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




elementArray = [ele01 ele02 ele03 ele04 ele05 ele06 ele07 ele08 ele09];


elementArray.setPropertyValue('CROSS_SECTION',1);
elementArray(1:7).setPropertyValue('YOUNGS_MODULUS',1000);
elementArray(8).setPropertyValue('YOUNGS_MODULUS',500);
elementArray(9).setPropertyValue('YOUNGS_MODULUS',500);

elementIds = elementArray.getId;

node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');

node03.fixDof('DISPLACEMENT_X');



model = FemModel(nodeArray, elementArray);
assembling = SimpleAssembler2(model);
[stiffnessMatrix,reducedStiffnessMatrix] = assembling.assembleGlobalStiffnessMatrix(model);
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
