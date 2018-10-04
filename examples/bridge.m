clear;

node01 = Node(1,0,0,0);
node02 = Node(2,10,5,0);
node03 = Node(3,10,0,0);
node04 = Node(4,20,8,0);
node05 = Node(5,20,0,0);
node06 = Node(6,30,9,0);
node07 = Node(7,30,0,0);
node08 = Node(8,40,8,0);
node09 = Node(9,40,0,0);
node10 = Node(10,50,5,0);
node11 = Node(11,50,0,0);
node12 = Node(12,60,0,0);

nodeArray = [node01 node02 node03 node04 node05 node06 node07 node08 node09 ...
    node10 node11 node12];

nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});

ele01 = BarElement3d2n(1,[node01 node03]);
ele02 = BarElement3d2n(2,[node03 node05]);
ele03 = BarElement3d2n(3,[node05 node07]);
ele04 = BarElement3d2n(4,[node07 node09]);
ele05 = BarElement3d2n(5,[node09 node11]);
ele06 = BarElement3d2n(6,[node11 node12]);
ele07 = BarElement3d2n(7,[node01 node02]);
ele08 = BarElement3d2n(8,[node02 node04]);
ele09 = BarElement3d2n(9,[node04 node06]);
ele10 = BarElement3d2n(10,[node06 node08]);
ele11 = BarElement3d2n(11,[node08 node10]);
ele12 = BarElement3d2n(12,[node10 node12]);
ele13 = BarElement3d2n(13,[node02 node03]);
ele14 = BarElement3d2n(14,[node04 node05]);
ele15 = BarElement3d2n(15,[node06 node07]);
ele16 = BarElement3d2n(16,[node08 node09]);
ele17 = BarElement3d2n(17,[node10 node11]);
ele18 = BarElement3d2n(18,[node02 node05]);
ele19 = BarElement3d2n(19,[node04 node07]);
ele20 = BarElement3d2n(20,[node07 node08]);
ele21 = BarElement3d2n(21,[node09 node10]);

elementArray = [ele01 ele02 ele03 ele04 ele05 ele06 ele07 ele08 ele09 ...
    ele10 ele11 ele12 ele13 ele14 ele15 ele16 ele17 ele18 ele19 ...
    ele20 ele21];
tmp = [ele18 ele19 ele20 ele21];
tmp.setPropertyValue('CROSS_SECTION',1);
tmp = [ele01 ele02 ele03 ele04 ele05 ele06];
tmp.setPropertyValue('CROSS_SECTION',2);
tmp = [ele07 ele08 ele09 ele10 ele11 ele12];
tmp.setPropertyValue('CROSS_SECTION',10);
tmp = [ele13 ele14 ele15 ele16 ele17];
tmp.setPropertyValue('CROSS_SECTION',3);

elementArray.setPropertyValue('YOUNGS_MODULUS',1000);


elementIds = elementArray.getId;

node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');
node12.fixDof('DISPLACEMENT_Y');
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);

addPointLoad([node03 node05 node09 node11],10,[0 -1 0]);
addPointLoad(node07,16,[0 -1 0]);

model = FemModel(nodeArray, elementArray);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
forceVector = assembling.applyExternalForces(model);
% reducedForceVector = assembling.reducedForceVector;

solver = SimpleSolvingStrategy(model);
x = solver.solve();

step = 1;

stressVector = computeElementStress(elementArray, step)';  %% besser in Model? sodass für alle elemente


VerschiebungDofs = model.getDofArray.getValue(step);

nodalForces = solver.getNodalForces(step);




n02d = node02.getDofArray;
n02d(1).getValueType;
n02d(1).getValue(step);

n09d = node09.getDofArray;              
n09d(2).getValueType;
n09d(2).getValue(step);

n01d = node01.getDofArray;
n01d(1).getValueType;
n01d(1).getValue(step);


