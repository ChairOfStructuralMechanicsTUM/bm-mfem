clear;

node01 = Node(1,0,0);
node02 = Node(2,2,0);
node03 = Node(3,4,0);
node04 = Node(4,0,1);
node05 = Node(5,2,1);
node06 = Node(6,4,1);
node07 = Node(7,0,2);
node08 = Node(8,2,2);
node09 = Node(9,4,2);

nodeArray = [node01 node02 node03 node04 node05 ...
    node06 node07 node08 node09];

nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

leftDofIDs=[1,2,3,4,5,6];
rightDofIDs=[13,14,15,16,17,18];

leftNodes = [node01 node04 node07];
rightNodes = [node03 node06 node 09];

ele01=[node01 node02 node03]; %alternativ über Quadril.El.2d4n?
ele02=[node04 node05 node06]; 
ele03=[node07 node08 node09]; 

%es fehlt: Angaben zu Masse und Steifigkeit
elementArray=[ele01 ele02 ele03];

model = FemModel(nodeArray, elementArray);

blochInverse1D=BlochInverse1D(model,leftNodes,rightNodes);
assembling=SimpleAssembler(model);

stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
            
massMatrix = assembling.assembleGlobalMassMatrix(model);

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
