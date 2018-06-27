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

leftNodes = [node01 node04 node07];
rightNodes = [node03 node06 node 09];

ele01=[node01 node02 node03];
ele02=[node04 node05 node06];
ele03=[node07 node08 node09];

elementArray=[ele01 ele02 ele03];

model = FemModel(nodeArray, elementArray);

blochInverse1D=BlochInverse1D(model,leftNodes,rightNodes);
assembling