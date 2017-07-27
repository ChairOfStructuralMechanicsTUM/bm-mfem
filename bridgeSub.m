clear;

% Bridge with Substructuring 
node01 = Node(1,0,0,0);
node02 = Node(2,10,5,0);
node03 = Node(3,10,0,0);
node04 = Node(4,20,8,0);
node05 = Node(5,20,0,0);
node06 = Node(6,30,9,0);
node07 = Node(7,30,0,0);

node13 = Node(13,30,9,0);
node14 = Node(14,30,0,0);
node08 = Node(8,40,8,0);
node09 = Node(9,40,0,0);
node10 = Node(10,50,5,0);
node11 = Node(11,50,0,0);
node12 = Node(12,60,0,0);

nodeArray01 = [node01 node02 node03 node04 node05 node06 node07 ];

mat = Material('test');
mat.addParameter('YOUNGS_MODULUS', 1000);

ele01 = BarElement3d2n(1,[node01 node03], mat, 2);
ele02 = BarElement3d2n(2,[node03 node05], mat, 2);
ele03 = BarElement3d2n(3,[node05 node07], mat, 2);
ele04 = BarElement3d2n(4,[node07 node09], mat, 2);
ele05 = BarElement3d2n(5,[node09 node11], mat, 2);
ele06 = BarElement3d2n(6,[node11 node12], mat, 2);
ele07 = BarElement3d2n(7,[node01 node02], mat, 10);
ele08 = BarElement3d2n(8,[node02 node04], mat, 10);
ele09 = BarElement3d2n(9,[node04 node06], mat, 10);
ele10 = BarElement3d2n(10,[node06 node08], mat, 10);
ele11 = BarElement3d2n(11,[node08 node10], mat, 10);
ele12 = BarElement3d2n(12,[node10 node12], mat, 10);
ele13 = BarElement3d2n(13,[node02 node03], mat, 3);
ele14 = BarElement3d2n(14,[node04 node05], mat, 3);
ele15 = BarElement3d2n(15,[node06 node07], mat, 1.5);  % Boundary
ele16 = BarElement3d2n(16,[node08 node09], mat, 3);
ele17 = BarElement3d2n(17,[node10 node11], mat, 3);
ele18 = BarElement3d2n(18,[node02 node05], mat, 1);
ele19 = BarElement3d2n(19,[node04 node07], mat, 1);
ele20 = BarElement3d2n(20,[node07 node08], mat, 1);
ele21 = BarElement3d2n(21,[node09 node10], mat, 1);

ele22 = BarElement3d2n(22,[node13 node14], mat, 1.5);  % Boundary
ele23 = BarElement3d2n(23,[node14 node09], mat, 2);
ele24 = BarElement3d2n(24,[node13 node08], mat, 10);
ele25 = BarElement3d2n(25,[node14 node08], mat, 1);
elementArray01 = [ele01 ele02 ele03 ele07 ele08 ele09 ele13...
    ele14 ele15 ele18 ele19];


node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');
node05.fixDof('DISPLACEMENT_Y');

arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray01);

addPointLoad([node03 node05 ],10,[0 -1 0]);
Substructure01 = FemModel(nodeArray01, elementArray01);

% Sub2:


nodeArray02=[node08 node09 node10 node11 node12 node13 node14];

elementArray02=[ ele05 ele06 ele11 ele12 ele16 ele17...
    ele21 ele22 ele23 ele24 ele25];

mat = Material('test');
mat.addParameter('YOUNGS_MODULUS', 1000);

addPointLoad([node09 node11],10,[0 -1 0]);
%addPointLoad(node14,8,[0 -1 0]);
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray02);

Substructure02 = FemModel(nodeArray02, elementArray02);

[u1, u2]=CallPCG(Substructure01,Substructure02,[6 7]);
