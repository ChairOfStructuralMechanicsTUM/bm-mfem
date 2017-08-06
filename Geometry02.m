
% Manuelle Erzeugung zweier Substructuren:
% Vertikal geschnitten 

clear;


% Sub1
node01 = Node(1,0,5,0);
node02 = Node(2,0,0,0);
node03 = Node(3,5,5,0);
node04 = Node(4,5,0,0);


nodeArray = [node01 node02 node03 node04 ];

mat = Material('test');
mat.addParameter('YOUNGS_MODULUS', 1000);

ele01 = BarElement3d2n(1,[node02 node04], mat, 10);
ele02 = BarElement3d2n(2,[node01 node03], mat, 10);
ele03 = BarElement3d2n(3,[node03 node04], mat, 5);
ele04 = BarElement3d2n(4,[node01 node02], mat, 10);
ele05 = BarElement3d2n(5,[node02 node03], mat, 10);


elementArray = [ele01 ele02 ele03 ele04 ele05];


arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);
node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');  
node02.fixDof('DISPLACEMENT_X');
addPointLoad(node03,50,[0 -1 0]);
Sub1 = FemModel(nodeArray, elementArray);
%Sub1.interfaceNodes=[3 4];
%SimpleSolvingStrategy.solve(Sub1);


% Sub2: [interface-B ele03=ele0 mit node3/4 = node7/8]
node07 = Node(7,10,5,0);
node08 = Node(8,10,0,0);
node05 = Node(5,5,5,0);
node06 = Node(6,5,0,0); 
nodeArray=[node05 node06 node07 node08];

mat = Material('test');
mat.addParameter('YOUNGS_MODULUS', 1000);

ele06 = BarElement3d2n(6,[node06 node08], mat, 10);
ele07 = BarElement3d2n(7,[node05 node07], mat, 10);
ele08 = BarElement3d2n(8,[node05 node06], mat, 5);
ele09 = BarElement3d2n(9,[node06 node07], mat, 10);
ele10 = BarElement3d2n(10,[node07 node08], mat, 10);
elementArray = [ele06 ele07 ele08 ele09 ele10];


arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);
addPointLoad(node05,50,[0 -1 0]);
Sub2 = FemModel(nodeArray, elementArray);
%Sub2interfaceNodes=[3 4];
%[RidgedBodyModes]=computeRidgedBodyModes(Sub2,interfaceNodes);

[u]=CallPCG(Sub1,Sub2);


% Kontroll Ergebnisse aus femModel:
%  1  -0.0500
%  2   0.0500
%  3  -0.2414
%  4  -0.0000
%  5  -0.2414
%  6   0.0500
%  7  -0.2914
%  8  -0.0000
%  9  -0.2914
U1=[-0.0500 0.0500 -0.2414  -0.0000 -0.2414]';
 U2=[ 0.0500 -0.2914 -0.0000 -0.2914 0.0500 -0.2414 -0.0000 -0.2414 ]';