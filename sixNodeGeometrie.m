clear;

%nodes
% node01 = Node(1,0,0,0);
% node02 = Node(2,0,10,0);
% node03 = Node(3,10,0,0);
% node04 = Node(4,10,10,0);
% node05 = Node(5,20,0,0);
% node06 = Node(6,20,10,0);

node01 = Node(1,0,0,0);
node02 = Node(2,0,10,0);
node03 = Node(3,0,20,0);
node04 = Node(4,10,0,0);
node05 = Node(5,10,10,0);
node06 = Node(6,10,20,0);
node07 = Node(7,20,0,0);
node08 = Node(8,20,10,0);
node09 = Node(9,20,20,0);

%%%SUBSTRUCTURE01
% node01 = Node(1,0,0,0);
% node02 = Node(2,0,10,0);
% node03 = Node(3,0,20,0);
% node04 = Node(4,10,0,0);
% node05 = Node(5,10,10,0);
% node06 = Node(6,10,20,0);
% 
% nodeArray = [node01 node02 node03 node04 node05 node06];

%%%SUBSTRUCTURE02
% node04 = Node(1,10,0,0);
% node05 = Node(2,10,10,0);
% node06 = Node(3,10,20,0);
% node07 = Node(4,20,0,0);
% node08 = Node(5,20,10,0);
% node09 = Node(6,20,20,0);
% 
% nodeArray = [node04 node05 node06 node07 node08 node09];



nodeArray = [node01 node02 node03 node04 node05 node06 node07 node08 node09];
%nodeArray = [node01 node02 node03 node04 node05 node06];


mat = Material('test');
mat.addParameter('YOUNGS_MODULUS', 1000);

%elements
% ele01 = BarElement3d2n(1,[node01 node03], mat, 2);
% ele02 = BarElement3d2n(2,[node03 node05], mat, 2);
% ele03 = BarElement3d2n(3,[node05 node06], mat, 2);
% ele04 = BarElement3d2n(4,[node06 node04], mat, 2);
% ele05 = BarElement3d2n(5,[node04 node02], mat, 2);
% ele06 = BarElement3d2n(6,[node02 node01], mat, 2);
% ele07 = BarElement3d2n(7,[node03 node04], mat, 2);
% ele08 = BarElement3d2n(8,[node01 node04], mat, 2);
% ele09 = BarElement3d2n(9,[node03 node06], mat, 2);


ele01 = BarElement3d2n(1,[node01 node04], mat, 2);
ele02 = BarElement3d2n(2,[node04 node07], mat, 2);
ele03 = BarElement3d2n(3,[node07 node08], mat, 2);
ele04 = BarElement3d2n(4,[node08 node09], mat, 2);
ele05 = BarElement3d2n(5,[node09 node06], mat, 2);
ele06 = BarElement3d2n(6,[node06 node03], mat, 2);
ele07 = BarElement3d2n(7,[node03 node02], mat, 2);
ele08 = BarElement3d2n(8,[node02 node01], mat, 2);
ele09 = BarElement3d2n(9,[node01 node05], mat, 2);
ele10 = BarElement3d2n(10,[node04 node08], mat, 2);
ele11 = BarElement3d2n(11,[node05 node08], mat, 2);
ele12 = BarElement3d2n(12,[node02 node05], mat, 2);
ele13 = BarElement3d2n(13,[node02 node06], mat, 2);
ele14 = BarElement3d2n(14,[node05 node09], mat, 2);
ele15 = BarElement3d2n(15,[node04 node05], mat, 2);
ele16 = BarElement3d2n(16,[node05 node06], mat, 2);

%%%SBSTRUCTURE01
% ele01 = BarElement3d2n(1,[node01 node04], mat, 2);
% ele06 = BarElement3d2n(2,[node06 node03], mat, 2);
% ele07 = BarElement3d2n(3,[node03 node02], mat, 2);
% ele08 = BarElement3d2n(4,[node02 node01], mat, 2);
% ele09 = BarElement3d2n(5,[node01 node05], mat, 2);
% ele12 = BarElement3d2n(6,[node02 node05], mat, 2);
% ele13 = BarElement3d2n(7,[node02 node06], mat, 2);
% ele15 = BarElement3d2n(8,[node04 node05], mat, 2);
% ele16 = BarElement3d2n(9,[node05 node06], mat, 2);

%elementArray = [ele01 ele06 ele07 ele08 ele09 ele12 ele13 ele15 ele16];

%%%SUBSTRUCTURE02
% ele02 = BarElement3d2n(1,[node04 node07], mat, 2);
% ele03 = BarElement3d2n(2,[node07 node08], mat, 2);
% ele04 = BarElement3d2n(3,[node08 node09], mat, 2);
% ele05 = BarElement3d2n(4,[node09 node06], mat, 2);
% ele10 = BarElement3d2n(5,[node04 node08], mat, 2);
% ele11 = BarElement3d2n(6,[node05 node08], mat, 2);
% ele14 = BarElement3d2n(7,[node05 node09], mat, 2);
% ele15 = BarElement3d2n(8,[node04 node05], mat, 2);
% ele16 = BarElement3d2n(9,[node05 node06], mat, 2);
% 
% elementArray = [ele02 ele03 ele04 ele05 ele10 ele11 ele14 ele15 ele16];


elementArray = [ele01 ele02 ele03 ele04 ele05 ele06 ele07 ele08 ele09 ele10 ele11 ele12 ele13 ele14 ele15 ele16];
%elementArray = [ele01 ele02 ele03 ele04 ele05 ele06 ele07 ele08 ele09];


%boundary conditions
node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');
node07.fixDof('DISPLACEMENT_Y');
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);

%Load
%addPointLoad(node02,10,[0 1 0]);
%addPointLoad(node03,16,[0 -1 0]);
%addPointLoad([node05 node06],10,[1 -1 0]);

addPointLoad(node02,10,[0 1 0]);
%addPointLoad(node04,16,[0 -1 0]);
addPointLoad(node08,10,[0 1 0]);

%create FemModel
modelNine = FemModel(nodeArray, elementArray);

%Substructure 
%eleIntf = elements at interface
%eleIntf = [ele07];
eleIntf = [ele15, ele16]; 
[substructure01, substructure02] = modelNine.divide(eleIntf);

%solve
%SimpleSolvingStrategy.solve(modelNine);

%solve indicidual models using FETI
SimpleSolvingStrategy.solve(substructure01, substructure02);

%Visualize Substructures
%substructure01 = Visualization(substructure01);
%substructure02 = Visualization(substructure02);
originalSystem = Visualization(modelNine);


% figure
% plotUndeformed(substructure01);
% plotDeformed(substructure01);
% figure
% plotUndeformed(substructure02);
% plotDeformed(substructure02);
% figure
% plotUndeformed(originalSystem);
% plotDeformed(originalSystem);



