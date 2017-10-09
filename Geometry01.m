% Example Geometry

% 1      3       5
% °------°-------°
% |     /|     / |
% |   /  |   /   |
% | /    | /     |
% °------°-------°
% 2      4       6

clear;
node01 = Node(1,0,5,0);
node02 = Node(2,0,0,0);
node03 = Node(3,5,5,0);
node04 = Node(4,5,0,0);
node05 = Node(5,10,5,0);
node06 = Node(6,10,0,0);


nodeArray = [node01 node02 node03 node04 node05 node06 ];

mat = Material('test');
mat.addParameter('YOUNGS_MODULUS', 1000);

ele01 = BarElement3d2n(8,[node02 node04], mat, 10);
ele02 = BarElement3d2n(9,[node04 node06], mat, 10);
ele03 = BarElement3d2n(1,[node01 node03], mat, 10);
ele04 = BarElement3d2n(2,[node03 node05], mat, 10);
ele05 = BarElement3d2n(3,[node01 node02], mat, 10);
ele06 = BarElement3d2n(5,[node03 node04], mat, 10);
ele07 = BarElement3d2n(7,[node05 node06], mat, 10);
ele08 = BarElement3d2n(4,[node02 node03], mat, 10);
ele09 = BarElement3d2n(6,[node04 node05], mat, 10);


elementArray = [ele01 ele02 ele03 ele04 ele05 ele06 ele07 ele08 ele09];

node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');  
node02.fixDof('DISPLACEMENT_X');
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);


addPointLoad(node03,100,[0 -1 0]);

model = FemModel(nodeArray, elementArray);
x=SimpleSolvingStrategy.solve(model);
% SUBSTRRUCTURE

 % dim:{1=X;2=Y;3=Z}, Boundary= InterfaceCoordinate
  Substructure=createSubstructure(model,1,5);  
  Substructure01=Substructure(1);
  Substructure02=Substructure(2);
 
 updateDofs(model.nodeArray,Substructure01,Substructure02)
 u=CallPCG(Substructure01,Substructure02);

