% 
clear
node01 = Node(1,0,10,0);
node02 = Node(2,0,5,0);
node03 = Node(3,0,0,0);

node04 = Node(4,5,10,0);
node05 = Node(5,5,5,0);
node06 = Node(6,5,0,0);


nodeArray = [node01 node02 node03 node04 node05 node06];

mat = Material('test');
mat.addParameter('YOUNGS_MODULUS', 1000);

%horizontal
ele01 = BarElement3d2n(1,[node01 node04], mat, 10);
ele02 = BarElement3d2n(2,[node02 node05], mat, 10);
ele03 = BarElement3d2n(3,[node03 node06], mat, 10);
%vertikal
ele04 = BarElement3d2n(4,[node01 node02], mat, 10);
ele05 = BarElement3d2n(5,[node04 node05], mat, 10);
ele06 = BarElement3d2n(6,[node02 node03], mat, 10);
ele07 = BarElement3d2n(7,[node05 node06], mat, 10);
%diagonal
ele08 = BarElement3d2n(8,[node02 node04], mat, 10);
ele09 = BarElement3d2n(9,[node03 node05], mat, 10);

elementArray=[ele01 ele02 ele03 ele04 ele05 ele06 ele07 ele08...
              ele09];

node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');
node04.fixDof('DISPLACEMENT_Y');
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);
addPointLoad([node02 node05],100000000,[-1 -1 0]);
%addPointLoad([node11],10,[-1 0 0]);
model = FemModel(nodeArray, elementArray);

x = SimpleSolvingStrategy.solve(model);





[Substructure]=createSubstructure(model,2,5);  
Substructure02=Substructure(1);
Substructure01=Substructure(2);
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), Substructure02.nodeArray);
addPointLoad(Substructure01.nodeArray(end-1:end),50000000,[-1 -1 0]);
addPointLoad(Substructure02.nodeArray(end-1:end),50000000,[-1 -1 0]);
%updateDofs(model.nodeArray,Substructure01,Substructure02)
[u]=CallPCG(Substructure01,Substructure02);
  
 Fabs=norm(x-u);
  Frel=norm(x-u)/norm(x);