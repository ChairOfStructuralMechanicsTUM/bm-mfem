clear;

node01 = Node(1,0,0,0);
node02 = Node(2,2,3,6);
node03 = Node(3,1,1,0);
node04 = Node(4,0,1,0);
nodeArray = [node01 node02 node03 node04];

steel = Material('Steel');
steel.addParameter('YOUNGS_MODULUS', 210000);

mat = Material('test');
mat.addParameter('YOUNGS_MODULUS', 343);

area = 10;

ele01 = BarElement3d2n(1,[node01 node02], mat, area);
ele02 = BarElement3d2n(2,[node02 node03], mat, area);
ele03 = BarElement3d2n(3,[node03 node04], mat, area);
ele04 = BarElement3d2n(4,[node04 node01], mat, area);
elementArray = [ele01 ele02 ele03 ele04];

node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');
node02.fixDof('DISPLACEMENT_Y');

model = FemModel(nodeArray, elementArray);

