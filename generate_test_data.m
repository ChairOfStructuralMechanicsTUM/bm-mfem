node01 = Node(1,0,0);
node02 = Node(2,1,0);
node03 = Node(3,1,1);
node04 = Node(4,0,1);
nodeArray = [node01 node02 node03 node04];

steel = Material('Steel');
steel.addParameter('YOUNGS_MODULUS', 210000);

ele01 = BarElement(1,[node01 node02], steel, 1);
ele02 = BarElement(2,[node02 node03], steel, 1);
ele03 = BarElement(3,[node03 node04], steel, 1);
ele04 = BarElement(4,[node04 node01], steel, 1);
elementArray = [ele01 ele02 ele03 ele04];

model = FemModel;
model.loadFemModel(nodeArray, elementArray);

