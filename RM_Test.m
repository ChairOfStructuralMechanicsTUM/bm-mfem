clear;
%example taken from Introduction to FEM (Felippa) 23-5
a=5;
b=2*a;


%Initialize Nodes
node01 = Node(1,0,0,0);
node02 = Node(2,a,0,0);
node03 = Node(3,a,b,0);
node04 = Node(4,0,b,0);
% nodeArray = [node01 node02 node03 node04 ];
nodeArray = [node01 node03 node02 node04 ];
nodeArray.addDof({'DISPLACEMENT_Z', 'ROTATION_X', 'ROTATION_Y'});

ele = ReissnerMindlinElement3d4n(1,[node01 node02 node03 node04]);

ele.setPropertyValue('THICKNESS', 1);
ele.setPropertyValue('YOUNGS_MODULUS', 8000);
ele.setPropertyValue('SHEAR_MODULUS', 4000);
ele.setPropertyValue('POISSON_RATIO', 0.3);


