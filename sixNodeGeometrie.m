clear;

%nodes
node01 = Node(1,0,0,0);
node02 = Node(2,10,0,0);
node03 = Node(3,20,0,0);
node04 = Node(4,20,10,0);
node05 = Node(5,10,10,0);
node06 = Node(6,0,10,0);

nodeArray = [node01 node02 node03 node04 node05 node06];

mat = Material('test');
mat.addParameter('YOUNGS_MODULUS', 1000);

%elements
ele01 = BarElement3d2n(1,[node01 node02], mat, 2);
ele02 = BarElement3d2n(2,[node02 node03], mat, 2);
ele03 = BarElement3d2n(3,[node03 node04], mat, 2);
ele04 = BarElement3d2n(4,[node04 node05], mat, 2);
ele05 = BarElement3d2n(5,[node05 node06], mat, 2);
ele06 = BarElement3d2n(6,[node06 node01], mat, 2);
ele07 = BarElement3d2n(7,[node01 node05], mat, 2);
ele08 = BarElement3d2n(8,[node05 node02], mat, 2);
ele09 = BarElement3d2n(9,[node02 node04], mat, 2);

elementArray = [ele01 ele02 ele03 ele04 ele05 ele06 ele07 ele08 ele09];

%boundary conditions
node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');
node03.fixDof('DISPLACEMENT_Y');
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);

%Load
addPointLoad([node02 node05],10,[0 -1 0]);
addPointLoad(node02,16,[0 -1 0]);

%create FemModel
modelSix = FemModel(nodeArray, elementArray);

%solve
SimpleSolvingStrategy.solve(modelSix);

%Substructure
n = 2; %number of substructures wanted
idsNodes = []; %nxm matrix with all the names of nodes at which we want to substructure (m=number of cut nodes)
idsElements = []; %nxl matrix with all the names of elements which are at the interface(l=number of interface elements)
substructure = zeros(1,n); %matrix containing all the substructures
%ERROR bei mehrfacher Teilung müssen die Restgebiete und nicht wieder
%modelSix ins substructuring gegeben werden.
for subs = 1:1:n
    substructure(1,subs) = Substructure.divide(modelSix, idsNodes(subs,:), idsElements(subs,:), substructure, subs);
end





