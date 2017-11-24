clear;
clc;

%create Nodes
node01 = Node(1,0,0,0);
node02 = Node(2,0,10,0);
node03 = Node(3,0,20,0);
node04 = Node(4,0,30,0);
node05 = Node(5,0,40,0);
node06 = Node(6,10,0,0);
node07 = Node(7,10,10,0);
node08 = Node(8,10,20,0);
node09 = Node(9,10,30,0);
node10 = Node(10,10,40,0);
node11 = Node(11,20,0,0);
node12 = Node(12,20,10,0);
node13 = Node(13,20,20,0);
node14 = Node(14,20,30,0);
node15 = Node(15,20,40,0);
node16 = Node(16,30,0,0);
node17 = Node(17,30,10,0);
node18 = Node(18,30,20,0);
node19 = Node(19,30,30,0);
node20 = Node(20,30,40,0);
node21 = Node(21,40,0,0);
node22 = Node(22,40,10,0);
node23 = Node(23,40,20,0);
node24 = Node(24,40,30,0);
node25 = Node(25,40,40,0);


nodeArray = [node01 node02 node03 node04 node05 node06 node07 node08... 
             node09 node10 node11 node12 node13 node14 node15 node16...
             node17 node18 node19 node20 node21 node22 node23 node24...
             node25];

mat = Material('test');
mat.addParameter('YOUNGS_MODULUS', 1000);

%create Elements
ele01 = BarElement3d2n(1,[node01 node02], mat, 2);
ele02 = BarElement3d2n(2,[node01 node06], mat, 2);
ele03 = BarElement3d2n(3,[node01 node07], mat, 2);
ele04 = BarElement3d2n(4,[node02 node07], mat, 2);
ele05 = BarElement3d2n(5,[node02 node08], mat, 2);
ele06 = BarElement3d2n(6,[node02 node03], mat, 2);
ele07 = BarElement3d2n(7,[node03 node08], mat, 2);
ele08 = BarElement3d2n(8,[node03 node09], mat, 2);
ele09 = BarElement3d2n(9,[node03 node04], mat, 2);
ele10 = BarElement3d2n(10,[node04 node09], mat, 2);
ele11 = BarElement3d2n(11,[node04 node10], mat, 2);
ele12 = BarElement3d2n(12,[node04 node05], mat, 2);
ele14 = BarElement3d2n(14,[node05 node10], mat, 2);
ele15 = BarElement3d2n(15,[node06 node11], mat, 2);
ele16 = BarElement3d2n(16,[node06 node12], mat, 2);
ele17 = BarElement3d2n(17,[node06 node07], mat, 2);
ele18 = BarElement3d2n(18,[node07 node12], mat, 2);
ele19 = BarElement3d2n(19,[node07 node13], mat, 2);
ele20 = BarElement3d2n(20,[node07 node08], mat, 2);
ele21 = BarElement3d2n(21,[node08 node13], mat, 2);
ele22 = BarElement3d2n(22,[node08 node14], mat, 2);
ele23 = BarElement3d2n(23,[node08 node09], mat, 2);
ele24 = BarElement3d2n(24,[node09 node14], mat, 2);
ele25 = BarElement3d2n(25,[node09 node15], mat, 2);
ele26 = BarElement3d2n(26,[node09 node10], mat, 2);
ele27 = BarElement3d2n(27,[node10 node15], mat, 2);
ele28 = BarElement3d2n(28,[node11 node16], mat, 2);
ele29 = BarElement3d2n(29,[node11 node17], mat, 2);
ele30 = BarElement3d2n(30,[node11 node12], mat, 2);
ele31 = BarElement3d2n(31,[node12 node17], mat, 2);
ele32 = BarElement3d2n(32,[node12 node18], mat, 2);
ele33 = BarElement3d2n(33,[node12 node13], mat, 2);
ele34 = BarElement3d2n(34,[node13 node18], mat, 2);
ele35 = BarElement3d2n(35,[node13 node19], mat, 2);
ele36 = BarElement3d2n(36,[node13 node14], mat, 2);
ele37 = BarElement3d2n(37,[node14 node19], mat, 2);
ele38 = BarElement3d2n(38,[node14 node20], mat, 2);
ele39 = BarElement3d2n(39,[node14 node15], mat, 2);
ele40 = BarElement3d2n(40,[node15 node20], mat, 2);
ele41 = BarElement3d2n(41,[node16 node21], mat, 2);
ele42 = BarElement3d2n(42,[node16 node22], mat, 2);
ele43 = BarElement3d2n(43,[node16 node17], mat, 2);
ele44 = BarElement3d2n(44,[node17 node22], mat, 2);
ele45 = BarElement3d2n(45,[node17 node23], mat, 2);
ele46 = BarElement3d2n(46,[node17 node18], mat, 2);
ele47 = BarElement3d2n(47,[node18 node23], mat, 2);
ele48 = BarElement3d2n(48,[node18 node24], mat, 2);
ele49 = BarElement3d2n(49,[node18 node19], mat, 2);
ele50 = BarElement3d2n(50,[node19 node24], mat, 2);
ele51 = BarElement3d2n(51,[node19 node25], mat, 2);
ele52 = BarElement3d2n(52,[node19 node20], mat, 2);
ele53 = BarElement3d2n(53,[node20 node25], mat, 2);
ele54 = BarElement3d2n(54,[node21 node22], mat, 2);
ele55 = BarElement3d2n(55,[node22 node23], mat, 2);
ele56 = BarElement3d2n(56,[node23 node24], mat, 2);
ele57 = BarElement3d2n(57,[node24 node25], mat, 2);


elementArray = [ele01 ele02 ele03 ele04 ele05 ele06 ele07 ele08 ele09... 
                ele10 ele11 ele12 ele14 ele15 ele16 ele17 ele18...
                ele19 ele20 ele21 ele22 ele23 ele24 ele25 ele26 ele27...
                ele28 ele29 ele30 ele31 ele32 ele33 ele34 ele35 ele36...
                ele37 ele38 ele39 ele40 ele41 ele42 ele43 ele44 ele45...
                ele46 ele47 ele48 ele49 ele50 ele51 ele52 ele53 ele54...
                ele55 ele56 ele57];

%boundary conditions
node01.fixDof('DISPLACEMENT_X');
node01.fixDof('DISPLACEMENT_Y');
node02.fixDof('DISPLACEMENT_X');
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);

%Point Loads
addPointLoad([node03, node07, node12, node13],10,[1 -1 0]);
addPointLoad([node04],10,[1 -1 0]);
addPointLoad([node20, node22, node23],10,[1 -1 0]);

%set IdTracker with current max NodeId and ElementId
idt = IdTracker(max(nodeArray.getId), max(elementArray.getId));

%create FemModel
modelNine = FemModel(nodeArray, elementArray);

%solve FEM
u = SimpleSolvingStrategy.solve(modelNine);

%Order displacements
[displacements1, displacementsIntf1] = orderDisplacements(modelNine);
%Output of displacements
displacements1;

%Substructure: divide structure into multiple substructures by iteratively
%diving the given structure into two parts
%Function divide performs substructuring. It needs elements along which one
%wants to substructure (eleIntf) and the current IdTracker (idt).

%first elements along which one wants to substructure are collected by Id
eleIntf = [ele30, ele33, ele36, ele39]; 

%call to substructure
[substructure01, substructure02] = modelNine.divide(eleIntf,idt);


%find other bar elements with certain Ids
elements = substructure01.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 7 || elements(ii).getId == 21
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure03, substructure04] = substructure01.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure02.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 34 || elements(ii).getId == 47
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure05, substructure06] = substructure02.divide(eleIntf, idt);

%gather all final substructures in a matrix (substructures)
substructures = [substructure03 substructure04 substructure05 substructure06];

%solve FETI
u = SimpleSolvingStrategy.solve(substructures);

%Order displacements
[displacements, displacementsIntf] = orderDisplacements(substructures);
%Output of displacements
displacements;
