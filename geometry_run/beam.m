clear;
clc;

%import gmsh model from directory
io = ModelIO('C:\users\Frederik Schaal\Desktop\Bachelorarbeit\Matlab Code\geometry_msh_files\beam.msh');
model = io.readModel;

%boundary conditions
nodeArray = model.getAllNodes;
elementArray = model.getAllElements;
nodeArray(1).fixDof('DISPLACEMENT_X');
nodeArray(1).fixDof('DISPLACEMENT_Y');
nodeArray(2).fixDof('DISPLACEMENT_X');
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);

%Point Loads, select points with loads by Id
points = [5 10 15 20 25 30 35 40 45 50];
for ii = 1:length(points)
    addPointLoad(nodeArray(points(ii)),1,[1 -1 0]);
end

%set IdTracker with current max NodeId and ElementId
idt = IdTracker(max(nodeArray.getId), max(elementArray.getId));

%solve FEM
u = SimpleSolvingStrategy.solve(model);

%Order displacements
[displacements1, displacementsIntf1] = orderDisplacements(model);
%Output of displacements
displacements1;

%Substructure: divide structure into multiple substructures by iteratively
%diving the given structure into two parts
%Function divide performs substructuring. It needs elements along which one
%wants to substructure (eleIntf) and the current IdTracker (idt).

%first elements along which one wants to substructure are collected by Id
elements = model.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 18 || elements(ii).getId == 19 || ...
            elements(ii).getId == 20 || elements(ii).getId == 21
        eleIntf = [eleIntf elements(ii)];
    end
end

%call to substructure
[substructure01, substructure02] = model.divide(eleIntf,idt);

%find other bar elements with certain Ids
elements = substructure02.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 36 || elements(ii).getId == 37 || elements(ii).getId == 38 ...
            || elements(ii).getId == 39
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure03, substructure04] = substructure02.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure01.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 68 || elements(ii).getId == 69 ...
            || elements(ii).getId == 70 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure05, substructure06] = substructure01.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure03.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 71 || elements(ii).getId == 72 ...
            || elements(ii).getId == 73 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure07, substructure08] = substructure03.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure04.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 74 || elements(ii).getId == 75 ...
            || elements(ii).getId == 76 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure09, substructure10] = substructure04.divide(eleIntf, idt);

%gather all final substructures in a matrix (substructures)
substructures = [substructure05 substructure06 substructure07 ...
                    substructure08 substructure09 substructure10];

%solve FETI
u = SimpleSolvingStrategy.solve(substructures);

%Order displacements
[displacements, displacementsIntf] = orderDisplacements(substructures);
%Output of displacements
displacements;
