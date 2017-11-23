clear;
clc;

%import gmsh model from directory
io = ModelIO('C:\users\Frederik Schaal\Desktop\Bachelorarbeit\Matlab Code\geometry_msh_files\beamMed5.msh');
model = io.readModel;

%boundary conditions
nodeArray = model.getAllNodes;
elementArray = model.getAllElements;
nodeArray(1).fixDof('DISPLACEMENT_X');
nodeArray(1).fixDof('DISPLACEMENT_Y');
nodeArray(2).fixDof('DISPLACEMENT_X');
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);

%Point Loads, select points with loads by Id
%vertical loads
pointsVer = [5 15 25 35 45 55 65];
% pointsVer = [1 11 21 31 41 51 61];
%horizontal loads
pointsHor = [8 18 28 38 48 58];

for ii = 1:length(pointsVer)
    addPointLoad(nodeArray(pointsVer(ii)),1,[0 -1 0]);
end
for ii = 1:length(pointsHor)
    addPointLoad(nodeArray(pointsHor(ii)),1,[1 0 0]);
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
    if elements(ii).getId == 77 || elements(ii).getId == 78 || ...
            elements(ii).getId == 79 || elements(ii).getId == 80
        eleIntf = [eleIntf elements(ii)];
    end
end

%call to substructure
[substructure01, substructure02] = model.divide(eleIntf,idt);

%find other bar elements with certain Ids
elements = substructure02.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 93 || elements(ii).getId == 94 || elements(ii).getId == 95 ...
            || elements(ii).getId == 96
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure03, substructure04] = substructure02.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure01.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 7 || elements(ii).getId == 12 ...
            || elements(ii).getId == 17 || elements(ii).getId == 22 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure05, substructure06] = substructure01.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure03.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 27 || elements(ii).getId == 32 ...
            || elements(ii).getId == 37 || elements(ii).getId == 42
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure07, substructure08] = substructure03.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure04.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 47 || elements(ii).getId == 52 ...
            || elements(ii).getId == 57 || elements(ii).getId == 62 
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
