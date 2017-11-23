clear;
clc;

%import gmsh model from directory
io = ModelIO('C:\users\Frederik Schaal\Desktop\Bachelorarbeit\Matlab Code\geometry_msh_files\beamMed6.msh');
model = io.readModel;

%boundary conditions
nodeArray = model.getAllNodes;
elementArray = model.getAllElements;
nodeArray(1).fixDof('DISPLACEMENT_X');
nodeArray(1).fixDof('DISPLACEMENT_Y');
nodeArray(11).fixDof('DISPLACEMENT_Y');
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);

%Point Loads, select points with loads by Id
%vertical loads
pointsVer = [5 15 25 35 45 55 65 75];
%horizontal loads
pointsHor = [8 18 28 38 48 58 68 78];

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
    if elements(ii).getId == 96 || elements(ii).getId == 97 || ...
            elements(ii).getId == 98 || elements(ii).getId == 99
        eleIntf = [eleIntf elements(ii)];
    end
end

%call to substructure
[substructure01, substructure02] = model.divide(eleIntf,idt);

%find other bar elements with certain Ids
elements = substructure02.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 116 || elements(ii).getId == 117 || elements(ii).getId == 118 ...
            || elements(ii).getId == 119
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure03, substructure04] = substructure02.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure01.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 7 || elements(ii).getId == 12 ...
            || elements(ii).getId == 17 || elements(ii).getId == 22 ...
            || elements(ii).getId == 27
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure05, substructure06] = substructure01.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure03.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 32 || elements(ii).getId == 37 ...
            || elements(ii).getId == 42 || elements(ii).getId == 47 ...
            || elements(ii).getId == 52
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure07, substructure08] = substructure03.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure04.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 57 || elements(ii).getId == 62 ...
            || elements(ii).getId == 67 || elements(ii).getId == 72 ...
            || elements(ii).getId == 77
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
