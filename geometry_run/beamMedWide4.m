clear;
clc;

%import gmsh model from directory
io = ModelIO('C:\users\Frederik Schaal\Desktop\Bachelorarbeit\Matlab Code\geometry_msh_files\beamMedWide.msh');
model = io.readModel;

%boundary conditions
nodeArray = model.getAllNodes;
elementArray = model.getAllElements;
nodeArray(1).fixDof('DISPLACEMENT_X');
nodeArray(1).fixDof('DISPLACEMENT_Y');
nodeArray(14).fixDof('DISPLACEMENT_Y');
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);

%Point Loads, select points with loads by Id
%vertical loads
pointsVer = [7 15 29 43 57 71 85];
%horizontal loads
pointsHor = [11 25 39 53 67 81];

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
    if elements(ii).getId == 103 || elements(ii).getId == 104 || ...
            elements(ii).getId == 105 || elements(ii).getId == 106 || ...
            elements(ii).getId == 107 || elements(ii).getId == 108
        eleIntf = [eleIntf elements(ii)];
    end
end

%call to substructure
[substructure01, substructure02] = model.divide(eleIntf,idt);

%find other bar elements with certain Ids
elements = substructure02.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 121 || elements(ii).getId == 122 ...
            || elements(ii).getId == 123 || elements(ii).getId == 124 ...
            || elements(ii).getId == 125 || elements(ii).getId == 126
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure03, substructure04] = substructure02.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure04.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 139 || elements(ii).getId == 140 ...
            || elements(ii).getId == 141 || elements(ii).getId == 142 ...
            || elements(ii).getId == 143 || elements(ii).getId == 144
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure05, substructure06] = substructure04.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure01.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 10 || elements(ii).getId == 17 ...
            || elements(ii).getId == 24 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure07, substructure08] = substructure01.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure03.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 31 || elements(ii).getId == 38 ...
            || elements(ii).getId == 45 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure09, substructure10] = substructure03.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure05.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 52 || elements(ii).getId == 59 ...
            || elements(ii).getId == 66 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure11, substructure12] = substructure05.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure06.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 73 || elements(ii).getId == 80 ...
            || elements(ii).getId == 87 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure13, substructure14] = substructure06.divide(eleIntf, idt);

%gather all final substructures in a matrix (substructures)
substructures = [substructure07 substructure08 substructure09 ...
                    substructure10 substructure11 substructure12 ...
                    substructure13 substructure14];

%solve FETI
u = SimpleSolvingStrategy.solve(substructures);

%Order displacements
[displacements, displacementsIntf] = orderDisplacements(substructures);
%Output of displacements
displacements;
