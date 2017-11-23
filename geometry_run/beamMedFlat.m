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
    if elements(ii).getId == 109 || elements(ii).getId == 110 || ...
            elements(ii).getId == 111 || elements(ii).getId == 112 || ...
            elements(ii).getId == 113 || elements(ii).getId == 114
        eleIntf = [eleIntf elements(ii)];
    end
end

%call to substructure
[substructure01, substructure02] = model.divide(eleIntf,idt);


%find other bar elements with certain Ids
elements = substructure02.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 133 || elements(ii).getId == 134 ...
            || elements(ii).getId == 135 || elements(ii).getId == 136 ...
            || elements(ii).getId == 137 || elements(ii).getId == 138 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure03, substructure04] = substructure02.divide(eleIntf,idt);

%find other bar elements with certain Ids
elements = substructure01.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 11 || elements(ii).getId == 18 ...
            || elements(ii).getId == 25 || elements(ii).getId == 32
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure05, substructure06] = substructure01.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure06.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 9 || elements(ii).getId == 16 ...
            || elements(ii).getId == 23 || elements(ii).getId == 30 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure07, substructure08] = substructure06.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure03.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 39 || elements(ii).getId == 46 ...
            || elements(ii).getId == 53 || elements(ii).getId == 60
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure09, substructure10] = substructure03.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure10.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 37 || elements(ii).getId == 44 ...
            || elements(ii).getId == 51 || elements(ii).getId == 58
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure11, substructure12] = substructure10.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure04.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 67 || elements(ii).getId == 74 ...
            || elements(ii).getId == 81 || elements(ii).getId == 88
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure13, substructure14] = substructure04.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure14.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 65 || elements(ii).getId == 72 ...
            || elements(ii).getId == 79 || elements(ii).getId == 86
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure15, substructure16] = substructure14.divide(eleIntf, idt);

%gather all final substructures in a matrix (substructures)
substructures = [substructure05 substructure07 substructure08 ...
                    substructure09 substructure11 substructure12 ...
                    substructure13 substructure15 substructure16];

%solve FETI                
u = SimpleSolvingStrategy.solve(substructures);

%Order displacements
[displacements, displacementsIntf] = orderDisplacements(substructures);
%Output of displacements
displacements;
