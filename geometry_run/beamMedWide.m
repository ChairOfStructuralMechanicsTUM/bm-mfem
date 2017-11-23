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
[substructure01, substructure02] = model.divide(eleIntf,idt);


%find new bar element with certain Ids
elements = substructure02.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 133 || elements(ii).getId == 134 ...
            || elements(ii).getId == 135 || elements(ii).getId == 136 ...
            || elements(ii).getId == 137 || elements(ii).getId == 138 
        eleIntf = [eleIntf elements(ii)];
    end
end

%call to substructure
[substructure03, substructure04] = substructure02.divide(eleIntf,idt);

%find other bar elements with certain Ids
elements = substructure01.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 10 || elements(ii).getId == 17 ...
            || elements(ii).getId == 24 || elements(ii).getId == 31
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure05, substructure06] = substructure01.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure03.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 38 || elements(ii).getId == 45 ...
            || elements(ii).getId == 52 || elements(ii).getId == 59 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure07, substructure08] = substructure03.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure04.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 66 || elements(ii).getId == 73 ...
            || elements(ii).getId == 80 || elements(ii).getId == 87
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
