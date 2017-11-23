clear;
clc;

%import gmsh model from directory
io = ModelIO('C:\users\Frederik Schaal\Desktop\Bachelorarbeit\Matlab Code\geometry_msh_files\beamLong.msh');
model = io.readModel;

%boundary conditions
nodeArray = model.getAllNodes;
elementArray = model.getAllElements;
nodeArray(1).fixDof('DISPLACEMENT_X');
nodeArray(1).fixDof('DISPLACEMENT_Y');
nodeArray(6).fixDof('DISPLACEMENT_Y');
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);

%Point Loads, select points with loads by Id
%vertical loads
pointsVer = [5 15 25 35 45 55 65 75 85 95];
%horizontal loads
pointsHor = [8 18 28 38 48 58 68 78 88];
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
    if elements(ii).getId == 36 || elements(ii).getId == 37 || ...
            elements(ii).getId == 38 || elements(ii).getId == 39
        eleIntf = [eleIntf elements(ii)];
    end
end

%call to substructure
[substructure01, substructure02] = model.divide(eleIntf,idt);

%find other bar elements with certain Ids
elements = substructure02.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 135 || elements(ii).getId == 136 || elements(ii).getId == 137 ...
            || elements(ii).getId == 138
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure03, substructure04] = substructure02.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure01.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 68 || elements(ii).getId == 69 ...
            || elements(ii).getId == 70 || elements(ii).getId == 71 ...
            || elements(ii).getId == 72 || elements(ii).getId == 73 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure05, substructure06] = substructure01.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure03.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 74 || elements(ii).getId == 75 ...
            || elements(ii).getId == 76 || elements(ii).getId == 185 ...
            || elements(ii).getId == 186 || elements(ii).getId == 187
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure07, substructure08] = substructure03.divide(eleIntf, idt);

%find other bar elements with certain Ids
elements = substructure04.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 188 || elements(ii).getId == 189 ...
            || elements(ii).getId == 190 || elements(ii).getId == 191 ...
            || elements(ii).getId == 192 || elements(ii).getId == 193 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure09, substructure10] = substructure04.divide(eleIntf, idt);

%gather all final substructures in a matrix (substructures)
substructures = [substructure05 substructure06 substructure07 ...
                    substructure08 substructure09 substructure10];
%solve FETI
u = SimpleSolvingStrategy.solve(substructures);

%Order displacements by coordinates
[displacements, displacementsIntf] = orderDisplacements(substructures);
%Output displacements
displacements;
