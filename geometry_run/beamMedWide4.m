clear;
clc;
io = ModelIO('geometry_msh_files/beamMedWide.msh');
model = io.readModel;

%boundary conditions
nodeArray = model.getAllNodes;
elementArray = model.getAllElements;
nodeArray(1).fixDof('DISPLACEMENT_X');
nodeArray(1).fixDof('DISPLACEMENT_Y');
nodeArray(14).fixDof('DISPLACEMENT_Y');
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);

%Point Loads
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

%set IdTracker
idt = IdTracker(max(nodeArray.getId), max(elementArray.getId));

%solve FEM
u = SimpleSolvingStrategy.solve(model);
[displacements1, displacementsIntf1] = orderDisplacements(model);
displacements1


%Substructure 
% eleIntf = elements at interface 
elements = model.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 103 || elements(ii).getId == 104 || ...
            elements(ii).getId == 105 || elements(ii).getId == 106 || ...
            elements(ii).getId == 107 || elements(ii).getId == 108
        eleIntf = [eleIntf elements(ii)];
    end
end
[substructure01, substructure02] = model.divide(eleIntf,idt);

%find new bar element with certain Ids
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

%find new bar element with certain Ids
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

%find new bar element with certain Ids
elements = substructure01.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 10 || elements(ii).getId == 17 ...
            || elements(ii).getId == 24 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure07, substructure08] = substructure01.divide(eleIntf, idt);

%find new bar element with certain Ids
elements = substructure03.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 31 || elements(ii).getId == 38 ...
            || elements(ii).getId == 45 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure09, substructure10] = substructure03.divide(eleIntf, idt);

%find new bar element with certain Ids
elements = substructure05.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 52 || elements(ii).getId == 59 ...
            || elements(ii).getId == 66 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure11, substructure12] = substructure05.divide(eleIntf, idt);

%find new bar element with certain Ids
elements = substructure06.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 73 || elements(ii).getId == 80 ...
            || elements(ii).getId == 87 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure13, substructure14] = substructure06.divide(eleIntf, idt);

substructures = [substructure07 substructure08 substructure09 ...
                    substructure10 substructure11 substructure12 ...
                    substructure13 substructure14];

%solve individual models using FETI
u = SimpleSolvingStrategy.solve(substructures);

%Order displacements
[displacements, displacementsIntf] = orderDisplacements(substructures);
displacements

%Plots
% v = Visualization(model);
% substructure01 = Visualization(substructure07);
% substructure02 = Visualization(substructure08);
% substructure03 = Visualization(substructure09);
% substructure04 = Visualization(substructure10);
% substructure05 = Visualization(substructure11);
% substructure06 = Visualization(substructure12);
% substructure07 = Visualization(substructure13);
% substructure08 = Visualization(substructure14);
% 
% 
% % figure
% % plotUndeformed(v);
% figure
% plotUndeformed(substructure01);
% figure
% plotUndeformed(substructure02);
% figure
% plotUndeformed(substructure03);
% % plotDeformed(substructure03);
% figure
% plotUndeformed(substructure04);
% % plotDeformed(substructure04);
% figure
% plotUndeformed(substructure05);
% % plotDeformed(substructure05);
% figure
% plotUndeformed(substructure06);
% % plotDeformed(substructure06);
% figure
% plotUndeformed(substructure07);
% figure
% plotUndeformed(substructure08);
