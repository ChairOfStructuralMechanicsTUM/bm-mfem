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
[displacements, displacementsIntf] = orderDisplacements(model);
displacements


%Substructure 
% eleIntf = elements at interface 
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

[substructure03, substructure04] = substructure02.divide(eleIntf,idt);

%find new bar element with certain Ids
elements = substructure01.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 10 || elements(ii).getId == 17 ...
            || elements(ii).getId == 24 || elements(ii).getId == 31
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure05, substructure06] = substructure01.divide(eleIntf, idt);
% 
%find new bar element with certain Ids
elements = substructure03.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 38 || elements(ii).getId == 45 ...
            || elements(ii).getId == 52 || elements(ii).getId == 59 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure07, substructure08] = substructure03.divide(eleIntf, idt);

%find new bar element with certain Ids
elements = substructure04.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 66 || elements(ii).getId == 73 ...
            || elements(ii).getId == 80 || elements(ii).getId == 87
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure09, substructure10] = substructure04.divide(eleIntf, idt);

substructures = [substructure05 substructure06 substructure07 ...
                    substructure08 substructure09 substructure10];
% %solve individual models using FETI
%substructures = [substructure01 substructure02];
u = SimpleSolvingStrategy.solve(substructures);

%Order displacements
[displacements, displacementsIntf] = orderDisplacements(substructures);
displacements

%Plots
% v = Visualization(model);
% substructure01 = Visualization(substructure01);
% substructure02 = Visualization(substructure02);
% substructure03 = Visualization(substructure03);
% substructure04 = Visualization(substructure04);
% substructure05 = Visualization(substructure05);
% substructure06 = Visualization(substructure06);
% substructure07 = Visualization(substructure07);
% substructure08 = Visualization(substructure08);
% substructure09 = Visualization(substructure09);
% substructure10 = Visualization(substructure10);

% figure
% plotUndeformed(v);
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
% figure
% plotUndeformed(substructure09);
% figure
% plotUndeformed(substructure10);