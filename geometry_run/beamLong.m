clear;
clc;

io = ModelIO('geometry_msh_files/beamLong.msh');
model = io.readModel;

%boundary conditions
nodeArray = model.getAllNodes;
elementArray = model.getAllElements;
nodeArray(1).fixDof('DISPLACEMENT_X');
nodeArray(1).fixDof('DISPLACEMENT_Y');
nodeArray(6).fixDof('DISPLACEMENT_Y');
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);

%Point Loads
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

%set IdTracker
idt = IdTracker(max(nodeArray.getId), max(elementArray.getId));

%solve FEM
u = SimpleSolvingStrategy.solve(model);

%output ordered displacements
[displacements1, displacementsIntf1] = orderDisplacements(model);
displacements1

%Substructure 
%eleIntf = elements at interface 
elements = model.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 36 || elements(ii).getId == 37 || ...
            elements(ii).getId == 38 || elements(ii).getId == 39
        eleIntf = [eleIntf elements(ii)];
    end
end
[substructure01, substructure02] = model.divide(eleIntf,idt);


%find new bar element with certain Ids
elements = substructure02.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 135 || elements(ii).getId == 136 || elements(ii).getId == 137 ...
            || elements(ii).getId == 138
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure03, substructure04] = substructure02.divide(eleIntf, idt);

%find new bar element with certain Ids
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

%find new bar element with certain Ids
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

%find new bar element with certain Ids
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

substructures = [substructure05 substructure06 substructure07 ...
                    substructure08 substructure09 substructure10];
%solve individual models using FETI
u = SimpleSolvingStrategy.solve(substructures);

%Order displacements by coordinates
[displacements, displacementsIntf] = orderDisplacements(substructures);
%Output displacements
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

