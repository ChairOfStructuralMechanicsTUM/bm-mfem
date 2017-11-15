clear;
clc;
io = ModelIO('beam.msh');
model = io.readModel;

%boundary conditions
nodeArray = model.getAllNodes;
elementArray = model.getAllElements;
nodeArray(1).fixDof('DISPLACEMENT_X');
nodeArray(1).fixDof('DISPLACEMENT_Y');
nodeArray(2).fixDof('DISPLACEMENT_X');
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), nodeArray);

%Point Loads
points = [5 10 15 20 25 30 35 40 45 50];
for ii = 1:length(points)
    addPointLoad(nodeArray(points(ii)),1,[1 -1 0]);
end

%set IdTracker
idt = IdTracker(max(nodeArray.getId), max(elementArray.getId));

%solve FEM
u = SimpleSolvingStrategy.solve(model);

%Order displacements
[displacements1, displacementsIntf1] = orderDisplacements(model);
%Output of displacements
displacements1

%Substructure 
% eleIntf = elements at interface 
elements = model.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 18 || elements(ii).getId == 19 || ...
            elements(ii).getId == 20 || elements(ii).getId == 21
        eleIntf = [eleIntf elements(ii)];
    end
end
[substructure01, substructure02] = model.divide(eleIntf,idt);
% 
% 
% %find new bar element with certain Ids
elements = substructure02.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 36 || elements(ii).getId == 37 || elements(ii).getId == 38 ...
            || elements(ii).getId == 39
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure03, substructure04] = substructure02.divide(eleIntf, idt);
% 
% %find new bar element with certain Ids
elements = substructure01.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 68 || elements(ii).getId == 69 ...
            || elements(ii).getId == 70 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure05, substructure06] = substructure01.divide(eleIntf, idt);
% 
%find new bar element with certain Ids
elements = substructure03.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 71 || elements(ii).getId == 72 ...
            || elements(ii).getId == 73 
        eleIntf = [eleIntf elements(ii)];
    end
end

[substructure07, substructure08] = substructure03.divide(eleIntf, idt);

%find new bar element with certain Ids
elements = substructure04.getAllElements;
eleIntf = [];
for ii = 1:length(elements)
    if elements(ii).getId == 74 || elements(ii).getId == 75 ...
            || elements(ii).getId == 76 
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

