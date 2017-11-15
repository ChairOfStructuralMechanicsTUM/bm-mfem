clear;
clc;
%model & substructures
io = ModelIO('board3D.msh');
model = io.readModel;
io = ModelIO('board4_3D_1.msh');
substructure01 = io.readModel;
io = ModelIO('board4_3D_2.msh');
substructure02 = io.readModel;
io = ModelIO('board4_3D_3.msh');
substructure03 = io.readModel;
io = ModelIO('board4_3D_4.msh');
substructure04 = io.readModel;
substructures = [substructure01 substructure02 substructure03 ...
                    substructure04];
                
%boundary conditions
nodeArray = model.getAllNodes;
elementArray = model.getAllElements;

nodeArrays = cell(1,length(substructures));
elementArrays = cell(1,length(substructures));
for numSubs = 1:length(substructures)
    nodeArrays{1,numSubs} = substructures(1,numSubs).getAllNodes;
    elementArrays{1,numSubs} = substructures(1,numSubs).getAllElements;
end

nodeArray(1).fixDof('DISPLACEMENT_X');
nodeArray(1).fixDof('DISPLACEMENT_Y');
nodeArray(1).fixDof('DISPLACEMENT_Z');
nodeArray(2).fixDof('DISPLACEMENT_Y');
nodeArray(2).fixDof('DISPLACEMENT_Z');
nodeArray(191).fixDof('DISPLACEMENT_Y');

nodeArrays{1,2}(1).fixDof('DISPLACEMENT_X');
nodeArrays{1,2}(1).fixDof('DISPLACEMENT_Y');
nodeArrays{1,2}(1).fixDof('DISPLACEMENT_Z');
nodeArrays{1,2}(2).fixDof('DISPLACEMENT_Y');
nodeArrays{1,2}(2).fixDof('DISPLACEMENT_Z');
nodeArrays{1,1}(51).fixDof('DISPLACEMENT_Y');

%Point Loads 
points = [11 12 13 14 15 16 17 18 19 106 107 108 109 110 111 112 113 114 ...
          201 202 203 204 205 206 207 208 209];
for ii = 1:length(points)
    addPointLoad(nodeArray(points(ii)),0.5,[1 -1 0]);
end
% 
% points = [3 4 5 6 13 14 15 16 23 24 25 26];
% for ii = 1:length(points)
%     if ii == 9 || ii == 10 || ii == 11 || ii == 12
%         addPointLoad(nodeArrays{1,1}(points(ii)),0.25,[1 -1 0]);
%     else
%         addPointLoad(nodeArrays{1,1}(points(ii)),0.5,[1 -1 0]);
%     end
% end
% 
% points = [3 4 5 6 13 14 15 16 23 24 25 26];
% for ii = 1:length(points)
%     if ii == 1 || ii == 2 || ii == 3 || ii == 4
%         addPointLoad(nodeArrays{1,2}(points(ii)),0.25,[1 -1 0]);
%     else
%         addPointLoad(nodeArrays{1,2}(points(ii)),0.5,[1 -1 0]);
%     end
% end

points = [2 3 4 5 6 7 8 9 10 ...
          52 53 54 55 56 57 58 59 60];
for ii = 1:length(points)
    if ii == 1 || ii == 2 || ii == 3 || ii == 4 || ii == 5 ||...
            ii == 6 || ii == 7 || ii == 8 || ii == 9
        addPointLoad(nodeArrays{1,3}(points(ii)),0.25,[1 -1 0]);
    else
        addPointLoad(nodeArrays{1,3}(points(ii)),0.5,[1 -1 0]);
    end
end

points = [2 3 4 5 6 7 8 9 10 ...
          52 53 54 55 56 57 58 59 60];
for ii = 1:length(points)
    if ii == 10 || ii == 11 || ii == 12 || ii == 13 || ii == 14 ||...
            ii == 15 || ii == 16 || ii == 17 || ii == 18
        addPointLoad(nodeArrays{1,4}(points(ii)),0.25,[1 -1 0]);
    else
        addPointLoad(nodeArrays{1,4}(points(ii)),0.5,[1 -1 0]);
    end
end

%set IdTracker & change nodeIds, so that they are non repeating
idt = IdTracker(max(nodeArrays{1,1}.getId), ...
            max(elementArrays{1,1}.getId));
for numSubs = 2:length(substructures)
        id = idt.getNodeId;
        for nodes = 1:length(nodeArrays{1,numSubs})
            nodeArrays{1,numSubs}(nodes).setId(id+nodes);
        end
        idt.setNodeId(id+nodes);
end

%solve FEM
u = SimpleSolvingStrategy.solve(model)
                
%solve individual models using FETI
u = SimpleSolvingStrategy.solve(substructures);

% for numSubs = 1:length(substructures)
%     u{1,numSubs}
% end

%Order displacements
[displacements, displacementsIntf] = orderDisplacements(substructures);
displacements
                
%Plots
% v = Visualization(model);
% substructure01 = Visualization(substructure01);
% substructure02 = Visualization(substructure02);
% substructure03 = Visualization(substructure03);
% substructure04 = Visualization(substructure04);
% 
% figure
% plotUndeformed(v);
% line([0 10],[0 0],[0 0], 'Color', 'red', 'LineWidth',5);
% line([0 0],[0 10],[0 0], 'Color', 'green', 'LineWidth',5);
% line([0 0],[0 0],[0 10], 'Color', 'yellow', 'LineWidth',5);
% 
% figure
% plotUndeformed(substructure01);
% line([0 10],[0 0],[0 0], 'Color', 'red', 'LineWidth',5);
% line([0 0],[0 10],[0 0], 'Color', 'green', 'LineWidth',5);
% line([0 0],[0 0],[0 10], 'Color', 'yellow', 'LineWidth',5);
% 
% figure
% plotUndeformed(substructure02);
% line([0 10],[0 0],[0 0], 'Color', 'red', 'LineWidth',5);
% line([0 0],[0 10],[0 0], 'Color', 'green', 'LineWidth',5);
% line([0 0],[0 0],[0 10], 'Color', 'yellow', 'LineWidth',5);
% 
% figure
% plotUndeformed(substructure03);
% line([90 100],[0 0],[0 0], 'Color', 'red', 'LineWidth',5);
% line([90 90],[0 10],[0 0], 'Color', 'green', 'LineWidth',5);
% line([90 90],[0 0],[0 10], 'Color', 'yellow', 'LineWidth',5);
% % % plotDeformed(substructure03);
% figure
% plotUndeformed(substructure04);
% line([90 100],[0 0],[0 0], 'Color', 'red', 'LineWidth',5);
% line([90 90],[0 10],[0 0], 'Color', 'green', 'LineWidth',5);
% line([90 90],[0 0],[0 10], 'Color', 'yellow', 'LineWidth',5);
% % plotDeformed(substructure04);