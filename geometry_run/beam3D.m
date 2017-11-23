clear;
clc;

%import gmsh model from directory
io = ModelIO('C:\users\Frederik Schaal\Desktop\Bachelorarbeit\Matlab Code\geometry_msh_files\beam3D.msh');
model = io.readModel;

%import substructures which are created in gmsh before 
io = ModelIO('C:\users\Frederik Schaal\Desktop\Bachelorarbeit\Matlab Code\geometry_msh_files\beam4_3D_1.msh');
substructure01 = io.readModel;
io = ModelIO('C:\users\Frederik Schaal\Desktop\Bachelorarbeit\Matlab Code\geometry_msh_files\beam4_3D_2.msh');
substructure02 = io.readModel;
io = ModelIO('C:\users\Frederik Schaal\Desktop\Bachelorarbeit\Matlab Code\geometry_msh_files\beam4_3D_3.msh');
substructure03 = io.readModel;
io = ModelIO('C:\users\Frederik Schaal\Desktop\Bachelorarbeit\Matlab Code\geometry_msh_files\beam4_3D_4.msh');
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

%boundary conditions for FEM
nodeArray(1).fixDof('DISPLACEMENT_X');
nodeArray(1).fixDof('DISPLACEMENT_Y');
nodeArray(1).fixDof('DISPLACEMENT_Z');
nodeArray(2).fixDof('DISPLACEMENT_Z');
nodeArray(2).fixDof('DISPLACEMENT_Y');
nodeArray(20).fixDof('DISPLACEMENT_Z');

%boundary conditions for FETI
nodeArrays{1,1}(1).fixDof('DISPLACEMENT_X');
nodeArrays{1,1}(1).fixDof('DISPLACEMENT_Y');
nodeArrays{1,1}(1).fixDof('DISPLACEMENT_Z');
nodeArrays{1,1}(2).fixDof('DISPLACEMENT_Z');
nodeArrays{1,1}(2).fixDof('DISPLACEMENT_Y');
nodeArrays{1,1}(11).fixDof('DISPLACEMENT_Z');

%Point Loads, select points with loads by Id 
points = [3 4 5 6 12 13 14 18 19 22 23 24 25 31 32 33 37 38 ...
          41 42 43 44 50 51 52 56 57 60 61 62 63 69 70 71 75 76 ...
          79 80 81 82 88 89 90 94 95];
%Point Loads for FEM
for ii = 1:length(points)
    addPointLoad(nodeArray(points(ii)),0.5,[2 0.5 0.5]);
end

%Point Loads for FETI (half point loads at interface nodes)
points = [3 4 5 6 13 14 15 16 23 24 25 26];
for ii = 1:length(points)
    if ii == 9 || ii == 10 || ii == 11 || ii == 12
        addPointLoad(nodeArrays{1,1}(points(ii)),0.25,[2 0.5 0.5]);
    else
        addPointLoad(nodeArrays{1,1}(points(ii)),0.5,[2 0.5 0.5]);
    end
end
points = [3 4 5 6 13 14 15 16 23 24 25 26];
for ii = 1:length(points)
    if ii == 1 || ii == 2 || ii == 3 || ii == 4
        addPointLoad(nodeArrays{1,2}(points(ii)),0.25,[2 0.5 0.5]);
    else
        addPointLoad(nodeArrays{1,2}(points(ii)),0.5,[2 0.5 0.5]);
    end
end
points = [3 4 5 9 10 13 14 15 19 20 23 24 25 29 30];
for ii = 1:length(points)
    if ii == 11 || ii == 12 || ii == 13 || ii == 14 || ii == 15
        addPointLoad(nodeArrays{1,3}(points(ii)),0.25,[2 0.5 0.5]);
    else
        addPointLoad(nodeArrays{1,3}(points(ii)),0.5,[2 0.5 0.5]);
    end
end
points = [3 4 5 9 10 13 14 15 19 20 23 24 25 29 30];
for ii = 1:length(points)
    if ii == 1 || ii == 2 || ii == 3 || ii == 4 || ii == 5
        addPointLoad(nodeArrays{1,4}(points(ii)),0.25,[2 0.5 0.5]);
    else
        addPointLoad(nodeArrays{1,4}(points(ii)),0.5,[2 0.5 0.5]);
    end
end

%set IdTracker with current max NodeId and ElementId & change Ids so that
%they are non repeating
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
u = SimpleSolvingStrategy.solve(model);

%order displacements by coordinates
[displacements1, displacementsIntf1] = orderDisplacements(model);
%Output of displacements
displacements1;
                
%solve FETI
u = SimpleSolvingStrategy.solve(substructures);

%Order displacements
[displacements, displacementsIntf] = orderDisplacements(substructures);
%Output of displacements
displacements;
                