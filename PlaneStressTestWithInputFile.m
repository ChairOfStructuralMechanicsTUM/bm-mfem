%% PlaneStressElement Test with Gmsh Input File
%% Import file

io = ModelIO('tests/Beam10m.msh');
model = io.readModel;
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

nodeArray = model.getAllNodes();
elementArray = model.getAllElements();

%% Line BC

idLowerLCorner = model.getModelPart('LowerLCorner').getId();
idUpperLCorner = model.getModelPart('UpperLCorner').getId();

addLineBC(idLowerLCorner,idUpperLCorner,nodeArray);

%% Set Properties

elementArray = model.getAllElements();

elementArray.setPropertyValue('THICKNESS', 1);
elementArray.setPropertyValue('YOUNGS_MODULUS', 100000);
elementArray.setPropertyValue('POISSON_RATIO', 0.3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT', 4);
elementArray.setPropertyValue('DENSITY', 1);

%% Point load
model.getModelPart('RMidPoint').setDofLoad('DISPLACEMENT_Y',-10)
%% Line load
% startNodeId = model.getModelPart('UpperLCorner').getId();
% endNodeId = model.getModelPart('UpperRCorner').getId();
% 
% addConstLineLoad(startNodeId,endNodeId,nodeArray,10,[0 -1 0]);


%% Solving system
solver = SimpleSolvingStrategy(model);
solver.solve();
actualDisplacementY = model.getAllNodes.getDofValue('DISPLACEMENT_Y');

lower_midpoint = model.getModelPart('RMidPoint').getId();
fprintf("Displacement at Node %i : %f\n", lower_midpoint, nodeArray(lower_midpoint).getDofValue('DISPLACEMENT_Y'));

%% Plot 

vis = Visualization(model);
vis.plotUndeformed();
vis.plotDeformed();