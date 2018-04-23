%% PlaneStressElement Test with Gmsh Input File
%% Import file

io = ModelIO('tests/Trig2ndOrd10x40.msh');
model = io.readModel;
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

nodeArray = model.getAllNodes();
elementArray = model.getAllElements();

%% BC
model.getModelPart('LLCorner').fixDof('DISPLACEMENT_X');
model.getModelPart('LLCorner').fixDof('DISPLACEMENT_Y');
model.getModelPart('LRCorner').fixDof('DISPLACEMENT_Y');

%% Line BC

% idLowerLCorner = model.getModelPart('LLCorner').getId();
% idUpperLCorner = model.getModelPart('ULCorner').getId();
% 
% addLineBC(idLowerLCorner,idUpperLCorner,nodeArray);

%% Set Properties
E = 100000;
elementArray.setPropertyValue('THICKNESS', 1);
elementArray.setPropertyValue('YOUNGS_MODULUS', E);
elementArray.setPropertyValue('POISSON_RATIO', 0.3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT', 3);
elementArray.setPropertyValue('DENSITY', 1);

%% Point load
% model.getModelPart('RMidPoint').setDofLoad('DISPLACEMENT_Y',-10)
%% Line load
q = 1;
startNodeId = model.getModelPart('ULCorner').getId();
endNodeId = model.getModelPart('URCorner').getId();

addConstLineLoad(startNodeId,endNodeId,nodeArray,q,[0 -1 0]);


%% Solving system
solver = SimpleSolvingStrategy(model);
solver.solve();
actualDisplacementY = model.getAllNodes.getDofValue('DISPLACEMENT_Y');

evaluation_point = model.getModelPart('LMPoint').getId();
dy = nodeArray(evaluation_point).getDofValue('DISPLACEMENT_Y');
fprintf("Displacement at Node %i: %f\n", evaluation_point, dy);

%% Analytical Solution
l = 10;
I = 1/12;
w = -(q*l^4)/(76.8*E*I);
err = (w-dy)*100/dy;
fprintf("Analytical Solution for Simply Supported Beam: %f\n",w);
fprintf("Error: %3.2f %%\n", err);
%% Plot 

vis = Visualization(model);
% vis.plotUndeformed();
% vis.plotDeformed();
% vis.plotField('sigma_xx');
