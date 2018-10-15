%% PlaneStressElement Test with Gmsh Input File
%% Import file

io = ModelIO('tests/VM16_Trig6.msh');
model = io.readModel;
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

nodes = model.getAllNodes();
elements = model.getAllElements();

%% BC (Simply Supported)
nodes(1).fixDof('DISPLACEMENT_X');
nodes(1).fixDof('DISPLACEMENT_Y');
nodes(4).fixDof('DISPLACEMENT_X');
nodes(4).fixDof('DISPLACEMENT_Y');
% nodes(4).fixDof('DISPLACEMENT_Y');

%% Line BC (Cantiliver)

% idLowerLCorner = model.getModelPart('LLCorner').getId();
% idUpperLCorner = model.getModelPart('ULCorner').getId();
% 
% addLineBC(idLowerLCorner,idUpperLCorner,nodes);
% 
% 
% idLowerLCorner = model.getModelPart('LRCorner').getId();
% idUpperLCorner = model.getModelPart('URCorner').getId();
% 
% addLineBC(idLowerLCorner,idUpperLCorner,nodes);
%% Set Properties
E = 2e+11;
elements.setPropertyValue('THICKNESS', 0.1);
elements.setPropertyValue('YOUNGS_MODULUS', E);
elements.setPropertyValue('POISSON_RATIO', 0.3);
elements.setPropertyValue('NUMBER_GAUSS_POINT', 3);
elements.setPropertyValue('DENSITY', 7000);

%% Point load
% nodes(3).setDofLoad('DISPLACEMENT_X',-1000)
% nodes(2).setDofLoad('DISPLACEMENT_X',1000)
nodes(3).setDofLoad('DISPLACEMENT_Y',150)
nodes(2).setDofLoad('DISPLACEMENT_Y',150)
%% Line load
% q = 1000000;
% startNodeId = model.getModelPart('ULCorner').getId();
% endNodeId = model.getModelPart('URCorner').getId();
% 
% addConstLineLoad(startNodeId,endNodeId,nodes,q,[0 -1 0]);


%% Solving system
solver = SimpleSolvingStrategy(model);
solver.solve();
% actualDisplacementY = model.getAllNodes.getDofValue('DISPLACEMENT_Y');
% [maxDisplacement, idx] = max(actualDisplacementY);
% fprintf("Max Displacement: %e\n", maxDisplacement);

% step = 1;
% [stressValue, element_connect] = computeElementStress(elements,nodes,step);
% field = stressValue(1,:)';
% evaluation_point = model.getModelPart('LMPoint').getId();
% stress = field(12);
% dy = nodeArray(evaluation_point).getDofValue('DISPLACEMENT_Y');
% fprintf("Stress: %f\n", stress);

% %% Analytical Solution
% l = 10;
% I = 1/12;
% w = -(q*l^4)/(76.8*E*I);
% err = (w-dy)*100/dy;
% fprintf("Analytical Solution for Simply Supported Beam: %f\n",w);
% fprintf("Error: %3.2f %%\n", err);
%% Eigenfrequencies
% 
% eigensolver = EigensolverStrategy(model);
% eigensolver.solve(5);
% eigensolver.assignModeShapes();
% modes = eigensolver.getModalMatrix;
% eigenfrequencies = eigensolver.getEigenfrequencies();
% 
% fprintf("Eigenfrequencies: %f\n" , eigenfrequencies);
%% Plot 

vis = Visualization(model);
% vis.plotField('displacement_absolute',3);

vis.plotDeformed();
vis.plotField('sigma_xx');
vis.plotConstrain(1);
vis.plotLoad(1);
% vis.plotUndeformed();