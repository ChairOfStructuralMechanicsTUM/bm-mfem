%% Averaged and Unaveraged Stresses
clear;

%% Import file

io = ModelIO('tests/VM16_Trig3.msh');
model = io.readModel;
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

nodes = model.getAllNodes();
elements = model.getAllElements();

%% BC (Simply Supported)
nodes(1).fixDof('DISPLACEMENT_X');
nodes(1).fixDof('DISPLACEMENT_Y');
nodes(4).fixDof('DISPLACEMENT_X');
nodes(4).fixDof('DISPLACEMENT_Y');

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

%% Plot 

vis = Visualization(model);

% Plot Averaged Stress
vis.plotField('sigma_xx');

% Plot Unaveraged Stress
figure('Name','Unaveraged')
hold on
axis equal

for i = 1:length(elements)
    
X = elements(i).getNodes.getX();
Y = elements(i).getNodes.getY();
ord = elements(i).drawOrder();
field = sigma_xx(i,:);

fill(X(ord),Y(ord),field(ord))

end

hold off
colormap parula(9)
colorbar

