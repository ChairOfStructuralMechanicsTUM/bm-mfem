clc;
clear all;

% EXAMPLES for FEM-Calculations with 3D Tetrahedron Elements
% TetrahedronElement3d4n: Tetraeder with 4 nodes in 3D (linear)


% CHOICE of Input-Method %

method = 'manuell';
% method = 'mdpa-input';


if strcmp(method, 'manuell')


%  ========================= %
% GEOMETRY BY MANUELL INPUT %
%  ========================= %

% Assignment of Nodes
node07 = Node(7,0,0,0);
node08 = Node(8,1,0,0);
node06 = Node(6,1,1,0);
node03 = Node(3,0,1,0);
node02 = Node(2,0,0,1);
node05 = Node(5,1,0,1);
node04 = Node(4,1,1,1);
node01 = Node(1,0,1,1);


nodeArray = [node01 node02 node03 node04 node05 node06 node07 node08];


% Assignment DOFs
nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});


% Assignment of Elements
ele01 = TetrahedronElement3d4n(1,[node07 node03 node01 node06]);
ele02 = TetrahedronElement3d4n(2,[node06 node04 node01 node07]);
ele03 = TetrahedronElement3d4n(3,[node02 node01 node04 node07]);
ele04 = TetrahedronElement3d4n(4,[node07 node02 node05 node04]);
ele05 = TetrahedronElement3d4n(5,[node07 node08 node06 node04]);
ele06 = TetrahedronElement3d4n(6,[node08 node05 node04 node07]);

elementArray = [ele01 ele02 ele03 ele04 ele05 ele06];


% Assignment of Material Properties
elementArray.setPropertyValue('YOUNGS_MODULUS',7e8);
elementArray.setPropertyValue('POISSON_RATIO',0.34);
% elementArray.setPropertyValue('NUMBER_GAUSS_POINT',1); 
% --> No GaussPointIntegrationScheme implemented up to now
elementArray.setPropertyValue('DENSITY',2699);


% Assignment of BCs
node03.fixAllDofs();
node06.fixAllDofs();
node07.fixAllDofs();
node08.fixAllDofs();


% Definition of Loads
addPointLoad([node02 node05],1000,[0 0 -1]);


% Definition of the femModel
model = FemModel(nodeArray, elementArray);



elseif strcmp(method, 'mdpa-input')

% ====================================== %
% GEOMETRY INPUT FROM GiD (.mdpa-input) %
% ====================================== %


% MDPA_INPUT Read a large mdpa model and plot it
io=MdpaInput('Cube_Tetrahedron.mdpa'); %specify input file
model = io.readModel(); %read the model


% Assignment DOFs
model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y", "DISPLACEMENT_Z"]);


% Assignment of Material Properties
model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e8);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
% model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',2699);


% Set BCs for specified elements:
% supportedNodes = [1249, 1259, 1263, 1275, 1285, 1295, 1311, 1319, 1325, 1329, 1331, 644, 661, 684, 731, 795, 862, 955, 1037, 1112, 1194, 1251];
model.getModelPart('GENERIC_FixedNodes').getNodes().fixAllDofs();
% model.getNode(supportedNodes).fixAllDofs();
% barycenter(model.getElement(1)); %is not working for arbitrary hexahedra elements
% 
% % % % % % Set BCs for all elements in 'left_support' and 'right_support' (modelParts):
% % % % % %set boundary conditions for all elements in 'left_support' and 'right_support':
% % % % % model.getModelPart('GENERIC_left_support').getNodes().fixAllDofs();
% % % % % model.getModelPart('GENERIC_right_support').getNodes().fixAllDofs();
% 
% 
% Set load for specified elements:
% loads = model.getModelPart('GENERIC_load').getNodes();
% addPointLoad(loads,1000,[0 0 -1]);
model.getNodes([5,2]).setDofLoad('DISPLACEMENT_Z',1000);
% 
% %Important Remark: Visualization with Label is not possible.
% %ToDo: Calculate System with Hole
% %ToDo: Visualization of Stresses

end


% ====================================== %
%% FE-CALCULATIONS (for all input dates) %
% ====================================== %

% Determination of Global Matrices
assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
forceVector = assembling.applyExternalForces(model);
% reducedForceVector = assembling.reducedForceVector;
% massMatrix = assembling.assembleGlobalMassMatrix(model);

% Solving
solver = SimpleSolvingStrategy(model);
% solver = SimpleHarmonicSolvingStrategy(model,100);
x = solver.solve();


% step = 1;           % timestep

% VerschiebungDofs = model.getDofArray.getValue(step);
verschiebung_z = model.getAllNodes().getDofValue('DISPLACEMENT_Z');
% nodalForces = solver.getNodalForces(step);

% Visualisierung der Lösung
v = Visualization(model);
v.setScaling(10000);
v.plotUndeformed(1)    %Insert 0: NoLabel / 1: WithLabels
v.plotDeformed