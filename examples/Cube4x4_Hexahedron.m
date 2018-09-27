clc;
clear all;

% %  ========================= %
% %% GEOMETRY BY MANUELL INPUT %
% %  ========================= %
% 
% % Assignment of Nodes
% node01 = Node(1,0,0,0);
% node02 = Node(2,1,0,0);
% node03 = Node(3,2,0,0);
% node04 = Node(4,3,0,0);
% node05 = Node(5,4,0,0);
% node06 = Node(6,0,1,0);
% node07 = Node(7,1,1,0);
% node08 = Node(8,2,1,0);
% 
% nodeArray = [node01 node02 node03 node04 node05 node06 node07 node08];
% 
% 
% % Assignment DOFs
% nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});
% 
% 
% % Assignment of Elements
% ele01 = HexahedronElement3d8n(1,nodeArray);
% 
% 
% % Assignment of Material Properties
% ele01.setPropertyValue('YOUNGS_MODULUS',7e10);
% ele01.setPropertyValue('POISSON_RATIO',0.34);
% ele01.setPropertyValue('NUMBER_GAUSS_POINT',2);
% ele01.setPropertyValue('DENSITY',2699);
% 
% 
% % Assignment of BCs
% node01.fixAllDofs();
% node02.fixAllDofs();
% node03.fixAllDofs();
% node04.fixAllDofs();
% 
% 
% % Definition of Loads
% addPointLoad([node05 node06],100,[0 -1 0]);
% 
% 
% % Definition of the femModel
% model = FemModel(nodeArray, ele01);



% % ================================== %
% %% GEOMETRY FOR CUBOIDS AUTOMATICALLY %
% % ================================== %

% Dimension of the structure
Lx=1;
Ly=1;
Lz=1;

% Number of elements in specific directions
nx=2;
ny=2;
nz=2;

% Calculation of the dimension of the elements (defined via L and n)
dx=Lx/nx;
dy=Ly/ny;
dz=Lz/nz;

numnodes=(nx + 1)*(ny + 1)*(nz + 1);
numele=nx*ny*nz;

model = FemModel();

% Generation of Nodes
id=0;

for k=1:(ny + 1)
    for j=1:(nz + 1)
        for i=1:(nx + 1)
            id=id+1;
            model.addNewNode(id,(i-1)*dx,(k-1)*dy,(j-1)*dz);
        end
    end
end

%Knotennummerierung: Erst in x-Richtung, dann in z-Richtung 1 hoch, wieder
%x-Richtung, wenn alle z-Nummeriert, y+1 und wie gehabt bei z=0 startend in
%x-Richtung.

% Assignment DOFs
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z'});

% Generation of Elements
id = 0;

for k=1:ny
    for j=1:nz
        for i=1:nx
            id=id+1;
            a = i + (j-1)*(nx+1) + (k-1)*(nx+1)*(nz+1);
            model.addNewElement('HexahedronElement3d8n',id,[a, a+1, a+1+(nx+1)*(nz+1), a+(nx+1)*(nz+1), a+(nx+1), a+1+(nx+1), a+1+(nx+1)*(nz+1)+(nx+1), a+(nx+1)*(nz+1)+(nx+1)]);
        end
    end
end

% Assignment of Material Properties
model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e8);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
model.getAllElements.setPropertyValue('DENSITY',2699);


% Definition of BCs
% for i=1:(nx+1)*2
%     model.getNode(i).fixDof('DISPLACEMENT_X');
%     model.getNode(i).fixDof('DISPLACEMENT_Y');
%     model.getNode(i).fixDof('DISPLACEMENT_Z');
% end
% 
% for i=((nx+1)*(nz+1)*ny)+1 :numnodes
%     model.getNode(i).fixDof('DISPLACEMENT_X');
%     model.getNode(i).fixDof('DISPLACEMENT_Y');
%     model.getNode(i).fixDof('DISPLACEMENT_Z');
% end
model.getNodes([1 2 3 12 21 20 19 10]).fixAllDofs();

% Definition of Loads
a=(nx+1)*(nz+1)*ny/2 + (nx+1)*nz + nx/2 + 1;
addPointLoad(model.getNodes([7 25]),1000,[0 0 -1]);




% ====================================== %
%% GEOMETRY INPUT FROM GiD (.mdpa-input) %
% ====================================== %


% % MDPA_INPUT Read a large mdpa model and plot it
% io=MdpaInput('Plate_neu.mdpa'); %specify input file
% model = io.readModel(); %read the model
% 
% 
% % Assignment DOFs
% model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y", "DISPLACEMENT_Z"]);
% 
% 
% % Assignment of Material Properties
% model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e8);
% model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
% model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
% model.getAllElements.setPropertyValue('DENSITY',2699);
% 
% 
% % Set BCs for specified elements:
% % supportedNodes = [1249, 1259, 1263, 1275, 1285, 1295, 1311, 1319, 1325, 1329, 1331, 644, 661, 684, 731, 795, 862, 955, 1037, 1112, 1194, 1251];
% model.getModelPart('GENERIC_support').getNodes().fixAllDofs();
% % model.getNode(supportedNodes).fixAllDofs();
% % barycenter(model.getElement(1)); is not working for arbitrary hexahedra elements
% 
% % % % % % Set BCs for all elements in 'left_support' and 'right_support' (modelParts):
% % % % % %set boundary conditions for all elements in 'left_support' and 'right_support':
% % % % % model.getModelPart('GENERIC_left_support').getNodes().fixAllDofs();
% % % % % model.getModelPart('GENERIC_right_support').getNodes().fixAllDofs();
% 
% 
% % Set load for specified elements:
% loads = model.getModelPart('GENERIC_load').getNodes();
% addPointLoad(loads,1000,[0 0 -1]);
% % model.getNode(252).setDofLoad('DISPLACEMENT_X',100000);
% 
% %Important Remark: Visualization with Label is not possible.
% %ToDo: Calculate System with Hole




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

% Visualisierung der L�sung
v = Visualization(model);
v.setScaling(10000);
v.plotUndeformed(1)    %Insert 0: NoLabel / 1: WithLabels
v.plotDeformed


