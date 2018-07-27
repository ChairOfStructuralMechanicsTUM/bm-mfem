%% Create bigger Model

clc;
clear all;

% Dimension of the structure
Lx=5;
Ly=5;
Lz=0.25;

% Number of elements in specific directions
nx=20;
ny=20;
nz=1;

% Calculation of the dimension of the elements (defined via L and n)
dx=Lx/nx;
dy=Ly/ny;
dz=Lz/nz;

numnodes=(nx + 1)*(ny + 1)*(nz + 1);
numele=nx*ny*nz;

model = FemModel();

% Generation of nodes
id=0;

for k=1:(ny + 1)
    for j=1:(nz + 1)
        for i=1:(nx + 1)
            id=id+1;
            model.addNewNode(id,(i-1)*dx,(k-1)*dy,(j-1)*dz);
        end
    end
end

% Assignment DOFs
model.getAllNodes.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y', 'DISPLACEMENT_SOLID_Z', 'PORE_PRESSURE'});

% Generation of elements
id = 0;

for k=1:ny
    for j=1:nz
        for i=1:nx
            id=id+1;
            a = i + (j-1)*(nx+1) + (k-1)*(nx+1)*(nz+1);
            model.addNewElement('MixedPorousElement3d8n',id,[a, a+1, a+1+(nx+1)*(nz+1), a+(nx+1)*(nz+1), a+(nx+1), a+1+(nx+1), a+1+(nx+1)*(nz+1)+(nx+1), a+(nx+1)*(nz+1)+(nx+1)]);
        end
    end
end


% assignment of material properties
model.getAllElements.setPropertyValue('DENSITY_SOLID',30);
model.getAllElements.setPropertyValue('LAMBDA_SOLID',905357);
model.getAllElements.setPropertyValue('MUE_SOLID',264062);
model.getAllElements.setPropertyValue('DAMPING_SOLID',0);

model.getAllElements.setPropertyValue('DENSITY_FLUID',1.21);
model.getAllElements.setPropertyValue('VISCOSITY_FLUID',1.84e-5);
model.getAllElements.setPropertyValue('STANDARD_PRESSURE_FLUID',101);
model.getAllElements.setPropertyValue('HEAT_CAPACITY_FLUID',1.4);
model.getAllElements.setPropertyValue('PRANDTL_NUMBER_FLUID',0.71);

model.getAllElements.setPropertyValue('POROSITY',0.96);
model.getAllElements.setPropertyValue('TORTUOSITY',1.7);
model.getAllElements.setPropertyValue('FLOW_RESISTIVITY',32e3);
model.getAllElements.setPropertyValue('VISCOUS_LENGHT',90);
model.getAllElements.setPropertyValue('THERMAL_LENGTH',165);

model.getAllElements.setPropertyValue('FREQUENCY',100);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);

% Definition of BCs
for i=1:(nx+1)*2
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_X');
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_Y');
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_Z');
end

for i=((nx+1)*(nz+1)*ny)+1 :numnodes
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_X');
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_Y');
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_Z');
end

% Definition of loading
a=(nx+1)*(nz+1)*ny/2 + (nx+1)*nz + nx/2 + 1;
addPointLoadPorous(model.getNode(a),1,[0 0 -1]);

% Determination of global matrices
assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);          
%massMatrix = assembling.assembleGlobalMassMatrix(model);

% Solving
%solver = SimpleSolvingStrategy(model);
solver = SimpleHarmonicSolvingStrategy(model,100);
x = solver.solve();
step = 1;
VerschiebungDofs = model.getDofArray.getValue(step);
nodalForces = solver.getNodalForces(step);

% Visualisierung der Lösung
v = Visualization(model);
v.setScaling(100000);
% v.plotUndeformed
v.plotDeformed
