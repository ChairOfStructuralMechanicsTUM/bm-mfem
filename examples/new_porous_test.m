%% Create bigger Model

clc;
clear all;

% Abmessungen des Körpers
Lx=5;
Ly=5;
Lz=0.25;

% Anzahl der Elemente in jew. Richtung
nx=20;
ny=20;
nz=1;

% Abmessungen der Elemente Berechnen (durch L und n definiert - keine
% Eingabe)
dx=Lx/nx;
dy=Ly/ny;
dz=Lz/nz;

numnodes=(nx + 1)*(ny + 1)*(nz + 1);
numele=nx*ny*nz;

model = FemModel();

% Erzeugen der Knoten
id=0;

for k=1:(ny + 1)
    for j=1:(nz + 1)
        for i=1:(nx + 1)
            id=id+1;
            model.addNewNode(id,(i-1)*dx,(k-1)*dy,(j-1)*dz);
        end
    end
end

% Zuweisung der DOF
model.getAllNodes.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y', 'DISPLACEMENT_SOLID_Z','DISPLACEMENT_FLUID_X', 'DISPLACEMENT_FLUID_Y', 'DISPLACEMENT_FLUID_Z'});

% Erzeugen der Elemente
id = 0;

for k=1:ny
    for j=1:nz
        for i=1:nx
            id=id+1;
            a = i + (j-1)*(nx+1) + (k-1)*(nx+1)*(nz+1);
            model.addNewElement('PorousElement3d8n',id,[a, a+1, a+1+(nx+1)*(nz+1),  a+(nx+1)*(nz+1), a+(nx+1), a+1+(nx+1), a+1+(nx+1)*(nz+1)+(nx+1), a+(nx+1)*(nz+1)+(nx+1) ]);
        end
    end
end

% Zuweisung der Materialparameter
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

model.getAllElements.setPropertyValue('FREQUENCY',10);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);

% Definition der Lagerungen
for i=1:(nx+1)*2
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_X');
    model.getNode(i).fixDof('DISPLACEMENT_FLUID_X');
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_Y');
    model.getNode(i).fixDof('DISPLACEMENT_FLUID_Y');
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_Z');
    model.getNode(i).fixDof('DISPLACEMENT_FLUID_Z');
end

% Definition der Lagerungen
for i=((nx+1)*(nz+1)*ny)+1 :numnodes
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_X');
    model.getNode(i).fixDof('DISPLACEMENT_FLUID_X');
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_Y');
    model.getNode(i).fixDof('DISPLACEMENT_FLUID_Y');
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_Z');
    model.getNode(i).fixDof('DISPLACEMENT_FLUID_Z');
end

% Definition der Last
a=(nx+1)*(nz+1)*ny/2 + (nx+1)*nz + nx/2 + 1;
addPointLoadPorous(model.getNode(a),1,[0 0 -1]);

% Bestimmen der globalen Systemmatrizen
assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);          
massMatrix = assembling.assembleGlobalMassMatrix(model);

% Lösen
%solver = SimpleSolvingStrategy(model);
solver = SimpleHarmonicSolvingStrategy(model,10);
x = solver.solve();
step = 1;
VerschiebungDofs = model.getDofArray.getValue(step);
nodalForces = solver.getNodalForces(step);

% Visualisierung der Lösung
v = Visualization(model);
v.setScaling(10000);
%v.plotUndeformed
v.plotDeformed
