%% Clasical
clc;
clear all;

model = FemModel();

% Dimension of the structure
Lx=2;
Ly=1;

% Number of elements in specific directions
nx=20;
ny=10;

% Calculation of the dimension of the elements (defined via L and n)
dx=Lx/nx;
dy=Ly/ny;

numnodes=(nx + 1)*(ny + 1);
numele=nx*ny;

% Generation of nodes
id=0;

for j=1:(ny + 1)
    for i=1:(nx + 1)
        id=id+1;
        model.addNewNode(id,(i-1)*dx,(j-1)*dy);
    end
end

% Assignment DOFs
model.getAllNodes.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y', 'DISPLACEMENT_FLUID_X', 'DISPLACEMENT_FLUID_Y'});

% Generation of elements
id = 0;


for j=1:ny
    for i=1:nx
        id=id+1;
        a = i + (j-1)*(nx+1);
        model.addNewElement('ClassicalPorousElement2d4n',id,[a, a+1, a+1+(nx+1), a+(nx+1)]);
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

model.getAllElements.setPropertyValue('FREQUENCY',510);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);

% Definition of BCs
for i=1:(ny+1)
    model.getNode(1+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_X');
    model.getNode(1+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_Y');
    model.getNode(1+(i-1)*(nx+1)).fixDof('DISPLACEMENT_FLUID_X');
    model.getNode(1+(i-1)*(nx+1)).fixDof('DISPLACEMENT_FLUID_Y');
end

% Definition of loading
addPointLoadPorous(model.getNode(numnodes),1,[0 -1]);

% for k=1:200
%     model.getAllElements.setPropertyValue('FREQUENCY',k)
%     assembling = SimpleAssembler(model);
%     [~, Kred] = SimpleAssembler.assembleGlobalStiffnessMatrix(model);
%     [~, Mred] = SimpleAssembler.assembleGlobalMassMatrix(model);
%     stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
%     massMatrix = assembling.assembleGlobalMassMatrix(model);
%     if k <= 30
%         scale=1e-5;
% %     elseif (4 <= k) && (k <= 30)
% %         scale=1e-5;
%     else
%         scale=1e-6;
%     end
%     a(k,1)=det(scale*(Kred-k^2*Mred));
% end

% Determination of global matrices
assembling = SimpleAssembler(model);
[stiffnessMatrix, Kred] = assembling.assembleGlobalStiffnessMatrix(model);           
[massMatrix, Mred] = assembling.assembleGlobalMassMatrix(model);
  
% Solving
solver = SimpleHarmonicSolvingStrategy(model,510);
x = solver.solve();

step = 1;

VerschiebungDofs = model.getDofArray.getValue(step);

for i=1:numnodes
    Displacement_Solid (1+(i-1)*2,1) = VerschiebungDofs(1+(i-1)*4);
    Displacement_Solid (2+(i-1)*2,1) = VerschiebungDofs(2+(i-1)*4);
end

for i=1:2*numnodes
    Betrag(i,1)=sqrt(real(Displacement_Solid(i))^2 + imag(Displacement_Solid(i))^2);
    Phase(i,1)=atan(imag(Displacement_Solid(i))/real(Displacement_Solid(i)));
end

nodalForces = solver.getNodalForces(step);

v = Visualization(model);
v.setScaling(3e4);
%v.plotUndeformed
v.plotDeformed
    
%% Total

clc;
clear all;

model = FemModel();

% Dimension of the structure
Lx=5;
Ly=5;

% Number of elements in specific directions
nx=10;
ny=10;

% Calculation of the dimension of the elements (defined via L and n)
dx=Lx/nx;
dy=Ly/ny;

numnodes=(nx + 1)*(ny + 1);
numele=nx*ny;

% Generation of nodes
id=0;

for j=1:(ny + 1)
    for i=1:(nx + 1)
        id=id+1;
        model.addNewNode(id,(i-1)*dx,(j-1)*dy);
    end
end

% Assignment DOFs
model.getAllNodes.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y', 'DISPLACEMENT_TOTAL_X', 'DISPLACEMENT_TOTAL_Y'});

% Generation of elements
id = 0;


for j=1:ny
    for i=1:nx
        id=id+1;
        a = i + (j-1)*(nx+1);
        model.addNewElement('TotalPorousElement2d4n',id,[a, a+1, a+1+(nx+1), a+(nx+1)]);
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

% model.getAllElements.setPropertyValue('FREQUENCY',200);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);

% Definition of BCs
for i=1:(ny+1)
    model.getNode(1+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_X');
    model.getNode(1+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_Y');
    model.getNode(1+(i-1)*(nx+1)).fixDof('DISPLACEMENT_TOTAL_X');
    model.getNode(1+(i-1)*(nx+1)).fixDof('DISPLACEMENT_TOTAL_Y');
end

% Definition of loading
addPointLoadPorous(model.getNode(numnodes),1,[0 -1]);

for k=1:200
    model.getAllElements.setPropertyValue('FREQUENCY',k)
    assembling = SimpleAssembler(model);
    [~, Kred] = SimpleAssembler.assembleGlobalStiffnessMatrix(model);
    [~, Mred] = SimpleAssembler.assembleGlobalMassMatrix(model);
    stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
    massMatrix = assembling.assembleGlobalMassMatrix(model);
    if k <= 1
        scale=1e-4;
    elseif (2 <= k) && (k <= 38)
        scale=1e-5;
    else
        scale=1e-6;
    end
    a(k,1)=det(scale*(Kred-k^2*Mred));
end

% % Determination of global matrices
% assembling = SimpleAssembler(model);
% stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);           
% massMatrix = assembling.assembleGlobalMassMatrix(model);
% 
% % Solving
% solver = SimpleHarmonicSolvingStrategy(model,200);
% x = solver.solve();
% 
% step = 1;
% 
% VerschiebungDofs = model.getDofArray.getValue(step);
% 
% for i=1:numnodes
%     Displacement_Solid (1+(i-1)*2,1) = VerschiebungDofs(1+(i-1)*4);
%     Displacement_Solid (2+(i-1)*2,1) = VerschiebungDofs(2+(i-1)*4);
% end
% 
% for i=1:2*numnodes
%     Betrag(i,1)=sqrt(real(Displacement_Solid(i))^2 + imag(Displacement_Solid(i))^2);
%     Phase(i,1)=atan(imag(Displacement_Solid(i))/real(Displacement_Solid(i)));
% end
% 
% nodalForces = solver.getNodalForces(step);
% 
% v = Visualization(model);
% v.setScaling(1e5);
% %v.plotUndeformed
% v.plotDeformed

%% Mixed

clc;
clear all;

model = FemModel();

% Dimension of the structure
Lx=5;
Ly=5;

% Number of elements in specific directions
nx=10;
ny=10;

% Calculation of the dimension of the elements (defined via L and n)
dx=Lx/nx;
dy=Ly/ny;

numnodes=(nx + 1)*(ny + 1);
numele=nx*ny;

% Generation of nodes
id=0;

for j=1:(ny + 1)
    for i=1:(nx + 1)
        id=id+1;
        model.addNewNode(id,(i-1)*dx,(j-1)*dy);
    end
end

% Assignment DOFs
model.getAllNodes.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y', 'PORE_PRESSURE'});

% Generation of elements
id = 0;


for j=1:ny
    for i=1:nx
        id=id+1;
        a = i + (j-1)*(nx+1);
        model.addNewElement('MixedPorousElement2d4n',id,[a, a+1, a+1+(nx+1), a+(nx+1)]);
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

% model.getAllElements.setPropertyValue('FREQUENCY',200);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);

% Definition of BCs
for i=1:(ny+1)
    model.getNode(1+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_X');
    model.getNode(1+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_Y');
end

for i=1:(nx+1)
    model.getNode(i).fixDof('PORE_PRESSURE');
end

for i=1:(ny+1)
    model.getNode(i*(nx+1)).fixDof('PORE_PRESSURE');
end

a=(nx+1)*ny+1;
b=numnodes;
for i=a:b
    model.getNode(i).fixDof('PORE_PRESSURE');
end

% Definition of loading
addPointLoadPorous(model.getNode(numnodes),1,[0 -1]);

for k=1:200
    model.getAllElements.setPropertyValue('FREQUENCY',k)
    assembling = SimpleAssembler(model);
    [~, Kred] = SimpleAssembler.assembleGlobalStiffnessMatrix(model);
    [~, Mred] = SimpleAssembler.assembleGlobalMassMatrix(model);
    stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
    massMatrix = assembling.assembleGlobalMassMatrix(model);
    a(k,1)=det(1e-4*(Kred-k^2*Mred));
end

% % Determination of global matrices
% assembling = SimpleAssembler(model);
% stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);           
% massMatrix = assembling.assembleGlobalMassMatrix(model);
% 
% % Solving
% solver = SimpleHarmonicSolvingStrategy(model,200);
% x = solver.solve();
% 
% step = 1;
% 
% VerschiebungDofs = model.getDofArray.getValue(step);
% 
% for i=1:numnodes
%     a=1+(i-1)*2
%     b=2+(i-1)*2
%     c=1+(i-1)*3
%     d=2+(i-1)*3
%     Displacement_Solid (1+(i-1)*2,1) = VerschiebungDofs(1+(i-1)*3);
%     Displacement_Solid (2+(i-1)*2,1) = VerschiebungDofs(2+(i-1)*3);
% end
% 
% for i=1:2*numnodes
%     Betrag(i,1)=sqrt(real(Displacement_Solid(i))^2 + imag(Displacement_Solid(i))^2);
%     Phase(i,1)=atan(imag(Displacement_Solid(i))/real(Displacement_Solid(i)));
% end
% 
% nodalForces = solver.getNodalForces(step);
% 
% v = Visualization(model);
% v.setScaling(1e5);
% %v.plotUndeformed
% v.plotDeformed
