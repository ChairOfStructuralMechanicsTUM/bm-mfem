%% SET PARAMETERS
clc 
clear all

Lx = 5;
Ly = 0.5;

nx = 20;
ny = 4;

% If SolvingMetod = 1 -> Harmonic solver
% If SolvingMetod = 2 -> Newmark solver
SolvingMethod = 1;

LoadValue = 1;
LoadDirection = [0 -1];

%% Standard (u_s,u_f) displacement


model = FemModel();

dx = Lx/nx;
dy = Ly/ny;

numnodes = (nx + 1) * (ny+1);
numele = nx*ny;

id = 0;

for j = 1:(ny+1)
    for i= 1:(nx+1)
        id = id+1;
        model.addNewNode(id,(i-1)*dx,(j-1)*dy);
    end
end


model.getAllNodes.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y', 'DISPLACEMENT_FLUID_X', 'DISPLACEMENT_FLUID_Y'});

id = 0;

for j = 1:ny
    for i = 1:nx
        id = id+1;
        a=i + (j-1)*(nx+1);
        model.addNewElement('BiotAllardElement2d4n',id,[a,a+1,a+1+(nx+1),a+(nx+1)]);
    end
end

model.getAllElements.setPropertyValue('DENSITY_S',30);
model.getAllElements.setPropertyValue('LAMBDA_S',905357);
model.getAllElements.setPropertyValue('MUE_S',264062);
model.getAllElements.setPropertyValue('ETA_S',0);
% Fluid Parameters:
model.getAllElements.setPropertyValue('DENSITY_F',1.21);
model.getAllElements.setPropertyValue('ETA_F',1.84e-5);
model.getAllElements.setPropertyValue('PRESSURE_0_F',101);
model.getAllElements.setPropertyValue('HEAT_CAPACITY_RATIO_F',1.4);
model.getAllElements.setPropertyValue('PRANDL_NUMBER_F',0.71);
% Porous Parameters:
model.getAllElements.setPropertyValue('POROSITY',0.96);
model.getAllElements.setPropertyValue('TORTUOSITY',1.7);
model.getAllElements.setPropertyValue('STATIC_FLOW_RESISTIVITY',32e3);
model.getAllElements.setPropertyValue('VISCOUS_CHARACT_LENGTH',90);
model.getAllElements.setPropertyValue('THERMAL_CHARACT_LENGTH',165);
% General Parameters:
model.getAllElements.setPropertyValue('OMEGA',100);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);

for i = 1:(ny+1)
    model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_X')
    model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_Y')
    model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_FLUID_X')
    model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_FLUID_Y')
end


addPointLoadPorous(model.getNode(numnodes),LoadValue,LoadDirection);
 
assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);     
massMatrix = assembling.assembleGlobalMassMatrix(model);

SolveAndPlot(SolvingMethod, model)




%% Mixed (u_s,p) displacement
clc
clear model ans solver stiffnessMatrix massMatrix 

model = FemModel();

dx = Lx/nx;
dy = Ly/ny;

numnodes = (nx + 1) * (ny+1);
numele = nx*ny;

id = 0;

for j = 1:(ny+1)
    for i= 1:(nx+1)
        id = id+1;
        model.addNewNode(id,(i-1)*dx,(j-1)*dy);
    end
end


model.getAllNodes.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y', 'PRESSURE_FLUID'});

id = 0;

for j = 1:ny
    for i = 1:nx
        id = id+1;
        a=i + (j-1)*(nx+1);
        model.addNewElement('AtallaElement2d4n',id,[a,a+1,a+1+(nx+1),a+(nx+1)]);
    end
end

model.getAllElements.setPropertyValue('DENSITY_S',30);
model.getAllElements.setPropertyValue('LAMBDA_S',905357);
model.getAllElements.setPropertyValue('MUE_S',264062);
model.getAllElements.setPropertyValue('ETA_S',0);
% Fluid Parameters:
model.getAllElements.setPropertyValue('DENSITY_F',1.21);
model.getAllElements.setPropertyValue('ETA_F',1.84e-5);
model.getAllElements.setPropertyValue('PRESSURE_0_F',101);
model.getAllElements.setPropertyValue('HEAT_CAPACITY_RATIO_F',1.4);
model.getAllElements.setPropertyValue('PRANDL_NUMBER_F',0.71);
% Porous Parameters:
model.getAllElements.setPropertyValue('POROSITY',0.96);
model.getAllElements.setPropertyValue('TORTUOSITY',1.7);
model.getAllElements.setPropertyValue('STATIC_FLOW_RESISTIVITY',32e3);
model.getAllElements.setPropertyValue('VISCOUS_CHARACT_LENGTH',90);
model.getAllElements.setPropertyValue('THERMAL_CHARACT_LENGTH',165);
% General Parameters:
model.getAllElements.setPropertyValue('OMEGA',100);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);

for i = 1:(ny+1)
    model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_X')
    model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_Y')
    model.getNode(i+(i-1)*(nx+1)).fixDof('PRESSURE_FLUID')
   
end


addPointLoadPorous(model.getNode(numnodes),LoadValue,LoadDirection);
 
assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);     
massMatrix = assembling.assembleGlobalMassMatrix(model);

SolveAndPlot(SolvingMethod, model)



%% Total (u_s,u_t) displacement

clc
clear model ans solver stiffnessMatrix massMatrix 


model = FemModel();

dx = Lx/nx;
dy = Ly/ny;

numnodes = (nx + 1) * (ny+1);
numele = nx*ny;

id = 0;

for j = 1:(ny+1)
    for i= 1:(nx+1)
        id = id+1;
        model.addNewNode(id,(i-1)*dx,(j-1)*dy);
    end
end


model.getAllNodes.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y', 'DISPLACEMENT_TOTAL_X', 'DISPLACEMENT_TOTAL_Y'});

id = 0;

for j = 1:ny
    for i = 1:nx
        id = id+1;
        a=i + (j-1)*(nx+1);
        model.addNewElement('DazelElement2d4n',id,[a,a+1,a+1+(nx+1),a+(nx+1)]);
    end
end

model.getAllElements.setPropertyValue('DENSITY_S',30);
model.getAllElements.setPropertyValue('LAMBDA_S',905357);
model.getAllElements.setPropertyValue('MUE_S',264062);
model.getAllElements.setPropertyValue('ETA_S',0);
% Fluid Parameters:
model.getAllElements.setPropertyValue('DENSITY_F',1.21);
model.getAllElements.setPropertyValue('ETA_F',1.84e-5);
model.getAllElements.setPropertyValue('PRESSURE_0_F',101);
model.getAllElements.setPropertyValue('HEAT_CAPACITY_RATIO_F',1.4);
model.getAllElements.setPropertyValue('PRANDL_NUMBER_F',0.71);
% Porous Parameters:
model.getAllElements.setPropertyValue('POROSITY',0.96);
model.getAllElements.setPropertyValue('TORTUOSITY',1.7);
model.getAllElements.setPropertyValue('STATIC_FLOW_RESISTIVITY',32e3);
model.getAllElements.setPropertyValue('VISCOUS_CHARACT_LENGTH',90);
model.getAllElements.setPropertyValue('THERMAL_CHARACT_LENGTH',165);
% General Parameters:
model.getAllElements.setPropertyValue('OMEGA',100);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);

for i = 1:(ny+1)
    model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_X')
    model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_Y')
    model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_TOTAL_X')
    model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_TOTAL_Y')
end


addPointLoadPorous(model.getNode(numnodes),LoadValue,LoadDirection);
 
assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);     
massMatrix = assembling.assembleGlobalMassMatrix(model);

SolveAndPlot(SolvingMethod, model)

%% EVALUATION


% Welche ids entsprechen den solid displacements??

SolidDisplacements(:,1) = VerschiebungDofs_biot;
SolidDisplacements(:,2) = VerschiebungDofs_mixed;
SolidDisplacements(:,3) = VerschiebungDofs_total;

