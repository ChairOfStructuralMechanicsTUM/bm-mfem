 %% SET PARAMETERS


%function [Answer,lambda1,lambda2,lambda3] = ModelTest(POROSITY, DENSITY_S, DENSITY_F, OMEGA, NUMBER_GAUSS_POINT)

function [DisplacementPorousNode1, DisplacementPorousNode2, DisplacementPorousNode3, DisplacementPorousNode4, DisplacementPorousTime] = ModelTestForGridCovergence(Lx,Ly,nx,ny)
clc

LoadValue = 10;
LoadDirection = [0 -1];


DENSITY_S = 30;
DENSITY_F = 1.21;
POROSITY = 0.96;
OMEGA = 100;
NUMBER_GAUSS_POINT = 2;
%Properties
p.DENSITY_S = DENSITY_S;
p.LAMBDA_S = 905357;
p.MUE_S = 264062;
p.ETA_S = 0;
p.DENSITY_F = DENSITY_F;
p.ETA_F = 1.84e-5;
p.PRESSURE_0_F = 101;
p.HEAT_CAPACITY_RATIO_F = 1.4;
p.PRANDL_NUMBER_F = 0.71;
p.POROSITY = POROSITY;
p.TORTUOSITY = 1.7;
p.STATIC_FLOW_RESISTIVITY = 32e3;
p.VISCOUS_CHARACT_LENGTH = 90;
p.THERMAL_CHARACT_LENGTH = 165;
p.OMEGA = omega;
p.NUMBER_GAUSS_POINT = NUMBER_GAUSS_POINT;


%% Standard (u_s,u_f) displacement
clc
clear model ans solver stiffnessMatrix massMatrix DisplacementAllard DisplacementAllardSolid
tic
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

model.getAllElements.setPropertyValue('DENSITY_S',p.DENSITY_S);
model.getAllElements.setPropertyValue('LAMBDA_S',p.LAMBDA_S);
model.getAllElements.setPropertyValue('MUE_S',p.MUE_S);
model.getAllElements.setPropertyValue('ETA_S',p.ETA_S);
% Fluid Parameters:
model.getAllElements.setPropertyValue('DENSITY_F',p.DENSITY_F);
model.getAllElements.setPropertyValue('ETA_F',p.ETA_F);
model.getAllElements.setPropertyValue('PRESSURE_0_F',p.PRESSURE_0_F);
model.getAllElements.setPropertyValue('HEAT_CAPACITY_RATIO_F',p.HEAT_CAPACITY_RATIO_F);
model.getAllElements.setPropertyValue('PRANDL_NUMBER_F',p.PRANDL_NUMBER_F);
% Porous Parameters:
model.getAllElements.setPropertyValue('POROSITY',p.POROSITY);
model.getAllElements.setPropertyValue('TORTUOSITY',p.TORTUOSITY);
model.getAllElements.setPropertyValue('STATIC_FLOW_RESISTIVITY',p.STATIC_FLOW_RESISTIVITY);
model.getAllElements.setPropertyValue('VISCOUS_CHARACT_LENGTH',p.VISCOUS_CHARACT_LENGTH);
model.getAllElements.setPropertyValue('THERMAL_CHARACT_LENGTH',p.THERMAL_CHARACT_LENGTH);
% General Parameters:
model.getAllElements.setPropertyValue('OMEGA',p.OMEGA);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',p.NUMBER_GAUSS_POINT);

for i = 1:1+nx:numnodes
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_X')
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_Y')
    model.getNode(i).fixDof('DISPLACEMENT_FLUID_X')
    model.getNode(i).fixDof('DISPLACEMENT_FLUID_Y') 
end

addPointLoadPorous(model.getNode(numnodes),LoadValue,LoadDirection); 
assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
massMatrix = assembling.assembleGlobalMassMatrix(model);

solver = SimpleHarmonicSolvingStrategy(model,p.OMEGA);
x = solver.solve();
step = 1;
Allardtime=toc;

 %v = Visualization(model);
 %v.setScaling(1);
 %v.plotUndeformed
 %v.plotDeformed;

%
%% Mixed (u_s,p) displacement
clc
clear model ans solver stiffnessMatrix massMatrix 
tic
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


model.getAllNodes.addDof({'PRESSURE_FLUID', 'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y'});

id = 0;

for j = 1:ny
    for i = 1:nx
        id = id+1;
        a=i + (j-1)*(nx+1);
        model.addNewElement('AtallaElement2d4n',id,[a,a+1,a+1+(nx+1),a+(nx+1)]);
    end
end

model.getAllElements.setPropertyValue('DENSITY_S',p.DENSITY_S);
model.getAllElements.setPropertyValue('LAMBDA_S',p.LAMBDA_S);
model.getAllElements.setPropertyValue('MUE_S',p.MUE_S);
model.getAllElements.setPropertyValue('ETA_S',p.ETA_S);
% Fluid Parameters:
model.getAllElements.setPropertyValue('DENSITY_F',p.DENSITY_F);
model.getAllElements.setPropertyValue('ETA_F',p.ETA_F);
model.getAllElements.setPropertyValue('PRESSURE_0_F',p.PRESSURE_0_F);
model.getAllElements.setPropertyValue('HEAT_CAPACITY_RATIO_F',p.HEAT_CAPACITY_RATIO_F);
model.getAllElements.setPropertyValue('PRANDL_NUMBER_F',p.PRANDL_NUMBER_F);
% Porous Parameters:
model.getAllElements.setPropertyValue('POROSITY',p.POROSITY);
model.getAllElements.setPropertyValue('TORTUOSITY',p.TORTUOSITY);
model.getAllElements.setPropertyValue('STATIC_FLOW_RESISTIVITY',p.STATIC_FLOW_RESISTIVITY);
model.getAllElements.setPropertyValue('VISCOUS_CHARACT_LENGTH',p.VISCOUS_CHARACT_LENGTH);
model.getAllElements.setPropertyValue('THERMAL_CHARACT_LENGTH',p.THERMAL_CHARACT_LENGTH);
% General Parameters:
model.getAllElements.setPropertyValue('OMEGA',p.OMEGA);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',p.NUMBER_GAUSS_POINT);

for i = 1:1+nx:numnodes
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_X')
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_Y') 
end

for j = 2:1+nx
   model.getNode(j).fixDof('PRESSURE_FLUID')
end

for j = 1+nx:1+nx:numnodes
   model.getNode(j).fixDof('PRESSURE_FLUID')
end

for j = numnodes-(nx-1):1:numnodes
   model.getNode(j).fixDof('PRESSURE_FLUID')
end

addPointLoadPorous(model.getNode(numnodes),LoadValue,LoadDirection);
assembling = SimpleAssembler(model);    
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);     
massMatrix = assembling.assembleGlobalMassMatrix(model);

solver = SimpleHarmonicSolvingStrategy(model,p.OMEGA);
x = solver.solve();
  
step = 1;
Atallatime = toc;

% 
 %v2 = Visualization(model);
 %v2.setScaling(1);
 %%v2.plotUndeformed
 %v2.plotDeformed;
% 
% 

%% Total (u_s,u_t) displacement

clc
clear model ans solver stiffnessMatrix massMatrix 
close all
tic

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

model.getAllElements.setPropertyValue('DENSITY_S',p.DENSITY_S);
model.getAllElements.setPropertyValue('LAMBDA_S',p.LAMBDA_S);
model.getAllElements.setPropertyValue('MUE_S',p.MUE_S);
model.getAllElements.setPropertyValue('ETA_S',p.ETA_S);
% Fluid Parameters:
model.getAllElements.setPropertyValue('DENSITY_F',p.DENSITY_F);
model.getAllElements.setPropertyValue('ETA_F',p.ETA_F);
model.getAllElements.setPropertyValue('PRESSURE_0_F',p.PRESSURE_0_F);
model.getAllElements.setPropertyValue('HEAT_CAPACITY_RATIO_F',p.HEAT_CAPACITY_RATIO_F);
model.getAllElements.setPropertyValue('PRANDL_NUMBER_F',p.PRANDL_NUMBER_F);
% Porous Parameters:
model.getAllElements.setPropertyValue('POROSITY',p.POROSITY);
model.getAllElements.setPropertyValue('TORTUOSITY',p.TORTUOSITY);
model.getAllElements.setPropertyValue('STATIC_FLOW_RESISTIVITY',p.STATIC_FLOW_RESISTIVITY);
model.getAllElements.setPropertyValue('VISCOUS_CHARACT_LENGTH',p.VISCOUS_CHARACT_LENGTH);
model.getAllElements.setPropertyValue('THERMAL_CHARACT_LENGTH',p.THERMAL_CHARACT_LENGTH);
% General Parameters:
model.getAllElements.setPropertyValue('OMEGA',p.OMEGA);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',p.NUMBER_GAUSS_POINT);


for i = 1:1+nx:numnodes
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_X')
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_Y')
    model.getNode(i).fixDof('DISPLACEMENT_TOTAL_X')
    model.getNode(i).fixDof('DISPLACEMENT_TOTAL_Y') 
end


addPointLoadPorous(model.getNode(numnodes),LoadValue,LoadDirection);
 
assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);     
massMatrix = assembling.assembleGlobalMassMatrix(model);

solver = SimpleHarmonicSolvingStrategy(model,p.OMEGA);
x = solver.solve(); 
Dazeltime=toc;
step = 1;

 %v3 = Visualization(model);
 %v3.setScaling(1);
 %%v3.plotUndeformed
 %v3.plotDeformed;

%% EVALUATION

uxAllard = modelAllard.getAllNodes.getDofValue('DISPLACEMENT_SOLID_X');

uyAllard = modelAllard.getAllNodes.getDofValue('DISPLACEMENT_SOLID_Y');

uxAtalla = modelAtalla.getAllNodes.getDofValue('DISPLACEMENT_SOLID_X');

uyAtalla = modelAtalla.getAllNodes.getDofValue('DISPLACEMENT_SOLID_Y');

uxDazel = modelDazel.getAllNodes.getDofValue('DISPLACEMENT_SOLID_X');

uyDazel = modelDazel.getAllNodes.getDofValue('DISPLACEMENT_SOLID_Y');

DisplacementPorousNode4(1,1)=uxAllard(nx+1);
DisplacementPorousNode4(1,2)=uyAllard(nx+1);
DisplacementPorousNode4(1,3)=uxAtalla(nx+1);
DisplacementPorousNode4(1,4)=uyAtalla(nx+1);
DisplacementPorousNode4(1,5)=uxDazel(nx+1);
DisplacementPorousNode4(1,6)=uyDazel(nx+1);

DisplacementPorousNode3(1,1)=uxAllard(round((nx+1)*3/4));
DisplacementPorousNode3(1,2)=uyAllard(round((nx+1)*3/4));
DisplacementPorousNode3(1,3)=uxAtalla(round((nx+1)*3/4));
DisplacementPorousNode3(1,4)=uyAtalla(round((nx+1)*3/4));
DisplacementPorousNode3(1,5)=uxDazel(round((nx+1)*3/4));
DisplacementPorousNode3(1,6)=uyDazel(round((nx+1)*3/4));

DisplacementPorousNode2(1,1)=uxAllard(round((nx+1)*2/4));
DisplacementPorousNode2(1,2)=uyAllard(round((nx+1)*2/4));
DisplacementPorousNode2(1,3)=uxAtalla(round((nx+1)*2/4));
DisplacementPorousNode2(1,4)=uyAtalla(round((nx+1)*2/4));
DisplacementPorousNode2(1,5)=uxDazel(round((nx+1)*2/4));
DisplacementPorousNode2(1,6)=uyDazel(round((nx+1)*2/4));

DisplacementPorousNode1(1,1)=uxAllard(round((nx+1)*1/4));
DisplacementPorousNode1(1,2)=uyAllard(round((nx+1)*1/4));
DisplacementPorousNode1(1,3)=uxAtalla(round((nx+1)*1/4));
DisplacementPorousNode1(1,4)=uyAtalla(round((nx+1)*1/4));
DisplacementPorousNode1(1,5)=uxDazel(round((nx+1)*1/4));
DisplacementPorousNode1(1,6)=uyDazel(round((nx+1)*1/4));

DisplacementPorousTime(1,1)=Allardtime;
DisplacementPorousTime(1,2)=Atallatime;
DisplacementPorousTime(1,3)=Dazeltime;



end
