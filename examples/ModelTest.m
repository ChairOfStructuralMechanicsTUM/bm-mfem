 %% SET PARAMETERS


%function [Answer,lambda1,lambda2,lambda3] = ModelTest(POROSITY, DENSITY_S, DENSITY_F, OMEGA, NUMBER_GAUSS_POINT)
%function
clc

Lx = 0.05;
Ly = 0.01;
ny = 10;
nx = 50;


%Ideal:
% Lx = 0.05;
% Ly = 0.01;
% ny = 20;
% nx = 100;

%Laenger(1)v1:(Error: nE=4000)
% Lx = 0.1;
% Ly = 0.01;
% ny = 20;
% nx = 200;

%Laenger(1)v2:(Working aber error wavelength)
% Lx = 0.1;
% Ly = 0.01; 
% ny = 15;
% nx = 150;

%Schmaler(1):
%  Lx = 0.05;
%  Ly = 0.005;
%  ny = 10;
%  nx = 100;

%Laenger(2):(nE = 4500 -> wsl Errormeldung)
% Lx = 0.2;
% Ly = 0.01;
% ny = 15;
% nx = 300;

%Schmaler(2):
%  Lx = 0.005;
%  Ly = 0.00025;
%  ny = 10;
%  nx = 200;

LoadValue = 1;
LoadDirection = [0 -1];
loadposition = (nx+1)*(ny+1);

DENSITY_S = 30;
DENSITY_F = 1.21;
POROSITY = 0.96;
OMEGA = 1000;
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
p.OMEGA = OMEGA;
p.NUMBER_GAUSS_POINT = NUMBER_GAUSS_POINT;

ratio=WavelengthCheck(p.OMEGA,p.TORTUOSITY,p.ETA_F,p.DENSITY_F,DENSITY_S,p.STATIC_FLOW_RESISTIVITY,p.VISCOUS_CHARACT_LENGTH,p.POROSITY,p.HEAT_CAPACITY_RATIO_F,p.PRESSURE_0_F,...
    p.PRANDL_NUMBER_F,p.THERMAL_CHARACT_LENGTH,p.ETA_S,p.LAMBDA_S,p.MUE_S,Lx,Ly,nx,ny);



%% Classical (u_s,u_f) displacement
clc
clear model ans solver stiffnessMatrix massMatrix DisplacementAllard DisplacementAllardSolid
fprintf("Processing classical.. \n")
model = FemModel();
tic
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

% for i = 1:(ny+1)
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_X')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_Y')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_FLUID_X')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_FLUID_Y')
% end


for i = 1:1+nx:numnodes
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_X')
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_Y')
    model.getNode(i).fixDof('DISPLACEMENT_FLUID_X')
    model.getNode(i).fixDof('DISPLACEMENT_FLUID_Y') 
end



addPointLoadPorous(model.getNode(loadposition),LoadValue,LoadDirection); 
assembling = SimpleAssembler(model);
fprintf('Assembling Stiffness \n')
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
fprintf('Assembling Mass \n')
massMatrix = assembling.assembleGlobalMassMatrix(model);

fprintf('Solving \n')
solver = SimpleHarmonicSolvingStrategy(model,p.OMEGA);
x = solver.solve();
step = 1;
fprintf('Solved \n')
%DisplacementAllard = sparse(nx*ny*4,1);
DisplacementAllard(:,1) = model.getDofArray.getValue(step);
DisplacementAllard(:,2) = real(model.getDofArray.getValue(step)); 
DisplacementAllard(:,3) = imag(model.getDofArray.getValue(step));
DisplacementAllard(:,4) = sqrt((DisplacementAllard(:,2)).^2 + (DisplacementAllard(:,3)).^2);
DisplacementAllard(:,5) = atan(DisplacementAllard(:,3)./DisplacementAllard(:,2));

%DisplacementAllardSolid = zeros(nx*ny*2,1);
u=0;
for i=1:4:size(DisplacementAllard,1)
    u=u+1;
    DisplacementAllardSolid(u,1)=DisplacementAllard(i,1); %Solid_x_Komplex
    DisplacementAllardSolid(u,2)=DisplacementAllard(i,4); %Solid_x_Betrag   
    DisplacementAllardSolid(u,3)=DisplacementAllard(i,5); %Solid_x_Phasenwinkel
end
u=0;

for i=2:4:size(DisplacementAllard,1)
    u=u+1;
    DisplacementAllardSolid(u,4)=DisplacementAllard(i,1); %Solid_y_Komplex
    DisplacementAllardSolid(u,5)=DisplacementAllard(i,4); %Solid_y_Betrag   
    DisplacementAllardSolid(u,6)=DisplacementAllard(i,5); %Solid_y_Phasenwinkel
end

clear u

%nodalForces_biot = (solver.getNodalForces(step));
 

 v = Visualization(model);
 v.setScaling(1);
 %v.plotUndeformed
 v.plotDeformed;
time = toc
fprintf("Finished")

%% Mixed (u_s,p) displacement
clc
clear model ans solver stiffnessMatrix massMatrix DisplacementAtalla DisplacementAtallaSolid
fprintf("Processing mixed.. \n")
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


% for i = 1:(ny+1)
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_X')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_Y') 
% end

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


addPointLoadPorous(model.getNode(loadposition),LoadValue,LoadDirection);
 
assembling = SimpleAssembler(model);    
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);     
massMatrix = assembling.assembleGlobalMassMatrix(model);

%SolveAndPlot(model)

solver = SimpleHarmonicSolvingStrategy(model,p.OMEGA);
x = solver.solve();
   
step = 1;

%nodalForces_atalla = solver.getNodalForces(step);
  
%VerschiebungDofs_atalla = (model.getDofArray.getValue(step));

DisplacementAtalla(:,1) = model.getDofArray.getValue(step);
DisplacementAtalla(:,2) = real(model.getDofArray.getValue(step)); 
DisplacementAtalla(:,3) = imag(model.getDofArray.getValue(step));
DisplacementAtalla(:,4) = sqrt((DisplacementAtalla(:,2)).^2 + (DisplacementAtalla(:,3)).^2);
DisplacementAtalla(:,5) = atan(DisplacementAtalla(:,3)./DisplacementAtalla(:,2));


u=0;
for i=1:3:size(DisplacementAtalla,1)
    u=u+1;
    DisplacementAtallaSolid(u,1)=DisplacementAtalla(i,1); %Solid_x_Komplex
    DisplacementAtallaSolid(u,2)=DisplacementAtalla(i,4); %Solid_x_Betrag   
    DisplacementAtallaSolid(u,3)=DisplacementAtalla(i,5); %Solid_x_Phasenwinkel
end
u=0;
for i=2:3:size(DisplacementAtalla,1)
    u=u+1;
    DisplacementAtallaSolid(u,4)=DisplacementAtalla(i,1); %Solid_y_Komplex
    DisplacementAtallaSolid(u,5)=DisplacementAtalla(i,4); %Solid_y_Betrag   
    DisplacementAtallaSolid(u,6)=DisplacementAtalla(i,5); %Solid_y_Phasenwinkel
end

clear u

% 
 v2 = Visualization(model);
 v2.setScaling(1e3);
 %v2.plotUndeformed
 v2.plotDeformed;
% 
% 

fprintf("Finished")


%% Total (u_s,u_t) displacement

clc
clear model ans solver stiffnessMatrix massMatrix  DisplacementDazel DisplacementDazelSolid
close all
fprintf("Processing total.. \n")

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

% for i = 1:(ny+1)
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_X')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_Y')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_TOTAL_X')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_TOTAL_Y')
% end


for i = 1:1+nx:numnodes
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_X')
    model.getNode(i).fixDof('DISPLACEMENT_SOLID_Y')
    model.getNode(i).fixDof('DISPLACEMENT_TOTAL_X')
    model.getNode(i).fixDof('DISPLACEMENT_TOTAL_Y') 
end




addPointLoadPorous(model.getNode(loadposition),LoadValue,LoadDirection);
 
assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);     
massMatrix = assembling.assembleGlobalMassMatrix(model);

%SolveAndPlot(model)

solver = SimpleHarmonicSolvingStrategy(model,p.OMEGA);
x = solver.solve();
   
step = 1;


DisplacementDazel(:,1) = model.getDofArray.getValue(step);
DisplacementDazel(:,2) = real(model.getDofArray.getValue(step)); 
DisplacementDazel(:,3) = imag(model.getDofArray.getValue(step));
DisplacementDazel(:,4) = sqrt((DisplacementDazel(:,2)).^2 + (DisplacementDazel(:,3)).^2);
DisplacementDazel(:,5) = atan(DisplacementDazel(:,3)./DisplacementDazel(:,2));


u=0;
for i=1:4:size(DisplacementDazel,1)
    u=u+1;
    DisplacementDazelSolid(u,1)=DisplacementDazel(i,1); %Solid_x_Komplex
    DisplacementDazelSolid(u,2)=DisplacementDazel(i,4); %Solid_x_Betrag   
    DisplacementDazelSolid(u,3)=DisplacementDazel(i,5); %Solid_x_Phasenwinkel
end
u=0;
for i=2:4:size(DisplacementDazel,1)
    u=u+1;
    DisplacementDazelSolid(u,4)=DisplacementDazel(i,1); %Solid_y_Komplex
    DisplacementDazelSolid(u,5)=DisplacementDazel(i,4); %Solid_y_Betrag   
    DisplacementDazelSolid(u,6)=DisplacementDazel(i,5); %Solid_y_Phasenwinkel
end

clear u

    

%nodalForces_dazel = solver.getNodalForces(step);
  


 v3 = Visualization(model);
 v3.setScaling(1e3);
 %v3.plotUndeformed
 v3.plotDeformed;

fprintf("Finished")


%% Homogenous Non-Porous Element

clc
clear model ans solver stiffnessMatrix massMatrix 
fprintf("Processing homogenous.. \n")

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


model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});

id = 0;

for j = 1:ny
    for i = 1:nx
        id = id+1;
        a=i + (j-1)*(nx+1);
        model.addNewElement('QuadrilateralElement2d4n',id,[a,a+1,a+1+(nx+1),a+(nx+1)]);
    end
end



POISSON_RATIO = p.LAMBDA_S/(2*(p.LAMBDA_S+p.MUE_S));
K = p.LAMBDA_S*(1+POISSON_RATIO)/(3*POISSON_RATIO);
E = 3*K*(1-2*POISSON_RATIO);

model.getAllElements.setPropertyValue('LAMBDA',p.LAMBDA_S);
model.getAllElements.setPropertyValue('MU',p.MUE_S);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',p.NUMBER_GAUSS_POINT);
model.getAllElements.setPropertyValue('DENSITY',p.DENSITY_S);

% for i = 1:(ny+1)
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_X')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_Y')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_FLUID_X')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_FLUID_Y')
% end


for i = 1:1+nx:numnodes
    model.getNode(i).fixDof('DISPLACEMENT_X')
    model.getNode(i).fixDof('DISPLACEMENT_Y')

end



addPointLoad(model.getNode(loadposition),LoadValue,LoadDirection);
 
assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);     
massMatrix = assembling.assembleGlobalMassMatrix(model);


solver = SimpleHarmonicSolvingStrategy(model,p.OMEGA);
x = solver.solve();
   
step = 1;
    

DisplacementHomogenous(:,1) = model.getDofArray.getValue(step);
DisplacementHomogenous(:,2) = real(model.getDofArray.getValue(step)); 
DisplacementHomogenous(:,3) = imag(model.getDofArray.getValue(step));
DisplacementHomogenous(:,4) = sqrt((DisplacementHomogenous(:,2)).^2 + (DisplacementHomogenous(:,3)).^2);
DisplacementHomogenous(:,5) = atan(DisplacementHomogenous(:,3)./DisplacementHomogenous(:,2));


u=0;
for i=1:2:size(DisplacementHomogenous,1)
    u=u+1;
    DisplacementHomogenousSolid(u,1)=DisplacementHomogenous(i,1); %Solid_x_Komplex
    DisplacementHomogenousSolid(u,2)=DisplacementHomogenous(i,4); %Solid_x_Betrag   
    DisplacementHomogenousSolid(u,3)=DisplacementHomogenous(i,5); %Solid_x_Phasenwinkel
end
u=0;

for i=2:2:size(DisplacementHomogenous,1)
    u=u+1;
    DisplacementHomogenousSolid(u,4)=DisplacementHomogenous(i,1); %Solid_y_Komplex
    DisplacementHomogenousSolid(u,5)=DisplacementHomogenous(i,4); %Solid_y_Betrag   
    DisplacementHomogenousSolid(u,6)=DisplacementHomogenous(i,5); %Solid_y_Phasenwinkel
end

clear u;

%nodalForces_homogenous = (solver.getNodalForces(step));
 

%AllardPlot = figure('Name','AllardPlot');

v = Visualization(model);
v.setScaling(1);
%v.plotUndeformed
v.plotDeformed;

fprintf("Finished")

%% EVALUATION

 DisplacementDifferencePorous(:,1) = DisplacementAllardSolid(:,5) - DisplacementAtallaSolid(:,5);
 DisplacementDifferencePorous(:,2) = DisplacementAllardSolid(:,5) - DisplacementDazelSolid(:,5);
 DisplacementDifferencePorous(:,3) = DisplacementDazelSolid(:,5) - DisplacementAtallaSolid(:,5);
 DisplacementPorous(:,1)= DisplacementAllardSolid(:,5);
 DisplacementPorous(:,2)= DisplacementAtallaSolid(:,5);
 DisplacementPorous(:,3)= DisplacementDazelSolid(:,5);
 DisplacementPorousSum = sum(DisplacementPorous);
 % Plot Displacement Differences

DisplacementPorousEndNode_y(:,1)=DisplacementAllardSolid([1:nx+1],4);
DisplacementPorousEndNode_y(:,2)=DisplacementAtallaSolid([1:nx+1],4);
DisplacementPorousEndNode_y(:,3)=DisplacementDazelSolid([1:nx+1],4);

DisplacementDifferenceNonPorous(:,1) = DisplacementAllardSolid(:,1) - DisplacementHomogenousSolid(:,1);
DisplacementDifferenceNonPorous(:,2) = DisplacementAtallaSolid(:,1) - DisplacementHomogenousSolid(:,1);
DisplacementDifferenceNonPorous(:,3) = DisplacementDazelSolid(:,1) - DisplacementHomogenousSolid(:,1);

%DisplacementDifferenceNonPorousEndNode(1,1) = DisplacementPorousEndNode(1,1) - DisplacementHomogenousSolid(nx,1);
%DisplacementDifferenceNonPorousEndNode(1,2) = DisplacementPorousEndNode(1,2) - DisplacementHomogenousSolid(nx,1);
%DisplacementDifferenceNonPorousEndNode(1,3) = DisplacementPorousEndNode(1,3) - DisplacementHomogenousSolid(nx,1);

%DisplacementDifferenceNonPorousSum = sum(DisplacementDifferenceNonPorous);
%DisplacementDifferenceNonPorousEndNodePercent = DisplacementDifferenceNonPorousEndNode./DisplacementPorousEndNode;

% MeanDisplacementAllard_x = mean(DisplacementAllardSolid(:,2));  %  Die Mean-Displacements sind 
% MeanDisplacementAtalla_x = mean(DisplacementAtallaSolid(:,2));  %  sich sehr ähnlich!
% MeanDisplacementDazel_x = mean(DisplacementDazelSolid(:,2));
% 
% MeanDisplacementAllard_y = mean(DisplacementAllardSolid(:,5)); 
% MeanDisplacementAtalla_y = mean(DisplacementAtallaSolid(:,5)); 
% MeanDisplacementDazel_y = mean(DisplacementDazelSolid(:,5));
% 
% VarDisplacementAllard_x = var(DisplacementAllardSolid(:,2));  
% VarDisplacementAtalla_x = var(DisplacementAtallaSolid(:,2)); 
% VarDisplacementDazel_x = var(DisplacementDazelSolid(:,2));
% VarDisplacementAllard_y = var(DisplacementAllardSolid(:,5));  
% VarDisplacementAtalla_y = var(DisplacementAtallaSolid(:,5));  
% VarDisplacementDazel_y = var(DisplacementDazelSolid(:,5));
% 
% SumDisplacementAllard_x = sum(DisplacementAllardSolid(:,2));
% SumDisplacementAtalla_x = sum(DisplacementAtallaSolid(:,2)); 
% SumDisplacementDazel_x = sum(DisplacementDazelSolid(:,2));
% SumDisplacementAllard_y = sum(DisplacementAllardSolid(:,5)); 
% SumDisplacementAtalla_y = sum(DisplacementAtallaSolid(:,5));  
% SumDisplacementDazel_y = sum(DisplacementDazelSolid(:,5));
% 
% SumDifferences = sum(DisplacementDifferencePorous)

 figure1 = figure("Name","Solid x-Displacements");
 hold on
 plot([1:1:size(DisplacementAllardSolid(:,2),1)],DisplacementAllardSolid(:,2),'xk')
 plot([1:1:size(DisplacementAtallaSolid(:,2),1)],DisplacementAtallaSolid(:,2),'xb')
 plot([1:1:size(DisplacementDazelSolid(:,2),1)],DisplacementDazelSolid(:,2),'xr')

 figure2 = figure("Name","Solid y-Displacements");
 hold on
 plot([1:1:size(DisplacementAllardSolid(:,5),1)],DisplacementAllardSolid(:,5),'xk')
 plot([1:1:size(DisplacementAtallaSolid(:,5),1)],DisplacementAtallaSolid(:,5),'xb')
 plot([1:1:size(DisplacementDazelSolid(:,5),1)],DisplacementDazelSolid(:,5),'xr')

figure3 = figure("Name","Solid x-Angle");
hold on
plot([1:1:size(DisplacementAllardSolid(:,3),1)],DisplacementAllardSolid(:,3),'xk')
plot([1:1:size(DisplacementAtallaSolid(:,3),1)],DisplacementAtallaSolid(:,3),'xb')
plot([1:1:size(DisplacementDazelSolid(:,3),1)],DisplacementDazelSolid(:,3),'xr')

figure4 = figure("Name","Solid y-Angle");
hold on
plot([1:1:size(DisplacementAllardSolid(:,6),1)],DisplacementAllardSolid(:,6),'xk')
plot([1:1:size(DisplacementAtallaSolid(:,6),1)],DisplacementAtallaSolid(:,6),'xb')
plot([1:1:size(DisplacementDazelSolid(:,6),1)],DisplacementDazelSolid(:,6),'xr')

figure3 = figure("Name","Solid x-Imaginaer");
hold on
plot([1:1:size(DisplacementAllardSolid(:,3),1)],imag(DisplacementAllardSolid(:,1)),'xk')
plot([1:1:size(DisplacementAtallaSolid(:,3),1)],imag(DisplacementAtallaSolid(:,1)),'xb')
plot([1:1:size(DisplacementDazelSolid(:,3),1)],imag(DisplacementDazelSolid(:,1)),'xr')

figure4 = figure("Name","Solid y-Imaginaer");
hold on
plot([1:1:size(DisplacementAllardSolid(:,6),1)],imag(DisplacementAllardSolid(:,3)),'xk')
plot([1:1:size(DisplacementAtallaSolid(:,6),1)],imag(DisplacementAtallaSolid(:,3)),'xb')
plot([1:1:size(DisplacementDazelSolid(:,6),1)],imag(DisplacementDazelSolid(:,3)),'xr')


figure5 = figure("Name","Solid y-Real-Bottom");
hold on
plot([1:1:nx],imag(DisplacementAllardSolid([1:nx],5)),'xk')
plot([1:1:nx],imag(DisplacementAtallaSolid([1:nx],5)),'xb')
plot([1:1:nx],imag(DisplacementDazelSolid([1:nx],5)),'xr')

%end
