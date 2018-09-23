 %% SET PARAMETERS


%function [Answer,lambda1,lambda2,lambda3] = ModelTest(POROSITY, DENSITY_S, DENSITY_F, OMEGA, NUMBER_GAUSS_POINT)
%function
clc
Lx = 0.05;
Ly = 0.01;

ny = 18;
nx = 5*ny;


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
p.OMEGA = OMEGA;
p.NUMBER_GAUSS_POINT = NUMBER_GAUSS_POINT;

WavelengthCheck(p.OMEGA,p.TORTUOSITY,p.ETA_F,p.DENSITY_F,DENSITY_S,p.STATIC_FLOW_RESISTIVITY,p.VISCOUS_CHARACT_LENGTH,p.POROSITY,p.HEAT_CAPACITY_RATIO_F,p.PRESSURE_0_F,...
    p.PRANDL_NUMBER_F,p.THERMAL_CHARACT_LENGTH,p.ETA_S,p.LAMBDA_S,p.MUE_S,Lx,Ly,nx,ny);



%% Standard (u_s,u_f) displacement
clc
clear model ans solver stiffnessMatrix massMatrix DisplacementAllard DisplacementAllardSolid
fprintf("Processing.. \n")
modelAllard = FemModel();

dx = Lx/nx;
dy = Ly/ny;

numnodes = (nx + 1) * (ny+1);
numele = nx*ny;

id = 0;

for j = 1:(ny+1)
    for i= 1:(nx+1)
        id = id+1;
        modelAllard.addNewNode(id,(i-1)*dx,(j-1)*dy);
    end
end


modelAllard.getAllNodes.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y', 'DISPLACEMENT_FLUID_X', 'DISPLACEMENT_FLUID_Y'});

id = 0;

for j = 1:ny
    for i = 1:nx
        id = id+1;
        a=i + (j-1)*(nx+1);
        modelAllard.addNewElement('BiotAllardElement2d4n',id,[a,a+1,a+1+(nx+1),a+(nx+1)]);
    end
end

modelAllard.getAllElements.setPropertyValue('DENSITY_S',p.DENSITY_S);
modelAllard.getAllElements.setPropertyValue('LAMBDA_S',p.LAMBDA_S);
modelAllard.getAllElements.setPropertyValue('MUE_S',p.MUE_S);
modelAllard.getAllElements.setPropertyValue('ETA_S',p.ETA_S);
% Fluid Parameters:
modelAllard.getAllElements.setPropertyValue('DENSITY_F',p.DENSITY_F);
modelAllard.getAllElements.setPropertyValue('ETA_F',p.ETA_F);
modelAllard.getAllElements.setPropertyValue('PRESSURE_0_F',p.PRESSURE_0_F);
modelAllard.getAllElements.setPropertyValue('HEAT_CAPACITY_RATIO_F',p.HEAT_CAPACITY_RATIO_F);
modelAllard.getAllElements.setPropertyValue('PRANDL_NUMBER_F',p.PRANDL_NUMBER_F);
% Porous Parameters:
modelAllard.getAllElements.setPropertyValue('POROSITY',p.POROSITY);
modelAllard.getAllElements.setPropertyValue('TORTUOSITY',p.TORTUOSITY);
modelAllard.getAllElements.setPropertyValue('STATIC_FLOW_RESISTIVITY',p.STATIC_FLOW_RESISTIVITY);
modelAllard.getAllElements.setPropertyValue('VISCOUS_CHARACT_LENGTH',p.VISCOUS_CHARACT_LENGTH);
modelAllard.getAllElements.setPropertyValue('THERMAL_CHARACT_LENGTH',p.THERMAL_CHARACT_LENGTH);
% General Parameters:
modelAllard.getAllElements.setPropertyValue('OMEGA',p.OMEGA);
modelAllard.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',p.NUMBER_GAUSS_POINT);

% for i = 1:(ny+1)
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_X')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_Y')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_FLUID_X')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_FLUID_Y')
% end


for i = 1:1+nx:numnodes
    modelAllard.getNode(i).fixDof('DISPLACEMENT_SOLID_X')
    modelAllard.getNode(i).fixDof('DISPLACEMENT_SOLID_Y')
    modelAllard.getNode(i).fixDof('DISPLACEMENT_FLUID_X')
    modelAllard.getNode(i).fixDof('DISPLACEMENT_FLUID_Y') 
end



addPointLoadPorous(modelAllard.getNode(numnodes),LoadValue,LoadDirection); 
assembling = SimpleAssembler(modelAllard);
fprintf('Assembling Stiffness \n')
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(modelAllard);
fprintf('Assembling Mass \n')
massMatrix = assembling.assembleGlobalMassMatrix(modelAllard);

fprintf('Solving \n')
solver = SimpleHarmonicSolvingStrategy(modelAllard,p.OMEGA);
x = solver.solve();
step = 1;
fprintf('Solved \n')
%DisplacementAllard = sparse(nx*ny*4,1);
DisplacementAllard(:,1) = modelAllard.getDofArray.getValue(step);
DisplacementAllard(:,2) = real(modelAllard.getDofArray.getValue(step)); 
DisplacementAllard(:,3) = imag(modelAllard.getDofArray.getValue(step));
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
 

 %v = Visualization(modelAllard);
 %v.setScaling(1);
 %v.plotUndeformed
 %v.plotDeformed;
 
fprintf("Finished")

%% Mixed (u_s,p) displacement
clc
clear model ans solver stiffnessMatrix massMatrix 
fprintf("Processing.. \n")
modelAtalla = FemModel();

dx = Lx/nx;
dy = Ly/ny;

numnodes = (nx + 1) * (ny+1);
numele = nx*ny;

id = 0;

for j = 1:(ny+1)
    for i= 1:(nx+1)
        id = id+1;
        modelAtalla.addNewNode(id,(i-1)*dx,(j-1)*dy);
    end
end


modelAtalla.getAllNodes.addDof({'PRESSURE_FLUID', 'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y'});

id = 0;

for j = 1:ny
    for i = 1:nx
        id = id+1;
        a=i + (j-1)*(nx+1);
        modelAtalla.addNewElement('AtallaElement2d4n',id,[a,a+1,a+1+(nx+1),a+(nx+1)]);
    end
end

modelAtalla.getAllElements.setPropertyValue('DENSITY_S',p.DENSITY_S);
modelAtalla.getAllElements.setPropertyValue('LAMBDA_S',p.LAMBDA_S);
modelAtalla.getAllElements.setPropertyValue('MUE_S',p.MUE_S);
modelAtalla.getAllElements.setPropertyValue('ETA_S',p.ETA_S);
% Fluid Parameters:
modelAtalla.getAllElements.setPropertyValue('DENSITY_F',p.DENSITY_F);
modelAtalla.getAllElements.setPropertyValue('ETA_F',p.ETA_F);
modelAtalla.getAllElements.setPropertyValue('PRESSURE_0_F',p.PRESSURE_0_F);
modelAtalla.getAllElements.setPropertyValue('HEAT_CAPACITY_RATIO_F',p.HEAT_CAPACITY_RATIO_F);
modelAtalla.getAllElements.setPropertyValue('PRANDL_NUMBER_F',p.PRANDL_NUMBER_F);
% Porous Parameters:
modelAtalla.getAllElements.setPropertyValue('POROSITY',p.POROSITY);
modelAtalla.getAllElements.setPropertyValue('TORTUOSITY',p.TORTUOSITY);
modelAtalla.getAllElements.setPropertyValue('STATIC_FLOW_RESISTIVITY',p.STATIC_FLOW_RESISTIVITY);
modelAtalla.getAllElements.setPropertyValue('VISCOUS_CHARACT_LENGTH',p.VISCOUS_CHARACT_LENGTH);
modelAtalla.getAllElements.setPropertyValue('THERMAL_CHARACT_LENGTH',p.THERMAL_CHARACT_LENGTH);
% General Parameters:
modelAtalla.getAllElements.setPropertyValue('OMEGA',p.OMEGA);
modelAtalla.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',p.NUMBER_GAUSS_POINT);


% for i = 1:(ny+1)
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_X')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_Y') 
% end

for i = 1:1+nx:numnodes
    modelAtalla.getNode(i).fixDof('DISPLACEMENT_SOLID_X')
    modelAtalla.getNode(i).fixDof('DISPLACEMENT_SOLID_Y') 
end

for j = 2:1+nx
   modelAtalla.getNode(j).fixDof('PRESSURE_FLUID')
end

for j = 1+nx:1+nx:numnodes
   modelAtalla.getNode(j).fixDof('PRESSURE_FLUID')
end

for j = numnodes-(nx-1):1:numnodes
   modelAtalla.getNode(j).fixDof('PRESSURE_FLUID')
end


addPointLoadPorous(modelAtalla.getNode(numnodes),LoadValue,LoadDirection);
 
assembling = SimpleAssembler(modelAtalla);    
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(modelAtalla);     
massMatrix = assembling.assembleGlobalMassMatrix(modelAtalla);

%SolveAndPlot(model)

solver = SimpleHarmonicSolvingStrategy(modelAtalla,p.OMEGA);
x = solver.solve();
   
step = 1;

%nodalForces_atalla = solver.getNodalForces(step);
  
%VerschiebungDofs_atalla = (model.getDofArray.getValue(step));

DisplacementAtalla(:,1) = modelAtalla.getDofArray.getValue(step);
DisplacementAtalla(:,2) = real(modelAtalla.getDofArray.getValue(step)); 
DisplacementAtalla(:,3) = imag(modelAtalla.getDofArray.getValue(step));
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
%  v2 = Visualization(modelAtalla);
%  v2.setScaling(1);
%  %v2.plotUndeformed
%  v2.plotDeformed;
% % 
% 

fprintf("Finished")


%% Total (u_s,u_t) displacement

clc
clear model ans solver stiffnessMatrix massMatrix 
close all
fprintf("Processing.. \n")

modelDazel = FemModel();

dx = Lx/nx;
dy = Ly/ny;

numnodes = (nx + 1) * (ny+1);
numele = nx*ny;

id = 0;

for j = 1:(ny+1)
    for i= 1:(nx+1)
        id = id+1;
        modelDazel.addNewNode(id,(i-1)*dx,(j-1)*dy);
    end
end


modelDazel.getAllNodes.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y', 'DISPLACEMENT_TOTAL_X', 'DISPLACEMENT_TOTAL_Y'});

id = 0;

for j = 1:ny
    for i = 1:nx
        id = id+1;
        a=i + (j-1)*(nx+1);
        modelDazel.addNewElement('DazelElement2d4n',id,[a,a+1,a+1+(nx+1),a+(nx+1)]);
    end
end

modelDazel.getAllElements.setPropertyValue('DENSITY_S',p.DENSITY_S);
modelDazel.getAllElements.setPropertyValue('LAMBDA_S',p.LAMBDA_S);
modelDazel.getAllElements.setPropertyValue('MUE_S',p.MUE_S);
modelDazel.getAllElements.setPropertyValue('ETA_S',p.ETA_S);
% Fluid Parameters:
modelDazel.getAllElements.setPropertyValue('DENSITY_F',p.DENSITY_F);
modelDazel.getAllElements.setPropertyValue('ETA_F',p.ETA_F);
modelDazel.getAllElements.setPropertyValue('PRESSURE_0_F',p.PRESSURE_0_F);
modelDazel.getAllElements.setPropertyValue('HEAT_CAPACITY_RATIO_F',p.HEAT_CAPACITY_RATIO_F);
modelDazel.getAllElements.setPropertyValue('PRANDL_NUMBER_F',p.PRANDL_NUMBER_F);
% Porous Parameters:
modelDazel.getAllElements.setPropertyValue('POROSITY',p.POROSITY);
modelDazel.getAllElements.setPropertyValue('TORTUOSITY',p.TORTUOSITY);
modelDazel.getAllElements.setPropertyValue('STATIC_FLOW_RESISTIVITY',p.STATIC_FLOW_RESISTIVITY);
modelDazel.getAllElements.setPropertyValue('VISCOUS_CHARACT_LENGTH',p.VISCOUS_CHARACT_LENGTH);
modelDazel.getAllElements.setPropertyValue('THERMAL_CHARACT_LENGTH',p.THERMAL_CHARACT_LENGTH);
% General Parameters:
modelDazel.getAllElements.setPropertyValue('OMEGA',p.OMEGA);
modelDazel.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',p.NUMBER_GAUSS_POINT);

% for i = 1:(ny+1)
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_X')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_SOLID_Y')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_TOTAL_X')
%     model.getNode(i+(i-1)*(nx+1)).fixDof('DISPLACEMENT_TOTAL_Y')
% end


for i = 1:1+nx:numnodes
    modelDazel.getNode(i).fixDof('DISPLACEMENT_SOLID_X')
    modelDazel.getNode(i).fixDof('DISPLACEMENT_SOLID_Y')
    modelDazel.getNode(i).fixDof('DISPLACEMENT_TOTAL_X')
    modelDazel.getNode(i).fixDof('DISPLACEMENT_TOTAL_Y') 
end




addPointLoadPorous(modelDazel.getNode(numnodes),LoadValue,LoadDirection);
 
assembling = SimpleAssembler(modelDazel);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(modelDazel);     
massMatrix = assembling.assembleGlobalMassMatrix(modelDazel);

%SolveAndPlot(model)

solver = SimpleHarmonicSolvingStrategy(modelDazel,p.OMEGA);
x = solver.solve();
   
step = 1;


DisplacementDazel(:,1) = modelDazel.getDofArray.getValue(step);
DisplacementDazel(:,2) = real(modelDazel.getDofArray.getValue(step)); 
DisplacementDazel(:,3) = imag(modelDazel.getDofArray.getValue(step));
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
  


%  v3 = Visualization(modelDazel);
%  v3.setScaling(1);
%  %v3.plotUndeformed
%  v3.plotDeformed;

fprintf("Finished")

%% Plots
%%Verschiebungsfiguren + Imaginärteil

%ALLARD:
%Auslesen von Knotenkoordinaten

xxAllard = modelAllard.getAllNodes.getX();

yyAllard = modelAllard.getAllNodes.getY();

%Auslesen von Knotenverschiebungen

uxAllard = modelAllard.getAllNodes.getDofValue('DISPLACEMENT_SOLID_X');

uyAllard = modelAllard.getAllNodes.getDofValue('DISPLACEMENT_SOLID_Y');

%Skalierung zur besseren Visualisierung der Ergebnisse (damit man was

%sieht)

scaling = 1;

%Berechnen der Phase bzgl. ux oder uy; 41 und 21 sind hier die Anzahl der

%Knoten in x- bzw. y-Richtung.

z = reshape(imag(uyAllard),91,19);

%z = reshape(angle(uy),41,21);

%Berechnen der Knotenkoordinaten im Verformten System

xxx = reshape(xxAllard+scaling*real(uxAllard.'),91,19);

yyy = reshape(yyAllard+scaling*real(uyAllard.'),91,19);

%Abbilden der Ergebnisse

figure()


subplot(2,1,1)

surf(xxx,yyy,z,'FaceColor','interp')

xlabel("x [m]")
ylabel("y [m]")
c = colorbar
% c.Limits = [0 4e-8] 
view(0,90)

%% ATALLA:

%Auslesen von Knotenkoordinaten

xxAtalla = modelAtalla.getAllNodes.getX();

yyAtalla = modelAtalla.getAllNodes.getY();

%Auslesen von Knotenverschiebungen

uxAtalla = modelAtalla.getAllNodes.getDofValue('DISPLACEMENT_SOLID_X');

uyAtalla = modelAtalla.getAllNodes.getDofValue('DISPLACEMENT_SOLID_Y');

%Skalierung zur besseren Visualisierung der Ergebnisse (damit man was

%sieht)

scaling = 1;

%Berechnen der Phase bzgl. ux oder uy; 41 und 21 sind hier die Anzahl der

%Knoten in x- bzw. y-Richtung.

z = reshape(imag(uyAtalla),91,19);

%z = reshape(angle(uy),41,21);

%Berechnen der Knotenkoordinaten im Verformten System

xxx = reshape(xxAtalla+scaling*real(uxAtalla.'),91,19);

yyy = reshape(yyAtalla+scaling*real(uyAtalla.'),91,19);

%Abbilden der Ergebnisse

figure()


subplot(2,1,1)

surf(xxx,yyy,z,'FaceColor','interp')

xlabel("x [m]")
ylabel("y [m]")
c = colorbar
% c.Limits = [0 4e-8] 
view(0,90)



%% Dazel:
%Auslesen von Knotenkoordinaten

xxDazel = modelDazel.getAllNodes.getX();

yyDazel = modelDazel.getAllNodes.getY();

%Auslesen von Knotenverschiebungen

uxDazel = modelDazel.getAllNodes.getDofValue('DISPLACEMENT_SOLID_X');

uyDazel = modelDazel.getAllNodes.getDofValue('DISPLACEMENT_SOLID_Y');

%Skalierung zur besseren Visualisierung der Ergebnisse (damit man was

%sieht)

scaling = 1;

%Berechnen der Phase bzgl. ux oder uy; 41 und 21 sind hier die Anzahl der

%Knoten in x- bzw. y-Richtung.

z = reshape(imag(uyDazel),91,19);

%z = reshape(angle(uy),41,21);

%Berechnen der Knotenkoordinaten im Verformten System

xxx = reshape(xxDazel+scaling*real(uxDazel.'),91,19);

yyy = reshape(yyDazel+scaling*real(uyDazel.'),91,19);

%Abbilden der Ergebnisse
figure()


subplot(2,1,1)

surf(xxx,yyy,z,'FaceColor','interp')

xlabel("x [m]")
ylabel("y [m]")
c = colorbar
% c.Limits = [0 4e-8] 
view(0,90)


%% Bottom Nodes Plot
figure1 = figure("Name","Solid y-Displacements");
hold on
plot([1:1:91],real(uyAllard([1:1:91],1)),'*k')
plot([1:1:91],real(uyAtalla([1:1:91],1)),'xb')
plot([1:1:91],real(uyDazel([1:1:91],1)),'.r')
xlabel("Knoten-Id")
ylabel("Re(u_s_y)")
legend("(u_s,u_f)-Formulierung","(u_s,p)-Formulierung","(u_s,u_t)-Formulierung")

% figure2 = figure("Name","Solid y-Displacements");
% hold on
% plot([1:1:91],imag(uyAllard([1:1:91],1)),'*k')
% % plot([1:1:91],imag(DisplacementAtallaSolid([1:1:91],1)),'xb')
% % plot([1:1:91],imag(DisplacementDazelSolid([1:1:91],1)),'.r')
% xlabel("Knoten-Id")
% ylabel("Im(u_s_y)")
% legend("(u_s,u_f)-Formulierung")

figure3 = figure("Name","Solid y-Displacements");
hold on
plot([1:1:91],imag(uyAllard([1:1:91],1)),'*k')
plot([1:1:91],imag(uyAtalla([1:1:91],1)),'xb')
plot([1:1:91],imag(uyDazel([1:1:91],1)),'.r')
xlabel("Knoten-Id")
ylabel("Im(u_s_y)")
legend("(u_s,u_f)-Formulierung","(u_s,p)-Formulierung","(u_s,u_t)-Formulierung")

figure4 = figure("Name","Solid y-Displacements");
hold on
plot([1:1:91],angle(uyAllard([1:1:91],1)),'*k')
plot([1:1:91],angle(uyAtalla([1:1:91],1)),'xb')
plot([1:1:91],angle(uyDazel([1:1:91],1)),'.r')
xlabel("Knoten-Id")
ylabel("\phi(u_s_y)")
legend("(u_s,u_f)-Formulierung","(u_s,p)-Formulierung","(u_s,u_t)-Formulierung")
















