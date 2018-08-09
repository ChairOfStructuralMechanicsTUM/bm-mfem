% Dieser Test verknüpft ein homogenes Element vom Typ QuadrilateralElement2d4n
% mit einem porösen Element vom Typ Porous2d4n
%
% Dieser Test ist eine Verbesserung des porous_homogen_test_2 durch:
% - Auslagerung der Matrizen-Assemblierung in eine eigene Funktion
% (Assembler_Porous_Homogen)
% - Implementierung der Berechnung der Massenmatrix
% - Implementierung der Berechnung der reduzierten Matrizen infolge von
% aufgebrachten Randbedingungen
%
% Einschränkungen:
% - der Rand der beiden Elemente ist parallel zur x-Achse
% - es wird angenommen, dass die Elemente fest verbunden sind (bounded)


tic;
clear;
close all;
 
% Creating a new model 

% Adding nodes
model = FemModel();
 
model.addNewNode(1,0,0,0);
model.addNewNode(2,1,0,0);
model.addNewNode(3,1,1,0);
model.addNewNode(4,0,1,0);
model.addNewNode(5,1,2,0);
model.addNewNode(6,0,2,0);
 
% Creating different types of elements
model.addNewElement('QuadrilateralElement2d4n',1,[1 2 3 4]);
model.addNewElement('Porous2d4n',2,[4 3 5 6]);
 
% Defining model parts
model.addNewModelPart('Homogen',[1 2 3 4],1);
model.addNewModelPart('Porous',[4 3 5 6],2);
 
% Setting property values and DOF for model parts
mp1 = model.getModelPart('Homogen');
mp1.getElements().setPropertyValue('YOUNGS_MODULUS',96);
mp1.getElements().setPropertyValue('POISSON_RATIO',1/3);
mp1.getElements().setPropertyValue('NUMBER_GAUSS_POINT',2);
mp1.getElements().setPropertyValue('DENSITY',7860);
 
mp2 = model.getModelPart('Porous');
%Values
mp2.getElements.setPropertyValue('FREQUENCY',100); %angenommen 100 Hz
mp2.getElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
%Frame
mp2.getElements.setPropertyValue('FRAME_LAME_PARAMETER_LAMBDA',114400);
mp2.getElements.setPropertyValue('FRAME_LAME_PARAMETER_MU',28600);
%Fluid
mp2.getElements.setPropertyValue('HEAT_CAPACITY_RATIO',1.4);
mp2.getElements.setPropertyValue('AMBIENT_FLUID_STANDARD_PRESSURE',101300); %angenommen 1013 hPa
mp2.getElements.setPropertyValue('AMBIENT_FLUID_VISCOSITY',1.84*10^-5);
mp2.getElements.setPropertyValue('AMBIENT_FLUID_DENSITY',1.21);
mp2.getElements.setPropertyValue('PRANDTL_NUMBER',0.71);
%Porous
mp2.getElements.setPropertyValue('THERMAL_CHAR_LENGTH',226*10^-6);
mp2.getElements.setPropertyValue('POROSITY',0.9);
%Mass Matrix
mp2.getElements.setPropertyValue('FRAME_DENSITY',300);
mp2.getElements.setPropertyValue('TORTUOSITY',7.8);
mp2.getElements.setPropertyValue('VISCOUS_CHAR_LENGTH',226*10^-6);
mp2.getElements.setPropertyValue('STATIC_FLOW_RESISTIVITY',2500);

% adding dofs to elements
mp1.getNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);
mp2.getNodes.addDof(["FRAME_DISPLACEMENT_X", "FRAME_DISPLACEMENT_Y", ...
    "FLUID_DISPLACEMENT_X", "FLUID_DISPLACEMENT_Y"]);

% Setting BC
model.getNode(1).fixAllDofs();
model.getNode(4).fixAllDofs();
model.getNode(6).fixAllDofs();

[K,Kred] = Assembler_Porous_Homogen.assembleGlobalStiffnessMatrix(model);
[M,Mred] = Assembler_Porous_Homogen.assembleGlobalMassMatrix(model);

time=toc;
