%% 
close all;
clear all; 
% clc; 
%% Initialization

a= 1; 
b= 1; 

node(1) = Node(1, 0,0); 
node(2) = Node(2, a,0); 
node(3) = Node(3, a,b); 
node(4) = Node(4, 0,b); 

nodeArray = node(:)';
nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z', 'ROTATION_X',  'ROTATION_Y', 'ROTATION_Z'});

ele(1) = ShellElement3d4n(1, [node(1) node(2) node(3) node(4)]); 
elementArray = ele(:)';


elementArray.setPropertyValue('THICKNESS', 1);
elementArray.setPropertyValue('YOUNGS_MODULUS', 8000);
elementArray.setPropertyValue('POISSON_RATIO', 1/3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT', 4);
elementArray.setPropertyValue('DENSITY', 1);
elementArray.setPropertyValue('SHEAR_CORRECTION_FACTOR', 5/6);

%% Solving system

model = FemModel(nodeArray,elementArray);
assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);

