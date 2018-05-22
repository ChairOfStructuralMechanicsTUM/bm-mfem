%% Preamble
close all;
clear all; 
clc; 
%% Initialization

a = 1; 
b = 1; 

node(1) = Node(1, 0, 0); 
node(2) = Node(2, a, 0); 
node(3) = Node(3, a, b); 
node(4) = Node(4, 0, b); 

nodeArray = node(:)';
% nodeArray.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z', 'ROTATION_X',  'ROTATION_Y', 'ROTATION_Z'});
nodeArray.addDof({'DISPLACEMENT_Z', 'ROTATION_X',  'ROTATION_Y'});
ele(1) = ShellElement3d4n(1, [node(1) node(2) node(3) node(4)]); 
elementArray = ele(:)';

elementArray.setPropertyValue('THICKNESS', 0.01);
elementArray.setPropertyValue('YOUNGS_MODULUS', 10920);
elementArray.setPropertyValue('POISSON_RATIO', 0.3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT', 4);
%% Solving system

model = FemModel(nodeArray,elementArray);
assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);


