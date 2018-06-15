%% Preamble
close all;
clear all; 
% clc; 
%% Initialization

nx = 10; 
ny = 10; 

[ model, x0, xl, y0, yl ] = createRectangularPlate( 1, 1, nx, ny,'elementType', 'ShellElement3d4n');
model.getAllNodes.addDof({'DISPLACEMENT_Z', 'ROTATION_X', 'ROTATION_Y','DISPLACEMENT_X', 'DISPLACEMENT_Y', 'ROTATION_Z'})

model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 10920);
model.getAllElements.setPropertyValue('POISSON_RATIO',.3);
model.getAllElements.setPropertyValue('THICKNESS', 0.005);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',4);
model.getAllElements.setPropertyValue('DENSITY',7860);
% model.getAllElements.setPropertyValue('SHEAR_CORRECTION_FACTOR',5/6);
% model.getAllElements.addProperty('FULL_INTEGRATION',false);

middle = fix((nx+1)*(ny+1)/2)+1; 

support = [x0 xl y0 yl];
% support.fixAllDofs();
support.fixDof('DISPLACEMENT_Z');
s = EigensolverStrategy(model); 
s.solve(10);
ef=s.getEigenfrequencies


% model.getNode(middle).setDofLoad('DISPLACEMENT_Z',   -.00025);
% solver = SimpleSolvingStrategy(model);
% solver.solve();
% v=Visualization(model);
% v.setScaling(50);
% v.plotUndeformed()
% v.plotDeformed()

model.getNode(middle).getDofValue('DISPLACEMENT_Z')
