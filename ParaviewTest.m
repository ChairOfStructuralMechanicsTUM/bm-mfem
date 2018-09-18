%% Preamble
close all;
clear; 
% clc; 
%% Initialization

nx = 20;
ny = 20;
[model, x0, xl, y0, yl] = createRectangularPlate(1, 1, nx, ny, 'elementType', 'ReissnerMindlinElement3d4n');
model.getAllNodes.addDof(["DISPLACEMENT_Z", "ROTATION_X", "ROTATION_Y"]);
support = [x0 xl y0 yl];
support.fixDof('DISPLACEMENT_Z');

% Properties
model.getAllElements.setPropertyValue('THICKNESS', .00125);
model.getAllElements.setPropertyValue('YOUNGS_MODULUS', 2.07e11);
model.getAllElements.setPropertyValue('POISSON_RATIO', 0.3);
model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT', 4);
model.getAllElements.setPropertyValue('DENSITY', 7850);
model.getAllElements.setPropertyValue('SHEAR_CORRECTION_FACTOR', 5/6);


middle = fix((nx+1)*(ny+1)/2)+1; 

% Solver
model.getNode(middle).setDofLoad('DISPLACEMENT_Z', 50);

dt = .01;
time = 0;
endTime = .2;
solver = NewmarkSolvingStrategy(model, dt);

while time < endTime
    solver.solve();
    time = time + dt;               
end
 
% v=Visualization(model);
% v.setScaling(50);
% v.plotUndeformed()
% v.plotDeformed()

tic
v =  VisualizationParaviewXML(model, 'VTU_test', 'DISPLACEMENT');
v.pvdWrite();
toc

actualDisplacementZ = model.getNode(middle).getDofValue('DISPLACEMENT_Z',1);

