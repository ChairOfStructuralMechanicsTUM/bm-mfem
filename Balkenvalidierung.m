    
clear;
 
%% BEAM ONLY
model = FemModel();
model.addNewNode(1,0,0,0);
model.addNewNode(2,0.25,0,0);
model.addNewNode(3,0.5,0,0);
model.addNewNode(4,0.75,0,0);
model.addNewNode(5,1,0,0);
model.getAllNodes.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y', 'DISPLACEMENT_Z', ...
    'ROTATION_X', 'ROTATION_Y', 'ROTATION_Z'});
model.addNewElement('BeamElement3d2n',1,[1 2]);
model.addNewElement('BeamElement3d2n',2,[2 3]);
model.addNewElement('BeamElement3d2n',3,[3 4]);
model.addNewElement('BeamElement3d2n',4,[4 5]);

model.getAllElements.setPropertyValue('IY',8.33e-10);
model.getAllElements.setPropertyValue('IZ',8.33e-10);
model.getAllElements.setPropertyValue('IT',8.33e-10);
model.getAllElements.setPropertyValue('YOUNGS_MODULUS',210000000000);
model.getAllElements.setPropertyValue('POISSON_RATIO',0.3);
model.getAllElements.setPropertyValue('CROSS_SECTION',1e-4);
model.getAllElements.setPropertyValue('DENSITY',7800);

%% ADDING SPRING
% model.addNewNode(6,0.5,0.1,0);
% model.getNode(6).addDof(["DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z"]);
% model.getNode(3).addDof(["DISPLACEMENT_Z"]);
% model.addNewElement('SpringDamperElement3d2n',5,[3 6]);
% model.addNewElement('ConcentratedMassElement3d1n',6, 6);
% model.getElement(5).setPropertyValue('ELEMENTAL_STIFFNESS',8.6359e+8);
% model.getElement(5).setPropertyValue('ELEMENTAL_DAMPING',0);
% model.getElement(6).setPropertyValue('ELEMENTAL_MASS',7e0);
% model.getElement(6).setPropertyValue('VOLUME_ACCELERATION',10);
% 
% model.getNode(6).fixDof('DISPLACEMENT_X');
% model.getNode(6).fixDof('DISPLACEMENT_Z');
% model.getNode(1).fixDof('DISPLACEMENT_Z');
% model.getNode(2).fixDof('DISPLACEMENT_Z');
% model.getNode(3).fixDof('DISPLACEMENT_Z');
% model.getNode(4).fixDof('DISPLACEMENT_Z');
% model.getNode(5).fixDof('DISPLACEMENT_Z');


%% VISUALIZATION

v=Visualization(model); 
v.plotUndeformed()  
 
assembling = SimpleAssembler(model);

%% Up to here everything is done as it would be set up of Standard FEM

% Create Bloch Solver
solver = BlochInverse1D_mm(model);

% define number of phases and number of bands
numberOfPhases = 20;
numberOfBands = 3;

% call the solve function of the solver
% phases: contains the discret values of the phase
% frequencies: contains the solution, each row of represents one band
[phases,frequencies]=solver.solve(numberOfPhases,numberOfBands);

%% Visualize results
figure()

for i=1:numberOfBands
    plot(phases,frequencies(i,:),'r')
    hold on
end

title('DC, Kreis 22.5Grad - mit Feder fixiert')
xlabel('Phase k')
ylabel('frequenzy f')
xlim([0 pi])
ylim([0 100])
