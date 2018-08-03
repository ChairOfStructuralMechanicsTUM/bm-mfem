clear;

node01 = Node(1,0,0);
node02 = Node(2,2,0);
node03 = Node(3,2,1);
node04 = Node(4,0,1);

nodeArray = [node01 node02 node03 node04];
nodeArray.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y','PRESSURE_FLUID'});

<<<<<<< HEAD:examples/MixedFormulation2d4n_test.m
ele01 = MixedDisplacementElement2d4n(1,[node01 node02 node03 node04]);
=======


ele01 = MixedDisplacementElement2d4n(1,[node01 node02 node03 node04]);

%ele01 = MixedFormulationElement2d4n(1,[node01 node02 node03 node04]);


>>>>>>> 7883e3470bdd45c2d8ea706f0ee3bf7cb639b5c7:examples/MixedDisplacement2d4n_test.m
elementArray = ele01;

% Solid Parameters:
elementArray.setPropertyValue('DENSITY_S',30);
elementArray.setPropertyValue('LAMBDA_S',905357);
elementArray.setPropertyValue('MUE_S',264062);
elementArray.setPropertyValue('ETA_S',0);
% Fluid Parameters:
elementArray.setPropertyValue('DENSITY_F',1.21);
elementArray.setPropertyValue('ETA_F',1.84e-5);
elementArray.setPropertyValue('PRESSURE_0_F',101);
elementArray.setPropertyValue('HEAT_CAPACITY_RATIO_F',1.4);
elementArray.setPropertyValue('PRANDL_NUMBER_F',0.71);
% Porous Parameters:
elementArray.setPropertyValue('POROSITY',0.96);
elementArray.setPropertyValue('TORTUOSITY',1.7);
elementArray.setPropertyValue('STATIC_FLOW_RESISTIVITY',32e3);
elementArray.setPropertyValue('VISCOUS_CHARACT_LENGTH',90);
elementArray.setPropertyValue('THERMAL_CHARACT_LENGTH',165);
% General Parameters:
elementArray.setPropertyValue('OMEGA',100);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);

% Boudary Conditions for solid phase:
node01.fixDof('DISPLACEMENT_SOLID_X');
node01.fixDof('DISPLACEMENT_SOLID_Y');
node04.fixDof('DISPLACEMENT_SOLID_X');
node04.fixDof('DISPLACEMENT_SOLID_Y');
% Boundary Condition for pressure: p = 0 at boundaries
node02.fixDof('PRESSURE_FLUID');
node03.fixDof('PRESSURE_FLUID');

addPointLoadPorous(node03,1,[0 -1]);

model = FemModel(nodeArray, elementArray);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);     
massMatrix = assembling.assembleGlobalMassMatrix(model);

<<<<<<< HEAD:examples/MixedFormulation2d4n_test.m
solver = SimpleHarmonicSolvingStrategy(model,100);
x = solver.solve();

step = 1;

VerschiebungDofs = model.getDofArray.getValue(step);
=======
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
         
massMatrix = assembling.assembleGlobalMassMatrix(model);
>>>>>>> 7883e3470bdd45c2d8ea706f0ee3bf7cb639b5c7:examples/MixedDisplacement2d4n_test.m

nodalForces = solver.getNodalForces(step);

v = Visualization(model);
v.setScaling(1e3);
v.plotUndeformed
v.plotDeformed

% Newmark-A
% v = Visualization(model);
% v.setScaling(1e3);
% v.plotUndeformed
% v.plotDeformed
% v = Visualization(model);
% v.setScaling(100);
% v.plotUndeformed
% %v.plotDeformed
% dt = .05;
% time = 0;
% endTime = 1.5;
% ls = linspace(time,endTime,endTime/dt+1);
% % disp = zeros(1,endTime/dt+1);
% solver = NewmarkSolvingStrategy(model, dt);
% for j = 1:30
%     solver.solve();
%     time = j*0.05;
%     hold on
%     v.plotDeformed();
%     view(3)
%     axis([-1.1 1.1 -1.1 1.1 -1.1 2.1])
%     F(j) = getframe(gcf);
% end
% fig = figure;
% movie(fig,F,1,2)
% step = 30;
% VerschiebungDofs = model.getDofArray.getValue(step);
% Newmark-B

