%% Classic
clc
clear;


node01 = Node(1,0,0);
node02 = Node(2,2,0);
node03 = Node(3,2,1);
node04 = Node(4,0,1);

nodeArray = [node01 node02 node03 node04];

nodeArray.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y', 'DISPLACEMENT_FLUID_X', 'DISPLACEMENT_FLUID_Y'});

ele01 = ClassicalPorousElement2d4n(1,[node01 node02 node03 node04]);

elementArray = [ele01];

% assignment of material properties
elementArray.setPropertyValue('DENSITY_SOLID',30);
elementArray.setPropertyValue('LAMBDA_SOLID',905357);
elementArray.setPropertyValue('MUE_SOLID',264062);
elementArray.setPropertyValue('DAMPING_SOLID',0);

elementArray.setPropertyValue('DENSITY_FLUID',1.21);
elementArray.setPropertyValue('VISCOSITY_FLUID',1.84e-5);
elementArray.setPropertyValue('STANDARD_PRESSURE_FLUID',101);
elementArray.setPropertyValue('HEAT_CAPACITY_FLUID',1.4);
elementArray.setPropertyValue('PRANDTL_NUMBER_FLUID',0.71);

elementArray.setPropertyValue('POROSITY',0.96);
elementArray.setPropertyValue('TORTUOSITY',1.7);
elementArray.setPropertyValue('FLOW_RESISTIVITY',32e3);
elementArray.setPropertyValue('VISCOUS_LENGHT',90);
elementArray.setPropertyValue('THERMAL_LENGTH',165);

elementArray.setPropertyValue('FREQUENCY',10);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);

node01.fixDof('DISPLACEMENT_SOLID_X');
node01.fixDof('DISPLACEMENT_SOLID_Y');
node04.fixDof('DISPLACEMENT_SOLID_X');
node04.fixDof('DISPLACEMENT_SOLID_Y');

node01.fixDof('DISPLACEMENT_FLUID_X');
node01.fixDof('DISPLACEMENT_FLUID_Y');
node04.fixDof('DISPLACEMENT_FLUID_X');
node04.fixDof('DISPLACEMENT_FLUID_Y');


addPointLoadPorous(node03,1,[0 -1]);

model = FemModel(nodeArray, elementArray);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);          
massMatrix = assembling.assembleGlobalMassMatrix(model);

solver = SimpleHarmonicSolvingStrategy(model,100);
x = solver.solve();

step = 1;

VerschiebungDofs = model.getDofArray.getValue(step);

nodalForces = solver.getNodalForces(step);

v = Visualization(model);
v.setScaling(1e3);
v.plotUndeformed
v.plotDeformed
    
    
%% Mixed
clc
clear;

node01 = Node(1,0,0);
node02 = Node(2,2,0);
node03 = Node(3,2,1);
node04 = Node(4,0,1);

nodeArray = [node01 node02 node03 node04];

nodeArray.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y', 'PORE_PRESSURE'});

ele01 = MixedPorousElement2d4n(1,[node01 node02 node03 node04]);

elementArray = [ele01];

% assignment of material properties
elementArray.setPropertyValue('DENSITY_SOLID',30);
elementArray.setPropertyValue('LAMBDA_SOLID',905357);
elementArray.setPropertyValue('MUE_SOLID',264062);
elementArray.setPropertyValue('DAMPING_SOLID',0);

elementArray.setPropertyValue('DENSITY_FLUID',1.21);
elementArray.setPropertyValue('VISCOSITY_FLUID',1.84e-5);
elementArray.setPropertyValue('STANDARD_PRESSURE_FLUID',101);
elementArray.setPropertyValue('HEAT_CAPACITY_FLUID',1.4);
elementArray.setPropertyValue('PRANDTL_NUMBER_FLUID',0.71);

elementArray.setPropertyValue('POROSITY',0.96);
elementArray.setPropertyValue('TORTUOSITY',1.7);
elementArray.setPropertyValue('FLOW_RESISTIVITY',32e3);
elementArray.setPropertyValue('VISCOUS_LENGHT',90);
elementArray.setPropertyValue('THERMAL_LENGTH',165);

elementArray.setPropertyValue('FREQUENCY',100);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);

node01.fixDof('DISPLACEMENT_SOLID_X');
node01.fixDof('DISPLACEMENT_SOLID_Y');
node04.fixDof('DISPLACEMENT_SOLID_X');
node04.fixDof('DISPLACEMENT_SOLID_Y');
node02.fixDof('PORE_PRESSURE');
node03.fixDof('PORE_PRESSURE');

addPointLoadPorous(node03,1,[0 -1]);

model = FemModel(nodeArray, elementArray);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);          
massMatrix = assembling.assembleGlobalMassMatrix(model);

solver = SimpleHarmonicSolvingStrategy(model,100);
x = solver.solve();

step = 1;

VerschiebungDofs = model.getDofArray.getValue(step);

nodalForces = solver.getNodalForces(step);

v = Visualization(model);
v.setScaling(1e3);
v.plotUndeformed
v.plotDeformed
    
%% Total
clc
clear;

node01 = Node(1,0,0);
node02 = Node(2,2,0);
node03 = Node(3,2,1);
node04 = Node(4,0,1);

nodeArray = [node01 node02 node03 node04];

nodeArray.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y', 'DISPLACEMENT_TOTAL_X', 'DISPLACEMENT_TOTAL_Y'});

ele01 = TotalPorousElement2d4n(1,[node01 node02 node03 node04]);

elementArray = [ele01];

% assignment of material properties
elementArray.setPropertyValue('DENSITY_SOLID',30);
elementArray.setPropertyValue('LAMBDA_SOLID',905357);
elementArray.setPropertyValue('MUE_SOLID',264062);
elementArray.setPropertyValue('DAMPING_SOLID',0);

elementArray.setPropertyValue('DENSITY_FLUID',1.21);
elementArray.setPropertyValue('VISCOSITY_FLUID',1.84e-5);
elementArray.setPropertyValue('STANDARD_PRESSURE_FLUID',101);
elementArray.setPropertyValue('HEAT_CAPACITY_FLUID',1.4);
elementArray.setPropertyValue('PRANDTL_NUMBER_FLUID',0.71);

elementArray.setPropertyValue('POROSITY',0.96);
elementArray.setPropertyValue('TORTUOSITY',1.7);
elementArray.setPropertyValue('FLOW_RESISTIVITY',32e3);
elementArray.setPropertyValue('VISCOUS_LENGHT',90);
elementArray.setPropertyValue('THERMAL_LENGTH',165);

elementArray.setPropertyValue('FREQUENCY',100);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);

node01.fixDof('DISPLACEMENT_SOLID_X');
node01.fixDof('DISPLACEMENT_SOLID_Y');
node04.fixDof('DISPLACEMENT_SOLID_X');
node04.fixDof('DISPLACEMENT_SOLID_Y');

node01.fixDof('DISPLACEMENT_TOTAL_X');
node01.fixDof('DISPLACEMENT_TOTAL_Y');
node04.fixDof('DISPLACEMENT_TOTAL_X');
node04.fixDof('DISPLACEMENT_TOTAL_Y');


addPointLoadPorous(node03,1,[0 -1]);

model = FemModel(nodeArray, elementArray);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);          
massMatrix = assembling.assembleGlobalMassMatrix(model);

solver = SimpleHarmonicSolvingStrategy(model,100);
x = solver.solve();

step = 1;

VerschiebungDofs = model.getDofArray.getValue(step);

nodalForces = solver.getNodalForces(step);

v = Visualization(model);
v.setScaling(1e3);
v.plotUndeformed
v.plotDeformed    
