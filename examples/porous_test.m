clear;

node01 = Node(1,-1,-1,-1);
node02 = Node(2,1,-1,-1);
node03 = Node(3,1,1,-1);
node04 = Node(4,-1,1,-1);
node05 = Node(5,-1,-1,1);
node06 = Node(6,1,-1,1);
node07 = Node(7,1,1,1);
node08 = Node(8,-1,1,1);
node09 = Node(9,-1,-1,2);
node10 = Node(10,1,-1,2);
node11 = Node(11,1,1,2);
node12 = Node(12,-1,1,2);



nodeArray = [node01 node02 node03 node04 node05 node06 node07 node08 node09 node10 node11 node12];

nodeArray.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y', 'DISPLACEMENT_SOLID_Z','DISPLACEMENT_FLUID_X', 'DISPLACEMENT_FLUID_Y', 'DISPLACEMENT_FLUID_Z'});

ele01 = PorousElement3d8n(1,[node01 node02 node03 node04 node05 node06 node07 node08]);
ele02 = PorousElement3d8n(2,[node05 node06 node07 node08 node09 node10 node11 node12]);

elementArray = [ele01 ele02];

elementArray.setPropertyValue('DENSITY_SOLID',9);
elementArray.setPropertyValue('LAMBDA_SOLID',9);
elementArray.setPropertyValue('MUE_SOLID',9);
elementArray.setPropertyValue('DAMPING_SOLID',9);

elementArray.setPropertyValue('DENSITY_FLUID',9);
elementArray.setPropertyValue('VISCOSITY_FLUID',9);
elementArray.setPropertyValue('STANDARD_PRESSURE_FLUID',9);
elementArray.setPropertyValue('HEAT_CAPACITY_FLUID',9);
elementArray.setPropertyValue('PRANDTL_NUMBER_FLUID',9);

elementArray.setPropertyValue('POROSITY',9);
elementArray.setPropertyValue('TORTUOSITY',9);
elementArray.setPropertyValue('FLOW_RESISTIVITY',9);
elementArray.setPropertyValue('VISCOUS_LENGHT',9);
elementArray.setPropertyValue('THERMAL_LENGTH',9);

elementArray.setPropertyValue('FREQUENCY',9);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);

% node01.fixDof('DISPLACEMENT_SOLID_X');
% node01.fixDof('DISPLACEMENT_SOLID_Y');
% node01.fixDof('DISPLACEMENT_SOLID_Z');
% node04.fixDof('DISPLACEMENT_SOLID_X');
% node04.fixDof('DISPLACEMENT_SOLID_Y');
% node04.fixDof('DISPLACEMENT_SOLID_Z');
% node05.fixDof('DISPLACEMENT_SOLID_X');
% node05.fixDof('DISPLACEMENT_SOLID_Y');
% node05.fixDof('DISPLACEMENT_SOLID_Z');
% node08.fixDof('DISPLACEMENT_SOLID_X');
% node08.fixDof('DISPLACEMENT_SOLID_Y');
% node08.fixDof('DISPLACEMENT_SOLID_Z');
% 
% node01.fixDof('DISPLACEMENT_FLUID_X');
% node01.fixDof('DISPLACEMENT_FLUID_Y');
% node01.fixDof('DISPLACEMENT_FLUID_Z');
% node04.fixDof('DISPLACEMENT_FLUID_X');
% node04.fixDof('DISPLACEMENT_FLUID_Y');
% node04.fixDof('DISPLACEMENT_FLUID_Z');
% node05.fixDof('DISPLACEMENT_FLUID_X');
% node05.fixDof('DISPLACEMENT_FLUID_Y');
% node05.fixDof('DISPLACEMENT_FLUID_Z');
% node08.fixDof('DISPLACEMENT_FLUID_X');
% node08.fixDof('DISPLACEMENT_FLUID_Y');
% node08.fixDof('DISPLACEMENT_FLUID_Z');

% addPointLoadPorous(node03,1,[0 -1]);
% addPointLoadPorous(node09,1,[1 0]);

model = FemModel(nodeArray, elementArray);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
            
massMatrix = assembling.assembleGlobalMassMatrix(model);

solver = SimpleSolvingStrategy(model);
x = solver.solve();

step = 1;

VerschiebungDofs = model.getDofArray.getValue(step);

nodalForces = solver.getNodalForces(step);

v = Visualization(model);
v.setScaling(1);
v.plotUndeformed
v.plotDeformed

