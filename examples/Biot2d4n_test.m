clear;

node01 = Node(1,-1,-1);
node02 = Node(2,1,-1);
node03 = Node(3,1,1);
node04 = Node(4,-1,1);


nodeArray = [node01 node02 node03 node04];

nodeArray.addDof({'DISPLACEMENT_SOLID_X', 'DISPLACEMENT_SOLID_Y','DISPLACEMENT_FLUID_X', 'DISPLACEMENT_FLUID_Y'});

ele01 = BiotElement2d4n(1,[node01 node02 node03 node04]);
% ele02 = PorousElement3d8n(2,[node05 node06 node07 node08 node09 node10 node11 node12]);

elementArray = ele01;

% Solid Parameters:
elementArray.setPropertyValue('DENSITY_S',30);
elementArray.setPropertyValue('LAMBDA',905357);
elementArray.setPropertyValue('MU',264062);
elementArray.setPropertyValue('ETA_S',0);
% Fluid Parameters:
elementArray.setPropertyValue('DENSITY_F',1.21);
elementArray.setPropertyValue('ETA_F',1.84e-5);
elementArray.setPropertyValue('PRESSURE_0',101);
elementArray.setPropertyValue('HEAT_CAPACITY_RATIO',1.4);
elementArray.setPropertyValue('PRANDL_NUMBER',0.71);
% Porous Parameters:
elementArray.setPropertyValue('POROSITY',0.96);
elementArray.setPropertyValue('TORTUOSITY',1.7);
elementArray.setPropertyValue('STATIC_FLOW_RESISTIVITY',32e3);
elementArray.setPropertyValue('VISCOUS_CHARACT_LENGTH',90);
elementArray.setPropertyValue('THERMAL_CHARACT_LENGTH',165);
%General Parameters:
elementArray.setPropertyValue('OMEGA',9);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT',2);

node01.fixDof('DISPLACEMENT_SOLID_X');
node01.fixDof('DISPLACEMENT_SOLID_Y');
node04.fixDof('DISPLACEMENT_SOLID_X');
node04.fixDof('DISPLACEMENT_SOLID_Y');


node01.fixDof('DISPLACEMENT_FLUID_X');
node01.fixDof('DISPLACEMENT_FLUID_Y');
node04.fixDof('DISPLACEMENT_FLUID_X');
node04.fixDof('DISPLACEMENT_FLUID_Y');


addPointLoadPorous(node03,10000,[0 -1]);


model = FemModel(nodeArray, elementArray);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.assembleGlobalStiffnessMatrix(model);
            
massMatrix = assembling.assembleGlobalMassMatrix(model);

solver = NewmarkSolvingStrategy(model);
x = solver.solve();

step = 1;

VerschiebungDofs = model.getDofArray.getValue(step);

nodalForces = solver.getNodalForces(step);

v = Visualization(model);
v.setScaling(1);
v.plotUndeformed
v.plotDeformed

