clear
close all

% Creating a new model and adding nodes
model = FemModel();

model.addNewNode(1,0,0,0);
model.addNewNode(2,1,0,0);
model.addNewNode(3,1,1,0);
model.addNewNode(4,0,1,0);
model.addNewNode(5,1,2,0);
model.addNewNode(6,0,2,0);

% Creating different types of elements
model.addNewElement('QuadrilateralElement2d4n',1,[1 2 3 4]);
model.addNewElement('ClassicalPorousElement2d4n',2,[4 3 5 6]);

% Defining model parts
model.addNewModelPart('Quadri',[1 2 3 4],1);
model.addNewModelPart('Porous',[4 3 5 6],2);

% Setting property values and DOF for model parts
mp1 = model.getModelPart('Quadri');
mp1.getElements().setPropertyValue('YOUNGS_MODULUS',96);
mp1.getElements().setPropertyValue('POISSON_RATIO',1/3);
mp1.getElements().setPropertyValue('NUMBER_GAUSS_POINT',2);
mp1.getElements().setPropertyValue('DENSITY',7860);

mp1.getNodes().addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);

mp2 = model.getModelPart('Porous');
mp2.getElements().setPropertyValue('DENSITY_SOLID',30);
mp2.getElements().setPropertyValue('LAMBDA_SOLID',905357);
mp2.getElements().setPropertyValue('MUE_SOLID',264062);
mp2.getElements().setPropertyValue('DAMPING_SOLID',0);

mp2.getElements().setPropertyValue('DENSITY_FLUID',1.21);
mp2.getElements().setPropertyValue('VISCOSITY_FLUID',1.84e-5);
mp2.getElements().setPropertyValue('STANDARD_PRESSURE_FLUID',101);
mp2.getElements().setPropertyValue('HEAT_CAPACITY_FLUID',1.4);
mp2.getElements().setPropertyValue('PRANDTL_NUMBER_FLUID',0.71);

mp2.getElements().setPropertyValue('POROSITY',0.96);
mp2.getElements().setPropertyValue('TORTUOSITY',1.7);
mp2.getElements().setPropertyValue('FLOW_RESISTIVITY',32e3);
mp2.getElements().setPropertyValue('VISCOUS_LENGHT',90);
mp2.getElements().setPropertyValue('THERMAL_LENGTH',165);

mp2.getElements().setPropertyValue('FREQUENCY',510);
mp2.getElements().setPropertyValue('NUMBER_GAUSS_POINT',2);

mp2.getNodes().addDof(["DISPLACEMENT_SOLID_X", "DISPLACEMENT_SOLID_Y", "DISPLACEMENT_FLUID_X", "DISPLACEMENT_FLUID_Y"]);

% Determining intersection (common nodes)
ids1 = mp1.getNodes().getId();
ids2 = model.getModelPart('Porous').getNodes().getId();
kopplung = intersect(ids1,ids2);

% Setting BC
model.getNode(1).fixAllDofs();
model.getNode(4).fixAllDofs();
model.getNode(6).fixAllDofs();

model.initialize();

% Assembling matrices; CAUTION: Nodes at intesection might own additional
% DOFS which produce zero entries (lines and colums) in the respective
% stiffness matrices of the model part
[k1, k1r] = SimpleAssembler.assembleGlobalStiffnessMatrix(model,'Quadri');
[k2, k2r] = SimpleAssembler.assembleGlobalStiffnessMatrix(model,'Porous');

[m1, m1r] = SimpleAssembler.assembleGlobalMassMatrix(model,'Quadri');
[m2, m2r] = SimpleAssembler.assembleGlobalMassMatrix(model,'Porous');