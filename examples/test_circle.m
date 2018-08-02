clear
close all

io=MdpaInput('test_circle.mdpa'); %specify input file
model = io.readModel(); %read the model
model.getAllNodes.addDof(["ROTATION_Y","DISPLACEMENT_Y", "DISPLACEMENT_X", "DISPLACEMENT_Z", ...
    "ROTATION_X",  "ROTATION_Z"]);
model.getAllElements.addProperty('NUMBER_GAUSS_POINT',2);


%set boundary conditions for all elements in 'left_support' and 'right_support':
model.getModelPart('GENERIC_left_support').getNodes().fixAllDofs();
model.getModelPart('GENERIC_right_support').getNodes().fixAllDofs();

%set load for all elements in 'inner_circle':
model.getModelPart('GENERIC_inner_circle').getNodes().setDofLoad('DISPLACEMENT_X',10000);


v=Visualization(model); %set up visualization
v.plotUndeformed  %visualize

