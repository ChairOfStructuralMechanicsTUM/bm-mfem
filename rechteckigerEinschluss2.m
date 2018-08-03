clear
close all

io=MdpaInput('Bsp_rechteckigerEinschluss_mK_3D.mdpa'); %specify input file
model = io.readModel(); %read the model

model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);
a=model.getAllModelParts

leftNodes = model.getModelPart('GENERIC_leftNodes').getNodes()
leftNodes2 = model.getModelPart('GENERIC_leftNodes2').getNodes()
rightNodes = model.getModelPart('GENERIC_rightNodes').getNodes()
rightNodes2 = model.getModelPart('GENERIC_rightNodes2').getNodes()
innerNodes = model.getModelPart('GENERIC_innerNodes').getNodes()
innerNodes2 = model.getModelPart('GENERIC_innerNodes2').getNodes()

% % %set boundary conditions for all elements in 'left_support' and 'right_support':
% % model.getModelPart('GENERIC_left_support').getNodes().fixAllDofs();
% % model.getModelPart('GENERIC_right_support').getNodes().fixAllDofs();
% % 
% % %set load for all elements in 'inner_circle':
% % model.getModelPart('GENERIC_inner_circle').getNodes().setDofLoad('DISPLACEMENT_X',10000);
% % 
% % 
% % v=Visualization(model); %set up visualization
% % v.plotUndeformed()  %visualize