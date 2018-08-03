
%%%%%%%%%%%%%%%%%%%%%%
clear
close all

io=MdpaInput('Bsp_rechteckigerEinschluss_mK.mdpa'); %specify input file
model = io.readModel(); %read the model

model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);
a=model.getAllModelParts

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