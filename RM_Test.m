% TO DO:
%       - Implement the ReissnerMindlin Test
%       - Include the possibility to choose between Reduced Integration and
%         Full Integration 
%% 
close all;
clear all; 
clc; 
%% Initialization

number_Elements = 10;
length_x = 1;

fprintf("%i x %i Elements\n", number_Elements,number_Elements);

a=linspace(0,length_x,number_Elements+1);
ii = 1;
for i=1:length(a)
    for j=1:length(a)
        node(ii) = Node(ii,a(j),a(i),0);
        ii = ii + 1; 
    end
end

ii = 1; 
for i = 1 : length(a) : (length(a)-1)^2
    for j= i:(i+length(a)-2)
        ele(ii) = ReissnerMindlinElement3d4n(ii, [node(j) node(j+1) node(j+length(a)+1) node(j+length(a))]);
        ii = ii + 1;
    end
end

nodeArray = node(:)';
nodeArray.addDof({'DISPLACEMENT_Z', 'ROTATION_X', 'ROTATION_Y'});

ii = 1; 
for i=1:length(node)
    if node(i).getX == 0 || node(i).getY == 0 || node(i).getX == length_x || node(i).getY == length_x
        boundary(ii) = node(i);
        ii= ii+1;
    end
end

boundary.fixDof('DISPLACEMENT_Z');
boundary.fixDof('ROTATION_X');
boundary.fixDof('ROTATION_Y');

elementArray = ele(:)';

elementArray.setPropertyValue('THICKNESS', .001);
elementArray.setPropertyValue('YOUNGS_MODULUS', 10920);
elementArray.setPropertyValue('POISSON_RATIO', 0.3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT', 4);
elementArray.setPropertyValue('DENSITY', 1);
elementArray.setPropertyValue('SHEAR_CORRECTION_FACTOR', 5/6);
%% Solving system
model = FemModel(nodeArray,elementArray);
% model.getAllNodes.setDofLoad('DISPLACEMENT_Z', -0.25);
model.getNode(61).setDofLoad('DISPLACEMENT_Z', -0.25);

solver = SimpleSolvingStrategy(model);
x = solver.solve();

node(61).getDofValue('DISPLACEMENT_Z')
%% Eigenfrequencies
% eigensolver = EigensolverStrategy(model);
% eigensolver.solve(10);
% eigenfrequencies = eigensolver.getEigenfrequencies();
% fprintf("Eigenfrequencies: %f\n" , eigenfrequencies);
%% Plot 
% vis = Visualization(model); 
% vis.plotUndeformed();
% vis.plotDeformed();



