close all;
clear all; 
clc; 
%% Initialization
number_Elements = 2;
fprintf("%i x %i Elements\n", number_Elements,number_Elements);

a=linspace(0,1,number_Elements+1);

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
    if node(i).getX == 0 || node(i).getY == 0 || node(i).getX == 1 || node(i).getY == 1
        boundary(ii) = node(i);
        ii= ii+1;
    end
end

nodeArray.fixDof('DISPLACEMENT_Z');
nodeArray.fixDof('ROTATION_X');
nodeArray.fixDof('ROTATION_Y');

elementArray = ele(:)';

elementArray.setPropertyValue('THICKNESS', 0.0125);
elementArray.setPropertyValue('YOUNGS_MODULUS', 200e9);
elementArray.setPropertyValue('SHEAR_MODULUS', 200e9/2.5);
elementArray.setPropertyValue('POISSON_RATIO', 0.3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT', 4);
elementArray.setPropertyValue('DENSITY', 1);
elementArray.setPropertyValue('SHEAR_CORRECTION_FACTOR', 5/6);

addPointLoad_Plate(node(5), 100, [1,0,0]);


%% Solving the system
model = FemModel(nodeArray,elementArray);

assembling = SimpleAssembler(model);
stiffnessMatrix = assembling.stiffnessMatrix;
forceVector = assembling.forceVector
reducedForceVector = assembling.reducedForceVector;

solver = SimpleSolvingStrategy(model);
x = solver.solve();

step = 1;

%% Plot 
vis = Visualization(model); 
vis.plotUndeformed();
vis.plotDeformed();
%% EigenFrquencies

% eigensolver = EigensolverStrategy(model);
% eigensolver.solve(5);
% eigenfrequencies = eigensolver.getEigenfrequencies();
% disp(eigenfrequencies);





