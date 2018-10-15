%% 
close all;
clear all;
clc;
%% Initialization

number_Elements = 20;
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
nodeArray = node(:)';
nodeArray.addDof({'DISPLACEMENT_X','DISPLACEMENT_Y'});

ii = 1; 
for i = 1 : length(a) : (length(a)-1)^2
    for j= i:(i+length(a)-2)
        ele(ii) = PlaneStressElement3d4n(ii, [node(j) node(j+1) node(j+length(a)+1) node(j+length(a))]);
        ii = ii + 1;
    end
end
elementArray = ele(:)';

ii = 1; 
for i=1:length(node)
    if node(i).getX == 0 %&& node(i).getY == 0 || node(i).getX == length_x && node(i).getY == 0
        node(i).fixDof('DISPLACEMENT_X');
        node(i).fixDof('DISPLACEMENT_Y');
        ii= ii+1;
    end
end



elementArray.setPropertyValue('THICKNESS', 0.2);
elementArray.setPropertyValue('YOUNGS_MODULUS', 1000);
elementArray.setPropertyValue('POISSON_RATIO', 0.3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT', 4);
elementArray.setPropertyValue('DENSITY', 1);

%% Line load
for i=1:length(node)
    if node(i).getX == 0 && node(i).getY == length_x
        IdUpperLCorner = node(i).getId;
    end
    if node(i).getX == length_x && node(i).getY == length_x
        IdUpperRCorner = node(i).getId;
    end
end
        
addConstLineLoad(IdUpperLCorner,IdUpperRCorner,nodeArray,10,[0 -1 0]);

%% Solving system
model = FemModel(nodeArray,elementArray);
% solver = SimpleSolvingStrategy(model);
% solver.solve();
% actualDisplacementY = model.getAllNodes.getDofValue('DISPLACEMENT_Y');
% 
% lower_midpoint = number_Elements/2 + 1;
% fprintf("Displacement at Node %i : %f\n", lower_midpoint, node(lower_midpoint).getDofValue('DISPLACEMENT_Y'));
%% Eigenfrequencies

eigensolver = EigensolverStrategy(model);
eigensolver.solve(5);
eigenfrequencies = sort(eigensolver.getEigenfrequencies());
fprintf("Eigenfrequencies: %f\n" , eigenfrequencies);

eigensolver = EigensolverStrategy(model);
eigensolver.solve(5);
modes = eigensolver.getModalMatrix;

%% Plot 

vis = Visualization(model); 
% vis.plotUndeformed();
% vis.plotDeformed();
vis.plotField('modes');


