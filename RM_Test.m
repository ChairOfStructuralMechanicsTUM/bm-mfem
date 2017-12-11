close all; 
clear all;  


%% example taken from Introduction to FEM (Felippa) 23-5
% a=5; 
% b=2*a;
% node01 = Node(1,0,0,0);
% node02 = Node(2,a,0,0);
% node03 = Node(3,a,b,0);
% node04 = Node(4,0,b,0);
% 
% nodeArray = [node01 node02 node03 node04]; 
% nodeArray.addDof({'DISPLACEMENT_Z', 'ROTATION_X', 'ROTATION_Y'});
% ele = ReissnerMindlinElement3d4n(1,[nodeArray]); 
% 
% 
% elementArray = ele;
% elementArray.setPropertyValue('THICKNESS', 1);
% elementArray.setPropertyValue('YOUNGS_MODULUS', 8000);
% elementArray.setPropertyValue('SHEAR_MODULUS', 0);
% elementArray.setPropertyValue('POISSON_RATIO', 1/3);
% elementArray.setPropertyValue('NUMBER_GAUSS_POINT', 4);
% elementArray.setPropertyValue('DENSITY', 1);
% 
% 
% model = FemModel(nodeArray,elementArray);
% 
% [K,~] = SimpleAssembler.assembleGlobalStiffnessMatrix(model);

%% 

number_Elements = 30;
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

arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), boundary);
% arrayfun(@(node) node.fixDof('ROTATION_X'), boundary);
% arrayfun(@(node) node.fixDof('ROTATION_Y'), boundary);

elementArray = ele(:)';

elementArray.setPropertyValue('THICKNESS', 0.01);
elementArray.setPropertyValue('YOUNGS_MODULUS', 10920);
elementArray.setPropertyValue('SHEAR_MODULUS', 4200);
elementArray.setPropertyValue('POISSON_RATIO', 0.3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT', 4);
elementArray.setPropertyValue('DENSITY', 1);
% elementArray.setPropertyValue('SHEAR_CORRECTION_FACTOR', 0.8601);

elementIds = elementArray.getId;

model = FemModel(nodeArray,elementArray);

vis = Visualization(model); 
vis.plotUndeformed(); 

[K,~] = SimpleAssembler.assembleGlobalStiffnessMatrix(model); 
M = SimpleAssembler.assembleGlobalMassMatrix(model);

eigensolver = EigensolverStrategy(model);
eigensolver.solve(5);
eigenfrequencies = eigensolver.getEigenfrequencies();
disp(eigenfrequencies);

