close all; 
clear all;  
clc;

%% example taken from Introduction to FEM (Felippa) 23-5
% a=5; 
% b=2*a;
% node01 = Node(1,0,0,0);
% node02 = Node(2,a,0,0);
% node03 = Node(3,a,b,0);
% node04 = Node(4,0,b,0);

% nodeArray = [node01 node02 node03 node04]; 
% nodeArray.addDof({'DISPLACEMENT_Z', 'ROTATION_X', 'ROTATION_Y'});
% ele = ReissnerMindlinElement3d4n(1,[nodeArray]); 
% 
% 
% elementArray = ele;
% elementArray.setPropertyValue('THICKNESS', 1);
% elementArray.setPropertyValue('YOUNGS_MODULUS', 8000);
% elementArray.setPropertyValue('SHEAR_MODULUS', 3.0769e+03);
% elementArray.setPropertyValue('POISSON_RATIO', 1/3);
% elementArray.setPropertyValue('NUMBER_GAUSS_POINT', 4);
% elementArray.setPropertyValue('DENSITY', 1);
% 
% elementIds = elementArray.getId;
% 
% model = FemModel(nodeArray,elementArray);
% 
% [K,~] = SimpleAssembler.assembleGlobalStiffnessMatrix(model);

%% 
a=linspace(0,1,11);

counter = 1;
for i=1:length(a)
    for j=1:length(a)
        node(counter) = Node(counter,a(j),a(i),0);
        counter = counter + 1; 
    end
end

ii = 1; 
for i = 1 : length(a) : (length(a)-1)^2
    for j= i:i+9
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
arrayfun(@(node) node.fixDof('ROTATION_X'), boundary);
arrayfun(@(node) node.fixDof('ROTATION_Y'), boundary);


% arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), [node(1) node(2) node(3) node(4) ...
%                                                  node(5) node(6) node(7) node(8) ...
%                                                  node(9) node(10) node(11) node(12)...
%                                                  node(22) node(23) node(33) node(34)...
%                                                  node(44) node(45) node(55) node(56)...
%                                                  node(66) node(67) node(77) node(78)...
%                                                  node(88) node(89) node(99) node(100)...
%                                                  node(110) node(111) node(112) node(113)...
%                                                  node(114) node(115) node(116) node(117)...
%                                                  node(118) node(119) node(120) node(121)]);
                                          
% arrayfun(@(node) node.fixDof('ROTATION_X'), [node(1) node(2) node(3) node(4) ...
%                                                  node(5) node(6) node(7) node(8) ...
%                                                  node(9) node(10) node(11) node(12)...
%                                                  node(22) node(23) node(33) node(34)...
%                                                  node(44) node(45) node(55) node(56)...
%                                                  node(66) node(67) node(77) node(78)...
%                                                  node(88) node(89) node(99) node(100)...
%                                                  node(110) node(111) node(112) node(113)...
%                                                  node(114) node(115) node(116) node(117)...
%                                                  node(118) node(119) node(120) node(121)]);
%                                              
% arrayfun(@(node) node.fixDof('ROTATION_Y'), [node(1) node(2) node(3) node(4) ...
%                                                  node(5) node(6) node(7) node(8) ...
%                                                  node(9) node(10) node(11) node(12)...
%                                                  node(22) node(23) node(33) node(34)...
%                                                  node(44) node(45) node(55) node(56)...
%                                                  node(66) node(67) node(77) node(78)...
%                                                  node(88) node(89) node(99) node(100)...
%                                                  node(110) node(111) node(112) node(113)...
%                                                  node(114) node(115) node(116) node(117)...
%                                                  node(118) node(119) node(120) node(121)]);

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

eigensolver = EigensolverStrategy(model)










