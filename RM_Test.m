
clear; clc;
%example taken from Introduction to FEM (Felippa) 23-5

%Initialize Nodes
% node01 = Node(1,0,0,0);
% node02 = Node(2,0.5,0,0);
% node03 = Node(3,1,0,0);
% node04 = Node(4,0,0.5,0);
% node05 = Node(5,0.5,0.5,0);
% node06 = Node(6,1,0.5,0);
% node07 = Node(7,0,1,0);
% node08 = Node(8,0.5,1,0);
% node09 = Node(9,1,1,0);

a=linspace(0,1,11);
counter = 1;

for i=1:length(a)
    for j=1:length(a)
        node(counter) = Node(counter,a(j),a(i),0);
        counter = counter + 1; 
    end
end

counter = 1;
for i=1:10
    ele(counter) = ReissnerMindlinElement3d4n(counter,[node(i) node(i+1) node(i+length(a)+1) node(i+length(a))]);
%     scatter(node(i).getX(), node(i).getY());
%     scatter(node(i+1).getX(), node(i+1).getY());
%     scatter(node(i+length(a)+1).getX(), node(i+length(a)+1).getY());
%     scatter(node(i+length(a)).getX(), node(i+length(a)).getY());
%     
%     text(node(i).getX(), node(i).getY(), num2str(counter));
%     text(node(i+1).getX(), node(i+1).getY(), num2str(counter));
%     text(node(i+length(a)+1).getX(), node(i+length(a)+1).getY(), num2str(counter));
%     text(node(i+length(a)).getX(), node(i+length(a)).getY(), num2str(counter));
    counter = counter + 1;
end

for i=12:21
    ele(counter) = ReissnerMindlinElement3d4n(counter,[node(i) node(i+1) node(i+length(a)+1) node(i+length(a))]);
    counter = counter + 1;
end


for i=23:32
    ele(counter) = ReissnerMindlinElement3d4n(counter,[node(i) node(i+1) node(i+length(a)+1) node(i+length(a))]);
    counter = counter + 1;
end

for i=34:43
    ele(counter) = ReissnerMindlinElement3d4n(counter,[node(i) node(i+1) node(i+length(a)+1) node(i+length(a))]);
    counter = counter + 1;
end
for i=45:54
    ele(counter) = ReissnerMindlinElement3d4n(counter,[node(i) node(i+1) node(i+length(a)+1) node(i+length(a))]);

    counter = counter + 1;
end
for i=56:65
    ele(counter) = ReissnerMindlinElement3d4n(counter,[node(i) node(i+1) node(i+length(a)+1) node(i+length(a))]);

    counter = counter + 1;
end

for i=67:76
    ele(counter) = ReissnerMindlinElement3d4n(counter,[node(i) node(i+1) node(i+length(a)+1) node(i+length(a))]);
    counter = counter + 1;
end

for i=78:87
    ele(counter) = ReissnerMindlinElement3d4n(counter,[node(i) node(i+1) node(i+length(a)+1) node(i+length(a))]);
    counter = counter + 1;
end
for i=89:98
    ele(counter) = ReissnerMindlinElement3d4n(counter,[node(i) node(i+1) node(i+length(a)+1) node(i+length(a))]);
    counter = counter + 1;
end
for i=100:109
    ele(counter) = ReissnerMindlinElement3d4n(counter,[node(i) node(i+1) node(i+length(a)+1) node(i+length(a))]);
    counter = counter + 1;
end

nodeArray = [node(:)'];
% % nodeArray = [node01 node02 node03 node04 node05 node06 node07 node08 node09];
% 
nodeArray.addDof({'DISPLACEMENT_Z', 'ROTATION_X', 'ROTATION_Y'});
% 
% 
% 
% % ele01 = ReissnerMindlinElement3d4n(1,[node01 node02 node05 node04]);
% % ele02 = ReissnerMindlinElement3d4n(2,[node02 node03 node06 node05]); 
% % ele03 = ReissnerMindlinElement3d4n(3,[node04 node05 node08 node07]);
% % ele04 = ReissnerMindlinElement3d4n(4,[node05 node06 node09 node08]);
% % ele01 = ReissnerMindlinElement3d4n(1,[node(1) node(2) node(5) node(4)]);
% % ele02 = ReissnerMindlinElement3d4n(2,[node(2) node(3) node(6) node(5)]); 
% % ele03 = ReissnerMindlinElement3d4n(3,[node(4) node(5) node(8) node(7)]);
% % ele04 = ReissnerMindlinElement3d4n(4,[node(5) node(6) node(9) node(8)]);
% 
% 
% 
arrayfun(@(node) node.fixDof('DISPLACEMENT_Z'), [node(1) node(2) node(3) node(4) ...
                                                 node(6) node(7) node(8) node(9) ...
                                                 node(10) node(11) node(12) node(23)...
                                                 node(34) node(45) node(56) node(67) ...
                                                 node(78) node(89) node(100) node(111)...
                                                 node(112) node(113) node(114) node(115)...
                                                 node(116) node(117) node(118) node(119)...
                                                 node(120) node(121) node(110) node(99) ...
                                                 node(88) node(77) node(66) node(55) ...
                                                 node(44) node(33) node(22)]);
                                             
arrayfun(@(node) node.fixDof('ROTATION_X'), [node(1) node(2) node(3) node(4) ...
                                                 node(6) node(7) node(8) node(9) ...
                                                 node(10) node(11) node(12) node(23)...
                                                 node(34) node(45) node(56) node(67) ...
                                                 node(78) node(89) node(100) node(111)...
                                                 node(112) node(113) node(114) node(115)...
                                                 node(116) node(117) node(118) node(119)...
                                                 node(120) node(121) node(110) node(99) ...
                                                 node(88) node(77) node(66) node(55) ...
                                                 node(44) node(33) node(22)]);
                                             
arrayfun(@(node) node.fixDof('ROTATION_Y'), [node(1) node(2) node(3) node(4) ...
                                                 node(6) node(7) node(8) node(9) ...
                                                 node(10) node(11) node(12) node(23)...
                                                 node(34) node(45) node(56) node(67) ...
                                                 node(78) node(89) node(100) node(111)...
                                                 node(112) node(113) node(114) node(115)...
                                                 node(116) node(117) node(118) node(119)...
                                                 node(120) node(121) node(110) node(99) ...
                                                 node(88) node(77) node(66) node(55) ...
                                                 node(44) node(33) node(22)]);

elementArray = [ele(:)'];


% % elementArray = [ele01 ele02 ele03 ele04];

elementArray.setPropertyValue('THICKNESS', 0.01);
elementArray.setPropertyValue('YOUNGS_MODULUS', 10920);
elementArray.setPropertyValue('SHEAR_MODULUS', 4200);
elementArray.setPropertyValue('POISSON_RATIO', 0.3);
elementArray.setPropertyValue('NUMBER_GAUSS_POINT', 4);
elementArray.setPropertyValue('DENSITY', 1);
% elementArray.setPropertyValue('SHEAR_CORRECTION_FACTOR', 0.8601);



elementIds = elementArray.getId;

model = FemModel(nodeArray,elementArray);

% assembling = SimpleAssembler(model);
[K,~] = SimpleAssembler.assembleGlobalStiffnessMatrix(model); 
M = SimpleAssembler.assembleGlobalMassMatrix(model);


[V,D] = eig(K,M);
D = diag(sqrt(D)*1*sqrt(1/4200))
% D = diag(sqrt(D))
[D,ii] = sort(D);
VV = V(:,ii);





