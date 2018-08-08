% Dieser Test verknüpft ein homogenes Element vom Typ QuadrilateralElement4d2n
% mit einem porösen Element vom Typ Porous2d4n
% 
% Einschränkungen:
% - der Rand der beiden Elemente ist parallel zur x-Achse
% - es wird angenommen, dass die Elemente fest verbunden sind (bounded)
% - händische Zuordnung der dofs ist nötig, da Reihenfolge im globalen
% dof-Array und in der K-Matrix nicht identisch sind



clear;

node01 = Node(1,0,0);
node02 = Node(2,1,0);
node03 = Node(3,2,0);
node04 = Node(4,2,1);
node05 = Node(5,1,1);
node06 = Node(6,0,1);

nodeArray = [node01 node02 node03 node04 node05 node06];
nodeArray_homogen = [node01 node02 node05 node06];
nodeArray_porous = [node02 node03 node04 node05];

nodeArray_homogen.addDof({'DISPLACEMENT_X', 'DISPLACEMENT_Y'});
nodeArray_porous.addDof({'FRAME_DISPLACEMENT_X', 'FRAME_DISPLACEMENT_Y',...
                'FLUID_DISPLACEMENT_X', 'FLUID_DISPLACEMENT_Y'});

element_homogen = QuadrilateralElement2d4n(1,nodeArray_homogen);
element_porous = Porous2d4n(1,nodeArray_porous);

% Parameter homogenes Element
element_homogen.setPropertyValue('YOUNGS_MODULUS',96);
element_homogen.setPropertyValue('POISSON_RATIO',1/3);
element_homogen.setPropertyValue('NUMBER_GAUSS_POINT',2);
element_homogen.setPropertyValue('DENSITY',7860);

% Parameter poröses Element
%Values
element_porous.setPropertyValue('FREQUENCY',100); %angenommen 100 Hz
element_porous.setPropertyValue('NUMBER_GAUSS_POINT',2);
%Frame
element_porous.setPropertyValue('FRAME_LAME_PARAMETER_LAMBDA',114400);
element_porous.setPropertyValue('FRAME_LAME_PARAMETER_MU',28600);
%Fluid
element_porous.setPropertyValue('HEAT_CAPACITY_RATIO',1.4);
element_porous.setPropertyValue('AMBIENT_FLUID_STANDARD_PRESSURE',101300); %angenommen 1013 hPa
element_porous.setPropertyValue('AMBIENT_FLUID_VISCOSITY',1.84*10^-5);
element_porous.setPropertyValue('AMBIENT_FLUID_DENSITY',1.21);
element_porous.setPropertyValue('PRANDTL_NUMBER',0.71);
%Porous
element_porous.setPropertyValue('THERMAL_CHAR_LENGTH',226*10^-6);
element_porous.setPropertyValue('POROSITY',0.9);
%Mass Matrix
element_porous.setPropertyValue('FRAME_DENSITY',300); % berechnet aus Porosität und Angabe
element_porous.setPropertyValue('TORTUOSITY',7.8);
element_porous.setPropertyValue('VISCOUS_CHAR_LENGTH',226*10^-6);
element_porous.setPropertyValue('STATIC_FLOW_RESISTIVITY',2500);

K_h = element_homogen.computeLocalStiffnessMatrix;
K_p = element_porous.computeLocalStiffnessMatrix;

% ID-Zuweisung an dofs
dofArray = arrayfun(@(node) node.getDofArray, nodeArray, 'UniformOutput', false);
dofArray = [dofArray{:}];
for iDofs=1:length(dofArray)
    dofArray(iDofs).setId(iDofs);
end

%Aufstellung Gesamtmatrix K aus K für poröses und homogenes Element
size_h = size(K_h,1);
size_p = size(K_p,1);
K = zeros(size_h+size_p);
K([1:size_h],[1:size_h])=K_h;
K([size_h+1:end],[size_h+1:end])=K_p;


%Aufbringung der Übergangsbedingungen und Reduzierung der Matrix
Zuordnung = [(1:24);[1 2 3 4 17 18 23 24 7 8 11 12 15 16 21 22 5 6 9 10 13 14 19 20]];
%Zuordnung ist notwendig, da Reihenfolge der dofs in globalem dof-Array
%nicht mit der Reihenfolge in der K-Matrix übereinstimmt
delete = zeros(6,1);
count = 0;
for iNode = 1: length(nodeArray)
    node = nodeArray(iNode);
    dofs_node = node.getDofArray;
    ndofs_node = length(dofs_node);
    if ndofs_node == 6
        count=count+1;
        id_homogen_x = dofs_node(1).getId;
        id_homogen_x = Zuordnung(2,id_homogen_x);
        id_homogen_y = dofs_node(2).getId;
        id_homogen_y = Zuordnung(2,id_homogen_y);
        id_fluid_x = dofs_node(3).getId;
        id_fluid_x = Zuordnung(2,id_fluid_x);
        id_fluid_y = dofs_node(4).getId;
        id_fluid_y = Zuordnung(2,id_fluid_y);
        id_frame_x = dofs_node(5).getId;
        id_frame_x = Zuordnung(2,id_frame_x);
        id_frame_y = dofs_node(6).getId;
        id_frame_y = Zuordnung(2,id_frame_y);
        

        K(id_homogen_x,:)=K(id_homogen_x,:)+K(id_frame_x,:);
        K(:,id_homogen_x)=K(:,id_homogen_x)+K(:,id_frame_x);
        K(id_homogen_x,:)=K(id_homogen_x,:)+K(id_fluid_x,:);
        K(:,id_homogen_x)=K(:,id_homogen_x)+K(:,id_fluid_x);
        K(id_homogen_y,:)=K(id_homogen_y,:)+K(id_frame_y,:);
        K(:,id_homogen_y)=K(:,id_homogen_y)+K(:,id_frame_y);
        
        delete([3*count-2,3*count-1,3*count]) = [id_frame_x,id_frame_y,id_fluid_x];
    end
end
K(delete,:) = [];
K(:,delete) = [];