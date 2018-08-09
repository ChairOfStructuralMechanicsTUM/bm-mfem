% Dieser Test verknüpft ein homogenes Element vom Typ QuadrilateralElement2d4n
% mit einem porösen Element vom Typ Porous2d4n
%
% Einschränkungen:
% - der Rand der beiden Elemente ist parallel zur x-Achse
% - es wird angenommen, dass die Elemente fest verbunden sind (bounded)

tic;
close all;
 
%% Creating a new model 

% Adding nodes
model = FemModel();
 
model.addNewNode(1,0,0,0);
model.addNewNode(2,1,0,0);
model.addNewNode(3,1,1,0);
model.addNewNode(4,0,1,0);
model.addNewNode(5,1,2,0);
model.addNewNode(6,0,2,0);
 
% Creating different types of elements
model.addNewElement('QuadrilateralElement2d4n',1,[1 2 3 4]);
model.addNewElement('Porous2d4n',2,[4 3 5 6]);
 
% Defining model parts
model.addNewModelPart('Homogen',[1 2 3 4],1);
model.addNewModelPart('Porous',[4 3 5 6],2);
 
% Setting property values and DOF for model parts
mp1 = model.getModelPart('Homogen');
mp1.getElements().setPropertyValue('YOUNGS_MODULUS',96);
mp1.getElements().setPropertyValue('POISSON_RATIO',1/3);
mp1.getElements().setPropertyValue('NUMBER_GAUSS_POINT',2);
mp1.getElements().setPropertyValue('DENSITY',7860);
 
mp2 = model.getModelPart('Porous');
%Values
mp2.getElements.setPropertyValue('FREQUENCY',100); %angenommen 100 Hz
mp2.getElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
%Frame
mp2.getElements.setPropertyValue('FRAME_LAME_PARAMETER_LAMBDA',114400);
mp2.getElements.setPropertyValue('FRAME_LAME_PARAMETER_MU',28600);
%Fluid
mp2.getElements.setPropertyValue('HEAT_CAPACITY_RATIO',1.4);
mp2.getElements.setPropertyValue('AMBIENT_FLUID_STANDARD_PRESSURE',101300); %angenommen 1013 hPa
mp2.getElements.setPropertyValue('AMBIENT_FLUID_VISCOSITY',1.84*10^-5);
mp2.getElements.setPropertyValue('AMBIENT_FLUID_DENSITY',1.21);
mp2.getElements.setPropertyValue('PRANDTL_NUMBER',0.71);
%Porous
mp2.getElements.setPropertyValue('THERMAL_CHAR_LENGTH',226*10^-6);
mp2.getElements.setPropertyValue('POROSITY',0.9);
%Mass Matrix
mp2.getElements.setPropertyValue('FRAME_DENSITY',300);
mp2.getElements.setPropertyValue('TORTUOSITY',7.8);
mp2.getElements.setPropertyValue('VISCOUS_CHAR_LENGTH',226*10^-6);
mp2.getElements.setPropertyValue('STATIC_FLOW_RESISTIVITY',2500);
 
%Add all dofs to all elements
mp1.getNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y",...
    "FRAME_DISPLACEMENT_X", "FRAME_DISPLACEMENT_Y", ...
    "FLUID_DISPLACEMENT_X", "FLUID_DISPLACEMENT_Y"]);
mp2.getNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y",...
    "FRAME_DISPLACEMENT_X", "FRAME_DISPLACEMENT_Y", ...
    "FLUID_DISPLACEMENT_X", "FLUID_DISPLACEMENT_Y"]);
           
 
%% Setting BC
% model.getNode(1).fixAllDofs();
% model.getNode(4).fixAllDofs();
% model.getNode(6).fixAllDofs();
 
model.initialize();

%% Computing elemental matrices and expanding them on full size

ndofs = length(model.getDofArray);

% matrix of homogenous part
K_h = mp1.getElements.computeLocalStiffnessMatrix;
K_h_full = zeros(ndofs);
id_vector_homogen = zeros(size(K_h,1),1);
node_array_mp1 = mp1.getElements.getNodeArray;

for i=1:length(node_array_mp1)
    node = node_array_mp1(i);
        id_vector_homogen(2*i-1) = node.getDof('DISPLACEMENT_X').getId;
        id_vector_homogen(2*i) = node.getDof('DISPLACEMENT_Y').getId;
end
K_h_full(id_vector_homogen,id_vector_homogen)= K_h;

% matrix of porous part
K_p = mp2.getElements.computeLocalStiffnessMatrix;
K_p_full = zeros(ndofs);
id_vector_porous = zeros(size(K_p,1),1);
node_array_mp2 = mp2.getElements.getNodeArray;

for i=1:length(node_array_mp2)
    node = node_array_mp2(i);
        id_vector_porous(2*i-1) = node.getDof('FRAME_DISPLACEMENT_X').getId;
        id_vector_porous(2*i) = node.getDof('FRAME_DISPLACEMENT_Y').getId;
        id_vector_porous((size(K_p,1)/2)+2*i-1) = node.getDof('FLUID_DISPLACEMENT_X').getId;
        id_vector_porous((size(K_p,1)/2)+2*i) = node.getDof('FLUID_DISPLACEMENT_Y').getId;
end
K_p_full(id_vector_porous,id_vector_porous)= K_p;
    
K_full=K_p_full+K_h_full;

%% Determining intersection (common nodes)
ids1 = mp1.getNodes().getId();
ids2 = mp2.getNodes().getId();
kopplung = intersect(ids1,ids2);

%% Applying coupling conditions and reducing matrix

% Applying coupling conditions
delete = zeros(length(kopplung),1);
node_array = model.getAllNodes;

for iNode = 1:length(kopplung)
    node_id = kopplung(iNode);
    node = node_array(node_id);
    id_homogen_x = node.getDof('DISPLACEMENT_X').getId;
    id_homogen_y = node.getDof('DISPLACEMENT_Y').getId;
    id_frame_x = node.getDof('FRAME_DISPLACEMENT_X').getId;
    id_frame_y = node.getDof('FRAME_DISPLACEMENT_Y').getId;
    id_fluid_y = node.getDof('FLUID_DISPLACEMENT_Y').getId;
        
    K_full(id_homogen_x,:)=K_full(id_homogen_x,:)+K_full(id_frame_x,:);
    K_full(:,id_homogen_x)=K_full(:,id_homogen_x)+K_full(:,id_frame_x);
    K_full(id_homogen_y,:)=K_full(id_homogen_y,:)+K_full(id_frame_y,:);
    K_full(:,id_homogen_y)=K_full(:,id_homogen_y)+K_full(:,id_frame_y);
    K_full(id_homogen_y,:)=K_full(id_homogen_y,:)+K_full(id_fluid_y,:);
    K_full(:,id_homogen_y)=K_full(:,id_homogen_y)+K_full(:,id_fluid_y);
    
    delete([3*iNode-2,3*iNode-1,3*iNode]) = [id_frame_x,id_frame_y,id_fluid_y];
end

%reducing matrix
for iID1=1:length(ids1)
    if ~ismember(ids1(iID1),kopplung)
        delete(end+1)=node_array_mp1(iID1).getDof('FRAME_DISPLACEMENT_X').getId;
        delete(end+1)=node_array_mp1(iID1).getDof('FRAME_DISPLACEMENT_Y').getId;
        delete(end+1)=node_array_mp1(iID1).getDof('FLUID_DISPLACEMENT_X').getId;
        delete(end+1)=node_array_mp1(iID1).getDof('FLUID_DISPLACEMENT_Y').getId;
    end
end
   
for iID2=1:length(ids2)
    if ~ismember(ids2(iID2),kopplung)
        delete(end+1)=node_array_mp2(iID2).getDof('DISPLACEMENT_X').getId;
        delete(end+1)=node_array_mp2(iID2).getDof('DISPLACEMENT_Y').getId;
    end
end

K_full(delete,:) = [];
K_full(:,delete) = [];

time = toc;

