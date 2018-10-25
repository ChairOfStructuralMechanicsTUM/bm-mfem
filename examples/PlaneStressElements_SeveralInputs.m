close all;
clear all;
clc;


% EXAMPLE: Frame with 2 floors (created out of splitted surfaced by GiD-Input-File)
% Elements: 2d-PlaneStress Elements



% CHOICE of Input-Method %

method = 'manuell';
% method = 'automized';
% method = 'mdpa-input';


if strcmp(method, 'manuell')
    
    %  =================================== %
    %% GEOMETRY BY AUTOMIZED MANUELL INPUT %
    %  =================================== %
    %  up to now: only rectangular plates meshed with quadrilaterals
    
    % Dimensions of the rectangualar plate
    length_x = 10;
    length_y = 2;
    
    % Mesh of the rectangualar plate
    number_Elements_x = 5;
    number_Elements_y = 1;
    
    fprintf("Rectangular Plate with %i x %i Elements\n", number_Elements_y,number_Elements_x);
    
    % Create Nodes
    xCoords=linspace(0,length_x,number_Elements_x+1);
    yCoords=linspace(0,length_y,number_Elements_y+1);
    
    ii = 1;
    
    for i = 1 : length(yCoords)
        for j=1:length(xCoords)
            node(ii) = Node(ii,xCoords(j),yCoords(i),0);
            ii = ii + 1;
        end
    end
    
    nodeArray = node(:)';
    nodeArray.addDof({'DISPLACEMENT_X','DISPLACEMENT_Y'});
    
    
    % Create Elements of the rectangular plate
    ii = 1;
    for j= 1:(length(yCoords)-1)
        for i = ((j-1)*length(xCoords)+1) : ((j-1)*length(xCoords)+length(xCoords)-1)
            ele(ii) = PlaneStressElement2d4n(ii, [node(i) node(i+1) node(i+length(xCoords)+1) node(i+length(xCoords))]);
            ii = ii + 1;
        end
    end
    elementArray = ele(:)';
    
    % Model SetUp
    model = FemModel(nodeArray,elementArray);
    
    % Assignment of Material Properties
    model.getAllElements.setPropertyValue('THICKNESS', 1);
    model.getAllElements.setPropertyValue('YOUNGS_MODULUS',30e6);
    model.getAllElements.setPropertyValue('POISSON_RATIO',0.0);
    model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',2);
    model.getAllElements.setPropertyValue('DENSITY',2699);
    
    
    % Define Boundary Conditions / Supports
    fixedNodes = [1,7];
    model.getNodes([1,7]).fixAllDofs();
    
    % ii = 1;
    % for i=1:length(node)
    %     if node(i).getX == 0 && node(i).getY == 0 || node(i).getX == length_x && node(i).getY == 0
    %         node(i).fixDof('DISPLACEMENT_X');
    %         node(i).fixDof('DISPLACEMENT_Y');
    %         ii= ii+1;
    %     end
    % end
    
    
    % Define Loads
    
    % Point Load
    F=300;
    model.getNodes([6,12]).setDofLoad('DISPLACEMENT_Y',F/2);
    
    
    % % % %% Line load
    % % % for i=1:length(node)
    % % %     if node(i).getX == 0 && node(i).getY == length_x
    % % %         IdUpperLCorner = node(i).getId;
    % % %     end
    % % %     if node(i).getX == length_x && node(i).getY == length_x
    % % %         IdUpperRCorner = node(i).getId;
    % % %     end
    % % % end
    % % %
    % % % addConstLineLoad(IdUpperLCorner,IdUpperRCorner,nodeArray,10,[0 -1 0]);
    
elseif strcmp(method, 'automized')
    
    
    %  =================================== %
    %% GEOMETRY BY AUTOMIZED MANUELL INPUT %
    %  =================================== %
    %  up to now: only rectangular plates meshed with triangles
    
    % Dimensions of the rectangualar plate
    length_x = 10;
    length_y = 2;
    
    % Mesh of the rectangualar plate
    number_Elements_x = 10;
    number_Elements_y = 1;
    
    % Create Nodes
    xCoords=linspace(0,length_x,(number_Elements_x/2+1));
    yCoords=linspace(0,length_y,number_Elements_y+1);
    
    ii = 1;
    
    for i = 1 : length(yCoords)
        for j=1:length(xCoords)
            node(ii) = Node(ii,xCoords(j),yCoords(i),0);
            ii = ii + 1;
        end
    end
    
    nodeArray = node(:)';
    nodeArray.addDof({'DISPLACEMENT_X','DISPLACEMENT_Y'});
    
    
    % Create Elements of the rectangular plate
    ii = 1;
    for j= 1:(length(yCoords)-1)
        for i = ((j-1)*length(xCoords)+1) : ((j-1)*length(xCoords)+length(xCoords)-1)
            ele(ii) = PlaneStressElement2d3n(ii, [node(i) node(i+length(xCoords)) node(i+length(xCoords)+1)]);
            ele(ii+1) = PlaneStressElement2d3n((ii+1), [node(i) node(i+1) node(i+length(xCoords)+1)]);
            ii = ii + 2;
        end
    end
    elementArray = ele(:)';
    
    % Model SetUp
    model = FemModel(nodeArray,elementArray);
    
    vis = Visualization(model);
    vis.setScaling(1000);
    vis.plotUndeformed(1);
    
    % Assignment of Material Properties
    model.getAllElements.setPropertyValue('THICKNESS', 1);
    model.getAllElements.setPropertyValue('YOUNGS_MODULUS',30e6);
    model.getAllElements.setPropertyValue('POISSON_RATIO',0.0);
    model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',1);
    model.getAllElements.setPropertyValue('DENSITY',2699);
    
    
    % Define Boundary Conditions / Supports
    fixedNodes = [1,7];
    model.getNodes([1,7]).fixAllDofs();
    
    % ii = 1;
    % for i=1:length(node)
    %     if node(i).getX == 0 && node(i).getY == 0 || node(i).getX == length_x && node(i).getY == 0
    %         node(i).fixDof('DISPLACEMENT_X');
    %         node(i).fixDof('DISPLACEMENT_Y');
    %         ii= ii+1;
    %     end
    % end
    
    
    % Define Loads
    
    % Point Load
    F=300;
    model.getNodes([6,12]).setDofLoad('DISPLACEMENT_Y',F/2);
    
    
    % % % %% Line load
    % % % for i=1:length(node)
    % % %     if node(i).getX == 0 && node(i).getY == length_x
    % % %         IdUpperLCorner = node(i).getId;
    % % %     end
    % % %     if node(i).getX == length_x && node(i).getY == length_x
    % % %         IdUpperRCorner = node(i).getId;
    % % %     end
    % % % end
    % % %
    % % % addConstLineLoad(IdUpperLCorner,IdUpperRCorner,nodeArray,10,[0 -1 0]);
    
    
elseif strcmp(method, 'mdpa-input')
    
    % ====================================== %
    %% GEOMETRY INPUT FROM GiD (.mdpa-input) %
    % ====================================== %
    
    % MDPA_INPUT Read a large mdpa model and plot it
    % io=MdpaInput('Frame_2floors_foundation_splitted.mdpa'); %specify input file
    io=MdpaInput('Frame_2floors_foundation_complete.mdpa'); %specify input file
    model = io.readModel(); %read the model
    
    
    % Assignment DOFs (here: geometry defined in x-y-plane)
    model.getAllNodes.addDof(["DISPLACEMENT_X","DISPLACEMENT_Y"]);
    
    
    % Assignment of Material Properties
    model.getAllElements.setPropertyValue('THICKNESS', 0.2);
    model.getAllElements.setPropertyValue('YOUNGS_MODULUS',7e8);
    model.getAllElements.setPropertyValue('POISSON_RATIO',0.34);
    model.getAllElements.setPropertyValue('NUMBER_GAUSS_POINT',3);
    model.getAllElements.setPropertyValue('DENSITY',2699);
    
    
    % Set BCs for specified elements:
    % supportedNodes = [370, 371, 372, 374, 378, 382, 384, 389, 393, 397, 403, 406, 413, 417, 424, 429, 435, 443, 450, 456, 461, 467, 474, 482, 490, 495, 501, 508, 513, 517, 520, 522, 523];
    supportedNodes = [349,350,352,357,362,368,374,381,389,399,407,417,426,436,445,457,465,473,478,481,484];
    % model.getModelPart('GENERIC_FixedNodes').getNodes().fixAllDofs();
    model.getNode(supportedNodes).fixAllDofs();
    % barycenter(model.getElement(1)); %is not working for arbitrary hexahedra elements
    
    % % % % % Set BCs for all elements in 'left_support' and 'right_support' (modelParts):
    % % % % %set boundary conditions for all elements in 'left_support' and 'right_support':
    % % % % model.getModelPart('GENERIC_left_support').getNodes().fixAllDofs();
    % % % % model.getModelPart('GENERIC_right_support').getNodes().fixAllDofs();
    
    
    % Set load for specified elements:
    % loads = model.getModelPart('GENERIC_load').getNodes();
    % addPointLoad(loads,1000,[0 0 -1]);
    model.getNodes([1,255]).setDofLoad('DISPLACEMENT_Y',-1000);
    
    
end



% ====================================== %
%% FE-CALCULATIONS (for all input dates) %
% ====================================== %

% Static Calculations

solver = SimpleSolvingStrategy(model);
solver.solve();
% actualDisplacementY = model.getAllNodes.getDofValue('DISPLACEMENT_Y');

% lower_midpoint = number_Elements/2 + 1;
% fprintf("Displacement at Node %i : %f\n", lower_midpoint, node(lower_midpoint).getDofValue('DISPLACEMENT_Y'));

verschiebung_y = model.getAllNodes().getDofValue('DISPLACEMENT_Y');

%% Eigenfrequencies

% eigensolver = EigensolverStrategy(model);
% eigensolver.solve(5);
% eigenfrequencies = sort(eigensolver.getEigenfrequencies());
% fprintf("Eigenfrequencies: %f\n" , eigenfrequencies);

%% Plot 

vis = Visualization(model); 
vis.setScaling(1000);
vis.plotUndeformed();
vis.plotDeformed();

computeElementStress( model.getAllElements, model.getAllNodes,1)

vis.plotField('sigma_xx')
% vis.plotField('sigma_yy')
% vis.plotField('sigma_xy')
