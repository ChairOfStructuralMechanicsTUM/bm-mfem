classdef MixedDisplacementElement2d4n < QuadrilateralElement
        
    properties (Access = private)

    end



    
    methods
%%         
    % Constructor
function mixed2d4n = MixedDisplacementElement2d4n(id,nodeArray)
            
    requiredPropertyNames = cellstr(["OMEGA","NUMBER_GAUSS_POINT", "LAMBDA", ...
                "MU","ETA_S","DENSITY_S","DENSITY_F","ETA_F","PRESSURE_0","HEAT_CAPACITY_RATIO",...
                 "PRANDL_NUMBER","POROSITY","TORTUOSITY","STATIC_FLOW_RESISTIVITY",...
                "VISCOUS_CHARACT_LENGTH","THERMAL_CHARACT_LENGTH"]);
            
    % Define the arguments for the super class constructor call
    if nargin == 0
        super_args = {};
        elseif nargin == 2
            if ~(length(nodeArray) == 4 && isa(nodeArray,'Node'))
               error('problem with the nodes in element %d', id);
            end
            super_args = {id, nodeArray, requiredPropertyNames};
    end
            
    % Call the super class contructor
    mixed2d4n@QuadrilateralElement(super_args{:});
    mixed2d4n.dofNames = cellstr(["DISPLACEMENT_SOLID_X", "DISPLACEMENT_SOLID_Y", ...
               "PRESSURE_FLUID"]);
end
%%   


    % Initialization
function initialize(mixed2d4n)
    mixed2d4n.lengthX = computeLength(mixed2d4n.nodeArray(1).getCoords, ...
    mixed2d4n.nodeArray(2).getCoords);        
    mixed2d4n.lengthY = computeLength(mixed2d4n.nodeArray(1).getCoords, ...
    mixed2d4n.nodeArray(4).getCoords);
            
    checkConvexity(mixed2d4n);
end
%%    


function responseDoF = getResponseDofArray(element, step)
            
    responseDoF = zeros(12,1);
    for itNodes = 1:1:4
        nodalDof = element.nodeArray(itNodes).getDofArray;
        nodalDof = nodalDof.';
                      
        for itDof = 3:(-1):1
            responseDoF(3*itNodes-(itDof-1),1) = nodalDof(4-itDof).getValue(step);
                
        end
        
    end
end          
%%    


    % Shape Function and Derivatives     
function [N_mat, N, B_s, B, J] = computeShapeFunction(mixed2d4n,xi,eta)

            N = [(1-xi)*(1-eta)/4    (1+xi)*(1-eta)/4    (1+xi)*(1+eta)/4    (1-xi)*(1+eta)/4];
            
            N_Diff_Par = [-(1-eta)/4    (1-eta)/4   (1+eta)/4   -(1+eta)/4 
                -(1-xi)/4     -(1+xi)/4   (1+xi)/4    (1-xi)/4];
            
            N_mat = sparse(2,8);
            N_mat(1,1:2:end) = N(:);
            N_mat(2,2:2:end) = N(:);
            
    % Coordinates of the nodes forming one element
            ele_coords = zeros(4,2);
            for i=1:4
                ele_coords(i,1) = mixed2d4n.nodeArray(i).getX;
                ele_coords(i,2) = mixed2d4n.nodeArray(i).getY;
            end
            
             
    % Jacobian
            J = N_Diff_Par * ele_coords;
            
    % Calculation of B-Matrix for Solid Phase
            B = J\N_Diff_Par; 
            Bx = B(1,1:4);
            By = B(2,1:4);
%             %%%% [8x8] = [8x3]*[3x3]*[3x8]
             B_s = [Bx(1),0,Bx(2),0,Bx(3),0,Bx(4),0;
             0,By(1),0,By(2),0,By(3),0,By(4);
             By(1),Bx(1),By(2),Bx(2),By(3),Bx(3),By(4),Bx(4)];
%     
%     %%%% Q:How to set up the other B-Matrices?
%             %%%% [8x4] = [8x3]*[3x3]*[3x4]
%             B_sf1 = zeros(8,3);
%             B_sf2 = zeros(3,4);
%             
%             %%%% [4x8] = [4x3]*[3x3]*[3x8]
%             B_fs1 = zeros(4,3);
%             B_fs2 = zeros(3,8);
% 
%             %%%% [4x4] = [4x3]*[3x3]*[3x4]
%             B_f = zeros(3,4);
%             
            
end
%%


function [totalElementStiffnessMatrix] = computeLocalStiffnessMatrix(mixed2d4n)
     

    ETA_S = mixed2d4n.getPropertyValue('ETA_S');
    LAMBDA = mixed2d4n.getPropertyValue('LAMBDA');
    MU = mixed2d4n.getPropertyValue('MU');
    OMEGA = mixed2d4n.getPropertyValue('OMEGA');
    POROSITY = mixed2d4n.getPropertyValue('POROSITY');
            
        
    HEAT_CAPACITY_RATIO= mixed2d4n.getPropertyValue('HEAT_CAPACITY_RATIO');    
    PRESSURE_0 = mixed2d4n.getPropertyValue('PRESSURE_0');
    ETA_F = mixed2d4n.getPropertyValue('ETA_F');
    PRANDL_NUMBER = mixed2d4n.getPropertyValue('PRANDL_NUMBER');
    THERMAL_CHARACT_LENGTH = mixed2d4n.getPropertyValue('THERMAL_CHARACT_LENGTH');
    DENSITY_F = mixed2d4n.getPropertyValue('DENSITY_F');
           
    p = mixed2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
         
    % Calculate Bulk-Modulus K_f:
    K_f = (HEAT_CAPACITY_RATIO*PRESSURE_0)/(HEAT_CAPACITY_RATIO-(HEAT_CAPACITY_RATIO-1)*...
    (1+(8*ETA_F)/(i*OMEGA*PRANDL_NUMBER*THERMAL_CHARACT_LENGTH^2*DENSITY_F)*...
    (1+(i*OMEGA*PRANDL_NUMBER*THERMAL_CHARACT_LENGTH^2*DENSITY_F)/(16*ETA_F))^(0.5))^(-1));
        
    % Hysteretic proportional damping model:
    LameCoeff = (1+i*ETA_S)*LAMBDA;
    ShearModulus = (1+i*ETA_S)*MU;
            
            
    % Biot's Elasticity Coefficients:
    R = POROSITY*K_f;
    Q = (1-POROSITY)*K_f;
    A = LameCoeff+(1-POROSITY)^2/POROSITY*K_f;
            
            
    % Elasticity Tensors:
    D_s = [LameCoeff+2*ShearModulus,LameCoeff,0;LameCoeff,LameCoeff+2*ShearModulus,0;0,0,ShearModulus];
    D_sf = -POROSITY*Q/R*eye(3);
    D_f = -POROSITY*eye(3);
    
    %%%% No coupling-term in fluid stress-strain relation
    D_fs = zeros(3);
            
 
    % Get gauss integration factors:
    [w,g]=returnGaussPoint(p);
           
            
    % Generate empty partial stiffness-matrices:
    stiffnessMatrix_s=sparse(8,8);
    stiffnessMatrix_f=sparse(4,4);
    stiffnessMatrix_sf = sparse(8,4);
    stiffnessMatrix_fs = sparse(4,8);
            
    % Execute full Gauss-Integration to obtain partial
    % stiffness-matrices for solid, fluid and coupling phases:
    for xi = 1 : p
        for eta = 1 : p
        [N_mat, N, B_s, B, J] = computeShapeFunction(mixed2d4n,g(xi),g(eta));

        stiffnessMatrix_s = stiffnessMatrix_s +...
        B_s' * D_s * B_s * det(J) * w(xi) * w(eta);
           
        stiffnessMatrix_f = stiffnessMatrix_f +             ...
        B'*B - N'*N * det(J) * w(xi) * w(eta);

         stiffnessMatrix_sf = stiffnessMatrix_sf +             ...
         B_sf1 * D_sf * B_sf2 * det(J) * w(xi) * w(eta); 
%                         
%     %%%% [4x8] = [4x3]*[3x3]*[3x8] = zeros(4,8)
%         stiffnessMatrix_fs = stiffnessMatrix_fs +             ...
%         B_fs1 * D_fs * B_fs2 * det(J) * w(xi) * w(eta); 
    
        end
    end
            
            
            
    % Generate empty total element stiffness matrix
    totalElementStiffnessMatrix = zeros(12,12);
            
    % Assemble total element stiffness matrix:
    totalElementStiffnessMatrix = [stiffnessMatrix_s, stiffnessMatrix_sf;...
    stiffnessMatrix_fs,stiffnessMatrix_f];
       
     
    
end    
%%


function [totalElementMassMatrix] = computeLocalMassMatrix(mixed2d4n)

    STATIC_FLOW_RESISTIVITY = mixed2d4n.getPropertyValue('STATIC_FLOW_RESISTIVITY');
    POROSITY = mixed2d4n.getPropertyValue('POROSITY');
    OMEGA = mixed2d4n.getPropertyValue('OMEGA');
    TORTUOSITY = mixed2d4n.getPropertyValue('TORTUOSITY');
    ETA_F = mixed2d4n.getPropertyValue('ETA_F');
    DENSITY_F = mixed2d4n.getPropertyValue('DENSITY_F');
    DENSITY_S = mixed2d4n.getPropertyValue('DENSITY_S');
    VISCOUS_CHARACT_LENGTH = mixed2d4n.getPropertyValue('VISCOUS_CHARACT_LENGTH');
    p = mixed2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            
    % G_J(OMEGA) = flow resistivity of air particles in the pores:
    AirFlowResistivity = (1+(4*1i*OMEGA*TORTUOSITY^2*ETA_F*DENSITY_F)/...
    (STATIC_FLOW_RESISTIVITY^2*VISCOUS_CHARACT_LENGTH^2*POROSITY^2))^(0.5);
            
            
    % Viscous Drag accounts for viscous body forces interacting between solid and fluid phases, 
    % proportional to the relative velocity:
    ViscousDrag = STATIC_FLOW_RESISTIVITY*POROSITY^2*AirFlowResistivity;
            

    % Inertial coupling term, related to the tortuosity, increases the fluid density  to model the motion
    % of fluid particles vibrating around the structural frame:
    Rho_a = POROSITY*DENSITY_F*(TORTUOSITY-1);
            
    % Equivalent densities for expressing the elastodynamic coupled equations in condensed form:
    Rho_sf = -Rho_a+1i*(ViscousDrag)/(OMEGA);
    Rho_s = (1-POROSITY)*DENSITY_S-Rho_sf;
    Rho_f = POROSITY*DENSITY_F-Rho_sf;
            
    %%%% Generate empty mass-matrices:
    M_s = zeros(8,8);
    M_sf = zeros(8,4);
    M_fs = zeros(4,8);
    M_f = zeros(4,4);
    
    totalElementMassMatrix = zeros(12,12);

    % Calculate geometric-only mass matrix:
    [w,g]=returnGaussPoint(p);
            
    for n=1:p
    xi=g(n);
        for m=1:p
        eta=g(m);            
        [N_mat, ~, ~, J] = computeShapeFunction(mixed2d4n,xi,eta);        
        M = M + (w(n)*w(m)*transpose(N_mat)*N_mat*det(J));
                    
        end
    end
            
            
    % Calculate solid-, fluid- and coupling-mass matrices:
    M_s = Rho_s*M;
    
    
    %%%%
    M_f = Rho_f*M;
    Msf = Rho_sf*M;
    Mfs = Rho_sf*M;
    
    % Assemble total element mass matrix:
    totalElementMassMatrix = [M_s,M_sf;M_fs,M_f];
         
             
end
%%  


    % Define possible DOF in the node-array:
function dofs = getDofList(element)
            dofs([1 3 5 7])     = element.nodeArray.getDof('DISPLACEMENT_SOLID_X');
            dofs([2 4 6 8])     = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y');
            dofs([9 11 13 15])  = element.nodeArray.getDof('DISPLACEMENT_FLUID_X');
            dofs([10 12 14 16]) = element.nodeArray.getDof('DISPLACEMENT_FLUID_Y');

end     
%%        


function vals = getValuesVector(element, step)
            vals = zeros(1,16);
            
            vals([1 3 5 7])     = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_X',step);
            vals([2 4 6 8])     = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_Y',step);
            vals([9 11 13 15])  = element.nodeArray.getDofValue('DISPLACEMENT_FLUID_X',step);
            vals([10 12 14 16]) = element.nodeArray.getDofValue('DISPLACEMENT_FLUID_Y',step);

end           
%%     


function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,16);
            
            [~, vals([1 3 5 7]), ~]     = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, vals([2 4 6 8]), ~]     = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, vals([9 11 13 15]), ~]  = element.nodeArray.getDof('DISPLACEMENT_FLUID_X').getAllValues(step);
            [~, vals([10 12 14 16]), ~] = element.nodeArray.getDof('DISPLACEMENT_FLUID_Y').getAllValues(step);

            
end
%%    


function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,16);
            
            [~, ~, vals([1 3 5 7])]     = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, ~, vals([2 4 6 8])]     = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, ~, vals([9 11 13 15])]  = element.nodeArray.getDof('DISPLACEMENT_FLUID_X').getAllValues(step);
            [~, ~, vals([10 12 14 16])] = element.nodeArray.getDof('DISPLACEMENT_FLUID_Y').getAllValues(step);
          
end
%% 


function [totalElementDampingMatrix] = computeLocalDampingMatrix(mixed2d4n)
            totalElementDampingMatrix = zeros(16,16);
end
%% 


function F = computeLocalForceVector(mixed2d4n)
            F = zeros(1,16);
end
        
    end
end