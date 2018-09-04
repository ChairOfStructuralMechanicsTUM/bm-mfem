classdef AtallaElement2d4n < PorousElement2d4n
        
    properties (Access = private)
    end
    
    
    methods
        
        % Constructor
        function mixed2d4n = AtallaElement2d4n(id,nodeArray)
            
            requiredPropertyNames = cellstr(["OMEGA","NUMBER_GAUSS_POINT", "LAMBDA_S", ...
                "MUE_S","ETA_S","DENSITY_S","DENSITY_F","ETA_F","PRESSURE_0_F","HEAT_CAPACITY_RATIO_F",...
                 "PRANDL_NUMBER_F","POROSITY","TORTUOSITY","STATIC_FLOW_RESISTIVITY",...
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
            mixed2d4n@PorousElement2d4n(super_args{:});
            mixed2d4n.dofNames = cellstr(["DISPLACEMENT_SOLID_X", "DISPLACEMENT_SOLID_Y", ...
               "PRESSURE_FLUID"]);
        end
 


        % Initialization
        function initialize(mixed2d4n)
            mixed2d4n.lengthX = computeLength(mixed2d4n.nodeArray(1).getCoords, ...
                mixed2d4n.nodeArray(2).getCoords);        
            
            mixed2d4n.lengthY = computeLength(mixed2d4n.nodeArray(1).getCoords, ...
                mixed2d4n.nodeArray(4).getCoords);
            
            checkConvexity(mixed2d4n);
        end
  

% 
%         function responseDoF = getResponseDofArray(element, step)
%             
%             responseDoF = zeros(12,1);
%             for itNodes = 1:1:4
%                 nodalDof = element.nodeArray(itNodes).getDofArray;
%                 nodalDof = nodalDof.';
%                       
%             for itDof = 3:(-1):1
%                 responseDoF(3*itNodes-(itDof-1),1) = nodalDof(4-itDof).getValue(step);
%                 
%             end
%         
%             end
%         end          
%   
% 
% 
%         Shape Function and Derivatives     
%     
%         function [N, Nf, B_b, B, J_det] = computeShapeFunction(mixed2d4n,xi,eta)
% 
%             Nf = [(1-xi)*(1-eta)/4    (1+xi)*(1-eta)/4    (1+xi)*(1+eta)/4    (1-xi)*(1+eta)/4];
%             
% 
%             N_Diff_Par = [-(1-eta)/4    (1-eta)/4   (1+eta)/4   -(1+eta)/4
%                 -(1-xi)/4     -(1+xi)/4   (1+xi)/4    (1-xi)/4];
%             
%             dN_xi  = [-(1-eta)/4    (1-eta)/4   (1+eta)/4   -(1+eta)/4]; 
%             dN_eta = [-(1-xi)/4     -(1+xi)/4   (1+xi)/4    (1-xi)/4]; 
%             
%             N = sparse(2,8);
%             N(1,1:2:end) = Nf(:);
%             N(2,2:2:end) = Nf(:);
%             
%             Coordinates of the nodes forming one element
%             ele_coords = zeros(4,2);
%             for i=1:4
%                 ele_coords(i,1) = mixed2d4n.nodeArray(i).getX;
%                 ele_coords(i,2) = mixed2d4n.nodeArray(i).getY;
%             end
%             
%             Jacobian
%             J = [dN_xi;dN_eta] * ele_coords;
%             J_det = det(J);
%             
%             Calculation of B-Matrix for Solid Phase     
%             B=J\N_Diff_Par;
%             Bx=B(1,1:4);
%             By=B(2,1:4);
%             B_b = [Bx(1),0,Bx(2),0,Bx(3),0,Bx(4),0;
%             0,By(1),0,By(2),0,By(3),0,By(4);
%             By(1),Bx(1),By(2),Bx(2),By(3),Bx(3),By(4),Bx(4)];
%        
%         end



        function [totalElementStiffnessMatrix] = computeLocalStiffnessMatrix(mixed2d4n)
     
            ETA_S = mixed2d4n.getPropertyValue('ETA_S');
            LAMBDA = mixed2d4n.getPropertyValue('LAMBDA_S');
            MUE = mixed2d4n.getPropertyValue('MUE_S');
            OMEGA = mixed2d4n.getPropertyValue('OMEGA');
            POROSITY = mixed2d4n.getPropertyValue('POROSITY');
            HEAT_CAPACITY_RATIO= mixed2d4n.getPropertyValue('HEAT_CAPACITY_RATIO_F');    
            PRESSURE_0 = mixed2d4n.getPropertyValue('PRESSURE_0_F');
            ETA_F = mixed2d4n.getPropertyValue('ETA_F');
            PRANDL_NUMBER = mixed2d4n.getPropertyValue('PRANDL_NUMBER_F');
            THERMAL_CHARACT_LENGTH = mixed2d4n.getPropertyValue('THERMAL_CHARACT_LENGTH');
            DENSITY_F = mixed2d4n.getPropertyValue('DENSITY_F');
           
            p = mixed2d4n.getPropertyValue('NUMBER_GAUSS_POINT');   
 
            STATIC_FLOW_RESISTIVITY = mixed2d4n.getPropertyValue('STATIC_FLOW_RESISTIVITY');
            TORTUOSITY = mixed2d4n.getPropertyValue('TORTUOSITY');
            VISCOUS_CHARACT_LENGTH = mixed2d4n.getPropertyValue('VISCOUS_CHARACT_LENGTH');
    
            % Determine relevant densities
            VISCOUS_DRAG = STATIC_FLOW_RESISTIVITY * POROSITY^2 * (1 + ...
            4 * 1i * TORTUOSITY^2 * ETA_F * DENSITY_F / (STATIC_FLOW_RESISTIVITY^2 * VISCOUS_CHARACT_LENGTH^2 * POROSITY^2));

            % Inertial coupling term, related to the tortuosity, increases the fluid density  to model the motion
            % of fluid particles vibrating around the structural frame:
            DENSITY_A = POROSITY * DENSITY_F * (TORTUOSITY -1);
            
            % Equivalent densities for expressing the elastodynamic coupled equations in condensed form:
            EQ_DENSITY_SF = (-1) * DENSITY_A + 1i * VISCOUS_DRAG / OMEGA;
            EQ_DENSITY_F = POROSITY * DENSITY_F - EQ_DENSITY_SF;
         
            % Determine Coefficients for EMat
            % Hysteretic proportional damping model:
            LAME_COEFF = (1+1i*ETA_S)*LAMBDA;
            SHEAR_MODULUS = (1+1i*ETA_S)*MUE;
                       
            % Calculate Bulk-Modulus K_f:
            K_f = (HEAT_CAPACITY_RATIO*PRESSURE_0)/(HEAT_CAPACITY_RATIO-(HEAT_CAPACITY_RATIO-1)*...
                (1+(8*ETA_F)/(1i*OMEGA*PRANDL_NUMBER*THERMAL_CHARACT_LENGTH^2*DENSITY_F)*...
                (1+(1i*OMEGA*PRANDL_NUMBER*THERMAL_CHARACT_LENGTH^2*DENSITY_F)/(16*ETA_F))^(0.5))^(-1));
      
            R = POROSITY * K_f;
            Q = (1 - POROSITY) * K_f;
            GAMMA = POROSITY * ((EQ_DENSITY_SF /EQ_DENSITY_F) - (Q/R)) ;
            
            % Stress-strain relation for solid phase
            D = [LAME_COEFF + 2 * SHEAR_MODULUS,LAME_COEFF,0,; ...
                LAME_COEFF, LAME_COEFF + 2 * SHEAR_MODULUS,0,; ...
                0,0,SHEAR_MODULUS];
            
            K_Matrix = zeros(8,8);
            C_Matrix = zeros(8,4);
            H_Matrix = zeros(4,4);
            
            [w,g]=returnGaussPoint(p);
            
    
            for zeta=1:p
                for eta=1:p
                    % N = 2x8 , Nf = 1x4, Be = 2x4, B = 3x8
                    [N, ~, Be, B, Jdet] = computeShapeFunction(mixed2d4n,g(zeta),g(eta));
                    % K = B'*B = 8x3 * 3x8 = 8x8
                    K_Matrix = K_Matrix + (w(zeta) * w(eta) * Jdet * B' * D * B);
                    % C = N'*Be = 8x2 * 2x4 = 8x4
                    C_Matrix = C_Matrix + (w(zeta) * w(eta) * Jdet * GAMMA * N' * Be);
                    % H = Be'*Be = 4x2 * 2x4 = 4x4
                    H_Matrix = H_Matrix + (w(zeta) * w(eta) * (POROSITY^2/EQ_DENSITY_F) * Jdet * Be' * Be);
            
                  
                end
            end
            
            empty= zeros(4,8);
            totalElementStiffnessMatrix = [K_Matrix, - C_Matrix; empty, H_Matrix];

        end    


        function [totalElementMassMatrix] = computeLocalMassMatrix(mixed2d4n)

            OMEGA = mixed2d4n.getPropertyValue('OMEGA');
            POROSITY = mixed2d4n.getPropertyValue('POROSITY');
            
            HEAT_CAPACITY_RATIO= mixed2d4n.getPropertyValue('HEAT_CAPACITY_RATIO_F');
            PRESSURE_0 = mixed2d4n.getPropertyValue('PRESSURE_0_F');
            ETA_F = mixed2d4n.getPropertyValue('ETA_F');
            PRANDL_NUMBER = mixed2d4n.getPropertyValue('PRANDL_NUMBER_F');
            THERMAL_CHARACT_LENGTH = mixed2d4n.getPropertyValue('THERMAL_CHARACT_LENGTH');
            DENSITY_F = mixed2d4n.getPropertyValue('DENSITY_F');
            DENSITY_S = mixed2d4n.getPropertyValue('DENSITY_S');
            p = mixed2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            
            STATIC_FLOW_RESISTIVITY = mixed2d4n.getPropertyValue('STATIC_FLOW_RESISTIVITY');
            TORTUOSITY = mixed2d4n.getPropertyValue('TORTUOSITY');
            VISCOUS_CHARACT_LENGTH = mixed2d4n.getPropertyValue('VISCOUS_CHARACT_LENGTH');
            
            % Determine relevant densities
            VISCOUS_DRAG = STATIC_FLOW_RESISTIVITY * POROSITY^2 * (1 + ...
                4 * 1i * TORTUOSITY^2 * ETA_F * DENSITY_F / (STATIC_FLOW_RESISTIVITY^2 * VISCOUS_CHARACT_LENGTH^2 * POROSITY^2));
            
            % Inertial coupling term, related to the tortuosity, increases the fluid density  to model the motion
            % of fluid particles vibrating around the structural frame:
            DENSITY_A = POROSITY * DENSITY_F * (TORTUOSITY -1);
            
            % Equivalent densities for expressing the elastodynamic coupled equations in condensed form:
            EQ_DENSITY_SF = (-1) * DENSITY_A + 1i * VISCOUS_DRAG / OMEGA;
            EQ_DENSITY_F = POROSITY * DENSITY_F - EQ_DENSITY_SF;
            EQ_DENSITY_S = (1 - POROSITY) * DENSITY_S - EQ_DENSITY_SF;
            EQ_DENSITY_TOTAL = EQ_DENSITY_S - EQ_DENSITY_SF^2/EQ_DENSITY_F;
            
            % Calculate Bulk-Modulus K_f:
            K_f = (HEAT_CAPACITY_RATIO*PRESSURE_0)/(HEAT_CAPACITY_RATIO-(HEAT_CAPACITY_RATIO-1)*...
                (1+(8*ETA_F)/(1i*OMEGA*PRANDL_NUMBER*THERMAL_CHARACT_LENGTH^2*DENSITY_F)*...
                (1+(1i*OMEGA*PRANDL_NUMBER*THERMAL_CHARACT_LENGTH^2*DENSITY_F)/(16*ETA_F))^(0.5))^(-1));
            
            R = POROSITY * K_f;
            Q = (1 - POROSITY) * K_f;
            GAMMA = POROSITY * ((EQ_DENSITY_SF /EQ_DENSITY_F) - (Q/R)) ;
            
            M_Matrix = zeros(8,8);
            Q_Matrix = zeros(4,4);
            C_Matrix = zeros(8,4);
            
            [w,g]=returnGaussPoint(p);
            
            
            for zeta=1:p
                for eta=1:p
                    % Nmat = 2x8 , Nf = 1x4, Be = 2x4, B = 3x8
                    [Nmat, Nf, Be, ~, Jdet] = computeShapeFunction(mixed2d4n,g(zeta),g(eta));
                    % M = N'*N = 8x2 * 2x8 = 8x8
                    M_Matrix = M_Matrix + (w(zeta) * w(eta) * Jdet * EQ_DENSITY_TOTAL * Nmat' * Nmat);
                    % Q = Nf'*Nf = 4x1 * 1x4 = 4x4
                    Q_Matrix = Q_Matrix + (w(zeta) * w(eta) * Jdet * (POROSITY^2/R) * Nf' * (Nf));
                    % C = N'*B = 8x2 * 2x4 = 8x4
                    C_Matrix = C_Matrix + (w(zeta) * w(eta) * Jdet * GAMMA * Nmat' * Be);
                    
                end
            end
            
            empty = zeros(8,4);
%             totalElementMassMatrix = [- OMEGA^2 * M_Matrix, empty'; - OMEGA^2 * C_Matrix',  - OMEGA^2 * Q_Matrix];
            totalElementMassMatrix = [M_Matrix, empty; C_Matrix', Q_Matrix]; 
            
        end
        
        
        % Define possible DOF in the node-array:
        function dofs = getDofList(element)
            dofs([1 3 5 7])     = element.nodeArray.getDof('DISPLACEMENT_SOLID_X');
            dofs([2 4 6 8])     = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y');
            dofs([9 10 11 12])  = element.nodeArray.getDof('PRESSURE_FLUID');      
        end
        
        
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,12);
            vals([1 3 5 7])     = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_X',step);
            vals([2 4 6 8])     = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_Y',step);
            vals([9 10 11 12])     = element.nodeArray.getDofValue('PRESSURE _FLUID',step);
        end
        
       
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,12);
            [~, vals([1 3 5 7]), ~]     = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, vals([2 4 6 8]), ~]     = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, vals([9 10 11 12]), ~]     = element.nodeArray.getDof('PRESSURE_FLUID').getAllValues(step);
          
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,12);       
            [~, ~, vals([1 3 5 7])]     = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, ~, vals([2 4 6 8])]     = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, ~, vals([9 10 11 12])]     = element.nodeArray.getDof('PRESSURE_FLUID').getAllValues(step);
          
        end
        
        function [totalElementDampingMatrix] = computeLocalDampingMatrix(mixed2d4n)
            totalElementDampingMatrix = zeros(12,12);
        end

        function F = computeLocalForceVector(mixed2d4n)
            F = zeros(1,12);
        end
        
    end
end