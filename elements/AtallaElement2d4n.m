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
            % Initialize properties:
            eta_s = mixed2d4n.getPropertyValue('ETA_S');
            lambda = mixed2d4n.getPropertyValue('LAMBDA_S');
            mu = mixed2d4n.getPropertyValue('MUE_S');
            phi = mixed2d4n.getPropertyValue('POROSITY');
            gamma = mixed2d4n.getPropertyValue('HEAT_CAPACITY_RATIO_F');    
            P_0 = mixed2d4n.getPropertyValue('PRESSURE_0_F');
            eta_f = mixed2d4n.getPropertyValue('ETA_F');
            Pr = mixed2d4n.getPropertyValue('PRANDL_NUMBER_F');
            Lambda_t = mixed2d4n.getPropertyValue('THERMAL_CHARACT_LENGTH');
            rho_f = mixed2d4n.getPropertyValue('DENSITY_F');
            sigma = mixed2d4n.getPropertyValue('STATIC_FLOW_RESISTIVITY');
            alpha_inf = mixed2d4n.getPropertyValue('TORTUOSITY');
            Lambda_v = mixed2d4n.getPropertyValue('VISCOUS_CHARACT_LENGTH');
            omega = mixed2d4n.getPropertyValue('OMEGA');
            p = mixed2d4n.getPropertyValue('NUMBER_GAUSS_POINT'); 
            % Calculate flow resistivity of air particles in the pores:
            G_J = (1+(4*1i*omega*alpha_inf^2*eta_f*rho_f)/...
            (sigma^2*Lambda_v^2*phi^2))^(0.5); 
            % Calculate viscous drag:
            bF = sigma * phi^2 * G_J;
            % Calculate inertial coupling term:
            Rho_a = phi * rho_f * (alpha_inf -1);
            % Calculate equivalent densities:
            Rho_sf = (-1) * Rho_a + 1i * bF / omega;
            Rho_ff = phi * rho_f - Rho_sf;
            % Calculate Lame coefficients:
            Lame_alpha = (1+1i*eta_s)*lambda;
            Lame_mu = (1+1i*eta_s)*mu;
            % Calculate fluid Bulk-Modulus:
            K_f = (gamma*P_0)/(gamma-(gamma-1)*...
                (1+(8*eta_f)/(1i*omega*Pr*Lambda_t^2*rho_f)*...
                (1+(1i*omega*Pr*Lambda_t^2*rho_f)/(16*eta_f))^(0.5))^(-1));
            % Calculate Biot coefficients:
            R = phi * K_f;
            Q = (1 - phi) * K_f;
            % Calculate Atallas gamma coefficient:
            GAMMA = phi * ((Rho_sf /Rho_ff) - (Q/R)) ;
            % Calculate solid phase elasticity tensor:
            D = [Lame_alpha + 2 * Lame_mu,Lame_alpha,0,; ...
                Lame_alpha, Lame_alpha + 2 * Lame_mu,0,; ...
                0,0,Lame_mu];
            % Generate empty partial stiffness-matrices:
            K_M = zeros(8,8);
            C_M = zeros(8,4);
            H_M = zeros(4,4);
            % Get gauss integration factors:
            [w,g]=returnGaussPoint(p);
            % Execute Gauss-Integration 
            for zeta=1:p
                for eta=1:p
                    % N = 2x8 , Nf = 1x4, Be = 2x4, B = 3x8
                    [N, ~, Be, B, Jdet] = computeShapeFunction(mixed2d4n,g(zeta),g(eta));
                    K_M = K_M + (w(zeta) * w(eta) * Jdet * B' * D * B);
                    C_M = C_M + (w(zeta) * w(eta) * Jdet * GAMMA * N' * Be);
                    H_M = H_M + (w(zeta) * w(eta) * (phi^2/Rho_ff) * Jdet * Be' * Be);
                end
            end
            % Assemble equivalent element stiffness matrix:
            empty= zeros(4,8);
            totalElementStiffnessMatrix = [K_M, - C_M; empty, H_M];
        end    


        function [totalElementMassMatrix] = computeLocalMassMatrix(mixed2d4n)
            % Initialize properties:
            phi = mixed2d4n.getPropertyValue('POROSITY');
            gamma= mixed2d4n.getPropertyValue('HEAT_CAPACITY_RATIO_F');
            P_0 = mixed2d4n.getPropertyValue('PRESSURE_0_F');
            eta_f = mixed2d4n.getPropertyValue('ETA_F');
            Pr = mixed2d4n.getPropertyValue('PRANDL_NUMBER_F');
            Lambda_t = mixed2d4n.getPropertyValue('THERMAL_CHARACT_LENGTH');
            rho_f = mixed2d4n.getPropertyValue('DENSITY_F');
            rho_s = mixed2d4n.getPropertyValue('DENSITY_S');
            Xi = mixed2d4n.getPropertyValue('STATIC_FLOW_RESISTIVITY');
            alpha_inf = mixed2d4n.getPropertyValue('TORTUOSITY');
            Lambda_v = mixed2d4n.getPropertyValue('VISCOUS_CHARACT_LENGTH');
            omega = mixed2d4n.getPropertyValue('OMEGA');            
            p = mixed2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            % Calculate flow resistivity of air particles in the pores:
            G_J = (1+(4*1i*omega*alpha_inf^2*eta_f*rho_f)/...
            (Xi^2*Lambda_v^2*phi^2))^(0.5); 
            % Calculate viscous drag:
            bF = Xi * phi^2 * G_J;
            % Calculate inertial coupling term:
            Rho_a = phi * rho_f * (alpha_inf -1);
            % Calculate equivalent densities:
            Rho_sf = (-1) * Rho_a + 1i * bF / omega;
            Rho_ff = phi * rho_f - Rho_sf;
            Rho_ss = (1 - phi) * rho_s - Rho_sf;
            Rho_total = Rho_ss - Rho_sf^2/Rho_ff;
            % Calculate fluid bulk modulus:
            K_f = (gamma*P_0)/(gamma-(gamma-1)*...
                (1+(8*eta_f)/(1i*omega*Pr*Lambda_t^2*rho_f)*...
                (1+(1i*omega*Pr*Lambda_t^2*rho_f)/(16*eta_f))^(0.5))^(-1));
            % Calculate Biot's coefficients:
            R = phi * K_f;
            Q = (1 - phi) * K_f;
            % Calculate Atallas gamma coefficient:
            GAMMA = phi * ((Rho_sf /Rho_ff) - (Q/R)) ;
            % Generate empty partial stiffness-matrices:
            M_M = zeros(8,8);
            Q_M = zeros(4,4);
            C_M = zeros(8,4);
            % Get gauss integration factors:
            [w,g]=returnGaussPoint(p);
            % Execute full Gauss-Integration 
            for zeta=1:p
                for eta=1:p
                    % Nmat = 2x8 , Nf = 1x4, Be = 2x4, B = 3x8
                    [Nmat, Nf, Be, ~, Jdet] = computeShapeFunction(mixed2d4n,g(zeta),g(eta));
                    M_M = M_M + (w(zeta) * w(eta) * Jdet * Rho_total * Nmat' * Nmat);
                    Q_M = Q_M + (w(zeta) * w(eta) * Jdet * (phi^2/R) * Nf' * (Nf));
                    C_M = C_M + (w(zeta) * w(eta) * Jdet * GAMMA * Nmat' * Be);
                end
            end
            % Assemble equivalent element mass matrix:
            empty = zeros(8,4);
            totalElementMassMatrix = [M_M, empty; C_M', Q_M]; 
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