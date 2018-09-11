classdef BiotAllardElement2d4n <  PorousElement2d4n
    
    properties (Access = private)
        
    end
    
    
    
    
    methods
        
        % Constructor
        function biot2d4n = BiotAllardElement2d4n(id,nodeArray)
            
            requiredPropertyNames = cellstr(["OMEGA","NUMBER_GAUSS_POINT", "LAMBDA_S", ...
                "MUE_S","ETA_S","DENSITY_S","DENSITY_F","ETA_F","PRESSURE_0_F","HEAT_CAPACITY_RATIO_F",...
                "PRANDL_NUMBER_F","POROSITY","TORTUOSITY","STATIC_FLOW_RESISTIVITY",...
                "VISCOUS_CHARACT_LENGTH","THERMAL_CHARACT_LENGTH"]);
            
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 4 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            % call the super class contructor
            biot2d4n@PorousElement2d4n(super_args{:});
            biot2d4n.dofNames = cellstr(["DISPLACEMENT_SOLID_X", "DISPLACEMENT_SOLID_Y", ...
                "DISPLACEMENT_FLUID_X", "DISPLACEMENT_FLUID_Y"]);
        end
        
        
        %Initialization
        function initialize(biot2d4n)
            biot2d4n.lengthX = computeLength(biot2d4n.nodeArray(1).getCoords, ...
                biot2d4n.nodeArray(2).getCoords);
            biot2d4n.lengthY = computeLength(biot2d4n.nodeArray(1).getCoords, ...
                biot2d4n.nodeArray(4).getCoords);
            
            checkConvexity(biot2d4n);
        end
        
        
        % function responseDoF = getResponseDofArray(element, step)
        %
        %     responseDoF = zeros(12,1);
        %     for itNodes = 1:1:4
        %         nodalDof = element.nodeArray(itNodes).getDofArray;
        %         nodalDof = nodalDof.';
        %
        %         for itDof = 3:(-1):1
        %             responseDoF(3*itNodes-(itDof-1),1) = nodalDof(4-itDof).getValue(step);
        %
        %         end
        %
        %     end
        % end
        %
        %
        %     % Shape Function and Derivatives
        % function [N_mat, N, B_b, J] = computeShapeFunction(biot2d4n,xi,eta)
        %
        %             N = [(1-xi)*(1-eta)/4    (1+xi)*(1-eta)/4    (1+xi)*(1+eta)/4    (1-xi)*(1+eta)/4];
        %
        %             N_Diff_Par = [-(1-eta)/4    (1-eta)/4   (1+eta)/4   -(1+eta)/4
        %                 -(1-xi)/4     -(1+xi)/4   (1+xi)/4    (1-xi)/4];
        %
        %             N_mat = sparse(2,8);
        %             N_mat(1,1:2:end) = N(:);
        %             N_mat(2,2:2:end) = N(:);
        %
        %     % Coordinates of the nodes forming one element
        %             ele_coords = zeros(4,2);
        %             for i=1:4
        %                 ele_coords(i,1) = biot2d4n.nodeArray(i).getX;
        %                 ele_coords(i,2) = biot2d4n.nodeArray(i).getY;
        %             end
        %
        %
        %     % Jacobian
        %             J = N_Diff_Par * ele_coords;
        %
        %     % Calculation of B-Matrix
        %             B=J\N_Diff_Par;%/Jdet;
        %             Bx=B(1,1:4);
        %             By=B(2,1:4);
        %
        %     % Calculation of B-Matrix
        %             %B_b = sparse(3,8);
        %             B_b=[Bx(1),0,Bx(2),0,Bx(3),0,Bx(4),0;
        %                 0,By(1),0,By(2),0,By(3),0,By(4);
        %                 By(1),Bx(1),By(2),Bx(2),By(3),Bx(3),By(4),Bx(4)];
        %
        %
        % end
        
        function [totalElementStiffnessMatrix] = computeLocalStiffnessMatrix(biot2d4n)
            % Initialize properties:
            eta_s = biot2d4n.getPropertyValue('ETA_S');
            lambda = biot2d4n.getPropertyValue('LAMBDA_S');
            mu = biot2d4n.getPropertyValue('MUE_S');
            omega = biot2d4n.getPropertyValue('OMEGA');
            phi = biot2d4n.getPropertyValue('POROSITY');
            gamma = biot2d4n.getPropertyValue('HEAT_CAPACITY_RATIO_F');
            P_0 = biot2d4n.getPropertyValue('PRESSURE_0_F');
            eta_f = biot2d4n.getPropertyValue('ETA_F');
            Pr = biot2d4n.getPropertyValue('PRANDL_NUMBER_F');
            Lambda_t = biot2d4n.getPropertyValue('THERMAL_CHARACT_LENGTH');
            rho_f = biot2d4n.getPropertyValue('DENSITY_F');
            p = biot2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            % Calculate Bulk-Modulus K_f:
            K_f = (gamma*P_0)/(gamma-(gamma-1)*...
                (1+(8*eta_f)/(1i*omega*Pr*Lambda_t^2*rho_f)*...
                (1+(1i*omega*Pr*Lambda_t^2*rho_f)/(16*eta_f))^(0.5))^(-1)); 
            % Calculate Lame coefficients:
            Lame_lambda = (1+1i*eta_s)*lambda;
            Lame_mu = (1+1i*eta_s)*mu;
            % Calculate Biot's elasticity coefficients:
            R = phi*K_f;
            Q = (1-phi)*K_f;
            A = Lame_lambda+(1-phi)^2/phi*K_f;
            % Elasticity Tensors:
            D_s = [A+2*Lame_mu,A,0;A,A+2*Lame_mu,0;0,0,Lame_mu];
            D_sf = [Q,Q,0;Q,Q,0;0,0,0];
            D_f = [R,R,0;R,R,0;0,0,0];
            % Get gauss integration factors:
            [w,g]=returnGaussPoint(p);
            % Generate empty partial stiffness-matrices:
            stiffnessMatrix_s=sparse(8,8);
            stiffnessMatrix_f=sparse(8,8);
            stiffnessMatrix_sf = sparse(8,8);
            % Execute gauss integration to obtain partial
            % stiffness-matrices for solid, fluid and coupling phases:
            for xi = 1 : p
                for eta = 1 : p
                    [~, ~, ~, B, Jdet] = computeShapeFunction(biot2d4n,g(xi),g(eta));
                    
                    stiffnessMatrix_s = stiffnessMatrix_s +...
                        B' * D_s * B * Jdet * w(xi) * w(eta);
                    
                    stiffnessMatrix_f = stiffnessMatrix_f +             ...
                        B' * D_f * B * Jdet * w(xi) * w(eta);
                    
                    stiffnessMatrix_sf = stiffnessMatrix_sf +             ...
                        B' * D_sf * B * Jdet * w(xi) * w(eta);
                    
                end
            end
            % Assemble total element stiffness matrix:
            totalElementStiffnessMatrix = [stiffnessMatrix_s, stiffnessMatrix_sf;...
                stiffnessMatrix_sf,stiffnessMatrix_f];   
        end
        
        
        function [totalElementMassMatrix] = computeLocalMassMatrix(biot2d4n)
            % Initialize properties:
            sigma = biot2d4n.getPropertyValue('STATIC_FLOW_RESISTIVITY');
            phi = biot2d4n.getPropertyValue('POROSITY');
            omega = biot2d4n.getPropertyValue('OMEGA');
            alpha_inf = biot2d4n.getPropertyValue('TORTUOSITY');
            eta_f = biot2d4n.getPropertyValue('ETA_F');
            rho_f = biot2d4n.getPropertyValue('DENSITY_F');
            rho_s = biot2d4n.getPropertyValue('DENSITY_S');
            Lambda_v = biot2d4n.getPropertyValue('VISCOUS_CHARACT_LENGTH');
            p = biot2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            % Calculate flow resistivity of air particles in the pores:
            G_J = (1+(4*1i*omega*alpha_inf^2*eta_f*rho_f)/...
                (sigma^2*Lambda_v^2*phi^2))^(0.5);
            % Calculate viscous drag:
            bF = sigma*phi^2*G_J;
            % Calculate inertial coupling term:
            Rho_a = phi*rho_f*(alpha_inf-1);
            % Calculate equivalent densities:
            Rho_sf = -Rho_a+1i*(bF)/(omega);
            Rho_ss = (1-phi)*rho_s-Rho_sf;
            Rho_ff = phi*rho_f-Rho_sf;
            % Generate empty mass-matrix:
            M=zeros(8,8);
            % Get gauss integration factors:
            [w,g]=returnGaussPoint(p);
            % Execute gauss integration to obtain geometric-only mass matrix:
            for n=1:p
                xi=g(n);
                for m=1:p
                    eta=g(m);
                    [Nmat, ~, ~, ~, Jdet] = computeShapeFunction(biot2d4n,xi,eta);
                    M = M + (w(n)*w(m)*transpose(Nmat)*Nmat*Jdet); 
                end
            end
            % Calculate solid-, fluid- and coupling partial mass matrices:
            M_s = Rho_ss*M;
            M_f = Rho_ff*M;
            M_sf = Rho_sf*M;
            % Assemble element mass matrix:
            totalElementMassMatrix = [M_s,M_sf;M_sf,M_f];
        end
        
        % Define possible DOF in the node-array:
        function dofs = getDofList(element)
            dofs([1 3 5 7])     = element.nodeArray.getDof('DISPLACEMENT_SOLID_X');
            dofs([2 4 6 8])     = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y');
            dofs([9 11 13 15])  = element.nodeArray.getDof('DISPLACEMENT_FLUID_X');
            dofs([10 12 14 16]) = element.nodeArray.getDof('DISPLACEMENT_FLUID_Y');
        end
        
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,16);
            
            vals([1 3 5 7])     = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_X',step);
            vals([2 4 6 8])     = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_Y',step);
            vals([9 11 13 15])  = element.nodeArray.getDofValue('DISPLACEMENT_FLUID_X',step);
            vals([10 12 14 16]) = element.nodeArray.getDofValue('DISPLACEMENT_FLUID_Y',step);
        end
        
        
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,16);
            
            [~, vals([1 3 5 7]), ~]     = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, vals([2 4 6 8]), ~]     = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, vals([9 11 13 15]), ~]  = element.nodeArray.getDof('DISPLACEMENT_FLUID_X').getAllValues(step);
            [~, vals([10 12 14 16]), ~] = element.nodeArray.getDof('DISPLACEMENT_FLUID_Y').getAllValues(step);
            
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,16);
            
            [~, ~, vals([1 3 5 7])]     = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, ~, vals([2 4 6 8])]     = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, ~, vals([9 11 13 15])]  = element.nodeArray.getDof('DISPLACEMENT_FLUID_X').getAllValues(step);
            [~, ~, vals([10 12 14 16])] = element.nodeArray.getDof('DISPLACEMENT_FLUID_Y').getAllValues(step);
            
        end
        
        function [totalElementDampingMatrix] = computeLocalDampingMatrix(biot2d4n)
                    totalElementDampingMatrix = zeros(16,16);
        end
        
        
        function update(biot2d4n)
        end
        
        function F = computeLocalForceVector(biot2d4n)
            F = zeros(1,16);
        end
    end
end