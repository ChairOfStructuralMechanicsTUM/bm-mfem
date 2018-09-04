classdef DazelElement2d4n < PorousElement2d4n

    
    properties (Access = private)
    end
    
    methods
        
        % constructor
        function totalporousElement2d4n = DazelElement2d4n(id,nodeArray)
            
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
            totalporousElement2d4n@PorousElement2d4n(super_args{:});
            totalporousElement2d4n.dofNames = cellstr(['DISPLACEMENT_SOLID_X'; 'DISPLACEMENT_SOLID_Y'; ...
                'DISPLACEMENT_TOTAL_X'; 'DISPLACEMENT_TOTAL_Y']);
        end
        
        %Initialization
        function initialize(totalporousElement2d4n)
            totalporousElement2d4n.lengthX = computeLength(totalporousElement2d4n.nodeArray(1).getCoords, ...
                totalporousElement2d4n.nodeArray(2).getCoords);
            
            totalporousElement2d4n.lengthY = computeLength(totalporousElement2d4n.nodeArray(1).getCoords, ...
                totalporousElement2d4n.nodeArray(4).getCoords);
            
            checkConvexity(totalporousElement2d4n);
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(totalporousElement2d4n)
            % Read Material Properties
            LAMBDA_SOLID = totalporousElement2d4n.getPropertyValue('LAMBDA_S');
            MUE_SOLID = totalporousElement2d4n.getPropertyValue('MUE_S');
            DAMPING_SOLID = totalporousElement2d4n.getPropertyValue('ETA_S');
            
            DENSITY_FLUID = totalporousElement2d4n.getPropertyValue('DENSITY_F');
            VISCOSITY_FLUID = totalporousElement2d4n.getPropertyValue('ETA_F');
            STANDARD_PRESSURE_FLUID = totalporousElement2d4n.getPropertyValue('PRESSURE_0_F');
            HEAT_CAPACITY_FLUID = totalporousElement2d4n.getPropertyValue('HEAT_CAPACITY_RATIO_F');
            PRANDTL_NUMBER_FLUID = totalporousElement2d4n.getPropertyValue('PRANDL_NUMBER_F');
            
            THERMAL_LENGTH = totalporousElement2d4n.getPropertyValue('THERMAL_CHARACT_LENGTH');
            
            omega = totalporousElement2d4n.getPropertyValue('OMEGA');
            p = totalporousElement2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            
            % Determine Coefficients needed
            K_EQ = HEAT_CAPACITY_FLUID * STANDARD_PRESSURE_FLUID / (HEAT_CAPACITY_FLUID - (HEAT_CAPACITY_FLUID -1)*...
                (1 + (8 * VISCOSITY_FLUID / ( 1i * THERMAL_LENGTH * PRANDTL_NUMBER_FLUID * omega * DENSITY_FLUID))*...
                (1 + (1i * DENSITY_FLUID * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2)/(16 * VISCOSITY_FLUID))^0.5)); 
            N = MUE_SOLID * ( 1 + 1i * DAMPING_SOLID);
            A = LAMBDA_SOLID * ( 1 + 1i * DAMPING_SOLID);
            P = A + 2 * N;

            D_0 = [1, 1, 0; ...
                1, 1, 0; ...
                0,0,0];
            D_1 = [0,-2,0; ...
                -2,0,0; ...
                0,0,1];
            
            K_0 = zeros(8,8);
            K_1 = zeros(8,8);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [~, ~, ~, B, Jdet] = computeShapeFunction(totalporousElement2d4n, xi, eta);
                    K_0=K_0+(w(i)*w(j)*Jdet*transpose(B)*(D_0*B));
                    K_1=K_1+(w(i)*w(j)*Jdet*transpose(B)*(D_1*B));
                end
            end
            
            stiffnessMatrix = [P * K_0 + N * K_1, zeros(8,8); zeros(8,8), K_EQ * K_0];
        end
        
        function massMatrix = computeLocalMassMatrix(totalporousElement2d4n)
            % Read Material Properties
            DENSITY_SOLID = totalporousElement2d4n.getPropertyValue('DENSITY_S');
            
            DENSITY_FLUID = totalporousElement2d4n.getPropertyValue('DENSITY_F');
            VISCOSITY_FLUID = totalporousElement2d4n.getPropertyValue('ETA_F');
            
            POROSITY = totalporousElement2d4n.getPropertyValue('POROSITY');
            TORTUOSITY = totalporousElement2d4n.getPropertyValue('TORTUOSITY');
            FLOW_RESISTIVITY = totalporousElement2d4n.getPropertyValue('STATIC_FLOW_RESISTIVITY');
            VISCOUS_LENGHT = totalporousElement2d4n.getPropertyValue('VISCOUS_CHARACT_LENGTH');
            
            omega = totalporousElement2d4n.getPropertyValue('OMEGA');
            
            % Determine relevant densities
            DENSITY_1 = (1 - POROSITY) * DENSITY_SOLID;
            DENSITY_2 = POROSITY * DENSITY_FLUID;
            DENSITY_12 = - POROSITY * DENSITY_FLUID * (TORTUOSITY - 1);
            
            alpha = 1 - 1i * POROSITY * FLOW_RESISTIVITY / (TORTUOSITY * DENSITY_FLUID * omega) * ...
                (1- (4 * 1i * TORTUOSITY ^2 * VISCOSITY_FLUID * DENSITY_FLUID * omega)/(FLOW_RESISTIVITY * VISCOUS_LENGHT * POROSITY))^0.5;
            b = 1i * omega * POROSITY * DENSITY_FLUID * (alpha - TORTUOSITY);
            
            COMPLEX_DENSITY_12 = DENSITY_12 - b/(1i * omega);
            COMPLEX_DENSITY_22 = DENSITY_2 -COMPLEX_DENSITY_12; 
            COMPLEX_DENSITY = DENSITY_1 - COMPLEX_DENSITY_12 - COMPLEX_DENSITY_12^2/(DENSITY_2 - COMPLEX_DENSITY_12);
            
            DENSITY_EQ = COMPLEX_DENSITY_22/(POROSITY^2);
            GAMMA = POROSITY * (COMPLEX_DENSITY_12/COMPLEX_DENSITY_22 - (1 - POROSITY)/POROSITY);
            DENSITY_S = COMPLEX_DENSITY + GAMMA^2 * DENSITY_EQ;
          
            ElementmassMatrix=zeros(8,8);
            
            p = totalporousElement2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            [w,g]=returnGaussPoint(p);
            
            for i=1:p                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [N_mat, ~, ~, ~, detJ] = computeShapeFunction(totalporousElement2d4n, xi, eta);
                    ElementmassMatrix = ElementmassMatrix + (w(i) * w(j) * transpose(N_mat) * N_mat * detJ);
                end
            end
            
            massMatrix = [DENSITY_S * ElementmassMatrix, GAMMA * DENSITY_EQ * ElementmassMatrix; ...
                GAMMA * DENSITY_EQ * ElementmassMatrix, DENSITY_EQ * ElementmassMatrix];
            
        end
        
        function dampingMatrix = computeLocalDampingMatrix(totalporousElement2d4n)
            dampingMatrix = zeros(16,16);
        end
        
         function dofs = getDofList(element)
            dofs([1 3 5 7]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_X');
            dofs([2 4 6 8]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y');
            dofs([9 11 13 15]) = element.nodeArray.getDof('DISPLACEMENT_TOTAL_X');
            dofs([10 12 14 16]) = element.nodeArray.getDof('DISPLACEMENT_TOTAL_Y');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,16); 
            vals([1 3 5 7]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_X',step);
            vals([2 4 6 8]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_Y',step);
            vals([9 11 13 15]) = element.nodeArray.getDofValue('DISPLACEMENT_TOTAL_X',step);
            vals([10 12 14 16]) = element.nodeArray.getDofValue('DISPLACEMENT_TOTAL_Y',step);
        end
        
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,16);
            [~, vals([1 3 5 7]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, vals([2 4 6 8]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, vals([9 11 13 15]), ~] = element.nodeArray.getDof('DISPLACEMENT_TOTAL_X').getAllValues(step);
            [~, vals([10 12 14 16]), ~] = element.nodeArray.getDof('DISPLACEMENT_TOTAL_Y').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,16);
            [~, ~, vals([1 3 5 7])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, ~, vals([2 4 6 8])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, ~, vals([9 11 13 15])] = element.nodeArray.getDof('DISPLACEMENT_TOTAL_X').getAllValues(step);
            [~, ~, vals([10 12 14 16])] = element.nodeArray.getDof('DISPLACEMENT_TOTAL_Y').getAllValues(step);
        end
        
        function update(totalporousElement2d4n)
        end
        
        function F = computeLocalForceVector(totalporousElement2d4n)
            F = zeros(1,16);
        end
        
    end
    
end
    

