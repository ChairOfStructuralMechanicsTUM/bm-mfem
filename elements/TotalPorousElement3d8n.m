classdef TotalPorousElement3d8n < PorousElement3d8n
    %   POROUSELEMENT3D8N Summary of this class goes here
    %   Detailed explanation goes here
    
    %   DOF: Solid-Displacement & Total-Displacement (see: DAZEL 2009)
    
    properties (Access = private)

    end
    
    methods
        % constructor
        function totalporousElement3d8n = TotalPorousElement3d8n(id, nodeArray)
            
            requiredPropertyNames = cellstr(["DENSITY_SOLID", "LAMBDA_SOLID", "MUE_SOLID", "DAMPING_SOLID", ...
                "DENSITY_FLUID", "VISCOSITY_FLUID", "STANDARD_PRESSURE_FLUID", "HEAT_CAPACITY_FLUID", "PRANDTL_NUMBER_FLUID", ...
                "POROSITY", "TORTUOSITY", "FLOW_RESISTIVITY", "VISCOUS_LENGHT", "THERMAL_LENGTH", ...
                "FREQUENCY", "NUMBER_GAUSS_POINT"]);
            
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 8 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            % call the super class constructor
            totalporousElement3d8n@PorousElement3d8n(super_args{:});
            totalporousElement3d8n.dofNames = cellstr(['DISPLACEMENT_SOLID_X'; 'DISPLACEMENT_SOLID_Y'; 'DISPLACEMENT_SOLID_Z'; ...
                'DISPLACEMENT_TOTAL_X'; 'DISPLACEMENT_TOTAL_Y'; 'DISPLACEMENT_TOTAL_Z']);
            
        end
        
        %Initialization
        function initialize(totalporousElement3d8n)
            totalporousElement3d8n.lengthX = computeLength(totalporousElement3d8n.nodeArray(1).getCoords, ...
                totalporousElement3d8n.nodeArray(2).getCoords);
            
            totalporousElement3d8n.lengthY = computeLength(totalporousElement3d8n.nodeArray(1).getCoords, ...
                totalporousElement3d8n.nodeArray(4).getCoords);
            
            totalporousElement3d8n.lengthY = computeLength(totalporousElement3d8n.nodeArray(1).getCoords, ...
                totalporousElement3d8n.nodeArray(5).getCoords);
            
            checkConvexity(totalporousElement3d8n);
        end
        
        % member functions
                function stiffnessMatrix = computeLocalStiffnessMatrix(totalporousElement3d8n)
            % Read Material Properties
            LAMBDA_SOLID = totalporousElement3d8n.getPropertyValue('LAMBDA_SOLID');
            MUE_SOLID = totalporousElement3d8n.getPropertyValue('MUE_SOLID');
            DAMPING_SOLID = totalporousElement3d8n.getPropertyValue('DAMPING_SOLID');
            
            DENSITY_FLUID = totalporousElement3d8n.getPropertyValue('DENSITY_FLUID');
            VISCOSITY_FLUID = totalporousElement3d8n.getPropertyValue('VISCOSITY_FLUID');
            STANDARD_PRESSURE_FLUID = totalporousElement3d8n.getPropertyValue('STANDARD_PRESSURE_FLUID');
            HEAT_CAPACITY_FLUID = totalporousElement3d8n.getPropertyValue('HEAT_CAPACITY_FLUID');
            PRANDTL_NUMBER_FLUID = totalporousElement3d8n.getPropertyValue('PRANDTL_NUMBER_FLUID');
            
            THERMAL_LENGTH = totalporousElement3d8n.getPropertyValue('THERMAL_LENGTH');
            
            omega = totalporousElement3d8n.getPropertyValue('FREQUENCY');
            p = totalporousElement3d8n.getPropertyValue('NUMBER_GAUSS_POINT');
            
            % Determine Coefficients needed
            K_EQ = HEAT_CAPACITY_FLUID * STANDARD_PRESSURE_FLUID / (HEAT_CAPACITY_FLUID - (HEAT_CAPACITY_FLUID -1)*...
                (1 + (8 * VISCOSITY_FLUID / ( 1i * THERMAL_LENGTH * PRANDTL_NUMBER_FLUID * omega * DENSITY_FLUID))*...
                (1 + (1i * DENSITY_FLUID * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2)/(16 * VISCOSITY_FLUID))^0.5)); 
            N = MUE_SOLID * ( 1 + 1i * DAMPING_SOLID);
            A = LAMBDA_SOLID * ( 1 + 1i * DAMPING_SOLID);
            P = A + 2 * N;

            D_0 = [1,1,1,0,0,0; ...
                1,1,1,0,0,0; ...
                1,1,1,0,0,0; ...
                0,0,0,0,0,0; ...
                0,0,0,0,0,0; ...
                0,0,0,0,0,0];
            D_1 = [0,-2,-2,0,0,0; ...
                -2,0,-2,0,0,0; ...
                -2,-2,0,0,0,0; ...
                0,0,0,1,0,0;...
                0,0,0,0,1,0; ...
                0,0,0,0,0,1];
            
            K_0 = zeros(24,24);
            K_1 = zeros(24,24);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                zeta=g(i);
                for j=1:p
                    eta=g(j);
                    for k=1:p
                        mue=g(k);
                        [~, ~, Be, ~, Jdet] = computeShapeFunction(totalporousElement3d8n, zeta, eta, mue);
                        K_0=K_0+(w(i)*w(j)*w(k)*Jdet*transpose(Be)*(D_0*Be));
                        K_1=K_1+(w(i)*w(j)*w(k)*Jdet*transpose(Be)*(D_1*Be));
                    end
                end
            end
            
            stiffnessMatrix = [P * K_0 + N * K_1, zeros(24,24); zeros(24,24), K_EQ * K_0];
        end
        
        function massMatrix = computeLocalMassMatrix(totalporousElement3d8n)
            % Read Material Properties
            DENSITY_SOLID = totalporousElement3d8n.getPropertyValue('DENSITY_SOLID');
            
            DENSITY_FLUID = totalporousElement3d8n.getPropertyValue('DENSITY_FLUID');
            VISCOSITY_FLUID = totalporousElement3d8n.getPropertyValue('VISCOSITY_FLUID');
            
            POROSITY = totalporousElement3d8n.getPropertyValue('POROSITY');
            TORTUOSITY = totalporousElement3d8n.getPropertyValue('TORTUOSITY');
            FLOW_RESISTIVITY = totalporousElement3d8n.getPropertyValue('FLOW_RESISTIVITY');
            VISCOUS_LENGHT = totalporousElement3d8n.getPropertyValue('VISCOUS_LENGHT');
            
            omega = totalporousElement3d8n.getPropertyValue('FREQUENCY');
            
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
          
            ElementmassMatrix=zeros(24,24);
            
            p = totalporousElement3d8n.getPropertyValue('NUMBER_GAUSS_POINT');
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                zeta=g(i);
                for j=1:p
                    eta=g(j);
                    for k=1:p
                        mue=g(k);
                        [N, ~, ~, ~, Jdet] = computeShapeFunction(totalporousElement3d8n, zeta, eta, mue);
                        ElementmassMatrix = ElementmassMatrix + (w(i) * w(j) * w(k)* transpose(N) * N * Jdet);
                    end
                end
            end
            
            massMatrix = [DENSITY_S * ElementmassMatrix, GAMMA * DENSITY_EQ * ElementmassMatrix; ...
                GAMMA * DENSITY_EQ * ElementmassMatrix, DENSITY_EQ * ElementmassMatrix];
            
        end
        
        function dampingMatrix = computeLocalDampingMatrix(totalporousElement3d8n)
            dampingMatrix = zeros(48,48);
        end
        
        function dofs = getDofList(element)
            dofs([1 4 7 10 13 16 19 22]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_X');
            dofs([2 5 8 11 14 17 20 23]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y');
            dofs([3 6 9 12 15 18 21 24]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_Z');
            dofs([25 28 31 34 37 40 43 46]) = element.nodeArray.getDof('DISPLACEMENT_TOTAL_X');
            dofs([26 29 32 35 38 41 44 47]) = element.nodeArray.getDof('DISPLACEMENT_TOTAL_Y');
            dofs([27 30 33 36 39 42 45 48]) = element.nodeArray.getDof('DISPLACEMENT_TOTAL_Z');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,48);
            
            vals([1 4 7 10 13 16 19 22]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_X',step);
            vals([2 5 8 11 14 17 20 23]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_Y',step);
            vals([3 6 9 12 15 18 21 24]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_Z',step);
            vals([25 28 31 34 37 40 43 46]) = element.nodeArray.getDofValue('DISPLACEMENT_TOTAL_X',step);
            vals([26 29 32 35 38 41 44 47]) = element.nodeArray.getDofValue('DISPLACEMENT_TOTAL_Y',step);
            vals([27 30 33 36 39 42 45 48]) = element.nodeArray.getDofValue('DISPLACEMENT_TOTAL_Z',step);
        end
        
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,48);
            
            [~, vals([1 4 7 10 13 16 19 22]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, vals([2 5 8 11 14 17 20 23]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, vals([3 6 9 12 15 18 21 24]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Z').getAllValues(step);
            [~, vals([25 28 31 34 37 40 43 46]), ~] = element.nodeArray.getDof('DISPLACEMENT_TOTAL_X').getAllValues(step);
            [~, vals([26 29 32 35 38 41 44 47]), ~] = element.nodeArray.getDof('DISPLACEMENT_TOTAL_Y').getAllValues(step);
            [~, vals([27 30 33 36 39 42 45 48]), ~] = element.nodeArray.getDof('DISPLACEMENT_TOTAL_Z').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,48);
            
            [~, ~, vals([1 4 7 10 13 16 19 22])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, ~, vals([2 5 8 11 14 17 20 23])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, ~, vals([3 6 9 12 15 18 21 24])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Z').getAllValues(step);
            [~, ~, vals([25 28 31 34 37 40 43 46])] = element.nodeArray.getDof('DISPLACEMENT_TOTAL_X').getAllValues(step);
            [~, ~, vals([26 29 32 35 38 41 44 47])] = element.nodeArray.getDof('DISPLACEMENT_TOTAL_Y').getAllValues(step);
            [~, ~, vals([27 30 33 36 39 42 45 48])] = element.nodeArray.getDof('DISPLACEMENT_TOTAL_Z').getAllValues(step);
        end
        
        function F = computeLocalForceVector(totalporousElement3d8n)
            F = zeros(1,48);
        end
        
    end
    
end
    

