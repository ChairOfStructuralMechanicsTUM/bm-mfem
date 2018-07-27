classdef MixedPorousElement2d4n < PorousElement2d4n
    %   MixedPorousElement2d4n Summary of this class goes here
    %   Detailed explanation goes here
    
    %   DOF: Solid-Displacement & Pore-Pressure (see: ALLARD 1998)
    
    properties (Access = private)
    end
    
    methods
        
        % constructor
        function mixedporousElement2d4n = MixedPorousElement2d4n(id,nodeArray)
            
            requiredPropertyNames = cellstr(["DENSITY_SOLID", "LAMBDA_SOLID", "MUE_SOLID", "DAMPING_SOLID", ...
                "DENSITY_FLUID", "VISCOSITY_FLUID", "STANDARD_PRESSURE_FLUID", "HEAT_CAPACITY_FLUID", "PRANDTL_NUMBER_FLUID", ...
                "POROSITY", "TORTUOSITY", "FLOW_RESISTIVITY", "VISCOUS_LENGHT", "THERMAL_LENGTH", ...
                "FREQUENCY", "NUMBER_GAUSS_POINT"]);
            
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
            mixedporousElement2d4n@PorousElement2d4n(super_args{:});
            mixedporousElement2d4n.dofNames = cellstr(["DISPLACEMENT_SOLID_X"; "DISPLACEMENT_SOLID_Y"; ...
                "PORE_PRESSURE"]);
        end
        
        %Initialization
        function initialize(mixedporousElement2d4n)
            mixedporousElement2d4n.lengthX = computeLength(mixedporousElement2d4n.nodeArray(1).getCoords, ...
                mixedporousElement2d4n.nodeArray(2).getCoords);
            
            mixedporousElement2d4n.lengthY = computeLength(mixedporousElement2d4n.nodeArray(1).getCoords, ...
                mixedporousElement2d4n.nodeArray(4).getCoords);
            
            checkConvexity(mixedporousElement2d4n);
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(mixedElementporous2d4n)
            % Read Material Properties
            LAMBDA_SOLID = mixedElementporous2d4n.getPropertyValue('LAMBDA_SOLID');
            MUE_SOLID = mixedElementporous2d4n.getPropertyValue('MUE_SOLID');
            DAMPING_SOLID = mixedElementporous2d4n.getPropertyValue('DAMPING_SOLID');
            
            DENSITY_FLUID = mixedElementporous2d4n.getPropertyValue('DENSITY_FLUID');
            VISCOSITY_FLUID = mixedElementporous2d4n.getPropertyValue('VISCOSITY_FLUID');
            STANDARD_PRESSURE_FLUID = mixedElementporous2d4n.getPropertyValue('STANDARD_PRESSURE_FLUID');
            HEAT_CAPACITY_FLUID = mixedElementporous2d4n.getPropertyValue('HEAT_CAPACITY_FLUID');
            PRANDTL_NUMBER_FLUID = mixedElementporous2d4n.getPropertyValue('PRANDTL_NUMBER_FLUID');
            
            POROSITY = mixedElementporous2d4n.getPropertyValue('POROSITY');
            THERMAL_LENGTH = mixedElementporous2d4n.getPropertyValue('THERMAL_LENGTH');
            
            omega = mixedElementporous2d4n.getPropertyValue('FREQUENCY');
            p = mixedElementporous2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            
            DENSITY_SOLID = mixedElementporous2d4n.getPropertyValue('DENSITY_SOLID');

            TORTUOSITY = mixedElementporous2d4n.getPropertyValue('TORTUOSITY');
            FLOW_RESISTIVITY = mixedElementporous2d4n.getPropertyValue('FLOW_RESISTIVITY');
            VISCOUS_LENGHT = mixedElementporous2d4n.getPropertyValue('VISCOUS_LENGHT');
            
            
            % Determine relevant densities
            VISCOUS_DRAG = FLOW_RESISTIVITY * POROSITY^2 * (1 + ...
                4 * 1i * TORTUOSITY^2 * VISCOSITY_FLUID * DENSITY_FLUID / (FLOW_RESISTIVITY^2 * VISCOUS_LENGHT^2 * POROSITY^2));
            APPARENT_MASS_DENSITY = POROSITY * DENSITY_FLUID * (TORTUOSITY -1);
            EQUIVALENT_COUPLING_DENSITY = (-1) * APPARENT_MASS_DENSITY + 1i * VISCOUS_DRAG / omega;
            EQUIVALENT_FLUID_DENSITY = POROSITY * DENSITY_FLUID - EQUIVALENT_COUPLING_DENSITY;
            EQUIVALENT_SOLID_DENSITY = (1 - POROSITY) * DENSITY_SOLID - EQUIVALENT_COUPLING_DENSITY;
            
            % Determine Coefficients for EMat
            LAMBDA = (1 + 1i * DAMPING_SOLID) * LAMBDA_SOLID;
            MUE = (1 + 1i * DAMPING_SOLID) * MUE_SOLID;
            
            K = HEAT_CAPACITY_FLUID * STANDARD_PRESSURE_FLUID /(HEAT_CAPACITY_FLUID - ...
                (1 - HEAT_CAPACITY_FLUID)*(1 + 8 * VISCOSITY_FLUID/(1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID)* ...
                (1 + 1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID /(16 * VISCOSITY_FLUID))^0.5));
            R = POROSITY * K;
            Q = (1 - POROSITY) * K;
            GAMMA = POROSITY * ((EQUIVALENT_COUPLING_DENSITY /EQUIVALENT_FLUID_DENSITY) - (Q/R)) ;
            
            % Matrix valid for isotropic material only
            D = [LAMBDA + 2 * MUE,LAMBDA,0; ...
                LAMBDA,LAMBDA + 2 * MUE,0; ...
                0,0,MUE];
            
            K_Matrix = zeros(8,8);
            C_Matrix = zeros(8,4);
            H_Matrix = zeros(4,4);
            
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                zeta=g(i);
                for j=1:p
                    eta=g(j);
                        [N_mat, ~, Be, B, J] = computeShapeFunction(mixedElementporous2d4n, zeta, eta);
                        K_Matrix=K_Matrix+(w(i)*w(j)*det(J)*transpose(Be)*(D*Be));
                        C_Matrix=C_Matrix+(w(i)*w(j)*det(J)*GAMMA*transpose(N_mat) * B);
                        H_Matrix=H_Matrix+(w(i)*w(j)*det(J)*(POROSITY^2/EQUIVALENT_FLUID_DENSITY)*transpose(B) * B);
                end
            end
            
            stiffnessMatrix = [K_Matrix, - C_Matrix; zeros(4,8), H_Matrix];
        end
            
        function massMatrix = computeLocalMassMatrix(mixedporousElement2d4n)
            % Read Material Properties
            LAMBDA_SOLID = mixedporousElement2d4n.getPropertyValue('LAMBDA_SOLID');
            MUE_SOLID = mixedporousElement2d4n.getPropertyValue('MUE_SOLID');
            DAMPING_SOLID = mixedporousElement2d4n.getPropertyValue('DAMPING_SOLID');
            
            DENSITY_FLUID = mixedporousElement2d4n.getPropertyValue('DENSITY_FLUID');
            VISCOSITY_FLUID = mixedporousElement2d4n.getPropertyValue('VISCOSITY_FLUID');
            STANDARD_PRESSURE_FLUID = mixedporousElement2d4n.getPropertyValue('STANDARD_PRESSURE_FLUID');
            HEAT_CAPACITY_FLUID = mixedporousElement2d4n.getPropertyValue('HEAT_CAPACITY_FLUID');
            PRANDTL_NUMBER_FLUID = mixedporousElement2d4n.getPropertyValue('PRANDTL_NUMBER_FLUID');
            
            POROSITY = mixedporousElement2d4n.getPropertyValue('POROSITY');
            THERMAL_LENGTH = mixedporousElement2d4n.getPropertyValue('THERMAL_LENGTH');
            
            omega = mixedporousElement2d4n.getPropertyValue('FREQUENCY');
            p = mixedporousElement2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            
            DENSITY_SOLID = mixedporousElement2d4n.getPropertyValue('DENSITY_SOLID');

            TORTUOSITY = mixedporousElement2d4n.getPropertyValue('TORTUOSITY');
            FLOW_RESISTIVITY = mixedporousElement2d4n.getPropertyValue('FLOW_RESISTIVITY');
            VISCOUS_LENGHT = mixedporousElement2d4n.getPropertyValue('VISCOUS_LENGHT');
            
            
            % Determine relevant densities
            VISCOUS_DRAG = FLOW_RESISTIVITY * POROSITY^2 * (1 + ...
                4 * 1i * TORTUOSITY^2 * VISCOSITY_FLUID * DENSITY_FLUID / (FLOW_RESISTIVITY^2 * VISCOUS_LENGHT^2 * POROSITY^2));
            APPARENT_MASS_DENSITY = POROSITY * DENSITY_FLUID * (TORTUOSITY -1);
            EQUIVALENT_COUPLING_DENSITY = (-1) * APPARENT_MASS_DENSITY + 1i * VISCOUS_DRAG / omega;
            EQUIVALENT_FLUID_DENSITY = POROSITY * DENSITY_FLUID - EQUIVALENT_COUPLING_DENSITY;
            EQUIVALENT_SOLID_DENSITY = (1 - POROSITY) * DENSITY_SOLID - EQUIVALENT_COUPLING_DENSITY;
            EQUIVALENT_DENSITY = EQUIVALENT_SOLID_DENSITY - EQUIVALENT_COUPLING_DENSITY^2/EQUIVALENT_FLUID_DENSITY; 
            
            % Determine Coefficients for EMat
            LAMBDA = (1 + 1i * DAMPING_SOLID) * LAMBDA_SOLID;
            MUE = (1 + 1i * DAMPING_SOLID) * MUE_SOLID;
            
            K = HEAT_CAPACITY_FLUID * STANDARD_PRESSURE_FLUID /(HEAT_CAPACITY_FLUID - ...
                (1 - HEAT_CAPACITY_FLUID)*(1 + 8 * VISCOSITY_FLUID/(1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID)* ...
                (1 + 1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID /(16 * VISCOSITY_FLUID))^0.5));
            R = POROSITY * K;
            Q = (1 - POROSITY) * K;
            GAMMA = POROSITY * ((EQUIVALENT_COUPLING_DENSITY /EQUIVALENT_FLUID_DENSITY) - (Q/R)) ;
            
            
            M_Matrix = zeros(8,8);
            C_Matrix = zeros(8,4);
            Q_Matrix = zeros(4,4);
            
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                zeta=g(i);
                for j=1:p
                    eta=g(j);
                        [N_mat, N, ~, B, J] = computeShapeFunction(mixedporousElement2d4n, zeta, eta);
                        M_Matrix=M_Matrix+(w(i)*w(j)*det(J)*EQUIVALENT_DENSITY*transpose(N_mat) * N_mat);
                        C_Matrix=C_Matrix+(w(i)*w(j)*det(J)*GAMMA*transpose(N_mat) * B);
                        Q_Matrix=Q_Matrix+(w(i)*w(j)*det(J)*(POROSITY^2/R)*transpose(N)*(N));
                end
            end
            
            massMatrix = [M_Matrix, zeros(8,4); transpose(C_Matrix), Q_Matrix];
                
        end
        
        function dampingMatrix = computeLocalDampingMatrix(mixedporousElement2d4n)
            dampingMatrix = zeros(12,12);
        end
        
         function dofs = getDofList(element)
            dofs([1 3 5 7]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_X');
            dofs([2 4 6 8]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y');
            dofs([9 10 11 12]) = element.nodeArray.getDof('PORE_PRESSURE');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,16); 
            vals([1 3 5 7]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_X',step);
            vals([2 4 6 8]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_Y',step);
            vals([9 10 11 12]) = element.nodeArray.getDofValue('PORE_PRESSURE',step);
        end
        
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,16);
            [~, vals([1 3 5 7]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, vals([2 4 6 8]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, vals([9 10 11 12]), ~] = element.nodeArray.getDof('PORE_PRESSURE').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,16);
            [~, ~, vals([1 3 5 7])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, ~, vals([2 4 6 8])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, ~, vals([9 10 11 12])] = element.nodeArray.getDof('PORE_PRESSURE').getAllValues(step);
        end
        
        function update(mixedporousElement2d4n)
        end
        
        function F = computeLocalForceVector(mixedporousElement2d4n)
            F = zeros(1,16);
        end
        
    end
    
end
    

