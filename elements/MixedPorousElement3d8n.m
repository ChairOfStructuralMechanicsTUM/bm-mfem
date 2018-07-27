classdef MixedPorousElement3d8n < PorousElement3d8n
    %   POROUSELEMENT3D8N Summary of this class goes here
    %   Detailed explanation goes here
    
    %   DOF: Solid-Displacement & Pore-Pressure (see: ALLARD 1998)
    
    properties (Access = private)

    end
    
    methods
        % constructor
        function mixedporousElement3d8n = MixedPorousElement3d8n(id, nodeArray)
            
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
            mixedporousElement3d8n@PorousElement3d8n(super_args{:});
            mixedporousElement3d8n.dofNames = cellstr(["DISPLACEMENT_SOLID_X"; "DISPLACEMENT_SOLID_Y"; "DISPLACEMENT_SOLID_Z"; "PORE_PRESSURE"]);
            
        end
        
        %Initialization
        function initialize(mixedporousElement3d8n)
            mixedporousElement3d8n.lengthX = computeLength(mixedporousElement3d8n.nodeArray(1).getCoords, ...
                mixedporousElement3d8n.nodeArray(2).getCoords);
            
            mixedporousElement3d8n.lengthY = computeLength(mixedporousElement3d8n.nodeArray(1).getCoords, ...
                mixedporousElement3d8n.nodeArray(4).getCoords);
            
            mixedporousElement3d8n.lengthY = computeLength(mixedporousElement3d8n.nodeArray(1).getCoords, ...
                mixedporousElement3d8n.nodeArray(5).getCoords);
            
            checkConvexity(mixedporousElement3d8n);
        end
        
        % member functions
        function stiffnessMatrix = computeLocalStiffnessMatrix(mixedporousElement3d8n)
            % Read Material Properties
            LAMBDA_SOLID = mixedporousElement3d8n.getPropertyValue('LAMBDA_SOLID');
            MUE_SOLID = mixedporousElement3d8n.getPropertyValue('MUE_SOLID');
            DAMPING_SOLID = mixedporousElement3d8n.getPropertyValue('DAMPING_SOLID');
            
            DENSITY_FLUID = mixedporousElement3d8n.getPropertyValue('DENSITY_FLUID');
            VISCOSITY_FLUID = mixedporousElement3d8n.getPropertyValue('VISCOSITY_FLUID');
            STANDARD_PRESSURE_FLUID = mixedporousElement3d8n.getPropertyValue('STANDARD_PRESSURE_FLUID');
            HEAT_CAPACITY_FLUID = mixedporousElement3d8n.getPropertyValue('HEAT_CAPACITY_FLUID');
            PRANDTL_NUMBER_FLUID = mixedporousElement3d8n.getPropertyValue('PRANDTL_NUMBER_FLUID');
            
            POROSITY = mixedporousElement3d8n.getPropertyValue('POROSITY');
            THERMAL_LENGTH = mixedporousElement3d8n.getPropertyValue('THERMAL_LENGTH');
            
            omega = mixedporousElement3d8n.getPropertyValue('FREQUENCY');
            p = mixedporousElement3d8n.getPropertyValue('NUMBER_GAUSS_POINT');

            TORTUOSITY = mixedporousElement3d8n.getPropertyValue('TORTUOSITY');
            FLOW_RESISTIVITY = mixedporousElement3d8n.getPropertyValue('FLOW_RESISTIVITY');
            VISCOUS_LENGHT = mixedporousElement3d8n.getPropertyValue('VISCOUS_LENGHT');
            
            % Determine relevant densities
            VISCOUS_DRAG = FLOW_RESISTIVITY * POROSITY^2 * (1 + ...
                4 * 1i * TORTUOSITY^2 * VISCOSITY_FLUID * DENSITY_FLUID / (FLOW_RESISTIVITY^2 * VISCOUS_LENGHT^2 * POROSITY^2));
            APPARENT_MASS_DENSITY = POROSITY * DENSITY_FLUID * (TORTUOSITY -1);
            EQUIVALENT_COUPLING_DENSITY = (-1) * APPARENT_MASS_DENSITY + 1i * VISCOUS_DRAG / omega;
            EQUIVALENT_FLUID_DENSITY = POROSITY * DENSITY_FLUID - EQUIVALENT_COUPLING_DENSITY;
            
            % Determine Coefficients further needed
            LAMBDA = (1 + 1i * DAMPING_SOLID) * LAMBDA_SOLID;
            MUE = (1 + 1i * DAMPING_SOLID) * MUE_SOLID;
            
            K = HEAT_CAPACITY_FLUID * STANDARD_PRESSURE_FLUID /(HEAT_CAPACITY_FLUID - ...
                (1 - HEAT_CAPACITY_FLUID)*(1 + 8 * VISCOSITY_FLUID/(1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID)* ...
                (1 + 1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID /(16 * VISCOSITY_FLUID))^0.5));
            R = POROSITY * K;
            Q = (1 - POROSITY) * K;
            GAMMA = POROSITY * ((EQUIVALENT_COUPLING_DENSITY /EQUIVALENT_FLUID_DENSITY) - (Q/R)) ;
            
            D = [LAMBDA + 2 * MUE,LAMBDA,LAMBDA,0,0,0,; ...
                LAMBDA,LAMBDA + 2 * MUE,LAMBDA,0,0,0,; ...
                LAMBDA,LAMBDA,LAMBDA + 2 * MUE,0,0,0,; ...
                0,0,0,MUE,0,0; ...
                0,0,0,0,MUE,0,; ...
                0,0,0,0,0,MUE];
            
            K_Matrix = zeros(24,24);
            C_Matrix = zeros(24,8);
            H_Matrix = zeros(8,8);
            
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                zeta=g(i);
                for j=1:p
                    eta=g(j);
                    for k=1:p
                        mue=g(k);
                        [N, ~, Be, B, Jdet] = computeShapeFunction(mixedporousElement3d8n, zeta, eta, mue);
                        K_Matrix=K_Matrix+(w(i)*w(j)*w(k)*Jdet*transpose(Be)*(D*Be));
                        C_Matrix=C_Matrix+(w(i)*w(j)*w(k)*Jdet*GAMMA*transpose(N) * B);
                        H_Matrix=H_Matrix+(w(i)*w(j)*w(k)*(POROSITY^2/EQUIVALENT_FLUID_DENSITY)*Jdet*transpose(B) * B);
                    end
                end
            end
            
            stiffnessMatrix = [K_Matrix, - C_Matrix; zeros(8,24), H_Matrix];
        end
        
        function massMatrix = computeLocalMassMatrix(mixedporousElement3d8n)
            % Read Material Properties
            DENSITY_FLUID = mixedporousElement3d8n.getPropertyValue('DENSITY_FLUID');
            VISCOSITY_FLUID = mixedporousElement3d8n.getPropertyValue('VISCOSITY_FLUID');
            STANDARD_PRESSURE_FLUID = mixedporousElement3d8n.getPropertyValue('STANDARD_PRESSURE_FLUID');
            HEAT_CAPACITY_FLUID = mixedporousElement3d8n.getPropertyValue('HEAT_CAPACITY_FLUID');
            PRANDTL_NUMBER_FLUID = mixedporousElement3d8n.getPropertyValue('PRANDTL_NUMBER_FLUID');
            
            POROSITY = mixedporousElement3d8n.getPropertyValue('POROSITY');
            THERMAL_LENGTH = mixedporousElement3d8n.getPropertyValue('THERMAL_LENGTH');
            
            omega = mixedporousElement3d8n.getPropertyValue('FREQUENCY');
            p = mixedporousElement3d8n.getPropertyValue('NUMBER_GAUSS_POINT');
            
            DENSITY_SOLID = mixedporousElement3d8n.getPropertyValue('DENSITY_SOLID');

            TORTUOSITY = mixedporousElement3d8n.getPropertyValue('TORTUOSITY');
            FLOW_RESISTIVITY = mixedporousElement3d8n.getPropertyValue('FLOW_RESISTIVITY');
            VISCOUS_LENGHT = mixedporousElement3d8n.getPropertyValue('VISCOUS_LENGHT');
            
            % Determine relevant densities
            VISCOUS_DRAG = FLOW_RESISTIVITY * POROSITY^2 * (1 + ...
                4 * 1i * TORTUOSITY^2 * VISCOSITY_FLUID * DENSITY_FLUID / (FLOW_RESISTIVITY^2 * VISCOUS_LENGHT^2 * POROSITY^2));
            APPARENT_MASS_DENSITY = POROSITY * DENSITY_FLUID * (TORTUOSITY -1);
            EQUIVALENT_COUPLING_DENSITY = (-1) * APPARENT_MASS_DENSITY + 1i * VISCOUS_DRAG / omega;
            EQUIVALENT_FLUID_DENSITY = POROSITY * DENSITY_FLUID - EQUIVALENT_COUPLING_DENSITY;
            EQUIVALENT_SOLID_DENSITY = (1 - POROSITY) * DENSITY_SOLID - EQUIVALENT_COUPLING_DENSITY;
            EQUIVALENT_DENSITY = EQUIVALENT_SOLID_DENSITY - EQUIVALENT_COUPLING_DENSITY^2/EQUIVALENT_FLUID_DENSITY; 
            
            % Determine Coefficients further needed
            K = HEAT_CAPACITY_FLUID * STANDARD_PRESSURE_FLUID /(HEAT_CAPACITY_FLUID - ...
                (1 - HEAT_CAPACITY_FLUID)*(1 + 8 * VISCOSITY_FLUID/(1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID)* ...
                (1 + 1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID /(16 * VISCOSITY_FLUID))^0.5));
            R = POROSITY * K;
            Q = (1 - POROSITY) * K;
            GAMMA = POROSITY * ((EQUIVALENT_COUPLING_DENSITY /EQUIVALENT_FLUID_DENSITY) - (Q/R)) ;
            
            M_Matrix = zeros(24,24);
            C_Matrix = zeros(24,8);
            Q_Matrix = zeros(8,8);
            
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                zeta=g(i);
                for j=1:p
                    eta=g(j);
                    for k=1:p
                        mue=g(k);
                        [N, Nf, ~, B, Jdet] = computeShapeFunction(mixedporousElement3d8n, zeta, eta, mue);
                        M_Matrix=M_Matrix+(w(i)*w(j)*w(k)*EQUIVALENT_DENSITY*Jdet*transpose(N) * N);
                        C_Matrix=C_Matrix+(w(i)*w(j)*w(k)*Jdet*GAMMA*transpose(N) * B);
                        Q_Matrix=Q_Matrix+(w(i)*w(j)*w(k)*(POROSITY^2/R)*Jdet*transpose(Nf)*(Nf));
                    end
                end
            end
            
            massMatrix = [M_Matrix, zeros(24,8); transpose(C_Matrix), Q_Matrix];
                
        end
        
        function dampingMatrix = computeLocalDampingMatrix(mixedporousElement3d8n)
            dampingMatrix = zeros(32,32);
        end
        
        function dofs = getDofList(element)
            dofs([1 4 7 10 13 16 19 22]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_X');
            dofs([2 5 8 11 14 17 20 23]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y');
            dofs([3 6 9 12 15 18 21 24]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_Z');
            dofs([25 26 27 28 29 30 31 32]) = element.nodeArray.getDof('PORE_PRESSURE');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,48);
            
            vals([1 4 7 10 13 16 19 22]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_X',step);
            vals([2 5 8 11 14 17 20 23]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_Y',step);
            vals([3 6 9 12 15 18 21 24]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_Z',step);
            vals([25 26 27 28 29 30 31 32]) = element.nodeArray.getDofValue('PORE_PRESSURE',step);
        end
        
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,48);
            
            [~, vals([1 4 7 10 13 16 19 22]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, vals([2 5 8 11 14 17 20 23]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, vals([3 6 9 12 15 18 21 24]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Z').getAllValues(step);
            [~, vals([25 26 27 28 29 30 31 32]), ~] = element.nodeArray.getDof('PORE_PRESSURE').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,48);
            
            [~, ~, vals([1 4 7 10 13 16 19 22])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, ~, vals([2 5 8 11 14 17 20 23])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, ~, vals([3 6 9 12 15 18 21 24])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Z').getAllValues(step);
            [~, ~, vals([25 26 27 28 29 30 31 32])] = element.nodeArray.getDof('PORE_PRESSURE').getAllValues(step);
        end
 
        function F = computeLocalForceVector(mixedporousElement3d8n)
            F = zeros(1,32);
        end
        
    end
    
end
    

