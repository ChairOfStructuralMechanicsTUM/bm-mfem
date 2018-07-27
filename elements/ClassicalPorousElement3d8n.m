classdef ClassicalPorousElement3d8n < PorousElement3d8n
    %   POROUSELEMENT3D8N Summary of this class goes here
    %   Detailed explanation goes here
    
    %   DOF: Solid-Displacement & Fluid-Displacement (see: KANG 1995 / RUMPLER 2012)
    
    properties (Access = private)

    end
    
    methods
        % constructor
        function classicalporousElement3d8n = ClassicalPorousElement3d8n(id, nodeArray)
            
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
            classicalporousElement3d8n@PorousElement3d8n(super_args{:});
            classicalporousElement3d8n.dofNames = cellstr(['DISPLACEMENT_SOLID_X'; 'DISPLACEMENT_SOLID_Y'; 'DISPLACEMENT_SOLID_Z'; ...
                'DISPLACEMENT_FLUID_X'; 'DISPLACEMENT_FLUID_Y'; 'DISPLACEMENT_FLUID_Z']);
            
        end
        
        %Initialization
        function initialize(classicalporousElement3d8n)
            classicalporousElement3d8n.lengthX = computeLength(classicalporousElement3d8n.nodeArray(1).getCoords, ...
                classicalporousElement3d8n.nodeArray(2).getCoords);
            
            classicalporousElement3d8n.lengthY = computeLength(classicalporousElement3d8n.nodeArray(1).getCoords, ...
                classicalporousElement3d8n.nodeArray(4).getCoords);
            
            classicalporousElement3d8n.lengthY = computeLength(classicalporousElement3d8n.nodeArray(1).getCoords, ...
                classicalporousElement3d8n.nodeArray(5).getCoords);
            
            checkConvexity(classicalporousElement3d8n);
        end
        
        % member functions
        function stiffnessMatrix = computeLocalStiffnessMatrix(classicalporousElement3d8n)
            % Read Material Properties
            LAMBDA_SOLID = classicalporousElement3d8n.getPropertyValue('LAMBDA_SOLID');
            MUE_SOLID = classicalporousElement3d8n.getPropertyValue('MUE_SOLID');
            DAMPING_SOLID = classicalporousElement3d8n.getPropertyValue('DAMPING_SOLID');
            
            DENSITY_FLUID = classicalporousElement3d8n.getPropertyValue('DENSITY_FLUID');
            VISCOSITY_FLUID = classicalporousElement3d8n.getPropertyValue('VISCOSITY_FLUID');
            STANDARD_PRESSURE_FLUID = classicalporousElement3d8n.getPropertyValue('STANDARD_PRESSURE_FLUID');
            HEAT_CAPACITY_FLUID = classicalporousElement3d8n.getPropertyValue('HEAT_CAPACITY_FLUID');
            PRANDTL_NUMBER_FLUID = classicalporousElement3d8n.getPropertyValue('PRANDTL_NUMBER_FLUID');
            
            POROSITY = classicalporousElement3d8n.getPropertyValue('POROSITY');
            THERMAL_LENGTH = classicalporousElement3d8n.getPropertyValue('THERMAL_LENGTH');
            
            omega = classicalporousElement3d8n.getPropertyValue('FREQUENCY');
            p = classicalporousElement3d8n.getPropertyValue('NUMBER_GAUSS_POINT');
            
            % Determine Coefficients for EMat
            LAMBDA = (1 + 1i * DAMPING_SOLID) * LAMBDA_SOLID;
            MUE = (1 + 1i * DAMPING_SOLID) * MUE_SOLID;
            K = HEAT_CAPACITY_FLUID * STANDARD_PRESSURE_FLUID /(HEAT_CAPACITY_FLUID - ...
                (1 - HEAT_CAPACITY_FLUID)*(1 + 8 * VISCOSITY_FLUID/(1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID)* ...
                (1 + 1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID /(16 * VISCOSITY_FLUID))^0.5));
            R = POROSITY * K;
            Q = (1 - POROSITY) * K;
            A = LAMBDA + ((1 - POROSITY)^2/ POROSITY) * K;
            
            % Matrix valid for isotropic material only
            D_S = [A + 2 * MUE,A,A,0,0,0,; ...
                A,A + 2 * MUE,A,0,0,0,; ...
                A,A,A + 2 * MUE,0,0,0,; ...
                0,0,0,MUE,0,0; ...
                0,0,0,0,MUE,0,; ...
                0,0,0,0,0,MUE];
            D_SF = [Q,Q,Q,0,0,0; ...
                Q,Q,Q,0,0,0; ...
                Q,Q,Q,0,0,0; ...
                0,0,0,0,0,0; ...
                0,0,0,0,0,0; ...
                0,0,0,0,0,0];
            D_F = [R,R,R,0,0,0; ...
                R,R,R,0,0,0; ...
                R,R,R,0,0,0; ...
                0,0,0,0,0,0; ...
                0,0,0,0,0,0; ...
                0,0,0,0,0,0];
            
            stiffnessMatrix_S = zeros(24,24);
            stiffnessMatrix_SF = zeros(24,24);
            stiffnessMatrix_F = zeros(24,24);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                zeta=g(i);
                for j=1:p
                    eta=g(j);
                    for k=1:p
                        mue=g(k);
                        [~, ~, Be, ~, Jdet] = computeShapeFunction(classicalporousElement3d8n, zeta, eta, mue);
                        stiffnessMatrix_S=stiffnessMatrix_S+(w(i)*w(j)*w(k)*Jdet*transpose(Be)*(D_S*Be));
                        stiffnessMatrix_SF=stiffnessMatrix_SF+(w(i)*w(j)*w(k)*Jdet*transpose(Be)*(D_SF*Be));
                        stiffnessMatrix_F=stiffnessMatrix_F+(w(i)*w(j)*w(k)*Jdet*transpose(Be)*(D_F*Be));
                    end
                end
            end
            
            stiffnessMatrix = [stiffnessMatrix_S stiffnessMatrix_SF; stiffnessMatrix_SF stiffnessMatrix_F];
        end
        
        function massMatrix = computeLocalMassMatrix(classicalporousElement3d8n)
            % Read Material Properties
            DENSITY_SOLID = classicalporousElement3d8n.getPropertyValue('DENSITY_SOLID');
            
            DENSITY_FLUID = classicalporousElement3d8n.getPropertyValue('DENSITY_FLUID');
            VISCOSITY_FLUID = classicalporousElement3d8n.getPropertyValue('VISCOSITY_FLUID');
            
            POROSITY = classicalporousElement3d8n.getPropertyValue('POROSITY');
            TORTUOSITY = classicalporousElement3d8n.getPropertyValue('TORTUOSITY');
            FLOW_RESISTIVITY = classicalporousElement3d8n.getPropertyValue('FLOW_RESISTIVITY');
            VISCOUS_LENGHT = classicalporousElement3d8n.getPropertyValue('VISCOUS_LENGHT');
            
            omega = classicalporousElement3d8n.getPropertyValue('FREQUENCY');
            p = classicalporousElement3d8n.getPropertyValue('NUMBER_GAUSS_POINT');
            
            % Determine relevant densities
            VISCOUS_DRAG = FLOW_RESISTIVITY * POROSITY^2 * (1 + ...
                4 * 1i * TORTUOSITY^2 * VISCOSITY_FLUID * DENSITY_FLUID / (FLOW_RESISTIVITY^2 * VISCOUS_LENGHT^2 * POROSITY^2));
            APPARENT_MASS_DENSITY = POROSITY * DENSITY_FLUID * (TORTUOSITY -1);
            EQUIVALENT_COUPLING_DENSITY = (-1) * APPARENT_MASS_DENSITY + 1i * VISCOUS_DRAG / omega;
            EQUIVALENT_FLUID_DENSITY = POROSITY * DENSITY_FLUID - EQUIVALENT_COUPLING_DENSITY;
            EQUIVALENT_SOLID_DENSITY = (1 - POROSITY) * DENSITY_SOLID - EQUIVALENT_COUPLING_DENSITY;
            
            ElementmassMatrix=zeros(24,24);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                zeta=g(i);
                for j=1:p
                    eta=g(j);
                    for k=1:p
                        mue=g(k);
                        [N, ~, ~, ~, Jdet] = computeShapeFunction(classicalporousElement3d8n, zeta, eta, mue);
                        ElementmassMatrix = ElementmassMatrix + (w(i) * w(j) * w(k) * transpose(N) * N * Jdet);
                    end
                end
            end
            
            massMatrix = [EQUIVALENT_SOLID_DENSITY * ElementmassMatrix, EQUIVALENT_COUPLING_DENSITY * ElementmassMatrix; ...
                EQUIVALENT_COUPLING_DENSITY * ElementmassMatrix, EQUIVALENT_FLUID_DENSITY * ElementmassMatrix];
            
        end
        
        function dampingMatrix = computeLocalDampingMatrix(classicalporousElement3d8n)
            dampingMatrix = zeros(48,48);
        end
        
        function dofs = getDofList(element)
            dofs([1 4 7 10 13 16 19 22]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_X');
            dofs([2 5 8 11 14 17 20 23]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y');
            dofs([3 6 9 12 15 18 21 24]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_Z');
            dofs([25 28 31 34 37 40 43 46]) = element.nodeArray.getDof('DISPLACEMENT_FLUID_X');
            dofs([26 29 32 35 38 41 44 47]) = element.nodeArray.getDof('DISPLACEMENT_FLUID_Y');
            dofs([27 30 33 36 39 42 45 48]) = element.nodeArray.getDof('DISPLACEMENT_FLUID_Z');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,48);
            
            vals([1 4 7 10 13 16 19 22]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_X',step);
            vals([2 5 8 11 14 17 20 23]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_Y',step);
            vals([3 6 9 12 15 18 21 24]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_Z',step);
            vals([25 28 31 34 37 40 43 46]) = element.nodeArray.getDofValue('DISPLACEMENT_FLUID_X',step);
            vals([26 29 32 35 38 41 44 47]) = element.nodeArray.getDofValue('DISPLACEMENT_FLUID_Y',step);
            vals([27 30 33 36 39 42 45 48]) = element.nodeArray.getDofValue('DISPLACEMENT_FLUID_Z',step);
        end
        
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,48);
            
            [~, vals([1 4 7 10 13 16 19 22]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, vals([2 5 8 11 14 17 20 23]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, vals([3 6 9 12 15 18 21 24]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Z').getAllValues(step);
            [~, vals([25 28 31 34 37 40 43 46]), ~] = element.nodeArray.getDof('DISPLACEMENT_FLUID_X').getAllValues(step);
            [~, vals([26 29 32 35 38 41 44 47]), ~] = element.nodeArray.getDof('DISPLACEMENT_FLUID_Y').getAllValues(step);
            [~, vals([27 30 33 36 39 42 45 48]), ~] = element.nodeArray.getDof('DISPLACEMENT_FLUID_Z').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,48);
            
            [~, ~, vals([1 4 7 10 13 16 19 22])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, ~, vals([2 5 8 11 14 17 20 23])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, ~, vals([3 6 9 12 15 18 21 24])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Z').getAllValues(step);
            [~, ~, vals([25 28 31 34 37 40 43 46])] = element.nodeArray.getDof('DISPLACEMENT_FLUID_X').getAllValues(step);
            [~, ~, vals([26 29 32 35 38 41 44 47])] = element.nodeArray.getDof('DISPLACEMENT_FLUID_Y').getAllValues(step);
            [~, ~, vals([27 30 33 36 39 42 45 48])] = element.nodeArray.getDof('DISPLACEMENT_FLUID_Z').getAllValues(step);
        end
        
        function F = computeLocalForceVector(classicalporousElement3d8n)
            F = zeros(1,48);
        end
        
    end
    
end
    

