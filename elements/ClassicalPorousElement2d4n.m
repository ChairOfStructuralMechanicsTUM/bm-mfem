classdef ClassicalPorousElement2d4n < PorousElement2d4n
    %   ClassicalPorousElement2d4n Summary of this class goes here
    %   Detailed explanation goes here
    
    %   DOF: Solid-Displacement & Fluid-Displacement (see: KANG 1995 / RUMPLER 2012)
    
    properties (Access = private)
    end
    
    methods
        
        % constructor
        function classicalporousElement2d4n = ClassicalPorousElement2d4n(id,nodeArray)
            
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
            classicalporousElement2d4n@PorousElement2d4n(super_args{:});
            classicalporousElement2d4n.dofNames = cellstr(['DISPLACEMENT_SOLID_X'; 'DISPLACEMENT_SOLID_Y'; ...
                'DISPLACEMENT_FLUID_X'; 'DISPLACEMENT_FLUID_Y']);
        end
        
        %Initialization
        function initialize(classicalporousElement2d4n)
            classicalporousElement2d4n.lengthX = computeLength(classicalporousElement2d4n.nodeArray(1).getCoords, ...
                classicalporousElement2d4n.nodeArray(2).getCoords);
            
            classicalporousElement2d4n.lengthY = computeLength(classicalporousElement2d4n.nodeArray(1).getCoords, ...
                classicalporousElement2d4n.nodeArray(4).getCoords);
            
            checkConvexity(classicalporousElement2d4n);
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(classicalporousElement2d4n)
            % Read Material Properties
            LAMBDA_SOLID = classicalporousElement2d4n.getPropertyValue('LAMBDA_SOLID');
            MUE_SOLID = classicalporousElement2d4n.getPropertyValue('MUE_SOLID');
            DAMPING_SOLID = classicalporousElement2d4n.getPropertyValue('DAMPING_SOLID');
            
            DENSITY_FLUID = classicalporousElement2d4n.getPropertyValue('DENSITY_FLUID');
            VISCOSITY_FLUID = classicalporousElement2d4n.getPropertyValue('VISCOSITY_FLUID');
            STANDARD_PRESSURE_FLUID = classicalporousElement2d4n.getPropertyValue('STANDARD_PRESSURE_FLUID');
            HEAT_CAPACITY_FLUID = classicalporousElement2d4n.getPropertyValue('HEAT_CAPACITY_FLUID');
            PRANDTL_NUMBER_FLUID = classicalporousElement2d4n.getPropertyValue('PRANDTL_NUMBER_FLUID');
            
            POROSITY = classicalporousElement2d4n.getPropertyValue('POROSITY');
            THERMAL_LENGTH = classicalporousElement2d4n.getPropertyValue('THERMAL_LENGTH');
            
            omega = classicalporousElement2d4n.getPropertyValue('FREQUENCY');
            p = classicalporousElement2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            
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
            D_S = [A + 2 * MUE,A,0; ...
                A,A + 2 * MUE,0; ...
                0,0,MUE];
            D_SF = [Q,Q,0; ...
                Q,Q,0; ...
                0,0,0];
            D_F = [R,R,0; ...
                R,R,0; ...
                0,0,0];
            
            stiffnessMatrix_S = zeros(8,8);
            stiffnessMatrix_F = zeros(8,8);
            stiffnessMatrix_SF = zeros(8,8);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [~, ~, Be, ~, J] = computeShapeFunction(classicalporousElement2d4n, xi, eta);
                    stiffnessMatrix_S=stiffnessMatrix_S+(w(i)*w(j)*det(J)*transpose(Be)*(D_S*Be));
                    stiffnessMatrix_F=stiffnessMatrix_F+(w(i)*w(j)*det(J)*transpose(Be)*(D_F*Be));
                    stiffnessMatrix_SF=stiffnessMatrix_SF+(w(i)*w(j)*det(J)*transpose(Be)*(D_SF*Be));
                end
            end
            
            stiffnessMatrix = [stiffnessMatrix_S stiffnessMatrix_SF; stiffnessMatrix_SF stiffnessMatrix_F];
        end
        
        function massMatrix = computeLocalMassMatrix(classicalporousElement2d4n)
            % Read Material Properties
            DENSITY_SOLID = classicalporousElement2d4n.getPropertyValue('DENSITY_SOLID');
            
            DENSITY_FLUID = classicalporousElement2d4n.getPropertyValue('DENSITY_FLUID');
            VISCOSITY_FLUID = classicalporousElement2d4n.getPropertyValue('VISCOSITY_FLUID');
            
            POROSITY = classicalporousElement2d4n.getPropertyValue('POROSITY');
            TORTUOSITY = classicalporousElement2d4n.getPropertyValue('TORTUOSITY');
            FLOW_RESISTIVITY = classicalporousElement2d4n.getPropertyValue('FLOW_RESISTIVITY');
            VISCOUS_LENGHT = classicalporousElement2d4n.getPropertyValue('VISCOUS_LENGHT');
            
            omega = classicalporousElement2d4n.getPropertyValue('FREQUENCY');
            p = classicalporousElement2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            
            % Determine relevant densities
            VISCOUS_DRAG = FLOW_RESISTIVITY * POROSITY^2 * (1 + ...
                4 * 1i * TORTUOSITY^2 * VISCOSITY_FLUID * DENSITY_FLUID / (FLOW_RESISTIVITY^2 * VISCOUS_LENGHT^2 * POROSITY^2));
            APPARENT_MASS_DENSITY = POROSITY * DENSITY_FLUID * (TORTUOSITY -1);
            EQUIVALENT_COUPLING_DENSITY = (-1) * APPARENT_MASS_DENSITY + 1i * VISCOUS_DRAG / omega;
            EQUIVALENT_FLUID_DENSITY = POROSITY * DENSITY_FLUID - EQUIVALENT_COUPLING_DENSITY;
            EQUIVALENT_SOLID_DENSITY = (1 - POROSITY) * DENSITY_SOLID - EQUIVALENT_COUPLING_DENSITY;
            
            ElementmassMatrix=zeros(8,8);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [N_mat, ~, ~, ~, J] = computeShapeFunction(classicalporousElement2d4n, xi, eta);
                    ElementmassMatrix = ElementmassMatrix + (w(i) * w(j) * transpose(N_mat) * N_mat * det(J));
                end
            end
            
            massMatrix = [EQUIVALENT_SOLID_DENSITY * ElementmassMatrix, EQUIVALENT_COUPLING_DENSITY * ElementmassMatrix; ...
                EQUIVALENT_COUPLING_DENSITY * ElementmassMatrix, EQUIVALENT_FLUID_DENSITY * ElementmassMatrix];
            
        end
        
        function dampingMatrix = computeLocalDampingMatrix(classicalporousElement2d4n)
            dampingMatrix = zeros(16,16);
        end
        
         function dofs = getDofList(element)
            dofs([1 3 5 7]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_X');
            dofs([2 4 6 8]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y');
            dofs([9 11 13 15]) = element.nodeArray.getDof('DISPLACEMENT_FLUID_X');
            dofs([10 12 14 16]) = element.nodeArray.getDof('DISPLACEMENT_FLUID_Y');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,16); 
            vals([1 3 5 7]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_X',step);
            vals([2 4 6 8]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_Y',step);
            vals([9 11 13 15]) = element.nodeArray.getDofValue('DISPLACEMENT_FLUID_X',step);
            vals([10 12 14 16]) = element.nodeArray.getDofValue('DISPLACEMENT_FLUID_Y',step);
        end
        
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,16);
            [~, vals([1 3 5 7]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, vals([2 4 6 8]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, vals([9 11 13 15]), ~] = element.nodeArray.getDof('DISPLACEMENT_FLUID_X').getAllValues(step);
            [~, vals([10 12 14 16]), ~] = element.nodeArray.getDof('DISPLACEMENT_FLUID_Y').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,16);
            [~, ~, vals([1 3 5 7])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, ~, vals([2 4 6 8])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, ~, vals([9 11 13 15])] = element.nodeArray.getDof('DISPLACEMENT_FLUID_X').getAllValues(step);
            [~, ~, vals([10 12 14 16])] = element.nodeArray.getDof('DISPLACEMENT_FLUID_Y').getAllValues(step);
        end
        
        function update(classicalporousElement2d4n)
        end
        
        function F = computeLocalForceVector(classicalporousElement2d4n)
            F = zeros(1,16);
        end
        
    end
    
end
    

