classdef BarElement3d2n < LinearElement
    %BARELEMENT3D2N Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    methods
        % constructor
        function barElement3d2n = BarElement3d2n(id, nodeArray)
            
            requiredPropertyNames = cellstr(["YOUNGS_MODULUS", "CROSS_SECTION", "DENSITY"]);
            
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 2 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            % call the super class constructor
            barElement3d2n@LinearElement(super_args{:});
            barElement3d2n.dofNames = cellstr(['DISPLACEMENT_X'; 'DISPLACEMENT_Y'; 'DISPLACEMENT_Z']);
                        
        end
        
        function initialize(element)
            element.localSystem = element.computeLocalSystem();
            element.length = computeLength(element.nodeArray(1).getCoords, ...
                    element.nodeArray(2).getCoords);
        end
        
        % getter functions
        
        function responseDoF = getResponseDofArray(barElement, step)
            
            responseDoF = zeros(6,1);
            
            for itNodes = 1:1:2
                nodalDof = barElement.nodeArray(itNodes).getDofArray;
                nodalDof = nodalDof.';
                
                for itDof = 3:(-1):1
                    responseDoF(3*itNodes-(itDof-1),1) = nodalDof(4-itDof).getValue(step);
                end
            end
            
        end
        
        % member functions
        function stiffnessMatrix = computeLocalStiffnessMatrix(barElement)
            dist = barElement.nodeArray(2).getCoords - barElement.nodeArray(1).getCoords;
            x21 = dist(1);
            y21 = dist(2);
            z21 = dist(3);
            stiffnessMatrix = [x21*x21 x21*y21 x21*z21 -x21*x21 -x21*y21 -x21*z21; ...
                x21*y21 y21*y21 y21*z21 -x21*y21 -y21*y21 -y21*z21; ...
                x21*z21 y21*z21 z21*z21 -x21*z21 -y21*z21 -z21*z21; ...
                -x21*x21 -x21*y21 -x21*z21 x21*x21 x21*y21 x21*z21; ...
                -x21*y21 -y21*y21 -y21*z21 x21*y21 y21*y21 y21*z21; ...
                -x21*z21 -y21*z21 -z21*z21 x21*z21 y21*z21 z21*z21];
            factor = (barElement.getProperties().getValue('YOUNGS_MODULUS') ...
                * barElement.getProperties().getValue('CROSS_SECTION')) ...
                / (barElement.length^3);
            stiffnessMatrix = factor * stiffnessMatrix;
        end
        
        function forceVector = computeLocalForceVector(element)
            forceVector = zeros(1,6);
        end
        
        function massMatrix = computeLocalMassMatrix(element)
            %COMPUTELOCALMASSMATRIX return the mass matrix
            %the mass is computed from density, length, and cross section
            %area and distributed equally among both nodes
            
            density = element.getProperties().getValue('DENSITY');
            length = element.length;
            area = element.getProperties().getValue('CROSS_SECTION');
            mass = density * length * area;
            massMatrix = mass * diag([.5 .5 .5 .5 .5 .5]);
        end
        
        function dampingMatrix = computeLocalDampingMatrix(element)
            dampingMatrix = zeros(6);
        end
        
        % Computation of the Internal Element Stresses
        function stressValue = computeElementStress(barElements, step)
            stressValue = zeros(length(barElements), 1);
            for ii = 1:length(barElements)
                barElement = barElements(ii);
                dist = barElement.nodeArray(2).getCoords - barElement.nodeArray(1).getCoords;
                CX = dist(1)/barElement.length;
                CY = dist(2)/barElement.length;
                CZ = dist(3)/barElement.length;
                nodalDisplacement = getResponseDofArray(barElement, step);
                stressValue(ii) = barElement.getProperties().getValue('YOUNGS_MODULUS')...
                    /barElement.length * [-CX  -CY  -CZ  CX  CY  CZ] * nodalDisplacement;  %Winkel überprüfen stets positiv
            end
        end
     
    end
    
end

