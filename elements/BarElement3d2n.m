classdef BarElement3d2n < BarElement
    %BARELEMENT3D2N Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    properties (Access = private, Constant = true)
%         dofNames = cellstr(['DISPLACEMENT_X'; 'DISPLACEMENT_Y'; 'DISPLACEMENT_Z']);
    end
    
    methods
        % constructor
        function barElement3d2n = BarElement3d2n(id, nodeArray, material, crossSectionArea)
            
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 4
                super_args = {id; material; crossSectionArea};
            end
            
            % call the super class constructor
            barElement3d2n@BarElement(super_args{:});
            
            barElement3d2n.dofNames = cellstr(['DISPLACEMENT_X'; 'DISPLACEMENT_Y'; 'DISPLACEMENT_Z']);
            barElement3d2n.requiredProperties = [ ];
            barElement3d2n.required3dProperties = [ ];
            
            % the constructor
            if nargin > 0
                if (length(nodeArray) == 2 && isa(nodeArray,'Node'))
                    barElement3d2n.nodeArray = nodeArray;
                else
                    error('problem with the nodes in element %d', id);
                end
                
%                 barElement3d2n.addDofs(barElement3d2n.dofNames);
                
                barElement3d2n.length = computeLength(barElement3d2n.nodeArray(1).getCoords, ...
                    barElement3d2n.nodeArray(2).getCoords);
            end
            
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
        
        function crossSectionArea = getCrossSectionArea(barElement)
            crossSectionArea = barElement.crossSectionArea;
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
            factor = (barElement.getMaterial().getValue('YOUNGS_MODULUS') ...
                * barElement.crossSectionArea) ...
                / (barElement.length^3);
            stiffnessMatrix = factor * stiffnessMatrix;
        end
        
        function forceVector = computeLocalForceVector(barElement)
            forceVector = zeros(1,6);
        end
        
        % Computation of the Internal Element Stresses
        function stressValue = computeElementStress(barElement, step)
            dist = barElement.nodeArray(2).getCoords - barElement.nodeArray(1).getCoords;
            CX = dist(1)/barElement.length;
            CY = dist(2)/barElement.length;
            CZ = dist(3)/barElement.length;
            nodalDisplacement = getResponseDofArray(barElement, step);
            stressValue = barElement.getMaterial().getValue('YOUNGS_MODULUS')...
                /barElement.length * [-CX  -CY  -CZ  CX  CY  CZ]*nodalDisplacement;  %Winkel überprüfen stets positiv
        end
        
     
    end
    
end

