classdef BarElement2d2n < BarElement
    %BARELEMENT2D2N: 2D Truss Element with only Normal Forces
    
    
    %  ==================================  %
    %% PROPERTIES OF THE 2D TRUSS ELEMENT %%
    %  ==================================  %
    
    properties (Access = private)
    end
    
    properties (Access = private, Constant = true)
%         dofNames = cellstr(['DISPLACEMENT_X'; 'DISPLACEMENT_Y']);
    end
    
    
    %  ==========================================  %
    %% METHODS OF APPLIED TO THE 2D TRUSS ELEMENT %%
    %  ==========================================  %
    
    
    methods
        
        % constructor
        
        function barElement2d2n = BarElement2d2n(id, nodeArray, material, crossSectionArea)
            
            if nargin == 0
                super_args = {};
            elseif nargin == 4
                super_args = {id; material; crossSectionArea};
            end
                
%             barElement2d2n@BarElement(id, material, crossSectionArea);
            barElement2d2n@BarElement(super_args{:});
            
            barElement2d2n.dofNames = cellstr(['DISPLACEMENT_X'; 'DISPLACEMENT_Y']);
            barElement2d2n.requiredProperties = [ ];
            barElement2d2n.required3dProperties = [ ];
            
            if nargin > 0
                
                if (length(nodeArray) == 2 && isa(nodeArray,'Node'))
                    barElement2d2n.nodeArray = nodeArray;
                else
                    error('problem with the nodes in element %d', id);
                end
                
%                 barElement2d2n.addDofs(barElement2d2n.dofNames);
                
                barElement2d2n.length = computeLength(barElement2d2n.nodeArray(1).getCoords, ...
                    barElement2d2n.nodeArray(2).getCoords);
            end
        end
        
        
        % getter functions
        
        function responseDoF = getResponseDofArray(barElement)
            
            responseDoF = zeros(4,1);
            
            for itNodes = 1:1:2
                nodalDof = barElement.nodeArray(itNodes).getDofArray;
                nodalDof = nodalDof.';
                
                for itDof = 2:(-1):1
                    responseDoF(2*itNodes-(itDof-1),1) = nodalDof(3-itDof).getValue;
                end
            end
            
        end
        
        
        function crossSectionArea = getCrossSectionArea(barElement)
            crossSectionArea = barElement.crossSectionArea;
        end
        
        
        % member functions
        
        % Computation of the Local Stiffness Matrix
        function stiffnessMatrix = computeLocalStiffnessMatrix(barElement)
            dist = barElement.nodeArray(2).getCoords - barElement.nodeArray(1).getCoords;
            x21 = dist(1);
            y21 = dist(2);
            stiffnessMatrix = [x21*x21 x21*y21 -x21*x21 -x21*y21; ...
                x21*y21 y21*y21 -x21*y21 -y21*y21; ...
                -x21*x21 -x21*y21 x21*x21 x21*y21; ...
                -x21*y21 -y21*y21 x21*y21 y21*y21];
            factor = (barElement.getMaterial().getValue('YOUNGS_MODULUS') ...
                * barElement.crossSectionArea) ...
                / (barElement.length^3);
            stiffnessMatrix = factor * stiffnessMatrix;
        end
        
        function forceVector = computeLocalForceVector(barElement)
            forceVector = zeros(1,6);
        end
        
        % Computation of the Element Stress
        function stressValue = computeElementStress(barElement)
            dist = barElement.nodeArray(2).getCoords - barElement.nodeArray(1).getCoords;
            cos = dist(1)/barElement.length;
            sin = dist(2)/barElement.length;
            nodalDisplacement = getResponseDofArray(barElement);
            stressValue = barElement.getMaterial().getValue('YOUNGS_MODULUS') ...
                /barElement.length* [-cos  -sin  cos  sin]*nodalDisplacement; %Winkel überprüfen stets positiv
            
        end
        
    end
    
end

