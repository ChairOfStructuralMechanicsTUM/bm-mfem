classdef BarElement2d2n < LinearElement
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
        
        function barElement2d2n = BarElement2d2n(id, nodeArray)
            
            requiredPropertyNames = cellstr(["YOUNGS_MODULUS", "CROSS_SECTION", "DENSITY"]);
            
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 2 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
                
%             barElement2d2n@BarElement(id, material, crossSectionArea);
            barElement2d2n@LinearElement(super_args{:});
            
            barElement2d2n.dofNames = cellstr(['DISPLACEMENT_X'; 'DISPLACEMENT_Y']);
%             barElement2d2n.requiredProperties = cellstr(["YOUNGS_MODULUS", "CROSS_SECTION", "DENSITY"]);
%             barElement2d2n.required3dProperties = [ ];
            
            if nargin > 0
                
%                 if (length(nodeArray) == 2 && isa(nodeArray,'Node'))
%                     barElement2d2n.nodeArray = nodeArray;
%                 else
%                     error('problem with the nodes in element %d', id);
%                 end
                
%                 barElement2d2n.addDofs(barElement2d2n.dofNames);
                
                barElement2d2n.length = computeLength(barElement2d2n.nodeArray(1).getCoords, ...
                    barElement2d2n.nodeArray(2).getCoords);
            end
        end
        
        function initialize(element) 
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
            factor = (barElement.getProperties().getValue('YOUNGS_MODULUS') ...
                * barElement.getProperties().getValue('CROSS_SECTION')) ...
                / (barElement.length^3);
            stiffnessMatrix = factor * stiffnessMatrix;
        end
        
        function forceVector = computeLocalForceVector(barElement)
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
            massMatrix = mass * diag([.5 .5 .5 .5]);
        end
        
        % Computation of the Element Stress
        function stressValue = computeElementStress(barElement)
            dist = barElement.nodeArray(2).getCoords - barElement.nodeArray(1).getCoords;
            cos = dist(1)/barElement.length;
            sin = dist(2)/barElement.length;
            nodalDisplacement = getResponseDofArray(barElement);
            stressValue = barElement.getProperties().getValue('YOUNGS_MODULUS') ...
                /barElement.length* [-cos  -sin  cos  sin]*nodalDisplacement; %Winkel überprüfen stets positiv
            
        end
        
    end
    
end

