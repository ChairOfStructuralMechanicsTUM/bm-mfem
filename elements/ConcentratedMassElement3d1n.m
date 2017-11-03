classdef ConcentratedMassElement3d1n < Element
    %CONCENTRATEDMASSELEMENT3D1N Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        % constructor
        function element = ConcentratedMassElement3d1n(id, nodeArray)
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                super_args = {id};
            end
            
            % call the super class constructor
            element@Element(super_args{:});
            element.dofNames = cellstr(['DISPLACEMENT_X'; 'DISPLACEMENT_Y'; 'DISPLACEMENT_Z']);
            element.requiredProperties = cellstr("ELEMENTAL_MASS");
            element.required3dProperties = cellstr("VOLUME_ACCELERATION");
            
            % the constructor
            if nargin > 0
                if (length(nodeArray) == 1 && isa(nodeArray,'Node'))
                    element.nodeArray = nodeArray;
                else
                    error('problem with the nodes in element %d', id);
                end
            end
            
        end
        
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(element)
            stiffnessMatrix = zeros(3);
        end
        
        function forceVector = computeLocalForceVector(ele)
           mass = ele.getPropertyValue('ELEMENTAL_MASS');
           volume_acceleration = ele.getPropertyValue('VOLUME_ACCELERATION');
           forceVector = volume_acceleration .* mass;
        end
        
        function massMatrix = computeLocalMassMatrix(element)
           mass = element.getPropertyValue('ELEMENTAL_MASS');
           massMatrix = zeros(3);
           massMatrix(1,1) = mass;
           massMatrix(2,2) = mass;
           massMatrix(3,3) = mass;
        end
        
        function dampingMatrix = computeLocalDampingMatrix(element)
            dampingMatrix = zeros(3);
        end
        
        
        
        function update(element)
            
        end
        
        function c = barycenter(element)
            nodes = element.getNodes;
            c = nodes(1).getCoords;
        end
        
        function plot = draw(element)
            plot = line(element.nodeArray.getX, element.nodeArray.getY, ...
                'color','b','Marker','d','MarkerFaceColor','b');
        end
        
        function plot = drawDeformed(element)
            plot = line(element.nodeArray.getX + ...
                element.nodeArray.getDofValue('DISPLACEMENT_X'), ...
                element.nodeArray.getY + ...
                element.nodeArray.getDofValue('DISPLACEMENT_Y'), ...
                'color','r','Marker','d','MarkerFaceColor','r');
        end
        
    end
    
end

