classdef SpringDamperElement3d2n < Element
    %SPRINGDAMPERELEMENT3D2N Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        length0 %initial length
    end
    
    methods
        % constructor
        function springDamperElement = SpringDamperElement3d2n(id, nodeArray, properties)
            if nargin == 0
                super_args = {};
            elseif nargin == 3
                super_args = {id; properties};
            end
            
            % call the super class constructor
            springDamperElement@Element(super_args{:});
            springDamperElement.dofNames = cellstr(['DISPLACEMENT_X'; 'DISPLACEMENT_Y'; 'DISPLACEMENT_Z']);
            springDamperElement.requiredProperties = [ ];
            springDamperElement.required3dProperties = cellstr(["ELEMENTAL_STIFFNESS", "ELEMENTAL_DAMPING"]);
            
            % the constructor
            if nargin > 0
                if (length(nodeArray) == 2 && isa(nodeArray,'Node'))
                    springDamperElement.nodeArray = nodeArray;
                else
                    error('problem with the nodes in element %d', id);
                end
                
                springDamperElement.length0 = computeLength(springDamperElement.nodeArray(1).getCoords, ...
                    springDamperElement.nodeArray(2).getCoords);
            end
            
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(springDamperElement)
           stiffness = springDamperElement.getPropertyValue('ELEMENTAL_STIFFNESS');
           stiffnessMatrix = zeros(6);
           for ii = 1:size(stiffnessMatrix,1)/2
               stiffnessMatrix(ii, ii) = stiffness(ii);
               stiffnessMatrix(ii+3, ii+3) = stiffness(ii);
               stiffnessMatrix(ii, ii+3) = -stiffness(ii);
               stiffnessMatrix(ii+3, ii) = -stiffness(ii);
           end
        end
        
        function dampingMatrix = computeLocalDampingMatrix(springDamperElement)
           stiffness = springDamperElement.getPropertyValue('ELEMENTAL_DAMPING');
           dampingMatrix = zeros(6);
           for ii = 1:size(dampingMatrix,1)/2
               dampingMatrix(ii, ii) = stiffness(ii);
               dampingMatrix(ii+3, ii+3) = stiffness(ii);
               dampingMatrix(ii, ii+3) = -stiffness(ii);
               dampingMatrix(ii+3, ii) = -stiffness(ii);
           end
        end
        
        function forceVector = computeLocalForceVector(ele)
           forceVector = zeros(1,6);
           stiffness = ele.getPropertyValue('ELEMENTAL_STIFFNESS');
           nodes = ele.getNodes;
           du = nodes(2).getDofValue('DISPLACEMENT_X') - nodes(1).getDofValue('DISPLACEMENT_X');
           dv = nodes(2).getDofValue('DISPLACEMENT_Y') - nodes(1).getDofValue('DISPLACEMENT_Y');
           dw = nodes(2).getDofValue('DISPLACEMENT_Z') - nodes(1).getDofValue('DISPLACEMENT_Z');
           
           forceVector(1) = + stiffness(1) * du;
           forceVector(2) = + stiffness(2) * dv;
           forceVector(3) = + stiffness(3) * dw;
           forceVector(4) = - stiffness(1) * du;
           forceVector(5) = - stiffness(2) * dv;
           forceVector(6) = - stiffness(3) * dw;
        end
        
        function update(springDamperElement)
            springDamperElement.length = computeLength(springDamperElement.nodeArray(1).getCoords, ...
                    springDamperElement.nodeArray(2).getCoords);
        end
        
        function c = barycenter(springDamperElement)
            nodes = springDamperElement.getNodes;
            c = (nodes(1).getCoords + nodes(2).getCoords) ./ 2;
        end
        
        function plot = draw(springDamperElement)
            plot = line(springDamperElement.nodeArray.getX, springDamperElement.nodeArray.getY);
        end
        
        function plot = drawDeformed(springDamperElement)
            plot = line(springDamperElement.nodeArray.getX + ...
                springDamperElement.nodeArray.getDofValue('DISPLACEMENT_X'), ...
                springDamperElement.nodeArray.getY + ...
                springDamperElement.nodeArray.getDofValue('DISPLACEMENT_Y'));
        end
    end
    
end

