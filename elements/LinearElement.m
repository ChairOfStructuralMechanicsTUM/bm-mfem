classdef (Abstract) LinearElement < Element
    %BARELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        crossSectionArea
        length
        localSystem
    end
    
    methods
        % constructor
        function e = LinearElement(id, nodeArray, requiredProperties)
            if nargin == 0
                super_args = {};
            elseif nargin == 3
                super_args = {id, nodeArray, requiredProperties};
            end
            
            e@Element(super_args{:});
        end
        
        % getter functions
        function len = getLength(element)
            if element.length == 0
                element.length = computeLength(element.nodeArray(1).getCoords, ...
                    element.nodeArray(2).getCoords);
            end
            len = element.length;
        end
        
        function lsystem = getLocalSystem(element)
            lsystem = element.localSystem;
        end
        
        % member functions
        function update(barElement)
            barElement.length = computeLength(barElement.nodeArray(1).getCoords, ...
                barElement.nodeArray(2).getCoords);
        end
        
        function c = barycenter(barElement)
            nodes = barElement.getNodes;
            c = (nodes(1).getCoords + nodes(2).getCoords) ./ 2;
        end
        
        function pl = draw(barElement)
            pl = line(barElement.nodeArray.getX, barElement.nodeArray.getY);
        end
        
        function pl = drawDeformed(barElement, step, scaling)
            pl = line(barElement.nodeArray.getX' + scaling * barElement.nodeArray.getDofValue('DISPLACEMENT_X', step), ...
                barElement.nodeArray.getY' + scaling * barElement.nodeArray.getDofValue('DISPLACEMENT_Y', step));
        end
        
    end
    
    methods (Access = protected)
        function lsystem = computeLocalSystem(ele, dirZ)
            %GETLOCALSYSTEM computes the local system of the 3d linear
            %element. The optional parameter DIRZ can be used to specify
            %the local z axis. By default, dirZ=[0 0 -1] and if the local x
            %axis and the global z axis coincide dirZ=[0 -1 0]
            
            dirX = ele.nodeArray(2).getCoords - ele.nodeArray(1).getCoords;
            dirX = dirX ./ norm(dirX);

            %change the local z, if the element is oriented in along the
            %global z axis
            if nargin == 1
                if abs(cos(10)) > dirX(3)
                    dirZ = [0 0 -1];
                else
                    dirZ = [0 -1 0];
                end
            end
            %project the local z axis onto the plane perpendicular to the
            %local x axis
            dirZ_proj = dirZ - (dot(dirZ, dirX) / norm(dirX)) .* dirX;            
            dirZ = dirZ_proj ./ norm(dirZ);
            
            dirY = cross(dirZ,dirX);
            
            lsystem=zeros(3);
            lsystem(1,:) = dirX;
            lsystem(2,:) = dirY;
            lsystem(3,:) = dirZ;
        end
        
        function tMat = getTransformationMatrix(ele)
            nDofs = length(ele.getDofs);
            lsystem = ele.getLocalSystem();
            tMat = zeros(nDofs);
            for ii = 1:3:nDofs
               tMat(ii:ii+2,ii:ii+2) = lsystem;
            end 
        end
        
        function cp = copyElement(obj)
            cp = copyElement@Element(obj);
            cp.crossSectionArea = obj.crossSectionArea;
            cp.length = obj.length;
        end
    end
    
end

