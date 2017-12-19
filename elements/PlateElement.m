classdef (Abstract) PlateElement < Element
    
    %PLATEELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        lengthX
        lengthY
    end
    
    methods
        % Constructor
        function plateElement = PlateElement(id, nodeArray, requiredProperties)
            if nargin == 0
                super_args = {};
            elseif nargin == 3
                super_args = {id, nodeArray, requiredProperties};
            end
            
            plateElement@Element(super_args{:});
        end
 
        function c = barycenter(plateElement)
            diag1X = [plateElement.nodeArray(1).getX() plateElement.nodeArray(3).getX()];
            diag1Y = [plateElement.nodeArray(1).getY() plateElement.nodeArray(3).getY()];
            diag2X = [plateElement.nodeArray(2).getX() plateElement.nodeArray(4).getX()];
            diag2Y = [plateElement.nodeArray(2).getY() plateElement.nodeArray(4).getY()];
            [c(1),c(2)] = polyxpoly(diag1X, diag1Y, diag2X, diag2Y);
        end
        % Check Convexity of quad
        function checkConvexity(plateElement)
            try 
                [~] = plateElement.barycenter();
            catch 
                error('Element %i is not convex', plateElement.getId());
            end
        end
        
        function pl = draw(plateElement)   
            x = [plateElement.nodeArray(1).getX, plateElement.nodeArray(2).getX, ... 
                 plateElement.nodeArray(3).getX, plateElement.nodeArray(4).getX,...
                 plateElement.nodeArray(1).getX];
             
            y = [plateElement.nodeArray(1).getY, plateElement.nodeArray(2).getY, ... 
                 plateElement.nodeArray(3).getY, plateElement.nodeArray(4).getY, ...
                 plateElement.nodeArray(1).getY];
             
            z = [plateElement.nodeArray(1).getZ, plateElement.nodeArray(2).getZ, ... 
                 plateElement.nodeArray(3).getZ, plateElement.nodeArray(4).getZ, ...
                 plateElement.nodeArray(1).getZ];
             
            pl = line(x,y,z);
        end
        
        function pl = drawDeformed(plateElement, step, scaling)
            x = [plateElement.nodeArray(1).getX, plateElement.nodeArray(2).getX, ...
                 plateElement.nodeArray(3).getX, plateElement.nodeArray(4).getX, ...
                 plateElement.nodeArray(1).getX];
             
            y = [plateElement.nodeArray(1).getY, plateElement.nodeArray(2).getY, ... 
                 plateElement.nodeArray(3).getY, plateElement.nodeArray(4).getY, ...
                 plateElement.nodeArray(1).getY];
             
            z = [plateElement.nodeArray(1).getZ + scaling * plateElement.nodeArray(1).getDofValue('DISPLACEMENT_Z', step), ... 
                 plateElement.nodeArray(2).getZ + scaling * plateElement.nodeArray(2).getDofValue('DISPLACEMENT_Z', step), ... 
                 plateElement.nodeArray(3).getZ + scaling * plateElement.nodeArray(3).getDofValue('DISPLACEMENT_Z', step), ...
                 plateElement.nodeArray(4).getZ + scaling * plateElement.nodeArray(4).getDofValue('DISPLACEMENT_Z', step), ...
                 plateElement.nodeArray(1).getZ + scaling * plateElement.nodeArray(1).getDofValue('DISPLACEMENT_Z', step)];
            
            pl = line(x,y,z);
        end

        function update(plateElement)
        end
        
        function F = computeLocalForceVector(plateElement)
            F = zeros(1,12);
        end

    end

end
