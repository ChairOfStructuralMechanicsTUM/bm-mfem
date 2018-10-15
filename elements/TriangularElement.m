classdef TriangularElement < Element
    %TRIANGULARELEMENT Base class for quadraliteral elements
    %   Detailed explanation goes here
    
    properties (Access = protected)
    end
    
    methods
        % Constructor
        function triangularElement = TriangularElement(id, nodeArray, requiredProperties)
            if nargin == 0
                super_args = {};
            elseif nargin == 3
                super_args = {id, nodeArray, requiredProperties};
            end
            
            triangularElement@Element(super_args{:});
        end
 
        function c = barycenter(triangularElement)
            for i = 1 : 3
                X(i) = triangularElement.nodeArray(i).getX();
                Y(i) = triangularElement.nodeArray(i).getY();
                Z(i) = triangularElement.nodeArray(i).getZ();
            end
            c(1) = mean(X);
            c(2) = mean(Y);
            c(3) = mean(Z);
        end
        % Check Convexity of triangle
        function checkConvexity(triangularElement)
            try 
                [~] = triangularElement.barycenter();
            catch 
                error('Element %i is not convex', triangularElement.getId());
            end
        end
        
        function pl = draw(triangularElement)   
            x = [triangularElement.nodeArray(1).getX, triangularElement.nodeArray(2).getX, ... 
                 triangularElement.nodeArray(3).getX, triangularElement.nodeArray(1).getX];
             
            y = [triangularElement.nodeArray(1).getY, triangularElement.nodeArray(2).getY, ... 
                 triangularElement.nodeArray(3).getY, triangularElement.nodeArray(1).getY];
             
            z = [triangularElement.nodeArray(1).getZ, triangularElement.nodeArray(2).getZ, ... 
                 triangularElement.nodeArray(3).getZ, triangularElement.nodeArray(1).getZ];
             
            pl = line(x,y,z);
        end

        function update(triangularElement)
        end
        
        function F = computeLocalForceVector(triangularElement)
            F = zeros(1,9);
        end
    end
end

