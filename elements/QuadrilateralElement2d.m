classdef (Abstract) QuadrilateralElement2d < Element
    %QUADRILATERALELEMENT Base class for quadrilateral elements
    %   Detailed explanation goes here
    
    properties (Access = protected)
        lengthX
        lengthY
    end
    
    methods
        % Constructor
        function quadrilateralElement = QuadrilateralElement2d(id, nodeArray, requiredProperties)
            if nargin == 0
                super_args = {};
            elseif nargin == 3
                super_args = {id, nodeArray, requiredProperties};
            end
            
            quadrilateralElement@Element(super_args{:});
        end
 
        function c = barycenter(quadrilateralElement)
            diag1X = [quadrilateralElement.nodeArray(1).getX() quadrilateralElement.nodeArray(3).getX()];
            diag1Y = [quadrilateralElement.nodeArray(1).getY() quadrilateralElement.nodeArray(3).getY()];
            diag2X = [quadrilateralElement.nodeArray(2).getX() quadrilateralElement.nodeArray(4).getX()];
            diag2Y = [quadrilateralElement.nodeArray(2).getY() quadrilateralElement.nodeArray(4).getY()];
            [c(1),c(2)] = polyxpoly(diag1X, diag1Y, diag2X, diag2Y);
        end
        % Check Convexity of quad
        function checkConvexity(quadrilateralElement)
            try 
                [~] = quadrilateralElement.barycenter();
            catch 
                error('Element %i is not convex', quadrilateralElement.getId());
            end
        end
        
        function pl = draw(quadrilateralElement)   
            x = [quadrilateralElement.nodeArray(1).getX, quadrilateralElement.nodeArray(2).getX, ... 
                 quadrilateralElement.nodeArray(3).getX, quadrilateralElement.nodeArray(4).getX,...
                 quadrilateralElement.nodeArray(1).getX];
             
            y = [quadrilateralElement.nodeArray(1).getY, quadrilateralElement.nodeArray(2).getY, ... 
                 quadrilateralElement.nodeArray(3).getY, quadrilateralElement.nodeArray(4).getY, ...
                 quadrilateralElement.nodeArray(1).getY];
             
            pl = line(x,y);
        end

        function update(quadrilateralElement)
        end
        
        function F = computeLocalForceVector(quadrilateralElement)
            F = zeros(1,12);
        end

    end

end

