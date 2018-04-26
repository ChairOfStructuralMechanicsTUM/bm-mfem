classdef (Abstract) QuadrilateralElement < Element
    %QUADRILATERALELEMENT Base class for quadrilateral elements
    %   Detailed explanation goes here
    
    properties (Access = protected)
        lengthX
        lengthY
    end
    
    methods
        % Constructor
        function quadrilateralElement = QuadrilateralElement(id, nodeArray, requiredProperties)
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
        
        function pl = draw(obj)
            x = [obj.nodeArray(1).getX, obj.nodeArray(2).getX, ...
                obj.nodeArray(3).getX, obj.nodeArray(4).getX,...
                obj.nodeArray(1).getX];
            
            y = [obj.nodeArray(1).getY, obj.nodeArray(2).getY, ...
                obj.nodeArray(3).getY, obj.nodeArray(4).getY, ...
                obj.nodeArray(1).getY];
            
            if(all(obj.getNodes().getDimension == 3))
                z = [obj.nodeArray(1).getZ, obj.nodeArray(2).getZ, ...
                    obj.nodeArray(3).getZ, obj.nodeArray(4).getZ, ...
                    obj.nodeArray(1).getZ];
                
                pl = line(x,y,z);
            else
                pl = line(x,y);
            end
            
        end
        
        function pl = drawDeformed(obj, step, scaling)
            x = [obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(2).getX + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(3).getX + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(4).getX + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_X', step),...
                obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step)];
            
            y = [obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(2).getY + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(3).getY + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(4).getY + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Y', step),...
                obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step)];
            
            if(all(obj.getNodes().getDimension == 3))
                z = [obj.nodeArray(1).getZ + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Z', step), ...
                    obj.nodeArray(2).getZ + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Z', step), ...
                    obj.nodeArray(3).getZ + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Z', step), ...
                    obj.nodeArray(4).getZ + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Z', step),...
                    obj.nodeArray(1).getZ + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Z', step)];
                
                pl = line(x,y,z);
            else
                pl = line(x,y);
            end
        end
        
        function update(quadrilateralElement)
        end
        
        function F = computeLocalForceVector(quadrilateralElement)
            F = zeros(1,12);
        end

    end

end
