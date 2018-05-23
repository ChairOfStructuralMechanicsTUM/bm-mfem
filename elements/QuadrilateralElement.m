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
        
        function checkConvexity(quadrilateralElement)
            % Check Convexity of quad
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
             
            z = [quadrilateralElement.nodeArray(1).getZ, quadrilateralElement.nodeArray(2).getZ, ... 
                 quadrilateralElement.nodeArray(3).getZ, quadrilateralElement.nodeArray(4).getZ, ...
                 quadrilateralElement.nodeArray(1).getZ];
             
            pl = line(x,y,z);
        end
        
        function pl = drawDeformed(obj, step, scaling)
            x = [obj.nodeArray(1).getX, ...
                    obj.nodeArray(2).getX, ...
                    obj.nodeArray(3).getX, ...
                    obj.nodeArray(4).getX,...
                    obj.nodeArray(1).getX];
                
            try obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step);
                x = x + scaling .* [obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step), ...
                        obj.nodeArray(2).getDofValue('DISPLACEMENT_X', step), ...
                        obj.nodeArray(3).getDofValue('DISPLACEMENT_X', step), ...
                        obj.nodeArray(4).getDofValue('DISPLACEMENT_X', step),...
                        obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step)];
            catch
            end
            
            y = [obj.nodeArray(1).getY, ...
                    obj.nodeArray(2).getY, ...
                    obj.nodeArray(3).getY, ...
                    obj.nodeArray(4).getY,...
                    obj.nodeArray(1).getY];
            
            try obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step);
                y = y + scaling .* [obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step), ...
                        obj.nodeArray(2).getDofValue('DISPLACEMENT_Y', step), ...
                        obj.nodeArray(3).getDofValue('DISPLACEMENT_Y', step), ...
                        obj.nodeArray(4).getDofValue('DISPLACEMENT_Y', step),...
                        obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step)];
            catch
            end
            
            z = [obj.nodeArray(1).getZ, ...
                    obj.nodeArray(2).getZ, ...
                    obj.nodeArray(3).getZ, ...
                    obj.nodeArray(4).getZ,...
                    obj.nodeArray(1).getZ];
            
            try obj.nodeArray(1).getDofValue('DISPLACEMENT_Z', step);
                z = z + scaling .* [obj.nodeArray(1).getDofValue('DISPLACEMENT_Z', step), ...
                    obj.nodeArray(2).getDofValue('DISPLACEMENT_Z', step), ...
                    obj.nodeArray(3).getDofValue('DISPLACEMENT_Z', step), ...
                    obj.nodeArray(4).getDofValue('DISPLACEMENT_Z', step),...
                    obj.nodeArray(1).getDofValue('DISPLACEMENT_Z', step)];
            catch
            end
            
            pl = line(x,y,z);
        end

        function update(quadrilateralElement)
        end
        
        function F = computeLocalForceVector(quadrilateralElement)
            F = zeros(1,12);
        end

    end

end
