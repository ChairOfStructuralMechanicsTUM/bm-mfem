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
 
        function c = barycenter(obj)
            %QUADRILATERALELEMENT.BARYCENTER returns the barycenter of the
            %element in the x-y plane
            diag1X = [obj.nodeArray(1).getX() obj.nodeArray(3).getX()];
            diag1Y = [obj.nodeArray(1).getY() obj.nodeArray(3).getY()];
            diag2X = [obj.nodeArray(2).getX() obj.nodeArray(4).getX()];
            diag2Y = [obj.nodeArray(2).getY() obj.nodeArray(4).getY()];
            
            try
                [c(1),c(2)] = polyxpoly(diag1X, diag1Y, diag2X, diag2Y);
            catch e
                if (strcmp(e.identifier,'MATLAB:subsassigndimmismatch'))
                    msg = ['Barycenter of element ', ...
                        num2str(obj.getId()), ...
                        ' could not be computed. ', ...
                        'Check the element for convexity.'];
                    causeException = MException('MATLAB:bm_mfem:barycenter',msg);
                    e = addCause(e,causeException);
                elseif (strcmp(e.identifier,'MATLAB:UndefinedFunction'))
                    msg = ['Function ''polyxpoly'' not found. ', ...
                        'Please install the Mapping Toolbox.'];
                    causeException = MException('MATLAB:bm_mfem:toolboxnotfound',msg);
                    e = addCause(e,causeException);
                end
                rethrow(e);
            end
        end

        function c = checkConvexity(obj)
            %QUADRILATERALELEMENT.CHECKCONVEXITY returns true, if the
            %element is convex
            diag1X = [obj.nodeArray(1).getX() obj.nodeArray(3).getX()];
            diag1Y = [obj.nodeArray(1).getY() obj.nodeArray(3).getY()];
            diag2X = [obj.nodeArray(2).getX() obj.nodeArray(4).getX()];
            diag2Y = [obj.nodeArray(2).getY() obj.nodeArray(4).getY()];
            
            c = 0;
            try
                res = polyxpoly(diag1X, diag1Y, diag2X, diag2Y);
            catch e
                if (strcmp(e.identifier,'MATLAB:subsassigndimmismatch'))
                    c = 0;
                elseif (strcmp(e.identifier,'MATLAB:UndefinedFunction'))
                    msg = ['Function ''polyxpoly'' not found. ', ...
                        'Please install the Mapping Toolbox.'];
                    causeException = MException('MATLAB:bm_mfem:toolboxnotfound',msg);
                    e = addCause(e,causeException);
                    rethrow(e);
                end

            end
            
            if ~isempty(res); c = 1; end
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
        
        function update(obj)
        end

    end

end
