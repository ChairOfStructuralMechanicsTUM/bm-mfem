classdef PorousElement2d4n < Element
    %   PorousElement2d4n Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        lengthX
        lengthY
    end
    
    methods
        % constructor
        function porousElement2d4n = PorousElement2d4n(id, nodeArray, requiredProperties)
            if nargin == 0
                super_args = {};
            elseif nargin == 3
                super_args = {id, nodeArray, requiredProperties};
            end
            porousElement2d4n@Element(super_args{:});
        end
        
        function c = barycenter(porousElement2d4n)
            diag1X = [porousElement2d4n.nodeArray(1).getX() porousElement2d4n.nodeArray(3).getX()];
            diag1Y = [porousElement2d4n.nodeArray(1).getY() porousElement2d4n.nodeArray(3).getY()];
            diag2X = [porousElement2d4n.nodeArray(2).getX() porousElement2d4n.nodeArray(4).getX()];
            diag2Y = [porousElement2d4n.nodeArray(2).getY() porousElement2d4n.nodeArray(4).getY()];
            [c(1),c(2)] = polyxpoly(diag1X, diag1Y, diag2X, diag2Y);
        end
        % Check Convexity of quad
        function checkConvexity(porousElement2d4n)
            try 
                [~] = porousElement2d4n.barycenter();
            catch 
                error('Element %i is not convex', porousElement2d4n.getId());
            end
        end
        
        % member functions
        function [Nmat, Nf, Be, B, Jdet] = computeShapeFunction(porousElement2d4n,xi,eta)
            % Shape Function and Derivatives
            Nf = [(1-xi)*(1-eta)/4    (1+xi)*(1-eta)/4    (1+xi)*(1+eta)/4    (1-xi)*(1+eta)/4];
            N_Diff_Par = [-(1-eta)/4    (1-eta)/4   (1+eta)/4   -(1+eta)/4
                -(1-xi)/4     -(1+xi)/4   (1+xi)/4    (1-xi)/4];   
            
            Nmat = sparse(2,8);
            Nmat(1,1:2:end) = Nf(:);
            Nmat(2,2:2:end) = Nf(:);
            
            % Coordinates of the nodes forming one element
            ele_coords = zeros(4,2);
            for i=1:4
                ele_coords(i,1) = porousElement2d4n.nodeArray(i).getX;
                ele_coords(i,2) = porousElement2d4n.nodeArray(i).getY;
            end
            
            % Jacobian
            J = N_Diff_Par * ele_coords;
            Jdet = det(J);
            
            % Calculation of B-Matrix
            Be=J\N_Diff_Par;
            Be_x=Be(1,1:4);
            Be_y=Be(2,1:4);
            
            B=[Be_x(1),0,Be_x(2),0,Be_x(3),0,Be_x(4),0;
                0,Be_y(1),0,Be_y(2),0,Be_y(3),0,Be_y(4);
                Be_y(1),Be_x(1),Be_y(2),Be_x(2),Be_y(3),Be_x(3),Be_y(4),Be_x(4)];
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
            x = [obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_SOLID_X', step), ...
                obj.nodeArray(2).getX + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_SOLID_X', step), ...
                obj.nodeArray(3).getX + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_SOLID_X', step), ...
                obj.nodeArray(4).getX + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_SOLID_X', step),...
                obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_SOLID_X', step)];
            
            y = [obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_SOLID_Y', step), ...
                obj.nodeArray(2).getY + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_SOLID_Y', step), ...
                obj.nodeArray(3).getY + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_SOLID_Y', step), ...
                obj.nodeArray(4).getY + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_SOLID_Y', step),...
                obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_SOLID_Y', step)];
            
            pl = line(x,y);
        end
        
        function update(porousElement2d4n)
            
        end
        
    end
    
end


