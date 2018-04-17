classdef PlaneStressElement3d3n < TriangularElement
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    methods
        %Constructor
        function planeStressElement3d3n = PlaneStressElement3d3n(id,nodeArray)
            
            requiredPropertyNames = cellstr(["YOUNGS_MODULUS", "POISSON_RATIO", ...
                                             "THICKNESS", "NUMBER_GAUSS_POINT", ...
                                             "DENSITY"]);
                                         
            %define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 3 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end

            %call the super class constructor
            planeStressElement3d3n@TriangularElement(super_args{:});
            planeStressElement3d3n.dofNames = cellstr(["DISPLACEMENT_X", "DISPLACEMENT_Y"]);
        end
    
        %Initialization
        function initialize(planeStressElement3d3n)
            checkConvexity(planeStressElement3d3n);
        end

        function responseDoF = getResponseDofArray(planeStressElement, step)

            responseDoF = zeros(8,1);
            for itNodes = 1:1:4
                nodalDof = planeStressElement.nodeArray(itNodes).getDofArray;
                nodalDof = nodalDof.';

                for itDof = 2:(-1):1
                    responseDoF(3*itNodes-(itDof-1),1) = nodalDof(4-itDof).getValue(step);
                end
            end
        end

        function [N_mat, N, B, J] = computeShapeFunction(planeStressElement3d3n,g)
            % Shape Function and Derivatives                    
            N = g;  
            
            N_mat = sparse(2,6);
            N_mat(1,1:2:end) = N(:);
            N_mat(2,2:2:end) = N(:);

            % Coordinates of the nodes forming one element 
            ele_coords = zeros(3,2); 
            for i=1:3
                ele_coords(1:3,1) = 1;
                ele_coords(i,2) = planeStressElement3d3n.nodeArray(i).getX;
                ele_coords(i,3) = planeStressElement3d3n.nodeArray(i).getY;
            end

            % Jacobian 
            J = 1/2 * det(ele_coords);
            
            N_Diff = [ele_coords(2,3)-ele_coords(3,3),ele_coords(3,3)-ele_coords(1,3),ele_coords(1,3)-ele_coords(2,3); ...
                        -ele_coords(2,2)+ele_coords(3,2),-ele_coords(3,2)+ele_coords(1,2),-ele_coords(1,2)+ele_coords(2,2)];
                    
            N_Diff = N_Diff / (2*J);
            
            % Assembling the B Matrix
            B = sparse(3,6);
            B(1,1:2:end) = N_Diff(1,:);
            B(3,2:2:end) = N_Diff(1,:);
            B(2,2:2:end) = N_Diff(2,:);
            B(3,1:2:end) = N_Diff(2,:);

        end

        function stiffnessMatrix = computeLocalStiffnessMatrix(planeStressElement3d3n)
            EModul = planeStressElement3d3n.getPropertyValue('YOUNGS_MODULUS');
            prxy = planeStressElement3d3n.getPropertyValue('POISSON_RATIO');
            nr_gauss_points = planeStressElement3d3n.getPropertyValue('NUMBER_GAUSS_POINT');
            thickness = planeStressElement3d3n.getPropertyValue('THICKNESS');

            % Moment-Curvature Equations
            D = [1    prxy    0; prxy     1   0; 0    0   (1-prxy)/2];
            % Material Matrix D
            D = D * EModul * thickness / (1-prxy^2);
            
            stiffnessMatrix = sparse(6,6);            
            for i = 1 : nr_gauss_points

                [w,g] = returnGaussPointTrig(nr_gauss_points,i);
                [~, ~, B, J] = computeShapeFunction(planeStressElement3d3n,g);
                stiffnessMatrix = stiffnessMatrix + ...
                w * B' * D * B * J;
            end
        end
      
        function pl = drawDeformed(planeStressElement3d3n, step, scaling)
    
            x = [planeStressElement3d3n.nodeArray(1).getX + scaling * planeStressElement3d3n.nodeArray(1).getDofValue('DISPLACEMENT_X', step), ... 
                 planeStressElement3d3n.nodeArray(2).getX + scaling * planeStressElement3d3n.nodeArray(2).getDofValue('DISPLACEMENT_X', step), ... 
                 planeStressElement3d3n.nodeArray(3).getX + scaling * planeStressElement3d3n.nodeArray(3).getDofValue('DISPLACEMENT_X', step), ...
                 planeStressElement3d3n.nodeArray(1).getX + scaling * planeStressElement3d3n.nodeArray(1).getDofValue('DISPLACEMENT_X', step)];
             
            y = [planeStressElement3d3n.nodeArray(1).getY + scaling * planeStressElement3d3n.nodeArray(1).getDofValue('DISPLACEMENT_Y', step), ... 
                 planeStressElement3d3n.nodeArray(2).getY + scaling * planeStressElement3d3n.nodeArray(2).getDofValue('DISPLACEMENT_Y', step), ... 
                 planeStressElement3d3n.nodeArray(3).getY + scaling * planeStressElement3d3n.nodeArray(3).getDofValue('DISPLACEMENT_Y', step), ...
                 planeStressElement3d3n.nodeArray(1).getY + scaling * planeStressElement3d3n.nodeArray(1).getDofValue('DISPLACEMENT_Y', step)];
            
            z = [planeStressElement3d3n.nodeArray(1).getZ, planeStressElement3d3n.nodeArray(2).getZ, ...
                 planeStressElement3d3n.nodeArray(3).getZ, planeStressElement3d3n.nodeArray(1).getZ];
             
            pl = line(x,y,z);
        end
        
        function dofs = getDofList(element)
            dofs([1 3 5]) = element.nodeArray.getDof('DISPLACEMENT_X'); 
           
            dofs([2 4 6]) = element.nodeArray.getDof('DISPLACEMENT_Y');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,6);
            
            vals([1 3 5]) = element.nodeArray.getDofValue('DISPLACEMENT_X',step);

            vals([2 4 6]) = element.nodeArray.getDofValue('DISPLACEMENT_Y',step);
        end
        
    end
end



