classdef PlaneStressElement3d4n < QuadrilateralElement
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    methods
        %Constructor
        function planeStressElement3d4n = PlaneStressElement3d4n(id,nodeArray)
            
            requiredPropertyNames = cellstr(["YOUNGS_MODULUS", "POISSON_RATIO", ...
                                             "THICKNESS", "NUMBER_GAUSS_POINT", ...
                                             "DENSITY"]);
                                         
            %define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 4 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end

            %call the super class constructor
            planeStressElement3d4n@QuadrilateralElement(super_args{:});
            planeStressElement3d4n.dofNames = cellstr(["DISPLACEMENT_X", "DISPLACEMENT_Y"]);
        end
    
        %Initialization
        function initialize(planeStressElement3d4n)
            planeStressElement3d4n.lengthY = computeLength(planeStressElement3d4n.nodeArray(1).getCoords, ...
                planeStressElement3d4n.nodeArray(2).getCoords);

            planeStressElement3d4n.lengthY = computeLength(planeStressElement3d4n.nodeArray(1).getCoords, ...
                planeStressElement3d4n.nodeArray(4).getCoords);

            checkConvexity(planeStressElement3d4n);
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

        function [N_mat, N, B, J] = computeShapeFunction(planeStressElemet3d4n,xi,eta)
            % Shape Function and Derivatives                    
            N = [(1-xi)*(1-eta)/4    (1+xi)*(1-eta)/4    (1+xi)*(1+eta)/4    (1-xi)*(1+eta)/4];  

            N_Diff_Par = [-(1-eta)/4    (1-eta)/4   (1+eta)/4   -(1+eta)/4
                          -(1-xi)/4     -(1+xi)/4   (1+xi)/4    (1-xi)/4];

            N_mat = sparse(2,8);
            N_mat(1,1:2:end) = N(:);
            N_mat(2,2:2:end) = N(:);

            % Coordinates of the nodes forming one element 
            ele_coords = zeros(4,2); 
            for i=1:4
                ele_coords(i,1) = planeStressElemet3d4n.nodeArray(i).getX;
                ele_coords(i,2) = planeStressElemet3d4n.nodeArray(i).getY;
            end

            % Jacobian 
            J = N_Diff_Par * ele_coords;

            N_Diff = J \ N_Diff_Par;

            % Assembling the B Matrix
            B = sparse(3,8);
            B(1,1:2:end) = N_Diff(1,:);     
            B(3,2:2:end) = N_Diff(1,:);
            B(2,2:2:end) = N_Diff(2,:);
            B(3,1:2:end) = N_Diff(2,:);

        end

        function stiffnessMatrix = computeLocalStiffnessMatrix(planeStressElement3d4n)
            EModul = planeStressElement3d4n.getPropertyValue('YOUNGS_MODULUS');
            prxy = planeStressElement3d4n.getPropertyValue('POISSON_RATIO');
            nr_gauss_points = planeStressElement3d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            thickness = planeStressElement3d4n.getPropertyValue('THICKNESS');

            % Moment-Curvature Equations
            D = [1    prxy    0; prxy     1   0; 0    0   (1-prxy)/2];
            % Material Bending Matrix D
            D = D * EModul * thickness / (1-prxy^2);

            [w,g] = returnGaussPoint(nr_gauss_points);
            stiffnessMatrix = sparse(8,8);
            for xi = 1 : nr_gauss_points
                for eta = 1 : nr_gauss_points
                    [~, ~, B, J] = computeShapeFunction(planeStressElement3d4n,g(xi),g(eta));

                    stiffnessMatrix = stiffnessMatrix + ...
                        B' * D *B * det(J) * w(xi) * w(eta);
                end
            end
        end
      
        function pl = drawDeformed(planeStressElement3d4n, step, scaling)
    
            x = [planeStressElement3d4n.nodeArray(1).getX + scaling * planeStressElement3d4n.nodeArray(1).getDofValue('DISPLACEMENT_X', step), ... 
                 planeStressElement3d4n.nodeArray(2).getX + scaling * planeStressElement3d4n.nodeArray(2).getDofValue('DISPLACEMENT_X', step), ... 
                 planeStressElement3d4n.nodeArray(3).getX + scaling * planeStressElement3d4n.nodeArray(3).getDofValue('DISPLACEMENT_X', step), ...
                 planeStressElement3d4n.nodeArray(4).getX + scaling * planeStressElement3d4n.nodeArray(4).getDofValue('DISPLACEMENT_X', step), ...
                 planeStressElement3d4n.nodeArray(1).getX + scaling * planeStressElement3d4n.nodeArray(1).getDofValue('DISPLACEMENT_X', step)];
             
            y = [planeStressElement3d4n.nodeArray(1).getY + scaling * planeStressElement3d4n.nodeArray(1).getDofValue('DISPLACEMENT_Y', step), ... 
                 planeStressElement3d4n.nodeArray(2).getY + scaling * planeStressElement3d4n.nodeArray(2).getDofValue('DISPLACEMENT_Y', step), ... 
                 planeStressElement3d4n.nodeArray(3).getY + scaling * planeStressElement3d4n.nodeArray(3).getDofValue('DISPLACEMENT_Y', step), ...
                 planeStressElement3d4n.nodeArray(4).getY + scaling * planeStressElement3d4n.nodeArray(4).getDofValue('DISPLACEMENT_Y', step), ...
                 planeStressElement3d4n.nodeArray(1).getY + scaling * planeStressElement3d4n.nodeArray(1).getDofValue('DISPLACEMENT_Y', step)];
            
            z = [planeStressElement3d4n.nodeArray(1).getZ, planeStressElement3d4n.nodeArray(2).getZ, ...
                 planeStressElement3d4n.nodeArray(3).getZ, planeStressElement3d4n.nodeArray(4).getZ, ...
                 planeStressElement3d4n.nodeArray(1).getZ];
             
            pl = line(x,y,z);
        end
        
        function dofs = getDofList(element)
            dofs([1 3 5 7]) = element.nodeArray.getDof('DISPLACEMENT_X'); 
           
            dofs([2 4 6 8]) = element.nodeArray.getDof('DISPLACEMENT_Y');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,8);
            
            vals([1 3 5 7]) = element.nodeArray.getDofValue('DISPLACEMENT_X',step);

            vals([2 4 6 8]) = element.nodeArray.getDofValue('DISPLACEMENT_Y',step);
        end
        
    end
end



