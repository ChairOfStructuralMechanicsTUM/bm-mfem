classdef PlaneStressElement3d3n < TriangularElement
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    methods
        %Constructor
        function obj = PlaneStressElement3d3n(id,nodeArray)
            
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
            obj@TriangularElement(super_args{:});
            obj.dofNames = cellstr(["DISPLACEMENT_X", "DISPLACEMENT_Y"]);
        end
    
        %Initialization
        function initialize(obj)
        end

        function [N_mat, N, B, J] = computeShapeFunction(obj,tCoord)
            
            % Shape Function and Derivatives   
             N = tCoord;           

             N_Diff_Par = [1 0 -1
                            0 1 -1];

            
            N_mat = sparse(2,6);
            N_mat(1,1:2:end) = N(:);
            N_mat(2,2:2:end) = N(:);
            
            % Coordinates of the nodes forming one element 
            ele_coords = zeros(3,2); 
            for i=1:3
                ele_coords(i,1) = obj.nodeArray(i).getX;
                ele_coords(i,2) = obj.nodeArray(i).getY;
            end

            
            % Jacobian 
            J = N_Diff_Par * ele_coords;
            
            Ntest_Diff = J \ N_Diff_Par;
                    
 
            % Assembling the B Matrix
            B = sparse(3,6);
            B(1,1:2:end) = Ntest_Diff(1,:);
            B(3,2:2:end) = Ntest_Diff(1,:);
            B(2,2:2:end) = Ntest_Diff(2,:);
            B(3,1:2:end) = Ntest_Diff(2,:);

        end

        function stiffnessMatrix = computeLocalStiffnessMatrix(obj)
            EModul = obj.getPropertyValue('YOUNGS_MODULUS');
            prxy = obj.getPropertyValue('POISSON_RATIO');
            nr_gauss_points = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            thickness = obj.getPropertyValue('THICKNESS');

            % Moment-Curvature Equations
            D = [1    prxy    0; prxy     1   0; 0    0   (1-prxy)/2];
            % Material Matrix D
            D = D * EModul * thickness / (1-prxy^2);
            
            stiffnessMatrix = sparse(6,6);            
            for i = 1 : nr_gauss_points

                [w,g] = returnGaussPointTrig(nr_gauss_points,i);
                [~, ~, B, J] = computeShapeFunction(obj,g);
                stiffnessMatrix = stiffnessMatrix + ...
                w * B' * D * B * 0.5 * det(J);
            end
        end
        
        function massMatrix = computeLocalMassMatrix(obj)
            %Formulation of the Massmatrix based on the Shape Functions
            
            density = obj.getPropertyValue('DENSITY');
            thickness = obj.getPropertyValue('THICKNESS');
            nr_gauss_points = obj.getPropertyValue('NUMBER_GAUSS_POINT');

            dens_mat = sparse(2,2);
            dens_mat(1,1) = density*thickness; 
            dens_mat(2,2) = dens_mat(1,1); 

            massMatrix = sparse(6,6);
            for i = 1 : nr_gauss_points
                
                [w,g] = returnGaussPointTrig(nr_gauss_points,i);
                [N_mat, ~, ~, J] = computeShapeFunction(obj,g);
  
                massMatrix = massMatrix + N_mat' * dens_mat * N_mat * 0.5 * det(J) * w;
            end
        end
        
        function dampingMatrix = computeLocalDampingMatrix(obj)
            eProperties = obj.getProperties;
            dampingMatrix = sparse(6,6);

            if (eProperties.hasValue('RAYLEIGH_ALPHA'))
                alpha = eProperties.getValue('RAYLEIGH_ALPHA');
                dampingMatrix = dampingMatrix + alpha * element.computeLocalMassMatrix;
            end

            if (eProperties.hasValue('RAYLEIGH_BETA'))
                beta = eProperties.getValue('RAYLEIGH_BETA');
                dampingMatrix = dampingMatrix + beta * element.computeLocalStiffnessMatrix;
            end
        end
            
        function pl = drawDeformed(obj, step, scaling)
    
            x = [obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step), ... 
                 obj.nodeArray(2).getX + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_X', step), ... 
                 obj.nodeArray(3).getX + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_X', step), ...
                 obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step)];
             
            y = [obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step), ... 
                 obj.nodeArray(2).getY + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Y', step), ... 
                 obj.nodeArray(3).getY + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Y', step), ...
                 obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step)];
            
            z = [obj.nodeArray(1).getZ, obj.nodeArray(2).getZ, ...
                 obj.nodeArray(3).getZ, obj.nodeArray(1).getZ];
             
            pl = line(x,y,z);
        end
        
        function dofs = getDofList(obj)
            dofs([1 3 5]) = obj.nodeArray.getDof('DISPLACEMENT_X'); 
            dofs([2 4 6]) = obj.nodeArray.getDof('DISPLACEMENT_Y');
        end
        
        function vals = getValuesVector(obj, step)
            vals = zeros(1,6);
            
            vals([1 3 5]) = obj.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals([2 4 6]) = obj.nodeArray.getDofValue('DISPLACEMENT_Y',step);
        end
    
    end
    
    methods (Static)
        function ord = drawOrder()
            ord = [1,2,3,1];
        end
        
        function p = stressPoints()
        %STRESSPOINTS returns locations where stresses are evaluated
            p = [1 0 0;0 1 0;0 0 1];
        end
    end
end



