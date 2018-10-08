classdef PlaneStressElement3d8n < QuadrilateralElement
    %PLANESTRESSELEMENT3D8N A 8-node quadrilateral plain stress element
    %   Details on the implementation can be found in ???
    %
    % Node numbering:
    %   4 - 7 - 3
    %   |       |
    %   8       6
    %   |       |
    %   1 - 5 - 2
    
    properties (Access = private)
    end
    
    methods
        %Constructor
        function obj = PlaneStressElement3d8n(id,nodeArray)
            
            requiredPropertyNames = cellstr(["YOUNGS_MODULUS", "POISSON_RATIO", ...
                                             "THICKNESS", "NUMBER_GAUSS_POINT", ...
                                             "DENSITY"]);
                                         
            %define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 8 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end

            %call the super class constructor
            obj@QuadrilateralElement(super_args{:});
            obj.dofNames = cellstr(["DISPLACEMENT_X", "DISPLACEMENT_Y"]);
        end
    
        %Initialization
        function initialize(obj)
            obj.lengthX = computeLength(obj.nodeArray(1).getCoords, ...
                obj.nodeArray(2).getCoords);

            obj.lengthY = computeLength(obj.nodeArray(1).getCoords, ...
                obj.nodeArray(4).getCoords);

            if ~ checkConvexity(obj)
                msg = [class(obj), ': Element ', ...
                    num2str(obj.getId), ' is not convex.'];
                e = MException('MATLAB:bm_mfem:elementNotConvex',msg);
                throw(e);
            end
            
            %check, if node numbering is correct
            n = obj.nodeArray;
            if ~ (isOnLineBetweenTwoPoints(n(1).getCoords, n(2).getCoords, n(5).getCoords) ...
                    && isOnLineBetweenTwoPoints(n(2).getCoords, n(3).getCoords, n(6).getCoords) ...
                    && isOnLineBetweenTwoPoints(n(3).getCoords, n(4).getCoords, n(7).getCoords) ...
                    && isOnLineBetweenTwoPoints(n(1).getCoords, n(4).getCoords, n(8).getCoords))
                msg = [class(obj),': Invalid node numbering ', ...
                    'in element ',num2str(obj.getId), '.'];
                e = MException('MATLAB:bm_mfem:invalidNodeNumbering',msg);
                throw(e);
            end
        end

        function [N_mat, N, B, J] = computeShapeFunction(obj,xi,eta)
            % Shape Function and Derivatives                    
            N = [-(1-xi)*(1-eta)*(1+xi+eta)/4,-(1+xi)*(1-eta)*(1-xi+eta)/4,-(1+xi)*(1+eta)*(1-xi-eta)/4,-(1-xi)*(1+eta)*(1+xi-eta)/4,...
                    (1-xi^2)*(1-eta)/2,(1+xi)*(1-eta^2)/2,(1-xi^2)*(1+eta)/2,(1-xi)*(1-eta^2)/2];

            N_Diff_Par = [(1-eta)*(2*xi+eta)/4,-(1-eta)*(eta-2*xi)/4,(1+eta)*(2*xi+eta)/4,-(1+eta)*(eta-2*xi)/4,...
                -xi*(1-eta),(1-eta^2)/2,-xi*(1+eta),-(1-eta^2)/2;...
                (1-xi)*(xi+2*eta)/4,(1+xi)*(2*eta-xi)/4,(1+xi)*(xi+2*eta)/4,(1-xi)*(2*eta-xi)/4,...
                -(1-xi^2)/2,-eta*(1+xi),(1-xi^2)/2,-eta*(1-xi)];
       
            N_mat = sparse(2,16);
            N_mat(1,1:2:end) = N(:);
            N_mat(2,2:2:end) = N(:);

            % Coordinates of the nodes forming one element 
            ele_coords = zeros(8,2); 
            for i=1:8
                ele_coords(i,1) = obj.nodeArray(i).getX;
                ele_coords(i,2) = obj.nodeArray(i).getY;
            end

            % Jacobian 
            J = N_Diff_Par * ele_coords;
            N_Diff = J \ N_Diff_Par;

            % Assembling the B Matrix
            B = zeros(3,16);
            B(1,1:2:end) = N_Diff(1,:);     
            B(3,2:2:end) = N_Diff(1,:);
            B(2,2:2:end) = N_Diff(2,:);
            B(3,1:2:end) = N_Diff(2,:);

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

            [w,g] = returnGaussPoint(nr_gauss_points);
            stiffnessMatrix = sparse(16,16);
            for xi = 1 : nr_gauss_points
                for eta = 1 : nr_gauss_points
                    [~, ~, B, J] = computeShapeFunction(obj,g(xi),g(eta));

                    stiffnessMatrix = stiffnessMatrix + ...
                        B' * D * B * det(J) * w(xi) * w(eta);
                end
            end
        end
        
        function massMatrix = computeLocalMassMatrix(obj)
            %Formulation of the Massmatrix based on the Shape Functions
            
            density = obj.getPropertyValue('DENSITY');
            thickness = obj.getPropertyValue('THICKNESS');
            nr_gauss_points = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            [w,g] = returnGaussPoint(nr_gauss_points);

            dens_mat = sparse(2,2);
            dens_mat(1,1) = density*thickness; 
            dens_mat(2,2) = dens_mat(1,1); 

            massMatrix = sparse(16,16);
            for xi = 1 : nr_gauss_points
                for eta = 1 : nr_gauss_points
                    [N_mat,~,~,J] = computeShapeFunction(obj,g(xi),g(eta));
                    
                    massMatrix = massMatrix + N_mat' * dens_mat * N_mat *det(J) * w(xi) * w(eta);
                end
            end
        end
        
        function pl = drawDeformed(obj, step, scaling)
    
            x = [obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step), ... 
                 obj.nodeArray(2).getX + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_X', step), ... 
                 obj.nodeArray(3).getX + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_X', step), ...
                 obj.nodeArray(4).getX + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_X', step), ...
                 obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step)];
             
            y = [obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step), ... 
                 obj.nodeArray(2).getY + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Y', step), ... 
                 obj.nodeArray(3).getY + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Y', step), ...
                 obj.nodeArray(4).getY + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Y', step), ...
                 obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step)];
            
            z = [obj.nodeArray(1).getZ, obj.nodeArray(2).getZ, ...
                 obj.nodeArray(3).getZ, obj.nodeArray(4).getZ, ...
                 obj.nodeArray(1).getZ];
             
            pl = line(x,y,z);
        end
        
        function dofs = getDofList(obj)
            dofs([1 3 5 7 9 11 13 15]) = obj.nodeArray.getDof('DISPLACEMENT_X'); 
            dofs([2 4 6 8 10 12 14 16]) = obj.nodeArray.getDof('DISPLACEMENT_Y');
        end
        
        function vals = getValuesVector(obj, step)
            vals = zeros(1,16);
            
            vals([1 3 5 7 9 11 13 15]) = obj.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals([2 4 6 8 10 12 14 16]) = obj.nodeArray.getDofValue('DISPLACEMENT_Y',step);
        end
        
    end
    
    methods (Static)
        function ord = drawOrder()
            ord = [1,5,2,6,3,7,4,8,1];
        end
        
        function p = stressPoints()
            %STRESSPOINTS returns locations where stresses are evaluated
            p = [-1 -1;1 -1;1 1;-1 1;0 -1;1 0;0 1;-1 0];
        end
    end
end



