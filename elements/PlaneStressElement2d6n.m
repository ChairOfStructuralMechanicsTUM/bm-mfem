classdef PlaneStressElement2d6n < TriangularElement
    %PLANESTRESSELEMENT2D6N A 6-node triangular plane stress element
    %   Details on the implementation can be found in ???
    %
    % Node numbering:
    %       3
    %      / \
    %     6   5
    %    /     \
    %   1 - 4 - 2
    
    properties (Access = private)
    end
    
    methods
        %Constructor
        function obj = PlaneStressElement2d6n(id,nodeArray)
            
            requiredPropertyNames = cellstr(["YOUNGS_MODULUS", "POISSON_RATIO", ...
                "THICKNESS", "NUMBER_GAUSS_POINT", ...
                "DENSITY"]);
            
            %define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 6 && isa(nodeArray,'Node'))
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
            %check, if node numbering is correct
            n = obj.nodeArray;
            if ~ (isOnLineBetweenTwoPoints(n(1).getCoords, n(2).getCoords, n(4).getCoords) ...
                    && isOnLineBetweenTwoPoints(n(2).getCoords, n(3).getCoords, n(5).getCoords) ...
                    && isOnLineBetweenTwoPoints(n(1).getCoords, n(3).getCoords, n(6).getCoords))
                msg = [class(obj),': Invalid node numbering ', ...
                    'in element ',num2str(obj.getId), '.'];
                e = MException('MATLAB:bm_mfem:invalidNodeNumbering',msg);
                throw(e);
            end
        end
        
        function [N_mat, N, B, J] = computeShapeFunction(obj,zeta)
            
            % Triangle Coords (Substituting zeta(3) = 1-xi-eta to form independent set)
            xi= zeta(1);
            eta = zeta(2);
            
            % Shape Function and Derivatives
            N = [xi*(2*xi-1) eta*(2*eta-1) (1-xi-eta)*(1-2*xi-2*eta) 4*xi*eta 4*(1-xi-eta)*eta 4*(1-eta-xi)*xi];
            
            N_Diff_Par = [4*xi-1 0 4*xi+4*eta-3 4*eta -4*eta -4*(2*xi+eta-1)
                0 4*eta-1 4*xi+4*eta-3 4*xi -4*(2*eta+xi-1) -4*xi];
            
            N_mat = sparse(2,12);
            N_mat(1,1:2:end) = N(:);
            N_mat(2,2:2:end) = N(:);
            
            % Coordinates of the nodes forming one element
            ele_coords = zeros(6,2);
            for i=1:6
                ele_coords(i,1) = obj.nodeArray(i).getX;
                ele_coords(i,2) = obj.nodeArray(i).getY;
            end
            
            % Jacobian
            J = N_Diff_Par * ele_coords;
            
            N_Diff = J \ N_Diff_Par;
            
            % Assembling the B Matrix
            B = sparse(3,12);
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
            
            stiffnessMatrix = sparse(12,12);
            for i = 1 : abs(nr_gauss_points)
                
                [w,g] = returnGaussPointTrig(nr_gauss_points,i);
                [~, ~, B, J] = computeShapeFunction(obj,g);
                stiffnessMatrix = stiffnessMatrix + ...
                    w * B' * D * B * 0.5 * det(J);
            end
        end
        
        function massMatrix = computeLocalMassMatrix(obj)
            %Formulation of the mass matrix based on the shape functions
            
            density = obj.getPropertyValue('DENSITY');
            thickness = obj.getPropertyValue('THICKNESS');
            nr_gauss_points = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            
            dens_mat = sparse(2,2);
            dens_mat(1,1) = density*thickness;
            dens_mat(2,2) = dens_mat(1,1);
            
            massMatrix = sparse(12,12);
            for i = 1 : nr_gauss_points
                
                [w,g] = returnGaussPointTrig(nr_gauss_points,i);
                [N_mat, ~, ~, J] = computeShapeFunction(obj,g);
                
                massMatrix = massMatrix + N_mat' * dens_mat * N_mat * 0.5 * det(J) * w;
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
            dofs([1 3 5 7 9 11]) = obj.nodeArray.getDof('DISPLACEMENT_X');
            dofs([2 4 6 8 10 12]) = obj.nodeArray.getDof('DISPLACEMENT_Y');
        end
        
        function vals = getValuesVector(obj, step)
            vals = zeros(1,6);
            
            vals([1 3 5 7 9 11]) = obj.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals([2 4 6 8 10 12]) = obj.nodeArray.getDofValue('DISPLACEMENT_Y',step);
        end
        
    end
    
    methods (Static)
        function ord = drawOrder()
            ord = [1,4,2,5,3,6,1];
        end
        
        function p = stressPoints()
            %STRESSPOINTS returns locations where stresses are evaluated
            p = [1 0 0;0 1 0;0 0 1;0.5 0.5 0;0 0.5 0.5;0.5 0 0.5];
        end
    end
end
