classdef QuadrilateralElement2d4n < QuadrilateralElement
    %QUADRILATERALELEMENT2D4N A basic quadrilateral element 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CAUTION SHEAR LOCKING NOT PEREVENTED YET      %
    % Element works only for certain configurations %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    methods
        % Constructor
        function quadrilateralElement2d4n = QuadrilateralElement2d4n(id,nodeArray)
            
            requiredPropertyNames = cellstr(["YOUNGS_MODULUS", "POISSON_RATIO", ...
                "NUMBER_GAUSS_POINT","DENSITY"]);
            
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 4 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            % call the super class contructor
            quadrilateralElement2d4n@QuadrilateralElement(super_args{:});
            quadrilateralElement2d4n.dofNames = cellstr(["DISPLACEMENT_X", ...
                "DISPLACEMENT_Y"]);
        end
        
        %Initialization
        function initialize(obj)
            obj.lengthX = computeLength(obj.nodeArray(1).getCoords, ...
                obj.nodeArray(2).getCoords);
            
            obj.lengthY = computeLength(obj.nodeArray(1).getCoords, ...
                obj.nodeArray(4).getCoords);
            
            if ~checkConvexity(obj)
                msg = ['QuaadrilateralElement2d4n: Element ', ...
                    num2str(obj.getId), ' is not convex.'];
                e = MException('MATLAB:bm_mfem:elementNotConvex',msg);
                throw(e);
            end
        end
        
        function [N_mat, N, Be, J] = computeShapeFunction(obj,xi,eta)
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
                ele_coords(i,1) = obj.nodeArray(i).getX;
                ele_coords(i,2) = obj.nodeArray(i).getY;
            end
            
            % Jacobian
            J = N_Diff_Par * ele_coords;
            
            % Calculation of B-Matrix
            B=J\N_Diff_Par;%/Jdet;
            Bx=B(1,1:4);
            By=B(2,1:4);
            
            Be=[Bx(1),0,Bx(2),0,Bx(3),0,Bx(4),0;
                0,By(1),0,By(2),0,By(3),0,By(4);
                By(1),Bx(1),By(2),Bx(2),By(3),Bx(3),By(4),Bx(4)];
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(obj)
            EModul = obj.getPropertyValue('YOUNGS_MODULUS');
            PoissonRatio = obj.getPropertyValue('POISSON_RATIO');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            % Calculate Materialmatrix
            Emat = EModul/(1-PoissonRatio^2)*[1 PoissonRatio 0; PoissonRatio 1 0; 0 0 (1 - PoissonRatio)/2];
            stiffnessMatrix=sparse(8,8);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [~, ~, B, J] = computeShapeFunction(obj, xi, eta);
                    stiffnessMatrix=stiffnessMatrix+(w(i)*w(j)*det(J)*transpose(B)*(Emat*B));
                end
            end
        end
        
        function massMatrix = computeLocalMassMatrix(obj)
            roh = obj.getPropertyValue('DENSITY');
            p = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            massMatrix=sparse(8,8);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [N_mat, ~, ~, J] = computeShapeFunction(obj,xi,eta);
                    massMatrix=massMatrix + (w(i)*w(j)*roh*transpose(N_mat)*N_mat*det(J));
                end
            end
        end
        
        
        function dampingMatrix = computeLocalDampingMatrix(obj)
            eProperties = obj.getProperties;
            dampingMatrix = sparse(8,8);
            
            if (eProperties.hasValue('RAYLEIGH_ALPHA'))
                alpha = eProperties.getValue('RAYLEIGH_ALPHA');
                dampingMatrix = dampingMatrix + alpha * obj.computeLocalMassMatrix;
            end
            
            if (eProperties.hasValue('RAYLEIGH_BETA'))
                beta = eProperties.getValue('RAYLEIGH_BETA');
                dampingMatrix = dampingMatrix + beta * obj.computeLocalStiffnessMatrix;
            end
            
        end
        
        function dofs = getDofList(obj)
            dofs([1 3 5 7]) = obj.nodeArray.getDof('DISPLACEMENT_X');
            dofs([2 4 6 8]) = obj.nodeArray.getDof('DISPLACEMENT_Y');
        end
        
        function vals = getValuesVector(obj, step)
            vals = zeros(1,8);
            
            vals([1 3 5 7]) = obj.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals([2 4 6 8]) = obj.nodeArray.getDofValue('DISPLACEMENT_Y',step);
        end
        
        function vals = getFirstDerivativesVector(obj, step)
            vals = zeros(1,8);
            
            [~, vals([1 3 5 7]), ~] = obj.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, vals([2 4 6 8]), ~] = obj.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(obj, step)
            vals = zeros(1,8);
            
            [~, ~, vals([1 3 5 7])] = obj.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, ~, vals([2 4 6 8])] = obj.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
        end      
        
        function F = computeLocalForceVector(obj)
            F = zeros(1,8);
        end
        
    end
end
