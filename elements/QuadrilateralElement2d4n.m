classdef QuadrilateralElement2d4n < QuadrilateralElement
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
                "NUMBER_GAUSS_POINT","DENSITY"]);  %Erstellen von Zeilenvektor mit 4 Schrifteinträgen
            
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
        function initialize(quadrilateralElement2d4n)
            quadrilateralElement2d4n.lengthX = computeLength(quadrilateralElement2d4n.nodeArray(1).getCoords, ...
                quadrilateralElement2d4n.nodeArray(2).getCoords);
            
            quadrilateralElement2d4n.lengthY = computeLength(quadrilateralElement2d4n.nodeArray(1).getCoords, ...
                quadrilateralElement2d4n.nodeArray(4).getCoords);
            
            checkConvexity(quadrilateralElement2d4n);
        end
        
        function [N_mat, N, Be, J] = computeShapeFunction(quadrilateralElement2d4n,xi,eta)
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
                ele_coords(i,1) = quadrilateralElement2d4n.nodeArray(i).getX;
                ele_coords(i,2) = quadrilateralElement2d4n.nodeArray(i).getY;
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
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(quadrilateralElement2d4n)
            EModul = quadrilateralElement2d4n.getPropertyValue('YOUNGS_MODULUS');
            PoissonRatio = quadrilateralElement2d4n.getPropertyValue('POISSON_RATIO');
            p = quadrilateralElement2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            % Calculate Materialmatrix
            Emat = EModul/(1-PoissonRatio^2)*[1 PoissonRatio 0; PoissonRatio 1 0; 0 0 (1 - PoissonRatio)/2];
            stiffnessMatrix=zeros(8,8);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [~, ~, B, J] = computeShapeFunction(quadrilateralElement2d4n, xi, eta);
                    stiffnessMatrix=stiffnessMatrix+(w(i)*w(j)*det(J)*transpose(B)*(Emat*B));
                end
            end
        end
        
        function massMatrix = computeLocalMassMatrix(quadrilateralElement2d4n)
            roh = quadrilateralElement2d4n.getPropertyValue('DENSITY');
            p = quadrilateralElement2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            massMatrix=zeros(8,8);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [N_mat, ~, ~, J] = computeShapeFunction(quadrilateralElement2d4n,xi,eta);
                    massMatrix=massMatrix + (w(i)*w(j)*roh*transpose(N_mat)*N_mat*det(J));
                end
            end
        end
        
        
        function dampingMatrix = computeLocalDampingMatrix(e)
            eProperties = e.getProperties;
            dampingMatrix = sparse(8,8);
            
            if (eProperties.hasValue('RAYLEIGH_ALPHA'))
                alpha = eProperties.getValue('RAYLEIGH_ALPHA');
                dampingMatrix = dampingMatrix + alpha * element.computeLocalMassMatrix;
            end
            
            if (eProperties.hasValue('RAYLEIGH_BETA'))
                beta = eProperties.getValue('RAYLEIGH_BETA');
                dampingMatrix = dampingMatrix + beta * element.computeLocalStiffnessMatrix;
            end
            
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
        
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,8);
            
            [~, vals([1 3 5 7]), ~] = element.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, vals([2 4 6 8]), ~] = element.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,8);
            
            [~, ~, vals([1 3 5 7])] = element.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, ~, vals([2 4 6 8])] = element.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
        end      
        
        function F = computeLocalForceVector(quadrilateralElement)
            F = zeros(1,8);
        end
        
    end
end
