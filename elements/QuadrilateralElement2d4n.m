classdef QuadrilateralElement2d4n < QuadrilateralElement
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CAUTION SHEAR LOCKING NOT PEREVENTED YET %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    methods
        % Constructor
        function quadrilateralElement2d4n = QuadrilateralElement2d4n(id,nodeArray)
            
            requiredPropertyNames = cellstr(["YOUNGS_MODULUS", "POISSON_RATIO", ...
                "NUMBER_GAUSS_POINT","DENSITY","SHEAR_CORRECTION_FACTOR"]);
            
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
        
        function responseDoF = getResponseDofArray(element, step)
            
            responseDoF = zeros(12,1);
            for itNodes = 1:1:4
                nodalDof = element.nodeArray(itNodes).getDofArray;
                nodalDof = nodalDof.';
                
                for itDof = 3:(-1):1
                    responseDoF(3*itNodes-(itDof-1),1) = nodalDof(4-itDof).getValue(step);
                end
            end
        end
        
        function [N_mat, N, B_b, B_s, J] = computeShapeFunction(quadrilateralElement2d4n,xi,eta)
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
            
            % Calculation of B_bending Matrix
            %B_b = sparse(3,8);
            B_b=[Bx(1),0,Bx(2),0,Bx(3),0,Bx(4),0;
                0,By(1),0,By(2),0,By(3),0,By(4);
                By(1),Bx(1),By(2),Bx(2),By(3),Bx(3),By(4),Bx(4)];
            
            % Calculation of B_shear Matrix
            %B_s = sparse(2,8);
            %If 3 DOF Element: B_s=[Bx(1),N(1),0,Bx(2),N(2),0,Bx(3),N(3),0,Bx(4),N(4),0; By(1),N(1),0,By(2),N(2),0,By(3),N(3),0,By(4),N(4),0];
            B_s=[N(1),0,N(2),0,N(3),0,N(4),0; 
                N(1),0,N(2),0,N(3),0,N(4),0];
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(quadrilateralElement2d4n)
            EModul = quadrilateralElement2d4n.getPropertyValue('YOUNGS_MODULUS');
            prxy = quadrilateralElement2d4n.getPropertyValue('POISSON_RATIO');
            p = quadrilateralElement2d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            alpha_shear = quadrilateralElement2d4n.getPropertyValue('SHEAR_CORRECTION_FACTOR');
            
            % Calculate Shear Modulus
            GModul = EModul/(2*(1+prxy));
            % Moment-Curvature Equations
            D_b = [1    prxy    0; prxy     1   0; 0    0   (1-prxy)/2];
            % Material Bending Matrix D_b
            D_b = D_b * (EModul) / ((1-prxy^2));
            % Material Shear Matrix D_s
            D_s = eye(2) * alpha_shear * GModul; 
            
            [w,g]=returnGaussPoint(p);
            
            stiffnessMatrix=sparse(8,8);
           
            if (quadrilateralElement2d4n.getProperties.hasValue('FULL_INTEGRATION')) && (quadrilateralElement2d4n.getPropertyValue('FULL_INTEGRATION'))
                for xi = 1 : p
                    for eta = 1 : p
                        [~, ~,B_b, B_s, J] = computeShapeFunction(quadrilateralElement2d4n,g(xi),g(eta));

                        stiffnessMatrix = stiffnessMatrix +             ...
                            B_b' * D_b * B_b * det(J) * w(xi) * w(eta) + ...
                            B_s' * D_s * B_s * det(J) * w(xi) * w(eta);                       
                    end
                end
                
            else
                for xi = 1 : p
                    for eta = 1 : p
                        [~, ~,B_b, ~, J] = computeShapeFunction(quadrilateralElement2d4n,g(xi),g(eta));

                        stiffnessMatrix = stiffnessMatrix +             ...
                            B_b' * D_b * B_b * det(J) * w(xi) * w(eta);
                    end
                end
                [w,g] = returnGaussPoint(1);
                for xi = 1 : 1
                    for eta = 1 : 1
                        [~, ~,~, B_s, J] = computeShapeFunction(quadrilateralElement2d4n,g(xi),g(eta));
                        
                        stiffnessMatrix = stiffnessMatrix +             ...
                            B_s' * D_s * B_s * det(J) * w(xi) * w(eta);                       
                    end
                end
            end
        end
        
        %The two Versions below work for some Elements, but lock for slender structures
        
        function stiffnessMatrix = computeLocalStiffnessMatrix_Option2(quadrilateralElement2d4n)
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
                    
                    [~, ~,Bx, By, J] = computeShapeFunction(quadrilateralElement2d4n, xi, eta);
                    
                    Be=[Bx(1),0,Bx(2),0,Bx(3),0,Bx(4),0;
                        0,By(1),0,By(2),0,By(3),0,By(4);
                        By(1),Bx(1),By(2),Bx(2),By(3),Bx(3),By(4),Bx(4)];
                    
                    stiffnessMatrix=stiffnessMatrix+(w(i)*w(j)*det(J)*transpose(Be)*(Emat*Be));
                    
                    
                end
            end
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix_Option3(quadrilateralElement2d4n,a,b)
            EModul = quadrilateralElement2d4n.getPropertyValue('YOUNGS_MODULUS');
            PoissonRatio = quadrilateralElement2d4n.getPropertyValue('POISSON_RATIO');
            s=a/b;
            width=1;
            
            k1 = (1 + PoissonRatio)*s;
            k2 = (1 - 3*PoissonRatio)*s;
            k3 = 2 + (1 - PoissonRatio)*s^2;
            k4 = 2*s^2 + (1 - PoissonRatio);
            k5 = (1 - PoissonRatio)*s^2 - 4;
            k6 = (1 - PoissonRatio)*s^2 - 1;
            k7 = 4*s^2 - (1 - PoissonRatio);
            k8 = s^2 - (1 -PoissonRatio);
            
            stiffnessMatrix = EModul*width/(24*s*(1-PoissonRatio^2)) * [  4*k3  3*k1  2*k5 -3*k2 -2*k3 -3*k1 -4*k6  3*k2;
                3*k1  4*k4  3*k2  4*k8 -3*k1 -2*k4 -3*k2 -2*k7;
                2*k5  3*k2  4*k3 -3*k1 -4*k6 -3*k2 -2*k3  3*k1;
                -3*k2  4*k8 -3*k1  4*k4  3*k2 -2*k7  3*k1 -2*k4;
                -2*k3 -3*k1 -4*k6  3*k2  4*k3  3*k1  2*k5 -3*k2;
                -3*k1 -2*k4 -3*k2 -2*k7  3*k1  4*k4  3*k2  4*k8;
                -4*k6 -3*k2 -2*k3  3*k1  2*k5  3*k2  4*k3 -3*k1;
                3*k2 -2*k7  3*k1 -2*k4 -3*k2  4*k8 -3*k1  4*k4];
            
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
                    
                    [N_mat, ~, ~, ~, J] = computeShapeFunction(quadrilateralElement2d4n,xi,eta);
                    
                    massMatrix=massMatrix + (w(i)*w(j)*roh*transpose(N_mat)*N_mat*det(J));
                    
                end
            end
        end
        
        
        %         function dampingMatrix = computeLocalDampingMatrix(e)
        %             eProperties = e.getProperties;
        %             dampingMatrix = sparse(12,12);
        %
        %             if (eProperties.hasValue('RAYLEIGH_ALPHA'))
        %                 alpha = eProperties.getValue('RAYLEIGH_ALPHA');
        %                 dampingMatrix = dampingMatrix + alpha * element.computeLocalMassMatrix;
        %             end
        %
        %             if (eProperties.hasValue('RAYLEIGH_BETA'))
        %                 beta = eProperties.getValue('RAYLEIGH_BETA');
        %                 dampingMatrix = dampingMatrix + beta * element.computeLocalStiffnessMatrix;
        %             end
        %
        %         end
        
        
        function dofs = getDofList(element)
            dofs([1 3 5 7]) = element.nodeArray.getDof('DISPLACEMENT_X');
            dofs([2 4 6 8]) = element.nodeArray.getDof('DISPLACEMENT_Y');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,8);
            
            vals([1 3 5 7]) = element.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals([2 4 6 8]) = element.nodeArray.getDofValue('DISPLACEMENT_Y',step);
        end
        
        %         function vals = getFirstDerivativesVector(element, step)
        %             vals = zeros(1,12);
        %
        %             [~, vals([1 4 7 10]), ~] = element.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
        %             [~, vals([2 5 8 11]), ~] = element.nodeArray.getDof('ROTATION_X').getAllValues(step);
        %             [~, vals([3 6 9 12]), ~] = element.nodeArray.getDof('ROTATION_Y').getAllValues(step);
        %         end
        
        %         function vals = getSecondDerivativesVector(element, step)
        %             vals = zeros(1,12);
        %
        %             [~, ~, vals([1 4 7 10])] = element.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
        %             [~, ~, vals([2 5 8 11])] = element.nodeArray.getDof('ROTATION_X').getAllValues(step);
        %             [~, ~, vals([3 6 9 12])] = element.nodeArray.getDof('ROTATION_Y').getAllValues(step);
        %         end
        
        
    end
end
