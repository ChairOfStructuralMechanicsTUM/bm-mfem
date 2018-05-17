classdef PlaneStressElement3d9n < QuadrilateralElement
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    methods
        %Constructor
        function planeStressElement3d9n = PlaneStressElement3d9n(id,nodeArray)
            
            requiredPropertyNames = cellstr(["YOUNGS_MODULUS", "POISSON_RATIO", ...
                                             "THICKNESS", "NUMBER_GAUSS_POINT", ...
                                             "DENSITY"]);
                                         
            %define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 9 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end

            %call the super class constructor
            planeStressElement3d9n@QuadrilateralElement(super_args{:});
            planeStressElement3d9n.dofNames = cellstr(["DISPLACEMENT_X", "DISPLACEMENT_Y"]);
        end
    
        %Initialization
        function initialize(planeStressElement3d9n)
            planeStressElement3d9n.lengthX = computeLength(planeStressElement3d9n.nodeArray(1).getCoords, ...
                planeStressElement3d9n.nodeArray(2).getCoords);

            planeStressElement3d9n.lengthY = computeLength(planeStressElement3d9n.nodeArray(1).getCoords, ...
                planeStressElement3d9n.nodeArray(4).getCoords);

            checkConvexity(planeStressElement3d9n);
        end

        function responseDoF = getResponseDofArray(planeStressElement, step)

            responseDoF = zeros(18,1);
            for itNodes = 1:1:4
                nodalDof = planeStressElement.nodeArray(itNodes).getDofArray;
                nodalDof = nodalDof.';

                for itDof = 2:(-1):1
                    responseDoF(3*itNodes-(itDof-1),1) = nodalDof(4-itDof).getValue(step);
                end
            end
        end

        function [N_mat, N, B, J] = computeShapeFunction(planeStressElement3d9n,xi,eta)
            % Shape Function and Derivatives                    
            N = [(1-xi)*(1-eta)*xi*eta/4,-(1+xi)*(1-eta)*xi*eta/4,(1+xi)*(1+eta)*xi*eta/4,-(1-xi)*(1+eta)*xi*eta/4,...
                    -(1-xi^2)*(1-eta)*eta/2,(1+xi)*(1-eta^2)*xi/2,(1-xi^2)*(1+eta)*eta/2,-(1-xi)*(1-eta)*xi*eta/2,(1-xi^2)*(1-eta^2)];  

            N_Diff_Par = [(1/4)*(-1+2*xi)*(-1+eta)*eta, (1/4)*(1+2*xi)*(-1+eta)*eta,(1/4)*(1+2*xi)*eta*(1+eta),...
                (1/4)*(-1+2*xi)*eta*(1+eta),-1*xi*(-1+eta)*eta,(-1/2)*(1+2*xi)*(-1+eta^2),...
                (-1)*xi*eta*(1+eta),(-1/2)*(-1+2*xi)*(-1+eta^2),2*xi*(-1+eta^2);...
                (1/4)*(-1+xi)*xi*(-1+2*eta),(1/4)*xi*(1+xi)*(-1+2*eta),(1/4)*xi*(1+xi)*(1+2*eta),...
                (1/4)*(-1+xi)*xi*(1+2*eta),(-1/2)*(-1+xi^2)*(-1+2*eta),-xi*(1+xi)*eta,...
                (-1/2)*(-1+xi^2)*(1+2*eta),-1*(-1+xi)*xi*eta,2*(-1+xi^2)*eta];



            N_mat = sparse(2,18);
            N_mat(1,1:2:end) = N(:);
            N_mat(2,2:2:end) = N(:);

            % Coordinates of the nodes forming one element 
            ele_coords = zeros(9,2); 
            for i=1:9
                ele_coords(i,1) = planeStressElement3d9n.nodeArray(i).getX;
                ele_coords(i,2) = planeStressElement3d9n.nodeArray(i).getY;
            end

            % Jacobian 
            J = N_Diff_Par * ele_coords;

            N_Diff = J \ N_Diff_Par;

            % Assembling the B Matrix
            B = zeros(3,18);
            B(1,1:2:end) = N_Diff(1,:);     
            B(3,2:2:end) = N_Diff(1,:);
            B(2,2:2:end) = N_Diff(2,:);
            B(3,1:2:end) = N_Diff(2,:);

        end

        function stiffnessMatrix = computeLocalStiffnessMatrix(planeStressElement3d9n)
            EModul = planeStressElement3d9n.getPropertyValue('YOUNGS_MODULUS');
            prxy = planeStressElement3d9n.getPropertyValue('POISSON_RATIO');
            nr_gauss_points = planeStressElement3d9n.getPropertyValue('NUMBER_GAUSS_POINT');
            thickness = planeStressElement3d9n.getPropertyValue('THICKNESS');

            % Moment-Curvature Equations
            D = [1    prxy    0; prxy     1   0; 0    0   (1-prxy)/2];
            % Material Matrix D
            D = D * EModul * thickness / (1-prxy^2);

            [w,g] = returnGaussPoint(nr_gauss_points);
            stiffnessMatrix = sparse(18,18);
            for xi = 1 : nr_gauss_points
                for eta = 1 : nr_gauss_points
                    [~, ~, B, J] = computeShapeFunction(planeStressElement3d9n,g(xi),g(eta));

                    stiffnessMatrix = stiffnessMatrix + ...
                        B' * D *B * det(J) * w(xi) * w(eta);
                end
            end
        end
        
        function massMatrix = computeLocalMassMatrix(planeStressElement3d9n)
            %Formulation of the Massmatrix based on the Shape Functions
            
            density = planeStressElement3d9n.getPropertyValue('DENSITY');
            thickness = planeStressElement3d9n.getPropertyValue('THICKNESS');
            nr_gauss_points = planeStressElement3d9n.getPropertyValue('NUMBER_GAUSS_POINT');
            [w,g] = returnGaussPoint(nr_gauss_points);

            dens_mat = sparse(2,2);
            dens_mat(1,1) = density*thickness; 
            dens_mat(2,2) = dens_mat(1,1); 

            massMatrix = sparse(18,18);
            for xi = 1 : nr_gauss_points
                for eta = 1 : nr_gauss_points
                    [N_mat,~,~,J] = computeShapeFunction(planeStressElement3d9n,g(xi),g(eta));
                    
                    massMatrix = massMatrix + N_mat' * dens_mat * N_mat *det(J) * w(xi) * w(eta);
                end
            end
        end        
        
        function dampingMatrix = computeLocalDampingMatrix(e)
            eProperties = e.getProperties;
            dampingMatrix = sparse(18,18);

            if (eProperties.hasValue('RAYLEIGH_ALPHA'))
                alpha = eProperties.getValue('RAYLEIGH_ALPHA');
                dampingMatrix = dampingMatrix + alpha * element.computeLocalMassMatrix;
            end

            if (eProperties.hasValue('RAYLEIGH_BETA'))
                beta = eProperties.getValue('RAYLEIGH_BETA');
                dampingMatrix = dampingMatrix + beta * element.computeLocalStiffnessMatrix;
            end
        end
        
        function pl = drawDeformed(planeStressElement3d9n, step, scaling)
    
            x = [planeStressElement3d9n.nodeArray(1).getX + scaling * planeStressElement3d9n.nodeArray(1).getDofValue('DISPLACEMENT_X', step), ... 
                 planeStressElement3d9n.nodeArray(2).getX + scaling * planeStressElement3d9n.nodeArray(2).getDofValue('DISPLACEMENT_X', step), ... 
                 planeStressElement3d9n.nodeArray(3).getX + scaling * planeStressElement3d9n.nodeArray(3).getDofValue('DISPLACEMENT_X', step), ...
                 planeStressElement3d9n.nodeArray(4).getX + scaling * planeStressElement3d9n.nodeArray(4).getDofValue('DISPLACEMENT_X', step), ...
                 planeStressElement3d9n.nodeArray(1).getX + scaling * planeStressElement3d9n.nodeArray(1).getDofValue('DISPLACEMENT_X', step)];
             
            y = [planeStressElement3d9n.nodeArray(1).getY + scaling * planeStressElement3d9n.nodeArray(1).getDofValue('DISPLACEMENT_Y', step), ... 
                 planeStressElement3d9n.nodeArray(2).getY + scaling * planeStressElement3d9n.nodeArray(2).getDofValue('DISPLACEMENT_Y', step), ... 
                 planeStressElement3d9n.nodeArray(3).getY + scaling * planeStressElement3d9n.nodeArray(3).getDofValue('DISPLACEMENT_Y', step), ...
                 planeStressElement3d9n.nodeArray(4).getY + scaling * planeStressElement3d9n.nodeArray(4).getDofValue('DISPLACEMENT_Y', step), ...
                 planeStressElement3d9n.nodeArray(1).getY + scaling * planeStressElement3d9n.nodeArray(1).getDofValue('DISPLACEMENT_Y', step)];
            
            z = [planeStressElement3d9n.nodeArray(1).getZ, planeStressElement3d9n.nodeArray(2).getZ, ...
                 planeStressElement3d9n.nodeArray(3).getZ, planeStressElement3d9n.nodeArray(4).getZ, ...
                 planeStressElement3d9n.nodeArray(1).getZ];
             
            pl = line(x,y,z);
        end
        
        function dofs = getDofList(element)
            dofs([1 3 5 7 9 11 13 15 17]) = element.nodeArray.getDof('DISPLACEMENT_X'); 
           
            dofs([2 4 6 8 10 12 14 16 18]) = element.nodeArray.getDof('DISPLACEMENT_Y');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,18);
            
            vals([1 3 5 7 9 11 13 15 17]) = element.nodeArray.getDofValue('DISPLACEMENT_X',step);

            vals([2 4 6 8 10 12 14 16 18]) = element.nodeArray.getDofValue('DISPLACEMENT_Y',step);
        end
        
        %Computation of the Internal Element Stresses
        function [stressValue, element_nodes] = computeElementStress(elementArray,nodeArray, step)

            element_nodes = zeros(length(elementArray),9);
            stressValue = zeros(3,length(nodeArray));

            
            for i = 1:length(elementArray)

                element_nodes(i,1:9) = elementArray(i).getNodes.getId();
%                 stressPoints = [-1 -1;1 -1;1 1;-1 1;-1 0;0 1;1 0;0 -1;0 0];
                stressPoints = [-1 -1;1 -1;1 1;-1 1;0 -1;1 0;0 1;-1 0;0 0];
                EModul = elementArray(i).getPropertyValue('YOUNGS_MODULUS');
                prxy = elementArray(i).getPropertyValue('POISSON_RATIO');
                % Moment-Curvature Equations
                D = [1    prxy    0; prxy     1   0; 0    0   (1-prxy)/2];
                % Material Matrix D
                D = D * EModul / (1-prxy^2);
                
                for j = 1:9
                    [~, ~, B, ~] = computeShapeFunction(elementArray(i),stressPoints(j,1),stressPoints(j,2));
                    displacement_e = getValuesVector(elementArray(i),1);
                    displacement_e = displacement_e';
                    strain_e = B * displacement_e;
                    stress_e = D * strain_e;
                    
                    % elementwise stress calculation
                    sigma_xx(i,j) = stress_e(1);
                    sigma_yy(i,j) = stress_e(2);
                    sigma_xy(i,j) = stress_e(3);

                end
            end
            
            for k = 1 : length(nodeArray)
                [I,J] = find(element_nodes == k);
                
                sum_sigma_xx = 0;
                sum_sigma_yy = 0;
                sum_sigma_xy = 0;
                
                for l = 1: length(I)
                    sum_sigma_xx = sum_sigma_xx + sigma_xx(I(l),J(l));
                    sum_sigma_yy = sum_sigma_yy + sigma_yy(I(l),J(l));
                    sum_sigma_xy = sum_sigma_xy + sigma_xy(I(l),J(l));
                end
                
                smooth_sigma_xx(k) = sum_sigma_xx/length(I);
                smooth_sigma_yy(k) = sum_sigma_yy/length(I);
                smooth_sigma_xy(k) = sum_sigma_xy/length(I);
            end
            
        stressValue(1,:) = smooth_sigma_xx;
        stressValue(2,:) = smooth_sigma_yy;
        stressValue(3,:) = smooth_sigma_xy;
        
        end   
    end
    
    methods (Static)
        function ord = drawOrder()

            ord = [1,5,2,6,3,7,4,8,1]; % Note: stress value at node 9 is ignored
        end
    end
end



