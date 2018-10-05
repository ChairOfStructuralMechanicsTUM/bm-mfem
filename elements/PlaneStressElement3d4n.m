classdef PlaneStressElement3d4n < QuadrilateralElement
    %PLANESTRESSELEMENT3D4N A 4-node quadrilateral plane stress element
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    methods
        %Constructor
        function obj = PlaneStressElement3d4n(id,nodeArray)
            
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
            obj@QuadrilateralElement(super_args{:});
            obj.dofNames = cellstr(["DISPLACEMENT_X", "DISPLACEMENT_Y"]);
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

        function responseDoF = getResponseDofArray(obj, step)

            responseDoF = zeros(8,1);
            for itNodes = 1:1:4
                nodalDof = obj.nodeArray(itNodes).getDofArray;
                nodalDof = nodalDof.';

                for itDof = 2:(-1):1
                    responseDoF(3*itNodes-(itDof-1),1) = nodalDof(4-itDof).getValue(step);
                end
            end
        end

        function [N_mat, N, B, J] = computeShapeFunction(obj,xi,eta)
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

            N_Diff = J \ N_Diff_Par;

            % Assembling the B Matrix
            B = sparse(3,8);
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
            stiffnessMatrix = sparse(8,8);
            for xi = 1 : nr_gauss_points
                for eta = 1 : nr_gauss_points
                    [~, ~, B, J] = computeShapeFunction(obj,g(xi),g(eta));

                    stiffnessMatrix = stiffnessMatrix + ...
                        B' * D *B * det(J) * w(xi) * w(eta);
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

            massMatrix = sparse(8,8);
            for xi = 1 : nr_gauss_points
                for eta = 1 : nr_gauss_points
                    [N_mat,~,~,J] = computeShapeFunction(obj,g(xi),g(eta));
                    
                    massMatrix = massMatrix + N_mat' * dens_mat * N_mat *det(J) * w(xi) * w(eta);
                end
            end
        end        
        
        function dampingMatrix = computeLocalDampingMatrix(obj)
            eProperties = obj.getProperties;
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
        
        function f = computeLocalForceVector(obj)
            f = zeros(1,8);
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
            [~, vals([1 3 5 7]), ~] = obj.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(obj, step)
            vals = zeros(1,8);
            
            [~, ~, vals([1 3 5 7])] = obj.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, ~, vals([1 3 5 7])] = obj.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
        end
        
        %Computation of Stresses
        function [stressValue, element_connect] = computeElementStress(obj,nodeArray,step)
        %COMPUTEELEMENTSTRESS stress in the element
        %   [STRESS, EC] = computeElementStress(ELEMENTS,NODES,STEP) returns
        %   element connectivity EC and the matrix STRESS containing:
        %       STRESS(1,:) = stress in xx
        %       STRESS(2,:) = stress in yy
        %       STRESS(3,:) = stress in xy
        %       STRESS(4,:) = first principal stress
        %       STRESS(5,:) = second principal stress
        %       STRESS(6,:) = von-Mises stress
        
            nElements = length(obj);
            nNodes = length(nodeArray);
            
            element_connect = zeros(nElements,4);
            stressValue = zeros(6,nNodes);
            sigma_xx = zeros(nElements,4);
            sigma_yy = zeros(nElements,4);
            sigma_xy = zeros(nElements,4);
            smooth_sigma_xx = zeros(1,nNodes);
            smooth_sigma_yy = zeros(1,nNodes);
            smooth_sigma_xy = zeros(1,nNodes);
            prin_I = zeros(1,nNodes);
            prin_II = zeros(1,nNodes);
            vm_stress = zeros(1,nNodes);
            
            for i = 1:nElements
                element_connect(i,1:4) = obj(i).getNodes.getId();
                stressPoints = [-1 -1;1 -1;1 1;-1 1];            
                EModul = obj(i).getPropertyValue('YOUNGS_MODULUS');
                prxy = obj(i).getPropertyValue('POISSON_RATIO');
                % Moment-Curvature Equations
                D = [1    prxy    0; prxy     1   0; 0    0   (1-prxy)/2];
                % Material Matrix D
                D = D * EModul / (1-prxy^2);
                
                for j = 1:4
                    [~, ~, B, ~] = computeShapeFunction(obj(i),stressPoints(j,1),stressPoints(j,2));
                    displacement_e = getValuesVector(obj(i),step);
                    displacement_e = displacement_e';
                    strain_e = B * displacement_e;
                    stress_e = D * strain_e;
                    
                    % elementwise stress calculation
                    sigma_xx(i,j) = stress_e(1);
                    sigma_yy(i,j) = stress_e(2);
                    sigma_xy(i,j) = stress_e(3);

                end
            end
            
            for k = 1 : nNodes
                [I,J] = find(element_connect == k);
                
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

                stress_ele = [smooth_sigma_xx(k) smooth_sigma_xy(k);
                smooth_sigma_xy(k) smooth_sigma_yy(k)];
                [~,lambda] = eig(stress_ele);

                prin_I(k) = lambda(1,1);
                prin_II(k) = lambda(2,2);
%                 prI_dir(k,:) = vec(:,1);
%                 prII_dir(k,:) = vec(:,2);
                
                vm_stress(k) = sqrt(prin_I(k).^2 + prin_II(k).^2 - prin_I(k) * prin_II(k));
            end
            
            stressValue(1,:) = smooth_sigma_xx;
            stressValue(2,:) = smooth_sigma_yy;
            stressValue(3,:) = smooth_sigma_xy;
            stressValue(4,:) = prin_I;
            stressValue(5,:) = prin_II;
            stressValue(6,:) = vm_stress;
            
        end
    end
    
    methods (Static)
        function ord = drawOrder()
            
            ord = [1,2,3,4,1];
        end
    end
end



