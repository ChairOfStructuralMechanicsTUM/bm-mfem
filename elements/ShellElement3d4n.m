classdef ShellElement3d4n < QuadrilateralElement 
    %SHELLELEMENT3D4N A quadrilateral shell element
    %   Detailed explanation goes here
    
    properties (Access = private)
    end 
    
    methods
        % Constructor
        function shellElement3d4n = ShellElement3d4n(id,nodeArray)
            
            requiredPropertyNames = cellstr(["YOUNGS_MODULUS", "POISSON_RATIO", ...
                                             "THICKNESS", "NUMBER_GAUSS_POINT", ...
                                             "DENSITY", "SHEAR_CORRECTION_FACTOR"]);
                                         
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
            shellElement3d4n@QuadrilateralElement(super_args{:});
            shellElement3d4n.dofNames = cellstr([ ...
                "DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", "ROTATION_X",  "ROTATION_Y",  "ROTATION_Z"]);
        end
        
        %Initialization
        function initialize(shellElement3d4n)
            shellElement3d4n.lengthX = computeLength(shellElement3d4n.nodeArray(1).getCoords, ...
                shellElement3d4n.nodeArray(2).getCoords);
            
            shellElement3d4n.lengthY = computeLength(shellElement3d4n.nodeArray(1).getCoords, ...
                shellElement3d4n.nodeArray(4).getCoords);
            
            checkConvexity(shellElement3d4n);
        end
        
        function responseDoF = getResponseDofArray(shellElement3d4n, step)
           
            responseDoF = zeros(24,1);
            for itNodes = 1:1:6
                nodalDof = shellElement3d4n.nodeArray(itNodes).getDofArray;
                nodalDof = nodalDof.';
                
                for itDof = 3:(-1):1
                    responseDoF(3*itNodes-(itDof-1),1) = nodalDof(4-itDof).getValue(step);
                end
            end
        end
        
        
        
        function [N,  J, B_mem, B_b, B_s] = computeShapeFunction(shellElement3d4n, xi, eta)
            % Shape Function and Derivatives                    
            N = [(1-xi)*(1-eta)/4    (1+xi)*(1-eta)/4    (1+xi)*(1+eta)/4    (1-xi)*(1+eta)/4];  
            
            N_Diff_Par = [-(1-eta)/4    (1-eta)/4   (1+eta)/4   -(1+eta)/4
                          -(1-xi)/4     -(1+xi)/4   (1+xi)/4    (1-xi)/4];

            
            % Coordinates of the nodes forming one element 
            ele_coords = zeros(4,2); 
            for i=1:4
                ele_coords(i,1) = shellElement3d4n.nodeArray(i).getX;
                ele_coords(i,2) = shellElement3d4n.nodeArray(i).getY;
            end
            
            % Jacobian and inverse Jacobian
            J = N_Diff_Par * ele_coords;
            inv_J = inv(J);
            inv_J_TR = [inv_J, zeros(2,2); zeros(2,2), inv_J]; 
            
            N_Diff = J \ N_Diff_Par;
 
            % Blended Displacement Interpolation
%             M1_xi = 1/2*(1-xi);
%             M2_xi = 1/2*(1+xi);
%             
%             M1_eta = 1/2*(1-eta); 
%             M2_eta = 1/2*(1+eta); 
%             
%             N1_xi = 1/2 - 3/4 * xi + (xi^3)/4;
%             N2_xi = 1/4 - xi/4 - (xi^2)/4 + (xi^3)/4; 
%             N3_xi = 1/2 + 3/4 * xi - (xi^3)/4;
%             N4_xi = -1/4 - xi/4 + (xi^2)/4 + (xi^3)/4;
%             
%             N1_eta = 1/2 - 3/4 * eta + (eta^3)/4;
%             N2_eta = 1/4 - eta/4 - (eta^2)/4 + (eta^3)/4; 
%             N3_eta = 1/2 + 3/4 * eta - (eta^3)/4;
%             N4_eta = -1/4 - eta/4 + (eta^2)/4 + (eta^3)/4;
%             
%             MN = [M1_xi*N1_eta    0    M1_xi*N2_eta    0    M2_xi*N1_eta    0    -M2_xi*N2_eta    0    M2_xi*N3_eta    0    -M2_xi*N4_eta    0    M1_xi*N3_eta    0    M1_xi*N4_eta    0
%                        0    M1_eta*N1_xi     0    M1_eta*N2_xi    0    M1_eta*N3_xi    0    M1_eta*N4_xi    0    M2_eta*N3_xi    0    M2_eta*N4_xi    0     M2_eta*N1_xi    0    M2_eta*N2_xi]; 
                   
                   
            MN_Diff_Par=  [ - eta^3/8 + (3*eta)/8 - 1/4,                                 0, - eta^3/8 + eta^2/8 + eta/8 - 1/8,                                       0, eta^3/8 - (3*eta)/8 + 1/4,                                0, - eta^3/8 + eta^2/8 + eta/8 - 1/8,                                        0, - eta^3/8 + (3*eta)/8 + 1/4,                                 0, - eta^3/8 - eta^2/8 + eta/8 + 1/8,                                       0, eta^3/8 - (3*eta)/8 - 1/4,                                0, - eta^3/8 - eta^2/8 + eta/8 + 1/8,                                        0;
                                      -(xi/2 - 1/2)*((3*eta^2)/4 - 3/4),                         0, (xi/2 - 1/2)*(eta/2 - (3*eta^2)/4 + 1/4),                              0, (xi/2 + 1/2)*((3*eta^2)/4 - 3/4),                       0, (xi/2 + 1/2)*(eta/2 - (3*eta^2)/4 + 1/4),                              0, -(xi/2 + 1/2)*((3*eta^2)/4 - 3/4),                         0, -(xi/2 + 1/2)*((3*eta^2)/4 + eta/2 - 1/4),                            0, (xi/2 - 1/2)*((3*eta^2)/4 - 3/4),                       0, -(xi/2 - 1/2)*((3*eta^2)/4 + eta/2 - 1/4),                            0;
                                      0, -(eta/2 - 1/2)*((3*xi^2)/4 - 3/4),                                 0, (eta/2 - 1/2)*(xi/2 - (3*xi^2)/4 + 1/4),                         0, (eta/2 - 1/2)*((3*xi^2)/4 - 3/4),                                 0, -(eta/2 - 1/2)*((3*xi^2)/4 + xi/2 - 1/4),                           0, -(eta/2 + 1/2)*((3*xi^2)/4 - 3/4),                                 0, (eta/2 + 1/2)*((3*xi^2)/4 + xi/2 - 1/4),                         0, (eta/2 + 1/2)*((3*xi^2)/4 - 3/4),                                 0, -(eta/2 + 1/2)*(xi/2 - (3*xi^2)/4 + 1/4);
                                      0, - xi^3/8 + (3*xi)/8 - 1/4,                                        0, - xi^3/8 + xi^2/8 + xi/8 - 1/8,                                0, xi^3/8 - (3*xi)/8 - 1/4,                                        0, - xi^3/8 - xi^2/8 + xi/8 + 1/8,                                 0, - xi^3/8 + (3*xi)/8 + 1/4,                                         0, xi^3/8 + xi^2/8 - xi/8 - 1/8,                                0, xi^3/8 - (3*xi)/8 + 1/4,                                         0, xi^3/8 - xi^2/8 - xi/8 + 1/8]; 
           
                            
           % Displacement Transformation Matrix
           Tr = sparse(16,12); 
           
           Tr(1,1) = 1; 
           Tr(2,2) = 1; 
           Tr(3,3) = (ele_coords(4,2) - ele_coords(1,2))/2; 
           Tr(4,3) = (ele_coords(2,1) - ele_coords(1,1))/2; 
           Tr(5,4) = 1; 
           Tr(6,5) = 1; 
           Tr(7,6) = (ele_coords(3,2) - ele_coords(2,2))/2; 
           Tr(8,6) = (ele_coords(2,1) - ele_coords(1,1))/2; 
           Tr(9,7) = 1; 
           Tr(10,8) = 1; 
           Tr(11,9) = (ele_coords(3,2) - ele_coords(2,2))/2; 
           Tr(12,9) = (ele_coords(3,1) - ele_coords(4,1))/2; 
           Tr(13,10) = 1; 
           Tr(14,11) = 1; 
           Tr(15,12) = (ele_coords(4,2) - ele_coords(1,2))/2;
           Tr(16,12) = (ele_coords(3,1) - ele_coords(4,1))/2; 
           
           % Computing the B_Membrane Matrix with a drilling dof
           B_mem = [1 0 0 0; 0 0 0 1; 0 1 1 0] * inv_J_TR * MN_Diff_Par * Tr; 
           
           % Assembling the B_bending Matrix
            B_b = sparse(3,12);
            B_b(1,2:3:end) = N_Diff(1,:);     
            B_b(3,3:3:end) = N_Diff(1,:);
            B_b(2,3:3:end) = N_Diff(2,:);
            B_b(3,2:3:end) = N_Diff(2,:);
            
            % Assembling the B_shear Matrix
            B_s = sparse(2,12);
            B_s(1,1:3:end) = N_Diff(1,:);
            B_s(1,2:3:end) = N(:);
            B_s(2,1:3:end) = N_Diff(2,:);
            B_s(2,3:3:end) = N(:);
           
           
        end
        
        
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(shellElement3d4n)
            EModul = shellElement3d4n.getPropertyValue('YOUNGS_MODULUS');
            prxy = shellElement3d4n.getPropertyValue('POISSON_RATIO');
            nr_gauss_points = shellElement3d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            thickness = shellElement3d4n.getPropertyValue('THICKNESS');
            alpha_shear = shellElement3d4n.getPropertyValue('SHEAR_CORRECTION_FACTOR');     % shear correction factor
        
            % Plane Stress Matrix
            D_mem = (EModul/(1-prxy^2))* [1  prxy  0; prxy  1  0 ; 0  0  (1-prxy)/2]; 
            
             % Calculate Shear Modulus
            GModul = EModul/(2*(1+prxy));
            % Moment-Curvature Equations
            D_b = [1    prxy    0; prxy     1   0; 0    0   (1-prxy)/2];
            % Material Bending Matrix D_b
            D_b = D_b * (EModul * thickness^3) / (12*(1-prxy^2));
            % Material Shear Matrix D_s
            D_s = eye(2) * alpha_shear * GModul * thickness; 
            
            [w,g] = returnGaussPoint(nr_gauss_points);
            
            stiffnessMatrixMemb = sparse (12,12);
            stiffnessMatrixBend = sparse (12,12);
            for xi = 1 : nr_gauss_points
                for eta = 1 : nr_gauss_points
                [~,  J, B_mem, ~, ~] = computeShapeFunction(shellElement3d4n,g(xi),g(eta));

                
                    stiffnessMatrixMemb = stiffnessMatrixMemb +             ...
                                                thickness* B_mem' * D_mem * B_mem * det(J) * w(xi) * w(eta);
                           
                end
            end
            
            for xi = 1 : nr_gauss_points
                    for eta = 1 : nr_gauss_points
                         [~,  J, ~, B_b, B_s] = computeShapeFunction(shellElement3d4n,g(xi),g(eta));

                        stiffnessMatrixBend = stiffnessMatrixBend +             ...
                            B_b' * D_b * B_b * det(J) * w(xi) * w(eta) + ...
                            B_s' * D_s * B_s * det(J) * w(xi) * w(eta);                       
                    end
            end
            
            stiffnessMatrix = sparse([stiffnessMatrixMemb, zeros(12,12); zeros(12,12), stiffnessMatrixBend]);
            
            
        end
        
        function dofs = getDofList(element)
            dofs([1  7  13 19]) = element.nodeArray.getDof('DISPLACEMENT_X');
            dofs([2  8  14 20]) = element.nodeArray.getDof('DISPLACEMENT_Y');
            dofs([4 10 16 22]) = element.nodeArray.getDof('DISPLACEMENT_Z');
            dofs([5 11 17 23]) = element.nodeArray.getDof('ROTATION_X');
            dofs([6 12 18 24]) = element.nodeArray.getDof('ROTATION_Y');
            dofs([3  9  15 21]) = element.nodeArray.getDof('ROTATION_Z');
%             dofs([1 4  7  10]) = element.nodeArray.getDof('DISPLACEMENT_X');
%             dofs([2 5  8  11]) = element.nodeArray.getDof('DISPLACEMENT_Y');
%             dofs([3 6  9  12]) = element.nodeArray.getDof('ROTATION_Z');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,24);
            vals([1  7  13 19]) = element.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals([2  8  14 20]) = element.nodeArray.getDofValue('DISPLACEMENT_Y',step);
            vals([4 10 16 22]) = element.nodeArray.getDofValue('DISPLACEMENT_Z',step);
            vals([5 11 17 23]) = element.nodeArray.getDofValue('ROTATION_X',step);
            vals([6 12 18 24]) = element.nodeArray.getDofValue('ROTATION_Y',step);
            vals([3  9  15 21]) = element.nodeArray.getDofValue('ROTATION_Z',step);
%             vals = zeros(1,12);
%             vals([1 4  7  10]) = element.nodeArray.getDofValue('DISPLACEMENT_X',step);
%             vals([2 5  8  11]) = element.nodeArray.getDofValue('DISPLACEMENT_Y',step);
%             vals([3 6  9  12]) = element.nodeArray.getDofValue('ROTATION_Z',step);
        end
    end
end




