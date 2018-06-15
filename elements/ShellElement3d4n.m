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
%                 "DISPLACEMENT_X", "DISPLACEMENT_Y",  "ROTATION_Z"]);
%                 "DISPLACEMENT_Z", "ROTATION_X",  "ROTATION_Y"]);
        
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
        
        
        
        function [J, B_mem, B_b] = computeShapeFunction(obj, xi, eta)
            
%             % Shape Function and Derivatives                    
%             N = [(1-xi)*(1-eta)/4    (1+xi)*(1-eta)/4    (1+xi)*(1+eta)/4    (1-xi)*(1+eta)/4];  

            N_Diff_Par = [-(1-eta)/4    (1-eta)/4   (1+eta)/4   -(1+eta)/4
                                  -(1-xi)/4     -(1+xi)/4   (1+xi)/4    (1-xi)/4];

            % Coordinates of the nodes forming one element 
            ele_coords = zeros(4,2); 
            for i=1:4
                ele_coords(i,1) = obj.nodeArray(i).getX;
                ele_coords(i,2) = obj.nodeArray(i).getY;
            end
            
            % Jacobian and inverse Jacobian
            J = N_Diff_Par * ele_coords;
            inv_J = inv(J);
            inv_J_TR = [inv_J, zeros(2,2); zeros(2,2), inv_J]; 
                   
            MN_Diff_Par= [ - eta^3/8 + (3*eta)/8 - 1/4,                                 0, - eta^3/8 + eta^2/8 + eta/8 - 1/8,                                       0, eta^3/8 - (3*eta)/8 + 1/4,                                0, - eta^3/8 + eta^2/8 + eta/8 - 1/8,                                        0, - eta^3/8 + (3*eta)/8 + 1/4,                                 0, - eta^3/8 - eta^2/8 + eta/8 + 1/8,                                       0, eta^3/8 - (3*eta)/8 - 1/4,                                0, - eta^3/8 - eta^2/8 + eta/8 + 1/8,                                        0;
                                    -(xi/2 - 1/2)*((3*eta^2)/4 - 3/4),                         0, (xi/2 - 1/2)*(eta/2 - (3*eta^2)/4 + 1/4),                              0, (xi/2 + 1/2)*((3*eta^2)/4 - 3/4),                       0, (xi/2 + 1/2)*(eta/2 - (3*eta^2)/4 + 1/4),                              0, -(xi/2 + 1/2)*((3*eta^2)/4 - 3/4),                         0, -(xi/2 + 1/2)*((3*eta^2)/4 + eta/2 - 1/4),                            0, (xi/2 - 1/2)*((3*eta^2)/4 - 3/4),                       0, -(xi/2 - 1/2)*((3*eta^2)/4 + eta/2 - 1/4),                            0;
                                    0, -(eta/2 - 1/2)*((3*xi^2)/4 - 3/4),                                 0, (eta/2 - 1/2)*(xi/2 - (3*xi^2)/4 + 1/4),                         0, (eta/2 - 1/2)*((3*xi^2)/4 - 3/4),                                 0, -(eta/2 - 1/2)*((3*xi^2)/4 + xi/2 - 1/4),                           0, -(eta/2 + 1/2)*((3*xi^2)/4 - 3/4),                                 0, (eta/2 + 1/2)*((3*xi^2)/4 + xi/2 - 1/4),                         0, (eta/2 + 1/2)*((3*xi^2)/4 - 3/4),                                 0, -(eta/2 + 1/2)*(xi/2 - (3*xi^2)/4 + 1/4);
                                    0, - xi^3/8 + (3*xi)/8 - 1/4,                                        0, - xi^3/8 + xi^2/8 + xi/8 - 1/8,                                0, xi^3/8 - (3*xi)/8 - 1/4,                                        0, - xi^3/8 - xi^2/8 + xi/8 + 1/8,                                 0, - xi^3/8 + (3*xi)/8 + 1/4,                                         0, xi^3/8 + xi^2/8 - xi/8 - 1/8,                                0, xi^3/8 - (3*xi)/8 + 1/4,                                         0, xi^3/8 - xi^2/8 - xi/8 + 1/8]; 
                                      
            MN_Diff_2 = zeros(4,16);
            MN_Diff_2(1,:) = [obj.dM_dx(1)*obj.Nx(1,eta) 0 ...
                obj.dM_dx(1) * -1*obj.Nx(2,eta) 0 ...
                obj.dM_dx(2) * obj.Nx(1,eta) 0 ...
                obj.dM_dx(2) * -1*obj.Nx(2,eta) 0 ...
                obj.dM_dx(2) * obj.Nx(3,eta) 0 ...
                obj.dM_dx(2) * -1*obj.Nx(4,eta) 0 ...
                obj.dM_dx(1) * obj.Nx(3,eta) 0 ...
                obj.dM_dx(1) * -1*obj.Nx(4,eta) 0 ];
            MN_Diff_2(2,:) = [obj.Mx(1,xi)*obj.dNx_dx(1,eta) 0 ...
                obj.Mx(1,xi) * -1*obj.dNx_dx(2,eta) 0 ...
                obj.Mx(2,xi) * obj.dNx_dx(1,eta) 0 ...
                obj.Mx(2,xi) * -1*obj.dNx_dx(2,eta) 0 ...
                obj.Mx(2,xi) * obj.dNx_dx(3,eta) 0 ...
                obj.Mx(2,xi) * -1*obj.dNx_dx(4,eta) 0 ...
                obj.Mx(1,xi) * obj.dNx_dx(3,eta) 0 ...
                obj.Mx(1,xi) * -1*obj.dNx_dx(4,eta) 0 ];
            MN_Diff_2(3,:) = [ 0 obj.Mx(1,eta) * obj.dNx_dx(1,xi) ...
                0 obj.Mx(1,eta) * obj.dNx_dx(2,xi) ...
                0 obj.Mx(1,eta) * obj.dNx_dx(3,xi) ...
                0 obj.Mx(1,eta) * obj.dNx_dx(4,xi) ...
                0 obj.Mx(2,eta) * obj.dNx_dx(3,xi) ...
                0 obj.Mx(2,eta) * obj.dNx_dx(4,xi) ...
                0 obj.Mx(2,eta) * obj.dNx_dx(1,xi) ...
                0 obj.Mx(2,eta) * obj.dNx_dx(2,xi) ...
            ];
            MN_Diff_2(4,:) = [ 0 obj.dM_dx(1) * obj.Nx(1,xi) ...
                0 obj.dM_dx(1) * obj.Nx(2,xi) ...
                0 obj.dM_dx(1) * obj.Nx(3,xi) ...
                0 obj.dM_dx(1) * obj.Nx(4,xi) ...
                0 obj.dM_dx(2) * obj.Nx(3,xi) ...
                0 obj.dM_dx(2) * obj.Nx(4,xi) ...
                0 obj.dM_dx(2) * obj.Nx(1,xi) ...
                0 obj.dM_dx(2) * obj.Nx(2,xi) ...
            ];
                                
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
           
            A = sparse(3,4); 
            A(1,1) = 1; 
            A(2,4) = 1; 
            A(3,2) = 1; 
            A(3,3) = 1;
            
            % Computing the B_Membrane Matrix with a drilling dof
            B_mem = A * inv_J_TR * MN_Diff_Par * Tr; 
            B_mem = A * inv_J_TR * MN_Diff_2 * Tr; 

            % Assembling the Discrete Kirchhoff Shape Function Matrix
            Psi_Diff_Par = zeros(4,12); 

            Psi_Diff_Par(1,1) = 0.75 * ((( 2 * xi * (1 - eta) * (ele_coords(1,1) - ele_coords(2,1))) / ( (ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)) ...
                                    - (((1 - eta^2) * (ele_coords(4,1) - ele_coords(1,1))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2  )));

            Psi_Diff_Par(1,2) = 0.375 * (((- 2 * xi * (1 - eta) * (ele_coords(1,1) - ele_coords(2,1)) * (ele_coords(1,2) - ele_coords(2,2))) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)) ...
                                    - (((1 - eta^2) * (ele_coords(4,1) - ele_coords(1,1)) * (ele_coords(4,2) - ele_coords(1,2))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)));

            Psi_Diff_Par(1,3) = 0.25 * ( - (xi - 1) * (eta - 1) - (eta - 1) * (eta + xi + 1)) ... 
                                        + ((xi * (1- eta) * (0.25 * (ele_coords(1,1) - ele_coords(2,1))^2 - 0.5 * (ele_coords(1,2) - ele_coords(2,2))^2  )) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)) ...
                                        + 0.5 * ( ((1- eta^2) * (0.25 * (ele_coords(1,1) - ele_coords(4,1))^2 - 0.5 * (ele_coords(1,2) - ele_coords(4,2))^2  )) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2));

            Psi_Diff_Par(1,4) = 0.75 * ( - (((1 - eta^2) * (ele_coords(2,1) - ele_coords(3,1))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2  )) ...
                                    - ((2 * xi * (1 - eta) * (ele_coords(1,1) - ele_coords(2,1))) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)));

            Psi_Diff_Par(1,5) = 0.375 * (((- 2 * xi * (1 - eta) * (ele_coords(1,1) - ele_coords(2,1)) * (ele_coords(1,2) - ele_coords(2,2))) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)) ...
                                    - (((1 - eta^2) * (ele_coords(2,1) - ele_coords(3,1)) * (ele_coords(2,2) - ele_coords(3,2))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)));

            Psi_Diff_Par(1,6) = 0.25 * ( - (xi + 1) * (eta - 1) + (eta - 1) * (eta - xi + 1)) ... 
                                        + ((xi * (1- eta) * (0.25 * (ele_coords(1,1) - ele_coords(2,1))^2 - 0.5 * (ele_coords(1,2) - ele_coords(2,2))^2  )) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)) ...
                                            - 0.5 * ( ((1- eta^2) * (0.25 * (ele_coords(2,1) - ele_coords(3,1))^2 - 0.5 * (ele_coords(2,2) - ele_coords(3,2))^2  )) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)); 

            Psi_Diff_Par(1,7) = 0.75 * ((((1 - eta^2) * (ele_coords(2,1) - ele_coords(3,1))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2  )) ...
                                         + ((2 * xi * (1 + eta) * (ele_coords(3,1) - ele_coords(4,1))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)));

            Psi_Diff_Par(1,8) = 0.375 * (((- 2 * xi * (1 + eta) * (ele_coords(3,1) - ele_coords(4,1)) * (ele_coords(3,2) - ele_coords(4,2))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)) ...
                                    + (((1 - eta^2) * (ele_coords(2,1) - ele_coords(3,1)) * (ele_coords(2,2) - ele_coords(3,2))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)));

            Psi_Diff_Par(1,9) = 0.25 * ((xi + 1) * (eta + 1) + (eta + 1) * (eta + xi - 1)) ... 
                                        + ((xi * (1 + eta) * (0.25 * (ele_coords(3,1) - ele_coords(4,1))^2 - 0.5 * (ele_coords(3,2) - ele_coords(4,2))^2  )) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)) ...
                                            - 0.5 * ( ((1- eta^2) * (0.25 * (ele_coords(2,1) - ele_coords(3,1))^2 - 0.5 * (ele_coords(2,2) - ele_coords(3,2))^2  )) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)); 

            Psi_Diff_Par(1,10) = 0.75 * ((((1 - eta^2) * (ele_coords(4,1) - ele_coords(1,1))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2  )) ...
                                         - ((2 * xi * (1 + eta) * (ele_coords(3,1) - ele_coords(4,1))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)));

            Psi_Diff_Par(1,11) = 0.375 * (((- 2 * xi * (1 + eta) * (ele_coords(3,1) - ele_coords(4,1)) * (ele_coords(3,2) - ele_coords(4,2))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)) ...
                                        - (((1 - eta^2) * (ele_coords(4,1) - ele_coords(1,1)) * (ele_coords(4,2) - ele_coords(1,2))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)));

            Psi_Diff_Par(1,12) = 0.25 * ((xi - 1) * (eta + 1) + (eta + 1) * (eta - xi + 1)) ... 
                                        + ((xi * (1 + eta) * (0.25 * (ele_coords(3,1) - ele_coords(4,1))^2 - 0.5 * (ele_coords(3,2) - ele_coords(4,2))^2  )) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)) ...
                                            + 0.5 * ( ((1- eta^2) * (0.25 * (ele_coords(1,1) - ele_coords(4,1))^2 - 0.5 * (ele_coords(1,2) - ele_coords(4,2))^2  )) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)); 

            Psi_Diff_Par(2,1) = 0.75 * ((((1 - xi^2) * (ele_coords(1,1) - ele_coords(2,1))) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2  )) ...
                                         - ((2 * eta * (1 - xi) * (ele_coords(4,1) - ele_coords(1,1))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)));

            Psi_Diff_Par(2,2) = 0.375 * (((- 2 * eta * (1 - xi) * (ele_coords(4,1) - ele_coords(1,1)) * (ele_coords(4,2) - ele_coords(1,2))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)) ...
                                        - (((1 - xi^2) * (ele_coords(1,1) - ele_coords(2,1)) * (ele_coords(1,2) - ele_coords(2,2))) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)));

            Psi_Diff_Par(2,3) =   0.25 * (- (xi - 1) * (eta - 1) - (xi - 1) * (eta + xi + 1)) ... 
                                        + ((eta * (1 - xi) * (0.25 * (ele_coords(1,1) - ele_coords(4,1))^2 - 0.5 * (ele_coords(1,2) - ele_coords(4,2))^2  )) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)) ...
                                        + 0.5 * (((1 - xi^2) * (0.25 * (ele_coords(1,1) - ele_coords(2,1))^2 - 0.5 * (ele_coords(1,2) - ele_coords(2,2))^2  )) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)); 

            Psi_Diff_Par(2,4) =   0.75 * ((( - (1 - xi^2) * (ele_coords(1,1) - ele_coords(2,1))) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2  )) ...
                                         + ((2 * eta * (1 + xi) * (ele_coords(2,1) - ele_coords(3,1))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)));

            Psi_Diff_Par(2,5) = 0.375 * (((- 2 * eta * (1 + xi) * (ele_coords(2,1) - ele_coords(3,1)) * (ele_coords(2,2) - ele_coords(3,2))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)) ...
                                        - (((1 - xi^2) * (ele_coords(1,1) - ele_coords(2,1)) * (ele_coords(1,2) - ele_coords(2,2))) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)));

            Psi_Diff_Par(2,6) =   0.25 * ((xi + 1) * (eta - 1) + (xi + 1) * (eta - xi + 1)) ... 
                                        + ((eta * (1 + xi) * (0.25 * (ele_coords(2,1) - ele_coords(3,1))^2 - 0.5 * (ele_coords(2,2) - ele_coords(3,2))^2  )) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)) ...
                                        + 0.5 * (((1 - xi^2) * (0.25 * (ele_coords(1,1) - ele_coords(2,1))^2 - 0.5 * (ele_coords(1,2) - ele_coords(2,2))^2  )) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)); 

            Psi_Diff_Par(2,7) =   0.75 * ((( - (1 - xi^2) * (ele_coords(3,1) - ele_coords(4,1))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2  )) ...
                                         - ((2 * eta * (1 + xi) * (ele_coords(2,1) - ele_coords(3,1))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)));

            Psi_Diff_Par(2,8) = 0.375 * (((- 2 * eta * (1 + xi) * (ele_coords(2,1) - ele_coords(3,1)) * (ele_coords(2,2) - ele_coords(3,2))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)) ...
                                        + (((1 - xi^2) * (ele_coords(3,1) - ele_coords(4,1)) * (ele_coords(3,2) - ele_coords(4,2))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)));

            Psi_Diff_Par(2,9) =   0.25 * ((xi + 1) * (eta + 1) + (xi + 1) * (eta + xi - 1)) ... 
                                        + ((eta * (1 + xi) * (0.25 * (ele_coords(2,1) - ele_coords(3,1))^2 - 0.5 * (ele_coords(2,2) - ele_coords(3,2))^2  )) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)) ...
                                        - 0.5 * (((1 - xi^2) * (0.25 * (ele_coords(3,1) - ele_coords(4,1))^2 - 0.5 * (ele_coords(3,2) - ele_coords(4,2))^2  )) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)); 

            Psi_Diff_Par(2,10) = 0.75 * ((((1 - xi^2) * (ele_coords(3,1) - ele_coords(4,1))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2  )) ...
                                         + ((2 * eta * (1 - xi) * (ele_coords(4,1) - ele_coords(1,1))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)));

            Psi_Diff_Par(2,11) = 0.375 * (((- 2 * eta * (1 - xi) * (ele_coords(4,1) - ele_coords(1,1)) * (ele_coords(4,2) - ele_coords(1,2))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)) ...
                                        + (((1 - xi^2) * (ele_coords(3,1) - ele_coords(4,1)) * (ele_coords(3,2) - ele_coords(4,2))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)));

            Psi_Diff_Par(2,12) =   0.25 * ( - (xi - 1) * (eta + 1) + (xi - 1) * (-eta + xi + 1)) ... 
                                        + ((eta * (1 - xi) * (0.25 * (ele_coords(1,1) - ele_coords(4,1))^2 - 0.5 * (ele_coords(1,2) - ele_coords(4,2))^2  )) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)) ...
                                        - 0.5 * (((1 - xi^2) * (0.25 * (ele_coords(3,1) - ele_coords(4,1))^2 - 0.5 * (ele_coords(3,2) - ele_coords(4,2))^2  )) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)); 

            Psi_Diff_Par(3,1) = 0.75 * ((( 2 * xi * (1 - eta) * (ele_coords(1,2) - ele_coords(2,2))) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)) ...
                                    - (((1 - eta^2) * (ele_coords(4,2) - ele_coords(1,2))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2  )));

            Psi_Diff_Par(3,2) = 0.25 * ((xi - 1) * (eta - 1) + (eta - 1) * (eta + xi + 1)) ... 
                                        - ((xi * (1- eta) * (- 0.5 * (ele_coords(1,1) - ele_coords(2,1))^2 + 0.25 * (ele_coords(1,2) - ele_coords(2,2))^2  )) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)) ...
                                        - 0.5 * (((1- eta^2) * (- 0.5 * (ele_coords(1,1) - ele_coords(4,1))^2 + 0.25 * (ele_coords(1,2) - ele_coords(4,2))^2  )) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)); 

            Psi_Diff_Par(3,3) = 0.375 * ((( 2 * xi * (1 - eta) * (ele_coords(1,1) - ele_coords(2,1)) * (ele_coords(1,2) - ele_coords(2,2))) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)) ...
                                    - (((1 - eta^2) * (ele_coords(4,1) - ele_coords(1,1)) * (ele_coords(4,2) - ele_coords(1,2))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)));

            Psi_Diff_Par(3,4) = 0.75 * ((( - 2 * xi * (1 - eta) * (ele_coords(1,2) - ele_coords(2,2))) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)) ...
                                    - (((1 - eta^2) * (ele_coords(2,2) - ele_coords(3,2))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2  )));

            Psi_Diff_Par(3,5) = 0.25 * ((xi + 1) * (eta - 1) - (eta - 1) * (eta - xi + 1)) ... 
                                        - ((xi * (1- eta) * (- 0.5 * (ele_coords(1,1) - ele_coords(2,1))^2 + 0.25 * (ele_coords(1,2) - ele_coords(2,2))^2  )) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)) ...
                                        + 0.5 * (((1 - eta^2) * (- 0.5 * (ele_coords(2,1) - ele_coords(3,1))^2 + 0.25 * (ele_coords(2,2) - ele_coords(3,2))^2  )) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)); 

            Psi_Diff_Par(3,6) = 0.375 * ((( 2 * xi * (1 - eta) * (ele_coords(1,1) - ele_coords(2,1)) * (ele_coords(1,2) - ele_coords(2,2))) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)) ...
                                    - (((1 - eta^2) * (ele_coords(2,1) - ele_coords(3,1)) * (ele_coords(2,2) - ele_coords(3,2))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)));

            Psi_Diff_Par(3,7) = 0.75 * (((2 * xi * (1 + eta) * (ele_coords(3,2) - ele_coords(4,2))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)) ...
                                    + (((1 - eta^2) * (ele_coords(2,2) - ele_coords(3,2))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2  )));

            Psi_Diff_Par(3,8) = 0.25 * ( - (xi + 1) * (eta + 1) - (eta + 1) * (eta + xi - 1)) ... 
                                        - ((xi * (1+ eta) * (- 0.5 * (ele_coords(3,1) - ele_coords(4,1))^2 + 0.25 * (ele_coords(3,2) - ele_coords(4,2))^2  )) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)) ...
                                        + 0.5 * (((1 - eta^2) * (- 0.5 * (ele_coords(2,1) - ele_coords(3,1))^2 + 0.25 * (ele_coords(2,2) - ele_coords(3,2))^2  )) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)); 

            Psi_Diff_Par(3,9) = 0.375 * ((( 2 * xi * (1 + eta) * (ele_coords(3,1) - ele_coords(4,1)) * (ele_coords(3,2) - ele_coords(4,2))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)) ...
                                    - (((1 - eta^2) * (ele_coords(2,1) - ele_coords(3,1)) * (ele_coords(2,2) - ele_coords(3,2))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)));                    

            Psi_Diff_Par(3,10) = 0.75 * (((- 2 * xi * (1 + eta) * (ele_coords(3,2) - ele_coords(4,2))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)) ...
                                    + (((1 - eta^2) * (ele_coords(4,2) - ele_coords(1,2))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2  )));

            Psi_Diff_Par(3,11) = 0.25 * ( - (xi - 1) * (eta + 1) - (eta + 1) * ( - eta + xi + 1)) ... 
                                        - ((xi * (1+ eta) * (- 0.5 * (ele_coords(3,1) - ele_coords(4,1))^2 + 0.25 * (ele_coords(3,2) - ele_coords(4,2))^2  )) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)) ...
                                        - 0.5 * (((1 - eta^2) * ( - 0.5 * (ele_coords(1,1) - ele_coords(4,1))^2 + 0.25 * (ele_coords(1,2) - ele_coords(4,2))^2  )) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)); 

            Psi_Diff_Par(3,12) = 0.375 * ((( 2 * xi * (1 + eta) * (ele_coords(3,1) - ele_coords(4,1)) * (ele_coords(3,2) - ele_coords(4,2))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)) ...
                                    + (((1 - eta^2) * (ele_coords(4,1) - ele_coords(1,1)) * (ele_coords(4,2) - ele_coords(1,2))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)));                    

            Psi_Diff_Par(4,1) = 0.75 * ((((1 - xi^2) * (ele_coords(1,2) - ele_coords(2,2))) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2  )) ...
                                         - ((2 * eta * (1 - xi) * (ele_coords(4,2) - ele_coords(1,2))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)));

            Psi_Diff_Par(4,2) = 0.25 * ((xi - 1) * (eta - 1) + (xi - 1) * (eta + xi + 1)) ... 
                                        - ((eta * (1 - xi) * ( - 0.5 * (ele_coords(1,1) - ele_coords(4,1))^2 + 0.25 * (ele_coords(1,2) - ele_coords(4,2))^2  )) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)) ...
                                        - 0.5 * (((1 - xi^2) * ( - 0.5 * (ele_coords(1,1) - ele_coords(2,1))^2 + 0.25 * (ele_coords(1,2) - ele_coords(2,2))^2  )) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)); 
            
            Psi_Diff_Par(4,3) = 0.375 * (((2 * eta * (1 - xi) * (ele_coords(4,1) - ele_coords(1,1)) * (ele_coords(4,2) - ele_coords(1,2))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)) ...
                                        + (((1 - xi^2) * (ele_coords(1,1) - ele_coords(2,1)) * (ele_coords(1,2) - ele_coords(2,2))) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)));   

            Psi_Diff_Par(4,4) =   0.75 * ((( - (1 - xi^2) * (ele_coords(1,2) - ele_coords(2,2))) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2  )) ...
                                         + ((2 * eta * (1 + xi) * (ele_coords(2,2) - ele_coords(3,2))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)));

            Psi_Diff_Par(4,5) = 0.25 * ( - (xi + 1) * (eta - 1) - (xi + 1) * (eta - xi + 1)) ... 
                                        - ((eta * (1 + xi) * ( - 0.5 * (ele_coords(2,1) - ele_coords(3,1))^2 + 0.25 * (ele_coords(2,2) - ele_coords(3,2))^2  )) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)) ...
                                        - 0.5 * (((1 - xi^2) * (-0.5 * (ele_coords(1,1) - ele_coords(2,1))^2 + 0.25 * (ele_coords(1,2) - ele_coords(2,2))^2  )) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)); 

            Psi_Diff_Par(4,6) = 0.375 * (((2 * eta * (1 + xi) * (ele_coords(2,1) - ele_coords(3,1)) * (ele_coords(2,2) - ele_coords(3,2))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)) ...
                                        + (((1 - xi^2) * (ele_coords(1,1) - ele_coords(2,1)) * (ele_coords(1,2) - ele_coords(2,2))) / ((ele_coords(1,1) - ele_coords(2,1))^2 + (ele_coords(1,2) - ele_coords(2,2))^2)));

            Psi_Diff_Par(4,7) = 0.75 * ((( - (1 - xi^2) * (ele_coords(3,2) - ele_coords(4,2))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2  )) ...
                                         - ((2 * eta * (1 + xi) * (ele_coords(2,2) - ele_coords(3,2))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)));

            Psi_Diff_Par(4,8) = 0.25 * ( - (xi + 1) * (eta + 1) - (xi + 1) * (eta + xi - 1)) ... 
                                        - ((eta * (1 + xi) * ( - 0.5 * (ele_coords(2,1) - ele_coords(3,1))^2 + 0.25 * (ele_coords(2,2) - ele_coords(3,2))^2  )) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)) ...
                                        + 0.5 * (((1 - xi^2) * ( - 0.5 * (ele_coords(3,1) - ele_coords(4,1))^2 + 0.25 * (ele_coords(3,2) - ele_coords(4,2))^2  )) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)); 
            
            Psi_Diff_Par(4,9) = 0.375 * (((2 * eta * (1 + xi) * (ele_coords(2,1) - ele_coords(3,1)) * (ele_coords(2,2) - ele_coords(3,2))) / ((ele_coords(2,1) - ele_coords(3,1))^2 + (ele_coords(2,2) - ele_coords(3,2))^2)) ...
                                        - (((1 - xi^2) * (ele_coords(3,1) - ele_coords(4,1)) * (ele_coords(3,2) - ele_coords(4,2))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)));  

            Psi_Diff_Par(4,10) = 0.75 * ((((1 - xi^2) * (ele_coords(3,2) - ele_coords(4,2))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2  )) ...
                                         + ((2 * eta * (1 - xi) * (ele_coords(4,2) - ele_coords(1,2))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)));

            Psi_Diff_Par(4,11) =   0.25 * ((xi - 1) * (eta + 1) - (xi - 1) * (-eta + xi + 1)) ... 
                                        - ((eta * (1 - xi) * ( - 0.5 * (ele_coords(1,1) - ele_coords(4,1))^2 + 0.25 * (ele_coords(1,2) - ele_coords(4,2))^2  )) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)) ...
                                        + 0.5 * (((1 - xi^2) * ( - 0.5 * (ele_coords(3,1) - ele_coords(4,1))^2 + 0.25 * (ele_coords(3,2) - ele_coords(4,2))^2  )) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)); 

            Psi_Diff_Par(4,12) = 0.375 * (((2 * eta * (1 - xi) * (ele_coords(4,1) - ele_coords(1,1)) * (ele_coords(4,2) - ele_coords(1,2))) / ((ele_coords(1,1) - ele_coords(4,1))^2 + (ele_coords(1,2) - ele_coords(4,2))^2)) ...
                                        - (((1 - xi^2) * (ele_coords(3,1) - ele_coords(4,1)) * (ele_coords(3,2) - ele_coords(4,2))) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)));
            
           % Computing the B_Bending Matrix
            B_b = A * inv_J_TR * Psi_Diff_Par; 
        end
        
        
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(shellElement3d4n)
            
            EModul = shellElement3d4n.getPropertyValue('YOUNGS_MODULUS');
            prxy = shellElement3d4n.getPropertyValue('POISSON_RATIO');
            nr_gauss_points = shellElement3d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            thickness = shellElement3d4n.getPropertyValue('THICKNESS');
        
            % Plane Stress Matrix
            D_mem = (EModul/(1-prxy^2))* [1  prxy  0; prxy  1  0 ; 0  0  (1-prxy)/2]; 
            
            % Material Bending Matrix D_b
            D_b = [1    prxy    0; prxy     1   0; 0    0   (1-prxy)/2] * (EModul * thickness^3) / (12*(1-prxy^2));
           
            [w,g] = returnGaussPoint(nr_gauss_points);
            
            stiffnessMatrixMemb = zeros (12,12);
            stiffnessMatrixBend = zeros (12,12);
            stiffnessMatrix = zeros(24,24);

            for xi = 1 : nr_gauss_points
                for eta = 1 : nr_gauss_points
                [J, B_mem, B_b] = computeShapeFunction(shellElement3d4n,g(xi),g(eta));

                
                    stiffnessMatrixMemb = stiffnessMatrixMemb +             ...
                                                thickness.* B_mem' * D_mem * B_mem * det(J) * w(xi) * w(eta);
                           
                    stiffnessMatrixBend = stiffnessMatrixBend +             ...
                                                B_b' * D_b * B_b * det(J) * w(xi) * w(eta);       
                
                end
            end
%             stiffnessMatrix = sparse([stiffnessMatrixMemb, zeros(12,12); zeros(12,12), stiffnessMatrixBend]);
            % assemble stiffness matrices
            stiffnessMatrix([1 2 6 7 8 12 13 14 18 19 20 24],[1 2 6 7 8 12 13 14 18 19 20 24]) = stiffnessMatrixMemb;
            stiffnessMatrix([3 4 5 9 10 11 15 16 17 21 22 23],[3 4 5 9 10 11 15 16 17 21 22 23]) = stiffnessMatrixBend;
%             stiffnessMatrix2(:,1) = stiffnessMatrixMemb(:,1);

%             stiffnessMatrix = stiffnessMatrixMemb;
%             for xi = 1 : nr_gauss_points
%                     for eta = 1 : nr_gauss_points
%                          [J, ~, B_b] = computeShapeFunction(shellElement3d4n,g(xi),g(eta));
% 
%                          stiffnessMatrixBend = stiffnessMatrixBend +             ...
%                                                     B_b' * D_b * B_b * det(J) * w(xi) * w(eta);               
%                     end
%             end
%               
%             stiffnessMatrix = stiffnessMatrixBend;
            
        end
        
        
        function massMatrix = computeLocalMassMatrix(shellElement3d4n)
            density = shellElement3d4n.getPropertyValue('DENSITY');
            thickness = shellElement3d4n.getPropertyValue('THICKNESS');
            nr_gauss_points = shellElement3d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            
            
            massMatrix = zeros(24,24); 
            
            % Coordinates of the nodes forming one element 
            ele_coords = zeros(4,2); 
            for i=1:4
                ele_coords(i,1) = shellElement3d4n.nodeArray(i).getX;
                ele_coords(i,2) = shellElement3d4n.nodeArray(i).getY;
            end
            
            
            lumpArea = computeLength(ele_coords(1,:), ele_coords(4,:)) * computeLength(ele_coords(1,:), ele_coords(2,:));
            nodalMass = 0.25 * density * thickness * lumpArea; 
            index = 1; 
            for i = 1 : 3
                
               massMatrix(index,index) = nodalMass; 
               massMatrix(index+1,index+1) = nodalMass; 
               massMatrix(index+2,index+2) = nodalMass; 
               index = i * 6;
            end
           
        end
        
        function dampingMatrix=computeLocalDampingMatrix(element)
        dampingMatrix = zeros(24,24); 
        end
        
        function dofs = getDofList(element)
            dofs([1  7  13 19]) = element.nodeArray.getDof('DISPLACEMENT_X');
            dofs([2  8  14 20]) = element.nodeArray.getDof('DISPLACEMENT_Y');
            dofs([3  9  15 21]) = element.nodeArray.getDof('DISPLACEMENT_Z');
            dofs([4 10 16 22]) = element.nodeArray.getDof('ROTATION_X');
            dofs([5 11 17 23]) = element.nodeArray.getDof('ROTATION_Y');
            dofs([6 12 18 24]) = element.nodeArray.getDof('ROTATION_Z');
            

        end
%         
%         function vals = getValuesVector(element, step)
%             vals = zeros(1,24);
%             vals([1  7  13 19]) = element.nodeArray.getDofValue('DISPLACEMENT_X',step);
%             vals([2  8  14 20]) = element.nodeArray.getDofValue('DISPLACEMENT_Y',step);
%             vals([4 10 16 22]) = element.nodeArray.getDofValue('DISPLACEMENT_Z',step);
%             vals([5 11 17 23]) = element.nodeArray.getDofValue('ROTATION_X',step);
%             vals([6 12 18 24]) = element.nodeArray.getDofValue('ROTATION_Y',step);
%             vals([3  9  15 21]) = element.nodeArray.getDofValue('ROTATION_Z',step);
%         end
        
%          function dofs = getDofList(element)
%              dofs([1 7 13 19]) = element.nodeArray.getDof('DISPLACEMENT_X');
%              dofs([2 8 14 20]) = element.nodeArray.getDof('DISPLACEMENT_Y');
%              dofs([3 9 15 21]) = element.nodeArray.getDof('ROTATION_Z');
%              dofs([4 10 16 22]) = element.nodeArray.getDof('DISPLACEMENT_Z');
%              dofs([5 11 17 23]) = element.nodeArray.getDof('ROTATION_Y');
%              dofs([6 12 18 24]) = element.nodeArray.getDof('ROTATION_X');
%          end
        
%          function dofs = getDofList(element)
%             dofs([1 4 7 10]) = element.nodeArray.getDof('DISPLACEMENT_Z');
%             dofs([2 5 8 11]) = element.nodeArray.getDof('ROTATION_X');
%             dofs([3 6 9 12]) = element.nodeArray.getDof('ROTATION_Y');
% end

%          function dofs = getDofList(element)
%             dofs([1 4 7 10]) = element.nodeArray.getDof('DISPLACEMENT_X');
%             dofs([2 5 8 11]) = element.nodeArray.getDof('DISPLACEMENT_Y');
%             dofs([3 6 9 12]) = element.nodeArray.getDof('ROTATION_Z');
% end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,12);
            
%             vals([1 4 7 10]) = element.nodeArray.getDofValue('DISPLACEMENT_Z',step);
%             vals([2 5 8 11]) = element.nodeArray.getDofValue('ROTATION_X',step);
%             vals([3 6 9 12]) = element.nodeArray.getDofValue('ROTATION_Y',step);
        end
    end
    
    methods (Access = private, Static = true)
        
%         function dMN_dxi(r,s,xi,eta)
%             if (r == 1) && (s == 1)
%                 
%             end
%         end
        
        function res = Mx(i,x)
            switch i
                case 1
                    res = 0.5 * (1-x);
                case 2
                    res = 0.5 * (1+x);
                otherwise
                    error('Mx only valid with 1 or 2, got %f',i)
            end
        end
        
        function res = dM_dx(i)
            switch i
                case 1
                    res = -0.5;
                case 2
                    res = 0.5;
                otherwise
                    error('dM_d only valid with 1 or 2, got %f',i)
            end
        end
        
        function res = Nx(i,x)
            switch i
                case 1
                    res = 0.5 - 0.75*x + 0.25*x^3;
                case 2
                    res = 0.25 - 0.25*x - 0.25*x^2 + 0.25*x^3;
                case 3
                    res = 0.5 + 0.75*x - 0.25*x^3;
                case 4
                    res = - 0.25 - 0.25*x + 0.25*x^2 + 0.25*x^3;
                otherwise
                    error('Nx only valid with 1,2,3 or 4, got %f',i)
            end
        end
        
        function res = dNx_dx(i,x)
            switch i
                case 1
                    res = -0.75 + 0.75 * x^2;
                case 2
                    res = -0.25 - 0.5*x + 0.75*x^2;
                case 3
                    res = 0.75 - 0.75*x^2;
                case 4
                    res = -0.25 + 0.5*x + 0.75*x^2;
                otherwise
                    error('dNx_dx only valid with 1,2,3 or 4, got %f',i)
            end
        end
        
    end
end




