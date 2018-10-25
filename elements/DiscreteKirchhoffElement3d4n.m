classdef DiscreteKirchhoffElement3d4n < QuadrilateralElement 
    %DISCRETEKIRCHHOFFELEMENT3D4N  A quadrilateral plate element based on
    %   the DKQ formulation
    %   
    %   see also SHELLELEMENT3D4N
    
    properties (Access = private)
    end
    
    methods
        eig
        % Constructor
        function obj = DiscreteKirchhoffElement3d4n(id,nodeArray)
            
            requiredPropertyNames = ["YOUNGS_MODULUS", "POISSON_RATIO", ...
                "THICKNESS", "NUMBER_GAUSS_POINT", "DENSITY", ...
                "SHEAR_CORRECTION_FACTOR"];
                                         
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
            obj@QuadrilateralElement(super_args{:});
            obj.dofNames = ["DISPLACEMENT_Z", ...
                "ROTATION_X", "ROTATION_Y"];
        end
        
        %Initialization
        function initialize(obj)
            obj.lengthX = computeLength(obj.nodeArray(1).getCoords, ...
                obj.nodeArray(2).getCoords);
            
            obj.lengthY = computeLength(obj.nodeArray(1).getCoords, ...
                obj.nodeArray(4).getCoords);
            
            checkConvexity(obj);
        end
        
        function check(obj)
            if obj.getPropertyValue('NUMBER_GAUSS_POINT') < 2
                obj.setPropertyValue('NUMBER_GAUSS_POINT', 2);
            end
            
            check@Element(obj);
        end
        
        function [N, N_mat, J] = computeShapeFunction(obj, xi, eta)
            % Shape Function and Derivatives
            N = [(1-xi)*(1-eta)/4    (1+xi)*(1-eta)/4    (1+xi)*(1+eta)/4    (1-xi)*(1+eta)/4];
            N_mat = sparse(3,12);
            N_mat(1,1:3:end) = N(:);
            N_mat(2,2:3:end) = N(:);
            N_mat(3,3:3:end) = N(:);

            N_Diff_Par = sparse([-(1-eta)/4    (1-eta)/4   (1+eta)/4   -(1+eta)/4
                                  -(1-xi)/4     -(1+xi)/4   (1+xi)/4    (1-xi)/4]);
            
            % Coordinates of the nodes forming one element
            ele_coords = zeros(4,2);
            for i=1:4
                ele_coords(i,1) = obj.nodeArray(i).getX;
                ele_coords(i,2) = obj.nodeArray(i).getY;
            end
            
            % Jacobian 
            J = N_Diff_Par * ele_coords;
        end

        function B_b = computeBMatrix(obj, xi, eta, J)

            % Coordinates of the nodes forming one element 
            ele_coords = zeros(4,2); 
            for i=1:4
                ele_coords(i,1) = obj.nodeArray(i).getX;
                ele_coords(i,2) = obj.nodeArray(i).getY;
            end
            
            % Inverse Jacobian
            inv_J = inv(J);
            inv_J_TR = sparse([1 1 2 2 3 3 4 4],[1 2 1 2 3 4 3 4],...
                            [inv_J(1,1) inv_J(1,2) inv_J(2,1) inv_J(2,2) inv_J(1,1) inv_J(1,2) inv_J(2,1) inv_J(2,2)],4,4);
            
            Psi_Diff_Par = sparse(4,12); 
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

            Psi_Diff_Par(1,12) = -0.25 * (eta - 2*eta*xi + eta^2 - 2*xi) + ...
                                        + (xi * (1 + eta) * (0.25 * (ele_coords(3,1) - ele_coords(4,1))^2 - 0.5 * (ele_coords(3,2) - ele_coords(4,2))^2 ) / ((ele_coords(3,1) - ele_coords(4,1))^2 + (ele_coords(3,2) - ele_coords(4,2))^2)) ...
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
                                    
            A = sparse([1 2 3 3],[1 4 2 3],[1 1 1 1],3,4);
            
            % Computing the B_Bending Matrix
            B_b = A * inv_J_TR * Psi_Diff_Par; 
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(obj)
            EModul = obj.getPropertyValue('YOUNGS_MODULUS');
            prxy = obj.getPropertyValue('POISSON_RATIO');
            nr_gauss_points = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            thickness = obj.getPropertyValue('THICKNESS');
            
            % Material Bending Matrix D_b
            D_b = sparse([1 1 2 2 3],[1 2 1 2 3],[1 prxy prxy 1 (1-prxy)/2],3,3);
            D_b = D_b * (EModul * thickness^3) / (12*(1-prxy^2));
            
            [w,g] = returnGaussPoint(nr_gauss_points);
            stiffnessMatrix = sparse(12,12);
            
            for xi = 1 : nr_gauss_points
                for eta = 1 : nr_gauss_points
                    [~,~,J] = obj.computeShapeFunction(g(xi), g(eta));
                    B_b = obj.computeBMatrix(g(xi), g(eta), J);
                    
                    stiffnessMatrix = stiffnessMatrix + ...
                        B_b' * D_b * B_b * det(J) * w(xi) * w(eta);
                end
            end
        end
        
        function massMatrix = computeLocalMassMatrix(obj)
            density = obj.getPropertyValue('DENSITY');
            thickness = obj.getPropertyValue('THICKNESS');
            nr_gauss_points = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            
            % Coordinates of the nodes forming one element
            ele_coords = zeros(4,2);
            for i=1:4
                ele_coords(i,1) = obj.nodeArray(i).getX;
                ele_coords(i,2) = obj.nodeArray(i).getY;
            end
            
            inertia = density*thickness;
            rot_inertia = density*thickness^3/12;
            dens_mat = sparse(1:3,1:3,[inertia rot_inertia*ones(1,2)]);
            
            [w,g] = returnGaussPoint(nr_gauss_points);
            massMatrix = sparse(12,12);
            for xi = 1 : nr_gauss_points
                for eta = 1 : nr_gauss_points
                    [~, N_mat, J] = computeShapeFunction(obj, g(xi), g(eta));
                    massMatrix = massMatrix + N_mat' * dens_mat * N_mat * det(J) * w(xi) * w(eta) ;
                end
            end

        end
                
        function pl = drawDeformed(obj, step, scaling)
            x = [obj.nodeArray(1).getX, obj.nodeArray(2).getX, ...
                 obj.nodeArray(3).getX, obj.nodeArray(4).getX, ...
                 obj.nodeArray(1).getX];
             
            y = [obj.nodeArray(1).getY, obj.nodeArray(2).getY, ... 
                 obj.nodeArray(3).getY, obj.nodeArray(4).getY, ...
                 obj.nodeArray(1).getY];
             
            z = [obj.nodeArray(1).getZ + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Z', step), ... 
                 obj.nodeArray(2).getZ + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Z', step), ... 
                 obj.nodeArray(3).getZ + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Z', step), ...
                 obj.nodeArray(4).getZ + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Z', step), ...
                 obj.nodeArray(1).getZ + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Z', step)];
            
            pl = line(x,y,z);
        end
        
        function dofs = getDofList(obj)
            dofs([1 4 7 10]) = obj.nodeArray.getDof('DISPLACEMENT_Z');
            dofs([2 5 8 11]) = obj.nodeArray.getDof('ROTATION_X');
            dofs([3 6 9 12]) = obj.nodeArray.getDof('ROTATION_Y');
        end
        
        function vals = getValuesVector(obj, step)
            vals = zeros(1,12);
            vals([1 4 7 10]) = obj.nodeArray.getDofValue('DISPLACEMENT_Z',step);
            vals([2 5 8 11]) = obj.nodeArray.getDofValue('ROTATION_X',step);
            vals([3 6 9 12]) = obj.nodeArray.getDofValue('ROTATION_Y',step);
        end
        
        function vals = getFirstDerivativesVector(obj, step)
            vals = zeros(1,12);
            [~, vals([1 4 7 10]), ~] = obj.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
            [~, vals([2 5 8 11]), ~] = obj.nodeArray.getDof('ROTATION_X').getAllValues(step);
            [~, vals([3 6 9 12]), ~] = obj.nodeArray.getDof('ROTATION_Y').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(obj, step)
            vals = zeros(1,12);
            [~, ~, vals([1 4 7 10])] = obj.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
            [~, ~, vals([2 5 8 11])] = obj.nodeArray.getDof('ROTATION_X').getAllValues(step);
            [~, ~, vals([3 6 9 12])] = obj.nodeArray.getDof('ROTATION_Y').getAllValues(step);
        end
    end
end
