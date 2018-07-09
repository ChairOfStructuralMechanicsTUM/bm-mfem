classdef DiscreteKirchhoffElement3d4n < QuadrilateralElement 
    %DISCRETEKIRCHHOFFELEMENT3D4N  A quadrilateral plate element
    %   Detailed explanation goes here
    
    properties (Access = private)
    end 
    
    methods
        % Constructor
        
        function discreteKirchhoffElement3d4n = DiscreteKirchhoffElement3d4n(id,nodeArray)
            
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
            discreteKirchhoffElement3d4n@QuadrilateralElement(super_args{:});
            discreteKirchhoffElement3d4n.dofNames = cellstr(["DISPLACEMENT_Z", ...
                                                "ROTATION_X", "ROTATION_Y"]);
        end
        
        %Initialization
        function initialize(discreteKirchhoffElement3d4n)
            discreteKirchhoffElement3d4n.lengthX = computeLength(discreteKirchhoffElement3d4n.nodeArray(1).getCoords, ...
                discreteKirchhoffElement3d4n.nodeArray(2).getCoords);
            
            discreteKirchhoffElement3d4n.lengthY = computeLength(discreteKirchhoffElement3d4n.nodeArray(1).getCoords, ...
                discreteKirchhoffElement3d4n.nodeArray(4).getCoords);
            
            checkConvexity(discreteKirchhoffElement3d4n);
        end
        
        function responseDoF = getResponseDofArray(discreteKirchhoffElement3d4n, step)
           
            responseDoF = zeros(12,1);
            for itNodes = 1:1:4
                nodalDof = discreteKirchhoffElement3d4n.nodeArray(itNodes).getDofArray;
                nodalDof = nodalDof.';
                
                for itDof = 3:(-1):1
                    responseDoF(3*itNodes-(itDof-1),1) = nodalDof(4-itDof).getValue(step);
                end
            end
        end

        function [B_b, J] = computeShapeFunction(discreteKirchhoffElement3d4n,xi,eta)
            N_Diff_Par = [-(1-eta)/4    (1-eta)/4   (1+eta)/4   -(1+eta)/4
                          -(1-xi)/4     -(1+xi)/4   (1+xi)/4    (1-xi)/4];

            % Coordinates of the nodes forming one element 
            ele_coords = zeros(4,2); 
            for i=1:4
                ele_coords(i,1) = discreteKirchhoffElement3d4n.nodeArray(i).getX;
                ele_coords(i,2) = discreteKirchhoffElement3d4n.nodeArray(i).getY;
            end
            
            % Jacobian and inverse Jacobian
            J = N_Diff_Par * ele_coords;
            inv_J = inv(J);
            inv_J_TR = [inv_J, zeros(2,2); zeros(2,2), inv_J]; 
            
            
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
                                    
            A = sparse(3,4); 
            A(1,1) = 1; 
            A(2,4) = 1; 
            A(3,2) = 1; 
            A(3,3) = 1;
            

            
            % Computing the B_Bending Matrix
            B_b = A * inv_J_TR * Psi_Diff_Par; 
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(discreteKirchhoffElement3d4n)            
            EModul = discreteKirchhoffElement3d4n.getPropertyValue('YOUNGS_MODULUS');
            prxy = discreteKirchhoffElement3d4n.getPropertyValue('POISSON_RATIO');
            nr_gauss_points = discreteKirchhoffElement3d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            thickness = discreteKirchhoffElement3d4n.getPropertyValue('THICKNESS');
            
            % Material Bending Matrix D_b
            D_b = [1    prxy    0; prxy     1   0; 0    0   (1-prxy)/2] * (EModul * thickness^3) / (12*(1-prxy^2));
           
            [w,g] = returnGaussPoint(nr_gauss_points);
            stiffnessMatrix = sparse(12,12);
            
            for xi = 1 : nr_gauss_points
                for eta = 1 : nr_gauss_points
                        [B_b, J] = computeShapeFunction(discreteKirchhoffElement3d4n,g(xi),g(eta));

                        stiffnessMatrix = stiffnessMatrix +             ...
                            B_b' * D_b * B_b * det(J) * w(xi) * w(eta);                       
                end
            end
        end
                
        function pl = drawDeformed(discreteKirchhoffElement3d4n, step, scaling)
            x = [discreteKirchhoffElement3d4n.nodeArray(1).getX, discreteKirchhoffElement3d4n.nodeArray(2).getX, ...
                 discreteKirchhoffElement3d4n.nodeArray(3).getX, discreteKirchhoffElement3d4n.nodeArray(4).getX, ...
                 discreteKirchhoffElement3d4n.nodeArray(1).getX];
             
            y = [discreteKirchhoffElement3d4n.nodeArray(1).getY, discreteKirchhoffElement3d4n.nodeArray(2).getY, ... 
                 discreteKirchhoffElement3d4n.nodeArray(3).getY, discreteKirchhoffElement3d4n.nodeArray(4).getY, ...
                 discreteKirchhoffElement3d4n.nodeArray(1).getY];
             
            z = [discreteKirchhoffElement3d4n.nodeArray(1).getZ + scaling * discreteKirchhoffElement3d4n.nodeArray(1).getDofValue('DISPLACEMENT_Z', step), ... 
                 discreteKirchhoffElement3d4n.nodeArray(2).getZ + scaling * discreteKirchhoffElement3d4n.nodeArray(2).getDofValue('DISPLACEMENT_Z', step), ... 
                 discreteKirchhoffElement3d4n.nodeArray(3).getZ + scaling * discreteKirchhoffElement3d4n.nodeArray(3).getDofValue('DISPLACEMENT_Z', step), ...
                 discreteKirchhoffElement3d4n.nodeArray(4).getZ + scaling * discreteKirchhoffElement3d4n.nodeArray(4).getDofValue('DISPLACEMENT_Z', step), ...
                 discreteKirchhoffElement3d4n.nodeArray(1).getZ + scaling * discreteKirchhoffElement3d4n.nodeArray(1).getDofValue('DISPLACEMENT_Z', step)];
            
            pl = line(x,y,z);
        end
        
        function dofs = getDofList(element)
            dofs([1 4 7 10]) = element.nodeArray.getDof('DISPLACEMENT_Z');
            dofs([2 5 8 11]) = element.nodeArray.getDof('ROTATION_X');
            dofs([3 6 9 12]) = element.nodeArray.getDof('ROTATION_Y');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,12);
            
            vals([1 4 7 10]) = element.nodeArray.getDofValue('DISPLACEMENT_Z',step);
            vals([2 5 8 11]) = element.nodeArray.getDofValue('ROTATION_X',step);
            vals([3 6 9 12]) = element.nodeArray.getDofValue('ROTATION_Y',step);
        end
        
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,12);
            
            [~, vals([1 4 7 10]), ~] = element.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
            [~, vals([2 5 8 11]), ~] = element.nodeArray.getDof('ROTATION_X').getAllValues(step);
            [~, vals([3 6 9 12]), ~] = element.nodeArray.getDof('ROTATION_Y').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,12);            
            
            [~, ~, vals([1 4 7 10])] = element.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
            [~, ~, vals([2 5 8 11])] = element.nodeArray.getDof('ROTATION_X').getAllValues(step);
            [~, ~, vals([3 6 9 12])] = element.nodeArray.getDof('ROTATION_Y').getAllValues(step);
        end
    end
end
