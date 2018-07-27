classdef ShellElement3d4n < QuadrilateralElement 
    %SHELLELEMENT3D4N A quadrilateral shell element with 6 dofs per node
    %   currently only rectangular shells in the x-y plane are supported
    %
    %   TODO: extend to shells other than in the x-y plane
    %   TODO: validate for skewed elements
    %   TODO: move the DKQ computation to functions for better readability
    %   TODO: add consistent mass matrix (dynamic calculations are not
    %       possible until then. Perhaps add a corotational formulation)
    %   TODO: lump area computation for all kinds of quads: 
    %       A=(1/2)|[(x3-x1)(y4-y2) +(x4-x2)(y1-y3)]|.
    %
    %   Based on Barrales, F. R. (2012): Development of a nonlinear quadrilateral 
    %       layered membrane element with drilling degrees of freedom and a 
    %       nonlinear quadrilateral thin flat layered shell element for the 
    %       modeling of reinforced concrete walls 
    %   Membrane part: Plane stress with rotational dof
    %   Bending part: Discrete Kirchhoff quadrilateral element
    
    properties (Access = private)
    end 
    
    methods
        % Constructor
        function obj = ShellElement3d4n(id,nodeArray)
            
            requiredPropertyNames = ["YOUNGS_MODULUS", "POISSON_RATIO", ...
                "THICKNESS", "NUMBER_GAUSS_POINT", ...
                "DENSITY", "SHEAR_CORRECTION_FACTOR"];
                                         
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
            obj.dofNames = ["DISPLACEMENT_X", "DISPLACEMENT_Y", "DISPLACEMENT_Z", ...
                "ROTATION_X",  "ROTATION_Y",  "ROTATION_Z"];
        
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
        
        function [J, B_mem, B_b] = computeShapeFunction(obj, xi, eta)
            
            % Shape Function and Derivatives
            %N = [(1-xi)*(1-eta)/4    (1+xi)*(1-eta)/4    (1+xi)*(1+eta)/4    (1-xi)*(1+eta)/4];

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
            
            MN_Diff_Par = zeros(4,16);
            MN_Diff_Par(1,:) = [obj.dM_dx(1)*obj.Nx(1,eta) 0 ...
                obj.dM_dx(1) * -1*obj.Nx(2,eta) 0 ...
                obj.dM_dx(2) * obj.Nx(1,eta) 0 ...
                obj.dM_dx(2) * -1*obj.Nx(2,eta) 0 ...
                obj.dM_dx(2) * obj.Nx(3,eta) 0 ...
                obj.dM_dx(2) * -1*obj.Nx(4,eta) 0 ...
                obj.dM_dx(1) * obj.Nx(3,eta) 0 ...
                obj.dM_dx(1) * -1*obj.Nx(4,eta) 0 ];
            MN_Diff_Par(2,:) = [obj.Mx(1,xi)*obj.dNx_dx(1,eta) 0 ...
                obj.Mx(1,xi) * -1*obj.dNx_dx(2,eta) 0 ...
                obj.Mx(2,xi) * obj.dNx_dx(1,eta) 0 ...
                obj.Mx(2,xi) * -1*obj.dNx_dx(2,eta) 0 ...
                obj.Mx(2,xi) * obj.dNx_dx(3,eta) 0 ...
                obj.Mx(2,xi) * -1*obj.dNx_dx(4,eta) 0 ...
                obj.Mx(1,xi) * obj.dNx_dx(3,eta) 0 ...
                obj.Mx(1,xi) * -1*obj.dNx_dx(4,eta) 0 ];
            MN_Diff_Par(3,:) = [ 0 obj.Mx(1,eta) * obj.dNx_dx(1,xi) ...
                0 obj.Mx(1,eta) * obj.dNx_dx(2,xi) ...
                0 obj.Mx(1,eta) * obj.dNx_dx(3,xi) ...
                0 obj.Mx(1,eta) * obj.dNx_dx(4,xi) ...
                0 obj.Mx(2,eta) * obj.dNx_dx(3,xi) ...
                0 obj.Mx(2,eta) * obj.dNx_dx(4,xi) ...
                0 obj.Mx(2,eta) * obj.dNx_dx(1,xi) ...
                0 obj.Mx(2,eta) * obj.dNx_dx(2,xi) ];
            MN_Diff_Par(4,:) = [ 0 obj.dM_dx(1) * obj.Nx(1,xi) ...
                0 obj.dM_dx(1) * obj.Nx(2,xi) ...
                0 obj.dM_dx(1) * obj.Nx(3,xi) ...
                0 obj.dM_dx(1) * obj.Nx(4,xi) ...
                0 obj.dM_dx(2) * obj.Nx(3,xi) ...
                0 obj.dM_dx(2) * obj.Nx(4,xi) ...
                0 obj.dM_dx(2) * obj.Nx(1,xi) ...
                0 obj.dM_dx(2) * obj.Nx(2,xi) ];
                                
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
            
           % Compute the B_Bending Matrix
            B_b = A * inv_J_TR * Psi_Diff_Par; 
        end
        
        
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(obj)
            
            EModul = obj.getPropertyValue('YOUNGS_MODULUS');
            prxy = obj.getPropertyValue('POISSON_RATIO');
            nr_gauss_points = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            thickness = obj.getPropertyValue('THICKNESS');
            
            % Plane Stress Matrix
            D_mem = (EModul/(1-prxy^2)) * ...
                [1  prxy  0; prxy  1  0 ; 0  0  (1-prxy)/2];
            
            % Material Bending Matrix D_b
            D_b = [1    prxy    0; ...
                prxy     1   0; ...
                0    0   (1-prxy)/2] * (EModul * thickness^3) / (12*(1-prxy^2));
            
            [w,g] = returnGaussPoint(nr_gauss_points);
            
            stiffnessMatrixMemb = zeros(12,12);
            stiffnessMatrixBend = zeros(12,12);
            stiffnessMatrix = sparse(24,24);
            
            for xi = 1 : nr_gauss_points
                for eta = 1 : nr_gauss_points
                    [J, B_mem, B_b] = computeShapeFunction(obj,g(xi),g(eta));
                    
                    stiffnessMatrixMemb = stiffnessMatrixMemb + ...
                        thickness.* B_mem' * D_mem * B_mem * det(J) * w(xi) * w(eta);
                    
                    stiffnessMatrixBend = stiffnessMatrixBend + ...
                        B_b' * D_b * B_b * det(J) * w(xi) * w(eta);
                    
                end
            end
            
            % assemble stiffness matrix:
            % membrane part: ux, uy, phiz
            % bending part: uz, phix, phiy
            stiffnessMatrix([1 2 6 7 8 12 13 14 18 19 20 24],...
                [1 2 6 7 8 12 13 14 18 19 20 24]) = stiffnessMatrixMemb;
            stiffnessMatrix([3 4 5 9 10 11 15 16 17 21 22 23],...
                [3 4 5 9 10 11 15 16 17 21 22 23]) = stiffnessMatrixBend;
            
        end
        
        
        function massMatrix = computeLocalMassMatrix(obj)
            density = obj.getPropertyValue('DENSITY');
            thickness = obj.getPropertyValue('THICKNESS');
            %nr_gauss_points = obj.getPropertyValue('NUMBER_GAUSS_POINT');
            
            % Coordinates of the nodes forming one element 
            ele_coords = zeros(4,2); 
            for i=1:4
                ele_coords(i,1) = obj.nodeArray(i).getX;
                ele_coords(i,2) = obj.nodeArray(i).getY;
            end
            
            lumpArea = computeLength(ele_coords(1,:), ele_coords(4,:)) * ...
                computeLength(ele_coords(1,:), ele_coords(2,:));
            nodalMass = 0.25 * density * thickness * lumpArea; 
            i = [1 2 3 7 8 9 13 14 15 19 20 21];
            
            massMatrix = sparse(i,i,nodalMass*ones(1,12),24,24);

        end
        
        function dampingMatrix = computeLocalDampingMatrix(obj)
            props = obj.getProperties;
            
            if props.hasValue('RAYLEIGH_ALPHA') && props.hasValue('RAYLEIGH_BETA')
                alpha = props.getValue('RAYLEIGH_ALPHA');
                beta = props.getValue('RAYLEIGH_BETA');
                dampingMatrix = alpha * obj.computeLocalMassMatrix + ...
                    beta * obj.computeLocalStiffnessMatrix;
            else
                dampingMatrix = sparse(24,24);
            end
            
        end
        
        function f=computeLocalForceVector(obj)
            f = zeros(1,24); 
        end
        
        function dofs = getDofList(obj)
            dofs([1  7  13 19]) = obj.nodeArray.getDof('DISPLACEMENT_X');
            dofs([2  8  14 20]) = obj.nodeArray.getDof('DISPLACEMENT_Y');
            dofs([3  9  15 21]) = obj.nodeArray.getDof('DISPLACEMENT_Z');
            dofs([4 10 16 22]) = obj.nodeArray.getDof('ROTATION_X');
            dofs([5 11 17 23]) = obj.nodeArray.getDof('ROTATION_Y');
            dofs([6 12 18 24]) = obj.nodeArray.getDof('ROTATION_Z');
        end
        
        function vals = getValuesVector(obj, step)
            vals = zeros(1,24);
            vals([1  7 13 19]) = obj.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals([2  8 14 20]) = obj.nodeArray.getDofValue('DISPLACEMENT_Y',step);
            vals([3  9 15 21]) = obj.nodeArray.getDofValue('DISPLACEMENT_Z',step);
            vals([4 10 16 22]) = obj.nodeArray.getDofValue('ROTATION_Z',step);
            vals([5 11 17 23]) = obj.nodeArray.getDofValue('ROTATION_X',step);
            vals([6 12 18 24]) = obj.nodeArray.getDofValue('ROTATION_Y',step);
        end
        
%         function vals = getFirstDerivativesVector(obj, step)
%             vals = zeros(1,24);
%             [~, vals([1  7 13 19]), ~] = obj.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
%             [~, vals([2  8 14 20]), ~] = obj.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
%             [~, vals([3  9 15 21]), ~] = obj.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
%             [~, vals([4 10 16 22]), ~] = obj.nodeArray.getDof('ROTATION_X').getAllValues(step);
%             [~, vals([5 11 17 23]), ~] = obj.nodeArray.getDof('ROTATION_Y').getAllValues(step);
%             [~, vals([6 12 18 24]), ~] = obj.nodeArray.getDof('ROTATION_Z').getAllValues(step);
%         end
%         
%         function vals = getSecondDerivativesVector(obj, step)
%             vals = zeros(1,24);
%             [~, ~, vals([1  7 13 19])] = obj.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
%             [~, ~, vals([2  8 14 20])] = obj.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
%             [~, ~, vals([3  9 15 21])] = obj.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
%             [~, ~, vals([4 10 16 22])] = obj.nodeArray.getDof('ROTATION_X').getAllValues(step);
%             [~, ~, vals([5 11 17 23])] = obj.nodeArray.getDof('ROTATION_Y').getAllValues(step);
%             [~, ~, vals([6 12 18 24])] = obj.nodeArray.getDof('ROTATION_Z').getAllValues(step);
%         end
        
    end
    
    methods (Access = private, Static = true)
        
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




