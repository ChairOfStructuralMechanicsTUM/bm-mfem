classdef ReissnerMindlinElement3d4n < PlateElement 
    %REISSNERMINDLINPLATE  Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
        thickness
    end 
    
    methods
        % Constructor
        function reissnerMindlinElement3d4n = ReissnerMindlinElement3d4n(id,nodeArray)
            
            requiredPropertyNames = cellstr(["YOUNGS_MODULUS", "SHEAR_MODULUS", ...
                                            "POISSON_RATIO", "THICKNESS"]); 
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
            reissnerMindlinElement3d4n@PlateElement(super_args{:});
            reissnerMindlinElement3d4n.dofNames = cellstr(['DISPLACEMENT_Z', 'ROTATION_X', 'ROTATION_Y']);
        end
        function initialize(element)
            element.lengthX = computeLength(element.nodeArray(1).getCoords, ...
                element.nodeArray(2).getCoords);
            element.lengthY = computeLength(element.nodeArray(1).getCoords, ...
                element.nodeArray(4).getCoords);
            fprintf("element.lengthX", element.lengthX);
            fprintf("element.lengthX", element.lengthY);
        end
                

        
        
        % getter functions
        
%         function thickness = getThickness(reissnerMindlinElement3d4n)
%             reissnerMindlinElement3d4n.thickness = reissnerMindlinElement3d4n.getPropertyValue('THICKNESS');
%             
%         end
%         
%         function lengthElementX = getNumElementsX(reissnerMindlinElement3d4n)
%             reissnerMindlinElement3d4n.numElementsX = reissnerMindlinElement3d4n.getPropertyValue('NUM_ELEMENTS_X');
%             lengthElementX = reissnerMindlinElement3d4n.lengthX / reissnerMindlinElement3d4n.numElementsX;
%         end
%         
%         function lengthElementY = getNumElementsY(reissnerMindlinElement3d4n)
%             reissnerMindlinElement3d4n.numElementsY = reissnerMindlinElement3d4n.getPropertyValue('NUM_ELEMENTS_Y');
%             lengthElementY = reissnerMindlinElement3d4n.lengthY / reissnerMindlinElement3d4n.numElementsY;
%         end
        



        function computeLocalStiffnessMatrix(reissnerMindlinElement3d4n)
    
            %Shape Functions
            syms xi eta;
            N1 = 1/4 * (1-xi) * (1-eta);
            N2 = 1/4 * (1+xi) * (1-eta);
            N3 = 1/4 * (1+xi) * (1+eta);
            N4 = 1/4 * (1-xi) * (1+eta);
            
            x = N1*x1 + N2*x2 + N3*x3 + N4*x4;
            y = N1*y1 + N2*y2 + N3*y3 + N4*y4;
            
            %Jacobian and Inverse Jacobian
            Jacobian = zeros(2);
            Jacobian(1,1) = diff(x, xi); 
            Jacobian(1,2) = diff(y, xi);
            Jacobian(1,2) = diff(x, eta);
            Jacobian(2,2) = diff(y, eta);
            
            invJacobian = inv(Jacobian);
            
            dxi_dx = invJacobian(1,1);
            deta_dx = invJacobian(1,2);
            dxi_dy = invJacobian(2,1);
            deta_dy = invJacobian(2,2);
            
            
            %Strain Displacement Matrix (B-Operator)
            dN1dx = diff(N1,xi) * dxi_dx + diff(N1,eta)* deta_dx;
            dN1dy = diff(N1,xi) * dxi_dy + diff(N1,eta) * deta_dy;
            dN2dx = diff(N2,xi) * dxi_dx + diff(N2,eta) * deta_dx;
            dN2dy = diff(N2,xi) * dxi_dy + diff(N2,eta) * deta_dy;
            dN3dx = diff(N3,xi) * dxi_dx + diff(N3,eta) * deta_dx;
            dN3dy = diff(N3,xi) * dxi_dy + diff(N3,eta) * deta_dy;
            dN4dx = diff(N4,xi) * dxi_dx + diff(N4,eta) * deta_dx;
            dN4dy = diff(N4,xi) * dxi_dy + diff(N4,eta) * deta_dy;
            
            % Moment-Curvature Equations
            D_b = sym(zeros(3,3));
            D_b(1,1) = 1;
            D_b(1,2) = reissnerMindlinElement3d4n.getPropertyValue('POISSON_RATIO');
            D_b(2,1) = D_b(1,2);
            D_b(2,2) = 1;
            D_b(3,3) = (1-reissnerMindlinElement3d4n.getPropertyValue('POISSON_RATIO'))/2;
            
            K = (reissnerMindlinElement3d4n.getPropertyValue('YOUNGS_MODULUS') * ... 
                 reissnerMindlinElement3d4n.getPropertyValue('THICKNESS')^3)/...
                 (12*(1-reissnerMindlinElement3d4n.getPropertyValue('POISSON_RATIO')^2));
             
            D_b = D_b* K; 

            % constructing the bending B matrix: B_b
            B_b = sym(zeros(3,12));
            
            B_b(1,2) = dN1dx;
            B_b(3,2) = dN1dy; 
            B_b(2,3) = dN1dy;
            B_b(3,3) = dN1dx;  
 
            B_b(1,5) = dN2dx;
            B_b(3,5) = dN2dy; 
            B_b(2,6) = dN2dy;
            B_b(3,6) = dN2dx;
            
            B_b(1,8) = dN3dx;
            B_b(3,8) = dN3dy; 
            B_b(2,9) = dN3dy;
            B_b(3,9) = dN3dx;

            B_b(1,11) = dN4dx;
            B_b(3,11) = dN4dy; 
            B_b(2,12) = dN4dy;
            B_b(3,12) = dN4dx;  

            BbT_Db_Bb = (transpose(B_b)*D_b)*B_b * det(Jacobian);

            %Shear Equation
            
            alpha = 5/6;     % shear correction factor
            D_s = sym(eye(2));
            D_s = D_s * alpha*reissnerMindlinElement3d4n.getPropertyValue('SHEAR_MODULUS')*...
                    reissnerMindlinElement3d4n.getPropertyValue('THICKNESS'); 
            
           % constructing the shear B matrix: B_s
           
            B_s = sym(zeros(2,12));
            B_s(1,1) = dN1dx;
            B_s(2,1) = dN1dy;
            B_s(1,2) = N1;
            B_s(2,2) = 0;
            B_s(1,3) = 0;
            B_s(2,3) = N1;

            B_s(1,4) = dN2dx;
            B_s(2,4) = dN2dy;
           
            B_s(1,5) = N2;
            B_s(2,5) = 0;
            B_s(1,6) = 0;
            B_s(2,6) = N2;

            B_s(1,7) = dN3dx;
            B_s(2,7) = dN3dy;
            B_s(1,8) = N3;
            B_s(2,8) = 0;
            B_s(1,9) = 0;
            B_s(2,9) = N3;

            B_s(1,10) = dN4dx;
            B_s(2,10) = dN4dy;
            B_s(1,11) = N4;
            B_s(2,11) = 0;
            B_s(1,12) = 0;
            B_s(2,12) = N4;
            
            BsT_Ds_Bs = (transpose(B_s)*D_s)*B_s *det(Jacobian);
            
            
            BbT_Db_Bb_xi = int(BbT_Db_Bb,xi,-1,1);
            BsT_Ds_Bs_xi = int(BsT_Ds_Bs,xi,-1,1);
            
            Ke = zeros(12,12);
            
            Ke = int(BbT_Db_Bb_xi,eta,-1,1) + int(BsT_Ds_Bs_xi,eta -1, 1)
            
            
        end
    end
end
