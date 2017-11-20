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
                                            "POISSON_RATIO", "THICKNESS",...
                                            "NUMBER_GAUSS_POINT"]); 
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
        end

%         getter functions
        
        function getThickness(reissnerMindlinElement3d4n)
            reissnerMindlinElement3d4n.thickness = reissnerMindlinElement3d4n.getPropertyValue('THICKNESS');
        end
        
%         function [N,dNdx, dNdy] = computeShapeFunction(reissnerMindlinElement3d4n, xi,eta)
%             %Shape Functions and derivatives
%             
%             N = {@(xi,eta) (1-xi)*(1-eta)/4, @(xi,eta) (1+xi)*(1-eta)/4, @(xi,eta) (1+xi)*(1+eta)/4, @(xi,eta) (1-xi)*(1+eta)/4};
%             dNdxi = {@(xi,eta) -(1-eta)/4,  @(xi,eta) (1-eta)/4, @(xi,eta) (1+eta)/4, @(xi,eta) -(1+eta)/4};
%             dNdeta = {@(xi,eta) -(1-xi)/4, @(xi,eta) -(1+xi)/4, @(xi,eta) (1+xi)/4, @(xi,eta) (1-xi)/4 };
%             
%             coords = zeros(2,4);
%             for i=1:4 
%                 coords(1,i) = reissnerMindlinElement3d4n.nodeArray(i).getX;
%                 coords(2,i) = reissnerMindlinElement3d4n.nodeArray(i).getY;
%             end
%             
%             for i=1:4
%             J{1,1} = @(xi,eta) dNdxi{i}*coords(1,i);
%             J{1,2} = @(xi,eta) dNdxi{i}*coords(2,i);
%             J{2,1} = @(xi,eta) dNdeta{i}*coords(1,i);
%             J{2,2} = @(xi,eta) dNdeta{i}*coords(2,i);
%             end
%             
%             dNdx = inv(J)*dNdxi;
%             dNdy = inv(J)*dNdeta; 
%             
%            
%         end
        
        function computeLocalStiffnessMatrix(reissnerMindlinElement3d4n)
            syms xi eta;
            N(1) = 1/4 * (1-xi) * (1-eta);
            N(2) = 1/4 * (1+xi) * (1-eta);
            N(3) = 1/4 * (1+xi) * (1+eta);
            N(4) = 1/4 * (1-xi) * (1+eta);
            
            coords = zeros(2,4);
            for i=1:4 
                coords(1,i) = reissnerMindlinElement3d4n.nodeArray(i).getX;
                coords(2,i) = reissnerMindlinElement3d4n.nodeArray(i).getY;
            end
            
            x = 0; y = 0;
            for i=1:4
            x = x + N(i) * coords(1,i);
            y = y + N(i) * coords(2,i);
            end

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
            dN1dx = diff(N(1),xi) * dxi_dx + diff(N(1),eta)* deta_dx;
            dN1dy = diff(N(1),xi) * dxi_dy + diff(N(1),eta) * deta_dy;
            dN2dx = diff(N(2),xi) * dxi_dx + diff(N(2),eta) * deta_dx;
            dN2dy = diff(N(2),xi) * dxi_dy + diff(N(2),eta) * deta_dy;
            dN3dx = diff(N(3),xi) * dxi_dx + diff(N(3),eta) * deta_dx;
            dN3dy = diff(N(3),xi) * dxi_dy + diff(N(3),eta) * deta_dy;
            dN4dx = diff(N(4),xi) * dxi_dx + diff(N(4),eta) * deta_dx;
            dN4dy = diff(N(4),xi) * dxi_dy + diff(N(4),eta) * deta_dy;
            
            poisson_ratio = reissnerMindlinElement3d4n.getPropertyValue('POISSON_RATIO');
            
            % Moment-Curvature Equations
            D_b = zeros(3,3);
            D_b(1,1) = 1;
            D_b(1,2) = poisson_ratio;
            D_b(2,1) = D_b(1,2);
            D_b(2,2) = 1;
            D_b(3,3) = (1-poisson_ratio)/2;
            
            K = (reissnerMindlinElement3d4n.getPropertyValue('YOUNGS_MODULUS') * ... 
                 reissnerMindlinElement3d4n.thickness^3)/...
                 (12*(1-poisson_ratio^2));
             
             
%             K = (reissnerMindlinElement3d4n.getPropertyValue('YOUNGS_MODULUS') * ... 
%                  reissnerMindlinElement3d4n.getPropertyValue('THICKNESS')^3)/...
%                  (1-reissnerMindlinElement3d4n.getPropertyValue('POISSON_RATIO')^2);
            D_b = D_b* K;

            % constructing the bending B matrix: B_b
            
            B_b = [0, dN1dx,  0,    0,  dN2dx,    0,      0,   dN3dx,    0,      0,  dN4dx,     0   
                   0,   0,   dN1dy, 0,    0,    dN2dy,    0,     0,    dN3dy,    0,    0,    dN4dy
                   0, dN1dy, dN1dx, 0,  dN2dy,  dN2dx,    0,   dN3dy,  dN3dx,    0   dN4dy,  dN4dx];
         
            
            BbT_Db_Bb = (B_b'*D_b)*B_b * det(Jacobian);

            %Shear Equation
            
            alpha = 5/6;     % shear correction factor
            D_s = sym(eye(2));
            D_s = D_s * alpha*reissnerMindlinElement3d4n.getPropertyValue('SHEAR_MODULUS')*...
                    reissnerMindlinElement3d4n.thickness;             
          
            % constructing the shear B matrix: B_s
           
            B_s = [dN1dx, N(1),  0,   dN2dx,  N(2),  0,   dN3dx,   N(3),   0,   dN4dx,   N(4),   0
                   dN1dy,  0,   N(1), dN2dy,   0,   N(2), dN3dy,    0,    N(3), dN4dy,    0,    N(4)];               
      
            BsT_Ds_Bs = (B_s'*D_s)*B_s *det(Jacobian);
            
            Ke = zeros(12,12);
            
            p = reissnerMindlinElement3d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            
            [w,g] = returnGaussPoint(p);
            for i=1:p
                for j=1:p
                    
                    Ke = Ke + eval(subs(BbT_Db_Bb,[xi,eta],[g(i),g(j)]) *w(i)*w(j));
                    Ke = Ke + eval(subs(BsT_Ds_Bs,[xi,eta],[g(i),g(j)]) *w(i)*w(j));
            
                end
            end
           
        end
        

        
    end
end
