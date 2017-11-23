classdef ReissnerMindlinElement3d4n < PlateElement 
    %REISSNERMINDLINPLATE  Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = protected)
    end 
    
    methods
        % Constructor
        function reissnerMindlinElement3d4n = ReissnerMindlinElement3d4n(id,nodeArray)
            
            requiredPropertyNames = cellstr(["YOUNGS_MODULUS", "SHEAR_MODULUS", ...
                                            "POISSON_RATIO", "THICKNESS",...
                                            "NUMBER_GAUSS_POINT", "DENSITY"]); 
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
            reissnerMindlinElement3d4n.dofNames = cellstr(["DISPLACEMENT_Z", ...
                                                "ROTATION_X", "ROTATION_Y"]);
        end
        function initialize(reissnerMindlinElement3d4n)
            reissnerMindlinElement3d4n.lengthX = computeLength(reissnerMindlinElement3d4n.nodeArray(1).getCoords, ...
                reissnerMindlinElement3d4n.nodeArray(2).getCoords);
            reissnerMindlinElement3d4n.lengthY = computeLength(reissnerMindlinElement3d4n.nodeArray(1).getCoords, ...
                reissnerMindlinElement3d4n.nodeArray(4).getCoords);
            
        end



        function responseDoF = getResponseDofArray(plateElement, step)
            responseDoF = zeros(12,1);
            
            for itNodes = 1:1:4
                nodalDof = plateElement.nodeArray(itNodes).getDofArray;
                nodalDof = nodalDof.';
                
                for itDof = 3:(-1):1
                    responseDoF(3*itNodes-(itDof-1),1) = nodalDof(4-itDof).getValue(step);
                end
            end
        end

        function [N_mass, N,B_b, B_s, J] = computeShapeFunction(reissnerMindlinElement3d4n,xi,eta)
            % Shape Function and Derivatives                    
            N(1) = (1-xi)*(1-eta)/4;
            N_Diff_Par(1,1) = -(1-eta)/4;
            N_Diff_Par(2,1) = -(1-xi)/4;

            N(2) = (1+xi)*(1-eta)/4;
            N_Diff_Par(1,2) = (1-eta)/4;
            N_Diff_Par(2,2) = -(1+xi)/4;

            N(3) = (1+xi)*(1+eta)/4;
            N_Diff_Par(1,3) = (1+eta)/4;
            N_Diff_Par(2,3) = (1+xi)/4;

            N(4) = (1-xi)*(1+eta)/4;
            N_Diff_Par(1,4) = -(1+eta)/4;
            N_Diff_Par(2,4) = (1-xi)/4;
       
            N_mass = sparse(3,12);
            N_mass(1,1:3:end) = N(:);
            N_mass(2,2:3:end) = N(:);
            N_mass(3,3:3:end) = N(:);
            
            coords = zeros(4,2); 
            for i=1:4
                coords(i,1) = reissnerMindlinElement3d4n.nodeArray(i).getX;
                coords(i,2) = reissnerMindlinElement3d4n.nodeArray(i).getY;
            end
            
            % Jacobian 
            J = N_Diff_Par * coords; 
            
            N_Diff = J \ N_Diff_Par;
            
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

        
        function stiffnessMatrix = computeLocalStiffnessMatrix(reissnerMindlinElement3d4n)
            
            Emodul = reissnerMindlinElement3d4n.getPropertyValue('YOUNGS_MODULUS');
            Gmodul = reissnerMindlinElement3d4n.getPropertyValue('SHEAR_MODULUS');
            poisson_ratio = reissnerMindlinElement3d4n.getPropertyValue('POISSON_RATIO');
            nr_gauss_points = reissnerMindlinElement3d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            thickness = reissnerMindlinElement3d4n.getPropertyValue('THICKNESS');
            alpha = 5/6;     % shear correction factor
            
            % Moment-Curvature Equations
            D_b = zeros(3,3);
            D_b(1,1) = 1;
            D_b(1,2) = poisson_ratio;
            D_b(2,1) = D_b(1,2);
            D_b(2,2) = 1;
            D_b(3,3) = (1-poisson_ratio)/2;
             
            K = (Emodul * thickness^3) / (12*(1-poisson_ratio^2));             
%              % K Matrix from Fellipa
%             K = (Emodul * thickness^3)/(1-poisson_ratio^2); 
            D_b = D_b* K;

            %Shear Equation
            D_s = eye(2) * alpha * Gmodul * thickness; 
                
            stiffnessMatrix = sparse(12,12);
            [w,g] = returnGaussPoint(nr_gauss_points);
            
            for xi=1:nr_gauss_points
                for eta=1:nr_gauss_points
                    [~, ~,B_b, B_s, J] = computeShapeFunction(reissnerMindlinElement3d4n,g(xi),g(eta));
                    
                    stiffnessMatrix = stiffnessMatrix + thickness * ...
                        B_b' * D_b * B_b *det(J) * w(xi) * w(eta) + ...
                        B_s' * D_s * B_s * det(J) * w(xi) * w(eta);

                end
            end
        end
        
        
        
        function massMatrix = computeLocalMassMatrix(reissnerMindlinElement3d4n)
            density = reissnerMindlinElement3d4n.getPropertyValue('DENSITY');
            thickness = reissnerMindlinElement3d4n.getPropertyValue('THICKNESS');
            nr_gauss_points = reissnerMindlinElement3d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            [w,g] = returnGaussPoint(nr_gauss_points);
            
%             a=reissnerMindlinElement3d4n.lengthY;
%             b=reissnerMindlinElement3d4n.lengthX;
%             a_w(1,1) = @(xi,eta)(1+2*xi)*(1-xi)^2*(1+2*eta)*(1-eta)^2;
%             a_w(1,2) = @(xi,eta)(1+2*xi)*(1-xi)^2*eta*(1-eta)^2*b;
%             a_w(1,3) = @(xi,eta)-xi*(1-xi)^2*(1+2*eta)*((1-eta)^2)*a;
%             a_w(1,4) = @(xi,eta)(1+2*xi)*(1-xi)^2*(3-2*eta)*eta^2;
%             a_w(1,5) = @(xi,eta)-(1+2*xi)*(1-xi)^2*(1-eta)*eta^2*b;
%             a_w(1,6) = @(xi,eta)-xi*(1-xi)^2*(3-2*eta)*eta^2*a;
%             a_w(1,7) = @(xi,eta)(3-2*xi)*xi^2*(3-2*eta)*eta^2;
%             a_w(1,8) = @(xi,eta)-(3-2*xi)*xi^2*(1-eta)*eta^2*b;
%             a_w(1,9) = @(xi,eta)(1-xi)*xi^2*(3-2*eta)*eta^2*a;
%             a_w(1,10) = @(xi,eta)(3-2*xi)*xi^2*(1+2*eta)*(1-eta)^2;
%             a_w(1,11) = @(xi,eta)(3-2*xi)*xi^2*eta*(1-eta)^2*b;
%             a_w(1,12) = @(xi,eta)(1-xi)*xi^2*(1+2*eta)*(1-eta)^2*a;

            
            dens_mat = zeros(3,3);
            dens_mat(1,1) = density*thickness; 
            dens_mat(2,2) = density*thickness^3/12; 
            dens_mat(3,3) = density*thickness^3/12; 
            
            massMatrix = zeros(12,12);
            
            for xi=1:nr_gauss_points
                for eta=1:nr_gauss_points
                    [N_mass, ~,~,~,J] = computeShapeFunction(reissnerMindlinElement3d4n,g(xi),g(eta));
                    
                    massMatrix = massMatrix + N_mass' * dens_mat * N_mass *det(J) * w(xi) * w(eta);
                end
            end
          
            disp(massMatrix);
            
        end
   
    end
end
