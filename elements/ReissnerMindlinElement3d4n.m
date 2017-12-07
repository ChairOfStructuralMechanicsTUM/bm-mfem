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

        function [N_mat, N,B_b, B_s, J] = computeShapeFunction(reissnerMindlinElement3d4n,xi,eta)
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
       
            N_mat = sparse(3,12);
            N_mat(1,1:3:end) = N(:);
            N_mat(2,2:3:end) = N(:);
            N_mat(3,3:3:end) = N(:);
            
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
            alpha = 0.8601;     % shear correction factor

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

            dens_mat = sparse(3,3);
            dens_mat(1,1) = density*thickness; 
            dens_mat(2,2) = density*thickness^3/12; 
            dens_mat(3,3) = dens_mat(2,2); 

            massMatrix = sparse( 12,12);

            for xi = 1 : nr_gauss_points
                for eta = 1 : nr_gauss_points
                    [N_mat, ~,~,~,J] = computeShapeFunction(reissnerMindlinElement3d4n,g(xi),g(eta));
                    massMatrix = massMatrix + N_mat' * dens_mat * N_mat *det(J) * w(xi) * w(eta);
                end
            end
        end
        
        
        function massMatrix = computeLocalMassMatrix_PRZEMIENIECKI(reissnerMindlinElement3d4n)
            density = reissnerMindlinElement3d4n.getPropertyValue('DENSITY');
            thickness = reissnerMindlinElement3d4n.getPropertyValue('THICKNESS');
            a = reissnerMindlinElement3d4n.getLengthX();
            b = reissnerMindlinElement3d4n.getLengthY();
            V = a * b * thickness; 
            
            massMatrix = zeros( 12,12);
            
            massMatrix(1,1) = 24336; 
            
            massMatrix(2,1) = 3432 * b; 
            massMatrix(2,2) = 624 * b^2; 
            
            massMatrix(3,1) = -3432 * a; 
            massMatrix(3,2) = -484 * a * b; 
            massMatrix(3,3) = 624 * a^2; 
            
            massMatrix(4,1) = 8424; 
            massMatrix(4,2) = 2028 * b; 
            massMatrix(4,3) = -1188 * a; 
            massMatrix(4,4) = 24336; 
            
            massMatrix(5,1) = -2028 * b; 
            massMatrix(5,2) = -468 * b^2; 
            massMatrix(5,3) = 286 * a * b; 
            massMatrix(5,4) = -3432 * b; 
            massMatrix(5,5) = 624 * b^2; 
            
            massMatrix(6,1) = -1188 * a; 
            massMatrix(6,2) = -286 * a * b; 
            massMatrix(6,3) = 216 * a^2; 
            massMatrix(6,4) = -3432 * a;
            massMatrix(6,5) = 484 * a * b; 
            massMatrix(6,6) = 624 * a^2; 
            
            massMatrix(7,1) = 2916; 
            massMatrix(7,2) = 702 * b;
            massMatrix(7,3) = -702 * a; 
            massMatrix(7,4) = 8424; 
            massMatrix(7,5) = -1188 * b;
            massMatrix(7,6) = -2028 * a; 
            massMatrix(7,7) = 24336; 
            
            massMatrix(8,1) = -702 * b; 
            massMatrix(8,2) = -162 * b^2; 
            massMatrix(8,3) = 169 * a * b; 
            massMatrix(8,4) = -1188 * b;
            massMatrix(8,5) = 216 * b^2;
            massMatrix(8,6) = 286 * a * b; 
            massMatrix(8,7) = -3432 * b; 
            massMatrix(8,8) = 624 * b^2; 
            
            massMatrix(9,1) = 702 * a; 
            massMatrix(9,2) = 169 * a * b; 
            massMatrix(9,3) = -162 * a^2; 
            massMatrix(9,4) = 2028 * a; 
            massMatrix(9,5) = -286 * a * b; 
            massMatrix(9,6) = -468 * a^2; 
            massMatrix(9,7) = 3432 * a; 
            massMatrix(9,8) = -484 * a * b; 
            massMatrix(9,9) = 624 * a^2; 
            
            massMatrix(10,1) = 8424; 
            massMatrix(10,2) = 1188 * b; 
            massMatrix(10,3) = -2028 * a; 
            massMatrix(10,4) = 2916; 
            massMatrix(10,5) = -702 * b; 
            massMatrix(10,6) = -702 * a; 
            massMatrix(10,7) = 8424;
            massMatrix(10,8) = -2028 * b; 
            massMatrix(10,9) = 1188 * a; 
            massMatrix(10,10) = 24336;
            
            massMatrix(11,1) = 1188 * b; 
            massMatrix(11,2) = 216 * b^2; 
            massMatrix(11,3) = -286 * a * b; 
            massMatrix(11,4) = 702 * b; 
            massMatrix(11,5) = -162 * b^2; 
            massMatrix(11,6) = -169 * a * b; 
            massMatrix(11,7) = 2028 * b; 
            massMatrix(11,8) = -468 * b^2; 
            massMatrix(11,9) = 286 * a * b; 
            massMatrix(11,10) = 3432 * b; 
            massMatrix(11,11) = 624 * b^2;
            
            massMatrix(12,1) = 2028 * a; 
            massMatrix(12,2) = 286 * a * b; 
            massMatrix(12,3) = -468 * a^2;
            massMatrix(12,4) = 702 * a; 
            massMatrix(12,5) = -169 * a * b; 
            massMatrix(12,6) = -162 * a^2; 
            massMatrix(12,7) = 1188 * a; 
            massMatrix(12,8) = -286 * a * b; 
            massMatrix(12,9) = 216 * a^2; 
            massMatrix(12,10) = 3432 * a; 
            massMatrix(12,11) = 484 * a * b;
            massMatrix(12,12) = 624 * a^2; 
            
            massMatrix(1,2) = massMatrix(2,1);
            massMatrix(1,3) = massMatrix(3,1);
            massMatrix(1,4) = massMatrix(4,1);
            massMatrix(1,5) = massMatrix(5,1);
            massMatrix(1,6) = massMatrix(6,1);
            massMatrix(1,7) = massMatrix(7,1);
            massMatrix(1,8) = massMatrix(8,1);
            massMatrix(1,9) = massMatrix(9,1);
            massMatrix(1,10) = massMatrix(10,1);
            massMatrix(1,11) = massMatrix(11,1);
            massMatrix(1,12) = massMatrix(12,1);
            
            massMatrix(2,3) = massMatrix(3,2);
            massMatrix(2,4) = massMatrix(4,2);
            massMatrix(2,5) = massMatrix(5,2);
            massMatrix(2,6) = massMatrix(6,2);
            massMatrix(2,7) = massMatrix(7,2);
            massMatrix(2,8) = massMatrix(8,2);
            massMatrix(2,9) = massMatrix(9,2);
            massMatrix(2,10) = massMatrix(10,2);
            massMatrix(2,11) = massMatrix(11,2);
            massMatrix(2,12) = massMatrix(12,2);
            
            massMatrix(3,4) = massMatrix(4,3);
            massMatrix(3,5) = massMatrix(5,3);
            massMatrix(3,6) = massMatrix(6,3);
            massMatrix(3,7) = massMatrix(7,3);
            massMatrix(3,8) = massMatrix(8,3);
            massMatrix(3,9) = massMatrix(9,3);
            massMatrix(3,10) = massMatrix(10,3);
            massMatrix(3,11) = massMatrix(11,3);
            massMatrix(3,12) = massMatrix(12,3);
            
            massMatrix(4,5) = massMatrix(5,4);
            massMatrix(4,6) = massMatrix(6,4);
            massMatrix(4,7) = massMatrix(7,4);
            massMatrix(4,8) = massMatrix(8,4);
            massMatrix(4,9) = massMatrix(9,4);
            massMatrix(4,10) = massMatrix(10,4);
            massMatrix(4,11) = massMatrix(11,4);
            massMatrix(4,12) = massMatrix(12,4);
            
            massMatrix(5,6) = massMatrix(6,5);
            massMatrix(5,7) = massMatrix(7,5);
            massMatrix(5,8) = massMatrix(8,5);
            massMatrix(5,9) = massMatrix(9,5);
            massMatrix(5,10) = massMatrix(10,5);
            massMatrix(5,11) = massMatrix(11,5);
            massMatrix(5,12) = massMatrix(12,5);
            
            massMatrix(6,7) = massMatrix(7,6);
            massMatrix(6,8) = massMatrix(8,6);
            massMatrix(6,9) = massMatrix(9,6);
            massMatrix(6,10) = massMatrix(10,6);
            massMatrix(6,11) = massMatrix(11,6);
            massMatrix(6,12) = massMatrix(12,6);
            
            massMatrix(7,8) = massMatrix(8,7);
            massMatrix(7,9) = massMatrix(9,7);
            massMatrix(7,10) = massMatrix(10,7);
            massMatrix(7,11) = massMatrix(11,7);
            massMatrix(7,12) = massMatrix(12,7);
            
            massMatrix(8,9) = massMatrix(9,8);
            massMatrix(8,10) = massMatrix(10,8);
            massMatrix(8,11) = massMatrix(11,8);
            massMatrix(8,12) = massMatrix(12,8);
            
            massMatrix(9,10) = massMatrix(10,9);
            massMatrix(9,11) = massMatrix(11,9);
            massMatrix(9,12) = massMatrix(12,9);
            
            massMatrix(10,11) = massMatrix(11,10);
            massMatrix(10,12) = massMatrix(12,10);
            
            massMatrix(11,12) = massMatrix(12,11);
            
            massMatrix = massMatrix * ((density * V)/176400);
            
            
        end
        
        function D = computeLocalDampingMatrix(e)
            D = zeros(12,12);
        end
   
    end
end
