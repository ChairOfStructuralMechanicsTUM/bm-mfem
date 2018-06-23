classdef Porous2d4n < QuadrilateralElement
    properties (access=private)
    end
   
    methods
        % Constructor
        function PorousElement = Porous2d4n(id,nodeArray)
            
            requiredPropertyNames = cellstr(["FRAME_DENSITY", "FRAME_LAME_PARAMETER_LAMBDA","FRAME_LAME_PARAMETER_MU" ...
                "AMBIENT_FLUID_DENSITY", "AMBIENT_FLUID_VISCOSITY", "AMBIENT_FLUID_PRESSURE", "HEAT_CAPACITY_RATIO", "PRANDTL_NUMBER" ...
                "POROSITY", "TORTUOSITY", "STATIC_FLOW_RESISTIVITY", "VISCOUS_CHAR_LENGTH", "THERMAL_CHAR_LENGTH",...
                "NUMBER_GAUSS_POINT","FREQUENCY"]);
            
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 4 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            % call the super class constructor
            PorousElement@QuadrilateralElement(super_args{:});
            PorousElement.dofNames = cellstr(["FRAME_DISPLACEMENT_X", "FRAME_DISPLACEMENT_Y",
                "FLUID_DISPLACEMENT_X", "FLUID_DISPLACEMENT_Y"]);
        end
        
        % Getter functions
        
        function dofs = getDofList(element)
%             dofs([1 5 9 13]) = element.nodeArray.getDof('FRAME_DISPLACEMENT_X');
%             dofs([2 6 10 14]) = element.nodeArray.getDof('FRAME_DISPLACEMENT_Y');
%             dofs([3 7 11 15]) = element.nodeArray.getDof('FLUID_DISPLACEMENT_X');
%             dofs([4 8 12 16]) = element.nodeArray.getDof('FLUID_DISPLACEMENT_Y');
        end
        
        function vals = getValuesVector(element, step)
%             vals = zeros(1,8);
%             
%             vals([1 3 5 7]) = element.nodeArray.getDofValue('DISPLACEMENT_X',step);
%             vals([2 4 6 8]) = element.nodeArray.getDofValue('DISPLACEMENT_Y',step);
        end
        
        %Initialization
        
        function initialize(Porous2d4n)
            Porous2d4n.lengthX = computeLength(Porous2d4n.nodeArray(1).getCoords, ...
                Porous2d4n.nodeArray(2).getCoords);
            
            Porous2d4n.lengthY = computeLength(Porous2d4n.nodeArray(1).getCoords, ...
                Porous2d4n.nodeArray(4).getCoords);
            
            checkConvexity(Porous2d4n);
        end
       
        %Computation functions
        
        function [N_mat, N, Be, J] = computeShapeFunction(porousElement,xi,eta)
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
                ele_coords(i,1) = porousElement.nodeArray(i).getX;
                ele_coords(i,2) = porousElement.nodeArray(i).getY;
            end
            
            % Jacobian
            J = N_Diff_Par * ele_coords;
            
            % Calculation of B-Matrix
            B=J\N_Diff_Par;%/Jdet;
            Bx=B(1,1:4);
            By=B(2,1:4);
            
            Be=[Bx(1),0,Bx(2),0,Bx(3),0,Bx(4),0;
                0,By(1),0,By(2),0,By(3),0,By(4);
                By(1),Bx(1),By(2),Bx(2),By(3),Bx(3),By(4),Bx(4)];
        end
        
        function [K_ss,K_sf,K_fs,K_ff] = computeLocalStiffnessMatrix(porousElement)
            omega = porousElement.getPropertyValue('FREQUENCY');
            %Frame
            lambda = porousElement.getPropertyValue('FRAME_LAME_PARAMETER_LAMBDA');
            mu = porousElement.getPropertyValue('FRAME_LAME_PARAMETER_MU');
            %Fluid
            gamma = porousElement.getPropertyValue('HEAT_CAPACITY_RATIO');
            P_0 = porousElement.getPropertyValue('AMBIENT_FLUID_PRESSURE');
            eta = porousElement.getPropertyValue('AMBIENT_FLUID_VISCOSITY');
            roh_f = porousElement.getPropertyValue('AMBIENT_FLUID_DENSITY');
            Pr = porousElement.getPropertyValue('PRANDTL_NUMBER');
            %Porous
            lambda_apostroph = porousElement.getPropertyValue('THERMAL_CHAR_LENGTH');
            phi = porousElement.getPropertyValue('POROSITY');
            
            
            temp1 = sqrt(1+((1i*omega*Pr*lambda_apostroph*lambda_apostroph*roh_f)/(16*eta)));
            temp2 = (1+((8*eta)/(1i*omega*Pr*lambda_apostroph*lambda_apostroph*roh_f))*temp1)^-1;
            K_f = (gamma*P_0)/(gamma-(gamma-1)*temp2);
            
            R = phi * K_f;
            Q = (1-phi) * K_f;
            A = lambda + (((1-phi)^2)/phi) * K_f;
            A_star = A + 2*mu;
            
            D_s = [ A_star, A, 0;
                    A, A_star, 0;
                    0, 0, mu];
            D_sf = [Q, Q, 0;
                    Q, Q, 0;
                    0, 0, 0];
            D_f = [ R, R, 0;
                    R, R, 0;
                    0, 0, 0];
                
            p = porousElement.getPropertyValue('NUMBER_GAUSS_POINT');    
            K_ss = zeros(8,8);
            K_sf = K_ss;
            K_fs = K_ss;
            K_ff = K_ss;
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [~, ~, B, J] = computeShapeFunction(porousElement, xi, eta);
                    K_ss = K_ss+(w(i)*w(j)*det(J)*transpose(B)*(D_s*B));
                    K_sf = K_sf+(w(i)*w(j)*det(J)*transpose(B)*(D_sf*B));
                    K_fs = K_sf;
                    K_ff = K_ff+(w(i)*w(j)*det(J)*transpose(B)*(D_f*B));
                end
            end     
        end
        
        function [M_ss,M_sf,M_fs,M_ff] = computeLocalMassMatrix(porousElement)
            omega = porousElement.getPropertyValue('FREQUENCY');
            roh_s = porousElement.getPropertyValue('FRAME_DENSITY');
            alpha_infinity  = porousElement.getPropertyValue('TORTUOSITY');
            eta = porousElement.getPropertyValue('AMBIENT_FLUID_VISCOSITY');
            roh_f = porousElement.getPropertyValue('AMBIENT_FLUID_DENSITY');
            sigma = porousElement.getPropertyValue('STATIC_FLOW_RESISTIVITY');
            lambda = porousElement.getPropertyValue('VISCOUS_CHAR_LENGTH');
            phi = porousElement.getPropertyValue('POROSITY');
            
            temp1 = ((4*1i*omega*alpha_infinity*alpha_infinity*eta*roh_f)/(sigma*sigma*lambda*lambda*phi*phi));
            b = sigma*phi*phi*(sqrt(1+temp1));
            roh_a = phi*roh_f*(alpha_infinity-1);
            
            roh_sf_complex = (-1*roh_a) + 1i * b/omega;
            roh_f_complex = phi*roh_f - roh_sf_complex;
            roh_s_complex = (1-phi)*roh_s - roh_sf_complex;
            
            p = porousElement.getPropertyValue('NUMBER_GAUSS_POINT');
            M_ss=zeros(8,8);
            M_sf=M_ss;
            M_fs=M_ss;
            M_ff=M_ss;
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                xi=g(i);
                for j=1:p
                    eta=g(j);
                    [N_mat, ~, ~, J] = computeShapeFunction(porousElement,xi,eta);
                    M_ss = M_ss + (w(i)*w(j)*roh_s_complex*transpose(N_mat)*N_mat*det(J));
                    M_sf = M_sf + (w(i)*w(j)*roh_sf_complex*transpose(N_mat)*N_mat*det(J));
                    M_fs = M_sf;
                    M_ff = M_ff + (w(i)*w(j)*roh_f_complex*transpose(N_mat)*N_mat*det(J));
                end
            end     
        end
        
        function F = computeLocalForceVector(quadrilateralElement)
        end
    end
end
