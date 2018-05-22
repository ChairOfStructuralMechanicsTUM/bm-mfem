classdef PorousElement3d8n < Element  
    %   POROUSELEMENT3D8N Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = private)
        lengthX
        lengthY
        lenghtZ
    end
    
    methods
        % constructor
        function porous3d8n = PorousElement3d8n(id, nodeArray)
            
            requiredPropertyNames = cellstr(["DENSITY_SOLID", "LAMBDA_SOLID", "MUE_SOLID", "DAMPING_SOLID", ...
                "DENSITY_FLUID", "VISCOSITY_FLUID", "STANDARD_PRESSURE_FLUID", "HEAT_CAPACITY_FLUID", "PRANDTL_NUMBER_FLUID", ...
                "POROSITY", "TORTUOSITY", "FLOW_RESISTIVITY", "VISCOUS_LENGHT", "THERMAL_LENGTH", ...
                "NUMBER_GAUSS_POINT"]);
            
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                if ~(length(nodeArray) == 8 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            % call the super class constructor
            porous3d8n@Element(super_args{:});
            porous3d8n.dofNames = cellstr(['DISPLACEMENT_SOLID_X'; 'DISPLACEMENT_SOLID_Y'; 'DISPLACEMENT_SOLID_Z'; ...
                'DISPLACEMENT_FLUID_X'; 'DISPLACEMENT_FLUID_Y'; 'DISPLACEMENT_FLUID_Z']);
            
        end
        
        function c = barycenter(hexahedron3d8n)
            diag1X = [porous3d8n.nodeArray(1).getX() porous3d8n.nodeArray(3).getX()];
            diag1Y = [porous3d8n.nodeArray(1).getY() porous3d8n.nodeArray(3).getY()];
            diag2X = [porous3d8n.nodeArray(2).getX() porous3d8n.nodeArray(4).getX()];
            diag2Y = [porous3d8n.nodeArray(2).getY() porous3d8n.nodeArray(4).getY()];
            [c(1),c(2)] = polyxpoly(diag1X, diag1Y, diag2X, diag2Y);
            c(3) = 0.5* (porous3d8n.nodeArray(1).getZ()+ porous3d8n.nodeArray(5).getZ());
        end
        
        
        % Check Convexity of quad
        function checkConvexity(porous3d8n)
            try
                [~] = porous3d8n.barycenter();
            catch
                error('Element %i is not convex', porous3d8n.getId());
            end
        end
        
        %Initialization
        function initialize(porous3d8n)
            porous3d8n.lengthX = computeLength(porous3d8n.nodeArray(1).getCoords, ...
                porous3d8n.nodeArray(2).getCoords);
            
            porous3d8n.lengthY = computeLength(porous3d8n.nodeArray(1).getCoords, ...
                porous3d8n.nodeArray(4).getCoords);
            
            porous3d8n.lengthY = computeLength(porous3d8n.nodeArray(1).getCoords, ...
                porous3d8n.nodeArray(5).getCoords);
            
            checkConvexity(porous3d8n);
        end

        % member functions
        function [N, Be, Jdet] = computeShapeFunction(porous3d8n, zeta, eta, mue)
            
            % Calculation of shape functions
            Nf = 1/8 * [(1-zeta)*(1-eta)*(1-mue);(1+zeta)*(1-eta)*(1-mue);(1+ zeta)*(1+eta)*(1-mue);(1-zeta)*(1+eta)*(1-mue);(1- zeta)*(1-eta)*(1+mue);(1+zeta)*(1-eta)*(1+mue);(1+ zeta)*(1+eta)*(1+mue);(1-zeta)*(1+eta)*(1+mue)];
            N=zeros(3,3*8);
            
            for l=1:8
                N(1,3*(l-1)+1)=Nf(l);
                N(2,3*(l-1)+2)=Nf(l);
                N(3,3*(l-1)+3)=Nf(l);
            end
            
            % Calculation of shape function derivatives
            dNzeta = 1/8 *[(-1)*(1-eta)*(1-mue),(1-eta)*(1-mue),(1+eta)*(1-mue),(-1)*(1+eta)*(1-mue),(-1)*(1-eta)*(1+mue),(1-eta)*(1+mue),(1+eta)*(1+mue),(-1)*(1+eta)*(1+mue)];
            dNeta = 1/8 *[(1- zeta)*(-1)*(1-mue),(1+zeta)*(-1)*(1-mue),(1+ zeta)*(1-mue),(1-zeta)*(1-mue),(1- zeta)*(-1)*(1+mue),(1+zeta)*(-1)*(1+mue),(1+ zeta)*(1+mue),(1- zeta)*(1+mue)];
            dNmue = 1/8 *[(1-zeta)*(1-eta)*(-1),(1+zeta)*(1-eta)*(-1),(1+zeta)*(1+eta)*(-1),(1-zeta)*(1+eta)*(-1),(1- zeta)*(1-eta),(1+zeta)*(1-eta),(1+ zeta)*(1+eta),(1-zeta)*(1+eta)];
            
            % Calculation of Jacobian matrix
            coords=zeros(8,3);
            for i=1:1:8
                coords(i,1:3)=porous3d8n.nodeArray(i).getCoords;
            end
            
            xn=coords(1:8,1);
            yn=coords(1:8,2);
            zn=coords(1:8,3);
            
            J11=dNzeta*xn;
            J12=dNzeta*yn;
            J13=dNzeta*zn;
            J21=dNeta*xn;
            J22=dNeta*yn;
            J23=dNeta*zn;
            J31=dNmue*xn;
            J32=dNmue*yn;
            J33=dNmue*zn;
            
            J=[J11, J12, J13; J21, J22, J23; J31, J32, J33];
            Jdet=det(J);
            Jinv=[J22*J33-J32*J23, J32*J13-J12*J33, J12*J23-J22*J13;
                J31*J23-J21*J33, J11*J33-J31*J13, J21*J13-J11*J23;
                J21*J32-J31*J22, J31*J12-J11*J32, J11*J22-J21*J12 ];
            
            % Calculation of B-Matrix
            B=Jinv*[dNzeta;dNeta;dNmue]/Jdet;
            Bx=B(1,1:8);
            By=B(2,1:8);
            Bz=B(3,1:8);
            
            Be=[Bx(1),0,0,Bx(2),0,0,Bx(3),0,0,Bx(4),0,0,Bx(5),0,0,Bx(6),0,0,Bx(7),0,0,Bx(8),0,0;
                0,By(1),0,0,By(2),0,0,By(3),0,0,By(4),0,0,By(5),0,0,By(6),0,0,By(7),0,0,By(8),0;
                0,0,Bz(1),0,0,Bz(2),0,0,Bz(3),0,0,Bz(4),0,0,Bz(5),0,0,Bz(6),0,0,Bz(7),0,0,Bz(8);
                By(1),Bx(1),0,By(2),Bx(2),0,By(3),Bx(3),0,By(4),Bx(4),0,By(5),Bx(5),0,By(6),Bx(6),0,By(7),Bx(7),0,By(8),Bx(8),0;
                0,Bz(1),By(1),0,Bz(2),By(2),0,Bz(3),By(3),0,Bz(4),By(4),0,Bz(5),By(5),0,Bz(6),By(6),0,Bz(7),By(7),0,Bz(8),By(8);
                Bz(1),0,Bx(1),Bz(2),0,Bx(2),Bz(3),0,Bx(3),Bz(4),0,Bx(4),Bz(5),0,Bx(5),Bz(6),0,Bx(6),Bz(7),0,Bx(7),Bz(8),0,Bx(8)];
            
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(porous3d8n, omega)
            % Read Material Properties
            DENSITY_SOLID = porous3d8n.getPropertyValue('DENSITY_SOLID');
            LAMBDA_SOLID = porous3d8n.getPropertyValue('LAMBDA_SOLID');
            MUE_SOLID = porous3d8n.getPropertyValue('MUE_SOLID');
            DAMPING_SOLID = porous3d8n.getPropertyValue('DAMPING_SOLID');
            
            DENSITY_FLUID = porous3d8n.getPropertyValue('DENSITY_FLUID');
            VISCOSITY_FLUID = porous3d8n.getPropertyValue('VISCOSITY_FLUID');
            STANDARD_PRESSURE_FLUID = porous3d8n.getPropertyValue('STANDARD_PRESSURE_FLUID');
            HEAT_CAPACITY_FLUID = porous3d8n.getPropertyValue('HEAT_CAPACITY_FLUID');
            PRANDTL_NUMBER_FLUID = porous3d8n.getPropertyValue('PRANDTL_NUMBER_FLUID');
            
            POROSITY = porous3d8n.getPropertyValue('POROSITY');
            TORTUOSITY = porous3d8n.getPropertyValue('TORTUOSITY');
            FLOW_RESISTIVITY = porous3d8n.getPropertyValue('FLOW_RESISTIVITY');
            VISCOUS_LENGHT = porous3d8n.getPropertyValue('VISCOUS_LENGHT');
            THERMAL_LENGTH = porous3d8n.getPropertyValue('THERMAL_LENGTH');
            
            p = porous3d8n.getPropertyValue('NUMBER_GAUSS_POINT');
            
            % Determine Coefficients for EMat
            LAMBDA = (1 + 1i * DAMPING_SOLID) * LAMBDA_SOLID;
            MUE = (1 + 1i * DAMPING_SOLID) * MUE_SOLID;
            K = HEAT_CAPACITY_FLUID * STANDARD_PRESSURE_FLUID /(HEAT_CAPACITY_FLUID - ...
                (1 - HEAT_CAPACITY_FLUID)*(1 + 8 * VISCOSITY_FLUID/(1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID)* ...
                (1 + 1i * omega * PRANDTL_NUMBER_FLUID * THERMAL_LENGTH^2 * DENSITY_FLUID /(16 * VISCOSITY_FLUID))^0.5));
            R = POROSITY * K;
            Q = (1 - POROSITY) * K;
            A = LAMBDA + ((1 - POROSITY)^2/ POROSITY) * K; 
            
            % Matrix valid for isotropic material only

            Emat=[A + 2 * MUE,A,A,0,0,0,Q,Q,Q,0,0,0; ...
                A,A + 2 * MUE,A,0,0,0,Q,Q,Q,0,0,0; ...
                A,A,A + 2 * MUE,0,0,0,Q,Q,Q,0,0,0; ...
                0,0,0,MUE,0,0,0,0,0,0,0,0; ...
                0,0,0,0,MUE,0,0,0,0,0,0,0; ...
                0,0,0,0,0,MUE,0,0,0,0,0,0; ...
                Q,Q,Q,0,0,0,R,R,R,0,0,0; ...
                Q,Q,Q,0,0,0,R,R,R,0,0,0; ...
                Q,Q,Q,0,0,0,R,R,R,0,0,0; ...
                0,0,0,0,0,0,0,0,0,0,0,0; ...
                0,0,0,0,0,0,0,0,0,0,0,0; ...
                0,0,0,0,0,0,0,0,0,0,0,0];
            
            
            stiffnessMatrix=zeros(48,48);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                zeta=g(i);
                for j=1:p
                    eta=g(j);
                    for k=1:p
                        mue=g(k);
                        [~, Be, Jdet] = computeShapeFunction(porous3d8n, zeta, eta, mue);
                        stiffnessMatrix=stiffnessMatrix+(w(i)*w(j)*w(k)*Jdet*transpose([Be; Be])*(Emat*[Be; Be]));
                    end
                end
            end
        end
          
        function massMatrix = computeLocalMassMatrix(porous3d8n)
            roh = porous3d8n.getPropertyValue('DENSITY');
            p = porous3d8n.getPropertyValue('NUMBER_GAUSS_POINT');
            massMatrix=zeros(24,24);
            [w,g]=returnGaussPoint(p);
            
%             for i=1:p
%                 zeta=g(i);
%                 for j=1:p
%                     eta=g(j);
%                     for k=1:p
%                         mue=g(k);
%                         [N, ~, Jdet] = computeShapeFunction(porous3d8n, zeta, eta, mue);
%                         massMatrix=massMatrix+(w(i)*w(j)*w(k)*roh*transpose(N)*N*Jdet);
%                     end
%                 end
%             end
%             
%         end
        
        function dofs = getDofList(element)
            dofs([1 7 13 19 25 31 37 43]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_X');
            dofs([2 8 14 20 26 32 38 44]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y');
            dofs([3 9 15 21 27 33 39 45]) = element.nodeArray.getDof('DISPLACEMENT_SOLID_Z');
            dofs([4 10 16 22 28 34 40 46]) = element.nodeArray.getDof('DISPLACEMENT_FLUID_X');
            dofs([5 11 17 23 29 35 41 47]) = element.nodeArray.getDof('DISPLACEMENT_FLUID_Y');
            dofs([6 12 18 24 30 36 42 48]) = element.nodeArray.getDof('DISPLACEMENT_FLUID_Z');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,48);
            
            vals([1 7 13 19 25 31 37 43]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_X',step);
            vals([2 8 14 20 26 32 38 44]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_Y',step);
            vals([3 9 15 21 27 33 39 45]) = element.nodeArray.getDofValue('DISPLACEMENT_SOLID_Z',step);
            vals([4 10 16 22 28 34 40 46]) = element.nodeArray.getDofValue('DISPLACEMENT_FLUID_X',step);
            vals([5 11 17 23 29 35 41 47]) = element.nodeArray.getDofValue('DISPLACEMENT_FLUID_Y',step);
            vals([6 12 18 24 30 36 42 48]) = element.nodeArray.getDofValue('DISPLACEMENT_FLUID_Z',step);
        end
        
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,48);
            
            [~, vals([1 7 13 19 25 31 37 43]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, vals([2 8 14 20 26 32 38 44]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, vals([3 9 15 21 27 33 39 45]), ~] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Z').getAllValues(step);
            [~, vals([4 10 16 22 28 34 40 46]), ~] = element.nodeArray.getDof('DISPLACEMENT_FLUID_X').getAllValues(step);
            [~, vals([5 11 17 23 29 35 41 47]), ~] = element.nodeArray.getDof('DISPLACEMENT_FLUID_Y').getAllValues(step);
            [~, vals([6 12 18 24 30 36 42 48]), ~] = element.nodeArray.getDof('DISPLACEMENT_FLUID_Z').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,48);
            
            [~, ~, vals([1 7 13 19 25 31 37 43])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_X').getAllValues(step);
            [~, ~, vals([2 8 14 20 26 32 38 44])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Y').getAllValues(step);
            [~, ~, vals([3 9 15 21 27 33 39 45])] = element.nodeArray.getDof('DISPLACEMENT_SOLID_Z').getAllValues(step);
            [~, ~, vals([4 10 16 22 28 34 40 46])] = element.nodeArray.getDof('DISPLACEMENT_FLUID_X').getAllValues(step);
            [~, ~, vals([5 11 17 23 29 35 41 47])] = element.nodeArray.getDof('DISPLACEMENT_FLUID_Y').getAllValues(step);
            [~, ~, vals([6 12 18 24 30 36 42 48])] = element.nodeArray.getDof('DISPLACEMENT_FLUID_Z').getAllValues(step);
        end
        
        function pl = draw(obj)
            x = [obj.nodeArray(1).getX, obj.nodeArray(2).getX, ...
                obj.nodeArray(3).getX, obj.nodeArray(4).getX,...
                obj.nodeArray(1).getX, obj.nodeArray(5).getX, ...
                obj.nodeArray(6).getX, obj.nodeArray(7).getX,...
                obj.nodeArray(8).getX, obj.nodeArray(5).getX,...
                obj.nodeArray(6).getX, obj.nodeArray(2).getX,...
                obj.nodeArray(3).getX, obj.nodeArray(7).getX,...
                obj.nodeArray(8).getX, obj.nodeArray(4).getX];
            
            y = [obj.nodeArray(1).getY, obj.nodeArray(2).getY, ...
                obj.nodeArray(3).getY, obj.nodeArray(4).getY,...
                obj.nodeArray(1).getY, obj.nodeArray(5).getY, ...
                obj.nodeArray(6).getY, obj.nodeArray(7).getY,...
                obj.nodeArray(8).getY, obj.nodeArray(5).getY,...
                obj.nodeArray(6).getY, obj.nodeArray(2).getY,...
                obj.nodeArray(3).getY, obj.nodeArray(7).getY,...
                obj.nodeArray(8).getY, obj.nodeArray(4).getY];
            
            z = [obj.nodeArray(1).getZ, obj.nodeArray(2).getZ, ...
                obj.nodeArray(3).getZ, obj.nodeArray(4).getZ,...
                obj.nodeArray(1).getZ, obj.nodeArray(5).getZ, ...
                obj.nodeArray(6).getZ, obj.nodeArray(7).getZ,...
                obj.nodeArray(8).getZ, obj.nodeArray(5).getZ,...
                obj.nodeArray(6).getZ, obj.nodeArray(2).getZ,...
                obj.nodeArray(3).getZ, obj.nodeArray(7).getZ,...
                obj.nodeArray(8).getZ, obj.nodeArray(4).getZ];
            
                pl = line(x,y,z);
            
        end
        
        function pl = drawDeformed(obj, step, scaling)
            % ONLY DRAWS DEFORMATION OF SOLID
            x = [obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_SOLID_X', step), ...
                obj.nodeArray(2).getX + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_SOLID_X', step), ...
                obj.nodeArray(3).getX + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_SOLID_X', step), ...
                obj.nodeArray(4).getX + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_SOLID_X', step),...
                obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_SOLID_X', step), ...
                obj.nodeArray(5).getX + scaling * obj.nodeArray(5).getDofValue('DISPLACEMENT_SOLID_X', step), ...
                obj.nodeArray(6).getX + scaling * obj.nodeArray(6).getDofValue('DISPLACEMENT_SOLID_X', step), ...
                obj.nodeArray(7).getX + scaling * obj.nodeArray(7).getDofValue('DISPLACEMENT_SOLID_X', step), ...
                obj.nodeArray(8).getX + scaling * obj.nodeArray(8).getDofValue('DISPLACEMENT_SOLID_X', step),...
                obj.nodeArray(5).getX + scaling * obj.nodeArray(5).getDofValue('DISPLACEMENT_SOLID_X', step), ...
                obj.nodeArray(6).getX + scaling * obj.nodeArray(6).getDofValue('DISPLACEMENT_SOLID_X', step), ...
                obj.nodeArray(2).getX + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_SOLID_X', step), ...
                obj.nodeArray(3).getX + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_SOLID_X', step), ...
                obj.nodeArray(7).getX + scaling * obj.nodeArray(7).getDofValue('DISPLACEMENT_SOLID_X', step),...
                obj.nodeArray(8).getX + scaling * obj.nodeArray(8).getDofValue('DISPLACEMENT_SOLID_X', step), ...
                obj.nodeArray(4).getX + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_SOLID_X', step)];
            
            y = [obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_SOLID_Y', step), ...
                obj.nodeArray(2).getY + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_SOLID_Y', step), ...
                obj.nodeArray(3).getY + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_SOLID_Y', step), ...
                obj.nodeArray(4).getY + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_SOLID_Y', step),...
                obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_SOLID_Y', step), ...
                obj.nodeArray(5).getY + scaling * obj.nodeArray(5).getDofValue('DISPLACEMENT_SOLID_Y', step), ...
                obj.nodeArray(6).getY + scaling * obj.nodeArray(6).getDofValue('DISPLACEMENT_SOLID_Y', step), ...
                obj.nodeArray(7).getY + scaling * obj.nodeArray(7).getDofValue('DISPLACEMENT_SOLID_Y', step),...
                obj.nodeArray(8).getY + scaling * obj.nodeArray(8).getDofValue('DISPLACEMENT_SOLID_Y', step),...
                obj.nodeArray(5).getY + scaling * obj.nodeArray(5).getDofValue('DISPLACEMENT_SOLID_Y', step), ...
                obj.nodeArray(6).getY + scaling * obj.nodeArray(6).getDofValue('DISPLACEMENT_SOLID_Y', step), ...
                obj.nodeArray(2).getY + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_SOLID_Y', step), ...
                obj.nodeArray(3).getY + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_SOLID_Y', step), ...
                obj.nodeArray(7).getY + scaling * obj.nodeArray(7).getDofValue('DISPLACEMENT_SOLID_Y', step),...
                obj.nodeArray(8).getY + scaling * obj.nodeArray(8).getDofValue('DISPLACEMENT_SOLID_Y', step), ...
                obj.nodeArray(4).getY + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_SOLID_Y', step)];
            
            z = [obj.nodeArray(1).getZ + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_SOLID_Z', step), ...
                obj.nodeArray(2).getZ + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_SOLID_Z', step), ...
                obj.nodeArray(3).getZ + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_SOLID_Z', step), ...
                obj.nodeArray(4).getZ + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_SOLID_Z', step),...
                obj.nodeArray(1).getZ + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_SOLID_Z', step), ...
                obj.nodeArray(5).getZ + scaling * obj.nodeArray(5).getDofValue('DISPLACEMENT_SOLID_Z', step), ...
                obj.nodeArray(6).getZ + scaling * obj.nodeArray(6).getDofValue('DISPLACEMENT_SOLID_Z', step), ...
                obj.nodeArray(7).getZ + scaling * obj.nodeArray(7).getDofValue('DISPLACEMENT_SOLID_Z', step),...
                obj.nodeArray(8).getZ + scaling * obj.nodeArray(8).getDofValue('DISPLACEMENT_SOLID_Z', step),...
                obj.nodeArray(5).getZ + scaling * obj.nodeArray(5).getDofValue('DISPLACEMENT_SOLID_Z', step), ...
                obj.nodeArray(6).getZ + scaling * obj.nodeArray(6).getDofValue('DISPLACEMENT_SOLID_Z', step), ...
                obj.nodeArray(2).getZ + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_SOLID_Z', step), ...
                obj.nodeArray(3).getZ + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_SOLID_Z', step), ...
                obj.nodeArray(7).getZ + scaling * obj.nodeArray(7).getDofValue('DISPLACEMENT_SOLID_Z', step),...
                obj.nodeArray(8).getZ + scaling * obj.nodeArray(8).getDofValue('DISPLACEMENT_SOLID_Z', step), ...
                obj.nodeArray(4).getZ + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_SOLID_Z', step)];
                
                pl = line(x,y,z);
                
        end
        
        function update(porous3d8n)
        end
        
        function F = computeLocalForceVector(porous3d8n)
            F = zeros(1,48);
        end

        
    end
      
end

