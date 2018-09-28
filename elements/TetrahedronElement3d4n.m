classdef TetrahedronElement3d4n < Element  %Class Tetrahedron to be implemented
    %TETRAHEDRON3D4N Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = private)
%         lengthX
%         lengthY
%         lenghtZ
    end
    
    methods
        % constructor
        function tetrahedron3d4n = TetrahedronElement3d4n(id, nodeArray)
            
            requiredPropertyNames = cellstr(["YOUNGS_MODULUS", "POISSON_RATIO","NUMBER_GAUSS_POINT", "DENSITY"]);
            
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
            tetrahedron3d4n@Element(super_args{:});
            tetrahedron3d4n.dofNames = cellstr(['DISPLACEMENT_X'; 'DISPLACEMENT_Y'; 'DISPLACEMENT_Z']);
             
        end
        
        
        % Compute Barycenter of a Tetrahedron
        function c = barycenter(tetrahedron3d4n)
            c = 1/4 * (tetrahedron3d4n.nodeArray(1).getCoords() + tetrahedron3d4n.nodeArray(2).getCoords() + tetrahedron3d4n.nodeArray(3).getCoords() + tetrahedron3d4n.nodeArray(4).getCoords());
        end
        
        
        % Check Convexity of Tetrahedron
        function checkConvexity(tetrahedron3d4n)
            try
                [~] = tetrahedron3d4n.barycenter();
            catch
                error('Element %i is not convex', tetrahedron3d4n.getId());
            end
        end
        
%         %Initialization
        function initialize(tetrahedron3d4n)
%             tetrahedron3d4n.lengthX = computeLength(tetrahedron3d4n.nodeArray(1).getCoords, ...
%                 tetrahedron3d4n.nodeArray(2).getCoords);
%             
%            tetrahedron3d4n.lengthY = computeLength(tetrahedron3d4n.nodeArray(1).getCoords, ...
%                 tetrahedron3d4n.nodeArray(3).getCoords);
%             
%            tetrahedron3d4n.lengthY = computeLength(tetrahedron3d4n.nodeArray(1).getCoords, ...
%                 tetrahedron3d4n.nodeArray(4).getCoords);
%             
%             checkConvexity(tetrahedron3d4n);
        end

        % member functions
        function [N, Be, Jdet] = computeShapeFunction(tetrahedron3d4n, theta_1, theta_2, theta_3, theta_4)
            
            % Calculation of shape functions
            Nf = [theta_1; theta_2; theta_3; theta_4];
            N=zeros(3,3*4);
            
            for l=1:4
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
                coords(i,1:3)=hexahedron3d8n.nodeArray(i).getCoords;
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
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(tetrahedron3d4n)
            Emodul = tetrahedron3d4n.getPropertyValue('YOUNGS_MODULUS');
            PoissonRatio = tetrahedron3d4n.getPropertyValue('POISSON_RATIO');
            p = tetrahedron3d4n.getPropertyValue('NUMBER_GAUSS_POINT');
            
            % Matrix valid for isotropic material only
            a=Emodul*(1-PoissonRatio)/( (1-2*PoissonRatio)*(1+PoissonRatio) );
            b=Emodul*PoissonRatio/( (1-2*PoissonRatio)*(1+PoissonRatio) );
            c=Emodul/(2*(1+PoissonRatio));
            Emat=[a,b,b,0,0,0;b,a,b,0,0,0;b,b,a,0,0,0;0,0,0,c,0,0;0,0,0,0,c,0;0,0,0,0,0,c];
            stiffnessMatrix=zeros(24,24);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                zeta=g(i);
                for j=1:p
                    eta=g(j);
                    for k=1:p
                        mue=g(k);
                        [~, Be, Jdet] = computeShapeFunction(hexahedron3d8n, zeta, eta, mue);
                        stiffnessMatrix=stiffnessMatrix+(w(i)*w(j)*w(k)*Jdet*transpose(Be)*(Emat*Be));
                    end
                end
            end
        end
          
        function massMatrix = computeLocalMassMatrix(tetrahedron3d4n)
            roh = tetrahedron3d4n.getPropertyValue('DENSITY');
            p = hexahedron3d8n.getPropertyValue('NUMBER_GAUSS_POINT');
            massMatrix=zeros(24,24);
            [w,g]=returnGaussPoint(p);
            
            for i=1:p
                zeta=g(i);
                for j=1:p
                    eta=g(j);
                    for k=1:p
                        mue=g(k);
                        [N, ~, Jdet] = computeShapeFunction(hexahedron3d8n, zeta, eta, mue);
                        massMatrix=massMatrix+(w(i)*w(j)*w(k)*roh*transpose(N)*N*Jdet);
                    end
                end
            end
            
        end
        
        function dampingMatrix = computeLocalDampingMatrix(e)
            eProperties = e.getProperties;
            dampingMatrix = sparse(12,12);
            
            if (eProperties.hasValue('RAYLEIGH_ALPHA'))
                alpha = eProperties.getValue('RAYLEIGH_ALPHA');
                dampingMatrix = dampingMatrix + alpha * element.computeLocalMassMatrix;
            end
            
            if (eProperties.hasValue('RAYLEIGH_BETA'))
                beta = eProperties.getValue('RAYLEIGH_BETA');
                dampingMatrix = dampingMatrix + beta * element.computeLocalStiffnessMatrix;
            end
            
        end
        
        function dofs = getDofList(element)
            dofs([1 4 7 10 13 16 19 22]) = element.nodeArray.getDof('DISPLACEMENT_X');
            dofs([2 5 8 11 14 17 20 23]) = element.nodeArray.getDof('DISPLACEMENT_Y');
            dofs([3 6 9 12 15 18 21 24]) = element.nodeArray.getDof('DISPLACEMENT_Z');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,24);
            
            vals([1 4 7 10 13 16 19 22]) = element.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals([2 5 8 11 14 17 20 23]) = element.nodeArray.getDofValue('DISPLACEMENT_Y',step);
            vals([3 6 9 12 15 18 21 24]) = element.nodeArray.getDofValue('DISPLACEMENT_Z',step);
        end
        
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,12);
            
            [~, vals([1 4 7 10 13 16 19 22]), ~] = element.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, vals([2 5 8 11 14 17 20 23]), ~] = element.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
            [~, vals([3 6 9 12 15 18 21 24]), ~] = element.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,12);            
            
            [~, ~, vals([1 4 7 10 13 16 19 22])] = element.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, ~, vals([2 5 8 11 14 17 20 23])] = element.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
            [~, ~, vals([3 6 9 12 15 18 21 24])] = element.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
        end
        
        function pl = draw(obj)
            x = [obj.nodeArray(1).getX, obj.nodeArray(2).getX, ...
                obj.nodeArray(3).getX, obj.nodeArray(1).getX,...
                obj.nodeArray(4).getX, obj.nodeArray(2).getX, ...
                obj.nodeArray(3).getX, obj.nodeArray(4).getX];
            
            y = [obj.nodeArray(1).getY, obj.nodeArray(2).getY, ...
                obj.nodeArray(3).getY, obj.nodeArray(1).getY,...
                obj.nodeArray(4).getY, obj.nodeArray(2).getY, ...
                obj.nodeArray(3).getY, obj.nodeArray(4).getY];
            
            z = [obj.nodeArray(1).getZ, obj.nodeArray(2).getZ, ...
                obj.nodeArray(3).getZ, obj.nodeArray(1).getZ,...
                obj.nodeArray(4).getZ, obj.nodeArray(2).getZ, ...
                obj.nodeArray(3).getZ, obj.nodeArray(4).getZ];
            
                pl = line(x,y,z);
        end
        
        function pl = drawDeformed(obj, step, scaling)
            x = [obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(2).getX + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(3).getX + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(4).getX + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_X', step),...
                obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(5).getX + scaling * obj.nodeArray(5).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(6).getX + scaling * obj.nodeArray(6).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(7).getX + scaling * obj.nodeArray(7).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(8).getX + scaling * obj.nodeArray(8).getDofValue('DISPLACEMENT_X', step),...
                obj.nodeArray(5).getX + scaling * obj.nodeArray(5).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(6).getX + scaling * obj.nodeArray(6).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(2).getX + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(3).getX + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(7).getX + scaling * obj.nodeArray(7).getDofValue('DISPLACEMENT_X', step),...
                obj.nodeArray(8).getX + scaling * obj.nodeArray(8).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(4).getX + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_X', step)];
            
            y = [obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(2).getY + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(3).getY + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(4).getY + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Y', step),...
                obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(5).getY + scaling * obj.nodeArray(5).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(6).getY + scaling * obj.nodeArray(6).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(7).getY + scaling * obj.nodeArray(7).getDofValue('DISPLACEMENT_Y', step),...
                obj.nodeArray(8).getY + scaling * obj.nodeArray(8).getDofValue('DISPLACEMENT_Y', step),...
                obj.nodeArray(5).getY + scaling * obj.nodeArray(5).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(6).getY + scaling * obj.nodeArray(6).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(2).getY + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(3).getY + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(7).getY + scaling * obj.nodeArray(7).getDofValue('DISPLACEMENT_Y', step),...
                obj.nodeArray(8).getY + scaling * obj.nodeArray(8).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(4).getY + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Y', step)];
            
            z = [obj.nodeArray(1).getZ + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(2).getZ + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(3).getZ + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(4).getZ + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Z', step),...
                obj.nodeArray(1).getZ + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(5).getZ + scaling * obj.nodeArray(5).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(6).getZ + scaling * obj.nodeArray(6).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(7).getZ + scaling * obj.nodeArray(7).getDofValue('DISPLACEMENT_Z', step),...
                obj.nodeArray(8).getZ + scaling * obj.nodeArray(8).getDofValue('DISPLACEMENT_Z', step),...
                obj.nodeArray(5).getZ + scaling * obj.nodeArray(5).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(6).getZ + scaling * obj.nodeArray(6).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(2).getZ + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(3).getZ + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(7).getZ + scaling * obj.nodeArray(7).getDofValue('DISPLACEMENT_Z', step),...
                obj.nodeArray(8).getZ + scaling * obj.nodeArray(8).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(4).getZ + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Z', step)];
                
                pl = line(x,y,z);
                
        end
        
        function update(hexahedron3d8n)
        end
        
        function F = computeLocalForceVector(hexahedron3d8n)
            F = zeros(1,24);
        end

        
    end
      
end

