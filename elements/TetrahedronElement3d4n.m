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
        function [NfDiff, Jdet] = computeShapeFunction(tetrahedron3d4n)
            
%             % Calculation of shape functions
%             Nf = [theta_1; theta_2; theta_3; theta_4];
%             N=zeros(3,3*4);
%             
%             for l=1:4
%                 N(1,3*(l-1)+1)=Nf(l);
%                 N(2,3*(l-1)+2)=Nf(l);
%                 N(3,3*(l-1)+3)=Nf(l);
%             end
            
            % Calculation of Determinant of Jacobian Matrix
            xDiff = zeros(4,4);
            yDiff = zeros(4,4);
            zDiff = zeros(4,4);
            
            for i = 1:4
                for j = 1:4
                    xDiff(i,j) = tetrahedron3d4n.nodeArray(i).getX - tetrahedron3d4n.nodeArray(j).getX;
                    yDiff(i,j) = tetrahedron3d4n.nodeArray(i).getY - tetrahedron3d4n.nodeArray(j).getY;
                    zDiff(i,j) = tetrahedron3d4n.nodeArray(i).getZ - tetrahedron3d4n.nodeArray(j).getZ;
                end
            end
                
%             J = [1 1 1 1;... 
%                  tetrahedron3d4n.nodeArray(1).getX   tetrahedron3d4n.nodeArray(2).getX   tetrahedron3d4n.nodeArray(3).getX   tetrahedron3d4n.nodeArray(4).getX;...  
%                  tetrahedron3d4n.nodeArray(1).getY   tetrahedron3d4n.nodeArray(2).getY   tetrahedron3d4n.nodeArray(3).getY   tetrahedron3d4n.nodeArray(4).getY;... 
%                  tetrahedron3d4n.nodeArray(1).getZ   tetrahedron3d4n.nodeArray(2).getZ   tetrahedron3d4n.nodeArray(3).getZ   tetrahedron3d4n.nodeArray(4).getZ];
%              
%              Jdet = det(J);
             Jdet = xDiff(2,1)*(yDiff(2,3)*zDiff(3,4)-yDiff(3,4)*zDiff(2,3))+xDiff(3,2)*(yDiff(3,4)*zDiff(1,2)-yDiff(1,2)*zDiff(3,4))+xDiff(4,3)*(yDiff(1,2)*zDiff(2,3)-yDiff(2,3)*zDiff(1,2));
            

            % Calculation of the derivatives of the Shape Functions
            NfDiffX = [yDiff(4,2)*zDiff(3,2)-yDiff(3,2)*zDiff(4,2)   yDiff(3,1)*zDiff(4,3)-yDiff(3,4)*zDiff(1,3)    yDiff(2,4)*zDiff(1,4)-yDiff(1,4)*zDiff(2,4)   yDiff(1,3)*zDiff(2,1)-yDiff(1,2)*zDiff(3,1)];  % Derivatives w.r.t. x
            NfDiffY = [xDiff(3,2)*zDiff(4,2)-xDiff(4,2)*zDiff(3,2)   xDiff(4,3)*zDiff(3,1)-xDiff(1,3)*zDiff(3,4)    xDiff(1,4)*zDiff(2,4)-xDiff(2,4)*zDiff(1,4)   xDiff(2,1)*zDiff(1,3)-xDiff(3,1)*zDiff(1,2)];  % Derivatives w.r.t. y
            NfDiffZ = [xDiff(4,2)*yDiff(3,2)-xDiff(3,2)*yDiff(4,2)   xDiff(3,1)*yDiff(4,3)-xDiff(3,4)*yDiff(1,3)    xDiff(2,4)*yDiff(1,4)-xDiff(1,4)*yDiff(2,4)   xDiff(1,3)*yDiff(2,1)-xDiff(1,2)*yDiff(3,1)];  % Derivatives w.r.t. z
            
            NfDiff = [NfDiffX; NfDiffY; NfDiffZ];    
            
            
           % Calculation of B-Matrix
% % %             Jinv=[J22*J33-J32*J23, J32*J13-J12*J33, J12*J23-J22*J13;
% % %             J31*J23-J21*J33, J11*J33-J31*J13, J21*J13-J11*J23;
% % %             J21*J32-J31*J22, J31*J12-J11*J32, J11*J22-J21*J12 ];
% % %             
% % %             B=Jinv*[dNzeta;dNeta;dNmue]/Jdet;
% % %             Bx=B(1,1:8);
% % %             By=B(2,1:8);
% % %             Bz=B(3,1:8);
% % %             
% % %             Be=[Bx(1),0,0,Bx(2),0,0,Bx(3),0,0,Bx(4),0,0,Bx(5),0,0,Bx(6),0,0,Bx(7),0,0,Bx(8),0,0;
% % %                 0,By(1),0,0,By(2),0,0,By(3),0,0,By(4),0,0,By(5),0,0,By(6),0,0,By(7),0,0,By(8),0;
% % %                 0,0,Bz(1),0,0,Bz(2),0,0,Bz(3),0,0,Bz(4),0,0,Bz(5),0,0,Bz(6),0,0,Bz(7),0,0,Bz(8);
% % %                 By(1),Bx(1),0,By(2),Bx(2),0,By(3),Bx(3),0,By(4),Bx(4),0,By(5),Bx(5),0,By(6),Bx(6),0,By(7),Bx(7),0,By(8),Bx(8),0;
% % %                 0,Bz(1),By(1),0,Bz(2),By(2),0,Bz(3),By(3),0,Bz(4),By(4),0,Bz(5),By(5),0,Bz(6),By(6),0,Bz(7),By(7),0,Bz(8),By(8);
% % %                 Bz(1),0,Bx(1),Bz(2),0,Bx(2),Bz(3),0,Bx(3),Bz(4),0,Bx(4),Bz(5),0,Bx(5),Bz(6),0,Bx(6),Bz(7),0,Bx(7),Bz(8),0,Bx(8)];
            
        end
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(tetrahedron3d4n)
            
            [NfDiff, Jdet] = computeShapeFunction(tetrahedron3d4n);
            
            % Compute B-Matrix
            Be = [NfDiff(1,1),0,0,            NfDiff(1,2),0,0,            NfDiff(1,3),0,0,            NfDiff(1,4),0,0;...
                  0,NfDiff(2,1),0,            0,NfDiff(2,2),0,            0,NfDiff(2,3),0,            0,NfDiff(2,4),0;...
                  0,0,NfDiff(3,1),            0,0,NfDiff(3,2),            0,0,NfDiff(3,3),            0,0,NfDiff(3,4);...
                  NfDiff(2,1),NfDiff(1,1),0,  NfDiff(2,2),NfDiff(1,2),0,  NfDiff(2,3),NfDiff(1,3),0,  NfDiff(2,4),NfDiff(1,4),0;...
                  0,NfDiff(3,1),NfDiff(2,1),  0,NfDiff(3,2),NfDiff(2,2),  0,NfDiff(3,3),NfDiff(2,3),  0,NfDiff(3,4),NfDiff(2,4);...
                  NfDiff(3,1),0,NfDiff(1,1),   NfDiff(3,2),0,NfDiff(1,2), NfDiff(3,3),0,NfDiff(1,3),  NfDiff(3,4),0,NfDiff(1,4)];
  
            % Compute Emat-Matrix (for isotropic material only)
            Emodul = tetrahedron3d4n.getPropertyValue('YOUNGS_MODULUS');
            PoissonRatio = tetrahedron3d4n.getPropertyValue('POISSON_RATIO');
            
            a=Emodul*(1-PoissonRatio)/( (1-2*PoissonRatio)*(1+PoissonRatio) );
            b=Emodul*PoissonRatio/( (1-2*PoissonRatio)*(1+PoissonRatio) );
            c=Emodul/(2*(1+PoissonRatio));
            
            Emat=[a,b,b,0,0,0;b,a,b,0,0,0;b,b,a,0,0,0;0,0,0,c,0,0;0,0,0,0,c,0;0,0,0,0,0,c];
                   
              
            % Compute Element Stiffness Matrix
            stiffnessMatrix=zeros(12,12);  
            stiffnessMatrix=(1/(6*Jdet))*(transpose(Be)*(Emat*Be));
             
           
% % %             [w,g]=returnGaussPoint(class(tetrahedron3d4n), p);
% % %             for i=1:p
% % %                 theta_1=g(i);
% % %                 for j=1:p
% % %                     theta_2=g(j);
% % %                     for k=1:p
% % %                         theta_3=g(k);
% % %                         for l=1:p
% % %                             theta_4=g(l);
% % %                                 [~, Be, Jdet] = computeShapeFunction(tetrahedron3d4n, theta_1, theta_2, theta_3, theta_4);
% % %                                 stiffnessMatrix=stiffnessMatrix+(w(i)*w(j)*w(k)*w(l)*Jdet*transpose(Be)*(Emat*Be));
% % %                         end
% % %                     end
% % %                 end
% % %             end
        end
          
        
        function massMatrix = computeLocalMassMatrix(tetrahedron3d4n)
           
        roh = tetrahedron3d4n.getPropertyValue('DENSITY');
        massMatrix = zeros(12,12);

         Jdet = xDiff(2,1)*(yDiff(2,3)*zDiff(3,4)-yDiff(3,4)*zDiff(2,3))+xDiff(3,2)*(yDiff(3,4)*zDiff(1,2)-yDiff(1,2)*zDiff(3,4))+xDiff(4,3)*(yDiff(1,2)*zDiff(2,3)-yDiff(2,3)*zDiff(1,2));
         V = 1/6*Jdet;
         massMatrix = roh*V/20 * [2 0 0 1 0 0 1 0 0 1 0 0;...
                                  0 2 0 0 1 0 0 1 0 0 1 0;...
                                  0 0 2 0 0 1 0 0 1 0 0 1;...
                                  1 0 0 2 0 0 1 0 0 1 0 0;...
                                  0 1 0 0 2 0 0 1 0 0 1 0;...
                                  0 0 1 0 0 2 0 0 1 0 0 1;...
                                  1 0 0 1 0 0 2 0 0 1 0 0;...
                                  0 1 0 0 1 0 0 2 0 0 1 0;...
                                  0 0 1 0 0 1 0 0 2 0 0 1;...
                                  1 0 0 1 0 0 1 0 0 2 0 0;...
                                  0 1 0 0 1 0 0 1 0 0 2 0;...
                                  0 0 1 0 0 1 0 0 1 0 0 2];

% % %             p = hexahedron3d8n.getPropertyValue('NUMBER_GAUSS_POINT');
% % %             [w,g]=returnGaussPoint(p);
% % %             
% % %             for i=1:p
% % %                 zeta=g(i);
% % %                 for j=1:p
% % %                     eta=g(j);
% % %                     for k=1:p
% % %                         mue=g(k);
% % %                         [N, ~, Jdet] = computeShapeFunction(hexahedron3d8n, zeta, eta, mue);
% % %                         massMatrix=massMatrix+(w(i)*w(j)*w(k)*roh*transpose(N)*N*Jdet);
% % %                     end
% % %                 end
% % %             end
% % %             
        end
        
        function dampingMatrix = computeLocalDampingMatrix(e)
            % up to now only Rayleigh-Damping
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
            dofs([1 4 7 10]) = element.nodeArray.getDof('DISPLACEMENT_X');
            dofs([2 5 8 11]) = element.nodeArray.getDof('DISPLACEMENT_Y');
            dofs([3 6 9 12]) = element.nodeArray.getDof('DISPLACEMENT_Z');
        end
        
        function vals = getValuesVector(element, step)
            vals = zeros(1,12);
            
            vals([1 4 7 10]) = element.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals([2 5 8 11]) = element.nodeArray.getDofValue('DISPLACEMENT_Y',step);
            vals([3 6 9 12]) = element.nodeArray.getDofValue('DISPLACEMENT_Z',step);
        end
        
        function vals = getFirstDerivativesVector(element, step)
            vals = zeros(1,12);
            
            [~, vals([1 4 7 10]), ~] = element.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, vals([2 5 8 11]), ~] = element.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
            [~, vals([3 6 9 12]), ~] = element.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
        end
        
        function vals = getSecondDerivativesVector(element, step)
            vals = zeros(1,12);            
            
            [~, ~, vals([1 4 7 10])] = element.nodeArray.getDof('DISPLACEMENT_X').getAllValues(step);
            [~, ~, vals([2 5 8 11])] = element.nodeArray.getDof('DISPLACEMENT_Y').getAllValues(step);
            [~, ~, vals([3 6 9 12])] = element.nodeArray.getDof('DISPLACEMENT_Z').getAllValues(step);
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
                obj.nodeArray(1).getX + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(4).getX + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_X', step),...
                obj.nodeArray(2).getX + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(3).getX + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_X', step), ...
                obj.nodeArray(4).getX + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_X', step)];
            
            y = [obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(2).getY + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(3).getY + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(1).getY + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(4).getY + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Y', step),...
                obj.nodeArray(2).getY + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(3).getY + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Y', step), ...
                obj.nodeArray(4).getY + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Y', step)];
            
            z = [obj.nodeArray(1).getZ + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(2).getZ + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(3).getZ + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(1).getZ + scaling * obj.nodeArray(1).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(4).getZ + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Z', step),...
                obj.nodeArray(2).getZ + scaling * obj.nodeArray(2).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(3).getZ + scaling * obj.nodeArray(3).getDofValue('DISPLACEMENT_Z', step), ...
                obj.nodeArray(4).getZ + scaling * obj.nodeArray(4).getDofValue('DISPLACEMENT_Z', step)];
                
                pl = line(x,y,z);
                
        end
        
        function update(tetrahedron3d4n)
        end
        
        function F = computeLocalForceVector(tetrahedron3d4n)
            F = zeros(1,12);
        end

        
    end
      
end

