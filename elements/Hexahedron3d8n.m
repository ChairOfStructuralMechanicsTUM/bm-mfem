classdef Hexahedron3d8n < Element  %Class Hexahedron to be implemented
    %HEXAHEDRON3D8N Summary of this class goes here
    %   Detailed explanation goes here
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Anpassen von Basisklasse   %
    %  Auslagern von Emat         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Access = private)
    end
    
    methods
        % constructor
        function hexahedron3d8n = Hexahedron3d8n(id, nodeArray)
            
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 2
                super_args = {id};
            end
            
            % call the super class constructor
            hexahedron3d8n@Element(super_args{:});
            
            hexahedron3d8n.dofNames = cellstr(['DISPLACEMENT_X'; 'DISPLACEMENT_Y'; 'DISPLACEMENT_Z']);
            hexahedron3d8n.requiredProperties = cellstr(["YOUNGS_MODULUS", "POISSON_RATIO","NUMBER_GAUSS_POINT"]);
            hexahedron3d8n.required3dProperties = [ ];
            
            % the constructor
            if nargin > 0
                if (length(nodeArray) == 8 && isa(nodeArray,'Node'))
                    hexahedron3d8n.nodeArray = nodeArray;
                else
                    error('problem with the nodes in element %d', id);
                end
                
            end
            
        end
        
        % getter functions
        
        function responseDoF = getResponseDofArray(hexahedron, step)
            
            responseDoF = zeros(24,1);
            
            for itNodes = 1:1:8
                nodalDof = hexahedron.nodeArray(itNodes).getDofArray;
                nodalDof = nodalDof.';
                
                for itDof = 3:(-1):1
                    responseDoF(3*itNodes-(itDof-1),1) = nodalDof(4-itDof).getValue(step);
                end
            end
            
        end
        
        % member functions
        
        function [Bx, By, Bz, Jdet] = computeShapeFunction(hexahedron3d8n, zeta, eta, mue)
            
            % Calculation of shape functions
            
            %Nf = 1/8 * [(1-zeta)*(1-eta)*(1-mue);(1+zeta)*(1-eta)*(1-mue);(1+ zeta)*(1+eta)*(1-mue);(1-zeta)*(1+eta)*(1-mue);(1- zeta)*(1-eta)*(1+mue);(1+zeta)*(1-eta)*(1+mue);(1+ zeta)*(1+eta)*(1+mue);(1-zeta)*(1+eta)*(1+mue)];
            
            % member functions
            
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
            % J11.*J22.*J33+J21.*J32.*J31+J31.*J12.*J23-J31.*J22.*J13-J11.*J32.*J23-J21.*J12.*J33;
            
            Jinv=[J22*J33-J32*J23, J32*J13-J12*J33, J12*J23-J22*J13;
                J31*J23-J21*J33, J11*J33-J31*J13, J21*J13-J11*J23;
                J21*J32-J31*J22, J31*J12-J11*J32, J11*J22-J21*J12 ];
            
            % Calculation of B-Matrix
            B=Jinv*[dNzeta;dNeta;dNmue]/Jdet;
            Bx=B(1,1:8);
            By=B(2,1:8);
            Bz=B(3,1:8);
            
        end
        
        
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(hexahedron3d8n)
            Emodul = hexahedron3d8n.getPropertyValue('YOUNGS_MODULUS');
            PoissonRatio = hexahedron3d8n.getPropertyValue('POISSON_RATIO');
            p = hexahedron3d8n.getPropertyValue('NUMBER_GAUSS_POINT');
            
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
                        
                        [Bx, By, Bz, Jdet] = computeShapeFunction(hexahedron3d8n, zeta, eta, mue);
                        %Be=[Bx, zeros(1,16); zeros(1,8), By, zeros(1,8); zeros(1,16), Bz; By, Bx, zeros(1,8); zeros(1,8), Bz, By; Bz, zeros(1,8), Bx];
                        Be=[Bx(1),0,0,Bx(2),0,0,Bx(3),0,0,Bx(4),0,0,Bx(5),0,0,Bx(6),0,0,Bx(7),0,0,Bx(8),0,0;
                            0,By(1),0,0,By(2),0,0,By(3),0,0,By(4),0,0,By(5),0,0,By(6),0,0,By(7),0,0,By(8),0;
                            0,0,Bz(1),0,0,Bz(2),0,0,Bz(3),0,0,Bz(4),0,0,Bz(5),0,0,Bz(6),0,0,Bz(7),0,0,Bz(8);
                            By(1),Bx(1),0,By(2),Bx(2),0,By(3),Bx(3),0,By(4),Bx(4),0,By(5),Bx(5),0,By(6),Bx(6),0,By(7),Bx(7),0,By(8),Bx(8),0;
                            0,Bz(1),By(1),0,Bz(2),By(2),0,Bz(3),By(3),0,Bz(4),By(4),0,Bz(5),By(5),0,Bz(6),By(6),0,Bz(7),By(7),0,Bz(8),By(8);
                            Bz(1),0,Bx(1),Bz(2),0,Bx(2),Bz(3),0,Bx(3),Bz(4),0,Bx(4),Bz(5),0,Bx(5),Bz(6),0,Bx(6),Bz(7),0,Bx(7),Bz(8),0,Bx(8)];
                        
                        stiffnessMatrix=stiffnessMatrix+(w(i)*w(j)*w(k)*Jdet*transpose(Be)*(Emat*Be));
                        
                    end
                end
            end
        end
        
        
    function computeLocalForceVector(e)
        
    end
    
    function barycenter(e)
        
    end
    
    function update(e)
        
    end
        
        
    end
    
    
    methods (Access = protected)
        
        function cp = copyElement(obj)
            cp = copyElement@Element(obj);
            % Extension might be necessary at this point
        end
        
    end
    
    
end

