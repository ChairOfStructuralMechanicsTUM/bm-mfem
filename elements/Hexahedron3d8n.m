classdef Hexahedron3d8n < Element %Class Hexahedron to be implemented
    %HEXAHEDRON3D8N Summary of this class goes here
    %   Detailed explanation goes here
    
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
            hexahedron3d8n.requiredProperties = cellstr(["YOUNGS_MODULUS", "POISSON_RATIO"]);
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
        
        function [Nf, Bx, By, Bz, Jdet] = computeShapeFunction(hexahedron3d8n)
            
            % Calculation of shape functions
            
            Nf = 1/8 * [(1-zeta)*(1-eta)*(1-mue);(1+zeta)*(1-eta)*(1-mue);(1+ zeta)*(1+eta)*(1-mue);(1-zeta)*(1+eta)*(1-mue);(1- zeta)*(1-eta)*(1+mue);(1+zeta)*(1-eta)*(1+mue);(1+ zeta)*(1+eta)*(1+mue);(1-zeta)*(1+eta)*(1+mue)];
            
            % member functions
            
            % Calculation of shape function derivatives
            
            dNzeta = [(-1)*(1-eta)*(1-mue),(1-eta)*(1-mue),(1+eta)*(1-mue),(-1)*(1+eta)*(1-mue),(-1)*(1-eta)*(1+mue),(1-eta)*(1+mue),(1+eta)*(1+mue),(-1)*(1+eta)*(1+mue)];
            dNeta = [(1- zeta)*(-1)*(1-mue),(1+zeta)*(-1)*(1-mue),(1+ zeta)*(1-mue),(1-zeta)*(1-mue),(1- zeta)*(-1)*(1+mue),(1+zeta)*(-1)*(1+mue),(1+ zeta)*(1+mue),(1- zeta)*(1+mue)];
            dNmue = [(1-zeta)*(1-eta)*(-1),(1+zeta)*(1-eta)*(-1),(1+zeta)*(1+eta)*(-1),(1-zeta)*(1+eta)*(-1),(1- zeta)*(1-eta),(1+zeta)*(1-eta),(1+ zeta)*(1+eta),(1-zeta)*(1+eta)];
            
            % Calculation of Jacobian matrix
            
            coords=zeros(8,3);
            for i=1:1:8
                coords(i,1:3)=hexahedron3d8n.Element.nodeArray(i).getCoords;
            end
            
            xn=coords(i,1);
            yn=coords(i,2);
            zn=coords(i,3);
            
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
            
            Jdet=J11*J22*J33+J21*J32*J31+J31*J12*J23-J31*J22*J13-J11*J32*J23-J21*J12*J33;
            
            Jinv=[J22*J33-J32*J23, J32*J13-J12*J33, J12*J23-J22*J13;
                J31*J23-J21*J33, J11*J33-J31*J13, J21*J13-J11*J23;
                J21*J32-J31*J22, J31*J12-J11*J32, J11*J22-J21*J12 ];
            
            % Calculation of B-Matrix
            Bx = Jinv*dNzeta/Jdet;
            By = Jinv*dNeta/Jdet;
            Bz = Jinv*dNmue/Jdet;
            
        end
        
        
        
        function stiffnessMatrix = computeLocalStiffnessMatrix(hexahedron3d8n)
            Emodul = hexahedron3d8n.getPropertyValue('YOUNGS_MODULUS');
        end
    end
    
    methods (Access = protected)
        
        function cp = copyElement(obj)
            cp = copyElement@Element(obj);
            % Extension might be necessary at this point
        end
        
    end
    
    
end

