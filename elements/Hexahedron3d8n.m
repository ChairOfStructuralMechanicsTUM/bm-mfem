classdef Hexahedron3d8n < Hexahedron %Class Hexahedron to be implemented
    %HEXAHEDRON3D8N Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
    end
    
    methods
        % constructor
        function hexahedron3d8n = Hexahedron3d8n(id, nodeArray, material)
            
            % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 3
                super_args = {id; material};
            end
            
            % call the super class constructor
            hexahedron3d8n@Hexahedron(super_args{:});
            
            hexahedron3d8n.dofNames = cellstr(['DISPLACEMENT_X'; 'DISPLACEMENT_Y'; 'DISPLACEMENT_Z']);
            hexahedron3d8n.requiredProperties = [ ];
            hexahedron3d8n.required3dProperties = [ ];
            
            % the constructor
            if nargin > 0
                if (length(nodeArray) == 8 && isa(nodeArray,'Node'))
                    hexahedron3d8n.nodeArray = nodeArray;
                else
                    error('problem with the nodes in element %d', id);
                end
                
                %In comparison to truss no length-calculation; volume
                %calculation might be necessary
                
            end
            
        end
        
        % getter functions
        
        function responseDoF = getResponseDofArray(hexahedron, step)
            
            responseDoF = zeros(24,1);
            
            for itNodes = 1:1:2
                nodalDof = hexahedron.nodeArray(itNodes).getDofArray;
                nodalDof = nodalDof.';
                
                for itDof = 3:(-1):1
                    responseDoF(3*itNodes-(itDof-1),1) = nodalDof(4-itDof).getValue(step);
                end
            end
            
        end
        
        

        % member functions
        
        % Calculation of shape functions
        % Nf = 1/8 * [(1-zeta)*(1-eta)*(1-mue);(1+zeta)*(1-eta)*(1-mue);(1+ zeta)*(1+eta)*(1-mue);(1-zeta)*(1+eta)*(1-mue);(1- zeta)*(1-eta)*(1+mue);(1+
        % zeta)*(1-eta)*(1+mue);(1+ zeta)*(1+eta)*(1+mue);(1-zeta)*(1+eta)*(1+mue)]
        
        % Calculation of Jacobian matrix
        % [dNzeta; dNeta; dNmue]= [(-1)*(1-eta)*(1-mue),(1-eta)*(1-mue),(1+eta)*(1-mue),(-1)*(1+eta)*(1-mue),(-1)*(1-eta)*(1+mue),(1-eta)*(1+mue),(1+eta)*(1+mue),(-1)*(1+eta)*(1+mue); 
        % (1- zeta)*(-1)*(1-mue),(1+zeta)*(-1)*(1-mue),(1+ zeta)*(1-mue),(1-zeta)*(1-mue),(1- zeta)*(-1)*(1+mue),(1+zeta)*(-1)*(1+mue),(1+ zeta)*(1+mue),(1- zeta)*(1+mue); 
        % (1-zeta)*(1-eta)*(-1),(1+zeta)*(1-eta)*(-1),(1+zeta)*(1+eta)*(-1),(1-zeta)*(1+eta)*(-1),(1- zeta)*(1-eta),(1+zeta)*(1-eta),(1+ zeta)*(1+eta),(1-zeta)*(1+eta)]
        
        % Numerical integration and calculation of mass and stiffnessmatrix
    
    end
    
end

