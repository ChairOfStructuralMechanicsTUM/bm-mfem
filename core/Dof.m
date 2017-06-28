classdef Dof < handle
    %DOF The degree of freedom class
    %   Detailed explanation goes here
    
    properties (Access = private)
        node
        value
        valueType
        dofLoad
        fixed = false
    end
    
    methods
        %constructor
        function dof = Dof(node, value, valueType)
            if nargin > 0
                if (isa(node,'Node'))
                    dof.node = node;
                else
                    error('invalid node')
                end
                dof.value = value;
                dof.valueType = valueType;
            end
        end
        
        % getter functions
        function node = getNode(dof)
           node = dof.node; 
        end
        
        function value = getValue(dofs)
            value = zeros;
            for ii = 1:length(dofs)
                value(ii) = dofs(ii).value;
            end
        end
        
        function valueType = getValueType(dof)
            valueType = dof.valueType;
        end
        
        function dofLoad = getDofLoad(dof)
            if isempty(dof.dofLoad)
                dofLoad = 0;
            else 
                dofLoad = dof.dofLoad;
            end
         end
        
        function fixed = isFixed(dof)
            fixed = dof.fixed;
        end
        
        % setter functions
        function fix(dof)
            dof.fixed = true;
        end
        
        function unfix(dof)
            dof.fixed = false;
        end
        
        function setValue(dof, value)
           dof.value = value; 
        end
        
        function setLoad(dof, dofLoad)
            dof.dofLoad = dofLoad;
        end
        
    end
    
end

