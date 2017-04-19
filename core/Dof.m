classdef Dof < handle
    %DOF The degree of freedom class
    %   Detailed explanation goes here
    
    properties (Access = private)
        node
        value
        valueType
        fixed = false
    end
    
    methods
        %constructor
        function dof = Dof(node, value, valueType)
            if (isa(node,'Node'))
                dof.node = node;
            else
                error('invalid node')
            end
            dof.value = value;
            dof.valueType = valueType;
        end
        
        % getter functions
        function node = getNode(dof)
           node = dof.node; 
        end
        
        function value = getValue(dof)
            value = dof.value;
        end
        
        function valueType = getValueType(dof)
            valueType = dof.valueType;
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
        
    end
    
end

