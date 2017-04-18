classdef Dof < handle
    %DOF The degree of freedom class
    %   Detailed explanation goes here
    
    properties (Access = private)
        node
        variable
        variableType
    end
    
    methods
        %constructor
        function dof = Dof(node, variable, variableType)
            if (isa(node,'node'))
                dof.node = node;
            else
                error('invalid node')
            end
            dof.variable = variable;
            dof.variableType = variableType;
        end
        
        % getter functions
        function node = getNode(dof)
           node = dof.node; 
        end
        
        function variable = getVariable(dof)
            variable = dof.variable;
        end
        
        function variableType = getVariableType(dof)
            variableType = dof.variableType;
        end
        
    end
    
end

