classdef Material < handle
    %MATERIAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        name
        parameters
    end
    
    methods
        % constructor
        function material = Material(name)
            material.name = name;
            material.parameters = containers.Map; 
        end
        
        % getter functions
        function name = getName(material)
            name = material.name;
        end
        
        function parameters = getParameters(material)
            parameters = material.parameters;
        end
        
        % member functions
        function addParameter(material, name, value)
            material.parameters(name) = value;
        end
        
        function value = getParameterValue(material, name)
            value = material.parameters(name);
        end
        
    end
    
end

