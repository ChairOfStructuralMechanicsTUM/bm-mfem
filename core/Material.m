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
            if nargin == 1
                material.name = name;
            end
            material.parameters = containers.Map; 
        end
        
        % getter functions
        function name = getName(material)
            name = material.name;
        end
        
        function parameters = getParameters(material)
            parameters = material.parameters;
        end
        
        function names = getValueNames(material)
            names = material.parameters.keys;
        end
        
        % member functions
        %         function setParameter(material, name, value)
        %             %alternative name for addValue for compatibility with the
        %             %PropertyContainer
        %             material.addValue(name, value);
        %         end
        
        function setValue(material, name, value)
            %alternative name for addValue for compatibility with the
            %PropertyContainer
            material.addValue(name, value);
        end
        
        function addValue(material, name, value)
            material.parameters(name) = value;
        end
        
        function value = getValue(material, name)
            value = material.parameters(name);
        end
        
        function value = getParameterValue(material, name)
            % alternative name for getValue for compatibility with the
            % PropertyContainer
            value = material.getValue(name);
        end
        
        function returnVal = hasValue(material, name)
            %HASVALUE returns true, if a value with NAME is specified in
            %the container. Returns false, if not
            returnVal = true;
            try
                material.parameters(name);
            catch
                returnVal = false;
            end
        end
        
    end
    
end

