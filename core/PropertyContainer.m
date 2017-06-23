classdef PropertyContainer < handle
    %PROPERTYCONTAINER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        propertyMap
    end
    
    methods
        % constructor
        function propertyContainer = PropertyContainer()
            propertyContainer.propertyMap = containers.Map; 
        end
        
        % getter functions
        
        function propertyMap = getPropertyMap(propertyContainer)
            propertyMap = propertyContainer.propertyMap;
        end
        
        % member functions
        function setValue(propertyContainer, name, value)
            propertyContainer.propertyMap(name) = value;
        end
        
        function value = getValue(propertyContainer, name)
            try
                value = propertyContainer.propertyMap(name);
            catch
                error('The value \"%s\" is not available in this container', name)
            end
        end
        
        function names = getValueNames(propertyContainer)
            names = propertyContainer.propertyMap.keys;
        end
        
    end
    
end

