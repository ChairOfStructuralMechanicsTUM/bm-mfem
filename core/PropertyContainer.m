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
            %SETVALUE sets the value NAME in the property container. An
            %already existing property is overwritten.
            propertyContainer.propertyMap(name) = value;
        end
        
        function setStepValue(propertyContainer, name, step, value)
            %ADDSTEPVALUE adds VALUE to the property container NAME at the
            %specified STEP. An already existing value at this step is
            %overwritten
            try
                vals = propertyContainer.propertyMap(name);
            catch
                error('The value \"%s\" is not available in this container', name)
            end
            vals(step) = value;
            propertyContainer.setValue(name, vals);
        end
        
        function value = getValue(propertyContainer, name, step)
            %GETVALUE returns the value stored under NAME. If the value is
            %an array, the step can be specified.
            
            try
                value = propertyContainer.propertyMap(name);
            catch
                error('The value \"%s\" is not available in this container', name)
            end
            
            if nargin == 3
                if mod(step,1) ~= 0
                    error('Please specify a valid step')
                end
                value = value(step);
            end
        end
        
        function names = getValueNames(propertyContainer)
            names = propertyContainer.propertyMap.keys;
        end
        
    end
    
end

