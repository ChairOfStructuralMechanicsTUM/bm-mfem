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
        
        function setStepValue(propertyContainer, name, value, step)
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
        
        function appendStepValue(propertyContainer, name, value)
           %APPENDSTEPVALUE appends the given VALUE to the property NAME.
           try
                vals = propertyContainer.propertyMap(name);
            catch
                error('The value \"%s\" is not available in this container', name)
            end
            vals(end+1) = value;
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
        
        function returnVal = hasValue(propertyContainer, name)
            %HASVALUE returns true, if a value with NAME is specified in
            %the container. Returns false, if not
            returnVal = true;
            try
                propertyContainer.propertyMap(name);
            catch
                returnVal = false;
            end
        end
        
        function names = getValueNames(propertyContainer)
            names = propertyContainer.propertyMap.keys;
        end
        
        function print(pc)
            %PRINT shows the keys and values of the property container
            names = pc.propertyMap.keys;
            fprintf('property container with %d values:\n', length(names));
            fprintf('name \t\t value\n');
            fprintf('----------------------\n');
            for ii = 1:length(names)
                name = cell2mat(names(ii));
                fprintf('%s \t\t %d\n', name, pc.propertyMap(name));
            end
            fprintf('\n\n');
        end
        
        function cp = copy(obj)
            %COPY performs a deep copy of a property container
            cp = PropertyContainer();
            map = obj.getPropertyMap();
            k = keys(map);
            v = values(map);
            for ii = 1:length(map)
               cp.setValue(k{ii},v{ii}); 
            end
        end
        
    end
    
end

