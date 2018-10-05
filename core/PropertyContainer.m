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
        function addValue(obj, name, value)
           %ADDVALUE adds a value NAME and initializes it with the given
           %VALUE, if specified. Otherwise it is initialized with 0 for 1d
           %variables or [0,0,0] for 3d variables.           
           if obj.hasValue(name)
               error('The property \"%s\" already exists in the container!', name)
           end
           
           type = checkPropertyName(name);
           
           if strcmp(type,'variable1d')
               if nargin == 2
                   value = 0;
               end
               if length(value) ~= 1
                   error('The value for property \"%s\" must be a single number.', name)
               end
               obj.propertyMap(name) = value;
               
           elseif strcmp(type,'variable3d')
               if nargin == 2
                   value = [0 0 0];
               end
               if length(value) ~= 3
                   error('The value for property \"%s\" must be a vector of length 3.', name)
               end
               obj.propertyMap(name) = value;
           
           elseif strcmp(type, 'flag')
              if nargin == 2
                  value = false; 
              end
              obj.propertyMap(name) = value;
              
           else
               error('A property with name \"%s\" is not defined in mfem!', name)
           end
           
        end
        
        function setValue(propertyContainer, name, value)
            %SETVALUE sets the value NAME in the property container. An
            %already existing property is overwritten.
            try
                propertyContainer.propertyMap(name);
            catch
                error('The value \"%s\" is not available in this container', name)
            end
            propertyContainer.propertyMap(name) = value;
        end
        
        function setStepValue(propertyContainer, name, value, step)
            %ADDSTEPVALUE adds VALUE to the property container NAME at the
            %specified STEP. An already existing value at this step is
            %overwritten
            direction = 0;
            if strcmp(name(end-1:end),'_X')
                direction = 1;
                name = name(1:end-2);
            elseif strcmp(name(end-1:end),'_Y')
                direction = 2;
                name = name(1:end-2);
            elseif strcmp(name(end-1:end),'_Z')
                direction = 3;
                name = name(1:end-2);
            end
            
            type = checkPropertyName(name);
            
            try
                vals = propertyContainer.propertyMap(name);
            catch
                error('The value \"%s\" is not available in this container', name)
            end
            
            if (direction == 0) && (strcmp(type,'variable1d'))
                vals(step) = value;
            elseif (direction == 0) && (strcmp(type,'variable3d'))
                vals(step,:) = value;
            elseif direction ~= 0
                tempVals = vals(step,:);
                tempVals(direction) = value;
                vals(step,:) = tempVals;
            end
            propertyContainer.setValue(name, vals);
        end
        
        function appendStepValue(propertyContainer, name, value)
            %APPENDSTEPVALUE appends the given VALUE to the property NAME.
            %If the property is 3d, VALUE has to be a vector of length 3.
            try
                vals = propertyContainer.propertyMap(name);
            catch
                error('The value \"%s\" is not available in this container', name)
            end
            
            type = checkPropertyName(name);
            if (strcmp(type,'variable1d'))
                vals(end+1) = value;
            elseif (strcmp(type,'variable3d'))
                if length(value) ~= 3
                    error('Values for property \"%s\" must be of length 3.', name)
                end
                vals(end+1,:) = value;
            end
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
                value = value(step,:);
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
               cp.addValue(k{ii},v{ii}); 
            end
        end
        
    end
    
end

