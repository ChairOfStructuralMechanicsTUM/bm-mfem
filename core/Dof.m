classdef Dof < handle
    %DOF The degree of freedom class
    %   Detailed explanation goes here
    
    properties (Access = private)
        id
        node
        value
        valueType
        dofLoad
        fixed = false
        
        firstDerivativeValue
        secondDerivativeValue
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
                dof.firstDerivativeValue = 0;
                dof.secondDerivativeValue = 0;
                dof.valueType = valueType;
                dof.id = -1; %uninitialized               
            end
        end
        
        % getter functions
        function node = getNode(dof)
           node = dof.node; 
        end
        
        function value = getValue(dofs, step)
            if (nargin == 1) || (strcmp(step, 'all'))
                nsteps = length(dofs(1).value);
                value = zeros(1,nsteps);
                for ii = 1:length(dofs)
                    values = dofs(ii).value;
                    if length(values) ~= nsteps
                        error('the dof value array of dof %d has a length of %d rather than the expected %d', ...
                            dofs(ii).getId, length(values), nsteps)
                    end
                    value(ii,:) = dofs(ii).value;
                end
                
            elseif nargin == 2
                if strcmp(step, 'end')
                    step = length(dofs(1).value);
                else
                    %check, if step exists
                    try
                        dofs(1).value(step);
                    catch
                        error('the specified step=%d does not exist for the dof', step)
                    end
                end
                value = zeros;
                for ii = 1:length(dofs)
                    value(ii,:) = dofs(ii).value(step);
                end
            end
        end
        
        function val = getFirstDerivativeValue(dof, step)
            [~, val, ~] = dof.getAllValues(step);
        end
        
        function val = getSecondDerivativeValue(dof, step)
            [~, ~, val] = dof.getAllValues(step);
        end
        
        function [value, firstDerivativeValue, secondDerivativeValue] = getAllValues(dofs, step)
            if (nargin == 1) || (strcmp(step, 'all'))
                nsteps = length(dofs(1).value);
                value = zeros(1,nsteps);
                firstDerivativeValue = zeros(1,nsteps);
                secondDerivativeValue = zeros(1,nsteps);
                
                for ii = 1:length(dofs)
                    values = dofs(ii).value;
                    if length(values) ~= nsteps
                        error('the dof value array of dof %d has a length of %d rather than the expected %d', ...
                            dofs(ii).getId, length(values), nsteps)
                    end
                    value(ii,:) = dofs(ii).value;
                    firstDerivativeValue(ii,:) = dofs(ii).firstDerivativeValue;
                    secondDerivativeValue(ii,:) = dofs(ii).secondDerivativeValue;
                end
                
            elseif nargin == 2
                if step == 'end'
                    step = length(dofs(1).value);
                else
                    %check, if step exists
                    try
                        dofs(1).value(step);
                    catch
                        error('the specified step=%d does not exist for the dof', step)
                    end
                end
                
                value = zeros;
                firstDerivativeValue = zeros;
                secondDerivativeValue = zeros;
                
                for ii = 1:length(dofs)
                    value(ii,:) = dofs(ii).value(step); 
                    firstDerivativeValue(ii,:) = dofs(ii).firstDerivativeValue(step);
                    secondDerivativeValue(ii,:) = dofs(ii).secondDerivativeValue(step);
                end
            end
        end
        
        
        function dofLoad = getDofLoad(dof)
            if isempty(dof.dofLoad)
                dofLoad = 0;
            else 
                dofLoad = dof.dofLoad;
            end
        end
 
        function id = getId(dofs)
           id = zeros;
           for ii = 1:length(dofs)
               id(ii) = dofs(ii).id;
           end
        end
        
        function valueType = getValueType(dof)
            valueType = dof.valueType;
        end
        
        function fixed = isFixed(dof)
            fixed = dof.fixed;
        end

        % setter functions
        function fix(dofs)
            for ii = 1:length(dofs)
                dofs(ii).fixed = true;
            end
        end

        function unfix(dofs)
            for ii = 1:length(dofs)
                dofs(ii).fixed = false;
            end
        end

        function setValue(dofs, values, step)
        %SETVALUE sets the value of the specified dof. If no step is
        %provided, the most recent step value is overwritten by value.
        %To append a new value, you can use DOF.APPENDVALUE
        %See also DOF.APPENDVALUE
            
            if length(dofs) ~= size(values,1)
                error('the arrays of dofs and values are not of the same size')
            end
            
            if nargin == 2
                for ii = 1:length(dofs)
                    dofs(ii).value(end) = values(ii);
                end
            elseif nargin == 3
                for ii = 1:length(dofs)
                   dofs(ii).value(step) = values(ii); 
                end
            end
        end
        
        function setFirstDerivativeValue(dofs, values, step)
        %SETFIRSTDERIVATIVEVALUE sets the value of the first derivative of
        %the specified dof.
        %See also DOF.SETVALUE
            
            if length(dofs) ~= size(values,1)
                error('the arrays of dofs and values are not of the same size')
            end
            
            if nargin == 2
                for ii = 1:length(dofs)
                    dofs(ii).firstDerivativeValue(end) = values(ii);
                end
            elseif nargin == 3
                for ii = 1:length(dofs)
                   dofs(ii).firstDerivativeValue(step) = values(ii); 
                end
            end
        end
        
        function setSecondDerivativeValue(dofs, values, step)
        %SETSECONDDERIVATIVEVALUE sets the value of the second derivative of
        %the specified dof.
        %See also DOF.SETVALUE
            
            if length(dofs) ~= size(values,1)
                error('the arrays of dofs and values are not of the same size')
            end
            
            if nargin == 2
                for ii = 1:length(dofs)
                    dofs(ii).secondDerivativeValue(end) = values(ii);
                end
            elseif nargin == 3
                for ii = 1:length(dofs)
                   dofs(ii).secondDerivativeValue(step) = values(ii); 
                end
            end
        end
        
        function appendValue(dofs, values)
        %APPENDVALUE appends the given values to the value array of the
        %dofs. No old step values are overwritten.
        %See also DOF.SETVALUE
            if length(dofs) ~= length(values)
                error('the arrays of dofs and values are not of the same size')
            end
            
            for ii = 1:length(dofs)
               dofs(ii).value(end+1) = values(ii); 
            end
        end
        
        function appendFirstDerivativeValue(dofs, values)
        %APPENDFIRSTDERIVATIVEVALUE appends the given values to the array 
        %of the first derivative values of the dofs. No old step values are 
        %overwritten.
        %See also DOF.SETVALUE
            if length(dofs) ~= length(values)
                error('the arrays of dofs and values are not of the same size')
            end
            
            for ii = 1:length(dofs)
               dofs(ii).firstDerivativeValue(end+1) = values(ii); 
            end
        end
        
        function appendSecondDerivativeValue(dofs, values)
        %APPENDSECONDDERIVATIVEVALUE appends the given values to the array 
        %of the second derivative values of the dofs. No old step values are 
        %overwritten.
        %See also DOF.SETVALUE
            if length(dofs) ~= length(values)
                error('the arrays of dofs and values are not of the same size')
            end
            
            for ii = 1:length(dofs)
               dofs(ii).secondDerivativeValue(end+1) = values(ii); 
            end
        end
        
        function removeValue(dofs, step)
           %REMOVEVALUE removes the value from the given step
           for ii = 1:length(dofs)
               dofs(ii).value(step) = [];
           end
        end
        
        function setLoad(dof, dofLoad)
            dof.dofLoad = dofLoad;
        end

        function setId(dof, id)
            dof.id = id;
        end

       
    end
    
end

