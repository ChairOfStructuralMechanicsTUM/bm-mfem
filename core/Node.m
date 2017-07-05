classdef Node < handle & matlab.mixin.Copyable
    %NODE The node class
    %   Parameters:
    %       id: unique identifier
    %       x, y(, z): coordinates
    
    properties (Access = private)
        id
        x
        y
        z
        dofArray
        valueMap
    end
    
    methods
        % constructor
        function node = Node(id, x, y, z)
            switch nargin
                case 0
                    % the empty constructor is needed in order to
                    % preallocate empty arrays of Nodes
                    node.id = -1;
                case 3
                    node.id = id;
                    node.x = x;
                    node.y = y;
                case 4
                    node.id = id;
                    node.x = x;
                    node.y = y;
                    node.z = z;
                otherwise
                    error('Wrong number of arguments')
            end
            node.valueMap = PropertyContainer();
        end
        
        
        % getter functions
        function id = getId(node)
            % Return the id of the node
            id = zeros;
            for ii = 1:length(node)
                id(ii) = node(ii).id;
            end
        end
        
        function coords = getCoords(node)
            % Return all coordinates in the form x, y(, z)
            coords = [node.x node.y node.z];
        end
        
        function x = getX(nodes)
            x = zeros;
            for ii = 1:length(nodes)
                x(ii) = nodes(ii).x;
            end
        end
        
        function y = getY(nodes)
            y = zeros;
            for ii = 1:length(nodes)
                y(ii) = nodes(ii).y;
            end
        end
        
        function z = getZ(nodes)
            z = zeros;
            for ii = 1:length(nodes)
                z(ii) = nodes(ii).z;
            end
        end
        
        function dofs = getDofArray(node)
            dofs = node.dofArray;
        end
        
       
        %setter functions
        function setDofArray(node, dofs)
            node.dofArray = dofs;
        end
        
        
        % member functions
        function addDof(nodes, dofNames)
            for itNode = 1:length(nodes)
                availableDofNames = arrayfun(@(dof) dof.getValueType, nodes(itNode).dofArray);
                for itDof = 1:length(dofNames)
                    if ~ isempty(availableDofNames)
                        if any(ismember(availableDofNames,dofNames))
                            continue
                        end
                    end
                    nodes(itNode).dofArray = [nodes(itNode).dofArray 
                        Dof(nodes(itNode),0.0,dofNames{itDof})];
                end
%                 nodes(itNode).dofArray = nodes(itNode).dofArray';
%                 test = nodes(itNode).dofArray;
            end
        end
        
        function fixDof(nodes, dof)
            for ii = 1:length(nodes)
                dofNames = arrayfun(@(ndof) ndof.getValueType, nodes(ii).dofArray,'UniformOutput',false);
                index = strfind(dofNames,dof,'ForceCellOutput',false);
                index = find(~cellfun(@isempty,index));
                nodes(ii).dofArray(index).fix;
            end
        end
        
        function setDofValue(nodes, dofName, load)  %not load initial displacement
            %SETDOFVALUE set the value of a specific dof
            % parameters: dof, load
            for ii = 1:length(nodes)
                dofNames = arrayfun(@(dof) dof.getValueType, nodes(ii).dofArray,'UniformOutput',false);
                index = strfind(dofNames,dofName,'ForceCellOutput',false);
                index = find(~cellfun(@isempty,index));
                nodes(ii).dofArray(index).setValue(load);
            end
        end
        
        
        
        function setDofLoad(nodes, dofName, load)
            %SETINITIALDOFLOAD set the load of a specific dof
            % parameters: dof, load
            for ii = 1:length(nodes)
                dofNames = arrayfun(@(dof) dof.getValueType, nodes(ii).dofArray,'UniformOutput',false);
                index = strfind(dofNames,dofName,'ForceCellOutput',false);
                index = find(~cellfun(@isempty,index));
                nodes(ii).dofArray(index).setLoad(load);
            end
        end
        
        
        
        function val = getDofValue(nodes, dof, varargin)
            numvarargs = length(varargin);
            optargs = { 1 };
            optargs(1:numvarargs) = varargin;
            step = optargs{:};
            
            val = zeros;
            for ii = 1:length(nodes)
                dofNames = arrayfun(@(dof) dof.getValueType, nodes(ii).dofArray,'UniformOutput',false);
                index = strfind(dofNames,dof,'ForceCellOutput',false);
                index = find(~cellfun(@isempty,index));
                val(ii) = nodes(ii).dofArray(index).getValue(step);
                %fprintf('node_%i: %s = %d\n',node.getId,dof,val);
            end
        end
         
        
        function addNewValue(nodes, valueNames)
            for itName = 1:length(valueNames)
                for itNode = 1:length(nodes)
                    try
                        nodes(itNode).valueMap.getValue(char(valueNames(itName)),1);
                    catch
                        nodes(itNode).valueMap.setValue(char(valueNames(itName)), 0.0);
                    end
                end
            end
        end
        
        function setStepValue(nodes, valueName, value, step)
           for ii = 1:length(nodes)
              nodes(ii).valueMap.setStepValue(valueName, value, step); 
           end
        end
        
        function appendStepValue(nodes, valueName, value)
           for ii = 1:length(nodes)
              nodes(ii).valueMap.appendStepValue(valueName, value); 
           end
        end
        
        function val = getValue(nodes, valueName, step)            
            if nargin == 2
                for ii = 1:length(nodes)
                    val = nodes(ii).valueMap.getValue(valueName);
                end
            elseif nargin == 3
                 for ii = 1:length(nodes)
                    val = nodes(ii).valueMap.getValue(valueName, step);
                end
            end
        end
        
%         function init = getInitialDofLoad(nodes, dof)
%             init = zeros;
%             for ii = 1:length(nodes)
%                 dofNames = arrayfun(@(dof) dof.getValueType, nodes(ii).dofArray);
%                 index = strfind(dofNames,dof,'ForceCellOutput',false);
%                 index = find(~cellfun(@isempty,index));
%                 val(ii) = nodes(ii).dofArray(index).getLoad;
%             end
%         end
        
        
        
        
        
    end
    
    methods (Access = protected)
        
        function cp = copyElement(obj)
           % copy constructor
           %cp = copyElement@matlab.mixin.Copyable(obj);
            coords = obj.getCoords;
            if (length(coords) == 2)
                cp = Node(obj.getId, coords(1), coords(2));
            else
                cp = Node(obj.getId, coords(1), coords(2), coords(3));
            end
        end
    end
    
end

