classdef Node < handle & matlab.mixin.Copyable
    %NODE The node class
    %   Parameters:
    %       id: unique identifier
    %       x, y(, z): coordinates
    
    properties (Access = public) %!!! was private, to public?
        id
        x
        y
        z
        dofArray
        responseDofArray
%         loadMap = containers.Map;
    end
    
    methods
        % constructor
        function node = Node(id, x, y, z)
            switch nargin
                case 0
                    % the empty constructor is needed in order to
                    % preallocate empty arrays of Nodes
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
        function fixDof(nodes, dof)
            for ii = 1:length(nodes)
                dofNames = arrayfun(@(dof) dof.getValueType, nodes(ii).dofArray);
                index = strfind(dofNames,dof,'ForceCellOutput',false);
                index = find(~cellfun(@isempty,index));
                nodes(ii).dofArray(index).fix;
            end
        end
        
        function setDofValue(nodes, dof, load)
            %SETDOFVALUE set the value of a specific dof
            % parameters: dof, load
            for ii = 1:length(nodes)
                dofNames = arrayfun(@(dof) dof.getValueType, nodes(ii).dofArray);
                index = strfind(dofNames,dof,'ForceCellOutput',false);
                index = find(~cellfun(@isempty,index));
                nodes(ii).dofArray(index).setValue(load);
            end
        end
        
        function val = getDofValue(nodes, dof)
            val = zeros;
            for ii = 1:length(nodes)
                dofNames = arrayfun(@(dof) dof.getValueType, nodes(ii).dofArray);
                index = strfind(dofNames,dof,'ForceCellOutput',false);
                index = find(~cellfun(@isempty,index));
                val(ii) = nodes(ii).dofArray(index).getValue;
                %fprintf('node_%i: %s = %d\n',node.getId,dof,val);
            end
        end
        
        %%% Start NEW
        %function that finds elements from nodes
        function  elements = findElements(nodeArray, allElements)
            elements = [];
            for jj = 1:length(nodeArray)
                for kk = 1:length(nodeArray) 
                    if kk == jj
                        continue
                    else
                    nodePair = [nodeArray(jj), nodeArray(kk)];
                    end
                    for ll = 1:length(allElements)
                        if ismember(nodePair, allElements(ll).getNodes)
                            elements = [elements, allElements(ll)];
                        end
                    end
                end
            end
        end   
        
        %function that finds all nodes along specific x-coordinate
        function nodesX = nodesAtX(allNodes, maxX)
            nodesX = [];
            for mm = 1:length(allNodes)
                if allNodes(mm).getX == maxX
                    nodesX = [nodesX, allNodes(mm)];
                end
            end
        
        end
        %%%END NEW
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

