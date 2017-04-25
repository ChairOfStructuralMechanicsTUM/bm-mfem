classdef Node < handle
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
            id = node.id;
        end
        
        function coords = getCoords(node)
           coords = [node.x node.y node.z]; 
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
         
        
    end
    
    methods (Access = private)
        
    end
    
end

