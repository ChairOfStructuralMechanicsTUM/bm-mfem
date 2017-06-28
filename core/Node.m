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
        
        %set Id Function
        function setId(node, id)
            node.id = id;
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
        
        function setDofLoad(nodes, dof, load)
            %SETINITIALDOFLOAD set the load of a specific dof
            % parameters: dof, load
            
            for ii = 1:length(nodes)
                dofNames = arrayfun(@(dof) dof.getValueType, nodes(ii).dofArray);
                index = strfind(dofNames,dof,'ForceCellOutput',false);
                index = find(~cellfun(@isempty,index));
                nodes(ii).dofArray(index).setLoad(load);
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
        
        %function to find elements which need to be copied for
        %substructuring. nodeIntf: nodes at Interface, nodes02: all nodes
        %in 2. substructure, totalElementArray: all elements
        function elementsToCopy = elementsForCopy(nodeIntf, nodes02, totalElementArray)
                nodes = [nodeIntf, nodes02];
                elementsToCopy = [];
                for kk = 1:length(nodeIntf)
                    for jj = 1:length(nodes)
                        if kk ~= jj
                            %all combinations of interface and other nodes
                            %to find elements for copy
                            nodePair = [nodeIntf(kk) nodes(jj)];
                            elementsToCopy = unique([elementsToCopy ... 
                                findElements(nodePair, totalElementArray)], 'stable');   
                        end
                    end
                end     
        end
       
        %function to half point loads at the interface. this implies that 
        %half the point load is at substructer01 the other at substructer02
        function halfPointLoads(nodes)
            for ii = 1:length(nodes)
                dof = nodes(ii).getDofArray;
                direction = zeros(1,length(dof));
                value = 0;
                for jj = 1:length(dof)
                    %val = dof(jj).getValue;
                    val = dof(jj).getDofLoad;
                    %if value of a direction not zero, safe that direction
                    if val ~= 0
                        direction(jj) = val/abs(val);
                        value = abs(val);
                    end
                end
                %add a new point load which overrides the old one
                addPointLoad(nodes(ii), 0.5*value, direction);
            end
        end


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
                            elements = unique([elements, allElements(ll)]);
                        end
                    end
                end
            end
        end
        
        %function that orders nodes into left and right part relativ to the
        %interface
        function [nodesLeft, nodesRight] = splitNodesX(nodeIntf, totalNodeArray)
            %interfaceNode furthest to right
            maxX = max(nodeIntf.getX());
            
            nodesLeft = [];
            nodesRight = [];
            
            for ii = 1:length(totalNodeArray)
                %see whether one node is more to right than interface
                if maxX >= totalNodeArray(ii).getX()
                    %all nodes left or at interface
                    nodesLeft = [nodesLeft totalNodeArray(ii)];
                else
                    %all nodes right of interface
                    nodesRight = [nodesRight totalNodeArray(ii)];
                end
            end
            
            %get max nodeId
            maxId = max(totalNodeArray.getId);
            %copy all Nodes at the Interface
            cp = copyElement(nodeIntf, maxId);
            nodesRight = [nodesRight, cp];
        end
        
        %function that orders the nodes in upper and lower part relative to
        %interface
        function [nodesDown, nodesUp] = splitNodesY(nodeIntf, totalNodeArray)
            %interface Node furthest up
            maxY = max(nodeIntf.getY());
            
            nodesUp = [];
            nodesDown = [];
            
            for ii = 1:length(totalNodeArray)
                %see whether one node is lower than interface
                if maxY >= totalNodeArray(ii).getY()
                    %all nodes below interface
                    nodesDown = [nodesDown totalNodeArray(ii)];
                else
                    %all nodes above interface
                    nodesUp = [nodesUp totalNodeArray(ii)];
                end
            end
            
            %get max nodeId
            maxId = max(totalNodeArray.getId);
            %copy all Nodes at the Interface
            cp = copyElement(nodeIntf, maxId);
            nodesUp = [nodesUp, cp];
        end
        
        %function that orders the elements in a right and a left half
        %relative to the interface
        function [elementsLeft, elementsRight] = splitElementsX(nodesLeft, nodesRight, totalElementArray)
            %all elements left or at interface
            elementsLeft = unique(findElements(nodesLeft, totalElementArray));
            %all elements right of interface
            elementsRight = unique(findElements(nodesRight, totalElementArray));
        end
        
        %function that orders the elements in an upper and a lower half
        %relative to the interface
        function [elementsUp, elementsDown] = splitElementsY(nodesUp, nodesDown, totalElementArray)
            %all elements left or at interface
            elementsUp = unique(findElements(nodesUp, totalElementArray));
            %all elements right of interface
            elementsDown = unique(findElements(nodesDown, totalElementArray));
        end
        
        %adds loads from boundary conditions to new nodes which are a
        %copy of an old interface node, and set up dofArray
        function setDof(obj, cp)
            dofs = getDofArray(obj);
            %find dofValues/Type/Load of old node which is copied
            for ii = 1:length(dofs)
                load(ii) = dofs(ii).getDofLoad;
                value(ii) = dofs(ii).getValue;
                type(ii) = dofs(ii).getValueType;
                fixed(ii) = dofs(ii).isFixed;
            end
            %set dofValue/Type of new node
            if length(dofs) == 2
                newDofs = [Dof(cp,value(1),type(1)) Dof(cp,value(2),type(2))];
                setDofArray(cp, newDofs);
            else
                newDofs = [Dof(cp,value(1),type(1)) Dof(cp,value(2),type(2)) Dof(cp,value(3),type(3))];
                setDofArray(cp, newDofs);
            end
            dofs = cp.getDofArray;
            %set dofLoad of new node
            if load ~= 0
                for jj = 1:length(dofs)
                    setLoad(dofs(jj), load(jj));
                end
            end
            %fix dofs if necessary
            for kk = 1:length(dofs)
                if fixed(kk) ~= 0
                    fix(dofs(kk));
                end
            end
        end
        
        %function to find copy of a node (same coordinates, different Id)
        function copy = findCopy(orig, nodes)
            for ii = 1:length(nodes)
                if orig.getX == nodes(ii).getX && orig.getY == nodes(ii).getY ...
                    && orig.getId ~= nodes(ii).getId && orig.getZ == nodes(ii).getZ
                    copy = nodes(ii);    
                    return
                else
                    copy = orig;
                end
            end
        end
    end
    
    methods (Access = protected) 
        %copy one/multiple nodes. also apply loads from original system to
        %new nodes
        function cp = copyElement(obj, maxId)
           % copy constructor
           % cp = copyElement@matlab.mixin.Copyable(obj);
           maxId = maxId+1;
           
           %only one Node
           if (length(obj) == 1)
               coords = obj.getCoords;
               if (length(coords) == 2)
                   cp = Node(maxId, coords(1), coords(2));
                   setDof(obj, cp);
               else
                   cp = Node(maxId, coords(1), coords(2), coords(3)); 
                   setDof(obj, cp);
               end
           
           %multiple nodes
           else
               cp =[];
               for jj = 1:length(obj)
                   coords = obj(jj).getCoords;
                   if (length(coords) == 2)
                       cp = [cp Node(maxId, coords(1), coords(2))];
                       setDof(obj(jj), cp(jj));
                   else
                       cp = [cp Node(maxId, coords(1), coords(2), coords(3))];
                       setDof(obj(jj), cp(jj));
                   end
                   %increase the Id
                   maxId = maxId+1;
               end
               
           end
        end
    end 
    %%%End NEW
end

