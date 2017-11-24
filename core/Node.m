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
        
        %find elements to copy from nodes
        function elementsToCopy = elementsForCopy(nodeIntf, nodes02, totalElementArray)
            %function to find elements which need to be copied for
            %substructuring. nodeIntf: nodes at Interface, nodes02: all nodes
            %in 2. substructure, totalElementArray: all elements
            
                nodes = [nodeIntf, nodes02];
                elementsToCopy = [];
                
                %iterate over all nodes at interface
                for kk = 1:length(nodeIntf)
                    %iterate over all nodes (2. substructure + interface)
                    for jj = 1:length(nodes)
                        if kk ~= jj
                            %all combinations of interface and other nodes
                            %to find elements for copy
                            nodePair = [nodeIntf(kk) nodes(jj)];
                            %Find elements that are at the interface and
                            %need to be copied. Filter out duplicates.
                            elementsToCopy = unique([elementsToCopy ... 
                                findElements(nodePair, totalElementArray)], 'stable');   
                        end
                    end
                end     
        end
        
        %sort nodes
        function nodeIntf = sortNodes(nodeIntfNoSort)
            %function sorting nodes in an array by their id
            
            nodeIntf = [];
            %iterate over unsorted nodes
            for ii = 1:length(nodeIntfNoSort)
                %find smallest node of "old" unsorted array
                mini = min(nodeIntfNoSort.getId);
                ll = length(nodeIntfNoSort);
                jj = 1;
                while jj <= ll
                    if nodeIntfNoSort(jj).getId == mini
                        %add smallest node to sorted array
                        nodeIntf = [nodeIntf nodeIntfNoSort(jj)];
                        %delete smallest node, now the smallest node in the
                        %unsorted array is the second smallest from before
                        nodeIntfNoSort(jj) = [];
                        ll = ll-1;
                        jj = jj-1;
                    end
                    jj = jj+1;
                end
            end
        end
               
        function halfPointLoads(nodes)
            %Function to half point loads at the interface. This implies that 
            %half the point load is at substructer01 the other at substructer02
            
            %find loads for each node
            for ii = 1:length(nodes)
                dofs = nodes(ii).getDofArray;
                for kk = 1:length(dofs)
                    loads(kk) = dofs(kk).getDofLoad;
                end
                temp = abs(loads);
                
                %smallest non-zero load at a node
                val = min(temp(temp>0));
                direction = zeros(1,length(loads));
                
                if (~isempty(find(loads)))
                    %Set new directions for loads. The directions imply the
                    %magnitude of each loads as well.
                    for jj = 1:length(loads)
                        direction(jj) = loads(jj)/val;
                    end
                    
                    %multiply with norm of direction because later in
                    %addPointLoad it is devided by this value
                    addPointLoad(nodes(ii), norm(direction)*0.5*val, direction);
                end
            end
        end

        %find element from node
        function  elements = findElements(nodeArray, allElements)
            %function that finds elements from nodes. nodeArray: are the
            %nodes for which the elements needs to be found
            
            elements = [];
            %iterate over nodes twice and try to find node Pairs
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
        
        %split nodes in x-direction
        function [nodesLeft, nodesRight] = splitNodesX(nodeIntf, totalNodeArray, idt)
            %function that orders nodes into left and right part relativ to the
            %interface. interfaceNode is node furthest to right
            
            %node furthest to the right
            maxX = max(nodeIntf.getX());
            
            nodesLeft = [];
            nodesRight = [];
            
            %iterate over all nodes
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

            %copy the nodes on the interface
            cp = copyElement(nodeIntf, idt);
            %add the copied interface nodes to the right half
            nodesRight = [nodesRight, cp];
        end
        
        %split nodes in y-direction
        function [nodesDown, nodesUp] = splitNodesY(nodeIntf, totalNodeArray, idt)
            %function that orders the nodes in upper and lower part relative to
            %interface. Interface Node is node furthest up.
            
            %node furthest up
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
            
            %copy the nodes on the interface
            cp = copyElement(nodeIntf, idt);
            %add the copied interface nodes to the right half
            nodesUp = [nodesUp, cp];
        end
        
        %split elements in x-direction
        function [elementsLeft, elementsRight] = splitElementsX(nodesLeft, nodesRight, totalElementArray)
            %function that orders the elements in a right and a left half
            %relative to the interface
            
            %all elements left or at interface
            elementsLeft = unique(findElements(nodesLeft, totalElementArray));
            %all elements right of interface
            elementsRight = unique(findElements(nodesRight, totalElementArray));
        end
        
        %split elements in y-direction
        function [elementsUp, elementsDown] = splitElementsY(nodesUp, nodesDown, totalElementArray)
            %function that orders the elements in an upper and a lower half
            %relative to the interface
            
            %all elements left or at interface
            elementsUp = unique(findElements(nodesUp, totalElementArray));
            %all elements right of interface
            elementsDown = unique(findElements(nodesDown, totalElementArray));
        end
        
        %update loads in copied nodes
        function setDof(obj, cp)
            %adds loads from boundary conditions to new nodes which are a
            %copy of an old interface node, and set up dofArray
            
            dofs = getDofArray(obj);
            %find dofValues/Type/Load of old node which is copied
            for ii = 1:length(dofs)
                load(ii) = dofs(ii).getDofLoad;
                value(ii) = dofs(ii).getValue;
                type(ii) = dofs(ii).getValueType;
                fixed(ii) = dofs(ii).isFixed;
            end

            %set dofValue/Type of new node and distinguish 2D, 3D case
            if length(dofs) == 2
                newDofs = [Dof(cp,value(1),type(1)) Dof(cp,value(2),type(2))];
                setDofArray(cp, newDofs);
            else
                newDofs = [Dof(cp,value(1),type(1)) Dof(cp,value(2),type(2)) Dof(cp,value(3),type(3))];
                setDofArray(cp, newDofs);
            end
            dofs = cp.getDofArray;
            
            %set dofLoad of new node
            if (~isempty(find(load)))   
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
        
        %find copies
        function copy = findCopy(orig, nodes)
            %function to find copy of a node (same coordinates, different Id)
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
        
        %find interface nodes
        function  nodeIntf = findIntfNode(orig, nodes)
            for ii = 1:length(nodes)
                if orig.getX == nodes(ii).getX && orig.getY == nodes(ii).getY ...
                    && orig.getZ == nodes(ii).getZ
                    nodeIntf = orig;
                    return;
                else
                    nodeIntf = Node(0,0,0,0);
                end
            end
        end
    end
    
    methods (Access = protected) 
        
        %copy one/multiple nodes. also apply loads from original system to
        %new nodes
        function cp = copyElement(obj, idt)
           % copy constructor
           % cp = copyElement@matlab.mixin.Copyable(obj);
           maxId = idt.getNodeId+1;
           
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
                   %increase maxId, done here because otherwise we have one
                   %Id too many in the tracker later
                   if jj > 1
                       maxId = maxId+1;
                   end
                   coords = obj(jj).getCoords;
                   if (length(coords) == 2)
                       cp = [cp Node(maxId, coords(1), coords(2))];
                       setDof(obj(jj), cp(jj));
                   else
                       cp = [cp Node(maxId, coords(1), coords(2), coords(3))];
                       setDof(obj(jj), cp(jj));
                   end
               end
           end
           %update IdTracker
           idt.setNodeId(maxId);
        end
    end 
end
