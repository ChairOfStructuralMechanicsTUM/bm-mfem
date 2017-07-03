classdef (Abstract) Element < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %ELEMENT The element class
    %   Abstract base class for all element implementations
    
    properties (Access = private)
        id
        material
    end
    properties (Access = protected)
        nodeArray
        dofNames
        requiredProperties
        required3dProperties
    end
    
    methods
        % constructor
        function element = Element(id, material)
            if (nargin > 0)
                element.id = id;
                if (isa(material,'Material') || isa(material,'PropertyContainer'))
                    element.material = material;
                else
                    error('problem with the material in element %d', id);
                end
                element.nodeArray = {};
            end
        end
    end
    
    methods (Abstract)
        update(element)     % update properties after e.g. nodes changed
        barycenter(element)
        computeLocalStiffnessMatrix(element)
        computeLocalForceVector(element)
    end
    
    methods (Sealed)
        % getter functions
        function id = getId(element)
            id = zeros;
            for ii = 1:length(element)
                id(ii) = element(ii).id;
            end
        end
        
        function material = getMaterial(element)
            material = element.material;
        end
        
        function value = getPropertyValue(element, valueName)
           value = element.material.getValue(valueName); 
        end
        
        function nodes = getNodes(elements)
            nodes = Node.empty;
            for ii = 1:length(elements)
                nodes(ii,:) = elements(ii).nodeArray;
            end
%             nodes = element.nodeArray;
        end
        
        function dofs = getDofs(element)
           dofs = Dof.empty;
           nodes = element.nodeArray;
           for ii = 1:length(nodes)
              dofs = [dofs nodes(ii).getDofArray'];
           end
        end
        
        function setPropertyValue(elements, valueName, value)
            for ii = 1:length(elements)
                elements(ii).material.setValue(valueName, value);
            end
        end
            
        
        function check(element)
        %CHECK checks, if all dofs and properties required by the element are available
            
            %check the dofs
            for iNode = 1:length(element.nodeArray)
                cNode = element.nodeArray(iNode);
                availableDofNames = arrayfun(@(dof) dof.getValueType, cNode.getDofArray, 'UniformOutput',false);
                diff = setxor(element.dofNames', availableDofNames);
                if ~ isempty(diff)
                    missingDofs = setdiff(element.dofNames', availableDofNames);
                    unknownDofs = setdiff(availableDofNames, element.dofNames');
                    
                    error('the following dofs are missing at node %d: %s\nthe following dofs at node %d are not defined for the element %s: %s\n', ...
                        cNode.getId, ...
                        strjoin(missingDofs,', '), ...
                        cNode.getId, ...
                        class(element), ...
                        strjoin(unknownDofs,', '))
                end
            end
            
            %check the properties
            valsToCheck = element.requiredProperties;
            properties = element.getMaterial;
            availableValueNames = properties.getValueNames;
            for ii = 1:length(valsToCheck)
                if ~ any(ismember(valsToCheck(ii), availableValueNames))
                    error('error in element %d: property %s is missing', element.id, cell2mat(valsToCheck(ii)))
                end
            end
            
            valsToCheck3d = element.required3dProperties;
            properties = element.getMaterial;
            availableValueNames = properties.getValueNames;
            for ii = 1:length(valsToCheck3d)
                if any(ismember(valsToCheck3d(ii), availableValueNames))
                    val = properties.getValue(cell2mat(valsToCheck3d(ii)));
                    if ~ (length(val) == 3)
                        error('error in element %d: property %s must have 3 values', element.id, cell2mat(valsToCheck3d(ii)))
                    end
                else
                    error('error in element %d: property %s is missing', element.id, cell2mat(valsToCheck3d(ii)))
                end
            end
        end
        
    end
    
    methods (Access = protected)
        
        function cp = copyElement(obj)
           cp = copyElement@matlab.mixin.Copyable(obj);
           obj.id = obj.id + 100;
        end
        
%         function addDofs(element, dofNames)
%             for itNode = 1:length(element.nodeArray)
%                 nodalDofs(1, length(dofNames)) = Dof;
%                 for itDof = 1:length(dofNames)
%                     newDof = Dof(element.nodeArray(itNode),0.0,dofNames(itDof));
%                     nodalDofs(itDof) = newDof;
%                 end
%                 element.nodeArray(itNode).setDofArray(nodalDofs);
%             end
%         end
        
        function addDofsToSingleNode(element, node)
            nodalDofs(1, length(element.dofNames)) = Dof;
            for itDof = 1:length(element.dofNames)
                newDof = Dof(node,0.0,element.dofNames(itDof));
                nodalDofs(itDof) = newDof;
            end
            node.setDofArray(nodalDofs);
        end
        
    end
    
    methods
        
        function overwriteNode(element, oldNode, newNode)
            if ~isa(newNode,'Node')
                error('invalid node')
            end
            
            for itNode = 1:length(element.nodeArray)
                if (element.nodeArray(itNode) == oldNode)
                    
                    element.nodeArray(itNode) = newNode;
                    element.addDofsToSingleNode(newNode);
                    element.update;
                    
%                  barElement3d2n.addDofs(barElement3d2n.dofNames);
%                 
%                 barElement3d2n.length = computeLength(barElement3d2n.nodeArray(1).getCoords, ...
%                     barElement3d2n.nodeArray(2).getCoords);
                end
            end
        end
        
    end
    
end

