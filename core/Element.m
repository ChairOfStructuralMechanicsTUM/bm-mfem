classdef (Abstract) Element < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %ELEMENT The element class
    %   Abstract base class for all element implementations
    
    properties (Access = public) % currently changed from private to public
        id
        material
        
    end
    properties (Access = public) % currently changed from private to public
        nodeArray
        dofNames
    end
    
    methods
        % constructor
        function element = Element(id, material)
            if (nargin > 0)
                element.id = id;
                if (isa(material,'Material'))
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
        
        function nodes = getNodes(element)
            nodes = element.nodeArray;
        end
        
    end
    
    methods (Access = protected)
        
        function cp = copyElement(obj)
           cp = copyElement@matlab.mixin.Copyable(obj);
           obj.id = obj.id + 100;
        end
        
        function addDofs(element, dofNames)
            for itNode = 1:length(element.nodeArray)
                nodalDofs(1, length(dofNames)) = Dof;
                for itDof = 1:length(dofNames)
                    newDof = Dof(element.nodeArray(itNode),0.0,dofNames(itDof));
                    nodalDofs(itDof) = newDof;
                end
                element.nodeArray(itNode).setDofArray(nodalDofs);
            end
        end
        
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
 
   %%% START -- Substructure_1
    methods (Access = public)
        function [elementArrayLeft, elementArrayRight]= divideElements(elementArray,xBoundary)
          elementArrayRight=[];
          elementArrayLeft=[];
          
        
        for ii=1:length(elementArray)
              currentElement=elementArray(ii);
              currentNodes = currentElement.getNodes;
              currentXCoords = currentNodes.getX;

              if currentXCoords(2) >= currentXCoords(1)    % get highest X-coordinat of the current 
                  xHigh= currentXCoords(2);                % element in orde to determ wether element is
                  xLow = currentXCoords(1);                % left or right of X-Boundary 
                  else
                      xHigh = currentXCoords(1);
                      xLow  = currentXCoords(2);
              end

                  if xHigh <= xBoundary       % if current element is left or on Boundary
                    elementArrayLeft = [elementArrayLeft currentElement]; % add element to left  elementArray                 
                  elseif xHigh == xBoundary && xLow == xBoundary
                     % setCrossSectionArea(currentElement,0.5*getCrossSectionArea(currentElement));
                      elmentArrayLeft   = [elementArrayLeft currentElement];
                      
                  elseif xHigh > xBoundary                                      
                     elementArrayRight = [elementArrayRight currentElement];
                  end
                  % maybe copy function should be used here! 
        end
        end
    end
%%% END---Substructure_1    
    
end

