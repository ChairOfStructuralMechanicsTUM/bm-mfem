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
        function [elementArrayLeft, elementArrayRight, InterfaceNodes]= divideElements(elementArray,dim,Boundary)
          elementArrayRight=[];
          elementArrayLeft=[];
          InterfaceNodes=[];
          
          for ii=1:length(elementArray)
              
              currentNodes = elementArray(ii).getNodes;
              currentCoords = currentNodes.getCoords;           %[ x1 x2 y1 y2 z1 z2]
              currentCoords = [currentCoords(2*dim-1),currentCoords(2*dim)];
              
              if currentCoords(2) >= currentCoords(1)    % get highest coordinat of the current element
                  highCoord= currentCoords(2);
                  lowCoord  = currentCoords(1);           % to compare element position                %
                  lowCoordPosition=1;                   %lowCoordPosition is needed to overwrite the Node in Case of Boundary position
              else
                  highCoord = currentCoords(1);
                  lowCoord  = currentCoords(2);
                  lowCoordPosition=2;
              end
              
              if highCoord <= Boundary  && lowCoord < Boundary         % if current element is left or on Boundary
                  elementArrayLeft = [elementArrayLeft elementArray(ii)]; % add element to left  elementArray
                  
              elseif highCoord == Boundary && lowCoord == Boundary
                                           
                      elementArrayLeft = [elementArrayLeft copyElement(elementArray(ii))];
                      elementArrayRight = [elementArrayRight copyElement(elementArray(ii))];
                      % Set Cross-Section-Area
                      setCrossSectionArea(elementArrayRight(length(elementArrayRight)),...
                          0.5*getCrossSectionArea(elementArray(ii)));
                      setCrossSectionArea(elementArrayLeft(length(elementArrayLeft)),...
                          0.5*getCrossSectionArea(elementArray(ii)));
                      
                      % Store Interface Information 
                      InterfaceNodes=[InterfaceNodes currentNodes];
                      
                      %overwrite  both Nodes of Copied Element: old Node = current
                      %Node ; new Node = copy of current Node
                                         
                      for j=1:2
                          overwriteNode(elementArrayRight(length(elementArrayRight)),...
                              currentNodes(j),copyElement(currentNodes(j)));
                      end
                      
                  
                  
              elseif highCoord > Boundary && lowCoord ==Boundary
                  elementArrayRight = [elementArrayRight elementArray(ii)];
                    
                % overwrite Boundary node
                     overwriteNode(elementArrayRight(length(elementArrayRight)),...
                     currentNodes(lowCoordPosition),copyElement(currentNodes(lowCoordPosition)));
              
                 % Store Interface Information       
                 InterfaceNodes=[InterfaceNodes currentNodes(lowCoordPosition)];
              
              else
                 elementArrayRight = [elementArrayRight elementArray(ii)]; 
              end
          end
        end
    end
    
%%% END---Substructure_1    
    
end

