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
       function [elementArrayLeft, elementArrayRight, interfaceNodes]= divideElements(elementArray,dim,Boundary)
           elementArrayRight=[];
           elementArrayLeft=[];
           interfaceNodes=[];
           
           for ii=1:length(elementArray)
               currentElement=elementArray(ii);
               currentNodes = currentElement.getNodes;
               currentCoords = currentNodes.getCoords;           %[ x1 x2 y1 y2 z1 z2]
               currentCoords = [currentCoords(2*dim-1),currentCoords(2*dim)];
               
               % get highest coordinat of the current element
               % to compare element position
               if currentCoords(2) >= currentCoords(1)    
                   highCoord = currentCoords(2);
                   lowCoord  = currentCoords(1);          
                   lowCoordPosition=1;                     
               else
                   highCoord = currentCoords(1);
                   lowCoord  = currentCoords(2);
                   lowCoordPosition=2;
               end
               % Case 1 one Node is Left or on of Boundary 
               if highCoord <= Boundary  && lowCoord < Boundary         % if current element is left or on Boundary
                   elementArrayLeft = [elementArrayLeft elementArray(ii)]; % add element to left  elementArray
                  
                % Case 2 element on Boundary    
               elseif highCoord == Boundary && lowCoord == Boundary
                    cp=copyElement(elementArray(ii));
                   elementArrayLeft = [elementArrayLeft cp];
                   elementArrayRight = [elementArrayRight cp ];

                   
                   % Set Cross-Section-Area
                   setCrossSectionArea(elementArrayRight(length(elementArrayRight)),...
                       0.5*getCrossSectionArea(elementArray(ii)));
                   setCrossSectionArea(elementArrayLeft(length(elementArrayLeft)),...
                       0.5*getCrossSectionArea(elementArray(ii)));
          
                   % Store Interface Information
                   interfaceNodes=[interfaceNodes currentNodes.getId];
                   
                   %overwrite  both Nodes of copied Element                 
                   for j=1:2
                       oldNode=currentNodes(j);
                       newNode=copyElement(currentNodes(j));
                       
                      % fix dofs that were fixed before
%                             for jj=1:length(oldNode.dofArray)
%                                 if  isFixed(oldNode.dofArray(jj))== 1
%                                     
%                                    newNode.dofArray(jj)=1;
%                                 else
%                                     newNode.dofArray(jj)=0;
%                                 end
%                             end
                            
                           overwriteNode(elementArrayRight(length(elementArrayRight)),...
                           oldNode,newNode);                     
                   end
                   
                   
                   % Case 3 one Node is right, the other on Boundary 
               elseif highCoord > Boundary && lowCoord ==Boundary
                   elementArrayRight =[elementArrayRight elementArray(ii)];
                   
                   % overwrite Boundary node
                   oldNode=currentNodes(lowCoordPosition);
                   newNode=copyElement(currentNodes(lowCoordPosition));
                   overwriteNode(elementArrayRight(length(elementArrayRight)),...
                       oldNode,newNode);
                   
                   % Fix Dofs that were fixed before
%                           for jj =1:length(currentNodes(lowCoordPosition).dofArray)
%                               if  isFixed(currentNodes(lowCoordPosition).dofArray(jj))== 'true'
%                                   fix(newNode.dofArray(jj));
%                               end
%                           end
            
                    % Store Interface Information
                   interfaceNodes=[interfaceNodes currentNodes(lowCoordPosition).getId];
               else
                   % Case 4 both nodes right of Boundary 
                   elementArrayRight = [elementArrayRight elementArray(ii)];
               end
           end
           interfaceNodes=unique(interfaceNodes);
           % fix Dofs
           %nodes=femM
       end
       
       %elementArrayLeft(find(elemenentArrayLeft==0))=[];
   end
   
   %%% END---Substructure_1
   
end

