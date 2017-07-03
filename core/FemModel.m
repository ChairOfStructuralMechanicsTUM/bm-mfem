classdef FemModel < handle
    %FEMMODEL The class of the complete model
    %   femModel = FemModel(nodeArray, elementArray, femModelParts)
    %   This class keeps track over all entities in the model
    
    properties (Access = private) 
        nodeArray
        elementArray
        %         dofArray = {}
        femModelParts = containers.Map
        
    end
    
    properties (Access = public)
        dofArray = {}
    end
    
    methods
        % constructor
        function femModel = FemModel(nodeArray, elementArray, femModelParts)
           
            nodeIds = arrayfun(@(node) node.getId, nodeArray);
            if (hasDuplicates(nodeIds))
                error('problem with node ids');
            else
                femModel.nodeArray = nodeArray;
            end
            
            femModel.dofArray = arrayfun(@(node) node.getDofArray, nodeArray, 'UniformOutput', false);
            femModel.dofArray = [femModel.dofArray{:}];
            
            elementIds = arrayfun(@(element) element.getId, elementArray);
            if (hasDuplicates(elementIds))
                error('problem with element ids');
            else
                femModel.elementArray = elementArray;
            end
            
            if nargin == 3
                femModel.femModelParts = femModelParts;
            end
            
        end
        
        % getter functions
        function nodeArray = getAllNodes(femModel)
           nodeArray = femModel.nodeArray; 
        end
        
        function elementArray = getAllElements(femModel)
           elementArray = femModel.elementArray; 
        end
        
        function dofArray = getDofArray(femModel)
           dofArray = femModel.dofArray; 
        end
        
        function node = getNode(femModel, id)
           node = femModel.nodeArray(id); 
        end
        
        function nodes = getNodes(femModel, ids)
            nodes = Node.empty;
            for ii = 1:length(ids)
                nodes(ii) = femModel.nodeArray(ids(ii));
            end
        end
        
        function element = getElement(femModel, id)
            element = femModel.elementArray(id);
        end
        
        function elements = getElements(femModel, ids)
           elements = Element.empty;
           for ii = 1:length(ids)
              elements(ii) = femModel.elementArray(ids(ii)); 
           end
        end
        
        function femModelParts = getAllModelParts(femModel)
            femModelParts = femModel.femModelParts;
        end
        
        function modelPart = getModelPart(femModel, name)
            modelPart = femModel.femModelParts(name);
        end
        
        % setter functions
        function addModelPart(femModel, name, entityArray)
            femModel.femModelParts(name) = entityArray;
        end
        
        %%%Start New
        %function to divide FemModel into Substructures.
        %FemModel: structure to be divided
        %eleIntf: elements that are on the Interface
        function [substructure01, substructure02] = divide(femModel, eleIntf)
            %all nodes/elements of Geomtrie
            totalNodeArray = femModel.getAllNodes;
            totalElementArray = femModel.getAllElements;
            maxEleId = max(totalElementArray.getId)+1;
            
            %half the crossSectionArea 
            halfCrossSectionArea(eleIntf);
            
            %all nodes at Interface 
            nodeIntfNoSort = findNodes(eleIntf);

            %sort nodes at interface by Id
            nodeIntf = sortNodes(nodeIntfNoSort);

            %half the point loads at interface nodes before copying
            halfPointLoads(nodeIntf);
            
            %see whether split is in x- or y-direction
            orientationX = unique(nodeIntf.getX());
            orientationY = unique(nodeIntf.getY());
            
            if size(orientationX) == 1
                [nodes01, nodes02] = splitNodesX(nodeIntf, totalNodeArray);
        
                %all elements that need to be copied
                elementsToCopy = elementsForCopy(nodeIntf, nodes02, totalElementArray);
                %all elements that have been copied
                copiedElements = callToCopy(elementsToCopy, nodes01, nodes02, maxEleId, eleIntf);
                %copied elements are added
                totalElementArray = [totalElementArray, copiedElements];
                [elements01, elements02] = splitElementsX(nodes01, nodes02, totalElementArray);
               
            elseif size(orientationY) == 1
                [nodes01, nodes02]= splitNodesY(nodeIntf, totalNodeArray);
                
                %all elements that need to be copied
                elementsToCopy = elementsForCopy(nodeIntf, nodes02, totalElementArray);
                %all elements that have been copied
                copiedElements = callToCopy(elementsToCopy, nodes01, nodes02, maxEleId, eleIntf);
                %copied elements are added
                totalElementArray = [totalElementArray, copiedElements];
                [elements01, elements02] = splitElementsY(nodes01, nodes02, totalElementArray);
            
            else
                disp('Chosen Elements are not in a vertical or horizontal line');
                return;
            end
            
            %give interface nodes to substructures
            nodeIntf01 = nodeIntf;
            nodeIntf02 = [];
            for ii = 1:length(nodeIntf01)
                nodeIntf02 = [nodeIntf02, findCopy(nodeIntf01(ii), nodes02)];
            end
            
            %create Substructure from NodeArray and ElementArray
            substructure01 = Substructure(nodes01, elements01, nodeIntf01);
            substructure02 = Substructure(nodes02, elements02, nodeIntf02);
        end 
        %%%End NEW
    end
end

