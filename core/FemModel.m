classdef FemModel < handle
    %FEMMODEL The class of the complete model
    %   femModel = FemModel(nodeArray, elementArray, femModelParts)
    %   This class keeps track over all entities in the model
    
    properties (Access = private) %changed from private to public
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
        
        function [substructure1, substructure2] = divide(femModel, eleIntf)
            %all nodes/elements of Geomtrie
            totalNodeArray = femModel.getAllNodes;
            totalElementArray = femModel.getAllElements;
            %all nodes at Interface
            nodeIntf = eleIntf.getNodes();
            %node furthest to the right of Interface elements
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
            
            nodesRight = [nodesRight, nodesAtX(totalNodeArray, maxX)];
            
            %all elements left or at interface
            elementsLeft = unique(findElements(nodesLeft, totalElementArray));
           
            %all elements right of interface
            elementsRight = unique(findElements(nodesRight, totalElementArray));
            
            %make area of interface elements half
            
            %create Substructure from NodeArray and ElementArray
            substructure1 = Substructure(nodesLeft, elementsLeft);
            substructure2 = Substructure(nodesRight, elementsRight);
            
            %%%End New
        end        
    end
end

