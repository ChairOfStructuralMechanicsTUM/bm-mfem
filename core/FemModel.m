classdef FemModel < handle
    %FEMMODEL The class of the complete model
    %   This class keeps track over all entities in the model 
    
    properties (Access = private)
        nodeArray
        elementArray
        dofArray = {}
        
    end
    
    methods
        % constructor
        function loadFemModel(femModel, nodeArray, elementArray)
           
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
            
        end
        
        % getter functions
        function nodeArray = getNodeArray(femModel)
           nodeArray = femModel.nodeArray; 
        end
        
        function elementArray = getElementArray(femModel)
           elementArray = femModel.elementArray; 
        end
        
        function dofArray = getDofArray(femModel)
           dofArray = femModel.dofArray; 
        end
        
    end
    
end

